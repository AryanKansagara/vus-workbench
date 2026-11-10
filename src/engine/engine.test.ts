import { test } from 'node:test'
import assert from 'node:assert/strict'
import { formatC, formatP, lengthChange, parseC, parseP } from './hgvs.ts'
import { revcomp, translate } from './seq.ts'
import { annotate, annotateHgvs, cToG, gToC, gToHgvsC, hgvsCToG, proteinEffect, type GVariant } from './transcript.ts'
import { autoMap, detectKind, labVariantKey, LAB_FIELDS, normalize, parseCsv, parseGtf, parseVcf, spliceAiMax, splitGt, validateLabCsv, validateVariantCsv, VARIANT_FIELDS } from './io.ts'
import { classify, combine, DEFAULT_GENE, EMPTY_EVIDENCE, oddsPathStrength, parseClass, parseCriteria, pvs1, revelStrength, suggest, vusTier } from './acmg.ts'
import { analyse, evaluateAll, phaseOf } from './pipeline.ts'
import { reconcile, severity } from './discrepancy.ts'
import { acceptorScore, donorScore, DONOR_MAX } from './splice.ts'
import { BUNDLED_GENES, BUNDLED_GENOME as G, BUNDLED_TRANSCRIPTS, LAB_CSV, REAL_VARIANTS_CSV, SAMPLE_GTF, sampleVcf, synthRef } from './sample.ts'

const [T1, T2] = BUNDLED_TRANSCRIPTS
let n = 0
const id = () => `id${n++}`
const hc = (s: string) => { const h = parseC(s); assert.ok(h, s); return h! }
const evAll = () => {
  const r = parseVcf(sampleVcf(), 's.vcf', G, id)
  return evaluateAll(r.variants.map(v => analyse(v, BUNDLED_TRANSCRIPTS, G)), {}, BUNDLED_GENES)
}
const byC = (list: ReturnType<typeof evAll>, c: string) => list.find(e => e.an.ann?.hgvsC.endsWith(':' + c))!

test('HGVS c. parser covers substitutions, intronic/UTR positions and indels, and round-trips', () => {
  const cases = ['NM_000059.4(BRCA2):c.68-7T>A', 'NM_007294.4(BRCA1):c.5266dup', 'NM_000492.4(CFTR):c.1521_1523del', 'SYN_000101.1(SYNG1):c.*35C>A', 'SYN_000101.1(SYNG1):c.-12G>T', 'NM_1.1:c.100_101insTT', 'NM_1.1:c.100_102delinsAG', 'NM_1.1:c.306+5G>A']
  for (const c of cases) assert.equal(formatC(hc(c)), c)
  assert.equal(formatC(hc('c.68-7T>A')), 'c.68-7T>A')
  const h = hc('BRCA2(NM_000059.4):c.68-7T>A')
  assert.equal(h.transcript, 'NM_000059.4'); assert.equal(h.gene, 'BRCA2')
  assert.deepEqual(hc('BRCA2:c.9976A>T').gene, 'BRCA2')
  assert.equal(formatC(hc('NM_1.1:c.76_78delACT')), 'NM_1.1:c.76_78del')
  assert.equal(parseC('c.100_101ins'), null)
  assert.equal(parseC('c.12A>Z'), null)
  assert.equal(parseC('g.12345A>G'), null)
})

test('length change counts across the c.-1/c.1 boundary and intronic spans', () => {
  assert.equal(lengthChange(hc('c.-1_2del')), -3)
  assert.equal(lengthChange(hc('c.5266dup')), 1)
  assert.equal(lengthChange(hc('c.1521_1523del')), -3)
  assert.equal(lengthChange(hc('c.100+3_100+6del')), -4)
  assert.equal(lengthChange(hc('c.100_102delinsAG')), -1)
})

test('HGVS p. parser accepts 3-letter, 1-letter, X and fs forms', () => {
  assert.deepEqual(parseP('p.(Gln35Ter)')?.kind, 'nonsense')
  assert.equal(parseP('p.R123X')?.kind, 'nonsense')
  assert.equal(parseP('p.Arg1699Gln')?.kind, 'missense')
  assert.equal(parseP('p.(Ala90=)')?.kind, 'synonymous')
  const fs = parseP('p.(Gln1756ProfsTer74)')!
  assert.equal(fs.kind, 'frameshift'); assert.equal(fs.pos, 1756); assert.equal(fs.fsTer, 74)
  assert.equal(parseP('p.Phe508del')?.kind, 'inframe_del')
  assert.equal(parseP('p.(Met1?)')?.kind, 'start_lost')
  assert.equal(formatP(parseP('p.R123*')!), 'p.Arg123Ter')
})

test('translation and reverse complement', () => {
  assert.equal(translate('ATGTTTTAAGGG'), 'MF*G')
  assert.equal(translate('ATGTTTTAAGGG', true), 'MF*')
  assert.equal(revcomp('AACGTN'), 'NACGTT')
})

test('genomic <-> c. mapping round-trips on plus and minus strand models', () => {
  for (const [t, cs] of [[T1, ['c.156+7', 'c.157-12', 'c.306+250']], [T2, ['c.50+7', 'c.51-12', 'c.162+3']]] as const) {
    for (const c of ['c.1', 'c.60', 'c.61', 'c.-40', 'c.*10', ...cs]) {
      const h = hc(`${t.id}:${c}A>G`)
      const g = cToG(t, h.start)
      assert.deepEqual(gToC(t, g).cpos, h.start, `${t.gene} ${c}`)
    }
  }
  // Minus strand: c. increases as genomic coordinates decrease.
  assert.equal(cToG(T2, { pos: 2, offset: 0 }), cToG(T2, { pos: 1, offset: 0 }) - 1)
})

test('intronic positions anchor to the nearer exon, ties to the donor side', () => {
  // SYNG1 intron 3 is 500 nt: position 250 from both ends is a tie.
  const donor = cToG(T1, { pos: 306, offset: 0 })
  assert.deepEqual(gToC(T1, donor + 250).cpos, { pos: 306, offset: 250 })
  assert.deepEqual(gToC(T1, donor + 251).cpos, { pos: 307, offset: -250 })
})

test('VCF normalisation left-aligns an indel written mid-homopolymer, and HGVS shifts it 3\' in the transcript', () => {
  const g = cToG(T1, { pos: 150, offset: 0 })
  const v = normalize(G, { chrom: 'chrS1', pos: g, ref: 'AA', alt: 'A' })
  assert.equal(v.leftAligned, true)
  assert.equal(v.ref.length, 2); assert.equal(v.alt.length, 1)
  assert.equal(v.pos, cToG(T1, { pos: 147, offset: 0 }))
  assert.equal(formatC(gToHgvsC(T1, G, v), false), 'c.152del')
  // Without reference only trimming happens.
  const u = normalize({}, { chrom: '7', pos: 10, ref: 'CAT', alt: 'CGT' })
  assert.deepEqual([u.pos, u.ref, u.alt, u.leftAligned], [11, 'A', 'G', false])
})

test('insertions that repeat the preceding bases are described as dup', () => {
  const v = synthRef('SYN_000101.1:c.560dup')!
  assert.equal(formatC(gToHgvsC(T1, G, v), false), 'c.560dup')
  assert.match(proteinEffect(T1, G, v)!.hgvsP, /^p\.\(Ser188ValfsTer48\)$/)
})

test('VCF import splits multi-allelic sites, re-expresses genotypes, and reads INFO per allele', () => {
  const r = parseVcf(sampleVcf(), 's.vcf', G, id)
  const site = r.variants.filter(v => v.source.line === r.variants.find(x => x.gt === '0/1' && x.evidence.revel === 0.52)?.source.line)
  assert.equal(site.length, 2)
  assert.deepEqual(site.map(s => s.gt), ['1/0', '0/1'])
  assert.deepEqual(site.map(s => s.evidence.revel), [0.81, 0.52])
  assert.equal(splitGt('1|2', 2), '0|1')
  assert.ok(r.issues.some(i => /Multi-allelic/.test(i.message)))
  assert.ok(r.issues.some(i => /left-aligned/.test(i.message)))
})

test('VCF validation reports problems with line numbers and skips bad records', () => {
  const pos = cToG(T1, { pos: 103, offset: 0 })
  const text = ['##fileformat=VCFv4.2', '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', `chrS1\tabc\t.\tC\tT\t.\tPASS\t.`, `chrS1\t${pos}\t.\tG\tT\t.\tPASS\t.`, 'chrS1\t5', `chrS1\t${pos}\t.\tC\t<DEL>\t.\tPASS\t.`, `chrS1\t${pos}\t.\tC\tT\t.\tLowQual\t.`].join('\n')
  const r = parseVcf(text, 'bad.vcf', G, id)
  const at = (l: number) => r.issues.filter(i => i.line === l).map(i => i.level)
  assert.deepEqual(at(3), ['error'])
  assert.deepEqual(at(4), ['error'])
  assert.deepEqual(at(5), ['error'])
  assert.deepEqual(at(6), ['warning'])
  assert.equal(r.variants.length, 1)
  assert.equal(r.variants[0].filter, 'LowQual')
})

test('SpliceAI INFO is reduced to the maximum delta score for the allele', () => {
  assert.equal(spliceAiMax('T|BRCA2|0.01|0.00|0.43|0.02|1|2|3|4,G|BRCA2|0.9|0|0|0|0|0|0|0', 'T'), 0.43)
  assert.ok(Number.isNaN(spliceAiMax('A|X|0.1|0|0|0|0|0|0|0', 'T')))
})

test('protein effects come from translating the edited transcript', () => {
  const e = (c: string) => proteinEffect(c.startsWith('SYN_0002') ? T2 : T1, G, synthRef(c)!)!
  assert.equal(e('SYN_000101.1:c.103C>T').hgvsP, 'p.(Gln35Ter)')
  assert.equal(e('SYN_000101.1:c.208C>T').hgvsP, 'p.(Arg70Trp)')
  assert.equal(e('SYN_000101.1:c.270C>T').hgvsP, 'p.(Ala90=)')
  assert.equal(e('SYN_000101.1:c.130_132del').hgvsP, 'p.(Asp44del)')
  // Minus-strand gene.
  assert.equal(e('SYN_000202.1:c.100G>A').hgvsP, 'p.(Gly34Ser)')
})

test('NMD: an early stop triggers decay, a stop in the last exon escapes it', () => {
  const list = evAll()
  assert.equal(byC(list, 'c.103C>T').an.ann!.nmd, true)
  assert.equal(byC(list, 'c.560dup').an.ann!.nmd, false)
  assert.equal(byC(list, 'c.103C>T').an.ann!.ptcAa, 35)
})

test('PVS1 decision tree: NMD -> very strong, last exon -> by protein removed, in-frame skip -> strong', () => {
  const list = evAll()
  const strength = (c: string) => byC(list, c).criteria.find(x => x.code === 'PVS1')?.strength
  assert.equal(strength('c.103C>T'), 'VeryStrong')
  assert.equal(strength('c.560dup'), 'Moderate')
  assert.equal(strength('c.306+1G>A'), 'Strong')
  assert.equal(byC(list, 'c.306+1G>A').an.ann!.skipInFrame, true)
  const a = byC(list, 'c.103C>T').an.ann!
  assert.equal(pvs1(a, { ...DEFAULT_GENE, lofMechanism: false }).strength, null)
  assert.equal(pvs1({ ...a, consequence: 'start_lost' }, DEFAULT_GENE).strength, 'Moderate')
})

test('Bayesian point thresholds and the stand-alone BA1 rule', () => {
  assert.deepEqual([10, 9, 6, 5, 0, -1, -6, -7].map(classify), ['Pathogenic', 'Likely pathogenic', 'Likely pathogenic', 'VUS', 'VUS', 'Likely benign', 'Likely benign', 'Benign'])
  assert.equal(combine([{ code: 'PVS1', strength: 'VeryStrong' }, { code: 'PM2', strength: 'Supporting' }, { code: 'PP3', strength: 'Supporting' }]).points, 9)
  assert.equal(combine([{ code: 'PVS1', strength: 'VeryStrong' }, { code: 'PS3', strength: 'Strong' }, { code: 'BA1', strength: 'StandAlone' }]).cls, 'Benign')
  assert.equal(combine([{ code: 'PS3', strength: 'Strong' }, { code: 'PM2', strength: 'Supporting' }, { code: 'PM4', strength: 'Moderate' }, { code: 'BS3', strength: 'Supporting' }]).points, 6)
  assert.equal(vusTier(4), 'hot'); assert.equal(vusTier(6), null)
})

test('double counting: PP3 is not combined with PVS1, and PP3/BP4 conflicts are excluded', () => {
  const c = combine([{ code: 'PVS1', strength: 'VeryStrong' }, { code: 'PP3', strength: 'Supporting' }])
  assert.equal(c.points, 8)
  assert.equal(c.excluded[0].code, 'PP3')
  assert.equal(combine([{ code: 'PP3', strength: 'Supporting' }, { code: 'BP4', strength: 'Supporting' }]).points, 0)
})

test('calibrated REVEL and OddsPath thresholds', () => {
  assert.deepEqual(revelStrength(0.95), { code: 'PP3', strength: 'Strong' })
  assert.deepEqual(revelStrength(0.8), { code: 'PP3', strength: 'Moderate' })
  assert.deepEqual(revelStrength(0.65), { code: 'PP3', strength: 'Supporting' })
  assert.equal(revelStrength(0.5), null)
  assert.deepEqual(revelStrength(0.1), { code: 'BP4', strength: 'Supporting' })
  assert.deepEqual(oddsPathStrength(20), { code: 'PS3', strength: 'Strong' })
  assert.deepEqual(oddsPathStrength(5), { code: 'PS3', strength: 'Moderate' })
  assert.deepEqual(oddsPathStrength(0.04), { code: 'BS3', strength: 'Strong' })
  assert.equal(oddsPathStrength(1), null)
})

test('population frequency drives PM2 (Supporting), BS1 and BA1 from gene thresholds', () => {
  const a = annotateHgvs(hc('BRCA2:c.7522G>A'), 'p.(Gly2508Ser)')
  const codes = (af: number) => suggest(a, { ...EMPTY_EVIDENCE, popmaxAf: af }, DEFAULT_GENE, null).map(s => `${s.code}_${s.strength}`)
  assert.deepEqual(codes(0), ['PM2_Supporting'])
  assert.deepEqual(codes(0.02), ['BS1_Strong'])
  assert.deepEqual(codes(0.2), ['BA1_StandAlone'])
  assert.deepEqual(codes(0.001), [])
})

test('functional evidence: strongest PS3 entry by OddsPath, Supporting when no OddsPath', () => {
  const a = annotateHgvs(hc('BRCA2:c.7522G>A'), 'p.(Gly2508Ser)')
  const f = (oddsPath: number | null, result: 'abnormal' | 'normal' = 'abnormal') => ({ id: id(), assay: 'HDR assay', result, oddsPath, note: '', editor: 'qa', at: 0 })
  const s = suggest(a, { ...EMPTY_EVIDENCE, functional: [f(null), f(25)] }, DEFAULT_GENE, null)
  assert.equal(s.find(x => x.code === 'PS3')?.strength, 'Strong')
  assert.equal(suggest(a, { ...EMPTY_EVIDENCE, functional: [f(null, 'normal')] }, DEFAULT_GENE, null)[0].strength, 'Supporting')
})

test('phase from genotypes and PM3 / BP2 for two variants in one recessive gene', () => {
  const v = (gt: string, ps?: string) => ({ id: id(), key: '', gt, ps, source: { file: '', line: 0, kind: 'vcf' as const }, evidence: EMPTY_EVIDENCE })
  assert.equal(phaseOf(v('0|1'), v('1|0')), 'trans')
  assert.equal(phaseOf(v('0|1'), v('0|1')), 'cis')
  assert.equal(phaseOf(v('0/1'), v('0|1')), 'unknown')
  assert.equal(phaseOf(v('0|1', '5'), v('1|0', '9')), 'unknown')
  const list = evAll()
  assert.equal(byC(list, 'c.100G>A').criteria.find(c => c.code === 'PM3')?.strength, 'Moderate')
  assert.ok(byC(list, 'c.250C>T').criteria.some(c => c.code === 'BP2'))
  assert.equal(byC(list, 'c.100G>A').result.points, 4)
})

test('local splice model: consensus donor scores highest, +5G>A weakens and +1G>A abolishes the site', () => {
  assert.equal(donorScore('CAGGTAAGT') > donorScore('CAGGTAAAT'), true)
  assert.ok(Math.abs(donorScore('CAGGTAAGT') - DONOR_MAX) < 1)
  assert.ok(acceptorScore('TTTCTCTTTCCTCAGG') - acceptorScore('TTTCTCTTTCCTCAAG') > 8)
  assert.ok(acceptorScore('TTTCTCTTTCCTCAGG') > acceptorScore('AGAGAGAGAGAGCAGG'))
  const list = evAll()
  assert.equal(byC(list, 'c.306+5G>A').an.splice?.verdict, 'weakens')
  assert.equal(byC(list, 'c.306+1G>A').an.splice?.verdict, 'disrupts')
  assert.ok(byC(list, 'c.306+5G>A').criteria.some(c => c.code === 'PP3' && /PWM/.test(c.evidence)))
})

test('HGVS-only annotation (no transcript model) still classifies position and type', () => {
  assert.equal(annotateHgvs(hc('BRCA2:c.68-7T>A')).consequence, 'splice_region')
  assert.equal(annotateHgvs(hc('BRCA2:c.68-2A>G')).consequence, 'splice_acceptor')
  assert.equal(annotateHgvs(hc('BRCA1:c.5266dup'), 'p.(Gln1756ProfsTer74)').ptcAa, 1829)
  assert.equal(annotateHgvs(hc('BRCA1:c.5266dup')).consequence, 'frameshift')
  assert.equal(annotateHgvs(hc('BRCA1:c.100A>G')).consequence, 'coding_sequence')
})

test('a last-exon stop without a model is sized from the protein length and lands at PVS1_Moderate', () => {
  const csv = parseCsv(REAL_VARIANTS_CSV)
  const r = validateVariantCsv(csv, autoMap(csv.rows[0], VARIANT_FIELDS), 'x.csv', G, id)
  assert.equal(r.issues.length, 0)
  const ev = evaluateAll(r.variants.map(v => analyse(v, BUNDLED_TRANSCRIPTS, G)), {}, BUNDLED_GENES)
  const k = ev.find(e => e.an.rec.hgvsC === 'c.9976A>T')!
  assert.equal(k.criteria.find(c => c.code === 'PVS1')?.strength, 'Moderate')
})

test('variant CSV validation flags bad HGVS, unknown genotypes and out-of-range scores by line', () => {
  const csv = parseCsv('gene,hgvs_c,revel,zygosity\nBRCA2,c.68-7T>A,1.7,het\nBRCA2,c.68-7T>Q,,\nBRCA2,c.100A>G,,maybe\n')
  const r = validateVariantCsv(csv, autoMap(csv.rows[0], VARIANT_FIELDS), 'x.csv', G, id)
  assert.equal(r.variants.length, 2)
  assert.deepEqual(r.issues.map(i => [i.line, i.level]), [[2, 'warning'], [3, 'error'], [4, 'warning']])
  assert.equal(r.variants[0].gt, '0/1')
})

test('lab criteria and class strings parse in the shorthand labs actually use', () => {
  const p = parseCriteria('PVS1, PM2_Supporting; PP3_Mod PS3_P BS3_Str foo')
  assert.deepEqual(p.applied.map(a => `${a.code}_${a.strength}`), ['PVS1_VeryStrong', 'PM2_Supporting', 'PP3_Moderate', 'PS3_Supporting', 'BS3_Strong'])
  assert.deepEqual(p.unknown, ['foo'])
  assert.equal(parseClass('Uncertain significance'), 'VUS')
  assert.equal(parseClass('likely_pathogenic'), 'Likely pathogenic')
  assert.equal(parseClass('maybe'), null)
})

test('discrepancy reconciliation: severity, criteria spread, and labs whose criteria do not add up', () => {
  assert.equal(severity(['Pathogenic', 'Benign']), 'major')
  assert.equal(severity(['Likely pathogenic', 'VUS']), 'clinically-significant')
  assert.equal(severity(['Pathogenic', 'Likely pathogenic']), 'confidence')
  const csv = parseCsv(LAB_CSV)
  const r = validateLabCsv(csv, autoMap(csv.rows[0], LAB_FIELDS), 'labs.csv', id)
  const d = reconcile(r.rows)
  assert.equal(d[0].variant, 'BRCA2:c.9976A>T')
  assert.equal(d[0].severity, 'major')
  const k = d.find(x => x.variant === 'BRCA2:c.7522G>A')!
  assert.equal(k.calls.find(c => c.lab === 'Synthetic Lab A')!.mismatch, true)
  assert.equal(k.diffs[0].code, 'PS3')
  assert.equal(labVariantKey('NM_000059.4:c.68-7T>A', 'BRCA2'), 'BRCA2:c.68-7T>A')
  assert.equal(labVariantKey('c.68-7T>A', 'BRCA2'), 'BRCA2:c.68-7T>A')
})

test('GTF import rebuilds the bundled transcript models exactly', () => {
  const r = parseGtf(SAMPLE_GTF, 'test')
  assert.equal(r.issues.length, 0)
  for (const t of BUNDLED_TRANSCRIPTS) {
    const g = r.transcripts.find(x => x.id === t.id)!
    assert.deepEqual([g.gene, g.strand, g.exons, g.cdsStart, g.cdsEnd], [t.gene, t.strand, t.exons, t.cdsStart, t.cdsEnd])
  }
  assert.equal(detectKind(SAMPLE_GTF), 'gtf')
  assert.equal(detectKind(sampleVcf()), 'vcf')
  assert.equal(detectKind(LAB_CSV), 'lab-csv')
  assert.equal(detectKind(REAL_VARIANTS_CSV), 'variant-csv')
})

test('c. to genomic rejects a reference mismatch (wrong transcript version)', () => {
  const ok = hgvsCToG(T1, G, hc('SYN_000101.1:c.103C>T'))
  assert.ok(ok)
  assert.equal(hgvsCToG(T1, G, hc('SYN_000101.1:c.103G>T')), null)
  const v: GVariant = ok!
  const a = annotate(T1, G, gToHgvsC(T1, G, v), v)
  assert.equal(a.consequence, 'stop_gained')
  assert.equal(a.exon, 2)
})
