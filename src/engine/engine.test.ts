import { test } from 'node:test'
import assert from 'node:assert/strict'
import { formatC, formatP, lengthChange, parseC, parseP } from './hgvs.ts'
import { revcomp, translate } from './seq.ts'
import { annotate, annotateHgvs, cdsBounds, cToG, gToC, gToHgvsC, hasReference, hgvsCToG, proteinEffect, refAt, type Genome, type GVariant, type Transcript } from './transcript.ts'
import { autoMap, detectKind, identities, inIndex, labVariantKey, LAB_FIELDS, normalize, parseCsv, parseFasta, parseGtf, parseVcf, spliceAiMax, splitGt, validateLabCsv, validateVariantCsv, VARIANT_FIELDS, variantIndex, type VariantRec } from './io.ts'
import { classify, combine, DEFAULT_GENE, EMPTY_EVIDENCE, oddsPathStrength, parseClass, parseCriteria, posterior, pvs1, revelStrength, suggest, vusTier, type Applied, type Strength } from './acmg.ts'
import { readFileSync } from 'node:fs'
import { clinvarLabRows, clinvarTerms, criteriaFromComment, parseVcvXml, titleMatches } from './clinvar.ts'
import { transcriptFromEnsembl } from './ensembl.ts'
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
  assert.equal(vusTier(5), 'hot'); assert.equal(vusTier(6), null)
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
  assert.deepEqual(revelStrength(0.25), { code: 'BP4', strength: 'Supporting' })
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

// ---------- Published references and real public data ----------

const fixture = (name: string) => readFileSync(new URL(`./fixtures/${name}`, import.meta.url), 'utf8')
const pts = (spec: Partial<Record<'vs' | 's' | 'm' | 'p', number>>, dir: 'P' | 'B' = 'P'): Applied[] => {
  const codes = dir === 'P' ? { vs: ['PVS1'], s: ['PS1', 'PS2', 'PS3', 'PS4'], m: ['PM1', 'PM2', 'PM3', 'PM5', 'PM6'], p: ['PP1', 'PP2', 'PP4', 'PP5'] } : { vs: [], s: ['BS1', 'BS2', 'BS3', 'BS4'], m: [], p: ['BP1', 'BP2', 'BP5', 'BP7'] }
  const str: Record<string, Strength> = { vs: 'VeryStrong', s: 'Strong', m: 'Moderate', p: 'Supporting' }
  return (Object.entries(spec) as [keyof typeof codes, number][]).flatMap(([k, n]) => codes[k].slice(0, n).map(code => ({ code, strength: str[k] })))
}

test('point thresholds are the Bayesian posteriors: prior 0.10, very strong odds 350 (Tavtigian 2018, 2020)', () => {
  const near = (a: number, b: number, tol = 5e-4) => assert.ok(Math.abs(a - b) < tol, `${a} vs ${b}`)
  near(posterior(0), 0.10)
  near(posterior(8), 350 * 0.1 / (349 * 0.1 + 1))
  // Pathogenic needs > 0.99, likely pathogenic >= 0.90 (6 points gives 0.8999, the value the framework accepts as 0.90).
  assert.ok(posterior(10) > 0.99 && posterior(9) < 0.99)
  near(posterior(6), 0.90)
  assert.ok(posterior(5) < 0.895)
  // Likely benign needs < 0.10, benign < 0.001.
  assert.ok(posterior(-1) < 0.10 && posterior(-6) > 0.001 && posterior(-7) < 0.001)
  for (let p = -12; p <= 14; p++) {
    const post = posterior(p)
    const expected = post > 0.99 ? 'Pathogenic' : post >= 0.8995 ? 'Likely pathogenic' : post < 0.001 ? 'Benign' : post < 0.10 - 1e-9 ? 'Likely benign' : 'VUS'
    assert.equal(classify(p), expected, `${p} points`)
  }
})

test('the 2015 ACMG/AMP combining rules land in the same class under points', () => {
  const cls = (a: Applied[]) => combine(a).cls
  // Pathogenic: PVS1 with one strong, or with a moderate and a supporting; one strong with three moderate.
  assert.equal(cls(pts({ vs: 1, s: 1 })), 'Pathogenic')
  assert.equal(cls(pts({ vs: 1, m: 1, p: 1 })), 'Pathogenic')
  assert.equal(cls(pts({ s: 1, m: 3 })), 'Pathogenic')
  assert.equal(cls(pts({ s: 1, m: 2, p: 2 })), 'Pathogenic')
  // Likely pathogenic: strong + moderate, strong + two supporting, three moderate, two moderate + two supporting, one moderate + four supporting.
  for (const spec of [{ s: 1, m: 1 }, { s: 1, p: 2 }, { m: 3 }, { m: 2, p: 2 }, { m: 1, p: 4 }]) assert.equal(cls(pts(spec)), 'Likely pathogenic', JSON.stringify(spec))
  // One short of each stays uncertain.
  for (const spec of [{ s: 1, p: 1 }, { m: 2, p: 1 }, { m: 1, p: 3 }]) assert.equal(cls(pts(spec)), 'VUS', JSON.stringify(spec))
  // Benign side: strong + supporting or two supporting is likely benign; two strong is benign.
  assert.equal(cls(pts({ s: 1, p: 1 }, 'B')), 'Likely benign')
  assert.equal(cls(pts({ p: 2 }, 'B')), 'Likely benign')
  assert.equal(cls(pts({ s: 2 }, 'B')), 'Benign')
  // Conflicting evidence nets out: one strong each way is 0 points.
  assert.equal(combine([...pts({ s: 1 }), ...pts({ s: 1 }, 'B')]).points, 0)
})

test('VUS temperature is one level per point, matching the ACGS posterior bands', () => {
  assert.deepEqual([0, 1, 2, 3, 4, 5].map(vusTier), ['ice cold', 'cold', 'cool', 'tepid', 'warm', 'hot'])
  assert.equal(vusTier(-1), null)
  // Band edges 10 / 18.8 / 32.5 / 50 / 67.5 / 81.2 / 90 %.
  assert.deepEqual([0, 1, 2, 3, 4, 5].map(p => Math.round(posterior(p) * 1000) / 10), [10, 18.8, 32.5, 50, 67.5, 81.2])
})

test('REVEL intervals follow Pejaver 2022 Table 2 at every boundary', () => {
  const at = (x: number) => { const r = revelStrength(x); return r ? `${r.code}_${r.strength}` : null }
  assert.deepEqual([0.003, 0.0031, 0.016, 0.0161, 0.183, 0.1831, 0.290, 0.2901].map(at), ['BP4_VeryStrong', 'BP4_Strong', 'BP4_Strong', 'BP4_Moderate', 'BP4_Moderate', 'BP4_Supporting', 'BP4_Supporting', null])
  assert.deepEqual([0.6439, 0.644, 0.7729, 0.773, 0.9319, 0.932].map(at), [null, 'PP3_Supporting', 'PP3_Supporting', 'PP3_Moderate', 'PP3_Moderate', 'PP3_Strong'])
  // A low score now carries its calibrated weight in the total.
  const a = annotateHgvs(hc('BRCA2:c.7522G>A'), 'p.(Gly2508Ser)')
  const bp4 = suggest(a, { ...EMPTY_EVIDENCE, revel: 0.1 }, DEFAULT_GENE, null).find(x => x.code === 'BP4')!
  assert.equal(bp4.strength, 'Moderate')
  assert.equal(combine([{ code: 'BP4', strength: bp4.strength }]).points, -2)
})

test('PP3 and PM1 together are limited to 4 points', () => {
  const c = combine([{ code: 'PP3', strength: 'Strong' }, { code: 'PM1', strength: 'Moderate' }, { code: 'PM2', strength: 'Supporting' }])
  assert.equal(c.points, 5)
  assert.equal(c.counted.find(x => x.code === 'PP3')!.points, 2)
  assert.equal(c.notes.length, 1)
  const ok = combine([{ code: 'PP3', strength: 'Moderate' }, { code: 'PM1', strength: 'Moderate' }])
  assert.equal(ok.points, 4)
  assert.equal(ok.notes.length, 0)
})

test('ClinVar classification terms and criteria spellings', () => {
  assert.equal(parseClass('Uncertain_significance'), 'VUS')
  assert.equal(parseClass('Pathogenic/Likely_pathogenic'), 'Likely pathogenic')
  assert.equal(parseClass('Benign/Likely benign'), 'Likely benign')
  assert.equal(parseClass('Likely pathogenic, low penetrance'), 'Likely pathogenic')
  assert.equal(parseClass('Pathogenic; drug response'), 'Pathogenic')
  assert.equal(parseClass('Pathogenic|risk_factor'), 'Pathogenic')
  for (const t of ['Conflicting classifications of pathogenicity', 'not provided', 'drug response', 'risk factor', 'Pathogenic/Likely risk allele']) assert.equal(parseClass(t), null, t)
  const p = parseCriteria('PVS1_Strong, PM2_Supp; BP5_Very Strong, BA1_Stand-alone (PS3_VeryStrong) PVS1_Moderate PP1_M.')
  assert.deepEqual(p.applied.map(a => `${a.code}_${a.strength}`), ['PVS1_Strong', 'PM2_Supporting', 'BP5_VeryStrong', 'BA1_StandAlone', 'PS3_VeryStrong', 'PP1_Moderate'])
  assert.deepEqual(p.unknown, [])
})

test('VCF 4.3 example from the specification (hts-specs simple.vcf)', () => {
  const r = parseVcf(fixture('hts-specs-simple.vcf'), 'simple.vcf', {}, id)
  assert.equal(r.meta.fileformat, 'VCFv4.3')
  assert.deepEqual(r.samples, ['NA00001', 'NA00002', 'NA00003'])
  assert.equal(r.records, 5)
  // NA00001 is 0|0 at the first two sites and at the microsatellite; it carries 1|2 at the multi-allelic site and the first (deletion) allele of the microsatellite.
  assert.deepEqual(r.variants.map(v => [v.key, v.gt]), [['20:1110696:A:G', '1|0'], ['20:1110696:A:T', '0|1'], ['20:1234567:GTC:G', '0/1']])
  assert.equal(r.issues.filter(i => i.level === 'error').length, 0)
  assert.ok(r.issues.some(i => /3 alleles skipped: NA00001/.test(i.message)))
})

test('structural, breakend and gVCF alleles are skipped without errors', () => {
  const vcf = ['##fileformat=VCFv4.4', '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1',
    'chrA\t2\t.\tT\t<DEL>\t.\tPASS\tSVLEN=-3\tGT\t0/1', 'chrA\t5\t.\tG\tG[chrA:9[\t.\tPASS\t.\tGT\t0/1', 'chrA\t14\t.\tT\t.CCCCCCG\t.\tPASS\t.\tGT\t0/1',
    'chrA\t20\t.\tA\tG,<NON_REF>\t.\t.\t.\tGT\t0/1', 'chrA\t21\t.\tC\t<NON_REF>\t.\t.\tEND=90\tGT\t0/0', ''].join('\n')
  const r = parseVcf(vcf, 'sv.vcf', {}, id)
  assert.deepEqual(r.variants.map(v => v.key), ['A:20:A:G'])
  assert.equal(r.issues.filter(i => i.level === 'error').length, 0)
  assert.equal(r.issues.filter(i => /skipped; this workbench handles sequence variants/.test(i.message)).length, 3)
})

test('ClinVar GRCh38 VCF: gene from GENEINFO, aggregate and conflicting classifications become lab rows', () => {
  const r = parseVcf(fixture('clinvar-grch38-excerpt.vcf'), 'clinvar.vcf', {}, id)
  assert.equal(r.meta.annotation, 'ClinVar')
  assert.equal(r.variants.length, 8)
  assert.equal(r.issues.filter(i => i.level === 'error').length, 0)
  assert.equal(r.variants[0].gene, 'OR4F5')
  assert.equal(r.variants.find(v => v.key === '1:930200:G:A')!.gene, 'SAMD11')
  const rows = (k: string) => r.labRows.filter(x => x.variant === k).map(x => [x.lab, x.cls])
  assert.deepEqual(rows('1:69134:A:G'), [['ClinVar (criteria provided, single submitter)', 'Likely benign']])
  // CLNSIGCONF=Uncertain_significance_(1)|Likely_benign_(2)
  assert.deepEqual(rows('1:930200:G:A'), [['ClinVar: 1 submitter (VUS)', 'VUS'], ['ClinVar: 2 submitters (Likely benign)', 'Likely benign']])
  // A record with no classification is still a variant, with no lab row.
  assert.deepEqual(rows('1:139320:C:G'), [])
  assert.equal(reconcile(r.labRows).find(d => d.variant === '1:930200:G:A')!.severity, 'clinically-significant')
})

test('VEP CSQ: the MANE Select entry supplies gene, transcript, HGVS, MAX_AF, REVEL and SpliceAI per allele', () => {
  const r = parseVcf(fixture('vep-csq-excerpt.vcf'), 'vep.vcf', {}, id)
  assert.equal(r.meta.annotation, 'VEP CSQ')
  const [a, b, c] = r.variants
  assert.deepEqual([a.gene, a.transcript, a.hgvsC, a.hgvsP], ['BRCA2', 'ENST00000380152.8', 'c.7522G>A', 'p.Gly2508Ser'])
  assert.deepEqual(a.evidence, { popmaxAf: 0.0001307, revel: 0.68, spliceAi: 0.01, functional: [] })
  // Two ALT alleles, sample is 0/1: only the first is carried, and its empty score fields stay empty rather than zero.
  assert.equal(r.variants.length, 3)
  assert.deepEqual([b.hgvsC, b.hgvsP, b.evidence.popmaxAf, b.evidence.revel, b.evidence.spliceAi], ['c.9976A>T', 'p.Lys3326Ter', 0.01393, null, null])
  assert.deepEqual([c.gene, c.hgvsC, c.evidence.revel, c.evidence.spliceAi], ['BRCA1', 'c.5096G>A', 0.94, 0.02])
  // No transcript model is loaded, yet the worklist is annotated from the HGVS the file carries.
  const ev = evaluateAll(r.variants.map(v => analyse(v, BUNDLED_TRANSCRIPTS, G)), {}, BUNDLED_GENES)
  assert.deepEqual(ev.map(e => e.an.ann?.consequence), ['missense', 'stop_gained', 'missense'])
  assert.equal(ev[1].result.cls, 'Benign')
  assert.ok(ev[2].result.counted.some(x => x.code === 'PP3' && x.strength === 'Strong'))
})

test('SnpEff ANN: first coding annotation, with dupC / delCTT / p.Lys3326* spellings', () => {
  const r = parseVcf(fixture('snpeff-ann-excerpt.vcf'), 'snpeff.vcf', {}, id)
  assert.equal(r.meta.annotation, 'SnpEff ANN')
  assert.deepEqual(r.variants.map(v => [v.gene, v.transcript, v.hgvsC, v.hgvsP, v.gt]), [
    ['BRCA2', 'NM_000059.4', 'c.9976A>T', 'p.Lys3326Ter', '0/1'],
    ['BRCA1', 'NM_007294.4', 'c.5266dupC', 'p.Gln1756fs', '0/1'],
    ['CFTR', 'NM_000492.4', 'c.1521_1523delCTT', 'p.Phe508del', '1/1'],
  ])
  assert.ok(r.variants.every(v => v.evidence.popmaxAf === null && v.evidence.revel === null))
  const ev = r.variants.map(v => analyse(v, BUNDLED_TRANSCRIPTS, G))
  assert.deepEqual(ev.map(e => e.ann?.consequence), ['stop_gained', 'frameshift', 'inframe_deletion'])
})

// Real MANE Select models (Ensembl release 116, GRCh38). Exon sequence with 25 nt of flank is stored; introns are filled with N.
function realModel(gene: 'BRCA1' | 'BRCA2'): { tx: Transcript; genome: Genome } {
  const f = JSON.parse(fixture('ensembl-brca.json'))[gene] as { transcript: Transcript; pad: number; exonSeq: string[] }
  const t = f.transcript
  const start = t.exons[0][0] - f.pad, end = t.exons[t.exons.length - 1][1] + f.pad
  const seq = Array.from<string>({ length: end - start + 1 }).fill('N')
  t.exons.forEach(([a], i) => [...f.exonSeq[i]].forEach((b, k) => { seq[a - f.pad - start + k] = b }))
  return { tx: t, genome: { [`${t.chrom}:${start}-${end}`]: { chrom: t.chrom, start, seq: seq.join('') } } }
}
const hgvsRec = (gene: string, hgvsC: string, evidence = EMPTY_EVIDENCE): VariantRec => ({ id: id(), key: `${gene}:${hgvsC}`, gene, hgvsC, source: { file: 't', line: 1, kind: 'csv' }, evidence })

test('Ensembl gene lookup becomes a transcript model named by its MANE RefSeq match', () => {
  const t = transcriptFromEnsembl(JSON.parse(fixture('ensembl-lookup-brca2.json')), 116)
  assert.equal(t.id, 'NM_000059.4')
  assert.deepEqual([t.gene, t.chrom, t.strand, t.exons.length], ['BRCA2', '13', 1, 27])
  assert.equal(cdsBounds(t).aaLength, 3418)
  assert.match(t.source, /Ensembl 116 MANE Select ENST00000380152\.8 \(GRCh38\)/)
  assert.throws(() => transcriptFromEnsembl({ display_name: 'XIST', Transcript: [] }), /no protein-coding transcript/)
})

test('real BRCA2 model (plus strand): public ClinVar descriptions map to their GRCh38 coordinates and protein changes', () => {
  const { tx, genome } = realModel('BRCA2')
  assert.ok(hasReference(genome, tx))
  const list = evaluateAll(['c.9976A>T', 'c.7522G>A', 'c.68-7T>A'].map(c => analyse(hgvsRec('BRCA2', c), [tx], genome)), {}, BUNDLED_GENES)
  const [stop, mis, intr] = list
  // ClinVar VCV000038266: NC_000013.11:g.32398489A>T, p.Lys3326Ter.
  assert.deepEqual(stop.an.g, { chrom: '13', pos: 32398489, ref: 'A', alt: 'T' })
  assert.equal(stop.an.ann!.hgvsP, 'p.(Lys3326Ter)')
  // The stop is in the last exon (27 of 27), escapes NMD and removes under 10% of the protein outside the DNA-binding domain: PVS1 at Moderate only.
  assert.deepEqual([stop.an.ann!.exon, stop.an.ann!.nmd], [27, false])
  assert.equal(stop.criteria.find(c => c.code === 'PVS1')!.strength, 'Moderate')
  // ClinVar VCV000052347: g.32356514G>A, p.Gly2508Ser.
  assert.deepEqual([mis.an.g!.pos, mis.an.ann!.hgvsP, mis.an.ann!.consequence], [32356514, 'p.(Gly2508Ser)', 'missense'])
  assert.deepEqual([intr.an.g!.pos, intr.an.ann!.consequence, intr.an.ann!.intron], [32319070, 'splice_region', 2])
  // And back: the genomic variant from a VCF gets the same c. description.
  assert.equal(gToHgvsC(tx, genome, { chrom: '13', pos: 32398489, ref: 'A', alt: 'T' }).start.pos, 9976)
})

test('real BRCA1 model (minus strand): c.5266dup and c.5096G>A, left-aligned as a VCF writes them', () => {
  const { tx, genome } = realModel('BRCA1')
  assert.deepEqual([tx.id, tx.strand, tx.exons.length, cdsBounds(tx).aaLength], ['NM_007294.4', -1, 23, 1863])
  const [dup, mis] = evaluateAll(['c.5266dup', 'c.5096G>A'].map(c => analyse(hgvsRec('BRCA1', c), [tx], genome)), {}, BUNDLED_GENES)
  // ClinVar VCV000017677 in VCF form: 17:43057062 T>TG.
  assert.deepEqual(dup.an.g, { chrom: '17', pos: 43057062, ref: 'T', alt: 'TG' })
  assert.equal(dup.an.ann!.hgvsP, 'p.(Gln1756ProfsTer74)')
  // The new stop (codon 1829) is in the last exon, so NMD is not predicted; the truncation cuts into the BRCT domain: PVS1 at Strong.
  assert.equal(dup.an.ann!.nmd, false)
  assert.equal(dup.criteria.find(c => c.code === 'PVS1')!.strength, 'Strong')
  assert.deepEqual([mis.an.g, mis.an.ann!.hgvsP], [{ chrom: '17', pos: 43063930, ref: 'C', alt: 'T' }, 'p.(Arg1699Gln)'])
  // A VCF record and the HGVS description of the same duplication resolve to one c. description.
  const fromVcf = analyse({ ...hgvsRec('BRCA1', ''), hgvsC: undefined, g: { chrom: 'chr17', pos: 43057062, ref: 'T', alt: 'TG' } }, [tx], genome)
  assert.equal(fromVcf.ann!.hgvsC, 'NM_007294.4(BRCA1):c.5266dup')
})

test('several reference regions on one chromosome; variants outside them are trimmed, not corrupted', () => {
  const g: Genome = { '7:101-110': { chrom: '7', start: 101, seq: 'ACGTTTTACG' }, '7:5001-5008': { chrom: 'chr7', start: 5001, seq: 'GGGCCCAA' } }
  assert.equal(refAt(g, 'chr7', 104, 4), 'TTTT')
  assert.equal(refAt(g, '7', 5004, 3), 'CCC')
  assert.equal(refAt(g, '7', 109, 5), null)
  assert.equal(refAt(g, '7', 3000, 1), null)
  // Inside a region the deletion is left-aligned through the T run; between regions it is only trimmed.
  assert.deepEqual(normalize(g, { chrom: '7', pos: 106, ref: 'TT', alt: 'T' }), { chrom: '7', pos: 103, ref: 'GT', alt: 'G', leftAligned: true })
  assert.deepEqual(normalize(g, { chrom: '7', pos: 3000, ref: 'CAA', alt: 'CA' }), { chrom: '7', pos: 3000, ref: 'CA', alt: 'C', leftAligned: false })
  const fa = parseFasta('>chr7:101-110\nACGTTTTACG\n>chr7:5001-5008\nGGGCCCAA\n>chrS9\nACGT\n')
  assert.deepEqual(Object.keys(fa.genome), ['chr7:101-110', 'chr7:5001-5008', 'chrS9'])
  assert.equal(refAt(fa.genome, '7', 5001, 3), 'GGG')
})

test('ClinVar record XML (VCV000038266): submitters, calls, dates and criteria written in comments', () => {
  const rec = parseVcvXml(fixture('clinvar-vcv38266-excerpt.xml'))!
  assert.deepEqual([rec.accession, rec.name, rec.gene, rec.aggregate, rec.aggregateReview], ['VCV000038266', 'NM_000059.4(BRCA2):c.9976A>T (p.Lys3326Ter)', 'BRCA2', 'Benign', 'reviewed by expert panel'])
  assert.equal(rec.submissions.length, 6)
  const prev = rec.submissions.find(s => s.scv === 'SCV000301783')!
  assert.deepEqual([prev.submitter, prev.classification, prev.reviewStatus, prev.method], ['PreventionGenetics, part of Exact Sciences', 'Benign', 'criteria provided, single submitter', 'ACMG Guidelines, 2015'])
  assert.equal(rec.submissions.find(s => /Mayo/.test(s.submitter))!.date, '2017-02-02')
  const asturias = rec.submissions.find(s => /Asturias/.test(s.submitter))!
  assert.deepEqual(asturias.criteria.map(c => `${c.code}_${c.strength}`), ['BS1_Strong', 'BS2_Strong', 'PVS1_Moderate'])
  const { rows, skipped } = clinvarLabRows(rec, 'BRCA2:c.9976A>T', id)
  // "not provided" is not a call.
  assert.deepEqual([rows.length, skipped], [5, 1])
  const d = reconcile(rows)[0]
  assert.equal(d.severity, 'clinically-significant')
  // The lab that wrote BS1+BS2+PVS1_Moderate: -4 -4 +2 = -6 points, likely benign, which is what it called.
  const call = d.calls.find(c => /Asturias/.test(c.lab))!
  assert.deepEqual([call.points, call.implied, call.mismatch], [-6, 'Likely benign', false])
  assert.equal(parseVcvXml('<html>busy</html>'), null)
})

test('criteria in free text: codes said not to apply are left out', () => {
  const codes = (t: string) => criteriaFromComment(t).map(c => `${c.code}_${c.strength}`)
  assert.deepEqual(codes('Classified as likely pathogenic (PVS1_Strong, PM2_Supporting; PP4).'), ['PVS1_Strong', 'PM2_Supporting', 'PP4_Supporting'])
  assert.deepEqual(codes('PM2 is not met. PS3 was not applied because the assay lacks controls. BS1 applies.'), ['BS1_Strong'])
  assert.deepEqual(codes('Functional data do not support PS3; BP5_Very Strong and BA1 met.'), ['BP5_VeryStrong', 'BA1_StandAlone'])
  assert.deepEqual(codes('The variant was observed in TP53BP1 and reported in PMID 12345.'), [])
  assert.deepEqual(clinvarTerms({ transcript: 'NM_000059.4', gene: 'BRCA2', hgvsC: 'c.9976A>T' }), ['"NM_000059" AND "c.9976A>T"', 'BRCA2[gene] AND "c.9976A>T"'])
  assert.ok(titleMatches('NM_000059.4(BRCA2):c.9976A>T (p.Lys3326Ter)', 'c.9976A>T'))
  assert.ok(!titleMatches('NM_000059.4(BRCA2):c.9976A>G (p.Lys3326Glu)', 'c.9976A>T'))
})

test('ClinVar tabular download layout: Name becomes HGVS c. and p., classifications go to the lab comparison', () => {
  const text = fixture('clinvar-result-excerpt.txt')
  assert.equal(detectKind(text, 'clinvar_result.txt'), 'variant-csv')
  const csv = parseCsv(text)
  const map = autoMap(csv.rows[0], VARIANT_FIELDS)
  assert.deepEqual(['hgvsC', 'gene', 'classification', 'reviewStatus', 'classDate'].map(k => csv.rows[0][map[k]]), ['Name', 'Gene(s)', 'Germline classification', 'Germline review status', 'Germline date last evaluated'])
  const r = validateVariantCsv(csv, map, 'clinvar_result.txt', G, id)
  assert.equal(r.issues.filter(i => i.level === 'error').length, 0)
  assert.deepEqual(r.variants.map(v => [v.gene, v.transcript, v.hgvsC, v.hgvsP]), [
    ['BRCA2', 'NM_000059.4', 'c.9976A>T', 'p.Lys3326Ter'], ['BRCA2', 'NM_000059.4', 'c.7522G>A', 'p.Gly2508Ser'],
    ['BRCA1', 'NM_007294.4', 'c.5266dup', 'p.Gln1756fs'], ['CFTR', 'NM_000492.4', 'c.1521_1523del', 'p.Phe508del'],
  ])
  assert.deepEqual(r.labRows.map(x => [x.variant, x.lab, x.cls]), [
    ['BRCA2:c.9976A>T', 'ClinVar (reviewed by expert panel)', 'Benign'], ['BRCA1:c.5266dup', 'ClinVar (reviewed by expert panel)', 'Pathogenic'], ['CFTR:c.1521_1523del', 'ClinVar (practice guideline)', 'Pathogenic'],
  ])
  assert.ok(r.issues.some(i => i.line === 3 && /Conflicting classifications/.test(i.message)))
})

test('one variant, one worklist entry: a VCF record and an HGVS row share an identity', () => {
  const vcf = parseVcf(fixture('snpeff-ann-excerpt.vcf'), 'a.vcf', {}, id).variants
  const csv = parseCsv(fixture('clinvar-result-excerpt.txt'))
  const list = validateVariantCsv(csv, autoMap(csv.rows[0], VARIANT_FIELDS), 'b.txt', G, id).variants
  const have = variantIndex(vcf)
  // BRCA2 c.9976A>T, BRCA1 c.5266dup (written c.5266dupC by SnpEff) and CFTR c.1521_1523del are already there; c.7522G>A is new.
  assert.deepEqual(list.filter(v => !inIndex(have, v)).map(v => v.hgvsC), ['c.7522G>A'])
  assert.deepEqual(identities(vcf[1]), ['17:43057062:T:TG', 'BRCA1:c.5266dup'])
  // Another sample carrying the same variant keeps its own entry; the same sample does not.
  assert.ok(!inIndex(have, { ...list[0], sample: 'CASE-2' }))
  assert.ok(inIndex(have, { ...list[0], sample: 'SAMPLE1' }))
})
