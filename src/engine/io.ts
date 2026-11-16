// Importers: VCF 4.x, CSV/TSV variant lists (HGVS), lab classification tables, GTF/GFF3 transcript subsets.
import { norm, parseCsv, parseNum } from './csv.ts'
import { formatC, formatP, parseC, parseP } from './hgvs.ts'
import { parseClass, parseCriteria, type Applied, type ClassName, type Evidence } from './acmg.ts'
import { refAt, sameChrom, type Genome, type GVariant, type Transcript } from './transcript.ts'

export interface Issue { line: number; level: 'error' | 'warning' | 'info'; message: string }

export interface VariantRec {
  id: string
  /** Stable identity used for de-duplication and snapshot comparison. */
  key: string
  gene?: string
  transcript?: string
  g?: GVariant
  hgvsC?: string
  hgvsP?: string
  /** Genotype for the first sample, relative to this ALT after multi-allelic splitting, e.g. "0|1". */
  gt?: string
  ps?: string
  sample?: string
  filter?: string
  source: { file: string; line: number; kind: 'vcf' | 'csv' }
  evidence: Evidence
}

export type FileKind = 'vcf' | 'variant-csv' | 'lab-csv' | 'gtf' | 'gff3' | 'fasta' | 'unknown'
export const FILE_KIND_LABEL: Record<FileKind, string> = { vcf: 'VCF', 'variant-csv': 'Variant list (HGVS)', 'lab-csv': 'Lab classification table', gtf: 'GTF transcript subset', gff3: 'GFF3 transcript subset', fasta: 'FASTA reference region', unknown: 'Unrecognised' }

export function detectKind(text: string, name = ''): FileKind {
  const head = text.slice(0, 4000)
  if (/^##fileformat=VCF/m.test(head) || /^#CHROM\tPOS/m.test(head) || /\.vcf$/i.test(name)) return 'vcf'
  if (/^>/.test(head) || /\.(fa|fasta)$/i.test(name)) return 'fasta'
  if (/^##gff-version\s*3/m.test(head) || /\.gff3?$/i.test(name)) return 'gff3'
  const firstData = head.split(/\r?\n/).find(l => l && !l.startsWith('#')) ?? ''
  const cols = firstData.split('\t')
  if (cols.length === 9 && /^\d+$/.test(cols[3]) && /\b(gene_id|transcript_id)\s+"/.test(cols[8])) return 'gtf'
  if (cols.length === 9 && /^\d+$/.test(cols[3]) && /(ID|Parent)=/.test(cols[8])) return 'gff3'
  const hdr = new Set(parseCsv(firstData || head.split(/\r?\n/)[0]).rows[0]?.map(norm) ?? [])
  if ([...hdr].some(h => LAB_FIELDS[1].aliases.includes(h)) && [...hdr].some(h => LAB_FIELDS[2].aliases.includes(h))) return 'lab-csv'
  if ([...hdr].some(h => /hgvs|c_dot|cdna|variant|^pos|^name$/.test(h))) return 'variant-csv'
  return 'unknown'
}

/** The contig name the loaded reference uses for this chromosome (with or without "chr"), if any region of it is loaded. */
export const findContig = (g: Genome, chrom: string) => Object.values(g).find(r => sameChrom(r.chrom, chrom))?.chrom

// ---------- Normalisation ----------

/** Left-aligns and trims a VCF-style variant (Tan et al. 2015). Without reference sequence only trimming is possible. */
export function normalize(g: Genome, v: GVariant): GVariant & { leftAligned: boolean } {
  let { pos, ref, alt } = v
  // Left-alignment needs reference bases around this position, not just somewhere on the chromosome.
  const contig = refAt(g, v.chrom, v.pos, 1) !== null ? v.chrom : undefined
  let leftAligned = !!contig
  for (let guard = 0; guard < 100000; guard++) {
    let changed = false
    const before = { ref, alt }
    if (ref.length && alt.length && ref[ref.length - 1] === alt[alt.length - 1] && (contig || (ref.length > 1 && alt.length > 1))) {
      ref = ref.slice(0, -1); alt = alt.slice(0, -1); changed = true
    }
    if (!ref.length || !alt.length) {
      const b = contig ? refAt(g, contig, pos - 1, 1) : null
      // Ran off the loaded region: keep the last valid representation.
      if (!b) { ({ ref, alt } = before); leftAligned = false; break }
      ref = b + ref; alt = b + alt; pos--; changed = true
    }
    if (!changed) break
  }
  while (ref.length > 1 && alt.length > 1 && ref[0] === alt[0]) { ref = ref.slice(1); alt = alt.slice(1); pos++ }
  return { chrom: v.chrom, pos, ref, alt, leftAligned }
}

/**
 * Every name a variant can be recognised by: its key, and gene plus canonical c. description
 * ("c.5266dupC" and "c.5266dup" are one variant), so a VCF record and an HGVS row do not become two entries.
 */
export function identities(v: VariantRec): string[] {
  const h = v.hgvsC ? parseC(v.hgvsC) : null
  return [v.key, ...(v.gene && h ? [`${v.gene}:${formatC(h, false)}`] : [])]
}

/** Identity to the samples it has been seen in ("" when the source had no sample, such as a reference list). */
export function variantIndex(list: VariantRec[]): Map<string, Set<string>> {
  const ix = new Map<string, Set<string>>()
  for (const v of list) for (const k of identities(v)) ix.set(k, (ix.get(k) ?? new Set()).add(v.sample ?? ''))
  return ix
}

/** Already present when the same variant exists for the same sample, or either side carries no sample. */
export const inIndex = (ix: Map<string, Set<string>>, v: VariantRec) => identities(v).some(k => { const s = ix.get(k); return !!s && (!v.sample || s.has(v.sample) || s.has('')) })

export const gKey = (v: GVariant) => `${v.chrom.replace(/^chr/i, '')}:${v.pos}:${v.ref}:${v.alt}`

// ---------- VCF ----------

export interface VcfResult { variants: VariantRec[]; labRows: LabRow[]; issues: Issue[]; samples: string[]; records: number; meta: { fileformat?: string; info: string[]; contigs: string[]; annotation?: 'VEP CSQ' | 'SnpEff ANN' | 'ClinVar' } }

const AF_KEYS = ['AF_popmax', 'AF_grpmax', 'popmax_AF', 'gnomAD_AF_popmax', 'MAX_AF', 'gnomADe_AF', 'gnomADg_AF', 'gnomAD_AF', 'AF']

/** SnpEff ANN sub-fields are fixed by the annotation standard (Cingolani et al., VCF annotation format 1.0). */
const ANN_FIELDS = ['Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'Gene_ID', 'Feature_Type', 'Feature_ID', 'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p']

interface TxAnnotation { gene?: string; transcript?: string; hgvsC?: string; hgvsP?: string; af?: number; revel?: number; spliceAi?: number }

const pct = (s: string) => { try { return decodeURIComponent(s) } catch { return s } }
const maxOf = (vals: (string | undefined)[]) => { const n = vals.flatMap(v => (v ?? '').split('&')).filter(v => v.trim() !== '').map(Number).filter(x => !Number.isNaN(x) && x >= 0); return n.length ? Math.max(...n) : undefined }

/**
 * Picks one transcript annotation for an ALT allele from VEP's CSQ (field order from the header) or SnpEff's ANN.
 * Preference: MANE Select, then VEP canonical, then the first coding entry with an HGVS c. description.
 */
export function pickAnnotation(raw: string, fields: string[], alt: string, alleleNum: number, nAlts: number): TxAnnotation | null {
  const ix = (k: string) => fields.indexOf(k)
  const entries = raw.split(',').map(e => e.split('|'))
  const num = ix('ALLELE_NUM')
  const mine = entries.filter(e => num >= 0 && e[num] ? Number(e[num]) === alleleNum : nAlts === 1 || e[0] === alt)
  if (!mine.length) return null
  const get = (e: string[], ...keys: string[]) => { for (const k of keys) { const v = ix(k) >= 0 ? e[ix(k)] : ''; if (v) return pct(v) } return '' }
  const score = (e: string[]) => (get(e, 'MANE_SELECT') ? 8 : 0) + (get(e, 'CANONICAL') === 'YES' ? 4 : 0) + (/c\./.test(get(e, 'HGVSc', 'HGVS.c')) ? 2 : 0) + (get(e, 'BIOTYPE', 'Transcript_BioType') === 'protein_coding' ? 1 : 0)
  const best = mine.reduce((a, b) => score(b) > score(a) ? b : a)
  const out: TxAnnotation = {}
  const gene = get(best, 'SYMBOL', 'Gene_Name'), c = get(best, 'HGVSc', 'HGVS.c'), p = get(best, 'HGVSp', 'HGVS.p')
  if (gene) out.gene = gene
  const feature = c.includes(':') ? c.slice(0, c.indexOf(':')) : get(best, 'Feature', 'Feature_ID')
  if (/^(N[MR]_|X[MR]_|ENST)/.test(feature)) out.transcript = feature
  if (/c\./.test(c) && parseC(c.replace(/^.*?(c\.)/, '$1'))) out.hgvsC = c.replace(/^.*?(c\.)/, '$1')
  const pp = /p\./.test(p) ? parseP(p.replace(/^.*?(p\.)/, '$1')) : null
  // Substitutions are rewritten in three-letter form ("p.Lys3326*" becomes "p.Lys3326Ter"); other types keep the tool's wording.
  if (pp) out.hgvsP = ['missense', 'nonsense', 'synonymous'].includes(pp.kind) && pp.pos ? formatP(pp) : p.replace(/^.*?(p\.)/, '$1')
  const af = maxOf([get(best, 'MAX_AF') || get(best, 'gnomADe_AF', 'gnomADg_AF', 'gnomAD_AF')])
  if (af !== undefined && af <= 1) out.af = af
  const revel = maxOf([get(best, 'REVEL', 'REVEL_score')])
  if (revel !== undefined && revel <= 1) out.revel = revel
  const sai = maxOf(['SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL'].map(k => get(best, k)))
  if (sai !== undefined && sai <= 1) out.spliceAi = sai
  return out
}

/** ClinVar VCF: aggregate CLNSIG, or one row per class from CLNSIGCONF ("Pathogenic(1)|Uncertain_significance(3)") when submitters conflict. */
export function clinvarRows(info: Record<string, string>, variant: string, gene: string | undefined, file: string, line: number, idGen: () => string): LabRow[] {
  const base = { variant, ...(gene ? { gene } : {}), criteria: [], line, file }
  const rows: LabRow[] = []
  if (info.CLNSIGCONF) {
    for (const part of info.CLNSIGCONF.split(/[|,]/)) {
      const m = /^(.*?)\((\d+)\)$/.exec(part.trim())
      const cls = m && parseClass(m[1])
      if (m && cls) rows.push({ ...base, id: idGen(), lab: `ClinVar: ${m[2]} submitter${m[2] === '1' ? '' : 's'} (${cls})`, cls })
    }
    if (rows.length) return rows
  }
  const cls = info.CLNSIG ? parseClass(info.CLNSIG) : null
  if (cls) rows.push({ ...base, id: idGen(), lab: `ClinVar (${(info.CLNREVSTAT ?? 'aggregate').replace(/_/g, ' ').replace(/ ,/g, ',')})`, cls })
  return rows
}

export function parseVcf(text: string, file: string, genome: Genome, idGen: () => string, maxRecords = 20000): VcfResult {
  const issues: Issue[] = []
  const meta: VcfResult['meta'] = { info: [], contigs: [] }
  const infoNumber: Record<string, string> = {}
  let samples: string[] = []
  let header: string[] | null = null
  const out: VariantRec[] = []
  const labRows: LabRow[] = []
  const seen = new Set<string>()
  let records = 0
  let csqFields: string[] | null = null
  let noCarrier = 0
  const unaligned = new Map<string, { n: number; line: number }>()
  const lines = text.split(/\r?\n/)
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i]
    const ln = i + 1
    if (!line.trim()) continue
    if (line.startsWith('##')) {
      const ff = /^##fileformat=(VCFv[\d.]+)/.exec(line)
      if (ff) meta.fileformat = ff[1]
      const inf = /^##INFO=<ID=([^,>]+),Number=([^,>]+)/.exec(line)
      if (inf) { meta.info.push(inf[1]); infoNumber[inf[1]] = inf[2] }
      const csq = /^##INFO=<ID=CSQ,.*Format: ([^"'>]+)/.exec(line)
      if (csq) csqFields = csq[1].trim().split('|').map(s => s.trim())
      const ct = /^##contig=<ID=([^,>]+)/.exec(line)
      if (ct) meta.contigs.push(ct[1])
      continue
    }
    if (line.startsWith('#')) {
      header = line.slice(1).split('\t')
      if (header.length < 8 || header[0] !== 'CHROM') issues.push({ line: ln, level: 'error', message: 'Header line must be tab-separated: #CHROM POS ID REF ALT QUAL FILTER INFO' })
      samples = header.slice(9)
      continue
    }
    if (!header) { issues.push({ line: ln, level: 'error', message: 'Data before the #CHROM header line' }); header = ['CHROM'] }
    const f = line.split('\t')
    if (f.length < 8) { issues.push({ line: ln, level: 'error', message: `Expected at least 8 tab-separated columns, found ${f.length}` }); continue }
    records++
    if (out.length >= maxRecords) { if (out.length === maxRecords) issues.push({ line: ln, level: 'warning', message: `Stopped at ${maxRecords.toLocaleString()} variants; filter the VCF to your gene panel first (see Docs)` }); continue }
    const [chrom, posS, id, refRaw, altRaw, , filter, info] = f
    const pos = Number(posS)
    if (!Number.isInteger(pos) || pos < 1) { issues.push({ line: ln, level: 'error', message: `POS "${posS}" is not a positive integer` }); continue }
    const ref = refRaw.toUpperCase()
    if (!/^[ACGTN]+$/.test(ref)) { issues.push({ line: ln, level: 'error', message: `REF "${refRaw}" must be A/C/G/T/N bases` }); continue }
    const contig = findContig(genome, chrom)
    if (contig) {
      const r = refAt(genome, contig, pos, ref.length)
      if (r !== null && r !== ref) { issues.push({ line: ln, level: 'error', message: `REF ${ref} does not match the reference (${r}) at ${chrom}:${pos}; wrong genome build?` }); continue }
    }
    const infoMap: Record<string, string> = {}
    if (info && info !== '.') for (const kv of info.split(';')) { const k = kv.indexOf('='); infoMap[k < 0 ? kv : kv.slice(0, k)] = k < 0 ? 'true' : kv.slice(k + 1) }
    const fmt = f[8]?.split(':') ?? []
    const s0 = f[9]?.split(':') ?? []
    const gtRaw = fmt.indexOf('GT') >= 0 ? s0[fmt.indexOf('GT')] : undefined
    const ps = fmt.indexOf('PS') >= 0 ? s0[fmt.indexOf('PS')] : undefined
    const alts = altRaw.split(',')
    if (alts.length > 1) issues.push({ line: ln, level: 'info', message: `Multi-allelic site split into ${alts.length} records` })
    alts.forEach((altU, ai) => {
      const alt = altU.toUpperCase()
      if (alt === '.' ) return
      if (alt === '*') { issues.push({ line: ln, level: 'info', message: 'ALT "*" (overlapping deletion) skipped' }); return }
      if (alt === '<NON_REF>' || alt === '<*>') return
      if (/^</.test(alt) || /[[\]]/.test(alt) || /^\.[ACGTN]|[ACGTN]\.$/.test(alt)) { issues.push({ line: ln, level: 'warning', message: `Symbolic/structural ALT ${altU} skipped; this workbench handles sequence variants` }); return }
      if (!/^[ACGTN]+$/.test(alt)) { issues.push({ line: ln, level: 'error', message: `ALT "${altU}" is not a base sequence` }); return }
      const n = normalize(genome, { chrom, pos, ref, alt })
      if (!n.ref && !n.alt) return
      if (n.pos !== pos || n.ref !== ref || n.alt !== alt) issues.push({ line: ln, level: 'info', message: `Normalised ${chrom}:${pos} ${ref}>${alt} to ${n.pos} ${n.ref}>${n.alt}${n.leftAligned ? ' (left-aligned)' : ''}` })
      if (!n.leftAligned && (ref.length !== 1 || alt.length !== 1)) { const u = unaligned.get(chrom) ?? { n: 0, line: ln }; u.n++; unaligned.set(chrom, u) }
      const g = { chrom: n.chrom, pos: n.pos, ref: n.ref, alt: n.alt }
      const key = gKey(g)
      if (seen.has(key)) { issues.push({ line: ln, level: 'warning', message: `Duplicate variant ${key} skipped` }); return }
      seen.add(key)
      const pick = (k: string) => {
        const raw = infoMap[k]
        if (raw === undefined) return undefined
        const parts = raw.split(',')
        return infoNumber[k] === 'A' || parts.length === alts.length ? parts[ai] ?? parts[0] : parts[0]
      }
      // The first sample must carry this allele; 0/0 and no-call rows of a multi-sample or gVCF file are not worklist variants.
      const gt = gtRaw ? splitGt(gtRaw, ai + 1) : undefined
      if (gt && !gt.includes('1')) { noCarrier++; seen.delete(key); return }
      const tx = infoMap.CSQ && csqFields ? pickAnnotation(infoMap.CSQ, csqFields, alt, ai + 1, alts.length) : infoMap.ANN ? pickAnnotation(infoMap.ANN, ANN_FIELDS, alt, ai + 1, alts.length) : null
      if (tx) meta.annotation = infoMap.CSQ ? 'VEP CSQ' : 'SnpEff ANN'
      const afKey = AF_KEYS.find(k => pick(k) !== undefined && k !== 'AF')
      const af = afKey ? parseNum(pick(afKey)!) : tx?.af ?? NaN
      const revel = parseNum(pick('REVEL') ?? pick('REVEL_score') ?? String(tx?.revel ?? ''))
      const sai = infoMap.SpliceAI ? spliceAiMax(infoMap.SpliceAI, alt) : tx?.spliceAi ?? NaN
      const rec: VariantRec = {
        id: idGen(), key, g, source: { file, line: ln, kind: 'vcf' },
        evidence: { popmaxAf: Number.isNaN(af) ? null : af, revel: Number.isNaN(revel) ? null : revel, spliceAi: Number.isNaN(sai) ? null : sai, functional: [] },
      }
      // ClinVar writes GENEINFO=SYMBOL:id|SYMBOL2:id2; the first gene is used.
      const gene = infoMap.GENE ?? infoMap.Gene ?? infoMap.SYMBOL ?? tx?.gene ?? infoMap.GENEINFO?.split(/[:|]/)[0]
      if (gene && gene !== '.') rec.gene = gene
      if (tx?.transcript) rec.transcript = tx.transcript
      if (tx?.hgvsC) rec.hgvsC = tx.hgvsC
      if (tx?.hgvsP) rec.hgvsP = tx.hgvsP
      if (infoMap.CLNSIG || infoMap.CLNSIGCONF) { meta.annotation = 'ClinVar'; labRows.push(...clinvarRows(infoMap, key, rec.gene, file, ln, idGen)) }
      if (filter && filter !== 'PASS' && filter !== '.') { rec.filter = filter; issues.push({ line: ln, level: 'warning', message: `FILTER is ${filter}; kept but flagged` }) }
      if (gt) rec.gt = gt
      if (ps && ps !== '.') rec.ps = ps
      if (samples[0]) rec.sample = samples[0]
      out.push(rec)
    })
  }
  if (!meta.fileformat) issues.unshift({ line: 1, level: 'warning', message: 'Missing ##fileformat line; assuming VCFv4.2' })
  if (!header) issues.unshift({ line: 1, level: 'error', message: 'No #CHROM header line found' })
  for (const [chrom, u] of unaligned) issues.push({ line: u.line, level: 'warning', message: `No reference loaded around ${u.n} indel or multi-base allele${u.n === 1 ? '' : 's'} on ${chrom}; trimmed but not left-aligned (first at this line)` })
  if (samples.length > 1) issues.push({ line: 1, level: 'info', message: `${samples.length} samples; genotypes read from the first (${samples[0]})` })
  if (noCarrier) issues.push({ line: 1, level: 'info', message: `${noCarrier} allele${noCarrier === 1 ? '' : 's'} skipped: ${samples[0] ?? 'the first sample'} is homozygous reference or not called there` })
  if (meta.annotation) issues.push({ line: 1, level: 'info', message: meta.annotation === 'ClinVar' ? `ClinVar annotations found: gene from GENEINFO, ${labRows.length} classification${labRows.length === 1 ? '' : 's'} added to the lab comparison` : `${meta.annotation} annotations found: gene, transcript, HGVS and scores read from the preferred transcript (MANE Select, then canonical)` })
  return { variants: out, labRows, issues, samples, records, meta }
}

/** Re-expresses a genotype relative to one ALT allele: "1|2" for allele 2 becomes "0|1". */
export function splitGt(gt: string, allele: number): string {
  return gt.replace(/\d+|\./g, a => a === '.' ? '.' : Number(a) === allele ? '1' : '0')
}

/** SpliceAI INFO: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL, comma-separated per allele/gene. */
export function spliceAiMax(raw: string, alt: string): number {
  let best = NaN
  for (const e of raw.split(',')) {
    const p = e.split('|')
    if (p.length < 6 || (p[0] !== alt && p[0] !== '.')) continue
    const m = Math.max(...p.slice(2, 6).map(Number).filter(x => !Number.isNaN(x)))
    if (!Number.isNaN(m) && (Number.isNaN(best) || m > best)) best = m
  }
  return best
}

// ---------- Column-mapped CSV ----------

export interface FieldSpec { key: string; label: string; required?: boolean; aliases: string[] }
export type Mapping = Record<string, number>

export const VARIANT_FIELDS: FieldSpec[] = [
  { key: 'hgvsC', label: 'HGVS c. (required unless genomic columns are mapped)', aliases: ['hgvs_c', 'hgvsc', 'c_dot', 'cdna', 'cdna_change', 'hgvs', 'variant', 'nucleotide_change', 'c_hgvs', 'hgvs_cdna', 'name'] },
  { key: 'hgvsP', label: 'HGVS p.', aliases: ['hgvs_p', 'hgvsp', 'p_dot', 'protein', 'protein_change', 'amino_acid_change', 'p_hgvs', 'hgvs_protein'] },
  { key: 'gene', label: 'Gene symbol', aliases: ['gene', 'gene_symbol', 'symbol', 'hgnc_symbol', 'gene(s)', 'genes'] },
  { key: 'transcript', label: 'Transcript', aliases: ['transcript', 'transcript_id', 'refseq', 'nm', 'feature'] },
  { key: 'chrom', label: 'Chromosome', aliases: ['chrom', 'chr', 'chromosome', 'contig', '#chrom'] },
  { key: 'pos', label: 'Position', aliases: ['pos', 'position', 'start', 'genomic_position'] },
  { key: 'ref', label: 'REF', aliases: ['ref', 'reference', 'ref_allele'] },
  { key: 'alt', label: 'ALT', aliases: ['alt', 'alternate', 'alt_allele'] },
  { key: 'af', label: 'Popmax allele frequency', aliases: ['popmax_af', 'af_popmax', 'gnomad_popmax', 'grpmax_af', 'gnomad_af', 'gnomad', 'af', 'allele_frequency', 'max_af'] },
  { key: 'revel', label: 'REVEL score', aliases: ['revel', 'revel_score'] },
  { key: 'spliceAi', label: 'SpliceAI max delta', aliases: ['spliceai', 'spliceai_max', 'spliceai_ds', 'splice_ai', 'spliceai_delta'] },
  { key: 'gt', label: 'Genotype / zygosity', aliases: ['gt', 'genotype', 'zygosity', 'zyg'] },
  { key: 'ps', label: 'Phase set', aliases: ['ps', 'phase_set', 'phase', 'haplotype'] },
  { key: 'sample', label: 'Sample / case ID', aliases: ['sample', 'sample_id', 'case', 'case_id', 'proband'] },
  // Present in ClinVar's tabular download; the calls go to the lab comparison.
  { key: 'classification', label: 'Reported classification (e.g. ClinVar)', aliases: ['germline_classification', 'clinical_significance_(last_reviewed)', 'clinical_significance', 'clinvar_classification'] },
  { key: 'reviewStatus', label: 'Review status', aliases: ['germline_review_status', 'review_status'] },
  { key: 'classDate', label: 'Date last evaluated', aliases: ['germline_date_last_evaluated', 'date_last_evaluated', 'last_evaluated'] },
]

export const LAB_FIELDS: FieldSpec[] = [
  { key: 'variant', label: 'Variant (HGVS c. or chrom:pos:ref:alt)', required: true, aliases: ['variant', 'hgvs', 'hgvs_c', 'hgvsc', 'c_dot', 'nucleotide_change', 'variant_name'] },
  { key: 'lab', label: 'Laboratory', required: true, aliases: ['lab', 'laboratory', 'submitter', 'lab_name', 'source'] },
  { key: 'classification', label: 'Classification', required: true, aliases: ['classification', 'class', 'interpretation', 'clinical_significance', 'call', 'significance'] },
  { key: 'criteria', label: 'ACMG criteria', aliases: ['criteria', 'acmg_criteria', 'evidence', 'acmg', 'codes'] },
  { key: 'gene', label: 'Gene', aliases: ['gene', 'gene_symbol', 'symbol'] },
  { key: 'date', label: 'Date evaluated', aliases: ['date', 'date_evaluated', 'evaluated', 'last_evaluated', 'report_date'] },
]

export function autoMap(headers: string[], fields: FieldSpec[]): Mapping {
  const hs = headers.map(norm)
  const map: Mapping = {}
  const used = new Set<number>()
  for (const f of fields) for (const a of f.aliases) {
    const i = hs.indexOf(a)
    if (i >= 0 && !used.has(i)) { map[f.key] = i; used.add(i); break }
  }
  return map
}

/** Normalises zygosity words to a genotype. */
export function toGt(raw: string): string | undefined {
  const s = raw.trim().toLowerCase()
  if (!s) return undefined
  if (/^[0-9.]+[|/][0-9.]+$/.test(s)) return s
  if (/^het/.test(s)) return '0/1'
  if (/^hom/.test(s)) return '1/1'
  if (/^hemi/.test(s)) return '1'
  return undefined
}

export interface CsvVariantResult { variants: VariantRec[]; labRows: LabRow[]; issues: Issue[]; total: number }

export function validateVariantCsv(csv: { rows: string[][]; lines: number[] }, map: Mapping, file: string, genome: Genome, idGen: () => string): CsvVariantResult {
  const issues: Issue[] = []
  const total = Math.max(0, csv.rows.length - 1)
  const genomic = ['chrom', 'pos', 'ref', 'alt'].every(k => k in map)
  if (!('hgvsC' in map) && !genomic) return { variants: [], labRows: [], issues: [{ line: 1, level: 'error', message: 'Map an HGVS c. column, or all four of chromosome, position, REF and ALT' }], total }
  if ('hgvsC' in map && !('gene' in map) && !('transcript' in map)) issues.push({ line: 1, level: 'warning', message: 'No gene or transcript column; HGVS values must carry a prefix such as NM_000059.4(BRCA2):c.68-7T>A' })
  const out: VariantRec[] = []
  const labRows: LabRow[] = []
  const seen = new Set<string>()
  const get = (rec: string[], k: string) => map[k] === undefined ? '' : (rec[map[k]] ?? '').trim()
  for (let r = 1; r < csv.rows.length; r++) {
    const rec = csv.rows[r], line = csv.lines[r]
    const num = (k: string, lo: number, hi: number) => {
      const s = get(rec, k)
      if (!s || s === '.' || /^(na|n\/a|-)$/i.test(s)) return null
      const n = parseNum(s)
      if (Number.isNaN(n) || n < lo || n > hi) { issues.push({ line, level: 'warning', message: `${k} "${s}" ignored (expected ${lo} to ${hi})` }); return null }
      return n
    }
    const ev: Evidence = { popmaxAf: num('af', 0, 1), revel: num('revel', 0, 1), spliceAi: num('spliceAi', 0, 1), functional: [] }
    const v: VariantRec = { id: idGen(), key: '', source: { file, line, kind: 'csv' }, evidence: ev }
    // ClinVar lists overlapping genes as "A|B"; the HGVS prefix, when present, names the right one.
    const gene = get(rec, 'gene').split('|')[0], tx = get(rec, 'transcript')
    if (gene) v.gene = gene
    if (tx) v.transcript = tx
    const c = get(rec, 'hgvsC')
    if (c) {
      const h = parseC(c.includes(':') || !(tx || gene) ? c : `${tx || gene}:${c.replace(/^.*?c\./, 'c.')}`)
      if (!h) { issues.push({ line, level: 'error', message: `HGVS "${c}" is not a supported c. description (substitution, del, dup, ins, delins)` }); continue }
      v.hgvsC = c.replace(/^.*?(c\.)/, '$1').replace(/\s+\(?p\..*$/, '')
      if (h.gene) v.gene = h.gene
      v.transcript ??= h.transcript
      if (!v.gene && !v.transcript) { issues.push({ line, level: 'error', message: 'No gene or transcript for this variant' }); continue }
    }
    // A protein change may ride along in the name: "... (p.Gln1756fs)".
    const pCol = get(rec, 'hgvsP')
    const p = /p\./.test(pCol) ? pCol : /\s\((p\.[^)]+)\)\s*$/.exec(c)?.[1] ?? pCol
    if (p && p !== '.' && p !== 'p.?') {
      if (!parseP(p.replace(/^.*?(p\.)/, '$1'))) issues.push({ line, level: 'warning', message: `HGVS p. "${p}" not recognised; ignored` })
      else v.hgvsP = p.replace(/^.*?(p\.)/, '$1')
    }
    if (genomic && get(rec, 'pos')) {
      const pos = Number(get(rec, 'pos')), ref = get(rec, 'ref').toUpperCase(), alt = get(rec, 'alt').toUpperCase()
      if (!Number.isInteger(pos) || !/^[ACGTN]*$/.test(ref) || !/^[ACGTN]*$/.test(alt)) issues.push({ line, level: 'warning', message: 'Genomic columns not valid; using HGVS only' })
      else {
        const chrom = get(rec, 'chrom')
        const contig = findContig(genome, chrom)
        const r = contig ? refAt(genome, contig, pos, ref.length) : null
        if (r !== null && ref && r !== ref) { issues.push({ line, level: 'error', message: `REF ${ref} does not match the reference (${r}) at ${chrom}:${pos}` }); continue }
        const n = normalize(genome, { chrom, pos, ref, alt })
        v.g = { chrom: n.chrom, pos: n.pos, ref: n.ref, alt: n.alt }
      }
    }
    if (!v.g && !v.hgvsC) { issues.push({ line, level: 'error', message: 'Neither HGVS c. nor genomic coordinates given' }); continue }
    const gt = toGt(get(rec, 'gt'))
    if (gt) v.gt = gt
    else if (get(rec, 'gt')) issues.push({ line, level: 'warning', message: `Genotype "${get(rec, 'gt')}" not recognised (use 0/1, 0|1, het, hom)` })
    if (get(rec, 'ps')) v.ps = get(rec, 'ps')
    if (get(rec, 'sample')) v.sample = get(rec, 'sample')
    v.key = v.g ? gKey(v.g) : `${(v.transcript || v.gene)!.replace(/\.\d+$/, '')}:${v.hgvsC}`
    if (seen.has(v.key + (v.sample ?? ''))) { issues.push({ line, level: 'warning', message: `Duplicate variant ${v.hgvsC ?? v.key} skipped` }); continue }
    seen.add(v.key + (v.sample ?? ''))
    out.push(v)
    const clsRaw = get(rec, 'classification').replace(/\(last reviewed.*$/i, '')
    if (clsRaw) {
      const cls = parseClass(clsRaw)
      const review = get(rec, 'reviewStatus'), date = get(rec, 'classDate')
      if (cls) labRows.push({ id: idGen(), variant: v.hgvsC ? labVariantKey(v.hgvsC, v.gene) : v.key, ...(v.gene ? { gene: v.gene } : {}), lab: review ? `ClinVar (${review})` : 'Reported classification', cls, criteria: [], ...(date ? { date } : {}), line, file })
      else issues.push({ line, level: 'info', message: `Classification "${clsRaw}" is not a single germline call; variant imported without it` })
    }
  }
  return { variants: out, labRows, issues, total }
}

// ---------- Lab classification tables ----------

export interface LabRow { id: string; variant: string; gene?: string; lab: string; cls: ClassName; criteria: Applied[]; date?: string; line: number; file: string }

export function validateLabCsv(csv: { rows: string[][]; lines: number[] }, map: Mapping, file: string, idGen: () => string): { rows: LabRow[]; issues: Issue[]; total: number } {
  const issues: Issue[] = []
  const total = Math.max(0, csv.rows.length - 1)
  for (const f of LAB_FIELDS) if (f.required && !(f.key in map)) issues.push({ line: 1, level: 'error', message: `Required column "${f.label}" is not mapped` })
  if (issues.length) return { rows: [], issues, total }
  const out: LabRow[] = []
  const get = (rec: string[], k: string) => map[k] === undefined ? '' : (rec[map[k]] ?? '').trim()
  for (let r = 1; r < csv.rows.length; r++) {
    const rec = csv.rows[r], line = csv.lines[r]
    const variant = get(rec, 'variant'), lab = get(rec, 'lab'), clsRaw = get(rec, 'classification')
    if (!variant) { issues.push({ line, level: 'error', message: 'Variant is empty' }); continue }
    if (!lab) { issues.push({ line, level: 'error', message: 'Laboratory is empty' }); continue }
    const cls = parseClass(clsRaw)
    if (!cls) { issues.push({ line, level: 'error', message: `Classification "${clsRaw}" not recognised (use P, LP, VUS, LB, B or the full names)` }); continue }
    const { applied, unknown } = parseCriteria(get(rec, 'criteria'))
    if (unknown.length) issues.push({ line, level: 'warning', message: `Unrecognised criteria ignored: ${unknown.join(', ')}` })
    const gene = get(rec, 'gene')
    out.push({ id: idGen(), variant: labVariantKey(variant, gene), ...(gene ? { gene } : {}), lab, cls, criteria: applied, ...(get(rec, 'date') ? { date: get(rec, 'date') } : {}), line, file })
  }
  return { rows: out, issues, total }
}

/** Groups lab rows by a normalised variant description: gene/transcript prefix kept, versions and p. dropped. */
export function labVariantKey(variant: string, gene?: string): string {
  const v = variant.trim().replace(/\s+\(?p\..*$/, '').replace(/\s+/g, '')
  const h = parseC(v.includes(':') || !gene ? v : `${gene}:${v}`)
  if (!h) return v
  const prefix = (h.gene ?? gene ?? h.transcript ?? '').replace(/\.\d+$/, '')
  return `${prefix}:${v.replace(/^.*?(c\.)/, '$1')}`
}

// ---------- GTF / GFF3 ----------

export function parseGtf(text: string, source: string): { transcripts: Transcript[]; issues: Issue[] } {
  const issues: Issue[] = []
  const tx = new Map<string, { gene: string; chrom: string; strand: 1 | -1; exons: [number, number][]; cds: [number, number][]; stop: [number, number][] }>()
  const geneNames = new Map<string, string>()
  const lines = text.split(/\r?\n/)
  for (let i = 0; i < lines.length; i++) {
    const l = lines[i]
    if (!l || l.startsWith('#')) continue
    const f = l.split('\t')
    if (f.length !== 9) { issues.push({ line: i + 1, level: 'error', message: `Expected 9 tab-separated columns, found ${f.length}` }); continue }
    const [chrom, , type, s, e, , strand, , attr] = f
    const start = Number(s), end = Number(e)
    if (!Number.isInteger(start) || !Number.isInteger(end) || start > end) { issues.push({ line: i + 1, level: 'error', message: 'Start/end are not valid coordinates' }); continue }
    const a: Record<string, string> = {}
    for (const m of attr.matchAll(/(\w+)\s+"([^"]*)"/g)) a[m[1]] = m[2]
    for (const kv of attr.split(';')) { const k = kv.indexOf('='); if (k > 0) a[kv.slice(0, k).trim()] = decodeURIComponent(kv.slice(k + 1).trim()) }
    if (type === 'gene' && a.ID) geneNames.set(a.ID, a.Name ?? a.gene_name ?? a.ID)
    if (type === 'mRNA' || type === 'transcript') {
      const id = (a.transcript_id ?? a.ID ?? '').replace(/^transcript:/, '')
      if (id && a.Parent) geneNames.set(id, geneNames.get(a.Parent) ?? a.Parent)
      continue
    }
    if (type !== 'exon' && type !== 'CDS' && type !== 'stop_codon') continue
    const id = (a.transcript_id ?? a.Parent ?? '').replace(/^transcript:/, '')
    if (!id) { issues.push({ line: i + 1, level: 'warning', message: `${type} without transcript_id/Parent skipped` }); continue }
    if (strand !== '+' && strand !== '-') { issues.push({ line: i + 1, level: 'error', message: `Strand "${strand}" must be + or -` }); continue }
    const t = tx.get(id) ?? { gene: a.gene_name ?? a.gene ?? '', chrom, strand: strand === '+' ? 1 : -1, exons: [], cds: [], stop: [] }
    if (!t.gene && (a.gene_name || a.gene)) t.gene = a.gene_name ?? a.gene
    ;(type === 'exon' ? t.exons : type === 'CDS' ? t.cds : t.stop).push([start, end])
    tx.set(id, t)
  }
  const transcripts: Transcript[] = []
  for (const [id, t] of tx) {
    if (!t.exons.length || !t.cds.length) { issues.push({ line: 0, level: 'warning', message: `${id}: needs both exon and CDS features; skipped` }); continue }
    const all = [...t.cds, ...t.stop]
    const cdsStart = Math.min(...all.map(x => x[0])), cdsEnd = Math.max(...all.map(x => x[1]))
    const exons = t.exons.sort((a, b) => a[0] - b[0])
    const len = [...t.cds, ...t.stop].reduce((s, x) => s + x[1] - x[0] + 1, 0)
    if (len % 3 !== 0) issues.push({ line: 0, level: 'warning', message: `${id}: CDS length ${len} is not a multiple of 3 (missing stop_codon rows?)` })
    transcripts.push({ id, gene: t.gene || geneNames.get(id) || id, chrom: t.chrom, strand: t.strand, exons, cdsStart, cdsEnd, source })
  }
  if (!transcripts.length) issues.push({ line: 0, level: 'error', message: 'No complete transcripts found (need exon and CDS rows with transcript_id or Parent)' })
  return { transcripts, issues }
}

// ---------- FASTA ----------

/** Reads FASTA records. Headers like ">chr13:32315000-32400000" (samtools faidx output) set the start coordinate; plain ">chr13" starts at 1. */
export function parseFasta(text: string): { genome: Genome; issues: Issue[] } {
  const genome: Genome = {}
  const issues: Issue[] = []
  let cur: { chrom: string; start: number; parts: string[]; line: number } | null = null
  const flush = () => {
    if (!cur) return
    const seq = cur.parts.join('').toUpperCase()
    if (/[^ACGTN]/.test(seq)) issues.push({ line: cur.line, level: 'warning', message: `${cur.chrom}: IUPAC codes other than N were read as N` })
    // Keyed by region so two genes on one chromosome can both be loaded.
    genome[cur.start === 1 ? cur.chrom : `${cur.chrom}:${cur.start}-${cur.start + seq.length - 1}`] = { chrom: cur.chrom, start: cur.start, seq: seq.replace(/[^ACGTN]/g, 'N') }
  }
  text.split(/\r?\n/).forEach((l, i) => {
    if (l.startsWith('>')) {
      flush()
      const h = l.slice(1).trim().split(/\s+/)[0]
      const m = /^(.+):(\d+)-(\d+)$/.exec(h.replace(/,/g, ''))
      cur = m ? { chrom: m[1], start: Number(m[2]), parts: [], line: i + 1 } : { chrom: h, start: 1, parts: [], line: i + 1 }
    } else if (l.trim()) {
      if (!cur) issues.push({ line: i + 1, level: 'error', message: 'Sequence before the first > header' })
      else cur.parts.push(l.trim())
    }
  })
  flush()
  if (!Object.keys(genome).length) issues.push({ line: 1, level: 'error', message: 'No FASTA records found' })
  return { genome, issues }
}

export { parseCsv }
