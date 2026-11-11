import { useMemo } from 'react'
import { CLASS_SHORT, critLabel, fmtAf, vusTier, type ClassName, type GeneCfg } from './engine/acmg.ts'
import { toCsv } from './engine/csv.ts'
import { autoMap, LAB_FIELDS, parseCsv, parseVcf, validateLabCsv, validateVariantCsv, VARIANT_FIELDS, type VariantRec } from './engine/io.ts'
import { analyse, evaluateAll, newCuration, type Curation, type Evaluated } from './engine/pipeline.ts'
import { BUNDLED_GENES, BUNDLED_GENOME, BUNDLED_TRANSCRIPTS, LAB_CSV, REAL_VARIANTS_CSV, SAMPLES, sampleVcf } from './engine/sample.ts'
import { CONSEQUENCE_LABEL } from './engine/transcript.ts'
import { defaultGene, editorName, getProjects, saveProject, useProjects, uuid, type Project, type Snapshot } from './store.ts'

export const context = (p: Project) => ({
  txs: [...p.transcripts, ...BUNDLED_TRANSCRIPTS.filter(b => !p.transcripts.some(t => t.id === b.id))],
  genome: { ...BUNDLED_GENOME, ...p.genome },
  genes: { ...BUNDLED_GENES, ...p.genes },
})

export function evaluateProject(p: Project): Evaluated[] {
  const c = context(p)
  return evaluateAll(p.variants.map(v => analyse(v, c.txs, c.genome)), p.curation, c.genes, defaultGene())
}

export function useProject(id: string) {
  const list = useProjects()
  const project = list?.find(p => p.id === id) ?? null
  const evaluated = useMemo(() => project ? evaluateProject(project) : [], [project])
  return { loading: list === null, project, evaluated }
}

export const geneCfg = (p: Project, gene: string): GeneCfg => p.genes[gene] ?? BUNDLED_GENES[gene] ?? defaultGene()

// ---------- Curation actions (every change is audited; class changes go into the reclassification history) ----------

export interface CurateOpts { at?: number; editor?: string }

export function applyCuration(p: Project, vid: string, action: string, detail: string, fn: (d: { rec: VariantRec; cur: Curation }) => void, o: CurateOpts = {}): Project {
  const at = o.at ?? Date.now(), editor = o.editor ?? editorName()
  const before = evaluateProject(p).find(e => e.an.rec.id === vid)!.result
  const variants = p.variants.map(v => v.id === vid ? structuredClone(v) : v)
  const rec = variants.find(v => v.id === vid)!
  const cur = structuredClone(p.curation[vid] ?? newCuration())
  if (!cur.history.length) cur.history.push({ at, editor: 'engine', from: null, to: before.cls, points: before.points, reason: 'Engine classification when curation started' })
  fn({ rec, cur })
  cur.audit.push({ at, editor, action, detail })
  if (cur.status === 'new') cur.status = 'in-review'
  const next = { ...p, variants, curation: { ...p.curation, [vid]: cur } }
  const after = evaluateProject(next).find(e => e.an.rec.id === vid)!.result
  if (after.cls !== before.cls) cur.history.push({ at, editor, from: before.cls, to: after.cls, points: after.points, reason: `${action}: ${detail}` })
  return next
}

export const curate = (p: Project, vid: string, action: string, detail: string, fn: (d: { rec: VariantRec; cur: Curation }) => void) =>
  saveProject(applyCuration(p, vid, action, detail, fn))

// ---------- Snapshots ----------

export function snapshotOf(evaluated: Evaluated[], label: string, at = Date.now()): Snapshot {
  return {
    id: uuid(), label, createdAt: at,
    rows: evaluated.map(e => ({ key: e.an.rec.key, label: e.an.label, cls: e.result.cls, points: e.result.points, criteria: e.result.counted.map(critLabel) })),
  }
}

export interface SnapDiff { reclassified: { key: string; label: string; from: ClassName; to: ClassName; dPoints: number; added: string[]; removed: string[] }[]; added: Snapshot['rows']; removed: Snapshot['rows']; unchanged: number }

export function compareSnapshots(a: Snapshot, b: Snapshot): SnapDiff {
  const A = new Map(a.rows.map(r => [r.key, r])), B = new Map(b.rows.map(r => [r.key, r]))
  const out: SnapDiff = { reclassified: [], added: [], removed: [], unchanged: 0 }
  for (const [k, r] of B) {
    const o = A.get(k)
    if (!o) { out.added.push(r); continue }
    const added = r.criteria.filter(c => !o.criteria.includes(c)), removed = o.criteria.filter(c => !r.criteria.includes(c))
    if (o.cls !== r.cls || added.length || removed.length) out.reclassified.push({ key: k, label: r.label, from: o.cls, to: r.cls, dPoints: r.points - o.points, added, removed })
    else out.unchanged++
  }
  for (const [k, r] of A) if (!B.has(k)) out.removed.push(r)
  return out
}

// ---------- Projects ----------

const blank = (name: string): Project => {
  const now = Date.now()
  return { id: uuid(), name, createdAt: now, updatedAt: now, variants: [], labRows: [], transcripts: [], genome: {}, genes: {}, curation: {}, snapshots: [], imports: [] }
}
export const createBlankProject = (name: string) => saveProject(blank(name))

const DAY = 86400000
const SYN_CURATOR = 'Synthetic curator (sample)'

export async function createSampleProject(sampleId: string): Promise<Project> {
  const existing = getProjects().find(p => p.sampleId === sampleId)
  if (existing) return existing
  const def = SAMPLES.find(s => s.id === sampleId)!
  const t0 = Date.UTC(2026, 6, 14, 9, 30)
  let p: Project = { ...blank(def.name), sampleId, createdAt: t0 }
  if (sampleId === 'synthetic-panel') {
    const r = parseVcf(sampleVcf(), 'SYN-PROBAND-01.synthetic.vcf', BUNDLED_GENOME, uuid)
    p = { ...p, variants: r.variants, imports: [{ id: uuid(), kind: 'vcf', fileName: 'SYN-PROBAND-01.synthetic.vcf', rows: r.variants.length, issues: r.issues.length, at: t0 }] }
    p.snapshots = [snapshotOf(evaluateProject(p), 'Initial import', t0)]
    const find = (c: string) => evaluateProject(p).find(e => e.an.ann?.hgvsC.endsWith(':' + c))!.an.rec.id
    const o = (d: number, h = 10): CurateOpts => ({ at: t0 + d * DAY + h * 3600000, editor: SYN_CURATOR })
    const splice = find('c.306+5G>A')
    p = applyCuration(p, splice, 'Note', 'Proband phenotype consistent with SYNG1 disorder; requested RNA study.', ({ cur }) => { cur.note = 'Near-splice variant; local PWM predicts donor weakening. RNA from blood requested.' }, o(1))
    p = applyCuration(p, splice, 'Functional evidence added', 'Minigene assay (synthetic): exon 3 skipping in 85% of transcripts, OddsPath 6.2', ({ rec }) => {
      rec.evidence.functional.push({ id: uuid(), assay: 'Minigene splicing assay (synthetic result)', result: 'abnormal', oddsPath: 6.2, note: 'Exon 3 skipping in 85% of transcripts vs 0% wild type.', editor: SYN_CURATOR, at: o(21).at! })
    }, o(21))
    p = applyCuration(p, splice, 'Criterion added', 'PP4_Supporting: phenotype highly specific (synthetic case notes)', ({ cur }) => {
      cur.overrides.PP4 = { code: 'PP4', met: true, strength: 'Supporting', evidence: 'Phenotype highly specific for SYNG1 disorder (synthetic case notes).', source: 'curator', editor: SYN_CURATOR, at: o(22).at! }
    }, o(22))
    p = applyCuration(p, splice, 'Status', 'Marked curated', ({ cur }) => { cur.status = 'curated' }, o(22, 15))
    const miss = find('c.208C>T')
    p = applyCuration(p, miss, 'Note', 'Second-tier review: segregation pending from two affected relatives.', ({ cur }) => { cur.note = 'Hot candidate: awaiting segregation (PP1) from relatives.' }, o(23))
    p.snapshots.push(snapshotOf(evaluateProject(p), 'Case review, RNA result in', o(22, 16).at))
  } else {
    const csv = parseCsv(REAL_VARIANTS_CSV)
    const v = validateVariantCsv(csv, autoMap(csv.rows[0], VARIANT_FIELDS), 'public-hgvs.synthetic-evidence.csv', BUNDLED_GENOME, uuid)
    const lc = parseCsv(LAB_CSV)
    const l = validateLabCsv(lc, autoMap(lc.rows[0], LAB_FIELDS), 'lab-calls.synthetic.csv', uuid)
    p = { ...p, variants: v.variants, labRows: l.rows, imports: [
      { id: uuid(), kind: 'variant-csv', fileName: 'public-hgvs.synthetic-evidence.csv', rows: v.variants.length, issues: v.issues.length, at: t0 },
      { id: uuid(), kind: 'lab-csv', fileName: 'lab-calls.synthetic.csv', rows: l.rows.length, issues: l.issues.length, at: t0 },
    ] }
    p.snapshots = [snapshotOf(evaluateProject(p), 'Initial import', t0)]
    const find = (c: string) => p.variants.find(x => x.hgvsC === c)!.id
    const o = (d: number): CurateOpts => ({ at: t0 + d * DAY, editor: SYN_CURATOR })
    p = applyCuration(p, find('c.5266dup'), 'Strength changed', 'PVS1 to Very strong: PTC at codon 1829 lies in exon 20 of 24, well upstream of the last junction (transcript review)', ({ cur }) => {
      cur.overrides.PVS1 = { code: 'PVS1', met: true, strength: 'VeryStrong', evidence: 'PTC at codon 1829 is in exon 20 of 24; predicted NMD (manual transcript review; no model loaded).', source: 'curator', editor: SYN_CURATOR, at: o(2).at! }
    }, o(2))
    p = applyCuration(p, find('c.5266dup'), 'Criterion added', 'PS4_Moderate from case-control data (illustrative)', ({ cur }) => {
      cur.overrides.PS4 = { code: 'PS4', met: true, strength: 'Moderate', evidence: 'Enriched in cases vs controls (illustrative value for the sample).', source: 'curator', editor: SYN_CURATOR, at: o(2).at! }
      cur.status = 'curated'
    }, o(2))
  }
  return saveProject(p)
}

// ---------- Exports ----------

export function worklistCsv(list: Evaluated[], p: Project): string {
  return toCsv([
    ['gene', 'transcript', 'hgvs_c', 'hgvs_p', 'consequence', 'splice_region', 'chrom', 'pos', 'ref', 'alt', 'genotype', 'popmax_af', 'revel', 'spliceai', 'criteria', 'points', 'classification', 'vus_tier', 'status', 'last_editor', 'source'],
    ...list.map(e => {
      const a = e.an.ann, r = e.an.rec, cur = p.curation[r.id], ev = r.evidence
      return [e.an.gene, e.an.tx?.id ?? r.transcript ?? '', a?.hgvsC.replace(/^.*:/, '') ?? r.hgvsC ?? '', a?.hgvsP ?? '', a ? CONSEQUENCE_LABEL[a.consequence] : '', a?.spliceRegion ? 'yes' : '',
        e.an.g?.chrom ?? '', e.an.g?.pos ?? '', e.an.g?.ref ?? '', e.an.g?.alt ?? '', r.gt ?? '', ev.popmaxAf === null ? '' : fmtAf(ev.popmaxAf), ev.revel ?? '', ev.spliceAi ?? '',
        e.result.counted.map(critLabel).join(' '), e.result.points, e.result.cls, vusTier(e.result.points) && e.result.cls === 'VUS' ? vusTier(e.result.points)! : '', cur?.status ?? 'new', cur?.audit.at(-1)?.editor ?? '', `${r.source.file}:${r.source.line}`]
    }),
  ])
}

/** Curated calls in the same layout the lab-table importer reads, so two groups can reconcile against each other. */
export function curatedCsv(list: Evaluated[], p: Project, lab: string): string {
  return toCsv([
    ['variant', 'gene', 'lab', 'classification', 'criteria', 'points', 'date_evaluated', 'curator'],
    ...list.filter(e => p.curation[e.an.rec.id]?.status === 'curated').map(e => {
      const cur = p.curation[e.an.rec.id]
      const last = cur.audit.at(-1)
      return [e.an.ann?.hgvsC ?? e.an.rec.key, e.an.gene, lab, e.result.cls, e.result.counted.map(critLabel).join(', '), e.result.points, last ? new Date(last.at).toISOString().slice(0, 10) : '', last?.editor ?? '']
    }),
  ])
}

export const classTone = (c: ClassName) => CLASS_SHORT[c].toLowerCase()
