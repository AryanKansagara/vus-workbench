// Ties the engine together: resolve a transcript, annotate, score splicing, work out phase, and evaluate ACMG criteria.
import { appliedOf, combine, DEFAULT_GENE, effective, suggest, type ClassName, type Combined, type CriterionState, type GeneCfg, type Phase, type PhasePartner, type Suggestion } from './acmg.ts'
import { parseC } from './hgvs.ts'
import type { VariantRec } from './io.ts'
import { spliceReport, type SpliceReport } from './splice.ts'
import { annotate, annotateHgvs, gToHgvsC, hgvsCToG, refAt, type Annotation, type Genome, type GVariant, type Transcript } from './transcript.ts'
import { findContig } from './io.ts'

export interface AuditEntry { at: number; editor: string; action: string; detail: string }
export interface ClassChange { at: number; editor: string; from: ClassName | null; to: ClassName; points: number; reason: string }
export interface Curation {
  overrides: Record<string, CriterionState>
  status: 'new' | 'in-review' | 'curated'
  note: string
  audit: AuditEntry[]
  history: ClassChange[]
}
export const newCuration = (): Curation => ({ overrides: {}, status: 'new', note: '', audit: [], history: [] })

export interface Analysed {
  rec: VariantRec
  tx: Transcript | null
  ann: Annotation | null
  g: GVariant | null
  splice: SpliceReport | null
  gene: string
  label: string
  problems: string[]
}

const stripVer = (s: string) => s.replace(/\.\d+$/, '')

export function resolveTranscript(rec: VariantRec, txs: Transcript[]): Transcript | null {
  if (rec.transcript) {
    const t = txs.find(t => stripVer(t.id) === stripVer(rec.transcript!))
    if (t) return t
  }
  if (rec.gene) {
    const t = txs.find(t => t.gene === rec.gene)
    if (t) return t
  }
  if (rec.g) {
    const pos = rec.g.pos
    return txs.find(t => t.chrom.replace(/^chr/i, '') === rec.g!.chrom.replace(/^chr/i, '') && pos >= t.exons[0][0] - 5000 && pos <= t.exons[t.exons.length - 1][1] + 5000) ?? null
  }
  return null
}

export function analyse(rec: VariantRec, txs: Transcript[], genome: Genome): Analysed {
  const problems: string[] = []
  const tx = resolveTranscript(rec, txs)
  let g = rec.g ?? null
  // Work in the contig name the reference uses (chr prefix or not).
  if (g && tx) g = { ...g, chrom: tx.chrom }
  let ann: Annotation | null = null
  if (tx) {
    const hasSeq = !!findContig(genome, tx.chrom)
    if (!g && rec.hgvsC) {
      const h = parseC(`${tx.id}:${rec.hgvsC}`)
      if (h && hasSeq) {
        g = hgvsCToG(tx, genome, h)
        if (!g) problems.push(`The reference base(s) do not match ${rec.hgvsC} on ${tx.id}; check the transcript version.`)
      }
      if (h) ann = annotate(tx, genome, { ...h, gene: tx.gene }, g, rec.hgvsP)
    } else if (g) {
      if (hasSeq && refAt(genome, tx.chrom, g.pos, g.ref.length) !== g.ref) problems.push('REF does not match the loaded reference.')
      const h = gToHgvsC(tx, genome, g)
      ann = annotate(tx, genome, h, hasSeq ? g : null, rec.hgvsP)
    }
  } else if (rec.hgvsC) {
    const h = parseC(`${rec.transcript ?? rec.gene}:${rec.hgvsC}`)
    if (h) ann = annotateHgvs({ ...h, ...(rec.gene ? { gene: rec.gene } : {}) }, rec.hgvsP)
  }
  if (!ann && !problems.length) problems.push(rec.g ? 'Outside every loaded transcript model. Import a GTF subset for this gene to annotate it.' : 'Could not annotate this variant.')
  const splice = tx && g && ann && findContig(genome, tx.chrom) ? spliceReport(tx, genome, g) : null
  const gene = tx?.gene ?? rec.gene ?? '—'
  const label = ann ? `${gene} ${ann.hgvsC.replace(/^.*:/, '')}${ann.hgvsP ? ' ' + ann.hgvsP : ''}` : rec.key
  return { rec, tx, ann, g, splice, gene, label, problems }
}

// ---------- Phase ----------

/** Phase of two variants from their genotypes: phased "0|1" vs "1|0" is trans; same haplotype is cis. */
export function phaseOf(a: VariantRec, b: VariantRec): Phase {
  const ga = a.gt ?? '', gb = b.gt ?? ''
  if (/^1[|/]1$/.test(ga) || /^1[|/]1$/.test(gb)) return 'homozygous'
  if (!ga.includes('|') || !gb.includes('|')) return 'unknown'
  if ((a.ps ?? '') !== (b.ps ?? '')) return 'unknown'
  const ha = ga.split('|').indexOf('1'), hb = gb.split('|').indexOf('1')
  if (ha < 0 || hb < 0) return 'unknown'
  return ha === hb ? 'cis' : 'trans'
}

// ---------- Evaluation ----------

export interface Evaluated { an: Analysed; suggestions: Suggestion[]; criteria: CriterionState[]; result: Combined; partners: PhasePartner[] }

export function evaluate(an: Analysed, cur: Curation | undefined, gene: GeneCfg, partners: PhasePartner[] = []): Evaluated {
  let ann = an.ann
  // Without a transcript model, size a truncation from the gene's protein length (gene setting).
  if (ann && ann.removedFraction == null && ann.ptcAa && gene.proteinLength) ann = { ...ann, removedFraction: Math.max(0, (gene.proteinLength - ann.ptcAa + 1) / gene.proteinLength) }
  const suggestions = ann ? suggest(ann, an.rec.evidence, gene, an.splice, partners) : []
  const criteria = effective(suggestions, cur?.overrides ?? {})
  return { an, suggestions, criteria, result: combine(appliedOf(criteria)), partners }
}

/** Evaluates a whole worklist. Phase partners are classified in a first pass without phase evidence, so PM3/BP2 never feed back on themselves. */
export function evaluateAll(list: Analysed[], cur: Record<string, Curation>, genes: Record<string, GeneCfg>, fallback: GeneCfg = DEFAULT_GENE): Evaluated[] {
  const cfg = (g: string) => genes[g] ?? fallback
  const first = list.map(an => evaluate(an, cur[an.rec.id], cfg(an.gene)))
  const byGene = new Map<string, Evaluated[]>()
  for (const e of first) {
    const k = `${e.an.gene}|${e.an.rec.sample ?? ''}`
    byGene.set(k, [...(byGene.get(k) ?? []), e])
  }
  return first.map(e => {
    const peers = (byGene.get(`${e.an.gene}|${e.an.rec.sample ?? ''}`) ?? []).filter(o => o !== e && (e.an.rec.gt || o.an.rec.gt))
    if (!peers.length) return e
    const partners = peers.map(o => ({ label: o.an.label, cls: o.result.cls, phase: phaseOf(e.an.rec, o.an.rec) }))
    return evaluate(e.an, cur[e.an.rec.id], cfg(e.an.gene), partners)
  })
}
