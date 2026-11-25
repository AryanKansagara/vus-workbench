// Ties the engine together: resolve a transcript, annotate, score splicing, work out phase, and evaluate ACMG criteria.
import { appliedOf, combine, DEFAULT_GENE, effective, suggest, type ClassName, type Combined, type CriterionState, type GeneCfg, type Phase, type PhasePartner, type Suggestion } from './acmg.ts'
import { parseC } from './hgvs.ts'
import { normalize, type VariantRec } from './io.ts'
import { spliceReport, type SpliceReport } from './splice.ts'
import { annotate, annotateHgvs, endsInStop, gToHgvsC, hasReference, hgvsCToG, orderProblem, positionProblem, refAt, withStop, type Annotation, type Genome, type GVariant, type Transcript } from './transcript.ts'

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

const isAccession = (s: string) => /^([NX][MR]_|ENST|LRG_|SYN_)/i.test(s)

/** By transcript accession in any version. A MANE Select model is named by its RefSeq accession and carries the identical Ensembl one in its source. */
const byAccession = (id: string, txs: Transcript[]) => txs.find(t => stripVer(t.id) === stripVer(id)) ?? txs.find(t => t.source.split(/\s+/).some(w => /^ENST/.test(w) && stripVer(w) === stripVer(id)))

export function resolveTranscript(rec: VariantRec, txs: Transcript[]): Transcript | null {
  if (rec.transcript) {
    const t = byAccession(rec.transcript, txs)
    if (t) return t
    // c. numbers belong to their transcript: without genomic coordinates another model of the gene cannot stand in for it.
    if (!rec.g && isAccession(rec.transcript)) return null
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
  const found = resolveTranscript(rec, txs)
  const tx = found && withStop(found, genome)
  let g = rec.g ?? null
  // Work in the contig name the reference uses (chr prefix or not).
  if (g && tx) g = { ...g, chrom: tx.chrom }
  let ann: Annotation | null = null
  if (tx) {
    const hasSeq = hasReference(genome, tx)
    if (endsInStop(tx, genome) === false) problems.push(`The coding sequence of ${tx.id} does not end in a stop codon in the loaded reference, so protein effects are not predicted. Check the model and the genome build.`)
    if (!g && rec.hgvsC) {
      const h = parseC(`${tx.id}:${rec.hgvsC}`)
      const bad = h && positionProblem(tx, h, hasSeq)
      if (bad) problems.push(`${bad} Not annotated; check the description and the transcript version.`)
      else if (h && hasSeq) {
        g = hgvsCToG(tx, genome, h)
        // Same left-aligned form a VCF caller writes, so one variant has one genomic identity.
        if (g) { const n = normalize(genome, g); g = { chrom: n.chrom, pos: n.pos, ref: n.ref, alt: n.alt } }
        if (!g) problems.push(`The reference base(s) do not match ${rec.hgvsC} on ${tx.id}; check the transcript version.`)
      }
      if (h && !bad) ann = annotate(tx, genome, { ...h, gene: tx.gene }, g, rec.hgvsP)
    } else if (g) {
      if (hasSeq && refAt(genome, tx.chrom, g.pos, g.ref.length) !== g.ref) problems.push('REF does not match the loaded reference.')
      const h = gToHgvsC(tx, genome, g)
      ann = annotate(tx, genome, h, hasSeq ? g : null, rec.hgvsP)
    }
  } else if (rec.hgvsC) {
    const h = parseC(`${rec.transcript ?? rec.gene}:${rec.hgvsC}`)
    const bad = h && orderProblem(h)
    if (bad) problems.push(`${bad} Not annotated.`)
    else if (h) ann = annotateHgvs({ ...h, ...(rec.gene ? { gene: rec.gene } : {}) }, rec.hgvsP)
    const other = rec.transcript && rec.gene ? txs.find(t => t.gene === rec.gene) : undefined
    if (other && ann) ann.notes = ann.notes.filter(n => !n.startsWith('No transcript model'))
    if (other) problems.push(`${rec.transcript} is not loaded. The loaded model for ${rec.gene} is ${other.id}, and c. numbers belong to one transcript, so this variant is read from its description alone. Import a GTF with ${rec.transcript}, or give genomic coordinates.`)
  }
  if (!ann && !problems.length) problems.push(rec.g ? 'No transcript model covers this position. Load the gene on the Genes tab (Ensembl, GRCh38) or import a GTF subset.' : 'Could not annotate this variant.')
  const splice = tx && g && ann && hasReference(genome, tx) ? spliceReport(tx, genome, g) : null
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
