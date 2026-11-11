// Transcript model: genomic <-> c. mapping, HGVS from a genomic variant, protein effect, consequence and NMD.
import { formatC, lengthChange, parseP, type CPos, type HgvsC, type PKind } from './hgvs.ts'
import { aa3, revcomp, translate } from './seq.ts'

export interface Transcript {
  id: string
  gene: string
  chrom: string
  strand: 1 | -1
  /** Genomic, 1-based inclusive, ascending order regardless of strand. */
  exons: [number, number][]
  /** Lowest and highest genomic coordinate of the CDS, stop codon included. */
  cdsStart: number
  cdsEnd: number
  source: string
  synthetic?: boolean
}

/** Reference sequence for a contig, + strand, starting at genomic position `start`. */
export interface RefSeq { chrom: string; start: number; seq: string }
export type Genome = Record<string, RefSeq>

/** A genomic variant in VCF style: 1-based `pos` of the first REF base. */
export interface GVariant { chrom: string; pos: number; ref: string; alt: string }

export function refAt(g: Genome, chrom: string, pos: number, len: number): string | null {
  const r = g[chrom]
  if (!r) return null
  const i = pos - r.start
  if (i < 0 || i + len > r.seq.length) return null
  return r.seq.slice(i, i + len)
}

// ---------- Structure ----------

/** Exons in transcript order (5' to 3'). */
export const txExons = (t: Transcript) => t.strand === 1 ? t.exons : [...t.exons].reverse()
const exLen = (e: [number, number]) => e[1] - e[0] + 1
export const txLength = (t: Transcript) => t.exons.reduce((s, e) => s + exLen(e), 0)

/** 1-based transcript (spliced) position of an exonic genomic coordinate, or null if intronic/outside. */
export function gToTx(t: Transcript, g: number): number | null {
  let acc = 0
  for (const e of txExons(t)) {
    if (g >= e[0] && g <= e[1]) return acc + (t.strand === 1 ? g - e[0] + 1 : e[1] - g + 1)
    acc += exLen(e)
  }
  return null
}

/** Genomic coordinate of a transcript position; positions beyond the ends are extrapolated linearly. */
export function txToG(t: Transcript, tx: number): number {
  const ex = txExons(t)
  if (tx < 1) return t.strand === 1 ? ex[0][0] + tx - 1 : ex[0][1] - (tx - 1)
  let acc = 0
  for (const e of ex) {
    if (tx <= acc + exLen(e)) { const k = tx - acc; return t.strand === 1 ? e[0] + k - 1 : e[1] - k + 1 }
    acc += exLen(e)
  }
  const last = ex[ex.length - 1]
  const over = tx - acc
  return t.strand === 1 ? last[1] + over : last[0] - over
}

export function cdsBounds(t: Transcript) {
  const a = gToTx(t, t.strand === 1 ? t.cdsStart : t.cdsEnd)!
  const b = gToTx(t, t.strand === 1 ? t.cdsEnd : t.cdsStart)!
  return { start: a, end: b, length: b - a + 1, aaLength: (b - a + 1) / 3 - 1 }
}

function txToC(t: Transcript, tx: number): CPos {
  const { start, end } = cdsBounds(t)
  if (tx < start) return { pos: tx - start, offset: 0 }
  if (tx > end) return { pos: tx - end, star: true, offset: 0 }
  return { pos: tx - start + 1, offset: 0 }
}

function cToTx(t: Transcript, c: CPos): number {
  const { start, end } = cdsBounds(t)
  if (c.star) return end + c.pos
  if (c.pos < 0) return start + c.pos
  return start + c.pos - 1
}

/** Where a genomic coordinate falls: exon index (transcript order), or the intron between two exons. */
export interface Locus { cpos: CPos; exon?: number; intron?: number; outside?: 'upstream' | 'downstream' }

export function gToC(t: Transcript, g: number): Locus {
  const ex = txExons(t)
  const tx = gToTx(t, g)
  if (tx !== null) return { cpos: txToC(t, tx), exon: ex.findIndex(e => g >= e[0] && g <= e[1]) }
  const lo = t.exons[0][0], hi = t.exons[t.exons.length - 1][1]
  if (g < lo || g > hi) {
    const up = t.strand === 1 ? g < lo : g > hi
    const d = up ? Math.abs(g - (t.strand === 1 ? lo : hi)) : Math.abs(g - (t.strand === 1 ? hi : lo))
    return { cpos: up ? txToC(t, 1 - d) : txToC(t, txLength(t) + d), outside: up ? 'upstream' : 'downstream' }
  }
  // Intronic: anchor to the nearer exon boundary; ties go to the donor side (c.N+d), as HGVS specifies.
  for (let i = 0; i < ex.length - 1; i++) {
    const a = ex[i], b = ex[i + 1]
    const donorG = t.strand === 1 ? a[1] : a[0]
    const accG = t.strand === 1 ? b[0] : b[1]
    const inside = t.strand === 1 ? g > donorG && g < accG : g < donorG && g > accG
    if (!inside) continue
    const d1 = Math.abs(g - donorG), d2 = Math.abs(accG - g)
    if (d1 <= d2) return { cpos: { ...txToC(t, gToTx(t, donorG)!), offset: d1 }, intron: i }
    return { cpos: { ...txToC(t, gToTx(t, accG)!), offset: -d2 }, intron: i }
  }
  throw new Error('unreachable')
}

export function cToG(t: Transcript, c: CPos): number {
  return txToG(t, cToTx(t, c)) + c.offset * t.strand
}

// ---------- Sequence ----------

/** Spliced mRNA (transcript orientation) or null if the reference is not loaded. */
export function mrna(t: Transcript, g: Genome): string | null {
  let s = ''
  for (const e of txExons(t)) {
    const r = refAt(g, t.chrom, e[0], exLen(e))
    if (r === null) return null
    s += t.strand === 1 ? r : revcomp(r)
  }
  return s
}

/** Shifts a trimmed indel as far as possible in one genomic direction (+1 right, -1 left). */
export function shiftIndel(g: Genome, chrom: string, start: number, ref: string, alt: string, dir: 1 | -1): { start: number; ref: string; alt: string } {
  if (ref && alt) return { start, ref, alt }
  let seq = ref || alt
  let s = start
  for (let guard = 0; guard < 10000; guard++) {
    if (dir === 1) {
      // Deletion occupies [s, s+len-1]; insertion sits before s.
      const next = refAt(g, chrom, ref ? s + seq.length : s, 1)
      if (!next || next !== seq[0]) break
      seq = seq.slice(1) + next
      s++
    } else {
      const prev = refAt(g, chrom, s - 1, 1)
      if (!prev || prev !== seq[seq.length - 1]) break
      seq = prev + seq.slice(0, -1)
      s--
    }
  }
  return ref ? { start: s, ref: seq, alt: '' } : { start: s, ref: '', alt: seq }
}

/** Removes shared leading then trailing bases. `start` is the genomic position of the first remaining REF base (for a pure insertion, the base after the insertion point). */
export function trim(v: GVariant): { start: number; ref: string; alt: string } {
  let { ref, alt } = v
  let start = v.pos
  while (ref && alt && ref[0] === alt[0]) { ref = ref.slice(1); alt = alt.slice(1); start++ }
  while (ref && alt && ref[ref.length - 1] === alt[alt.length - 1]) { ref = ref.slice(0, -1); alt = alt.slice(0, -1) }
  return { start, ref, alt }
}

/** HGVS c. for a genomic variant against a transcript, 3'-shifted in transcript orientation when sequence is available. */
export function gToHgvsC(t: Transcript, g: Genome, v: GVariant): HgvsC {
  let { start, ref, alt } = trim(v)
  if (!ref || !alt) ({ start, ref, alt } = shiftIndel(g, v.chrom, start, ref, alt, t.strand))
  const base = { transcript: t.id, gene: t.gene }
  const o = (s: string) => t.strand === 1 ? s : revcomp(s)
  const span = (a: number, b: number) => {
    const ca = gToC(t, a).cpos, cb = gToC(t, b).cpos
    return t.strand === 1 ? { start: ca, end: cb } : { start: cb, end: ca }
  }
  if (ref.length === 1 && alt.length === 1) return { ...base, start: gToC(t, start).cpos, kind: 'sub', ref: o(ref), alt: o(alt) }
  if (!alt) return { ...base, ...(ref.length === 1 ? { start: gToC(t, start).cpos } : span(start, start + ref.length - 1)), kind: 'del', ref: o(ref) }
  if (!ref) {
    // Duplication when the inserted bases repeat the bases just 5' of the insertion (in transcript orientation).
    const L = alt.length
    if (t.strand === 1 && refAt(g, v.chrom, start - L, L) === alt) return { ...base, ...(L === 1 ? { start: gToC(t, start - 1).cpos } : span(start - L, start - 1)), kind: 'dup', ref: alt }
    if (t.strand === -1 && refAt(g, v.chrom, start, L) === alt) return { ...base, ...(L === 1 ? { start: gToC(t, start).cpos } : span(start, start + L - 1)), kind: 'dup', ref: o(alt) }
    return { ...base, ...span(start - 1, start), kind: 'ins', alt: o(alt) }
  }
  return { ...base, ...(ref.length === 1 ? { start: gToC(t, start).cpos } : span(start, start + ref.length - 1)), kind: 'delins', ref: o(ref), alt: o(alt) }
}

/** Genomic variant (VCF style, with anchor base for indels) from HGVS c.; needs the reference. */
export function hgvsCToG(t: Transcript, g: Genome, h: HgvsC): GVariant | null {
  const a = cToG(t, h.start), b = h.end ? cToG(t, h.end) : a
  const lo = Math.min(a, b), hi = Math.max(a, b)
  const o = (s: string) => t.strand === 1 ? s : revcomp(s)
  const refSpan = refAt(g, t.chrom, lo, hi - lo + 1)
  if (refSpan === null) return null
  const anchor = (pos: number, ref: string, alt: string): GVariant | null => {
    const pad = refAt(g, t.chrom, pos - 1, 1)
    return pad ? { chrom: t.chrom, pos: pos - 1, ref: pad + ref, alt: pad + alt } : null
  }
  switch (h.kind) {
    case 'sub': return refSpan === o(h.ref!) ? { chrom: t.chrom, pos: lo, ref: refSpan, alt: o(h.alt!) } : null
    case 'del': return h.ref && o(h.ref) !== refSpan ? null : anchor(lo, refSpan, '')
    case 'dup': return anchor(hi + 1, '', refSpan)
    case 'ins': return hi - lo === 1 ? anchor(hi, '', o(h.alt!)) : null
    case 'delins': return { chrom: t.chrom, pos: lo, ref: refSpan, alt: o(h.alt!) }
  }
}

// ---------- Protein effect ----------

export interface ProteinEffect { hgvsP: string; ptcAa?: number }

/** Protein consequence for a variant fully inside exons, computed by translating the edited transcript. */
export function proteinEffect(t: Transcript, g: Genome, v: GVariant): ProteinEffect | null {
  const m = mrna(t, g)
  if (!m) return null
  const { start: s, ref, alt } = trim(v)
  const lo = ref ? s : s - 1, hi = ref ? s + ref.length - 1 : s
  const ta = gToTx(t, lo), tb = gToTx(t, hi)
  if (ta === null || tb === null || Math.abs(tb - ta) !== hi - lo) return null
  const cds = cdsBounds(t)
  const [x, y] = t.strand === 1 ? [ta, tb] : [tb, ta]
  const altO = t.strand === 1 ? alt : revcomp(alt)
  const edited = ref ? m.slice(0, x - 1) + altO + m.slice(y) : m.slice(0, x) + altO + m.slice(x)
  if (y < cds.start || x > cds.end) return null
  if (x < cds.start + 3 && y >= cds.start) return { hgvsP: 'p.(Met1?)' }
  const refP = translate(m.slice(cds.start - 1, cds.end), false)
  const altP = translate(edited.slice(cds.start - 1), true)
  let i = 0
  while (i < refP.length && i < altP.length && refP[i] === altP[i]) i++
  if (i >= refP.length || (i === refP.length - 1 && altP[i] === '*' && refP[i] === '*')) {
    const pos = Math.min(Math.floor(((ref ? x : x + 1) - cds.start) / 3) + 1, refP.length)
    return { hgvsP: `p.(${aa3(refP[pos - 1])}${pos}=)` }
  }
  const dl = altO.length - (ref ? y - x + 1 : 0)
  const at = (k: number) => aa3(refP[k]) + (k + 1)
  if (refP[i] === '*') return { hgvsP: `p.(Ter${i + 1}${aa3(altP[i] ?? 'X')}extTer?)` }
  if (altP[i] === '*') return { hgvsP: `p.(${at(i)}Ter)`, ptcAa: i + 1 }
  if (dl % 3 !== 0) {
    const ter = altP.indexOf('*', i)
    return { hgvsP: `p.(${at(i)}${aa3(altP[i])}fs${ter >= 0 ? `Ter${ter - i + 1}` : ''})`, ptcAa: ter >= 0 ? ter + 1 : undefined }
  }
  // In-frame: trim the common prefix (done) then suffix on the protein level.
  const rTail = refP.slice(i), aTail = altP.slice(0, altP.indexOf('*') + 1 || undefined).slice(i)
  let k = 0
  while (k < rTail.length && k < aTail.length && rTail[rTail.length - 1 - k] === aTail[aTail.length - 1 - k]) k++
  const delSeg = rTail.slice(0, rTail.length - k), insSeg = aTail.slice(0, aTail.length - k)
  const rng = (n: number) => n === 1 ? at(i) : `${at(i)}_${at(i + n - 1)}`
  if (delSeg.length === 1 && insSeg.length === 1) return { hgvsP: `p.(${at(i)}${aa3(insSeg)})` }
  if (!insSeg) return { hgvsP: `p.(${rng(delSeg.length)}del)` }
  if (!delSeg) return { hgvsP: `p.(${at(i - 1)}_${at(i)}ins${aa3(insSeg)})` }
  return { hgvsP: `p.(${rng(delSeg.length)}delins${aa3(insSeg)})`, ...(insSeg.includes('*') ? { ptcAa: i + insSeg.indexOf('*') + 1 } : {}) }
}

// ---------- Consequence ----------

export type Consequence =
  | 'frameshift' | 'stop_gained' | 'splice_donor' | 'splice_acceptor' | 'start_lost' | 'stop_lost'
  | 'inframe_deletion' | 'inframe_insertion' | 'inframe_delins' | 'missense' | 'splice_region'
  | 'synonymous' | 'coding_sequence' | 'utr5' | 'utr3' | 'intronic' | 'upstream' | 'downstream'

export const CONSEQUENCE_LABEL: Record<Consequence, string> = {
  frameshift: 'Frameshift', stop_gained: 'Nonsense', splice_donor: 'Splice donor', splice_acceptor: 'Splice acceptor', start_lost: 'Start lost', stop_lost: 'Stop lost',
  inframe_deletion: 'In-frame deletion', inframe_insertion: 'In-frame insertion', inframe_delins: 'In-frame delins', missense: 'Missense', splice_region: 'Splice region',
  synonymous: 'Synonymous', coding_sequence: 'Coding (effect unknown)', utr5: "5' UTR", utr3: "3' UTR", intronic: 'Intronic', upstream: 'Upstream', downstream: 'Downstream',
}

/** Null variant types that enter the PVS1 decision tree. */
export const NULL_TYPES: Consequence[] = ['frameshift', 'stop_gained', 'splice_donor', 'splice_acceptor', 'start_lost']

export interface Annotation {
  hgvsC: string
  hgvsP?: string
  consequence: Consequence
  /** Also within the splice region (exonic 1-3 bp or intronic 3-8 bp from a junction) though the primary consequence is different. */
  spliceRegion: boolean
  /** Distance into the intron (signed c. offset), when intronic. */
  intronOffset?: number
  exon?: number
  intron?: number
  exonCount?: number
  /** Premature stop position (amino acid), when known. */
  ptcAa?: number
  /** Null = cannot tell (no transcript model). */
  nmd?: boolean | null
  /** Fraction of the protein removed by truncation or in-frame skipping of the affected exon. */
  removedFraction?: number | null
  /** For canonical splice variants: whether skipping the adjacent exon keeps the frame. */
  skipInFrame?: boolean | null
  /** Variant lies in the last exon or within 50 nt of the last exon-exon junction. */
  lastExonRegion?: boolean
  model: 'transcript' | 'hgvs-only'
  notes: string[]
}

const isCanonical = (off: number) => off === 1 || off === 2 || off === -1 || off === -2

/** Consequence from HGVS notation alone (no transcript model): position class, length change and the p. description. */
export function annotateHgvs(h: HgvsC, hgvsP?: string): Annotation {
  const notes: string[] = []
  const pos = [h.start, ...(h.end ? [h.end] : [])]
  const off = pos.map(p => p.offset).find(o => o !== 0) ?? 0
  const allIntronic = pos.every(p => p.offset !== 0)
  const p = hgvsP ? parseP(hgvsP) : null
  const dl = lengthChange(h)
  let consequence: Consequence
  let spliceRegion = false
  const canonical = pos.some(q => isCanonical(q.offset))
  const inRegion = pos.some(q => Math.abs(q.offset) >= 3 && Math.abs(q.offset) <= 8)
  if (canonical) consequence = pos.find(q => isCanonical(q.offset))!.offset > 0 ? 'splice_donor' : 'splice_acceptor'
  else if (allIntronic) consequence = inRegion ? 'splice_region' : 'intronic'
  else if (pos.every(q => !q.star && q.pos < 0)) consequence = 'utr5'
  else if (pos.every(q => q.star)) consequence = 'utr3'
  else if (p && p.kind !== 'unknown') {
    consequence = pKindToConsequence(p.kind)
    if (consequence === 'frameshift' && dl !== null && dl % 3 === 0) notes.push('p. says frameshift but the c. length change is a multiple of 3.')
  } else if (dl !== null && dl % 3 !== 0) consequence = 'frameshift'
  else if (h.kind === 'sub') { consequence = 'coding_sequence'; notes.push('Add the p. description or a transcript model with sequence to resolve missense / synonymous / nonsense.') }
  else consequence = h.kind === 'del' ? 'inframe_deletion' : h.kind === 'delins' ? 'inframe_delins' : 'inframe_insertion'
  if (inRegion && consequence !== 'splice_region') spliceRegion = true
  const a: Annotation = { hgvsC: formatC(h), consequence, spliceRegion, model: 'hgvs-only', notes, nmd: null, removedFraction: null, skipInFrame: null }
  if (hgvsP) a.hgvsP = hgvsP.replace(/^.*:/, '')
  if (off) a.intronOffset = off
  if (p?.kind === 'nonsense' && p.pos) a.ptcAa = p.pos
  if (p?.kind === 'frameshift' && p.pos && p.fsTer) a.ptcAa = p.pos + p.fsTer - 1
  notes.push('No transcript model loaded for this gene: exon structure, NMD and splice-strength are not computed. Import a GTF subset to enable them.')
  return a
}

function pKindToConsequence(k: PKind): Consequence {
  switch (k) {
    case 'nonsense': return 'stop_gained'
    case 'frameshift': return 'frameshift'
    case 'synonymous': return 'synonymous'
    case 'missense': return 'missense'
    case 'start_lost': return 'start_lost'
    case 'stop_lost': return 'stop_lost'
    case 'inframe_del': return 'inframe_deletion'
    case 'inframe_ins': return 'inframe_insertion'
    case 'inframe_delins': return 'inframe_delins'
    default: return 'coding_sequence'
  }
}

/** Full annotation against a transcript model. `v` is optional for c.-only input; with it, protein effect is computed from sequence. */
export function annotate(t: Transcript, g: Genome, h: HgvsC, v?: GVariant | null, givenP?: string): Annotation {
  const effect = v ? proteinEffect(t, g, v) : null
  const hgvsP = effect?.hgvsP ?? givenP
  const a = annotateHgvs(h, hgvsP)
  a.model = 'transcript'
  a.notes = a.notes.filter(n => !n.startsWith('No transcript model') && !(effect && n.startsWith('Add the p.')))
  if (!effect && !givenP && !v) a.notes.push('Reference sequence is not loaded for this contig; protein effect is taken from the p. description if given.')
  if (effect?.ptcAa) a.ptcAa = effect.ptcAa
  const ex = txExons(t)
  a.exonCount = ex.length
  const cds = cdsBounds(t)
  const lastJunctionTx = txLength(t) - exLen(ex[ex.length - 1])
  // Exon / intron numbering and the exonic splice region (1-3 bp from a junction).
  const g0 = cToG(t, h.start)
  const loc = gToC(t, g0)
  if (loc.exon !== undefined) {
    a.exon = loc.exon + 1
    const e = ex[loc.exon]
    const fromStart = t.strand === 1 ? g0 - e[0] : e[1] - g0
    const fromEnd = t.strand === 1 ? e[1] - g0 : g0 - e[0]
    if ((loc.exon > 0 && fromStart < 3) || (loc.exon < ex.length - 1 && fromEnd < 3)) a.spliceRegion = a.consequence !== 'splice_region' && !a.consequence.startsWith('splice')
  } else if (loc.intron !== undefined) a.intron = loc.intron + 1
  // NMD: a stop in the last exon or within 50 nt of the last junction escapes decay (the 50-55 nt rule).
  if (a.ptcAa) {
    const ptcTx = cds.start + (a.ptcAa - 1) * 3
    a.lastExonRegion = ex.length === 1 || ptcTx > lastJunctionTx - 50
    a.nmd = !a.lastExonRegion
    a.removedFraction = Math.max(0, (cds.aaLength - a.ptcAa + 1) / cds.aaLength)
  } else if (a.consequence === 'frameshift' || a.consequence === 'stop_gained') {
    a.nmd = null
    a.notes.push('Premature stop position unknown; give p. with fsTer or load the reference sequence.')
  }
  // Canonical splice: model as skipping of the adjacent exon.
  if ((a.consequence === 'splice_donor' || a.consequence === 'splice_acceptor') && loc.intron !== undefined) {
    const idx = a.consequence === 'splice_donor' ? loc.intron : loc.intron + 1
    const e = ex[idx]
    const ea = gToTx(t, t.strand === 1 ? e[0] : e[1])!, eb = ea + exLen(e) - 1
    const coding = Math.max(0, Math.min(eb, cds.end) - Math.max(ea, cds.start) + 1)
    a.skipInFrame = coding % 3 === 0
    // Out-of-frame skipping of one of the last two exons puts the new stop in the NMD-escape zone.
    a.lastExonRegion = idx >= ex.length - 2
    a.nmd = a.skipInFrame ? false : !a.lastExonRegion
    a.removedFraction = a.skipInFrame ? coding / 3 / cds.aaLength : Math.max(0, (cds.end - Math.max(ea, cds.start) + 1) / 3 / cds.aaLength)
    a.notes.push(`Predicted to skip exon ${idx + 1} (${exLen(e)} nt, ${a.skipInFrame ? 'in frame' : 'out of frame'}).`)
  }
  if (effect && a.hgvsP && a.consequence === 'coding_sequence') a.consequence = 'missense'
  return a
}

export function hgvsLabel(a: Annotation) { return a.hgvsP ? `${a.hgvsC} ${a.hgvsP}` : a.hgvsC }
