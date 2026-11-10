// Local splice-site strength: position weight matrices for the 5' donor (9-mer, exon -3..intron +6) and the
// 3' acceptor (16-mer, intron -15..exon +1), scored as log2 odds against a uniform background. This is a
// MaxEnt-style approximation built from published consensus frequencies, not MaxEntScan itself.
import { refAt, trim, txExons, type Genome, type GVariant, type Transcript } from './transcript.ts'
import { revcomp } from './seq.ts'

type Col = [number, number, number, number] // A C G T
const DONOR: Col[] = [
  [.33, .36, .18, .13], [.60, .13, .14, .13], [.08, .04, .80, .08],
  [.001, .001, .997, .001], [.001, .001, .001, .997],
  [.59, .03, .35, .03], [.70, .08, .12, .10], [.07, .06, .80, .07], [.16, .15, .19, .50],
]
const PPT: Col = [.10, .33, .10, .47]
const ACCEPTOR: Col[] = [
  ...Array.from({ length: 10 }, () => PPT),
  [.22, .30, .10, .38], [.25, .30, .18, .27], [.05, .75, .01, .19],
  [.997, .001, .001, .001], [.001, .001, .997, .001],
  [.25, .14, .48, .13],
]
const IDX: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 }

function score(pwm: Col[], s: string): number {
  if (s.length !== pwm.length) return NaN
  let x = 0
  for (let i = 0; i < s.length; i++) {
    const k = IDX[s[i]]
    x += Math.log2((k === undefined ? .25 : pwm[i][k]) / .25)
  }
  return Math.round(x * 100) / 100
}
export const donorScore = (nine: string) => score(DONOR, nine)
export const acceptorScore = (sixteen: string) => score(ACCEPTOR, sixteen)
export const DONOR_MAX = Math.round(DONOR.reduce((s, c) => s + Math.log2(Math.max(...c) / .25), 0) * 100) / 100
export const ACCEPTOR_MAX = Math.round(ACCEPTOR.reduce((s, c) => s + Math.log2(Math.max(...c) / .25), 0) * 100) / 100

export interface SiteDelta { site: 'donor' | 'acceptor'; exon: number; refScore: number; altScore: number; delta: number; pct: number; refSeq: string; altSeq: string }
export interface SpliceReport {
  canonical: SiteDelta | null
  /** Strongest new or strengthened site of either type within the scanned window. */
  cryptic: (SiteDelta & { offset: number }) | null
  verdict: 'disrupts' | 'weakens' | 'creates' | 'none'
  summary: string
}

/** Edited transcript-orientation window around a genomic interval, with the variant applied. */
function windowSeqs(t: Transcript, g: Genome, v: GVariant, lo: number, hi: number) {
  const ref = refAt(g, t.chrom, lo, hi - lo + 1)
  if (ref === null) return null
  const { start, ref: r, alt } = trim(v)
  const i = start - lo
  if (i < 0 || i + r.length > ref.length) return null
  const altG = ref.slice(0, i) + alt + ref.slice(i + r.length)
  return t.strand === 1 ? { ref, alt: altG } : { ref: revcomp(ref), alt: revcomp(altG) }
}

/** Scores the nearest native splice site and scans +/-25 nt for created sites. Needs the reference sequence. */
export function spliceReport(t: Transcript, g: Genome, v: GVariant): SpliceReport | null {
  const ex = txExons(t)
  const pos = trim(v).start
  // Nearest exon boundary (in genomic coordinates) and its site type.
  let best: { site: 'donor' | 'acceptor'; exon: number; g: number } | null = null
  for (let i = 0; i < ex.length; i++) {
    const [s, e] = ex[i]
    const acc = t.strand === 1 ? s : e, don = t.strand === 1 ? e : s
    if (i > 0 && (!best || Math.abs(pos - acc) < Math.abs(pos - best.g))) best = { site: 'acceptor', exon: i + 1, g: acc }
    if (i < ex.length - 1 && (!best || Math.abs(pos - don) < Math.abs(pos - best.g))) best = { site: 'donor', exon: i + 1, g: don }
  }
  if (!best || !refAt(g, t.chrom, pos, 1)) return null
  const W = 40
  const lo = Math.min(pos, best.g) - W, hi = Math.max(pos, best.g) + W
  const w = windowSeqs(t, g, v, lo, hi)
  if (!w) return null
  const d = w.alt.length - w.ref.length
  // Transcript-oriented window index of the first changed base and of the base just after the REF span.
  const tIdx = (gp: number) => t.strand === 1 ? gp - lo : hi - gp
  const refLen = trim(v).ref.length
  const vStart = t.strand === 1 ? tIdx(pos) : refLen ? tIdx(pos + refLen - 1) : tIdx(pos) + 1
  const vEnd = vStart + refLen
  const toAlt = (k: number) => k >= vEnd ? k + d : k
  const nat = tIdx(best.g)
  const site = (kind: 'donor' | 'acceptor', start: number): SiteDelta => {
    const len = kind === 'donor' ? 9 : 16, fn = kind === 'donor' ? donorScore : acceptorScore
    const rs = w.ref.slice(start, start + len), as = w.alt.slice(toAlt(start), toAlt(start) + len)
    const r = fn(rs), al = fn(as)
    return { site: kind, exon: best!.exon, refScore: r, altScore: al, delta: round(al - r), pct: pct(r, al), refSeq: rs, altSeq: as }
  }
  const nativeStart = best.site === 'donor' ? nat - 2 : nat - 15
  const canonical = site(best.site, nativeStart)
  // Cryptic sites: alt windows overlapping the change that score well above their reference counterpart.
  let cryptic: SpliceReport['cryptic'] = null
  for (const kind of ['donor', 'acceptor'] as const) {
    const len = kind === 'donor' ? 9 : 16, fn = kind === 'donor' ? donorScore : acceptorScore
    for (let k = Math.max(0, vStart - len + 1); k < vStart + Math.max(alt(v), 1); k++) {
      if (kind === best.site && k === toAlt(nativeStart)) continue
      const as = w.alt.slice(k, k + len)
      const rk = k >= vStart + Math.max(alt(v), 0) ? k - d : k
      const rs = w.ref.slice(rk, rk + len)
      if (as.length !== len || rs.length !== len) continue
      const al = fn(as), r = fn(rs)
      if (al > 3 && al - r > 2 && (!cryptic || al - r > cryptic.delta)) cryptic = { site: kind, exon: best.exon, refScore: r, altScore: al, delta: round(al - r), pct: pct(r, al), refSeq: rs, altSeq: as, offset: k + (kind === 'donor' ? 2 : 15) - vStart }
    }
  }
  let verdict: SpliceReport['verdict'] = 'none'
  const c = canonical
  if (c.altScore < 0 || c.pct <= -60) verdict = 'disrupts'
  else if (c.pct <= -20) verdict = 'weakens'
  const cr = cryptic as SpliceReport['cryptic']
  if (verdict === 'none' && cr && cr.altScore >= c.altScore - 1) verdict = 'creates'
  const far = Math.abs(pos - best.g) > 50
  if (far && verdict !== 'creates') verdict = 'none'
  const summary = verdict === 'none' && far
    ? `${Math.abs(pos - best.g)} nt from the nearest splice site (${c.site} of exon ${c.exon}); no new site created.`
    : verdict === 'none'
    ? `Native ${c.site} of exon ${c.exon}: ${c.refScore} → ${c.altScore} bits (${c.pct >= 0 ? '+' : ''}${c.pct}%). No predicted effect.`
    : verdict === 'creates'
      ? `New ${cr!.site} site (${cr!.refScore} → ${cr!.altScore} bits) competes with the native ${c.site} (${c.altScore} bits).`
      : `Native ${c.site} of exon ${c.exon} ${verdict === 'disrupts' ? 'abolished' : 'weakened'}: ${c.refScore} → ${c.altScore} bits (${c.pct}%).`
  return { canonical, cryptic: cr, verdict, summary }
}

const round = (x: number) => Math.round(x * 100) / 100
const pct = (r: number, a: number) => r <= 0 ? 0 : Math.round(((a - r) / r) * 100)

const alt = (v: GVariant) => trim(v).alt.length
