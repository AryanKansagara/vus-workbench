// HGVS c. and p. parsing and formatting for the subset used in clinical reports.
import { AA1, AA3 } from './seq.ts'

/** A coding DNA position: `pos` is negative in the 5' UTR, `star` marks 3' UTR numbering, `offset` is the intronic distance. */
export interface CPos { pos: number; star?: boolean; offset: number }
export type CKind = 'sub' | 'del' | 'dup' | 'ins' | 'delins'
export interface HgvsC { transcript?: string; gene?: string; start: CPos; end?: CPos; kind: CKind; ref?: string; alt?: string }

const POS = String.raw`([*-]?)(\d+)(?:([+-])(\d+))?`
const RE_C = new RegExp(String.raw`^(?:([A-Za-z0-9_.]+)(?:\(([A-Za-z0-9-]+)\))?:)?c\.${POS}(?:_${POS})?(.*)$`)

function cpos(sign: string, n: string, osign?: string, off?: string): CPos {
  const p: CPos = { pos: sign === '-' ? -Number(n) : Number(n), offset: off ? (osign === '-' ? -Number(off) : Number(off)) : 0 }
  if (sign === '*') p.star = true
  return p
}

export function parseC(raw: string): HgvsC | null {
  const s = raw.trim().replace(/\s+/g, '')
  // Gene-only prefixes such as "BRCA2:c.68-7T>A" are accepted; "GENE(NM_...)" order is normalised too.
  const swapped = s.replace(/^([A-Za-z0-9-]+)\((N[MR]_[\d.]+|ENST[\d.]+)\):/, '$2($1):')
  const m = RE_C.exec(swapped)
  if (!m) return null
  let transcript = m[1], gene = m[2]
  if (transcript && !/^(N[MR]_|ENST|X[MR]_|LRG_|SYN)/i.test(transcript) && !gene) { gene = transcript; transcript = undefined }
  const start = cpos(m[3], m[4], m[5], m[6])
  const end = m[8] ? cpos(m[7], m[8], m[9], m[10]) : undefined
  const edit = m[11].toUpperCase()
  let e: RegExpExecArray | null
  const out = (o: Omit<HgvsC, 'start' | 'end' | 'transcript' | 'gene'>): HgvsC => ({ ...(transcript ? { transcript } : {}), ...(gene ? { gene } : {}), start, ...(end ? { end } : {}), ...o })
  if ((e = /^([ACGTN])>([ACGT])$/.exec(edit))) return end ? null : out({ kind: 'sub', ref: e[1], alt: e[2] })
  if ((e = /^DEL([ACGT]*)INS([ACGT]+)$/.exec(edit))) return out({ kind: 'delins', ...(e[1] ? { ref: e[1] } : {}), alt: e[2] })
  if ((e = /^DEL([ACGT]*|\d*)$/.exec(edit))) return out({ kind: 'del', ...(/^[ACGT]+$/.test(e[1]) ? { ref: e[1] } : {}) })
  if ((e = /^DUP([ACGT]*|\d*)$/.exec(edit))) return out({ kind: 'dup', ...(/^[ACGT]+$/.test(e[1]) ? { ref: e[1] } : {}) })
  if ((e = /^INS([ACGT]+)$/.exec(edit))) return end ? out({ kind: 'ins', alt: e[1] }) : null
  return null
}

export const formatPos = (p: CPos) => `${p.star ? '*' : ''}${p.pos}${p.offset > 0 ? '+' + p.offset : p.offset < 0 ? String(p.offset) : ''}`
export const samePos = (a: CPos, b: CPos) => a.pos === b.pos && !!a.star === !!b.star && a.offset === b.offset

/** Canonical form: deleted/duplicated bases are dropped, as the HGVS recommendations prefer. */
export function formatC(h: HgvsC, withPrefix = true): string {
  const prefix = withPrefix && (h.transcript || h.gene) ? `${h.transcript ?? h.gene}${h.transcript && h.gene ? `(${h.gene})` : ''}:` : ''
  const range = formatPos(h.start) + (h.end && !samePos(h.start, h.end) ? '_' + formatPos(h.end) : '')
  const edit = h.kind === 'sub' ? `${h.ref}>${h.alt}` : h.kind === 'del' ? 'del' : h.kind === 'dup' ? 'dup' : h.kind === 'ins' ? `ins${h.alt}` : `delins${h.alt}`
  return `${prefix}c.${range}${edit}`
}

/** Length in nucleotides of the reference span (introns included). Null when the span crosses UTR/intron numbering and cannot be counted from c. alone. */
export function spanLength(h: HgvsC): number | null {
  if (!h.end) return 1
  if (h.start.offset || h.end.offset || !!h.start.star !== !!h.end.star) return h.start.pos === h.end.pos && !!h.start.star === !!h.end.star ? h.end.offset - h.start.offset + 1 : null
  // There is no c.0: c.-1 is followed by c.1.
  const gap = h.start.pos < 0 && h.end.pos > 0 ? 1 : 0
  return h.end.pos - h.start.pos + 1 - gap
}

/** Net change in length (alt minus ref), or null when it cannot be told from the notation. */
export function lengthChange(h: HgvsC): number | null {
  const span = spanLength(h)
  switch (h.kind) {
    case 'sub': return 0
    case 'del': return span === null ? (h.ref ? -h.ref.length : null) : -span
    case 'dup': return span ?? (h.ref ? h.ref.length : null)
    case 'ins': return h.alt!.length
    case 'delins': return span === null ? null : h.alt!.length - span
  }
}

// ---------- Protein ----------

export type PKind = 'missense' | 'nonsense' | 'synonymous' | 'frameshift' | 'inframe_del' | 'inframe_ins' | 'inframe_delins' | 'start_lost' | 'stop_lost' | 'unknown'
export interface HgvsP { kind: PKind; pos?: number; ref?: string; alt?: string; fsTer?: number; predicted: boolean }

const AA = '(?:[A-Z][a-z]{2}|[A-Z*])'
const toOne = (s: string) => s.length === 1 ? (s === 'X' ? '*' : s) : (AA1[s] ?? 'X')

export function parseP(raw: string): HgvsP | null {
  let s = raw.trim().replace(/^.*:/, '').replace(/\s+/g, '')
  if (!s.startsWith('p.')) return null
  s = s.slice(2)
  const predicted = s.startsWith('(') && s.endsWith(')')
  if (predicted) s = s.slice(1, -1)
  if (s === '=' ) return { kind: 'synonymous', predicted }
  if (s === '?' || s === '0') return { kind: 'unknown', predicted }
  let m: RegExpExecArray | null
  if ((m = new RegExp(`^(${AA})(\\d+)=$`).exec(s))) return { kind: 'synonymous', pos: +m[2], ref: toOne(m[1]), alt: toOne(m[1]), predicted }
  if ((m = new RegExp(`^(${AA})(\\d+)(${AA})?fs(?:\\*|Ter|X)?(\\d+|\\?)?$`).exec(s))) {
    const p: HgvsP = { kind: 'frameshift', pos: +m[2], ref: toOne(m[1]), predicted }
    if (m[3]) p.alt = toOne(m[3])
    if (m[4] && m[4] !== '?') p.fsTer = +m[4]
    return p
  }
  if ((m = /^(Met|M)1(\?|[A-Z][a-z]{2}|[A-Z])$/.exec(s)) && m[2] !== 'Met' && m[2] !== 'M') return { kind: 'start_lost', pos: 1, ref: 'M', predicted }
  if ((m = new RegExp(`^(Ter|\\*|X)(\\d+)(${AA})ext`).exec(s))) return { kind: 'stop_lost', pos: +m[2], ref: '*', alt: toOne(m[3]), predicted }
  if ((m = new RegExp(`^(${AA})(\\d+)(${AA})$`).exec(s))) {
    const ref = toOne(m[1]), alt = toOne(m[3])
    const kind: PKind = alt === '*' ? 'nonsense' : ref === '*' ? 'stop_lost' : ref === alt ? 'synonymous' : +m[2] === 1 && ref === 'M' ? 'start_lost' : 'missense'
    return { kind, pos: +m[2], ref, alt, predicted }
  }
  if ((m = new RegExp(`^(${AA})(\\d+)(?:_(${AA})(\\d+))?(del|dup|ins|delins)([A-Za-z*]*)$`).exec(s))) {
    const op = m[5]
    const alt = m[6] ? (m[6].match(/[A-Z][a-z]{2}|[A-Z*]/g) ?? []).map(toOne).join('') : undefined
    const kind: PKind = op === 'del' ? 'inframe_del' : op === 'delins' ? (alt?.includes('*') ? 'nonsense' : 'inframe_delins') : 'inframe_ins'
    return { kind, pos: +m[2], ref: toOne(m[1]), ...(alt ? { alt } : {}), predicted }
  }
  return null
}

export function formatP(p: HgvsP): string {
  const w = (s: string) => p.predicted ? `p.(${s})` : `p.${s}`
  const a = (c?: string) => AA3[c ?? 'X'] ?? 'Xaa'
  switch (p.kind) {
    case 'synonymous': return w(p.pos ? `${a(p.ref)}${p.pos}=` : '=')
    case 'missense': case 'nonsense': return w(`${a(p.ref)}${p.pos}${a(p.alt)}`)
    case 'frameshift': return w(`${a(p.ref)}${p.pos}${p.alt ? a(p.alt) : ''}fs${p.fsTer ? `Ter${p.fsTer}` : ''}`)
    case 'start_lost': return w('Met1?')
    case 'stop_lost': return w(`Ter${p.pos}${a(p.alt)}extTer?`)
    case 'inframe_del': return w(`${a(p.ref)}${p.pos}del`)
    case 'inframe_ins': return w(`${a(p.ref)}${p.pos}ins${[...(p.alt ?? '')].map(a).join('')}`)
    case 'inframe_delins': return w(`${a(p.ref)}${p.pos}delins${[...(p.alt ?? '')].map(a).join('')}`)
    default: return w('?')
  }
}
