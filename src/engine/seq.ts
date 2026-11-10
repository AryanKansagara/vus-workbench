// Sequence helpers: genetic code, reverse complement, translation, amino acid names.

const BASES = 'TCAG'
const AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
export const CODON: Record<string, string> = {}
for (let i = 0; i < 64; i++) CODON[BASES[i >> 4] + BASES[(i >> 2) & 3] + BASES[i & 3]] = AAS[i]

const COMP: Record<string, string> = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' }
export const revcomp = (s: string) => [...s].reverse().map(b => COMP[b] ?? 'N').join('')
export const complement = (s: string) => [...s].map(b => COMP[b] ?? 'N').join('')

/** Translates in frame 0. Stops after the first stop codon when `toStop` is set; incomplete trailing codons are dropped. */
export function translate(nt: string, toStop = false): string {
  let out = ''
  for (let i = 0; i + 3 <= nt.length; i += 3) {
    const aa = CODON[nt.slice(i, i + 3)] ?? 'X'
    out += aa
    if (toStop && aa === '*') break
  }
  return out
}

export const AA3: Record<string, string> = {
  A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly', H: 'His', I: 'Ile', L: 'Leu', K: 'Lys',
  M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser', T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': 'Ter', X: 'Xaa', U: 'Sec',
}
export const AA1: Record<string, string> = Object.fromEntries(Object.entries(AA3).map(([k, v]) => [v, k]))
export const aa3 = (s: string) => [...s].map(c => AA3[c] ?? 'Xaa').join('')

/** Deterministic PRNG (mulberry32) for the synthetic gene models. */
export function rng(seed: number) {
  let a = seed >>> 0
  return () => {
    a = (a + 0x6d2b79f5) >>> 0
    let t = a
    t = Math.imul(t ^ (t >>> 15), t | 1)
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61)
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296
  }
}
