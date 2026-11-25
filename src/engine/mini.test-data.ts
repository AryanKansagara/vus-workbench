// The hand-made MINI gene of e2e/fixtures.ts for unit tests, on either strand. Exon 1 c.1-15, exon 2 c.16-75, exon 3 c.76-99.
import { EMPTY_EVIDENCE, type Evidence } from './acmg.ts'
import { analyse, type Analysed } from './pipeline.ts'
import { revcomp } from './seq.ts'
import type { Genome, Transcript } from './transcript.ts'

export const MINI_SEQ = [
  'TTTTTTTTTT', 'CCACC', 'ATGGCTGAATGGCAA',
  'GTAAGTATTTATTTATTTACTAACTTTTCTTTCTTTTCAG',
  'GATCTGAAAGGTACCCGTTCTGTTCATAACCCAATGGCTTATGAACTGGATAAAGGTTTC',
  'GTAAGTTTTAATTTAATTTACTAACTTTTCTCTCTTTCAG',
  'CGTGATTCTCTGGCTCCAGAATAA', 'TTGCCT', 'AAAAAAAAAA',
].join('')
export const MINI: Transcript = { id: 'NM_900001.1', gene: 'MINI', chrom: 'chr9', strand: 1, exons: [[1011, 1030], [1071, 1130], [1171, 1200]], cdsStart: 1016, cdsEnd: 1194, source: 'test' }
export const MINI_G: Genome = { 'chr9:1001-1210': { chrom: 'chr9', start: 1001, seq: MINI_SEQ } }

/** A plus-strand model and its sequence turned round, so every c. description must read the same. */
export function flip(t: Transcript, g: Genome): [Transcript, Genome] {
  const r = Object.values(g)[0], last = r.start + r.seq.length - 1, at = (p: number) => r.start + last - p
  return [
    { ...t, strand: -1, exons: t.exons.map(([a, b]) => [at(b), at(a)] as [number, number]).reverse(), cdsStart: at(t.cdsEnd), cdsEnd: at(t.cdsStart) },
    { x: { ...r, seq: revcomp(r.seq) } },
  ]
}
export const BOTH: [Transcript, Genome][] = [[MINI, MINI_G], flip(MINI, MINI_G)]

export const read = (hgvsC: string, [t, g]: [Transcript, Genome] = [MINI, MINI_G], evidence: Evidence = EMPTY_EVIDENCE, hgvsP?: string): Analysed =>
  analyse({ id: 'x', key: 'k', gene: t.gene, hgvsC, ...(hgvsP ? { hgvsP } : {}), source: { file: 'f', line: 2, kind: 'csv' }, evidence }, [t], g)
