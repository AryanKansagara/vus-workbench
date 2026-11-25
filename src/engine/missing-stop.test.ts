import { test } from 'node:test'
import assert from 'node:assert/strict'
import { EMPTY_EVIDENCE } from './acmg.ts'
import { analyse } from './pipeline.ts'
import type { VariantRec } from './io.ts'
import { BUNDLED_GENOME as G, BUNDLED_TRANSCRIPTS as TXS } from './sample.ts'
import { cdsBounds, type Genome, type Transcript } from './transcript.ts'
import { revcomp } from './seq.ts'

const rec = (gene: string, hgvsC: string): VariantRec => ({ id: 'x', key: 'k', gene, hgvsC, source: { file: 'f', line: 2, kind: 'csv' }, evidence: EMPTY_EVIDENCE })
const cut = (t: Transcript, n = 3): Transcript => t.strand === 1 ? { ...t, cdsEnd: t.cdsEnd - n } : { ...t, cdsStart: t.cdsStart + n }

test('a model without its stop gives the same annotation as the full one, on either strand', () => {
  for (const [t, list] of [[TXS[0], ['c.103C>T', 'c.130_132del', 'c.560dup', 'c.270C>T', 'c.306+1G>A', 'c.*35C>A', 'c.-12G>T']], [TXS[1], ['c.100G>A', 'c.7_9del', 'c.162+1G>A']]] as const)
    for (const c of list) {
      const full = analyse(rec(t.gene, c), [t], G), short = analyse(rec(t.gene, c), [cut(t)], G)
      assert.deepEqual(short.ann, full.ann, `${t.gene} ${c}`)
      assert.deepEqual(short.g, full.g)
      assert.deepEqual(short.problems, full.problems)
    }
})

test('the stop may lie across a splice junction', () => {
  // ATG GCT GAA T | AA: the CDS rows end at GAA and the stop is TA + A over the intron.
  for (const strand of [1, -1] as const) {
    const s = 'CCATGGCTGAATA' + 'GTAAGTTTTTCAG' + 'ACCCC', at = (p: number) => strand === 1 ? p : s.length + 1 - p
    const e = (a: number, b: number): [number, number] => strand === 1 ? [a, b] : [at(b), at(a)]
    const exons = [e(1, 13), e(27, 31)].sort((x, y) => x[0] - y[0])
    const t: Transcript = { id: 'T', gene: 'J', chrom: '1', strand, exons, cdsStart: Math.min(at(3), at(11)), cdsEnd: Math.max(at(3), at(11)), source: '' }
    const g: Genome = { 1: { chrom: '1', start: 1, seq: strand === 1 ? s : revcomp(s) } }
    const a = analyse(rec('J', 'c.4_6del'), [t], g)
    assert.equal(a.ann?.hgvsP, 'p.(Ala2del)')
    assert.equal(cdsBounds(a.tx!).length, 12)
  }
})

test('no stop in the sequence: no protein prediction, and a message', () => {
  for (const n of [1, 2, 4, 6]) {
    const a = analyse(rec('SYNG1', 'c.130_132del'), [cut(TXS[0], n)], G)
    assert.equal(a.ann?.hgvsP, undefined)
    assert.match(a.problems.join(), /does not end in a stop codon/)
  }
})

test('a model without sequence is left alone', () => {
  const a = analyse(rec('SYNG1', 'c.130_132del'), [cut(TXS[0])], {})
  assert.equal(a.ann?.consequence, 'inframe_deletion')
  assert.doesNotMatch(a.problems.join(), /stop codon/)
})
