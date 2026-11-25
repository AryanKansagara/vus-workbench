import { test } from 'node:test'
import assert from 'node:assert/strict'
import { parseC } from './hgvs.ts'
import { annotate, annotateHgvs, type Annotation } from './transcript.ts'
import { DEFAULT_GENE, EMPTY_EVIDENCE, pvs1, suggest } from './acmg.ts'
import { BUNDLED_GENOME as G, BUNDLED_TRANSCRIPTS } from './sample.ts'

const [T1, T2] = BUNDLED_TRANSCRIPTS
const h = (c: string) => parseC(`X:${c}`)!
const only = (c: string) => annotateHgvs(h(c)).consequence
const at = (c: string, t = T1): Annotation => annotate(t, G, h(c))

test('a deletion that removes +1/+2 or -1/-2 is a splice variant whatever its ends are', () => {
  assert.equal(only('c.300_306+5del'), 'splice_donor')
  assert.equal(only('c.300_306+5delinsTT'), 'splice_donor')
  assert.equal(only('c.157-10_160del'), 'splice_acceptor')
  assert.equal(only('c.157-25_306+25del'), 'splice_acceptor')
  assert.equal(only('c.157-25_400-30del'), 'splice_acceptor')
  assert.equal(only('c.306+1_306+10del'), 'splice_donor')
  assert.equal(only('c.300_307-30del'), 'splice_donor')
  assert.equal(only('c.306+30_310del'), 'splice_acceptor')
  assert.equal(only('c.306+30_394-30del'), 'splice_acceptor')
})

test('deletions that stop short of the canonical bases are unchanged', () => {
  assert.equal(only('c.306+3_307-3del'), 'splice_region')
  assert.equal(only('c.306+30_307-30del'), 'intronic')
  assert.equal(only('c.157-25_157-12del'), 'intronic')
  assert.equal(only('c.-1+30_1-30del'), 'intronic')
  assert.equal(only('c.306+20_306+40del'), 'intronic')
  assert.equal(only('c.100_102del'), 'inframe_deletion')
  assert.equal(only('c.100_103del'), 'frameshift')
})

test('an unknown length is not called in-frame', () => {
  assert.equal(only('c.600_*5del'), 'stop_lost')
  assert.equal(only('c.600_*5dup'), 'coding_sequence')
  assert.equal(only('c.-40-3_-40-3del'), 'splice_region')
})

test('frame and size follow the exons that are lost', () => {
  // Exon 3 is 150 nt, c.157-306: 50 of 200 residues.
  for (const c of ['c.300_306+5del', 'c.157-10_160del', 'c.157-25_306+25del']) {
    const a = at(c)
    assert.equal(a.skipInFrame, true, c)
    assert.equal(a.removedFraction, 0.25, c)
    assert.equal(pvs1(a, DEFAULT_GENE).strength, 'Strong', c)
    assert.ok(!suggest(a, { ...EMPTY_EVIDENCE, spliceAi: 0 }, DEFAULT_GENE, null).some(s => ['PM4', 'PP3', 'BP4', 'BP7'].includes(s.code)), c)
  }
  // Exons 3 and 4: 150 + 87 = 237 nt = 79 of 200 residues.
  const two = at('c.157-25_393+25del')
  assert.ok(Math.abs(two.removedFraction! - 79 / 200) < 1e-9)
  assert.match(two.notes.join(' '), /lose exons 3, 4 \(in frame\)/)
})

test('minus strand: exon 2 of SYNG2 is 112 nt, out of frame, with decay', () => {
  const a = at('c.51-9_162+9del', T2)
  assert.equal(a.skipInFrame, false)
  assert.equal(a.nmd, true)
  assert.equal(pvs1(a, DEFAULT_GENE).strength, 'VeryStrong')
})
