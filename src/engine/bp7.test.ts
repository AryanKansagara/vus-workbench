import { test } from 'node:test'
import assert from 'node:assert/strict'
import { DEFAULT_GENE, EMPTY_EVIDENCE, suggest } from './acmg.ts'
import { annotateHgvs } from './transcript.ts'
import { parseC } from './hgvs.ts'
import { BOTH, read } from './mini.test-data.ts'

const codes = (c: string, spliceAi: number | null = 0) => suggest(annotateHgvs(parseC(`X:${c}`)!), { ...EMPTY_EVIDENCE, spliceAi }, DEFAULT_GENE, null).map(s => s.code)

test('BP7 at or beyond +7 and -21 with no predicted splice effect (Walker 2023)', () => {
  for (const c of ['c.75+7T>C', 'c.75+8T>C', 'c.75+9T>C', 'c.75+50T>C', 'c.76-21T>C', 'c.76-22T>C', 'c.-5+7T>C', 'c.75+7_75+9del', 'c.75+7_76-21del', 'c.75+8dup']) assert.deepEqual(codes(c), ['BP4', 'BP7'], c)
})

test('no BP7 inside +6 or -20, also when one end of a range is there', () => {
  for (const c of ['c.75+6T>C', 'c.75+3T>C', 'c.76-20T>C', 'c.76-9T>C', 'c.76-8T>C', 'c.76-3T>C', 'c.75+6_75+7insA', 'c.75+5_75+9del', 'c.75+9_76-10del', 'c.76-25_76-12del', 'c.75+30_76-5del']) assert.ok(!codes(c).includes('BP7'), c)
  for (const c of ['c.75+1T>C', 'c.75+2del', 'c.70_75+9del', 'c.40A>G', 'c.*5A>G']) assert.ok(!codes(c).includes('BP7'), c)
})

test('a predicted splice effect or no prediction gives no BP7', () => {
  assert.deepEqual(codes('c.75+7T>C', 0.2), ['PP3'])
  assert.deepEqual(codes('c.75+7T>C', 0.15), [])
  assert.deepEqual(codes('c.75+7T>C', null), [])
})

test('the same on either strand with a model', () => {
  for (const [t, g] of BOTH) {
    const a = read('c.75+7T>C', [t, g], { ...EMPTY_EVIDENCE, spliceAi: 0 })
    assert.deepEqual(suggest(a.ann!, a.rec.evidence, DEFAULT_GENE, a.splice).map(s => s.code), ['BP4', 'BP7'])
    assert.equal(a.ann!.consequence, 'splice_region')
  }
})
