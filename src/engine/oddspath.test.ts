import { test } from 'node:test'
import assert from 'node:assert/strict'
import { annotateHgvs } from './transcript.ts'
import { parseC } from './hgvs.ts'
import { appliedOf, combine, DEFAULT_GENE, effective, EMPTY_EVIDENCE, suggest, type Strength } from './acmg.ts'

const ann = annotateHgvs(parseC('X:c.31C>T')!, 'p.(Arg11Cys)')
type R = 'abnormal' | 'normal'
const f = (oddsPath: number | null, result: R = 'abnormal') => ({ id: 'x', assay: 'Assay', result, oddsPath, note: '', editor: '', at: 0 })
const run = (...functional: ReturnType<typeof f>[]) => suggest(ann, { ...EMPTY_EVIDENCE, functional }, DEFAULT_GENE, null)
const met = (...l: ReturnType<typeof f>[]) => run(...l).filter(s => s.met).map(s => [s.code, s.strength])

test('OddsPath bands of Brnich 2019, at their boundaries', () => {
  const p: [number, Strength | null][] = [[351, 'VeryStrong'], [350, 'Strong'], [18.8, 'Strong'], [18.7, 'Moderate'], [4.3, 'Supporting'], [2.11, 'Supporting'], [2.1, null], [1, null], [0.48, null], [0.479, null], [0.04, null], [0.0001, null]]
  for (const [o, s] of p) assert.deepEqual(met(f(o)), s ? [['PS3', s]] : [], `damaging ${o}`)
  const b: [number, Strength | null][] = [[0.01, 'Strong'], [0.053, 'Moderate'], [0.23, 'Supporting'], [0.47, 'Supporting'], [0.48, null], [2.1, null], [3, null], [400, null]]
  for (const [o, s] of b) assert.deepEqual(met(f(o, 'normal')), s ? [['BS3', s]] : [], `normal ${o}`)
  assert.deepEqual(met(f(null)), [['PS3', 'Supporting']])
  assert.deepEqual(met(f(null, 'normal')), [['BS3', 'Supporting']])
})

test('a result that is not counted says why and gives no point', () => {
  const s = run(f(1), f(0.04))
  assert.equal(s.length, 1)
  assert.match(s[0].reason, /PS3 not applied.*OddsPath 1 is indeterminate.*OddsPath 0\.04 points to BS3/)
  assert.equal(combine(appliedOf(effective(s, {}))).points, 0)
})

test('a counted result stands beside one that is not', () => {
  const s = run(f(1), f(5), f(null))
  assert.deepEqual([s[0].met, s[0].strength], [true, 'Moderate'])
  assert.match(s[0].reason, /^2 functional damaging results.*Not counted: Assay \(OddsPath 1/)
})
