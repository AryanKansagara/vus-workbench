import { test } from 'node:test'
import assert from 'node:assert/strict'
import { annotateHgvs } from './transcript.ts'
import { parseC } from './hgvs.ts'
import { DEFAULT_GENE, EMPTY_EVIDENCE, suggest, type ClassName, type GeneCfg, type Phase, type PhasePartner } from './acmg.ts'

const ann = annotateHgvs(parseC('X:c.31C>T')!, 'p.(Arg11Cys)')
const q = (cls: ClassName, phase: Phase): PhasePartner => ({ label: `${cls} ${phase}`, cls, phase })
const P = 'Pathogenic', LP = 'Likely pathogenic'
const codes = (inheritance: GeneCfg['inheritance'], ...partners: PhasePartner[]) => suggest(ann, EMPTY_EVIDENCE, { ...DEFAULT_GENE, inheritance }, null, partners).filter(s => s.met).map(s => `${s.code} ${s.strength}`)

test('PM3 follows SVI v1.0 Table 1 and 2 for one proband', () => {
  for (const c of [P, LP] as const) assert.deepEqual(codes('AR', q(c, 'trans')), ['PM3 Moderate'])
  assert.deepEqual(codes('AR', q(P, 'unknown')), ['PM3 Supporting'])
  assert.deepEqual(codes('AR', q(LP, 'unknown')), [])
  assert.deepEqual(codes('AR', q('VUS', 'trans')), [])
  assert.deepEqual(codes('AR', q(P, 'homozygous')), [])
  assert.deepEqual(codes('AR'), [])
})

test('the best partner wins in every order', () => {
  const list = [q(LP, 'unknown'), q(LP, 'trans'), q(P, 'unknown'), q('Benign', 'cis')]
  for (let i = 0; i < list.length; i++) for (const l of [[...list.slice(i), ...list.slice(0, i)], [...list.slice(i), ...list.slice(0, i)].reverse()]) assert.deepEqual(codes('AR', ...l), ['PM3 Moderate'])
  // A partner on each haplotype: both observations are made.
  for (const l of [[q(P, 'cis'), q(P, 'trans')], [q(P, 'trans'), q(P, 'cis')]]) assert.deepEqual(codes('AR', ...l), ['PM3 Moderate', 'BP2 Supporting'])
})

test('BP2 in trans is for dominant genes only, in cis for any', () => {
  assert.deepEqual(codes('AD', q(P, 'trans')), ['BP2 Supporting'])
  assert.deepEqual(codes('XL', q(P, 'trans')), [])
  for (const i of ['AD', 'AR', 'XL'] as const) assert.deepEqual(codes(i, q(LP, 'cis')), ['BP2 Supporting'])
  for (const i of ['AD', 'XL'] as const) assert.deepEqual(codes(i, q(P, 'unknown')), [])
})
