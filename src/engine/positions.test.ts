import { test } from 'node:test'
import assert from 'node:assert/strict'
import { BOTH, MINI, read } from './mini.test-data.ts'
import { analyse } from './pipeline.ts'
import { parseC } from './hgvs.ts'
import { EMPTY_EVIDENCE } from './acmg.ts'

const refused = (c: string, why: RegExp) => { for (const m of BOTH) { const a = read(c, m); assert.equal(a.ann, null, c); assert.equal(a.g, null, c); assert.match(a.problems.join(), why, c) } }
const fine = (c: string) => { for (const m of BOTH) { const a = read(c, m); assert.ok(a.ann, c); assert.ok(a.g, c); assert.deepEqual(a.problems, [], c) } }

test('an intron offset must start at an exon boundary next to an intron', () => {
  for (const c of ['c.30+5C>T', 'c.30-5C>T', 'c.15-1G>A', 'c.16+1G>A', 'c.99+3A>G', 'c.1-3A>G', 'c.*6+1A>G', 'c.-5-1A>G', 'c.70_74+3del', 'c.200+1G>A']) refused(c, /not an intron position|past the end/)
  for (const c of ['c.15+1G>A', 'c.16-1G>A', 'c.75+5del', 'c.16-2_16-1del', 'c.70_75+3del', 'c.75+39_76-1del']) fine(c)
})

test('an offset must lie within the intron, which is 40 nt', () => {
  refused('c.15+41A>G', /outside the intron, which is 40 nt/)
  refused('c.16-41A>G', /outside the intron/)
  for (const m of BOTH) assert.ok(read('c.15+40G>A', m).ann)
})

test('a c. position past the stop codon is refused when the sequence is known', () => {
  refused('c.100A>G', /past the end of the coding sequence/)
  fine('c.*1T>G')
  // Without sequence the model may lack its stop codon.
  assert.ok(analyse({ id: 'x', key: 'k', gene: 'MINI', hgvsC: 'c.99A>G', source: { file: '', line: 1, kind: 'csv' }, evidence: EMPTY_EVIDENCE }, [{ ...MINI, cdsEnd: 1191 }], {}).ann)
})

test('a reversed range is refused, with and without a model', () => {
  for (const c of ['c.90_80del', 'c.35_31dup', 'c.75+9_75+3del', 'c.76-3_75+3del', 'c.*3_95del', 'c.3_-2del']) refused(c, /runs backwards/)
  const bare = (hgvsC: string) => analyse({ id: 'x', key: 'k', gene: 'NONE', hgvsC, source: { file: '', line: 1, kind: 'csv' }, evidence: EMPTY_EVIDENCE }, [], {})
  for (const c of ['c.90_80del', 'c.75+9_75+3del', 'c.*3_95del', 'c.3_-2del', 'c.100-5_80+3del']) { assert.equal(bare(c).ann, null, c); assert.match(bare(c).problems.join(), /runs backwards/) }
  for (const c of ['c.80_90del', 'c.75+3_75+9del', 'c.76-80_75+60del', 'c.1-80_-1+60del', 'c.-5_3del', 'c.90_*3del', 'c.80_80del', 'c.75_76-3del', 'c.76-3_76del']) { assert.ok(bare(c).ann, c); assert.deepEqual(bare(c).problems, [], c) }
})

test('c.0 and an offset of 0 are not positions', () => {
  for (const c of ['c.0A>G', 'c.-0A>G', 'c.*0A>G', 'c.0_3del', 'c.-3_0del', 'c.15+0A>G', 'c.15_16-0del']) assert.equal(parseC(`MINI:${c}`), null, c)
  for (const c of ['c.10A>G', 'c.-10A>G', 'c.100+10A>G']) assert.ok(parseC(`MINI:${c}`), c)
})

test('stated dup and delins bases are checked', () => {
  // c.11 is G, c.31-33 is CGT.
  for (const c of ['c.11dupT', 'c.31_33dupCGA', 'c.31_33dupCG', 'c.31_33delCGAinsTT', 'c.11delTinsAA']) for (const m of BOTH) { const a = read(c, m); assert.equal(a.g, null, c); assert.match(a.problems.join(), /do not match/, c) }
  for (const c of ['c.11dupG', 'c.11dup', 'c.31_33dupCGT', 'c.31_33delCGTinsTT', 'c.31_33delinsTT', 'c.11delGinsAA']) fine(c)
})
