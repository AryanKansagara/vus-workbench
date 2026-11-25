import { test } from 'node:test'
import assert from 'node:assert/strict'
import { BOTH, read } from './mini.test-data.ts'
import { annotateHgvs } from './transcript.ts'
import { parseC } from './hgvs.ts'

// 5' UTR ...CCACC, then ATG GCT GAA.
const of = (c: string) => BOTH.map(m => read(c, m).ann!).map(a => `${a.consequence} ${a.hgvsP ?? '-'}`)
const both = (c: string, want: string) => assert.deepEqual(of(c), [want, want], c)

test('an insertion after ATG keeps the start codon', () => {
  both('c.3_4insA', 'frameshift p.(Ala2SerfsTer2)')
  both('c.3_4insGTT', 'inframe_insertion p.(Met1_Ala2insVal)')
  assert.match(of('c.3dup')[0], /^frameshift p\.\(Ala2GlyfsTer/)
  assert.deepEqual(of('c.3dup')[0], of('c.3dup')[1])
})

test('an insertion before ATG is 5\' UTR, also when it is named inside the codon', () => {
  both('c.-1_1insA', 'utr5 -')
  both('c.-1_1insGCTTT', 'utr5 -')
  both('c.1_2dup', 'utr5 -') // AT ATG: the same as inserting AT before ATG
  both('c.-2_1del', 'utr5 -') // CCACC ATG > CC ATG: the same as deleting c.-4_-2
  // Without a model the position says it too.
  assert.equal(annotateHgvs(parseC('X:c.-1_1insA')!).consequence, 'utr5')
  assert.equal(annotateHgvs(parseC('X:c.3_4insA')!).consequence, 'frameshift')
})

test('a deletion that brings the same base back keeps the start codon', () => {
  // ATG GCT: deleting c.3 equals deleting c.4.
  both('c.3del', of('c.4del')[0])
  assert.match(of('c.3del')[0], /^frameshift/)
})

test('a first codon that does change is start-lost', () => {
  for (const c of ['c.1A>G', 'c.2T>C', 'c.3G>A', 'c.1_2insC', 'c.2_3insA', 'c.2del', 'c.1_3del', 'c.-1_1del', 'c.-1_2delinsG', 'c.2_4del', 'c.1delinsTTT']) both(c, 'start_lost p.(Met1?)')
})

test('the protein change given with a 5\' UTR insertion is not used over the sequence', () => {
  assert.equal(read('c.-1_1insA', BOTH[0], undefined, 'p.Met1?').ann?.consequence, 'utr5')
  // ATG GCT GAA without TGGC is ATG AAT: the codon is whole again.
  assert.match(read('c.2_5del').ann!.hgvsP!, /^p\.\(Ala2AsnfsTer/)
})
