import { test } from 'node:test'
import assert from 'node:assert/strict'
import { BOTH, read } from './mini.test-data.ts'

// Codons 1-5 M A E W Q, 6-25 D L K G T R S V H N P M A Y E L D K G F, 26-33 R D S L A P E *.
const p = (c: string, want: string, consequence?: string) => { for (const m of BOTH) { const a = read(c, m).ann!; assert.equal(a.hgvsP, want, c); if (consequence) assert.equal(a.consequence, consequence, c) } }

test('inserted residues that repeat those just before are dup', () => {
  p('c.31_33dup', 'p.(Arg11dup)', 'inframe_insertion')
  p('c.31_36dup', 'p.(Arg11_Ser12dup)', 'inframe_insertion')
  p('c.22_33dup', 'p.(Lys8_Arg11dup)')
  p('c.33_34insCGC', 'p.(Arg11dup)') // a synonymous copy
  p('c.4_6dup', 'p.(Ala2dup)')
  // Out of phase, CGT TCT becomes CGT TGT TCT: a duplicated base triplet that is not a duplicated residue.
  p('c.32_34dup', 'p.(Arg11_Ser12insCys)')
})

test('two copies are not a copy of the two residues before', () => {
  p('c.33_34insCGTCGT', 'p.(Arg11_Ser12insArgArg)')
})

test('other insertions stay insertions, and a new stop stays a stop', () => {
  p('c.33_34insGCT', 'p.(Arg11_Ser12insAla)', 'inframe_insertion')
  p('c.6_7insGATTAA', 'p.(Ala2_Glu3insAspTer)', 'stop_gained')
})

test('a copy of the last residue is not an extension; a lost stop still is', () => {
  p('c.94_96dup', 'p.(Glu32dup)', 'inframe_insertion')
  p('c.96_97insGCA', 'p.(Glu32_Ter33insAla)', 'inframe_insertion')
  for (const m of BOTH) for (const c of ['c.97T>C', 'c.98A>C', 'c.97_99del', 'c.97_98insC']) assert.equal(read(c, m).ann!.consequence, 'stop_lost', c)
})
