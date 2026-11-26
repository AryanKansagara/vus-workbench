import { test } from 'node:test'
import assert from 'node:assert/strict'
import { parseC, parseP } from './hgvs.ts'
import { annotate, annotateHgvs, gToHgvsC, hgvsCToG, proteinEffect } from './transcript.ts'
import { DEFAULT_GENE, pvs1 } from './acmg.ts'
import { pickAnnotation } from './io.ts'
import { BUNDLED_GENOME as G, BUNDLED_TRANSCRIPTS } from './sample.ts'

const [T1] = BUNDLED_TRANSCRIPTS
const run = (c: string) => {
  const v = hgvsCToG(T1, G, parseC(`${T1.id}:${c}`)!)!
  return { effect: proteinEffect(T1, G, v)!, ann: annotate(T1, G, gToHgvsC(T1, G, v), v) }
}

test('an in-frame insertion with a stop codon is nonsense with the stop in the description', () => {
  // SYNG1 codons 49-52: CTC AAA AAG CTG (Leu Lys Lys Leu). Codon 50 ends at c.150; GAT TAA is Asp Ter.
  const { effect, ann } = run('c.150_151insGATTAA')
  assert.equal(effect.hgvsP, 'p.(Lys50_Lys51insAspTer)')
  assert.equal(effect.ptcAa, 52)
  assert.equal(ann.consequence, 'stop_gained')
  assert.equal(ann.nmd, true)
  assert.equal(pvs1(ann, DEFAULT_GENE).strength, 'VeryStrong')
  // More bases after the stop change nothing.
  assert.equal(run('c.150_151insGATTAAGAT').effect.hgvsP, 'p.(Lys50_Lys51insAspTer)')
})

test('an in-frame delins that brings a stop in is nonsense', () => {
  // Codons 49-51 CTC AAA AAG replaced by GAT TGA CTC.
  const { effect, ann } = run('c.145_153delinsGATTGACTC')
  assert.equal(effect.hgvsP, 'p.(Leu49_Lys51delinsAspTer)')
  assert.equal(effect.ptcAa, 50)
  assert.equal(ann.consequence, 'stop_gained')
})

test('in-frame changes without a new stop are described as before', () => {
  assert.equal(run('c.130_132del').effect.hgvsP, 'p.(Asp44del)')
  assert.equal(run('c.150_151insGATGAT').effect.hgvsP, 'p.(Lys50_Lys51insAspAsp)')
  assert.equal(run('c.145_147delinsGATGAT').effect.hgvsP, 'p.(Leu49delinsAspAsp)')
  assert.equal(run('c.145_147delinsGATGAT').ann.consequence, 'inframe_delins')
})

test('p. delins and ins with Ter are read as nonsense at the right codon', () => {
  assert.deepEqual([parseP('p.(Lys51_Glu199delinsAspTer)')?.kind, parseP('p.(Lys51_Glu199delinsAspTer)')?.stopAt], ['nonsense', 52])
  assert.deepEqual([parseP('p.Lys51_Leu52insAspTer')?.kind, parseP('p.Lys51_Leu52insAspTer')?.stopAt], ['nonsense', 53])
  assert.equal(parseP('p.K51_L52delinsD*')?.stopAt, 52)
  assert.equal(parseP('p.Lys51_Leu52delinsArg')?.kind, 'inframe_delins')
  assert.equal(parseP('p.Lys51_Leu52delinsArg')?.alt, 'R')
  assert.equal(parseP('p.Lys51del')?.kind, 'inframe_del')
  assert.equal(parseP('p.Lys51dup')?.kind, 'inframe_ins')
  assert.equal(parseP('p.Arg51Ter')?.stopAt, undefined)
  assert.equal(annotateHgvs(parseC('X:c.151_156delinsGATTAA')!, 'p.(Lys51_Leu52delinsAspTer)').ptcAa, 52)
  assert.equal(annotateHgvs(parseC('X:c.103C>T')!, 'p.(Gln35Ter)').ptcAa, 35)
})

test('a VEP p. delins keeps its wording', () => {
  assert.equal(pickAnnotation('T|X|NP_1.1:p.Lys51_Leu52delinsAspTer', ['Allele', 'SYMBOL', 'HGVSp'], 'T', 1, 1)?.hgvsP, 'p.Lys51_Leu52delinsAspTer')
})
