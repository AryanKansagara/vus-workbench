import { test } from 'node:test'
import assert from 'node:assert/strict'
import { DEFAULT_GENE, pvs1 } from './acmg.ts'
import { flip, MINI_G, read } from './mini.test-data.ts'
import { annotate, type Genome, type Transcript } from './transcript.ts'
import { parseC } from './hgvs.ts'

// Coding 90 + 6 + 57 = 153 nt, 50 residues and the stop. Exon 2 is residues 31-32.
const ODD: Transcript = { id: 'NM_900003.1', gene: 'ODD', chrom: 'chr9', strand: 1, exons: [[1016, 1105], [1131, 1136], [1151, 1207]], cdsStart: 1016, cdsEnd: 1207, source: 'test' }
const at = (c: string, [t, g]: [Transcript, Genome], regions?: [number, number][]) => { const a = annotate(t, g, parseC(`X:${c}`)!); return { a, s: pvs1(a, { ...DEFAULT_GENE, criticalRegions: regions }).strength } }

test('skipping an exon in a critical region is PVS1_Strong even under 10% (Abou Tayoun 2018)', () => {
  for (const m of [[ODD, MINI_G], flip(ODD, MINI_G)] as [Transcript, Genome][]) for (const c of ['c.96+1T>G', 'c.91-1A>G', 'c.94_96+3del']) {
    assert.deepEqual(at(c, m).a.lostAa, [31, 32], c)
    for (const r of [[30, 35], [31, 31], [32, 40], [1, 31], [1, 50]] as [number, number][]) assert.equal(at(c, m, [r]).s, 'Strong', `${c} ${r}`)
    for (const r of [[1, 30], [33, 50]] as [number, number][]) assert.equal(at(c, m, [r]).s, 'Moderate', `${c} ${r}`)
    assert.equal(at(c, m, [[1, 30], [32, 33]]).s, 'Strong')
    assert.equal(at(c, m).s, 'Moderate')
  }
})

test('a lost frame takes everything to the end of the protein', () => {
  // Exon 2 of 7 nt: skipping it shifts the frame, and the new stop is in the last exon region.
  const t: Transcript = { ...ODD, exons: [[1016, 1105], [1131, 1137], [1152, 1207]] }
  const { a } = at('c.97+1T>G', [t, MINI_G])
  assert.equal(a.skipInFrame, false)
  assert.equal(a.lostAa![1], 50)
})

test('an exon without coding sequence loses no residue', () => {
  const t: Transcript = { ...ODD, cdsStart: 1131 }
  assert.equal(at('c.-1+1G>A', [t, MINI_G], [[1, 50]]).a.lostAa, undefined)
})

test('NMD still gives full strength and a truncation is read from its stop as before', () => {
  assert.equal(pvs1(read('c.11G>A').ann!, DEFAULT_GENE).strength, 'VeryStrong')
  assert.equal(pvs1(read('c.94G>T').ann!, { ...DEFAULT_GENE, criticalRegions: [[32, 32]] }).strength, 'Strong')
  assert.equal(pvs1(read('c.94G>T').ann!, { ...DEFAULT_GENE, criticalRegions: [[1, 31]] }).strength, 'Moderate')
})
