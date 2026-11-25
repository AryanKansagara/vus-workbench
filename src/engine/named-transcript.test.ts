import { test } from 'node:test'
import assert from 'node:assert/strict'
import { EMPTY_EVIDENCE } from './acmg.ts'
import { analyse, evaluate } from './pipeline.ts'
import { DEFAULT_GENE } from './acmg.ts'
import type { VariantRec } from './io.ts'
import { BUNDLED_GENOME as G, BUNDLED_TRANSCRIPTS as TXS } from './sample.ts'

const rec = (o: Partial<VariantRec>): VariantRec => ({ id: 'x', key: 'k', gene: 'SYNG1', source: { file: 'f', line: 2, kind: 'csv' }, evidence: EMPTY_EVIDENCE, ...o })

test('a c. description on a transcript that is not loaded is not placed on another one', () => {
  for (const hgvsC of ['c.103C>T', 'c.130_132del', 'c.560dup']) {
    const a = analyse(rec({ transcript: 'NM_999999.1', hgvsC }), TXS, G)
    assert.equal(a.tx, null)
    assert.equal(a.g, null)
    assert.equal(a.ann?.model, 'hgvs-only')
    assert.equal(a.ann?.hgvsP, undefined)
    assert.equal(a.gene, 'SYNG1')
    assert.match(a.problems.join(), /NM_999999\.1 is not loaded.*SYN_000101\.1/)
  }
  assert.equal(evaluate(analyse(rec({ transcript: 'NM_999999.1', hgvsC: 'c.103C>T' }), TXS, G), undefined, DEFAULT_GENE).result.points, 0)
})

test('the same accession in another version, no transcript, or a label that is not an accession use the model', () => {
  for (const transcript of ['SYN_000101.1', 'SYN_000101.7', 'SYN_000101', undefined, 'canonical', 'MANE']) assert.equal(analyse(rec({ transcript, hgvsC: 'c.103C>T' }), TXS, G).ann?.hgvsP, 'p.(Gln35Ter)', String(transcript))
})

test('the Ensembl twin of a MANE Select model matches', () => {
  const tx = { ...TXS[0], id: 'NM_900009.1', source: 'Ensembl 116 MANE Select ENST00000999999.8 (GRCh38)' }
  assert.equal(analyse(rec({ transcript: 'ENST00000999999.7', hgvsC: 'c.103C>T' }), [tx], G).tx?.id, 'NM_900009.1')
  assert.equal(analyse(rec({ transcript: 'ENST00000999998.1', hgvsC: 'c.103C>T' }), [tx], G).tx, null)
})

test('genomic coordinates still use the gene model', () => {
  const a = analyse(rec({ transcript: 'NM_999999.1', g: { chrom: 'chrS1', pos: 100863, ref: 'C', alt: 'T' } }), TXS, G)
  assert.equal(a.ann?.hgvsC, 'SYN_000101.1(SYNG1):c.103C>T')
})

test('with no model at all nothing changes', () => {
  const a = analyse(rec({ gene: 'BRCA2', transcript: 'NM_000059.4', hgvsC: 'c.68-7T>A' }), TXS, G)
  assert.equal(a.ann?.consequence, 'splice_region')
  assert.deepEqual(a.problems, [])
})
