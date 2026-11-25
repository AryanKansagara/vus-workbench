import { test } from 'node:test'
import assert from 'node:assert/strict'
import { parseFasta, parseVcf } from './io.ts'

const head = (samples = '\tFORMAT\tP1') => `##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO${samples}\n`
const vcf = (row: string, samples?: string) => { const r = parseVcf(head(samples) + row, 'f', {}, () => 'x'); return [r.variants.length, r.issues.filter(i => i.level === 'error').map(i => i.line)] }
const R = '1\t5\t.\tC\tT\t60\tPASS\tAF=0'

test('a VCF row with fewer columns than its header is an error', () => {
  for (const cut of [R, R + '\t', R + '\tGT', R + '\tGT\t', R + '\t\t0/1']) assert.deepEqual(vcf(cut), [0, [3]], JSON.stringify(cut))
  assert.deepEqual(vcf(R + '\tGT\t0/1', '\tFORMAT\tP1\tP2'), [0, [3]])
  assert.deepEqual(vcf(`${R}\n1\t6\t.\tC\tT\t60\tPASS\t.\tGT\t0/1`), [1, [3]])
})

test('complete rows import as before', () => {
  assert.deepEqual(vcf(R + '\tGT\t0/1'), [1, []])
  assert.deepEqual(vcf(R + '\tGT\t0/1\r'), [1, []])
  assert.deepEqual(vcf(R + '\tGT:AD\t0/1'), [1, []]) // trailing sample sub-fields may be dropped (VCF 4.3, 1.6.2)
  assert.deepEqual(vcf(R, ''), [1, []])
  assert.deepEqual(vcf(R, '\t'), [1, []]) // a stray tab at the end of the header line
  assert.deepEqual(vcf(R, '\tFORMAT'), [0, [3]])
  assert.deepEqual(vcf(R + '\tGT\t0/1\textra'), [1, []])
})

const fa = (t: string) => { const r = parseFasta(t); return [Object.keys(r.genome), r.issues.filter(i => i.level === 'error').map(i => i.line)] }
test('a FASTA record must have the length its header states', () => {
  assert.deepEqual(fa('>chr9:11-20\nACGTA\nCGTAC\n'), [['chr9:11-20'], []])
  for (const t of ['>chr9:11-20\nACGTACGTA\n', '>chr9:11-20\nACGTACGTACG\n', '>chr9:11-20\n', '>chr9:20-11\nACGT\n']) assert.deepEqual(fa(t), [[], [1]], t)
  assert.deepEqual(fa('>chr9:11-20\nACGT\n>chr8:1-4\nACGT\n'), [['chr8'], [1]])
  assert.deepEqual(fa('>chr9:1,011-1,014\nACGT'), [['chr9:1011-1014'], []])
})

test('an empty FASTA record is an error, whatever its header', () => {
  assert.deepEqual(fa('>chr9\n'), [[], [1]])
  assert.deepEqual(fa('>chr9\n\n>chr8\nACGT\n'), [['chr8'], [1]])
  assert.deepEqual(fa('>chr9\nACGT'), [['chr9'], []])
  assert.deepEqual(fa(''), [[], [1]])
})
