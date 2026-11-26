import { test } from 'node:test'
import assert from 'node:assert/strict'
import { parseNum, parseCsv } from './csv.ts'
import { autoMap, parseVcf, spliceAiMax, validateVariantCsv, VARIANT_FIELDS } from './io.ts'

let n = 0
const id = () => `id${n++}`
const vcf = (info: string, alt = 'T') => parseVcf(`##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t10\t.\tA\t${alt}\t.\tPASS\t${info}\n`, 'x.vcf', {}, id)

test('numbers: decimal comma after a leading zero, scientific notation, and nothing guessed', () => {
  for (const [s, v] of [['0,0001', 0.0001], ['0,85', 0.85], ['0.5', 0.5], ['1e-5', 1e-5], ['.5', 0.5], [' 1 ', 1], ['-1', -1], ['0', 0], ['5E-06', 5e-6]] as const) assert.equal(parseNum(s), v, s)
  for (const s of ['1,000', '1.234,5', '0x1', '5%', '', 'null', '$1', '1 000', '1,5', 'Infinity', '0,1,2', '<0.001']) assert.ok(Number.isNaN(parseNum(s)), s)
})

test('a decimal-comma list keeps its frequencies and an ambiguous one is refused with its line', () => {
  const csv = parseCsv('gene;hgvs_c;popmax_af;revel\nBRCA2;c.1A>G;0,0001;0,85\nBRCA2;c.2T>G;1,000;\n')
  const r = validateVariantCsv(csv, autoMap(csv.rows[0], VARIANT_FIELDS), 'x', {}, id)
  assert.deepEqual(r.variants.map(v => [v.evidence.popmaxAf, v.evidence.revel]), [[0.0001, 0.85], [null, null]])
  assert.deepEqual(r.issues.map(i => [i.line, i.level]), [[3, 'warning']])
})

test('SpliceAI: missing scores are not low scores', () => {
  assert.ok(Number.isNaN(spliceAiMax('T|G|.|.|.|.|.|.|.|.', 'T')))
  assert.ok(Number.isNaN(spliceAiMax('T|G||||', 'T')))
  assert.ok(Number.isNaN(spliceAiMax('T|G|0.01|.|0.00|0.00|1|2|3|4', 'T')))
  assert.equal(spliceAiMax('T|G|0.43|.|.|.|1|2|3|4', 'T'), 0.43)
  assert.equal(spliceAiMax('T|G|0.00|0.00|0.00|0.00|1|2|3|4', 'T'), 0)
  assert.equal(spliceAiMax('T|G|.|.|.|.|.|.|.|.,T|H|0.10|0.02|0.00|0.00|1|2|3|4', 'T'), 0.1)
  const r = vcf('SpliceAI=T|G|.|.|.|.|.|.|.|.')
  assert.equal(r.variants[0].evidence.spliceAi, null)
  assert.match(r.issues.map(i => i.message).join(), /no usable score/)
})

test('VCF frequencies and scores outside 0 to 1 are refused with the line', () => {
  const r = vcf('AF_popmax=-1;REVEL=7')
  assert.deepEqual([r.variants[0].evidence.popmaxAf, r.variants[0].evidence.revel], [null, null])
  assert.deepEqual(r.issues.filter(i => i.level === 'warning').map(i => [i.line, i.message]), [[3, 'AF_popmax "-1" ignored (expected 0 to 1)'], [3, 'REVEL "7" ignored (expected 0 to 1)']])
  const dot = vcf('AF_popmax=.;REVEL=.')
  assert.equal(dot.variants[0].evidence.popmaxAf, null)
  assert.equal(dot.issues.filter(i => i.level === 'warning').length, 0)
  assert.deepEqual([vcf('AF_popmax=0').variants[0].evidence.popmaxAf, vcf('REVEL_score=0.5').variants[0].evidence.revel], [0, 0.5])
})
