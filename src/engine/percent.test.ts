import { test } from 'node:test'
import assert from 'node:assert/strict'
import { parseGtf } from './io.ts'

const gff = (attr: string) => ['##gff-version 3', `1\tt\tgene\t1\t12\t.\t+\t.\tID=g1;${attr}`, '1\tt\tmRNA\t1\t12\t.\t+\t.\tID=T1;Parent=g1', '1\tt\texon\t1\t12\t.\t+\t.\tParent=T1', '1\tt\tCDS\t1\t12\t.\t+\t.\tParent=T1'].join('\n')
const gene = (attr: string) => parseGtf(gff(attr), 'x').transcripts[0]?.gene

test('a percent sign that is not an escape is kept as text', () => {
  for (const note of ['100% match', '%', '%%', '50%', '%ZZ', '%E0%A4%A', '%ff', 'a%2']) assert.equal(gene(`Note=${note};Name=MINI`), 'MINI', note)
  for (const name of ['100% match', '%', '%ZZ', '%ff', 'a%2']) assert.equal(gene(`Name=${name}`), name)
})

test('valid escapes are still decoded, also next to a stray percent sign', () => {
  assert.equal(gene('Name=A%3BB%2C%20C'), 'A;B, C')
  assert.equal(gene('Name=caf%C3%A9'), 'café')
  assert.equal(gene('Name=5% x%3By'), '5% x;y')
  assert.equal(gene('Name=%25'), '%')
})

test('a GTF with a percent sign and an equals sign in a quoted value imports', () => {
  const l = (t: string) => `1\tt\t${t}\t1\t12\t.\t+\t.\tgene_id "G"; transcript_id "T"; gene_name "G"; note "id=99% x";`
  assert.equal(parseGtf(`${l('exon')}\n${l('CDS')}`, 'x').transcripts.length, 1)
})
