import { test } from 'node:test'
import assert from 'node:assert/strict'
import { combine, parseCriteria } from './acmg.ts'

const read = (s: string) => parseCriteria(s).applied.map(a => `${a.code}_${a.strength}`).join(' ')

test('BP6 counts as benign supporting and PP5 as pathogenic supporting (Richards 2015)', () => {
  assert.equal(combine(parseCriteria('BP4, BP6').applied).points, -2)
  assert.equal(combine(parseCriteria('PP5').applied).points, 1)
  assert.equal(combine(parseCriteria('BA1 BS1 BP1').applied).points, -13)
})

test('a strength set apart from its code is kept', () => {
  for (const s of ['PS3 (Moderate)', 'PS3 - Moderate', 'PS3: Moderate', 'PS3 Moderate', 'PS3[moderate]', 'PS3_Moderate', 'PS3-Mod', 'PS3_M', 'PS3 – moderate']) assert.equal(read(s), 'PS3_Moderate', s)
  assert.equal(read('PP1 (Very Strong), PM2 supporting; BP4 strong'), 'PP1_VeryStrong PM2_Supporting BP4_Strong')
  assert.equal(read('PVS1 (Strong)\tPM2 (Supporting)'), 'PVS1_Strong PM2_Supporting')
})

test('a code stated as not met is left out and reported', () => {
  for (const s of ['PM2 not met', 'PM2 (not met)', 'PM2: not met', 'PM2_Supporting not met', 'PM2 not applied', 'PM2 N/A', 'not PM2', 'no PM2', 'PM2 unmet']) {
    const r = parseCriteria(`PVS1, ${s}, PP4`)
    assert.deepEqual(r.applied.map(a => a.code), ['PVS1', 'PP4'], s)
    assert.ok(r.unknown.some(u => u.includes('PM2')), s)
  }
})

test('lists read as before', () => {
  assert.equal(read('PVS1, PM2_Supporting; PP3_Mod PS3_P BS3_Str'), 'PVS1_VeryStrong PM2_Supporting PP3_Moderate PS3_Supporting BS3_Strong')
  assert.equal(read('PS3 PM2 PM1'), 'PS3_Strong PM2_Supporting PM1_Moderate')
  assert.equal(read('PS3, strong phenotype'), 'PS3_Strong')
  assert.equal(read('(PM2) (PP3)'), 'PM2_Supporting PP3_Supporting')
  assert.equal(read('PM2|PM2_Moderate'), 'PM2_Supporting')
  assert.deepEqual(parseCriteria('PS33 PVS12 xPM2').applied, [])
  assert.deepEqual(parseCriteria('').applied, [])
})
