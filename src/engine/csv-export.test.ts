import { test } from 'node:test'
import assert from 'node:assert/strict'
import { parseCsv, toCsv, uninert } from './csv.ts'
import { autoMap, LAB_FIELDS, validateLabCsv } from './io.ts'

const BAD = ['=HYPERLINK("http://x.test")', '=1+1', '+SUM(A1)', '-2+3', '@SUM(A1)', '\t=1', '\r=1', '=', '-', '--', '+1x', '-cmd|x']
const cell = (s: string) => parseCsv(toCsv([['h', 'x'], [s, 'x']])).rows[1][0]

test('text that starts like a formula is made inert', () => {
  for (const s of BAD) assert.equal(cell(s), `'${s}`, JSON.stringify(s))
})

test('numbers and ordinary text are unchanged', () => {
  for (const s of ['', 'a=b', 'c.31C>T', "'quoted", '1.0e-5', '-2', '+3', '-1.5e-3', 'Lab "A", B', 'two\nlines', ' =x', '0|1', '*']) assert.equal(cell(s), s)
  assert.equal(toCsv([['a', -2, 0.5, 1e-9, -10]]), 'a,-2,0.5,0,-10\n')
})

test('the apostrophe is taken off again, and only that one', () => {
  for (const s of BAD) assert.equal(uninert(cell(s)), s)
  for (const s of ["'quoted", "'", "''=x", "'-2", 'plain']) assert.equal(uninert(s), s)
})

test('an exported table re-imports under the same lab name', () => {
  const lab = '=HYPERLINK("http://x.test")'
  const csv = parseCsv(toCsv([['variant', 'gene', 'lab', 'classification', 'criteria', 'points', 'date_evaluated', 'curator'], ['BRCA2:c.1A>G', 'BRCA2', lab, 'VUS', 'PM2_Supporting', -3, '2026-01-02', '@me']]))
  const r = validateLabCsv(csv, autoMap(csv.rows[0], LAB_FIELDS), 'f', () => 'x')
  assert.deepEqual(r.issues, [])
  assert.equal(r.rows[0].lab, lab)
  assert.equal(csv.rows[1][5], '-3')
})
