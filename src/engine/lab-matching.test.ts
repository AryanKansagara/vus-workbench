import { test } from 'node:test'
import assert from 'node:assert/strict'
import { autoMap, dateKey, LAB_FIELDS, labVariantKey, parseCsv, validateLabCsv, type LabRow } from './io.ts'
import { reconcile } from './discrepancy.ts'

let n = 0
const read = (t: string) => { const c = parseCsv(t); return validateLabCsv(c, autoMap(c.rows[0], LAB_FIELDS), 'l.csv', () => `id${n++}`) }
const row = (o: Partial<LabRow>): LabRow => ({ id: `r${n++}`, variant: 'BRCA2:c.1A>G', lab: 'A', cls: 'VUS', criteria: [], line: 2, file: 'l.csv', ...o })

test('spellings of one change share a key', () => {
  for (const v of ['NM_007294.4:c.5266dupC', 'NM_007294.4:c.5266dup', 'NM_007294.3(BRCA1):c.5266dupc', 'c.5266dup', 'NM_007294.4(BRCA1):c.5266dup (p.Gln1756fs)']) assert.equal(labVariantKey(v, 'BRCA1'), 'BRCA1:c.5266dup', v)
  assert.equal(labVariantKey('CFTR:c.1521_1523delCTT'), 'CFTR:c.1521_1523del')
  assert.equal(labVariantKey('BRCA2:c.68-7t>a'), 'BRCA2:c.68-7T>A')
  assert.equal(labVariantKey('X:c.10delAinsTT'), 'X:c.10delinsTT')
  assert.notEqual(labVariantKey('X:c.10insA'), labVariantKey('X:c.10dup'))
  // Not a c. description: kept as written.
  assert.equal(labVariantKey('13:32340300:A:T'), '13:32340300:A:T')
})

test('rows saved before this change are grouped too', () => {
  const d = reconcile([row({ variant: 'BRCA1:c.5266dupC', cls: 'Pathogenic' }), row({ variant: 'BRCA1:c.5266dup', lab: 'B', cls: 'Benign' })])
  assert.deepEqual(d.map(x => [x.variant, x.severity]), [['BRCA1:c.5266dup', 'major']])
})

test('dates become sortable text and unclear ones say so', () => {
  const ok: [string, string][] = [['2024-01-09', '2024-01-09'], ['2024/1/9', '2024-01-09'], ['2024-01', '2024-01'], ['2024', '2024'], ['Jan 09, 2024', '2024-01-09'], ['September 18, 2019', '2019-09-18'], ['09 Jan 2024', '2024-01-09'], ['09-Jan-2024', '2024-01-09'], ['Jan 2024', '2024-01'], ['2024-02-29', '2024-02-29'], ['2024-01-09T10:00:00Z', '2024-01-09'], ['', ''], ['-', '']]
  for (const [a, b] of ok) assert.deepEqual(dateKey(a), { key: b, exact: true }, a)
  for (const [a, b] of [['03/04/2024', '2024'], ['2023-02-29', '2023'], ['2024-13-01', '2024'], ['45300', ''], ['soon', ''], ['Foo 12, 2024', '2024']]) assert.deepEqual(dateKey(a), { key: b, exact: false }, a)
})

test('the latest call wins whatever the date form', () => {
  const calls = (rows: Partial<LabRow>[]) => reconcile(rows.map(row))[0].calls.map(c => c.cls)
  assert.deepEqual(calls([{ cls: 'Benign', date: 'Jan 09, 2024' }, { cls: 'Pathogenic', date: 'Sep 18, 2019' }]), ['Benign'])
  assert.deepEqual(calls([{ cls: 'Benign', date: '2024-01-09' }, { cls: 'Pathogenic', date: 'Sep 18, 2019' }]), ['Benign'])
  assert.deepEqual(calls([{ cls: 'Benign', date: '31/01/2019' }, { cls: 'Pathogenic', date: '01/02/2024' }]), ['Pathogenic'])
  assert.deepEqual(calls([{ cls: 'Benign', date: '2024-01-09' }, { cls: 'Pathogenic' }]), ['Benign'])
})

test('an order that cannot be told is stated', () => {
  const d = reconcile([row({ cls: 'Benign', date: '03/04/2024' }), row({ cls: 'Pathogenic', date: '05/03/2024' })])[0]
  assert.match(d.summary, /dates do not say which is later; Pathogenic is shown/)
  assert.doesNotMatch(reconcile([row({ cls: 'Benign', date: '2023-01-01' }), row({ cls: 'Pathogenic', date: '2024-01-01' })])[0].summary, /dates do not say/)
  const r = read('variant,lab,classification,date\nBRCA2:c.1A>G,A,Benign,03/04/2024\nBRCA2:c.1A>G,B,Benign,-\n')
  assert.deepEqual(r.issues.map(i => [i.line, i.level]), [[2, 'warning']])
})
