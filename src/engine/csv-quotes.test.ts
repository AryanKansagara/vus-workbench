import { test } from 'node:test'
import assert from 'node:assert/strict'
import { parseCsv } from './csv.ts'
import { autoMap, LAB_FIELDS, validateLabCsv, validateVariantCsv, VARIANT_FIELDS } from './io.ts'

const H = 'variant,lab,classification,criteria\n'
const lab = (body: string) => { const c = parseCsv(H + body); return validateLabCsv(c, autoMap(c.rows[0] ?? [], LAB_FIELDS), 'f', () => 'x') }
const errors = (r: { issues: { level: string; line: number }[] }) => r.issues.filter(i => i.level === 'error').map(i => i.line)

test('a quote that never closes is an error at its line; the rows before it are kept', () => {
  const r = lab('V:c.1A>G,A,VUS,PM2\nV:c.1A>G,B,VUS,"PM2\nV:c.1A>G,C,Pathogenic,PVS1\n')
  assert.deepEqual(r.rows.map(x => x.lab), ['A'])
  assert.deepEqual(errors(r), [3])
  assert.match(r.issues[0].message, /never closed.*up to line/)
  // On the last line, at the end of the file, and just an opening quote.
  for (const tail of ['"PM2', '"PM2\n', '"', '"PM2""', '"PM2""\n']) assert.deepEqual([lab('V:c.1A>G,B,VUS,' + tail).rows.length, errors(lab('V:c.1A>G,B,VUS,' + tail))], [0, [2]], tail)
})

test('an unclosed quote that a later quote closes does not merge rows', () => {
  const r = lab('V:c.1A>G,A,VUS,"PM2\nV:c.1A>G,B,Pathogenic,"PVS1"\nV:c.1A>G,C,Benign,BA1\n')
  assert.deepEqual(r.rows.map(x => x.lab), ['C'])
  assert.deepEqual(errors(r), [2])
  assert.deepEqual(r.rows[0].criteria.map(c => c.code), ['BA1'])
  const s = lab('V:c.1A>G,A,VUS,"PM2\nV:c.1A>G,"Lab, B",Pathogenic,PVS1\n')
  assert.deepEqual([s.rows.length, errors(s)], [0, [2]])
})

test('a bad quote in the header refuses the file instead of reading the next row as the header', () => {
  for (const h of ['variant,lab,"classification,criteria\nV:c.1A>G,A,VUS,PM2\n', '"variant,"x,lab,classification\nvariant,lab,classification\nV:c.1A>G,A,VUS\n']) {
    const c = parseCsv(h)
    assert.deepEqual(c.rows, [])
    const r = validateLabCsv(c, autoMap([], LAB_FIELDS), 'f', () => 'x')
    assert.deepEqual([r.rows.length, errors(r)], [0, [1]])
    assert.deepEqual(errors(validateVariantCsv(c, autoMap([], VARIANT_FIELDS), 'f', {}, () => 'x')), [1])
  }
})

test('a variant list reports it too', () => {
  const c = parseCsv('gene,hgvs_c\nMINI,c.1A>G\nMINI,"c.2T>G\nMINI,c.3G>A\n')
  const r = validateVariantCsv(c, autoMap(c.rows[0], VARIANT_FIELDS), 'f', {}, () => 'x')
  assert.deepEqual([r.variants.length, errors(r)], [1, [3]])
})

test('well-formed quoting reads as before', () => {
  const rows = (s: string) => parseCsv(s).rows
  assert.deepEqual(rows('a,b\n"x, y","say ""hi"""\n'), [['a', 'b'], ['x, y', 'say "hi"']])
  assert.deepEqual(rows('a,b\n"two\nlines",z\r\n"","q"'), [['a', 'b'], ['two\nlines', 'z'], ['', 'q']])
  assert.deepEqual(rows('a,b\n5" probe,x"y\n'), [['a', 'b'], ['5" probe', 'x"y']])
  assert.deepEqual(rows('a\tb\n"Foo" syndrome\tz\n'), [['a', 'b'], ['Foo syndrome', 'z']])
  assert.deepEqual(rows('a,b\n"x" ,"y"  \n'), [['a', 'b'], ['x ', 'y  ']])
  assert.deepEqual(parseCsv('a,b\n"two\nlines",z\nq,r').lines, [1, 2, 4])
  assert.equal(parseCsv('a,b\n"x","y"').bad, undefined)
})
