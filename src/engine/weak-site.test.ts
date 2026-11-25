import { test } from 'node:test'
import assert from 'node:assert/strict'
import { MINI, MINI_G, read } from './mini.test-data.ts'
import { spliceReport } from './splice.ts'
import { cToG, refAt, type Transcript } from './transcript.ts'

// Exon 2 is 6 nt without real splice sites: the native scores are negative.
const ODD: Transcript = { id: 'NM_900003.1', gene: 'ODD', chrom: 'chr9', strand: 1, exons: [[1016, 1105], [1131, 1136], [1151, 1207]], cdsStart: 1016, cdsEnd: 1207, source: 'test' }
const rep = (t: Transcript, pos: number, alt: string) => spliceReport(t, MINI_G, { chrom: 'chr9', pos, ref: refAt(MINI_G, 'chr9', pos, 1)!, alt })!

test('a change that raises the score of a weak native site disrupts nothing', () => {
  const a = read('c.97-3T>C', [ODD, MINI_G])
  assert.ok(a.splice!.canonical!.refScore < 0)
  assert.ok(a.splice!.canonical!.delta > 0)
  assert.equal(a.splice!.verdict, 'none')
  assert.match(a.splice!.summary, /No predicted effect/)
})

test('over every single-base change near the sites: no drop, no disruption; every drop below 0 still disrupts', () => {
  let drops = 0, rises = 0
  for (const t of [ODD, MINI]) for (const [s, e] of t.exons) for (let pos = s - 16; pos <= e + 8; pos++) for (const alt of 'ACGT') {
    if (refAt(MINI_G, 'chr9', pos, 1) === alt || refAt(MINI_G, 'chr9', pos, 1) === null) continue
    const r = spliceReport(t, MINI_G, { chrom: 'chr9', pos, ref: refAt(MINI_G, 'chr9', pos, 1)!, alt })
    if (!r?.canonical) continue
    if (r.canonical.delta >= 0) { rises++; assert.notEqual(r.verdict, 'disrupts'); assert.notEqual(r.verdict, 'weakens') }
    else if (r.canonical.altScore < 0) { drops++; assert.equal(r.verdict, 'disrupts') }
  }
  assert.ok(drops > 50 && rises > 50, `${drops} ${rises}`)
})

test('losing a real site is still abolished', () => {
  assert.equal(rep(MINI, cToG(MINI, { pos: 75, offset: 1 }), 'A').verdict, 'disrupts')
  assert.equal(rep(MINI, cToG(MINI, { pos: 16, offset: -1 }), 'A').verdict, 'disrupts')
})
