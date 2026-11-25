import { test } from 'node:test'
import assert from 'node:assert/strict'
import { compareSnapshots, createBlankProject, evaluateProject, snapshotOf } from './projectData.ts'
import type { VariantRec } from './engine/io.ts'
import type { Snapshot, SnapRow } from './store.ts'

const r = (key: string, sample: string | undefined, points = 1): SnapRow => ({ key, ...(sample === undefined ? {} : { sample }), label: key, cls: points > 5 ? 'Likely pathogenic' : 'VUS', points, criteria: [`P${points}`] })
const s = (...rows: SnapRow[]): Snapshot => ({ id: 'x', label: 'x', createdAt: 0, rows })
const cmp = (a: Snapshot, b: Snapshot) => { const d = compareSnapshots(a, b); return [d.unchanged, d.reclassified.map(x => `${x.label} ${x.dPoints}`), d.added.map(x => x.sample ?? x.key), d.removed.map(x => x.sample ?? x.key)] }

test('each sample\'s row of one variant is compared', () => {
  assert.deepEqual(cmp(s(r('v', 'S1'), r('v', 'S2')), s(r('v', 'S1'), r('v', 'S2'))), [2, [], [], []])
  assert.deepEqual(cmp(s(r('v', 'S1'), r('v', 'S2')), s(r('v', 'S2'), r('v', 'S1', 9))), [1, ['v (S1) 8'], [], []])
  assert.deepEqual(cmp(s(r('v', 'S1')), s(r('v', 'S1'), r('v', 'S2'))), [1, [], ['S2'], []])
  assert.deepEqual(cmp(s(r('v', 'S1'), r('v', 'S2')), s(r('v', 'S2'))), [1, [], [], ['S1']])
  assert.deepEqual(cmp(s(r('v', '')), s(r('v', ''))), [1, [], [], []])
})

test('snapshots saved without samples still compare, with each other and with new ones', () => {
  const old = s(r('v', undefined), r('v', undefined, 2), r('w', undefined))
  assert.deepEqual(cmp(old, old), [3, [], [], []])
  assert.deepEqual(cmp(old, s(r('v', 'S1'), r('v', 'S2', 2), r('w', ''))), [3, [], [], []])
  assert.deepEqual(cmp(old, s(r('v', 'S1', 9), r('v', 'S2', 2), r('w', ''))), [2, ['v (S1) 8'], [], []])
  assert.deepEqual(cmp(s(r('v', undefined)), s(r('v', 'S1'), r('v', 'S2'))), [1, [], ['S2'], []])
  assert.deepEqual(cmp(s(), old)[2], ['v', 'v', 'w'])
})

test('a snapshot of a worklist keeps the sample', async () => {
  const p = await createBlankProject('t')
  const v = (id: string, sample?: string): VariantRec => ({ id, key: 'k', gene: 'SYNG1', hgvsC: 'c.103C>T', ...(sample ? { sample } : {}), source: { file: 'f', line: 1, kind: 'csv' }, evidence: { popmaxAf: null, revel: null, spliceAi: null, functional: [] } })
  const snap = snapshotOf(evaluateProject({ ...p, variants: [v('a', 'S1'), v('b')] }), 'x')
  assert.deepEqual(snap.rows.map(x => x.sample), ['S1', ''])
})
