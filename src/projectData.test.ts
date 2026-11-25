import { test } from 'node:test'
import assert from 'node:assert/strict'
import { createBlankProject, curate } from './projectData.ts'
import { getProject, saveProject, updateProject } from './store.ts'
import type { Strength } from './engine/acmg.ts'
import type { VariantRec } from './engine/io.ts'

const rec = (id: string): VariantRec => ({ id, key: id, gene: 'SYNG1', hgvsC: 'c.103C>T', source: { file: 'f', line: 1, kind: 'csv' }, evidence: { popmaxAf: null, revel: null, spliceAi: null, functional: [] } })
const setup = async () => { const p = await createBlankProject('t'); await saveProject({ ...p, variants: [rec('a'), rec('b')] }); return p.id }
const crit = (code: string, strength: Strength = 'Supporting') => ({ code, met: true, strength, evidence: '', source: 'curator' as const, editor: 'x', at: 1 })
const add = (id: string, code: string, vid = 'a') => curate(id, vid, 'Criterion applied', code, ({ cur }) => { cur.overrides[code] = crit(code) })

test('actions started in the same tick are all kept and logged', async () => {
  const id = await setup()
  await Promise.all([
    curate(id, 'a', 'Evidence updated', 'AF', ({ rec }) => { rec.evidence.popmaxAf = 0.2 }),
    add(id, 'PP4'), add(id, 'PP1', 'b'),
    updateProject(id, p => ({ ...p, name: 'renamed' })),
  ])
  const p = getProject(id)!
  assert.equal(p.variants[0].evidence.popmaxAf, 0.2)
  assert.deepEqual(p.curation.a.audit.map(x => x.detail), ['AF', 'PP4'])
  assert.ok(p.curation.b.overrides.PP1)
  assert.equal(p.name, 'renamed')
})

test('undo with nothing after it leaves no trace on the variant and keeps the rest of the worklist', async () => {
  const id = await setup()
  const undo = await add(id, 'PP4')
  await add(id, 'PP1', 'b')
  await updateProject(id, p => ({ ...p, name: 'renamed' }))
  await undo()
  const p = getProject(id)!
  assert.equal(p.curation.a, undefined)
  assert.ok(p.curation.b.overrides.PP1)
  assert.equal(p.name, 'renamed')
})

test('undo of one change keeps later changes to the same variant and is logged', async () => {
  const id = await setup()
  const undo = await add(id, 'PP4')
  await curate(id, 'a', 'Evidence updated', 'AF', ({ rec }) => { rec.evidence.popmaxAf = 0.2 })
  await curate(id, 'a', 'Note', 'n', ({ cur }) => { cur.note = 'n' })
  await undo()
  const p = getProject(id)!
  assert.equal(p.curation.a.overrides.PP4, undefined)
  assert.equal(p.curation.a.note, 'n')
  assert.equal(p.variants[0].evidence.popmaxAf, 0.2)
  assert.deepEqual(p.curation.a.audit.map(x => x.action), ['Criterion applied', 'Evidence updated', 'Note', 'Undone'])
})

test('undo leaves a criterion that was changed again since', async () => {
  const id = await setup()
  const undo = await add(id, 'PP4')
  await curate(id, 'a', 'Strength changed', 'PP4', ({ cur }) => { cur.overrides.PP4 = crit('PP4', 'Moderate') })
  await undo()
  assert.equal(getProject(id)!.curation.a.overrides.PP4.strength, 'Moderate')
  assert.equal(getProject(id)!.curation.a.audit.length, 2)
})

test('undo of a removed functional result puts back that result only', async () => {
  const id = await setup()
  const f = (fid: string) => ({ id: fid, assay: fid, result: 'abnormal' as const, oddsPath: null, note: '', editor: 'x', at: 1 })
  await curate(id, 'a', 'Functional evidence added', 'f1', ({ rec }) => { rec.evidence.functional.push(f('f1')) })
  const undo = await curate(id, 'a', 'Functional evidence removed', 'f1', ({ rec }) => { rec.evidence.functional = [] })
  await curate(id, 'a', 'Functional evidence added', 'f2', ({ rec }) => { rec.evidence.functional.push(f('f2')) })
  await undo()
  assert.deepEqual(getProject(id)!.variants[0].evidence.functional.map(x => x.id), ['f2', 'f1'])
})

test('an action on a deleted worklist or variant does not bring it back', async () => {
  const id = await setup()
  await curate('gone', 'a', 'x', 'y', () => {})
  await curate(id, 'gone', 'x', 'y', () => {})
  assert.equal(await updateProject('gone', p => p), null)
  assert.equal(getProject('gone'), undefined)
})
