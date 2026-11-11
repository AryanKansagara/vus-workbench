import { useMemo, useState } from 'react'
import type { Evaluated } from './engine/pipeline.ts'
import { compareSnapshots, snapshotOf } from './projectData.ts'
import { saveProject, type Project } from './store.ts'
import { deleteSnapshotCloud } from './sync.ts'
import { Empty, fmtDateTime, Link, toast } from './ui.tsx'
import { ClassBadge } from './viz.tsx'

const CURRENT = 'current'

export default function History({ project, evaluated }: { project: Project; evaluated: Evaluated[] }) {
  const snaps = [...project.snapshots].sort((a, b) => a.createdAt - b.createdAt)
  const [a, setA] = useState(snaps.length >= 2 ? snaps[snaps.length - 2].id : snaps[0]?.id ?? '')
  const [b, setB] = useState(snaps.length >= 2 ? snaps[snaps.length - 1].id : CURRENT)
  const current = useMemo(() => snapshotOf(evaluated, 'Current state'), [evaluated])
  const pick = (id: string) => id === CURRENT ? current : snaps.find(s => s.id === id)
  const A = pick(a), B = pick(b)
  const diff = A && B ? compareSnapshots(A, B) : null
  const events = useMemo(() => project.variants.flatMap(v => (project.curation[v.id]?.history ?? []).filter(h => h.from).map(h => ({ ...h, vid: v.id, label: evaluated.find(e => e.an.rec.id === v.id)?.an.label ?? v.key }))).sort((x, y) => y.at - x.at), [project, evaluated])

  if (!snaps.length) return <Empty title="No snapshots yet"><p className="muted">A snapshot is saved automatically on every import, or use Save snapshot above. Compare any two to see what was reclassified.</p></Empty>
  return (
    <div className="stack">
      <div className="card pad">
        <div className="row gap wrap end-items">
          <label className="field"><span>From</span><select value={a} onChange={e => setA(e.target.value)}>{snaps.map(s => <option key={s.id} value={s.id}>{s.label} · {fmtDateTime(s.createdAt)}</option>)}</select></label>
          <span aria-hidden="true" className="arrow-sep">→</span>
          <label className="field"><span>To</span><select value={b} onChange={e => setB(e.target.value)}>{snaps.map(s => <option key={s.id} value={s.id}>{s.label} · {fmtDateTime(s.createdAt)}</option>)}<option value={CURRENT}>Current state</option></select></label>
        </div>
        {diff && (
          <div className="stack-sm mt">
            <p className="small"><b>{diff.reclassified.filter(r => r.from !== r.to).length}</b> reclassified · <b>{diff.reclassified.filter(r => r.from === r.to).length}</b> criteria changed, same class · <b>{diff.added.length}</b> new · <b>{diff.removed.length}</b> removed · <b>{diff.unchanged}</b> unchanged</p>
            {diff.reclassified.length > 0 && (
              <div className="table-wrap">
                <table className="tbl">
                  <thead><tr><th>Variant</th><th>Class</th><th className="num">Δ points</th><th>Criteria added</th><th>Criteria removed</th></tr></thead>
                  <tbody>{diff.reclassified.map(r => (
                    <tr key={r.key}>
                      <td className="mono small">{r.label}</td>
                      <td>{r.from === r.to ? <ClassBadge cls={r.to} small /> : <><ClassBadge cls={r.from} small /> → <ClassBadge cls={r.to} small /></>}</td>
                      <td className="num mono">{r.dPoints > 0 ? '+' : ''}{r.dPoints}</td>
                      <td className="small">{r.added.join(' ') || '—'}</td>
                      <td className="small">{r.removed.join(' ') || '—'}</td>
                    </tr>
                  ))}</tbody>
                </table>
              </div>
            )}
            {diff.added.length > 0 && <p className="small">New: {diff.added.map(r => r.label).join(', ')}</p>}
            {diff.removed.length > 0 && <p className="small">Removed: {diff.removed.map(r => r.label).join(', ')}</p>}
          </div>
        )}
      </div>

      <div className="two-col">
        <div className="card pad">
          <h3 className="h5">Reclassifications</h3>
          {events.length ? <ol className="timeline">{events.map((h, i) => <li key={i}><span className="muted small">{fmtDateTime(h.at)} · {h.editor}</span><div><Link href={`/app/p/${project.id}/v/${h.vid}`} className="mono small">{h.label}</Link></div><div><ClassBadge cls={h.from!} small /> → <ClassBadge cls={h.to} small /></div><span className="small muted">{h.reason}</span></li>)}</ol> : <p className="muted small">No variant has changed class through curation yet.</p>}
        </div>
        <div className="card pad">
          <h3 className="h5">Snapshots</h3>
          <ul className="snap-list">
            {[...snaps].reverse().map(s => (
              <li key={s.id}>
                <div><b>{s.label}</b><div className="muted small">{fmtDateTime(s.createdAt)} · {s.rows.length} variants</div></div>
                <button type="button" className="btn danger-ghost sm" onClick={() => {
                  const prev = project
                  void saveProject({ ...project, snapshots: project.snapshots.filter(x => x.id !== s.id) }).then(() => {
                    toast(`Deleted snapshot ${s.label}.`, { undo: () => void saveProject(prev), ttl: 8000 })
                    // Undo re-saves the snapshot locally and the next sync uploads it again.
                    void deleteSnapshotCloud(s.id)
                  })
                }}>Delete</button>
              </li>
            ))}
          </ul>
        </div>
      </div>
    </div>
  )
}
