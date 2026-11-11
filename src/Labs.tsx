import { useMemo, useState } from 'react'
import { critLabel, STRENGTH_SHORT, type Applied } from './engine/acmg.ts'
import { toCsv } from './engine/csv.ts'
import { reconcile, SEVERITY_LABEL, type Severity } from './engine/discrepancy.ts'
import type { Evaluated } from './engine/pipeline.ts'
import type { Project } from './store.ts'
import { download, Empty, Link, slug } from './ui.tsx'
import { ClassBadge } from './viz.tsx'

const OURS = 'This worklist'

export default function Labs({ project, evaluated }: { project: Project; evaluated: Evaluated[] }) {
  const [sev, setSev] = useState<'all' | Severity>('all')
  const list = useMemo(() => reconcile(project.labRows), [project.labRows])
  const ours = useMemo(() => new Map(evaluated.filter(e => e.an.ann).map(e => [`${e.an.gene}:${e.an.ann!.hgvsC.replace(/^.*:/, '')}`, e])), [evaluated])
  if (!project.labRows.length) return (
    <Empty title="No lab classifications imported">
      <p className="muted">Import a table with one row per lab call: variant, lab, classification and (ideally) the ACMG criteria each lab applied. ClinVar submission exports and internal review spreadsheets both work.</p>
      <div className="row gap center-x"><Link href={`/app/p/${project.id}/import`} className="btn primary">Import lab calls</Link><Link href="/docs#labs" className="btn ghost">Table format</Link></div>
    </Empty>
  )
  const shown = list.filter(d => sev === 'all' || d.severity === sev)
  const counts = (s: Severity) => list.filter(d => d.severity === s).length
  const exportCsv = () => download(`${slug(project.name)}-lab-reconciliation.csv`, toCsv([
    ['variant', 'severity', 'lab', 'classification', 'criteria', 'criteria_points', 'criteria_imply', 'mismatch', 'date', 'summary'],
    ...list.flatMap(d => d.calls.map(c => [d.variant, SEVERITY_LABEL[d.severity], c.lab, c.cls, c.criteria.map(critLabel).join(' '), c.points ?? '', c.implied ?? '', c.mismatch ? 'yes' : '', c.date ?? '', d.summary])),
  ]))
  return (
    <div className="stack">
      <div className="row between wrap gap-sm">
        <div className="seg" role="group" aria-label="Filter by severity">
          <button type="button" aria-pressed={sev === 'all'} onClick={() => setSev('all')}>All {list.length}</button>
          {(['major', 'clinically-significant', 'confidence', 'concordant'] as Severity[]).map(s => <button key={s} type="button" aria-pressed={sev === s} onClick={() => setSev(s)}>{SEVERITY_LABEL[s]} {counts(s)}</button>)}
        </div>
        <button type="button" className="btn ghost sm" onClick={exportCsv}>Reconciliation CSV</button>
      </div>
      {shown.map(d => {
        const mine = ours.get(d.variant)
        const labs = [...d.calls.map(c => c.lab), ...(mine ? [OURS] : [])]
        const mineApplied: Applied[] = mine ? mine.result.counted.map(c => ({ code: c.code, strength: c.strength })) : []
        const codes = [...new Set([...d.diffs.map(x => x.code), ...d.calls.flatMap(c => c.criteria.map(a => a.code)), ...mineApplied.map(a => a.code)])]
        const strengthOf = (lab: string, code: string) => lab === OURS ? mineApplied.find(a => a.code === code)?.strength : d.calls.find(c => c.lab === lab)?.criteria.find(a => a.code === code)?.strength
        return (
          <article key={d.variant} className="card pad disc">
            <div className="row between wrap gap-sm">
              <h3 className="mono h5">{d.variant}</h3>
              <span className={`sev ${d.severity}`}>{SEVERITY_LABEL[d.severity]}</span>
            </div>
            <p className="small">{d.summary}</p>
            <div className="table-wrap">
              <table className="tbl disc-tbl">
                <thead>
                  <tr><th>Criterion</th>{labs.map(l => <th key={l}>{l === OURS ? <Link href={`/app/p/${project.id}/v/${mine!.an.rec.id}`}>{l}</Link> : l}</th>)}</tr>
                  <tr className="cls-row"><th>Call</th>{labs.map(l => {
                    const c = d.calls.find(x => x.lab === l)
                    return <th key={l}>{l === OURS ? <ClassBadge cls={mine!.result.cls} points={mine!.result.points} small /> : <><ClassBadge cls={c!.cls} small />{c!.date && <div className="muted small">{c!.date}</div>}</>}</th>
                  })}</tr>
                </thead>
                <tbody>
                  {codes.map(code => {
                    const diff = d.diffs.some(x => x.code === code)
                    return (
                      <tr key={code} className={diff ? 'diff' : ''}>
                        <td><b className={`crit ${code.startsWith('B') ? 'b' : 'p'}`}>{code}</b></td>
                        {labs.map(l => { const s = strengthOf(l, code); return <td key={l} className="small">{s ? STRENGTH_SHORT[s] : <span className="muted">—</span>}</td> })}
                      </tr>
                    )
                  })}
                  <tr className="sum">
                    <td>Points</td>
                    {labs.map(l => {
                      if (l === OURS) return <td key={l} className="mono">{mine!.result.points}</td>
                      const c = d.calls.find(x => x.lab === l)!
                      return <td key={l} className="mono">{c.points ?? '—'}{c.mismatch && <span className="flag" title="The lab's own criteria total a different class">= {c.implied}</span>}</td>
                    })}
                  </tr>
                </tbody>
              </table>
            </div>
          </article>
        )
      })}
    </div>
  )
}
