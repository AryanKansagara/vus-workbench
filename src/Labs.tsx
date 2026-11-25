import { useMemo, useState } from 'react'
import { clinvarLabRows, fetchClinVar } from './engine/clinvar.ts'
import { updateProject, uuid } from './store.ts'
import { toast } from './ui.tsx'
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
  // Lab rows name a variant by gene and c. description, or by genomic key when they came from a VCF.
  const ours = useMemo(() => new Map(evaluated.flatMap(e => [[e.an.rec.key, e] as const, ...(e.an.ann ? [[`${e.an.gene}:${e.an.ann.hgvsC.replace(/^.*:/, '')}`, e] as const] : [])])), [evaluated])
  const [busy, setBusy] = useState('')
  const lookups = useMemo(() => evaluated.filter(e => e.an.ann && !e.an.tx?.synthetic && !/^SYN/.test(e.an.gene) && e.an.gene !== '—'), [evaluated])
  const fetchAll = async () => {
    const sig = (r: { variant: string; lab: string; date?: string; cls: string }) => `${r.variant}|${r.lab}|${r.date ?? ''}|${r.cls}`
    const have = new Set(project.labRows.map(sig))
    const fresh: Project['labRows'] = []
    let found = 0, none = 0, failed = ''
    const todo = lookups.slice(0, 25)
    for (const [i, e] of todo.entries()) {
      const c = e.an.ann!.hgvsC.replace(/^.*:/, '')
      setBusy(`${i + 1} of ${todo.length}: ${e.an.gene} ${c}`)
      try {
        const rec = await fetchClinVar({ transcript: e.an.tx?.id ?? e.an.rec.transcript, gene: e.an.gene, hgvsC: c })
        if (!rec) { none++; continue }
        found++
        for (const r of clinvarLabRows(rec, `${e.an.gene}:${c}`, uuid).rows) if (!have.has(sig(r))) { have.add(sig(r)); fresh.push(r) }
      } catch (err) { failed = err instanceof Error ? err.message : 'lookup failed'; break }
    }
    setBusy('')
    if (fresh.length) await updateProject(project.id, p => ({ ...p, labRows: [...p.labRows, ...fresh], imports: [...p.imports, { id: uuid(), kind: 'lab-csv', fileName: `ClinVar submissions (${found} variants)`, rows: fresh.length, issues: 0, at: Date.now() }] }))
    if (failed) toast(`${failed}${fresh.length ? ` ${fresh.length} submissions were added before that.` : ''}`, { tone: 'error' })
    else toast(`${fresh.length} new ClinVar submissions for ${found} variant${found === 1 ? '' : 's'}${none ? `; ${none} not in ClinVar under that description` : ''}${lookups.length > todo.length ? `. First ${todo.length} variants only` : ''}.`, { tone: fresh.length ? 'ok' : undefined })
  }
  const clinvarBtn = lookups.length > 0 && (
    <button type="button" className="btn ghost sm" disabled={!!busy} onClick={() => void fetchAll()} title="Looks each worklist variant up in ClinVar by its HGVS description (NCBI E-utilities). Only the description is sent.">{busy ? `ClinVar ${busy}…` : `Fetch ClinVar submissions (${Math.min(25, lookups.length)})`}</button>
  )
  if (!project.labRows.length) return (
    <Empty title="No lab classifications imported">
      <p className="muted">Import a table with one row per lab call: variant, lab, classification and (ideally) the ACMG criteria each lab applied. ClinVar submission exports and internal review spreadsheets both work.</p>
      <div className="row gap center-x wrap"><Link href={`/app/p/${project.id}/import`} className="btn primary">Import lab calls</Link>{clinvarBtn}<Link href="/docs#labs" className="btn ghost">Table format</Link></div>
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
        <div className="row gap-sm wrap">{clinvarBtn}<button type="button" className="btn ghost sm" onClick={exportCsv}>Reconciliation CSV</button></div>
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
              <h3 className="mono h5">{mine && !d.variant.includes('c.') ? mine.an.label === d.variant ? `${mine.an.gene} ${d.variant}` : `${mine.an.label} (${d.variant})` : d.variant}</h3>
              <span className={`sev ${d.severity}`}>{SEVERITY_LABEL[d.severity]}</span>
            </div>
            <p className="small">{d.summary}</p>
            {labs.length > 8 ? (
              <div className="table-wrap">
                <table className="tbl small">
                  <thead><tr><th>Lab</th><th>Call</th><th>Evaluated</th><th>Criteria</th><th className="num">Points</th></tr></thead>
                  <tbody>
                    {mine && <tr><td><Link href={`/app/p/${project.id}/v/${mine.an.rec.id}`}>{OURS}</Link></td><td><ClassBadge cls={mine.result.cls} points={mine.result.points} small /></td><td /><td className="mono">{mineApplied.map(critLabel).join(' ')}</td><td className="num mono">{mine.result.points}</td></tr>}
                    {d.calls.map(c => <tr key={c.lab}><td>{c.lab}</td><td><ClassBadge cls={c.cls} small /></td><td className="muted">{c.date ?? ''}</td><td className="mono">{c.criteria.map(critLabel).join(' ') || <span className="muted">not stated</span>}</td><td className="num mono">{c.points ?? '—'}{c.mismatch && <span className="flag" title="The lab's own criteria total a different class">= {c.implied}</span>}</td></tr>)}
                  </tbody>
                </table>
              </div>
            ) : (
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
            )}
          </article>
        )
      })}
    </div>
  )
}
