import { useEffect } from 'react'
import { critLabel, fmtAf, STRENGTH_SHORT, vusTier } from './engine/acmg.ts'
import type { Evaluated } from './engine/pipeline.ts'
import { CONSEQUENCE_LABEL } from './engine/transcript.ts'
import { useProject } from './projectData.ts'
import { getDefaults, type Project } from './store.ts'
import { fmtDate, fmtDateTime, Link, Monogram } from './ui.tsx'
import { ClassBadge, Meter, Track } from './viz.tsx'

export default function Report({ id, vid }: { id: string; vid?: string }) {
  const { loading, project, evaluated } = useProject(id)
  useEffect(() => { if (project) document.title = `${vid ? 'Curation report' : 'Worklist report'} · ${project.name}` }, [project?.name, vid])
  if (loading) return <div className="boot" role="status"><span className="spinner" /> Loading</div>
  if (!project) return <div className="page narrow"><h1>Report not found</h1><Link href="/app">Back to worklists</Link></div>
  const e = vid ? evaluated.find(x => x.an.rec.id === vid) : null
  const d = getDefaults()
  return (
    <div className="report-page">
      <div className="report-bar no-print">
        <Link href={vid ? `/app/p/${id}/v/${vid}` : `/app/p/${id}`} className="btn ghost sm">Back</Link>
        <button type="button" className="btn primary sm" onClick={() => print()}>Print or save as PDF</button>
      </div>
      <article className="sheet">
        <header className="sheet-head">
          <div className="row gap-sm"><Monogram /><div><b>VUS Workbench</b><div className="small muted">{d.institution || 'Institution not set'} · {project.name}</div></div></div>
          <div className="small right">Generated {fmtDateTime(Date.now())}<br />Prepared by {d.editor || 'curator not set (Settings)'}</div>
        </header>
        <p className="ruo">Research use only. Not a diagnostic report. Classifications are curator-reviewed evidence summaries and must be confirmed under your laboratory's validated procedures before any clinical use.{project.sampleId ? ' This worklist is a synthetic sample: no patient data.' : ''}</p>
        {e ? <VariantSheet e={e} project={project} /> : <WorklistSheet list={evaluated} project={project} />}
        <footer className="sheet-foot small muted">Criteria per Richards et al. 2015 with ClinGen SVI refinements; points per Tavtigian et al. 2020 (P ≥ 10, LP 6 to 9, VUS 0 to 5, LB −1 to −6, B ≤ −7).</footer>
      </article>
    </div>
  )
}

function VariantSheet({ e, project }: { e: Evaluated; project: Project }) {
  const { an, result } = e
  const a = an.ann
  const cur = project.curation[an.rec.id]
  const ev = an.rec.evidence
  return (
    <>
      <h1 className="sheet-title">{an.gene} <span className="mono">{a?.hgvsC ?? an.rec.key}</span> {a?.hgvsP && <span className="mono">{a.hgvsP}</span>}</h1>
      <div className="sheet-grid">
        <dl className="facts">
          <div><dt>Classification</dt><dd><ClassBadge cls={result.cls} points={result.points} /> {result.points} points{result.cls === 'VUS' ? ` (${vusTier(result.points)})` : ''}</dd></div>
          {a && <div><dt>Consequence</dt><dd>{CONSEQUENCE_LABEL[a.consequence]}{a.spliceRegion ? ', splice region' : ''}{a.exon ? `, exon ${a.exon}` : a.intron ? `, intron ${a.intron}` : ''}</dd></div>}
          {an.g && <div><dt>Genomic</dt><dd className="mono">{an.g.chrom}:{an.g.pos} {an.g.ref}&gt;{an.g.alt}</dd></div>}
          {an.rec.gt && <div><dt>Genotype</dt><dd>{an.rec.gt}{an.rec.sample ? ` (${an.rec.sample})` : ''}</dd></div>}
          <div><dt>Popmax AF</dt><dd>{ev.popmaxAf === null ? 'not entered' : fmtAf(ev.popmaxAf)}</dd></div>
          <div><dt>REVEL / SpliceAI</dt><dd>{ev.revel ?? '—'} / {ev.spliceAi ?? '—'}</dd></div>
          <div><dt>Status</dt><dd>{cur?.status ?? 'new'}</dd></div>
        </dl>
        <Meter points={result.points} cls={result.cls} />
      </div>
      <Track tx={an.tx} ann={a} gpos={an.g?.pos} height={110} />
      {an.splice && <p className="small">Splice check: {an.splice.summary}</p>}
      <h2 className="h5">Criteria applied</h2>
      <table className="tbl small">
        <thead><tr><th>Criterion</th><th>Strength</th><th className="num">Points</th><th>Evidence</th><th>Source</th></tr></thead>
        <tbody>
          {e.criteria.filter(c => c.met).map(c => {
            const counted = result.counted.find(x => x.code === c.code)
            return <tr key={c.code}><td><b>{c.code}</b></td><td>{STRENGTH_SHORT[c.strength]}</td><td className="num">{counted ? counted.points : 'not counted'}</td><td>{c.evidence}</td><td>{c.source === 'engine' ? 'Engine rule' : `${c.editor}, ${fmtDate(c.at)}`}</td></tr>
          })}
        </tbody>
      </table>
      {result.excluded.map(x => <p key={x.code} className="small">{x.code}: {x.why}</p>)}
      {ev.functional.length > 0 && (<><h2 className="h5">Functional evidence</h2><ul className="small">{ev.functional.map(f => <li key={f.id}>{f.result === 'abnormal' ? 'Damaging' : 'Normal'}: {f.assay}{f.oddsPath !== null ? `, OddsPath ${f.oddsPath}` : ''}. {f.note} ({f.editor}, {fmtDate(f.at)})</li>)}</ul></>)}
      {cur?.note && <><h2 className="h5">Curator note</h2><p className="small">{cur.note}</p></>}
      {cur?.history.length ? (<><h2 className="h5">Classification history</h2><ul className="small">{cur.history.map((h, i) => <li key={i}>{fmtDate(h.at)}: {h.from ? `${h.from} → ` : ''}{h.to} ({h.points} points). {h.reason} [{h.editor}]</li>)}</ul></>) : null}
      {cur?.audit.length ? (<><h2 className="h5">Evidence log</h2><ul className="small audit-print">{cur.audit.map((x, i) => <li key={i}>{fmtDateTime(x.at)} · {x.editor} · {x.action}: {x.detail}</li>)}</ul></>) : null}
      <div className="sign"><div>Curated by</div><div>Reviewed by</div><div>Date</div></div>
    </>
  )
}

function WorklistSheet({ list, project }: { list: Evaluated[]; project: Project }) {
  return (
    <>
      <h1 className="sheet-title">{project.name}</h1>
      <p className="small">{list.length} variants · {Object.values(project.curation).filter(c => c.status === 'curated').length} curated</p>
      <table className="tbl small">
        <thead><tr><th>Variant</th><th>Consequence</th><th>Criteria</th><th className="num">Points</th><th>Class</th><th>Status</th></tr></thead>
        <tbody>{list.map(e => <tr key={e.an.rec.id}><td><b>{e.an.gene}</b> <span className="mono">{e.an.ann?.hgvsC.replace(/^.*:/, '') ?? e.an.rec.key} {e.an.ann?.hgvsP ?? ''}</span></td><td>{e.an.ann ? CONSEQUENCE_LABEL[e.an.ann.consequence] : '—'}</td><td>{e.result.counted.map(critLabel).join(' ')}</td><td className="num">{e.result.points}</td><td><ClassBadge cls={e.result.cls} small /></td><td>{project.curation[e.an.rec.id]?.status ?? 'new'}</td></tr>)}</tbody>
      </table>
    </>
  )
}
