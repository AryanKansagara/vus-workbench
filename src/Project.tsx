import { lazy, Suspense, useEffect, useMemo, useRef, useState } from 'react'
import { CLASSES, CLASS_SHORT, critLabel, fmtAf, vusTier, type ClassName } from './engine/acmg.ts'
import { phaseOf, type Evaluated } from './engine/pipeline.ts'
import { CONSEQUENCE_LABEL } from './engine/transcript.ts'
import { curatedCsv, geneCfg, snapshotOf, useProject, worklistCsv } from './projectData.ts'
import { getDefaults, saveProject, type Project } from './store.ts'
import { download, Empty, go, isTyping, Link, slug, SyntheticBadge, toast } from './ui.tsx'
import { ClassBadge } from './viz.tsx'

const VariantView = lazy(() => import('./Variant.tsx'))
const Labs = lazy(() => import('./Labs.tsx'))
const ImportTab = lazy(() => import('./Import.tsx'))
const History = lazy(() => import('./History.tsx'))
const Genes = lazy(() => import('./Genes.tsx'))

const TABS = [['worklist', 'Worklist'], ['labs', 'Lab discrepancies'], ['history', 'History'], ['genes', 'Genes & models'], ['import', 'Import']] as const

export default function ProjectView({ id, tab, vid }: { id: string; tab: string; vid?: string }) {
  const { loading, project, evaluated } = useProject(id)
  useEffect(() => { if (project) document.title = `${project.name} · VUS Workbench` }, [project?.name])
  if (loading) return <div className="boot" role="status"><span className="spinner" /> Loading</div>
  if (!project) return (
    <section className="page narrow"><Empty title="Worklist not found"><p className="muted">It may have been deleted, or it lives on another device. Sign in to sync it here.</p><Link href="/app" className="btn primary">All worklists</Link></Empty></section>
  )
  return (
    <section className="page">
      <div className="proj-head">
        <div>
          <p className="eyebrow"><Link href="/app">Worklists</Link> /</p>
          <h1 className="proj-title">
            <input aria-label="Worklist name" defaultValue={project.name} key={project.name} onBlur={e => { const v = e.target.value.trim(); if (v && v !== project.name) void saveProject({ ...project, name: v }) }} onKeyDown={e => { if (e.key === 'Enter') (e.target as HTMLInputElement).blur() }} />
            {project.sampleId && <SyntheticBadge />}
          </h1>
        </div>
        <div className="row gap-sm wrap">
          <button type="button" className="btn ghost sm" onClick={() => { void saveProject({ ...project, snapshots: [...project.snapshots, snapshotOf(evaluated, `Snapshot ${new Date().toLocaleDateString('en-GB')}`)] }); toast('Snapshot saved. Compare it in History.', { tone: 'ok' }) }}>Save snapshot</button>
          <Link href={`/app/p/${project.id}/report`} className="btn ghost sm">Worklist report</Link>
        </div>
      </div>
      {project.partial && <div className="banner warn">This device has the synced copy, which holds curated variants only. Re-import the VCF here to see the full worklist.</div>}
      <nav className="tabs" aria-label="Worklist sections">
        {TABS.map(([k, l]) => <Link key={k} href={`/app/p/${project.id}${k === 'worklist' ? '' : '/' + k}`} aria-current={tab === k || (k === 'worklist' && tab === 'variant') ? 'page' : undefined}>{l}{k === 'labs' && project.labRows.length ? <span className="count">{new Set(project.labRows.map(r => r.variant)).size}</span> : null}</Link>)}
      </nav>
      <Suspense fallback={<div className="boot" role="status"><span className="spinner" /> Loading</div>}>
        {tab === 'worklist' && <Worklist project={project} evaluated={evaluated} />}
        {tab === 'variant' && vid && <VariantView project={project} evaluated={evaluated} vid={vid} />}
        {tab === 'labs' && <Labs project={project} evaluated={evaluated} />}
        {tab === 'import' && <ImportTab project={project} />}
        {tab === 'history' && <History project={project} evaluated={evaluated} />}
        {tab === 'genes' && <Genes project={project} evaluated={evaluated} />}
      </Suspense>
    </section>
  )
}

type Filter = 'all' | ClassName
const STATUS_LABEL = { new: 'New', 'in-review': 'In review', curated: 'Curated' }

function Worklist({ project, evaluated }: { project: Project; evaluated: Evaluated[] }) {
  const [q, setQ] = useState('')
  const [cls, setCls] = useState<Filter>('all')
  const [gene, setGene] = useState('all')
  const [status, setStatus] = useState('all')
  const search = useRef<HTMLInputElement>(null)
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => { if (e.key === '/' && !isTyping(e)) { e.preventDefault(); search.current?.focus() } }
    addEventListener('keydown', onKey)
    return () => removeEventListener('keydown', onKey)
  }, [])
  const genes = useMemo(() => [...new Set(evaluated.map(e => e.an.gene))].sort(), [evaluated])
  const counts = useMemo(() => { const c: Partial<Record<ClassName, number>> = {}; for (const e of evaluated) c[e.result.cls] = (c[e.result.cls] ?? 0) + 1; return c }, [evaluated])
  const rows = evaluated.filter(e =>
    (cls === 'all' || e.result.cls === cls) && (gene === 'all' || e.an.gene === gene) &&
    (status === 'all' || (project.curation[e.an.rec.id]?.status ?? 'new') === status) &&
    (!q || `${e.an.label} ${e.an.rec.key} ${e.an.ann ? CONSEQUENCE_LABEL[e.an.ann.consequence] : ''}`.toLowerCase().includes(q.toLowerCase())))
    .sort((a, b) => CLASSES.indexOf(b.result.cls) - CLASSES.indexOf(a.result.cls) || b.result.points - a.result.points)
  const pairs = usePairs(evaluated).filter(p => p.phase !== 'unknown' || geneCfg(project, p.gene).inheritance === 'AR')

  if (!evaluated.length) return (
    <Empty title="No variants yet">
      <p className="muted">Import a VCF or an HGVS variant list to build the worklist. Lab tables and GTF models can come later.</p>
      <div className="row gap center-x"><Link href={`/app/p/${project.id}/import`} className="btn primary">Import variants</Link><Link href="/docs" className="btn ghost">How to export from your pipeline</Link></div>
    </Empty>
  )
  const name = slug(project.name)
  return (
    <div className="stack">
      <div className="dist card pad">
        <div className="row between wrap gap-sm">
          <strong>{evaluated.length} variants</strong>
          <div className="row gap-sm wrap">
            <button type="button" className="btn ghost sm" onClick={() => download(`${name}-worklist.csv`, worklistCsv(evaluated, project))}>Worklist CSV</button>
            <button type="button" className="btn ghost sm" onClick={() => { const csv = curatedCsv(evaluated, project, getDefaults().institution || 'My lab'); if (csv.split('\n').length <= 2) toast('Nothing is marked curated yet.', { tone: 'error' }); else download(`${name}-curated-classifications.csv`, csv) }}>Curated calls CSV</button>
          </div>
        </div>
        <div className="dist-bar" role="img" aria-label={CLASSES.map(c => `${c}: ${counts[c] ?? 0}`).join(', ')}>
          {[...CLASSES].reverse().map(c => counts[c] ? <button type="button" key={c} className={`seg-${CLASS_SHORT[c].toLowerCase()}`} style={{ flexGrow: counts[c] }} onClick={() => setCls(cls === c ? 'all' : c)} aria-pressed={cls === c} title={`${c}: ${counts[c]}`}>{CLASS_SHORT[c]} {counts[c]}</button> : null)}
        </div>
      </div>

      {pairs.length > 0 && (
        <div className="card pad pairs">
          <strong>Same gene, same sample</strong>
          <p className="muted small">Phase decides PM3 (recessive, in trans with pathogenic) and BP2 (in cis). Phased genotypes from the VCF are used; unphased pairs need parental testing or read-backed phasing.</p>
          <ul>
            {pairs.map(p => <li key={p.key}><b>{p.gene}</b> <span className="muted">{p.sample}</span>: {p.items.map((it, i) => <span key={it.an.rec.id}>{i > 0 && <span className={`phase ${p.phase}`}>{p.phase === 'trans' ? ' in trans with ' : p.phase === 'cis' ? ' in cis with ' : p.phase === 'homozygous' ? ' with ' : ' phase unknown with '}</span>}<Link href={`/app/p/${project.id}/v/${it.an.rec.id}`} className="mono">{it.an.ann?.hgvsC.replace(/^.*:/, '') ?? it.an.rec.key}</Link> <ClassBadge cls={it.result.cls} small /></span>)}</li>)}
          </ul>
        </div>
      )}

      <div className="filters row gap-sm wrap">
        <label className="field grow"><span>Search <kbd>/</kbd></span><input ref={search} type="search" value={q} onChange={e => setQ(e.target.value)} placeholder="Gene, c., p., consequence" /></label>
        <label className="field"><span>Gene</span><select value={gene} onChange={e => setGene(e.target.value)}><option value="all">All genes</option>{genes.map(g => <option key={g}>{g}</option>)}</select></label>
        <label className="field"><span>Class</span><select value={cls} onChange={e => setCls(e.target.value as Filter)}><option value="all">All classes</option>{[...CLASSES].reverse().map(c => <option key={c}>{c}</option>)}</select></label>
        <label className="field"><span>Status</span><select value={status} onChange={e => setStatus(e.target.value)}><option value="all">Any status</option><option value="new">New</option><option value="in-review">In review</option><option value="curated">Curated</option></select></label>
      </div>

      <div className="table-wrap">
        <table className="tbl worklist">
          <thead><tr><th>Variant</th><th>Consequence</th><th className="num">Popmax AF</th><th>Criteria</th><th className="num">Points</th><th>Class</th><th>Status</th></tr></thead>
          <tbody>
            {rows.map(e => {
              const a = e.an.ann, cur = project.curation[e.an.rec.id], href = `/app/p/${project.id}/v/${e.an.rec.id}`
              return (
                <tr key={e.an.rec.id} onClick={ev => { if (!(ev.target as HTMLElement).closest('a')) go(href) }}>
                  <td><Link href={href} className="var-link"><b>{e.an.gene}</b> <span className="mono">{a?.hgvsC.replace(/^.*:/, '') ?? e.an.rec.key}</span>{a?.hgvsP && <span className="mono muted"> {a.hgvsP}</span>}</Link>{e.an.rec.filter && <span className="flag">FILTER {e.an.rec.filter}</span>}{e.an.problems.length > 0 && <span className="flag" title={e.an.problems.join(' ')}>check</span>}</td>
                  <td>{a ? <span className={`csq csq-${a.consequence}`}>{CONSEQUENCE_LABEL[a.consequence]}</span> : <span className="muted small">not annotated</span>}{a?.spliceRegion && <span className="csq csq-splice_region">+ splice region</span>}</td>
                  <td className="num mono">{e.an.rec.evidence.popmaxAf === null ? <span className="muted">—</span> : fmtAf(e.an.rec.evidence.popmaxAf)}</td>
                  <td><div className="crit-list">{e.result.counted.map(c => <span key={c.code} className={`crit ${c.points < 0 ? 'b' : 'p'}`}>{critLabel(c)}</span>)}</div></td>
                  <td className="num mono">{e.result.points > 0 ? '+' : ''}{e.result.points}</td>
                  <td><ClassBadge cls={e.result.cls} points={e.result.points} small />{e.result.cls === 'VUS' && <span className="tier">{vusTier(e.result.points)}</span>}</td>
                  <td><span className={`status ${cur?.status ?? 'new'}`}>{STATUS_LABEL[cur?.status ?? 'new']}</span></td>
                </tr>
              )
            })}
          </tbody>
        </table>
        {!rows.length && <p className="muted center pad">No variants match these filters.</p>}
      </div>
    </div>
  )
}

function usePairs(evaluated: Evaluated[]) {
  return useMemo(() => {
    const m = new Map<string, Evaluated[]>()
    for (const e of evaluated) if (e.an.ann && e.an.rec.gt && !/^1[|/]1$/.test(e.an.rec.gt)) { const k = `${e.an.gene}|${e.an.rec.sample ?? ''}`; m.set(k, [...(m.get(k) ?? []), e]) }
    const out: { key: string; gene: string; sample: string; items: Evaluated[]; phase: string }[] = []
    for (const [k, items] of m) if (items.length >= 2) {
      const [gene, sample] = k.split('|')
      // List pathogenic first so the phase relation reads from it.
      const sorted = [...items].sort((a, b) => b.result.points - a.result.points)
      for (let i = 1; i < sorted.length; i++) out.push({ key: `${k}|${i}`, gene, sample, items: [sorted[0], sorted[i]], phase: phaseOf(sorted[0].an.rec, sorted[i].an.rec) })
    }
    return out
  }, [evaluated])
}
