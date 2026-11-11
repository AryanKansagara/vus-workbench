import { useEffect, useMemo, useState } from 'react'
import { SAMPLES } from './engine/sample.ts'
import { createBlankProject, createSampleProject, evaluateProject } from './projectData.ts'
import { deleteProject, useProjects, type Project } from './store.ts'
import { Empty, fmtDate, go, Link, Modal, SyntheticBadge, toast } from './ui.tsx'
import { ClassBadge } from './viz.tsx'
import { CLASSES } from './engine/acmg.ts'

export default function Workspace() {
  const list = useProjects()
  const [naming, setNaming] = useState(false)
  useEffect(() => {
    const q = new URLSearchParams(location.search)
    if (q.get('import') === '1') { history.replaceState({}, '', '/app'); setNaming(true) }
    const s = q.get('sample'), tab = q.get('tab')
    // Deep link into a sample variant, e.g. ?sample=synthetic-panel&variant=c.306%2B5G%3EA
    const variant = q.get('variant')
    if (s && SAMPLES.some(x => x.id === s)) void createSampleProject(s).then(p => {
      const v = variant ? evaluateProject(p).find(e => e.an.ann?.hgvsC.endsWith(':' + variant)) : null
      go(v ? `/app/p/${p.id}/v/${v.an.rec.id}` : `/app/p/${p.id}${tab ? '/' + tab : ''}`, { replace: true })
    })
  }, [])

  const create = async (name: string) => {
    const p = await createBlankProject(name || 'New worklist')
    go(`/app/p/${p.id}/import`)
  }
  const sample = async (id: string) => { const p = await createSampleProject(id); go(`/app/p/${p.id}`) }

  if (list === null) return <div className="boot" role="status"><span className="spinner" /> Loading</div>
  return (
    <section className="page">
      <div className="page-head">
        <div>
          <p className="eyebrow">Workspace</p>
          <h1>Worklists</h1>
          <p className="muted">One worklist per case, panel or reanalysis batch. Each keeps its variants, criteria, evidence log, lab calls and dated snapshots.</p>
        </div>
        <button type="button" className="btn primary" onClick={() => setNaming(true)}>New worklist</button>
      </div>

      {list.length === 0 ? (
        <Empty title="No worklists yet" icon={<EmptyArt />}>
          <p className="muted">Open the synthetic proband to see a curated splice variant and two VUS in one gene, or bring your own VCF or HGVS list.</p>
          <div className="row gap center-x wrap">
            <button type="button" className="btn primary" onClick={() => void sample('synthetic-panel')}>Try the sample</button>
            <button type="button" className="btn outline" onClick={() => setNaming(true)}>Import your data</button>
          </div>
        </Empty>
      ) : (
        <div className="proj-grid">{list.map(p => <ProjectCard key={p.id} p={p} />)}</div>
      )}

      <div className="samples">
        <h2 className="h4">Samples</h2>
        <div className="proj-grid">
          {SAMPLES.map(s => (
            <button key={s.id} type="button" className="card sample-card" onClick={() => void sample(s.id)}>
              <strong>{s.name}</strong> <SyntheticBadge />
              <span className="muted small">{s.blurb}</span>
            </button>
          ))}
        </div>
      </div>

      <Modal open={naming} onClose={() => setNaming(false)} title="New worklist">
        <form onSubmit={e => { e.preventDefault(); void create(String(new FormData(e.currentTarget).get('name') ?? '').trim()) }} className="stack">
          <label className="field"><span>Case, panel or batch name</span><input name="name" autoFocus placeholder="Hereditary cancer panel, case 2026-118" /></label>
          <p className="muted small">Next: drop a VCF, an HGVS variant list, a lab classification table or a GTF subset. Files are parsed in this browser.</p>
          <div className="row end gap"><button type="button" className="btn ghost" onClick={() => setNaming(false)}>Cancel</button><button className="btn primary">Create and import</button></div>
        </form>
      </Modal>
    </section>
  )
}

function ProjectCard({ p }: { p: Project }) {
  const counts = useMemo(() => {
    const c: Record<string, number> = {}
    for (const e of evaluateProject(p)) c[e.result.cls] = (c[e.result.cls] ?? 0) + 1
    return c
  }, [p])
  const curated = Object.values(p.curation).filter(c => c.status === 'curated').length
  return (
    <article className="card proj">
      <Link href={`/app/p/${p.id}`} className="proj-link">
        <h3>{p.name} {p.sampleId && <SyntheticBadge />}</h3>
        <p className="muted small">{p.variants.length} variants · {curated} curated · {new Set(p.labRows.map(r => r.lab)).size} labs · updated {fmtDate(p.updatedAt)}</p>
        <div className="row gap-sm wrap">
          {[...CLASSES].reverse().filter(c => counts[c]).map(c => <span key={c} className="row gap-xs"><ClassBadge cls={c} small /> <span className="small">{counts[c]}</span></span>)}
        </div>
      </Link>
      <div className="row gap-sm proj-actions">
        <Link href={`/app/p/${p.id}/labs`} className="btn ghost sm">Lab calls</Link>
        <Link href={`/app/p/${p.id}/history`} className="btn ghost sm">History</Link>
        <button type="button" className="btn danger-ghost sm" onClick={async () => {
          const undo = await deleteProject(p.id)
          toast(`Deleted ${p.name}.`, { undo: () => void undo() })
        }}>Delete</button>
      </div>
    </article>
  )
}

const EmptyArt = () => (
  <svg className="empty-art" viewBox="0 0 160 70" aria-hidden="true">
    <line x1="6" x2="154" y1="35" y2="35" className="ea-intron" />
    <rect x="14" y="22" width="30" height="26" rx="4" className="ea-exon" />
    <rect x="66" y="22" width="26" height="26" rx="4" className="ea-exon" />
    <rect x="116" y="26" width="30" height="18" rx="3" className="ea-utr" />
    <circle cx="100" cy="35" r="6" className="ea-var" />
  </svg>
)
