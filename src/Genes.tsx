import { useMemo } from 'react'
import { fmtAf, type GeneCfg } from './engine/acmg.ts'
import type { Evaluated } from './engine/pipeline.ts'
import { cdsBounds } from './engine/transcript.ts'
import { context, geneCfg } from './projectData.ts'
import { saveProject, type Project } from './store.ts'
import { Link, toast } from './ui.tsx'

const parseRegions = (s: string): [number, number][] => s.split(/[,;]/).map(x => x.trim()).filter(Boolean).map(x => x.split(/[-–]/).map(Number) as [number, number]).filter(([a, b]) => Number.isInteger(a) && Number.isInteger(b) && a <= b)

export default function Genes({ project, evaluated }: { project: Project; evaluated: Evaluated[] }) {
  const genes = useMemo(() => [...new Set(evaluated.map(e => e.an.gene).filter(g => g !== '—'))].sort(), [evaluated])
  const ctx = context(project)
  const set = (g: string, patch: Partial<GeneCfg>) => void saveProject({ ...project, genes: { ...project.genes, [g]: { ...geneCfg(project, g), ...patch } } })
  const num = (g: string, k: 'pm2Max' | 'bs1Min' | 'ba1Min', raw: string) => {
    const n = Number(raw)
    if (raw === '' || Number.isNaN(n) || n < 0 || n > 1) { toast('Frequencies must be between 0 and 1.', { tone: 'error' }); return }
    set(g, { [k]: n })
  }
  return (
    <div className="stack">
      <div className="card pad">
        <h3 className="h5">Gene settings</h3>
        <p className="muted small">These drive PVS1 (LoF mechanism, critical regions), PM2/BS1/BA1 thresholds, and how phase is read (PM3 vs BP2). Use your ClinGen VCEP specification where one exists.</p>
        {genes.length ? (
          <div className="table-wrap">
            <table className="tbl genes">
              <thead><tr><th>Gene</th><th>Inheritance</th><th>LoF mechanism</th><th>PM2 ≤</th><th>BS1 ≥</th><th>BA1 ≥</th><th>Protein (aa)</th><th>Critical regions (aa)</th></tr></thead>
              <tbody>{genes.map(g => {
                const c = geneCfg(project, g)
                return (
                  <tr key={g + JSON.stringify(c)}>
                    <td><b>{g}</b></td>
                    <td><select aria-label={`${g} inheritance`} value={c.inheritance} onChange={e => set(g, { inheritance: e.target.value as GeneCfg['inheritance'] })}><option value="AD">Dominant</option><option value="AR">Recessive</option><option value="XL">X-linked</option></select></td>
                    <td><input type="checkbox" aria-label={`${g} loss of function mechanism`} checked={c.lofMechanism} onChange={e => set(g, { lofMechanism: e.target.checked })} /></td>
                    {(['pm2Max', 'bs1Min', 'ba1Min'] as const).map(k => <td key={k}><input className="short" aria-label={`${g} ${k}`} defaultValue={fmtAf(c[k])} onBlur={e => { if (e.target.value !== fmtAf(c[k])) num(g, k, e.target.value) }} /></td>)}
                    <td><input className="short" type="number" min={1} aria-label={`${g} protein length`} defaultValue={c.proteinLength ?? ''} placeholder="from model" onBlur={e => { const n = Number(e.target.value); set(g, { proteinLength: e.target.value && n > 0 ? n : undefined }) }} /></td>
                    <td><input aria-label={`${g} critical regions`} defaultValue={(c.criticalRegions ?? []).map(r => r.join('-')).join(', ')} placeholder="e.g. 2479-3190" onBlur={e => set(g, { criticalRegions: parseRegions(e.target.value) })} /></td>
                  </tr>
                )
              })}</tbody>
            </table>
          </div>
        ) : <p className="muted small">Import variants to configure their genes.</p>}
      </div>

      <div className="card pad">
        <div className="row between wrap"><h3 className="h5">Transcript models</h3><Link href={`/app/p/${project.id}/import`} className="btn ghost sm">Import GTF / GFF3</Link></div>
        <div className="table-wrap">
          <table className="tbl">
            <thead><tr><th>Transcript</th><th>Gene</th><th>Contig</th><th>Strand</th><th className="num">Exons</th><th className="num">Protein</th><th>Reference</th><th>Source</th><th /></tr></thead>
            <tbody>{ctx.txs.map(t => {
              const imported = project.transcripts.some(x => x.id === t.id)
              return (
                <tr key={t.id}>
                  <td className="mono">{t.id}</td><td>{t.gene}</td><td className="mono">{t.chrom}</td><td>{t.strand === 1 ? '+' : '−'}</td><td className="num">{t.exons.length}</td>
                  <td className="num">{Number.isInteger(cdsBounds(t).aaLength) ? `${cdsBounds(t).aaLength} aa` : '?'}</td>
                  <td>{Object.keys(ctx.genome).some(k => k.replace(/^chr/i, '') === t.chrom.replace(/^chr/i, '')) ? 'loaded' : <span className="muted">not loaded</span>}</td>
                  <td className="small">{t.synthetic ? 'Bundled synthetic' : t.source}</td>
                  <td>{imported && <button type="button" className="btn danger-ghost sm" onClick={() => { const prev = project; void saveProject({ ...project, transcripts: project.transcripts.filter(x => x.id !== t.id) }).then(() => toast(`Removed ${t.id}.`, { undo: () => void saveProject(prev) })) }}>Remove</button>}</td>
                </tr>
              )
            })}</tbody>
          </table>
        </div>
        <p className="muted small">Genes without a model are annotated from HGVS alone: position class and type, but no exon structure, NMD or splice scoring. Add a GTF subset and a FASTA region (see <Link href="/docs#models">Docs</Link>) to enable them.</p>
      </div>
    </div>
  )
}
