import { useMemo, useRef, useState } from 'react'
import { autoMap, detectKind, FILE_KIND_LABEL, LAB_FIELDS, parseCsv, parseFasta, parseGtf, parseVcf, validateLabCsv, validateVariantCsv, VARIANT_FIELDS, type FieldSpec, type FileKind, type Issue, type LabRow, type Mapping, type VariantRec } from './engine/io.ts'
import type { Genome, Transcript } from './engine/transcript.ts'
import { LAB_CSV, REAL_VARIANTS_CSV, SAMPLE_GTF, sampleVcf } from './engine/sample.ts'
import { context, evaluateProject, snapshotOf } from './projectData.ts'
import { saveProject, uuid, type Project } from './store.ts'
import { download, go, Link, toast } from './ui.tsx'

interface Pending { id: string; name: string; kind: FileKind; text: string; map: Mapping }

const MAX = 50_000_000

export default function ImportTab({ project }: { project: Project }) {
  const [pending, setPending] = useState<Pending[]>([])
  const [drag, setDrag] = useState(false)
  const [paste, setPaste] = useState('')
  const input = useRef<HTMLInputElement>(null)

  const add = (name: string, text: string) => {
    const kind = detectKind(text, name)
    if (kind === 'unknown') { toast(`${name}: could not tell the format. Expected VCF, CSV/TSV with an HGVS or lab column, GTF/GFF3 or FASTA.`, { tone: 'error' }); return }
    const csv = kind === 'variant-csv' || kind === 'lab-csv' ? parseCsv(text) : null
    setPending(p => [...p, { id: uuid(), name, kind, text, map: csv ? autoMap(csv.rows[0] ?? [], kind === 'lab-csv' ? LAB_FIELDS : VARIANT_FIELDS) : {} }])
  }
  const readFiles = async (files: FileList | File[]) => {
    for (const f of [...files]) {
      if (f.size > MAX) { toast(`${f.name} is over 50 MB. Filter the VCF to your panel regions first (see Docs).`, { tone: 'error' }); continue }
      if (/\.gz$/i.test(f.name)) { toast(`${f.name} is compressed. Run bcftools view -Ov (or gunzip) first.`, { tone: 'error' }); continue }
      add(f.name, await f.text())
    }
  }

  return (
    <div className="stack">
      <div
        className={`drop${drag ? ' over' : ''}`}
        onDragOver={e => { e.preventDefault(); setDrag(true) }}
        onDragLeave={() => setDrag(false)}
        onDrop={e => { e.preventDefault(); setDrag(false); void readFiles(e.dataTransfer.files) }}
      >
        <div className="drop-ring" aria-hidden="true" />
        <h2 className="h4">Drop files here</h2>
        <p className="muted">VCF 4.x · CSV/TSV variant list with HGVS c. · lab classification table · GTF/GFF3 transcript subset · FASTA region. Several at once is fine.</p>
        <div className="row gap center-x wrap">
          <button type="button" className="btn primary" onClick={() => input.current?.click()}>Choose files</button>
          <Link href="/docs" className="btn ghost">How to export them</Link>
        </div>
        <input ref={input} type="file" multiple hidden accept=".vcf,.csv,.tsv,.txt,.gtf,.gff,.gff3,.fa,.fasta" onChange={e => { if (e.target.files) void readFiles(e.target.files); e.target.value = '' }} />
      </div>

      <details className="card pad">
        <summary><b>Paste instead</b> <span className="muted small">(a few HGVS lines, VCF records or a lab table)</span></summary>
        <textarea rows={6} value={paste} onChange={e => setPaste(e.target.value)} placeholder={'gene,hgvs_c,hgvs_p\nBRCA2,c.68-7T>A,'} aria-label="Pasted data" />
        <button type="button" className="btn ghost sm" disabled={!paste.trim()} onClick={() => { add('pasted text', paste); setPaste('') }}>Read pasted text</button>
      </details>

      <div className="row gap-sm wrap small">
        <span className="muted">Synthetic example files:</span>
        <button type="button" className="linkish" onClick={() => add('SYN-PROBAND-01.synthetic.vcf', sampleVcf())}>VCF</button>
        <button type="button" className="linkish" onClick={() => add('public-hgvs.synthetic-evidence.csv', REAL_VARIANTS_CSV)}>HGVS list</button>
        <button type="button" className="linkish" onClick={() => add('lab-calls.synthetic.csv', LAB_CSV)}>Lab table</button>
        <button type="button" className="linkish" onClick={() => add('synthetic-genes.gtf', SAMPLE_GTF)}>GTF</button>
        <span className="muted">·</span>
        <button type="button" className="linkish" onClick={() => download('SYN-PROBAND-01.synthetic.vcf', sampleVcf(), 'text/plain')}>download VCF</button>
      </div>

      {pending.map(p => <PendingFile key={p.id} p={p} project={project} onMap={map => setPending(x => x.map(y => y.id === p.id ? { ...y, map } : y))} onDone={() => setPending(x => x.filter(y => y.id !== p.id))} />)}

      {project.imports.length > 0 && (
        <div className="card pad">
          <h3 className="h5">Import history</h3>
          <table className="tbl small">
            <thead><tr><th>File</th><th>Type</th><th className="num">Records</th><th className="num">Issues</th><th>When</th></tr></thead>
            <tbody>{[...project.imports].reverse().map(i => <tr key={i.id}><td>{i.fileName}</td><td>{FILE_KIND_LABEL[i.kind as FileKind] ?? i.kind}</td><td className="num">{i.rows}</td><td className="num">{i.issues}</td><td>{new Date(i.at).toLocaleString('en-GB')}</td></tr>)}</tbody>
          </table>
        </div>
      )}
    </div>
  )
}

interface Parsed { variants?: VariantRec[]; labRows?: LabRow[]; transcripts?: Transcript[]; genome?: Genome; issues: Issue[]; total: number; headers?: string[] }

function PendingFile({ p, project, onMap, onDone }: { p: Pending; project: Project; onMap(m: Mapping): void; onDone(): void }) {
  const [showAll, setShowAll] = useState(false)
  const parsed: Parsed = useMemo(() => {
    const g = context(project).genome
    if (p.kind === 'vcf') { const r = parseVcf(p.text, p.name, g, uuid); return { variants: r.variants, issues: r.issues, total: r.records } }
    if (p.kind === 'gtf' || p.kind === 'gff3') { const r = parseGtf(p.text, p.name); return { transcripts: r.transcripts, issues: r.issues, total: r.transcripts.length } }
    if (p.kind === 'fasta') { const r = parseFasta(p.text); return { genome: r.genome, issues: r.issues, total: Object.keys(r.genome).length } }
    const csv = parseCsv(p.text)
    if (p.kind === 'lab-csv') { const r = validateLabCsv(csv, p.map, p.name, uuid); return { labRows: r.rows, issues: r.issues, total: r.total, headers: csv.rows[0] } }
    const r = validateVariantCsv(csv, p.map, p.name, g, uuid)
    return { variants: r.variants, issues: r.issues, total: r.total, headers: csv.rows[0] }
  }, [p, project.genome, project.transcripts])
  const errors = parsed.issues.filter(i => i.level === 'error').length
  const ok = (parsed.variants?.length ?? parsed.labRows?.length ?? parsed.transcripts?.length ?? Object.keys(parsed.genome ?? {}).length) || 0
  const fields: FieldSpec[] | null = p.kind === 'lab-csv' ? LAB_FIELDS : p.kind === 'variant-csv' ? VARIANT_FIELDS : null

  const commit = async () => {
    let next: Project = { ...project }
    let msg = ''
    if (parsed.variants) {
      const have = new Set(project.variants.map(v => v.key + (v.sample ?? '')))
      const fresh = parsed.variants.filter(v => !have.has(v.key + (v.sample ?? '')))
      next.variants = [...project.variants, ...fresh]
      msg = `Added ${fresh.length} variants${fresh.length < parsed.variants.length ? `, ${parsed.variants.length - fresh.length} already in the worklist` : ''}`
    }
    if (parsed.labRows) {
      const sig = (r: LabRow) => `${r.variant}|${r.lab}|${r.date ?? ''}|${r.cls}`
      const have = new Set(project.labRows.map(sig))
      const fresh = parsed.labRows.filter(r => !have.has(sig(r)))
      next.labRows = [...project.labRows, ...fresh]
      msg = `Added ${fresh.length} lab calls`
    }
    if (parsed.transcripts) {
      next.transcripts = [...project.transcripts.filter(t => !parsed.transcripts!.some(n => n.id === t.id)), ...parsed.transcripts]
      msg = `Loaded ${parsed.transcripts.length} transcript models`
    }
    if (parsed.genome) {
      next.genome = { ...project.genome, ...parsed.genome }
      msg = `Loaded reference for ${Object.keys(parsed.genome).join(', ')}`
    }
    next.imports = [...project.imports, { id: uuid(), kind: p.kind, fileName: p.name, rows: ok, issues: parsed.issues.length, at: Date.now() }]
    if (parsed.variants || parsed.transcripts || parsed.genome) next.snapshots = [...project.snapshots, snapshotOf(evaluateProject(next), `After importing ${p.name}`)]
    await saveProject(next)
    onDone()
    toast(msg, { tone: 'ok' })
    if (parsed.variants) go(`/app/p/${project.id}`)
    if (parsed.labRows) go(`/app/p/${project.id}/labs`)
  }

  const shownIssues = showAll ? parsed.issues : parsed.issues.slice(0, 12)
  return (
    <div className="card pad pending">
      <div className="row between wrap gap-sm">
        <div><strong>{p.name}</strong> <span className="tag">{FILE_KIND_LABEL[p.kind]}</span></div>
        <div className="row gap-sm"><button type="button" className="btn ghost sm" onClick={onDone}>Discard</button><button type="button" className="btn primary sm" disabled={!ok} onClick={() => void commit()}>Import {ok} {p.kind === 'lab-csv' ? 'calls' : p.kind === 'gtf' || p.kind === 'gff3' ? 'models' : p.kind === 'fasta' ? 'regions' : 'variants'}</button></div>
      </div>
      {fields && parsed.headers && (
        <div className="mapping">
          {fields.map(f => (
            <label key={f.key} className="field">
              <span>{f.label}{f.required && ' *'}</span>
              <select value={p.map[f.key] ?? ''} onChange={e => { const m = { ...p.map }; if (e.target.value === '') delete m[f.key]; else m[f.key] = Number(e.target.value); onMap(m) }}>
                <option value="">Not mapped</option>
                {parsed.headers!.map((h, i) => <option key={i} value={i}>{h}</option>)}
              </select>
            </label>
          ))}
        </div>
      )}
      <div className="report">
        <div className="row gap wrap small">
          <span><b>{ok}</b> ready</span>
          <span className={errors ? 'error' : 'muted'}><b>{errors}</b> errors (skipped)</span>
          <span className="muted"><b>{parsed.issues.filter(i => i.level === 'warning').length}</b> warnings</span>
          <span className="muted"><b>{parsed.issues.filter(i => i.level === 'info').length}</b> notes</span>
          <span className="muted">{parsed.total} {p.kind === 'vcf' ? 'records' : 'rows'} read</span>
        </div>
        {parsed.issues.length > 0 && (
          <ul className="issues">
            {shownIssues.map((i, k) => <li key={k} className={i.level}><span className="mono">{i.line ? `line ${i.line}` : 'file'}</span><span className="lvl">{i.level}</span><span>{i.message}</span></li>)}
          </ul>
        )}
        {parsed.issues.length > 12 && <button type="button" className="linkish small" onClick={() => setShowAll(s => !s)}>{showAll ? 'Show fewer' : `Show all ${parsed.issues.length}`}</button>}
      </div>
    </div>
  )
}
