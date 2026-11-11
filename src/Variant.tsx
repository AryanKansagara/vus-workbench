import { useEffect, useState } from 'react'
import { CRITERIA, DEF, fmtAf, pvs1, STRENGTH_SHORT, STRENGTHS, type CriterionState, type Strength } from './engine/acmg.ts'
import type { Evaluated } from './engine/pipeline.ts'
import { CONSEQUENCE_LABEL } from './engine/transcript.ts'
import { curate, geneCfg } from './projectData.ts'
import { editorName, saveProject, uuid, type Project } from './store.ts'
import { Empty, fmtDateTime, go, isTyping, Link, toast } from './ui.tsx'
import { ClassBadge, Meter, SiteSeq, Track } from './viz.tsx'

export default function VariantView({ project, evaluated, vid }: { project: Project; evaluated: Evaluated[]; vid: string }) {
  const idx = evaluated.findIndex(e => e.an.rec.id === vid)
  const e = evaluated[idx]
  useEffect(() => {
    const onKey = (k: KeyboardEvent) => {
      if (isTyping(k) || k.metaKey || k.ctrlKey) return
      const to = k.key === '[' ? evaluated[idx - 1] : k.key === ']' ? evaluated[idx + 1] : null
      if (to) go(`/app/p/${project.id}/v/${to.an.rec.id}`)
    }
    addEventListener('keydown', onKey)
    return () => removeEventListener('keydown', onKey)
  }, [idx, evaluated, project.id])
  if (!e) return <Empty title="Variant not found"><Link href={`/app/p/${project.id}`} className="btn primary">Back to the worklist</Link></Empty>
  const { an, result } = e
  const a = an.ann
  const cur = project.curation[vid]
  const ev = an.rec.evidence
  const gene = geneCfg(project, an.gene)
  const act = (action: string, detail: string, fn: Parameters<typeof curate>[4], undo = false) => {
    const prev = project
    void curate(project, vid, action, detail, fn).then(() => { if (undo) toast(`${action}: ${detail}`, { undo: () => void saveProject(prev) }) })
  }
  const setCrit = (c: CriterionState | undefined, code: string, patch: Partial<CriterionState>, action: string, detail: string) =>
    act(action, detail, ({ cur }) => {
      const base: CriterionState = c ?? { code, met: false, strength: DEF[code].strength, evidence: '', source: 'curator', editor: '', at: 0 }
      cur.overrides[code] = { ...base, ...patch, source: 'curator', editor: editorName(), at: Date.now() }
    }, true)
  const pv = a ? pvs1(a, gene) : null

  return (
    <div className="variant">
      <div className="var-head">
        <div>
          <div className="row gap-sm wrap muted small">
            <span>{an.tx ? `${an.tx.id} · ${an.tx.source}` : an.rec.transcript ?? 'No transcript model'}</span>
            {an.g && <span className="mono">{an.g.chrom}:{an.g.pos} {an.g.ref}&gt;{an.g.alt}</span>}
            {an.rec.gt && <span>GT {an.rec.gt}{an.rec.sample ? ` · ${an.rec.sample}` : ''}</span>}
            <span>from {an.rec.source.file}, line {an.rec.source.line}</span>
          </div>
          <h2 className="var-title"><span>{an.gene}</span> <span className="mono">{a?.hgvsC.replace(/^.*:/, '') ?? an.rec.key}</span> {a?.hgvsP && <span className="mono muted">{a.hgvsP}</span>}</h2>
        </div>
        <div className="row gap-sm wrap">
          <label className="field inline"><span>Status</span>
            <select value={cur?.status ?? 'new'} onChange={x => act('Status', `Marked ${x.target.value}`, ({ cur }) => { cur.status = x.target.value as 'new' })}>
              <option value="new">New</option><option value="in-review">In review</option><option value="curated">Curated</option>
            </select>
          </label>
          <Link href={`/app/p/${project.id}/report/${vid}`} className="btn primary sm">Curation report</Link>
          <div className="row gap-xs">
            <button type="button" className="icon-btn" aria-label="Previous variant ([)" disabled={idx <= 0} onClick={() => go(`/app/p/${project.id}/v/${evaluated[idx - 1].an.rec.id}`)}>‹</button>
            <span className="small muted">{idx + 1}/{evaluated.length}</span>
            <button type="button" className="icon-btn" aria-label="Next variant (])" disabled={idx >= evaluated.length - 1} onClick={() => go(`/app/p/${project.id}/v/${evaluated[idx + 1].an.rec.id}`)}>›</button>
          </div>
        </div>
      </div>

      <div className="var-grid">
        <div className="stack">
          <div className="card pad">
            <Track tx={an.tx} ann={a} gpos={an.g?.pos} />
            {a && (
              <dl className="facts">
                <div><dt>Consequence</dt><dd><span className={`csq csq-${a.consequence}`}>{CONSEQUENCE_LABEL[a.consequence]}</span>{a.spliceRegion && <span className="csq csq-splice_region">+ splice region</span>}</dd></div>
                {a.exon && <div><dt>Location</dt><dd>Exon {a.exon}{a.exonCount ? ` of ${a.exonCount}` : ''}</dd></div>}
                {a.intron && <div><dt>Location</dt><dd>Intron {a.intron}, {Math.abs(a.intronOffset ?? 0)} nt from the {a.intronOffset! > 0 ? 'donor' : 'acceptor'}</dd></div>}
                {a.ptcAa && <div><dt>Premature stop</dt><dd>Codon {a.ptcAa}</dd></div>}
                {a.nmd !== undefined && a.nmd !== null && <div><dt>NMD</dt><dd>{a.nmd ? 'Predicted' : 'Predicted to escape'}</dd></div>}
                {a.removedFraction != null && <div><dt>Protein affected</dt><dd>{(a.removedFraction * 100).toFixed(1)}%</dd></div>}
                <div><dt>Gene</dt><dd>{gene.inheritance}, LoF {gene.lofMechanism ? 'is' : 'is not'} a mechanism · <Link href={`/app/p/${project.id}/genes`}>edit</Link></dd></div>
              </dl>
            )}
            {[...(a?.notes ?? []), ...an.problems].map(n => <p key={n} className="note small">{n}</p>)}
          </div>

          <div className="card pad">
            <h3 className="h5">Splice-site check</h3>
            {an.splice?.canonical ? (
              <>
                <p className={`splice-verdict ${an.splice.verdict}`}>{an.splice.summary}</p>
                <div className="row gap wrap">
                  <div><div className="muted small">Native {an.splice.canonical.site}, exon {an.splice.canonical.exon}</div><SiteSeq refSeq={an.splice.canonical.refSeq} altSeq={an.splice.canonical.altSeq} site={an.splice.canonical.site} /></div>
                  {an.splice.cryptic && <div><div className="muted small">Candidate new {an.splice.cryptic.site} ({an.splice.cryptic.refScore} → {an.splice.cryptic.altScore} bits)</div><SiteSeq refSeq={an.splice.cryptic.refSeq} altSeq={an.splice.cryptic.altSeq} site={an.splice.cryptic.site} /></div>}
                </div>
                <p className="muted small">Local position weight matrix (MaxEnt-style approximation), log2 odds in bits. Use it as supporting context; SpliceAI or RNA evidence takes precedence.</p>
              </>
            ) : <p className="muted small">{an.tx ? 'Reference sequence for this contig is not loaded. Import a FASTA region (samtools faidx) on the Import tab to score splice sites.' : 'Needs a transcript model and reference sequence.'}</p>}
          </div>

          <div className="card pad">
            <h3 className="h5">Evidence inputs</h3>
            <div className="row gap wrap">
              <NumField label="Popmax AF" value={ev.popmaxAf} step="any" hint={`PM2 ≤ ${fmtAf(gene.pm2Max)}, BS1 ≥ ${fmtAf(gene.bs1Min)}, BA1 ≥ ${fmtAf(gene.ba1Min)}`} onCommit={v => act('Evidence updated', `Popmax AF ${ev.popmaxAf ?? 'blank'} → ${v ?? 'blank'}`, ({ rec }) => { rec.evidence.popmaxAf = v })} />
              <NumField label="REVEL" value={ev.revel} step="0.001" hint="Missense only, Pejaver 2022 calibration" onCommit={v => act('Evidence updated', `REVEL ${ev.revel ?? 'blank'} → ${v ?? 'blank'}`, ({ rec }) => { rec.evidence.revel = v })} />
              <NumField label="SpliceAI max Δ" value={ev.spliceAi} step="0.01" hint="PP3 ≥ 0.2, BP4 ≤ 0.1" onCommit={v => act('Evidence updated', `SpliceAI ${ev.spliceAi ?? 'blank'} → ${v ?? 'blank'}`, ({ rec }) => { rec.evidence.spliceAi = v })} />
            </div>
            <h4 className="h5 mt">Functional studies (PS3 / BS3)</h4>
            {ev.functional.length ? (
              <ul className="func-list">
                {ev.functional.map(f => (
                  <li key={f.id}>
                    <span className={`pill ${f.result === 'abnormal' ? 'p' : 'b'}`}>{f.result === 'abnormal' ? 'Damaging' : 'Normal'}</span>
                    <div><b>{f.assay}</b>{f.oddsPath !== null && <span className="mono small"> · OddsPath {f.oddsPath}</span>}<div className="muted small">{f.note} {f.editor}, {fmtDateTime(f.at)}</div></div>
                    <button type="button" className="btn danger-ghost sm" onClick={() => act('Functional evidence removed', f.assay, ({ rec }) => { rec.evidence.functional = rec.evidence.functional.filter(x => x.id !== f.id) }, true)}>Remove</button>
                  </li>
                ))}
              </ul>
            ) : <p className="muted small">None recorded. OddsPath maps to strength (Brnich 2019): over 18.7 Strong, over 4.3 Moderate, over 2.1 Supporting.</p>}
            <FunctionalForm onAdd={f => act('Functional evidence added', `${f.assay}${f.oddsPath !== null ? `, OddsPath ${f.oddsPath}` : ''}`, ({ rec }) => { rec.evidence.functional.push({ ...f, id: uuid(), editor: editorName(), at: Date.now() }) })} />
          </div>

          {e.partners.length > 0 && (
            <div className="card pad">
              <h3 className="h5">Other variants in {an.gene} ({an.rec.sample ?? 'same sample'})</h3>
              <ul className="plain">{e.partners.map(p => <li key={p.label}><span className={`phase ${p.phase}`}>{p.phase}</span> {p.label} <ClassBadge cls={p.cls} small /></li>)}</ul>
              <p className="muted small">{gene.inheritance === 'AR' ? 'Recessive gene: in trans with P/LP gives PM3; in cis gives BP2.' : 'Dominant gene: in trans or in cis with a pathogenic variant gives BP2.'}</p>
            </div>
          )}
        </div>

        <aside className="stack crit-col">
          <div className="card pad sticky">
            <div className="row between"><ClassBadge cls={result.cls} points={result.points} /><span className="mono big">{result.points > 0 ? '+' : ''}{result.points} pts</span></div>
            <Meter points={result.points} cls={result.cls} />
            {result.standAloneBenign && <p className="small note">BA1 is stand-alone: Benign regardless of other criteria.</p>}
            {result.excluded.map(x => <p key={x.code} className="small note">{x.code} not counted. {x.why}</p>)}
          </div>
          <div className="card pad">
            <div className="row between"><h3 className="h5">Criteria</h3><span className="muted small">Engine suggestions in teal, your decisions in orange</span></div>
            {pv?.strength && <details className="pvs1"><summary>PVS1 decision path</summary><ol>{pv.path.map(s => <li key={s}>{s}</li>)}</ol></details>}
            {(['P', 'B'] as const).map(dir => (
              <div key={dir} className="crit-group">
                <div className="crit-group-h">{dir === 'P' ? 'Pathogenic' : 'Benign'}</div>
                {CRITERIA.filter(c => c.dir === dir).map(d => {
                  const c = e.criteria.find(x => x.code === d.code)
                  const sug = e.suggestions.find(s => s.code === d.code)
                  const over = cur?.overrides[d.code]
                  return (
                    <div key={d.code} className={`crit-row${c?.met ? ' met' : ''}${over ? ' over' : sug ? ' sug' : ''}`}>
                      <label className="crit-toggle" title={d.label}>
                        <input type="checkbox" checked={!!c?.met} onChange={x => setCrit(c, d.code, { met: x.target.checked, ...(x.target.checked && !c ? { strength: d.strength } : {}) }, x.target.checked ? 'Criterion applied' : 'Criterion removed', `${d.code}${x.target.checked ? '' : ` (${sug ? 'overrides engine' : 'manual'})`}`)} />
                        <b>{d.code}</b>
                      </label>
                      {d.code !== 'BA1' ? (
                        <select aria-label={`${d.code} strength`} value={c?.strength ?? d.strength} disabled={!c?.met} onChange={x => setCrit(c, d.code, { strength: x.target.value as Strength }, 'Strength changed', `${d.code} ${STRENGTH_SHORT[c?.strength ?? d.strength]} → ${STRENGTH_SHORT[x.target.value as Strength]}`)}>
                          {STRENGTHS.map(s => <option key={s} value={s}>{STRENGTH_SHORT[s]}</option>)}
                        </select>
                      ) : <span className="small muted">Stand-alone</span>}
                      <div className="crit-ev">
                        <span className="small">{d.label}</span>
                        {c && c.evidence && <span className="small ev-text">{c.evidence} <i>{c.source === 'engine' ? 'engine' : `${c.editor}, ${fmtDateTime(c.at)}`}</i></span>}
                        {c?.met && <EvidenceEdit value={over?.evidence ?? ''} onCommit={t => setCrit(c, d.code, { evidence: t }, 'Evidence note', `${d.code}: ${t.slice(0, 80)}`)} />}
                        {over && sug && <button type="button" className="linkish small" onClick={() => act('Reset to engine', d.code, ({ cur }) => { delete cur.overrides[d.code] }, true)}>Reset to engine suggestion ({sug.met ? `${d.code}_${sug.strength}` : 'not met'})</button>}
                        {over && !sug && <button type="button" className="linkish small" onClick={() => act('Criterion cleared', d.code, ({ cur }) => { delete cur.overrides[d.code] }, true)}>Clear</button>}
                      </div>
                    </div>
                  )
                })}
              </div>
            ))}
          </div>
          <div className="card pad">
            <label className="field"><span>Curator note</span><NoteBox value={cur?.note ?? ''} onCommit={t => act('Note', t.slice(0, 120) || 'cleared', ({ cur }) => { cur.note = t })} /></label>
          </div>
        </aside>
      </div>

      <div className="var-logs">
        <div className="card pad">
          <h3 className="h5">Reclassification history</h3>
          {cur?.history.length ? (
            <ol className="timeline">
              {cur.history.map((h, i) => <li key={i}><span className="muted small">{fmtDateTime(h.at)} · {h.editor}</span><div>{h.from ? <><ClassBadge cls={h.from} small /> → </> : null}<ClassBadge cls={h.to} small /> <span className="mono small">{h.points > 0 ? '+' : ''}{h.points}</span></div><span className="small">{h.reason}</span></li>)}
            </ol>
          ) : <p className="muted small">Not reclassified. The engine call is {result.cls} at {result.points} points.</p>}
        </div>
        <div className="card pad">
          <h3 className="h5">Evidence log</h3>
          {cur?.audit.length ? (
            <ol className="audit">{[...cur.audit].reverse().map((x, i) => <li key={i}><span className="muted small">{fmtDateTime(x.at)}</span> <b>{x.editor}</b> <span className="tag">{x.action}</span> <span className="small">{x.detail}</span></li>)}</ol>
          ) : <p className="muted small">No curator actions yet. Every change you make here is logged with your name (set it in Settings) and time.</p>}
        </div>
      </div>
    </div>
  )
}

function NumField({ label, value, onCommit, step, hint }: { label: string; value: number | null; onCommit(v: number | null): void; step: string; hint: string }) {
  const [v, setV] = useState(value === null ? '' : String(value))
  useEffect(() => setV(value === null ? '' : String(value)), [value])
  const commit = () => {
    const n = v.trim() === '' ? null : Number(v)
    if (n !== null && (Number.isNaN(n) || n < 0 || n > 1)) { toast(`${label} must be between 0 and 1.`, { tone: 'error' }); setV(value === null ? '' : String(value)); return }
    if (n !== value) onCommit(n)
  }
  return <label className="field num-field"><span>{label}</span><input type="number" min={0} max={1} step={step} value={v} onChange={e => setV(e.target.value)} onBlur={commit} onKeyDown={e => { if (e.key === 'Enter') commit() }} placeholder="not entered" /><small className="muted">{hint}</small></label>
}

function EvidenceEdit({ value, onCommit }: { value: string; onCommit(t: string): void }) {
  const [v, setV] = useState(value)
  return <input className="ev-input" aria-label="Evidence note" value={v} placeholder="Add your evidence or citation" onChange={e => setV(e.target.value)} onBlur={() => { if (v.trim() !== value) onCommit(v.trim()) }} onKeyDown={e => { if (e.key === 'Enter') (e.target as HTMLInputElement).blur() }} />
}

function NoteBox({ value, onCommit }: { value: string; onCommit(t: string): void }) {
  const [v, setV] = useState(value)
  useEffect(() => setV(value), [value])
  return <textarea rows={3} value={v} onChange={e => setV(e.target.value)} onBlur={() => { if (v !== value) onCommit(v) }} placeholder="Case context, open questions, what would move this variant" />
}

function FunctionalForm({ onAdd }: { onAdd(f: { assay: string; result: 'abnormal' | 'normal'; oddsPath: number | null; note: string }): void }) {
  const [open, setOpen] = useState(false)
  if (!open) return <button type="button" className="btn ghost sm" onClick={() => setOpen(true)}>Add functional result</button>
  return (
    <form className="func-form" onSubmit={e => {
      e.preventDefault()
      const f = new FormData(e.currentTarget)
      const o = String(f.get('odds') ?? '').trim()
      const odds = o ? Number(o) : null
      if (odds !== null && (Number.isNaN(odds) || odds <= 0)) { toast('OddsPath must be a positive number.', { tone: 'error' }); return }
      onAdd({ assay: String(f.get('assay')).trim(), result: f.get('result') as 'abnormal', oddsPath: odds, note: String(f.get('note') ?? '').trim() })
      setOpen(false)
    }}>
      <label className="field"><span>Assay</span><input name="assay" required placeholder="Minigene, RNA-seq, saturation genome editing" /></label>
      <label className="field"><span>Result</span><select name="result"><option value="abnormal">Damaging (PS3)</option><option value="normal">Normal (BS3)</option></select></label>
      <label className="field"><span>OddsPath</span><input name="odds" type="number" step="any" min={0} placeholder="optional" /></label>
      <label className="field"><span>Note / citation</span><input name="note" placeholder="PMID, internal study ID" /></label>
      <div className="row gap-sm"><button className="btn primary sm">Add</button><button type="button" className="btn ghost sm" onClick={() => setOpen(false)}>Cancel</button></div>
    </form>
  )
}
