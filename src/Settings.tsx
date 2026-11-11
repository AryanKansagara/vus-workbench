import { useState } from 'react'
import { getDefaults, resetDefaults, setDefaults, type Defaults } from './store.ts'
import { setTheme, toast, useDark } from './ui.tsx'
import { useSync } from './sync.ts'

export default function Settings() {
  const [d, setD] = useState<Defaults>(getDefaults)
  const dark = useDark()
  const s = useSync()
  const save = (patch: Partial<Defaults>) => { const n = { ...d, ...patch }; setD(n); setDefaults(n) }
  const freq = (k: 'pm2Max' | 'bs1Min' | 'ba1Min', v: string) => { const n = Number(v); if (v !== '' && !Number.isNaN(n) && n >= 0 && n <= 1) save({ [k]: n }) }
  return (
    <section className="page narrow">
      <p className="eyebrow">Settings</p>
      <h1>Settings</h1>
      <div className="card pad stack">
        <h2 className="h4">Curator</h2>
        <label className="field"><span>Your name (recorded in every evidence log entry)</span><input value={d.editor} onChange={e => save({ editor: e.target.value })} placeholder="Dr Sam Patel, clinical scientist" /></label>
        <label className="field"><span>Laboratory or group (on reports and curated-call exports)</span><input value={d.institution} onChange={e => save({ institution: e.target.value })} placeholder="Molecular Genetics, Example Hospital" /></label>
      </div>
      <div className="card pad stack">
        <h2 className="h4">Defaults for genes without settings</h2>
        <div className="row gap wrap">
          <label className="field"><span>Inheritance</span><select value={d.inheritance} onChange={e => save({ inheritance: e.target.value as Defaults['inheritance'] })}><option value="AD">Dominant</option><option value="AR">Recessive</option><option value="XL">X-linked</option></select></label>
          <label className="field"><span>PM2 at or below</span><input className="short" type="number" step="any" min={0} max={1} defaultValue={d.pm2Max} onChange={e => freq('pm2Max', e.target.value)} /></label>
          <label className="field"><span>BS1 at or above</span><input className="short" type="number" step="any" min={0} max={1} defaultValue={d.bs1Min} onChange={e => freq('bs1Min', e.target.value)} /></label>
          <label className="field"><span>BA1 at or above</span><input className="short" type="number" step="any" min={0} max={1} defaultValue={d.ba1Min} onChange={e => freq('ba1Min', e.target.value)} /></label>
        </div>
        <p className="muted small">Saved as you type. Per-gene settings on each worklist's Genes tab take precedence.</p>
      </div>
      <div className="card pad stack">
        <h2 className="h4">Appearance</h2>
        <div className="seg" role="group" aria-label="Theme">
          <button type="button" aria-pressed={!dark} onClick={() => setTheme('light')}>Light</button>
          <button type="button" aria-pressed={dark} onClick={() => setTheme('dark')}>Dark</button>
        </div>
      </div>
      <div className="card pad stack">
        <h2 className="h4">Sync and data</h2>
        <p className="muted small">{s.user ? `Signed in as ${s.user.email}. Use the account menu in the header to sync, export, or delete cloud data.` : 'Signed out. Everything is stored in this browser. Sign in from the account menu to sync across devices.'}</p>
        <p className="muted small">Keyboard: <kbd>/</kbd> search the worklist, <kbd>[</kbd> <kbd>]</kbd> previous and next variant, <kbd>⌘Z</kbd> undo the last change while its toast is showing.</p>
        <button type="button" className="btn ghost sm" onClick={() => { resetDefaults(); setD(getDefaults()); toast('Defaults restored.') }}>Restore default settings</button>
      </div>
    </section>
  )
}
