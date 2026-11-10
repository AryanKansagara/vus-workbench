import { useState } from 'react'
import { getDefaults, setDefaults, type Defaults } from './store.ts'
import { setTheme, toast, useDark } from './ui.tsx'
import { useSync } from './sync.ts'

export default function Settings() {
  const [d, setD] = useState<Defaults>(getDefaults)
  const dark = useDark()
  const s = useSync()
  const save = (patch: Partial<Defaults>) => { const n = { ...d, ...patch }; setD(n); setDefaults(n) }
  return (
    <section className="page narrow">
      <p className="eyebrow">Settings</p>
      <h1>Settings</h1>
      <div className="card pad stack">
        <h2 className="h4">Defaults for new close packs</h2>
        <div className="row gap wrap">
          <label className="field"><span>Dollar threshold</span><input type="number" min={0} step={500} value={d.abs} onChange={e => save({ abs: Math.max(0, +e.target.value || 0) })} /></label>
          <label className="field"><span>Percent threshold</span><input type="number" min={0} step={1} value={Math.round(d.pct * 1000) / 10} onChange={e => save({ pct: Math.max(0, +e.target.value || 0) / 100 })} /></label>
          <label className="field"><span>Compare against</span><select value={d.basis} onChange={e => save({ basis: e.target.value as Defaults['basis'] })}><option value="budget">Budget</option><option value="mom">Prior month</option><option value="yoy">Same month last year</option></select></label>
          <label className="field"><span>Reporting currency</span><input className="short" value={d.reportingCurrency} onChange={e => save({ reportingCurrency: e.target.value.toUpperCase().slice(0, 3) })} /></label>
        </div>
        <label className="field"><span>Preparer name on board packs</span><input value={d.preparer} onChange={e => save({ preparer: e.target.value })} placeholder="Jordan Lee, FP&A" /></label>
        <p className="muted small">Saved as you type. Existing packs keep their own thresholds (Data tab).</p>
      </div>
      <div className="card pad stack">
        <h2 className="h4">Appearance</h2>
        <div className="seg" role="group" aria-label="Theme">
          <button type="button" aria-pressed={!dark} onClick={() => setTheme('light')}>Light</button>
          <button type="button" aria-pressed={dark} onClick={() => setTheme('dark')}>Dark</button>
        </div>
      </div>
      <div className="card pad stack">
        <h2 className="h4">Sync</h2>
        <p className="muted small">{s.user ? `Signed in as ${s.user.email}. Use the account menu in the header to sync, export, or delete cloud data.` : 'Signed out. Everything is stored in this browser. Sign in from the account menu to sync across devices.'}</p>
        <button type="button" className="btn ghost sm" onClick={() => { localStorage.removeItem('variance-pack.defaults'); setD(getDefaults()); toast('Defaults restored.') }}>Restore default settings</button>
      </div>
    </section>
  )
}
