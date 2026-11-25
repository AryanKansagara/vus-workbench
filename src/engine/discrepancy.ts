// Inter-laboratory discrepancy reconciliation: what each lab called a variant, what their own criteria add up to, and where they differ.
import { CLASSES, combine, critLabel, DEF, POINTS, type Applied, type ClassName } from './acmg.ts'
import { dateKey, labVariantKey, type LabRow } from './io.ts'

export type Severity = 'concordant' | 'confidence' | 'clinically-significant' | 'major'
export const SEVERITY_LABEL: Record<Severity, string> = {
  concordant: 'Concordant',
  confidence: 'Confidence difference',
  'clinically-significant': 'Clinically significant',
  major: 'Major (P/LP vs B/LB)',
}

const tier = (c: ClassName) => c === 'Pathogenic' || c === 'Likely pathogenic' ? 'P' : c === 'VUS' ? 'V' : 'B'

/** Severity of a set of calls: P vs LP is a confidence difference; P/LP vs VUS changes management; P/LP vs B/LB is major. */
export function severity(classes: ClassName[]): Severity {
  const t = new Set(classes.map(tier))
  if (t.has('P') && t.has('B')) return 'major'
  if (t.has('P') && t.has('V')) return 'clinically-significant'
  if (new Set(classes).size > 1) return t.has('V') && t.size > 1 ? 'clinically-significant' : 'confidence'
  return 'concordant'
}

export interface LabCall { lab: string; cls: ClassName; criteria: Applied[]; points: number | null; implied: ClassName | null; mismatch: boolean; date?: string; line: number; /** The lab has another, different call and the dates do not say which is later. */ unclear?: boolean }
export interface CriterionDiff { code: string; byLab: Record<string, Applied['strength'] | null>; spread: number }
export interface Discrepancy { variant: string; gene?: string; calls: LabCall[]; severity: Severity; diffs: CriterionDiff[]; summary: string }

export function reconcile(rows: LabRow[]): Discrepancy[] {
  const groups = new Map<string, LabRow[]>()
  for (const r of rows) { const k = labVariantKey(r.variant, r.gene); groups.set(k, [...(groups.get(k) ?? []), r]) }
  const out: Discrepancy[] = []
  for (const [variant, list] of groups) {
    // Latest call per lab wins.
    const byLab = new Map<string, LabRow>()
    for (const r of list) { const p = byLab.get(r.lab); if (!p || dateKey(r.date).key >= dateKey(p.date).key) byLab.set(r.lab, r) }
    // A missing or shorter date ("2024" against "2024-03-01") does not settle the order.
    const unclear = (r: LabRow) => list.some(o => {
      if (o.lab !== r.lab || o.cls === r.cls) return false
      const a = dateKey(r.date).key, b = dateKey(o.date).key
      return a.startsWith(b) || b.startsWith(a)
    })
    const calls: LabCall[] = [...byLab.values()].map(r => {
      if (!r.criteria.length) return { lab: r.lab, cls: r.cls, criteria: [], points: null, implied: null, mismatch: false, date: r.date, line: r.line, ...(unclear(r) ? { unclear: true } : {}) }
      const c = combine(r.criteria)
      return { lab: r.lab, cls: r.cls, criteria: r.criteria, points: c.points, implied: c.cls, mismatch: c.cls !== r.cls, date: r.date, line: r.line, ...(unclear(r) ? { unclear: true } : {}) }
    }).sort((a, b) => CLASSES.indexOf(b.cls) - CLASSES.indexOf(a.cls))
    const codes = [...new Set(calls.flatMap(c => c.criteria.map(a => a.code)))].sort((a, b) => Object.keys(DEF).indexOf(a) - Object.keys(DEF).indexOf(b))
    const diffs: CriterionDiff[] = codes.map(code => {
      const m: Record<string, Applied['strength'] | null> = {}
      for (const c of calls) m[c.lab] = c.criteria.find(a => a.code === code)?.strength ?? null
      const pts = calls.filter(c => c.criteria.length).map(c => m[c.lab] ? POINTS[m[c.lab]!] : 0)
      return { code, byLab: m, spread: pts.length ? Math.max(...pts) - Math.min(...pts) : 0 }
    }).filter(d => d.spread > 0).sort((a, b) => b.spread - a.spread)
    const sev = severity(calls.map(c => c.cls))
    out.push({ variant, gene: list.find(r => r.gene)?.gene, calls, severity: sev, diffs, summary: summarize(calls, diffs, sev) })
  }
  const rank: Severity[] = ['major', 'clinically-significant', 'confidence', 'concordant']
  return out.sort((a, b) => rank.indexOf(a.severity) - rank.indexOf(b.severity) || a.variant.localeCompare(b.variant))
}

function summarize(calls: LabCall[], diffs: CriterionDiff[], sev: Severity): string {
  const order = calls.filter(c => c.unclear).map(c => `${c.lab} has different calls whose dates do not say which is later; ${c.cls} is shown.`)
  if (calls.length < 2) return [`Only one lab (${calls[0].lab}) has classified this variant.`, ...order].join(' ')
  const parts: string[] = []
  if (sev === 'concordant') parts.push(`All ${calls.length} labs agree: ${calls[0].cls}.`)
  else if (calls.length > 6) {
    // Too many labs to name each one (a ClinVar record can have fifty): count by class, name the outliers.
    const by = CLASSES.map(k => [k, calls.filter(c => c.cls === k)] as const).filter(([, l]) => l.length).sort((a, b) => b[1].length - a[1].length)
    parts.push(`${calls.length} labs: ${by.map(([k, l]) => `${l.length} ${k}`).join(', ')}.`)
    const few = by.filter(([, l]) => l.length <= 2 && l.length < by[0][1].length).flatMap(([k, l]) => l.map(c => `${c.lab} (${k}${dateKey(c.date).key ? `, ${dateKey(c.date).key.slice(0, 4)}` : ''})`))
    if (few.length && few.length <= 4) parts.push(`Outliers: ${few.join('; ')}.`)
  } else parts.push(calls.map(c => `${c.lab}: ${c.cls}`).join('; ') + '.')
  parts.push(...order)
  const mis = calls.filter(c => c.mismatch)
  for (const c of mis) parts.push(`${c.lab}'s own criteria total ${c.points} points (${c.implied}), not ${c.cls}.`)
  if (diffs.length) {
    const top = diffs.slice(0, 3).map(d => {
      const have = Object.entries(d.byLab).filter(([, s]) => s).map(([l, s]) => `${l} ${critLabel({ code: d.code, strength: s! })}`)
      // Only labs that stated criteria can be said not to have used one.
      const stated = new Set(calls.filter(c => c.criteria.length).map(c => c.lab))
      const lack = Object.entries(d.byLab).filter(([l, s]) => !s && stated.has(l)).map(([l]) => l)
      const list = (x: string[]) => x.length > 3 ? `${x.slice(0, 2).join(', ')} and ${x.length - 2} more` : x.join(', ')
      return `${d.code} (${list(have)}${lack.length ? `; not used by ${list(lack)}` : ''})`
    })
    parts.push(`Largest evidence differences: ${top.join('; ')}.`)
  }
  return parts.join(' ')
}
