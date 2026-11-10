// Inter-laboratory discrepancy reconciliation: what each lab called a variant, what their own criteria add up to, and where they differ.
import { CLASSES, combine, critLabel, DEF, POINTS, type Applied, type ClassName } from './acmg.ts'
import type { LabRow } from './io.ts'

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

export interface LabCall { lab: string; cls: ClassName; criteria: Applied[]; points: number | null; implied: ClassName | null; mismatch: boolean; date?: string; line: number }
export interface CriterionDiff { code: string; byLab: Record<string, Applied['strength'] | null>; spread: number }
export interface Discrepancy { variant: string; gene?: string; calls: LabCall[]; severity: Severity; diffs: CriterionDiff[]; summary: string }

export function reconcile(rows: LabRow[]): Discrepancy[] {
  const groups = new Map<string, LabRow[]>()
  for (const r of rows) groups.set(r.variant, [...(groups.get(r.variant) ?? []), r])
  const out: Discrepancy[] = []
  for (const [variant, list] of groups) {
    // Latest call per lab wins.
    const byLab = new Map<string, LabRow>()
    for (const r of list) { const p = byLab.get(r.lab); if (!p || (r.date ?? '') >= (p.date ?? '')) byLab.set(r.lab, r) }
    const calls: LabCall[] = [...byLab.values()].map(r => {
      if (!r.criteria.length) return { lab: r.lab, cls: r.cls, criteria: [], points: null, implied: null, mismatch: false, date: r.date, line: r.line }
      const c = combine(r.criteria)
      return { lab: r.lab, cls: r.cls, criteria: r.criteria, points: c.points, implied: c.cls, mismatch: c.cls !== r.cls, date: r.date, line: r.line }
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
  if (calls.length < 2) return `Only one lab (${calls[0].lab}) has classified this variant.`
  const parts: string[] = []
  if (sev === 'concordant') parts.push(`All ${calls.length} labs agree: ${calls[0].cls}.`)
  else parts.push(calls.map(c => `${c.lab}: ${c.cls}`).join('; ') + '.')
  const mis = calls.filter(c => c.mismatch)
  for (const c of mis) parts.push(`${c.lab}'s own criteria total ${c.points} points (${c.implied}), not ${c.cls}.`)
  if (diffs.length) {
    const top = diffs.slice(0, 3).map(d => {
      const have = Object.entries(d.byLab).filter(([, s]) => s).map(([l, s]) => `${l} ${critLabel({ code: d.code, strength: s! })}`)
      const lack = Object.entries(d.byLab).filter(([, s]) => !s).map(([l]) => l)
      return `${d.code} (${have.join(', ')}${lack.length ? `; not used by ${lack.join(', ')}` : ''})`
    })
    parts.push(`Largest evidence differences: ${top.join('; ')}.`)
  }
  return parts.join(' ')
}
