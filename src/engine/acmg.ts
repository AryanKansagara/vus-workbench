// ACMG/AMP criteria with ClinGen SVI refinements, combined with the Bayesian points system (Tavtigian et al. 2020).
import type { Annotation } from './transcript.ts'
import type { SpliceReport } from './splice.ts'

export type Strength = 'StandAlone' | 'VeryStrong' | 'Strong' | 'Moderate' | 'Supporting'
export const STRENGTHS: Strength[] = ['VeryStrong', 'Strong', 'Moderate', 'Supporting']
export const STRENGTH_SHORT: Record<Strength, string> = { StandAlone: 'Stand-alone', VeryStrong: 'Very strong', Strong: 'Strong', Moderate: 'Moderate', Supporting: 'Supporting' }
export const POINTS: Record<Strength, number> = { StandAlone: 8, VeryStrong: 8, Strong: 4, Moderate: 2, Supporting: 1 }

export interface CriterionDef { code: string; dir: 'P' | 'B'; strength: Strength; label: string }
export const CRITERIA: CriterionDef[] = [
  { code: 'PVS1', dir: 'P', strength: 'VeryStrong', label: 'Null variant in a gene where loss of function causes disease' },
  { code: 'PS1', dir: 'P', strength: 'Strong', label: 'Same amino acid change as an established pathogenic variant' },
  { code: 'PS2', dir: 'P', strength: 'Strong', label: 'De novo, maternity and paternity confirmed' },
  { code: 'PS3', dir: 'P', strength: 'Strong', label: 'Well-established functional studies show a damaging effect' },
  { code: 'PS4', dir: 'P', strength: 'Strong', label: 'Prevalence in affected individuals significantly increased over controls' },
  { code: 'PM1', dir: 'P', strength: 'Moderate', label: 'Mutational hot spot or critical functional domain' },
  { code: 'PM2', dir: 'P', strength: 'Supporting', label: 'Absent or extremely rare in population databases' },
  { code: 'PM3', dir: 'P', strength: 'Moderate', label: 'Recessive: detected in trans with a pathogenic variant' },
  { code: 'PM4', dir: 'P', strength: 'Moderate', label: 'Protein length change (in-frame indel in a non-repeat region, stop-loss)' },
  { code: 'PM5', dir: 'P', strength: 'Moderate', label: 'Novel missense at a residue where a different pathogenic missense is known' },
  { code: 'PM6', dir: 'P', strength: 'Moderate', label: 'Assumed de novo, without confirmation of parentage' },
  { code: 'PP1', dir: 'P', strength: 'Supporting', label: 'Co-segregation with disease in multiple affected family members' },
  { code: 'PP2', dir: 'P', strength: 'Supporting', label: 'Missense in a gene with low benign missense variation' },
  { code: 'PP3', dir: 'P', strength: 'Supporting', label: 'Computational evidence supports a deleterious effect' },
  { code: 'PP4', dir: 'P', strength: 'Supporting', label: 'Phenotype or family history highly specific for the gene' },
  { code: 'BA1', dir: 'B', strength: 'StandAlone', label: 'Allele frequency above 5% (or the gene-specific BA1 threshold)' },
  { code: 'BS1', dir: 'B', strength: 'Strong', label: 'Allele frequency greater than expected for the disorder' },
  { code: 'BS2', dir: 'B', strength: 'Strong', label: 'Observed in a healthy adult with full penetrance expected' },
  { code: 'BS3', dir: 'B', strength: 'Strong', label: 'Well-established functional studies show no damaging effect' },
  { code: 'BS4', dir: 'B', strength: 'Strong', label: 'Lack of segregation in affected family members' },
  { code: 'BP1', dir: 'B', strength: 'Supporting', label: 'Missense in a gene where only truncating variants cause disease' },
  { code: 'BP2', dir: 'B', strength: 'Supporting', label: 'In trans with a pathogenic variant (dominant) or in cis (any inheritance)' },
  { code: 'BP3', dir: 'B', strength: 'Supporting', label: 'In-frame indel in a repetitive region without known function' },
  { code: 'BP4', dir: 'B', strength: 'Supporting', label: 'Computational evidence suggests no impact' },
  { code: 'BP5', dir: 'B', strength: 'Supporting', label: 'Found in a case with an alternate molecular basis for disease' },
  { code: 'BP7', dir: 'B', strength: 'Supporting', label: 'Synonymous (or deep intronic) with no predicted splice impact' },
]
export const DEF = Object.fromEntries(CRITERIA.map(c => [c.code, c])) as Record<string, CriterionDef>

export type ClassName = 'Pathogenic' | 'Likely pathogenic' | 'VUS' | 'Likely benign' | 'Benign'
export const CLASSES: ClassName[] = ['Benign', 'Likely benign', 'VUS', 'Likely pathogenic', 'Pathogenic']
export const CLASS_SHORT: Record<ClassName, string> = { Pathogenic: 'P', 'Likely pathogenic': 'LP', VUS: 'VUS', 'Likely benign': 'LB', Benign: 'B' }

/** Point thresholds from Tavtigian 2020: P >= 10, LP 6-9, VUS 0-5, LB -1 to -6, B <= -7. */
export function classify(points: number): ClassName {
  if (points >= 10) return 'Pathogenic'
  if (points >= 6) return 'Likely pathogenic'
  if (points >= 0) return 'VUS'
  if (points >= -6) return 'Likely benign'
  return 'Benign'
}

/** ClinGen-style VUS temperature, useful when deciding which VUS to revisit first. */
export function vusTier(points: number): 'hot' | 'warm' | 'tepid' | 'cool' | null {
  if (points < 0 || points > 5) return null
  return points >= 4 ? 'hot' : points >= 2 ? 'warm' : points >= 1 ? 'tepid' : 'cool'
}

export interface Applied { code: string; strength: Strength }
export interface Combined { points: number; cls: ClassName; counted: (Applied & { points: number })[]; excluded: (Applied & { why: string })[]; standAloneBenign: boolean }

/** Combines met criteria into points and a class, applying the double-counting rules. */
export function combine(applied: Applied[]): Combined {
  const has = (c: string) => applied.some(a => a.code === c)
  const excluded: Combined['excluded'] = []
  const counted: Combined['counted'] = []
  for (const a of applied) {
    let why = ''
    if (a.code === 'PP3' && has('PVS1')) why = 'Not combined with PVS1 (ClinGen SVI: the null effect is already counted).'
    else if (a.code === 'PM4' && has('PVS1')) why = 'PM4 and PVS1 describe the same protein-length change.'
    else if (a.code === 'BP4' && has('PP3')) why = 'PP3 and BP4 conflict; resolve the computational evidence first.'
    else if (a.code === 'PP3' && has('BP4')) why = 'PP3 and BP4 conflict; resolve the computational evidence first.'
    if (why) { excluded.push({ ...a, why }); continue }
    counted.push({ ...a, points: POINTS[a.strength] * (DEF[a.code]?.dir === 'B' ? -1 : 1) })
  }
  const points = counted.reduce((s, c) => s + c.points, 0)
  const standAloneBenign = counted.some(c => c.code === 'BA1')
  return { points, cls: standAloneBenign ? 'Benign' : classify(points), counted, excluded, standAloneBenign }
}

// ---------- Inputs ----------

export interface GeneCfg { inheritance: 'AD' | 'AR' | 'XL'; lofMechanism: boolean; pm2Max: number; bs1Min: number; ba1Min: number; criticalRegions?: [number, number][]; proteinLength?: number }
export const DEFAULT_GENE: GeneCfg = { inheritance: 'AD', lofMechanism: true, pm2Max: 0.0001, bs1Min: 0.01, ba1Min: 0.05 }

export interface FunctionalEntry { id: string; assay: string; result: 'abnormal' | 'normal'; oddsPath: number | null; note: string; editor: string; at: number }
export interface Evidence {
  /** Highest continental-population allele frequency (gnomAD popmax / grpmax). 0 = absent. */
  popmaxAf: number | null
  revel: number | null
  spliceAi: number | null
  functional: FunctionalEntry[]
}
export const EMPTY_EVIDENCE: Evidence = { popmaxAf: null, revel: null, spliceAi: null, functional: [] }

export type Phase = 'trans' | 'cis' | 'unknown' | 'homozygous'
export interface PhasePartner { label: string; cls: ClassName; phase: Phase }

export interface Suggestion { code: string; met: boolean; strength: Strength; reason: string }

// ---------- PVS1 (Abou Tayoun et al. 2018, simplified) ----------

export interface Pvs1Result { strength: Strength | null; path: string[] }

export function pvs1(a: Annotation, gene: GeneCfg): Pvs1Result {
  const path: string[] = []
  const c = a.consequence
  if (!['frameshift', 'stop_gained', 'splice_donor', 'splice_acceptor', 'start_lost'].includes(c)) return { strength: null, path: ['Not a null variant type.'] }
  if (!gene.lofMechanism) return { strength: null, path: ['Loss of function is not an established disease mechanism for this gene (gene setting).'] }
  path.push('Null variant; loss of function is a known mechanism.')
  if (c === 'start_lost') { path.push('Initiation codon: no evidence modelled for an alternative start; PVS1_Moderate.'); return { strength: 'Moderate', path } }
  const removed = a.removedFraction
  // A critical region downstream of the stop is lost with the truncated protein.
  const critical = a.ptcAa !== undefined && (gene.criticalRegions ?? []).some(([, hi]) => hi >= a.ptcAa!)
  const escape = (what: string): Pvs1Result => {
    path.push(what)
    if (critical) { path.push('Truncated or skipped region overlaps a critical region (gene setting): PVS1_Strong.'); return { strength: 'Strong', path } }
    if (removed == null) { path.push('Fraction of protein removed is unknown: PVS1_Moderate pending review.'); return { strength: 'Moderate', path } }
    if (removed > 0.1) { path.push(`Removes ${(removed * 100).toFixed(1)}% of the protein (>10%): PVS1_Strong.`); return { strength: 'Strong', path } }
    path.push(`Removes ${(removed * 100).toFixed(1)}% of the protein (<=10%): PVS1_Moderate.`)
    return { strength: 'Moderate', path }
  }
  if (c === 'frameshift' || c === 'stop_gained') {
    if (a.nmd === true) { path.push(`Premature stop at codon ${a.ptcAa} is upstream of the last 50 nt of the penultimate exon: predicted NMD. PVS1.`); return { strength: 'VeryStrong', path } }
    if (a.nmd === false) return escape(`Premature stop at codon ${a.ptcAa} is in the last exon or the last 50 nt of the penultimate exon: predicted to escape NMD.`)
    return escape('NMD cannot be predicted without a transcript model or stop position.')
  }
  // Canonical splice site.
  if (a.skipInFrame === null || a.skipInFrame === undefined) return escape('Exon-skipping outcome unknown without a transcript model.')
  if (!a.skipInFrame) {
    if (a.nmd) { path.push('Exon skipping is predicted to shift the frame and trigger NMD: PVS1.'); return { strength: 'VeryStrong', path } }
    return escape('Exon skipping shifts the frame but the new stop likely escapes NMD.')
  }
  return escape('Exon skipping is predicted to preserve the reading frame.')
}

// ---------- Computational evidence thresholds ----------

/** REVEL calibration (Pejaver et al. 2022). */
export function revelStrength(s: number): { code: 'PP3' | 'BP4'; strength: Strength } | null {
  if (s >= 0.932) return { code: 'PP3', strength: 'Strong' }
  if (s >= 0.773) return { code: 'PP3', strength: 'Moderate' }
  if (s >= 0.644) return { code: 'PP3', strength: 'Supporting' }
  if (s <= 0.003) return { code: 'BP4', strength: 'Strong' }
  if (s <= 0.016) return { code: 'BP4', strength: 'Moderate' }
  if (s <= 0.290) return { code: 'BP4', strength: 'Supporting' }
  return null
}

/** OddsPath to evidence strength for PS3/BS3 (Brnich et al. 2019). */
export function oddsPathStrength(o: number): { code: 'PS3' | 'BS3'; strength: Strength } | null {
  if (o > 350) return { code: 'PS3', strength: 'VeryStrong' }
  if (o > 18.7) return { code: 'PS3', strength: 'Strong' }
  if (o > 4.3) return { code: 'PS3', strength: 'Moderate' }
  if (o > 2.1) return { code: 'PS3', strength: 'Supporting' }
  if (o < 0.053) return { code: 'BS3', strength: 'Strong' }
  if (o < 0.23) return { code: 'BS3', strength: 'Moderate' }
  if (o < 0.48) return { code: 'BS3', strength: 'Supporting' }
  return null
}

const SPLICE_TYPES = ['splice_donor', 'splice_acceptor', 'splice_region', 'intronic', 'synonymous']

/** Engine suggestions; the curator can accept, change strength, or override any of them. */
export function suggest(a: Annotation, ev: Evidence, gene: GeneCfg, splice: SpliceReport | null, partners: PhasePartner[] = []): Suggestion[] {
  const out: Suggestion[] = []
  const add = (code: string, strength: Strength, reason: string) => out.push({ code, met: true, strength, reason })
  const p = pvs1(a, gene)
  if (p.strength) add('PVS1', p.strength, p.path.join(' '))

  const af = ev.popmaxAf
  if (af !== null) {
    if (af >= gene.ba1Min) add('BA1', 'StandAlone', `Popmax AF ${fmtAf(af)} is at or above the BA1 threshold ${fmtAf(gene.ba1Min)}.`)
    else if (af >= gene.bs1Min) add('BS1', 'Strong', `Popmax AF ${fmtAf(af)} is above the BS1 threshold ${fmtAf(gene.bs1Min)}.`)
    else if (af <= gene.pm2Max) add('PM2', 'Supporting', `${af === 0 ? 'Absent from' : `Popmax AF ${fmtAf(af)} in`} the population database (threshold ${fmtAf(gene.pm2Max)}). Applied at Supporting per ClinGen SVI.`)
  }

  // Computational: missense uses REVEL; splice uses SpliceAI, then the local PWM model. No double counting.
  const isNull = !!p.strength
  if (!isNull) {
    const splicey = SPLICE_TYPES.includes(a.consequence) || a.spliceRegion
    let comp: { code: 'PP3' | 'BP4'; strength: Strength; reason: string } | null = null
    if (a.consequence === 'missense' && ev.revel !== null) {
      const r = revelStrength(ev.revel)
      if (r) comp = { ...r, reason: `REVEL ${ev.revel} calibrated to ${r.code}_${r.strength} (Pejaver 2022).` }
    }
    if (ev.spliceAi !== null && ev.spliceAi >= 0.2 && (!comp || comp.code === 'BP4')) comp = { code: 'PP3', strength: 'Supporting', reason: `SpliceAI max delta ${ev.spliceAi} >= 0.2 (ClinGen SVI splicing).` }
    else if (comp?.code === 'BP4' && ev.spliceAi !== null && ev.spliceAi > 0.1) comp = null
    if (!comp && ev.spliceAi === null && splice && splicey && (splice.verdict === 'disrupts' || splice.verdict === 'weakens' || splice.verdict === 'creates'))
      comp = { code: 'PP3', strength: 'Supporting', reason: `Local PWM model: ${splice.summary}` }
    if (!comp && a.consequence !== 'missense' && splicey && ev.spliceAi !== null && ev.spliceAi <= 0.1)
      comp = { code: 'BP4', strength: 'Supporting', reason: `SpliceAI max delta ${ev.spliceAi} <= 0.1.` }
    if (comp) add(comp.code, comp.strength, comp.reason)

    const noSplice = ev.spliceAi !== null ? ev.spliceAi <= 0.1 : splice ? splice.verdict === 'none' : false
    const deepIntronic = a.consequence === 'intronic' && a.intronOffset !== undefined && (a.intronOffset >= 7 || a.intronOffset <= -21)
    if ((a.consequence === 'synonymous' && !a.spliceRegion || deepIntronic) && noSplice) add('BP7', 'Supporting', `${a.consequence === 'synonymous' ? 'Synonymous' : 'Deep intronic'} with no predicted splice impact (${ev.spliceAi !== null ? `SpliceAI ${ev.spliceAi}` : 'local PWM model'}).`)
  }

  if (['inframe_deletion', 'inframe_insertion', 'inframe_delins', 'stop_lost'].includes(a.consequence) && !isNull) add('PM4', 'Moderate', a.consequence === 'stop_lost' ? 'Stop-loss extends the protein.' : 'In-frame length change; confirm it is not in a repeat region (else BP3).')

  // Functional evidence: strongest entry in each direction.
  for (const dir of ['abnormal', 'normal'] as const) {
    const list = ev.functional.filter(f => f.result === dir)
    let bestF: { strength: Strength; text: string } | null = null
    for (const f of list) {
      const s = f.oddsPath !== null ? oddsPathStrength(f.oddsPath) : null
      const strength: Strength = s && s.code === (dir === 'abnormal' ? 'PS3' : 'BS3') ? s.strength : 'Supporting'
      if (!bestF || POINTS[strength] > POINTS[bestF.strength]) bestF = { strength, text: `${f.assay}${f.oddsPath !== null ? ` (OddsPath ${f.oddsPath})` : ' (no OddsPath: Supporting)'}` }
    }
    if (bestF) add(dir === 'abnormal' ? 'PS3' : 'BS3', bestF.strength, `${list.length} functional ${dir === 'abnormal' ? 'damaging' : 'normal'} result${list.length > 1 ? 's' : ''}; strongest: ${bestF.text}.`)
  }

  // Phase with other variants in the same gene.
  for (const q of partners) {
    const pathogenic = q.cls === 'Pathogenic' || q.cls === 'Likely pathogenic'
    if (!pathogenic) continue
    if (q.phase === 'cis') { add('BP2', 'Supporting', `In cis with ${q.cls.toLowerCase()} ${q.label}.`); break }
    if (gene.inheritance === 'AR' && q.phase === 'trans') { add('PM3', 'Moderate', `In trans with ${q.cls.toLowerCase()} ${q.label} in a recessive gene.`); break }
    if (gene.inheritance === 'AR' && q.phase === 'unknown') { add('PM3', 'Supporting', `With ${q.cls.toLowerCase()} ${q.label}, phase not confirmed (ClinGen SVI: Supporting).`); break }
    if (gene.inheritance !== 'AR' && q.phase === 'trans') { add('BP2', 'Supporting', `In trans with ${q.cls.toLowerCase()} ${q.label} in a dominant gene.`); break }
  }
  return out
}

export const fmtAf = (x: number) => x === 0 ? '0' : x >= 0.001 ? String(+x.toPrecision(3)) : x.toExponential(1)

// ---------- Curation state ----------

export interface CriterionState { code: string; met: boolean; strength: Strength; evidence: string; source: 'engine' | 'curator'; editor: string; at: number }

/** Effective criteria: curator decisions win over engine suggestions for the same code. */
export function effective(suggestions: Suggestion[], overrides: Record<string, CriterionState>, at = 0): CriterionState[] {
  const out = new Map<string, CriterionState>()
  for (const s of suggestions) if (!out.has(s.code)) out.set(s.code, { code: s.code, met: s.met, strength: s.strength, evidence: s.reason, source: 'engine', editor: 'engine', at })
  for (const o of Object.values(overrides)) out.set(o.code, o)
  return CRITERIA.map(c => out.get(c.code)).filter((x): x is CriterionState => !!x)
}

export const appliedOf = (list: CriterionState[]): Applied[] => list.filter(c => c.met).map(c => ({ code: c.code, strength: c.strength }))

/** "PM2_Supporting" or plain "PVS1" when used at the default strength. */
export const critLabel = (a: Applied) => a.strength === DEF[a.code]?.strength ? a.code : `${a.code}_${a.strength}`

/** Parses lab-report criteria such as "PVS1, PM2_Supporting; PP3_Mod, PS3_P". Unknown tokens are returned separately. */
export function parseCriteria(raw: string): { applied: Applied[]; unknown: string[] } {
  const applied: Applied[] = [], unknown: string[] = []
  for (const tok of raw.split(/[,;|\s]+/).map(s => s.trim()).filter(Boolean)) {
    const m = /^(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])(?:[_-]?(VeryStrong|Very_Strong|VStr|VS|Strong|Str|S|Moderate|Mod|M|Supporting|Sup|P))?$/i.exec(tok)
    if (!m) { unknown.push(tok); continue }
    const code = m[1].toUpperCase()
    const s = (m[2] ?? '').toLowerCase().replace('_', '')
    const strength: Strength = !s ? DEF[code]?.strength ?? 'Supporting' : s.startsWith('v') ? 'VeryStrong' : s.startsWith('s') && s !== 'sup' && s !== 'supporting' ? 'Strong' : s.startsWith('m') ? 'Moderate' : 'Supporting'
    applied.push({ code, strength })
  }
  return { applied, unknown }
}

/** Parses a classification label from a lab report. */
export function parseClass(raw: string): ClassName | null {
  const s = raw.trim().toLowerCase().replace(/[_-]+/g, ' ')
  if (/^(p|path|pathogenic)$/.test(s)) return 'Pathogenic'
  if (/^(lp|likely path(ogenic)?)$/.test(s)) return 'Likely pathogenic'
  if (/^(vus|vous|us|uncertain( significance)?|variant of uncertain significance)$/.test(s)) return 'VUS'
  if (/^(lb|likely benign)$/.test(s)) return 'Likely benign'
  if (/^(b|benign)$/.test(s)) return 'Benign'
  return null
}
