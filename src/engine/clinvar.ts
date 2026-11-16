// Per-submitter classifications from ClinVar (NCBI E-utilities, public domain data): the real lab-versus-lab picture for a variant.
import { parseClass, parseCriteria, type Applied } from './acmg.ts'
import type { LabRow } from './io.ts'

export const EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'

export interface ClinVarSubmission { submitter: string; scv: string; classification: string; reviewStatus: string; date?: string; method?: string; criteria: Applied[] }
export interface ClinVarRecord { variationId: string; accession: string; name: string; gene?: string; updated?: string; aggregate?: string; aggregateReview?: string; submissions: ClinVarSubmission[] }

const unesc = (s: string) => s.replace(/&gt;/g, '>').replace(/&lt;/g, '<').replace(/&quot;/g, '"').replace(/&apos;/g, "'").replace(/&#(\d+);/g, (_, n) => String.fromCodePoint(+n)).replace(/&amp;/g, '&')
const attr = (tag: string, name: string) => { const m = new RegExp(`\\b${name}="([^"]*)"`).exec(tag); return m ? unesc(m[1]) : undefined }
const text = (xml: string, tag: string) => { const m = new RegExp(`<${tag}(?:\\s[^>]*)?>([^<]*)</${tag}>`).exec(xml); return m ? unesc(m[1].trim()) : undefined }

const CODE = /\b(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])(?:_(Very[ _]?Strong|Stand[ _-]?alone|Strong|Moderate|Supporting|Supp|Sup|Mod|Str|VS|[MPS]))?\b/gi

/**
 * ACMG codes a submitter wrote in their free-text comment. Codes they say were not met are skipped.
 * This is a reading aid, not a structured field: ClinVar does not carry criteria separately.
 */
export function criteriaFromComment(comment: string): Applied[] {
  const found: string[] = []
  for (const m of comment.matchAll(CODE)) {
    const before = comment.slice(Math.max(0, m.index - 28), m.index).toLowerCase()
    const after = comment.slice(m.index + m[0].length, m.index + m[0].length + 32).toLowerCase()
    if (/\b(not|no|without|cannot|unable to)\b[^.;]*$/.test(before) || /^[\s:,(]*(is |was |were |are )?(not |n\/a|cannot |could not )/.test(after)) continue
    found.push(m[0])
  }
  return parseCriteria(found.join(',')).applied
}

/**
 * Reads the VCV XML that efetch returns (rettype=vcv). The format is machine-written and regular, so the few fields
 * needed are taken with patterns, which also keeps this usable outside a browser.
 */
export function parseVcvXml(xml: string): ClinVarRecord | null {
  const head = /<VariationArchive\b[^>]*>/.exec(xml)?.[0]
  if (!head) return null
  const agg = /<Classifications>[\s\S]*?<GermlineClassification\b[^>]*>([\s\S]*?)<\/GermlineClassification>/.exec(xml)?.[1] ?? ''
  const rec: ClinVarRecord = {
    variationId: attr(head, 'VariationID') ?? '', accession: attr(head, 'Accession') ?? '', name: attr(head, 'VariationName') ?? '',
    gene: attr(/<Gene\b[^>]*>/.exec(xml)?.[0] ?? '', 'Symbol'), updated: attr(head, 'DateLastUpdated'),
    aggregate: text(agg, 'Description'), aggregateReview: text(agg, 'ReviewStatus'), submissions: [],
  }
  for (const m of xml.matchAll(/<ClinicalAssertion\b[\s\S]*?<\/ClinicalAssertion>/g)) {
    const b = m[0]
    const acc = /<ClinVarAccession\b[^>]*>/.exec(b)?.[0] ?? ''
    const cls = /<Classification\b[^>]*>[\s\S]*?<\/Classification>/.exec(b)?.[0] ?? ''
    const classification = text(cls, 'GermlineClassification')
    const submitter = attr(acc, 'SubmitterName')
    if (!classification || !submitter) continue
    const comments = [...cls.matchAll(/<Comment\b[^>]*>([^<]*)<\/Comment>/g)].map(c => unesc(c[1])).join(' ')
    rec.submissions.push({
      submitter, scv: attr(acc, 'Accession') ?? '', classification, reviewStatus: text(cls, 'ReviewStatus') ?? '',
      date: attr(/<Classification\b[^>]*>/.exec(cls)![0], 'DateLastEvaluated') ?? attr(acc, 'DateUpdated'),
      method: /<Attribute Type="AssertionMethod">([^<]*)</.exec(b)?.[1], criteria: criteriaFromComment(comments),
    })
  }
  return rec
}

/** Lab rows for the discrepancy view. Submissions without a germline pathogenicity call ("not provided", "risk factor") are counted, not imported. */
export function clinvarLabRows(rec: ClinVarRecord, variant: string, idGen: () => string): { rows: LabRow[]; skipped: number } {
  const rows: LabRow[] = []
  let skipped = 0
  rec.submissions.forEach((s, i) => {
    const cls = parseClass(s.classification)
    if (!cls) { skipped++; return }
    rows.push({ id: idGen(), variant, ...(rec.gene ? { gene: rec.gene } : {}), lab: s.submitter, cls, criteria: s.criteria, ...(s.date ? { date: s.date } : {}), line: i + 1, file: `ClinVar ${rec.accession} ${s.scv}` })
  })
  return { rows, skipped }
}

/** Search terms tried in order: the full transcript description, then gene plus c. description. */
export function clinvarTerms(v: { transcript?: string; gene?: string; hgvsC: string }): string[] {
  const c = v.hgvsC.replace(/^.*?(c\.)/, '$1')
  const out: string[] = []
  if (v.transcript) out.push(`"${v.transcript.replace(/\.\d+$/, '')}" AND "${c}"`)
  if (v.gene) out.push(`${v.gene}[gene] AND "${c}"`)
  return out
}

/** True when a ClinVar title such as "NM_000059.4(BRCA2):c.9976A>T (p.Lys3326Ter)" describes this c. change. */
export const titleMatches = (title: string, hgvsC: string) => title.replace(/\s+\(p\..*$/, '').endsWith(':' + hgvsC.replace(/^.*?(c\.)/, '$1'))

const sleep = (ms: number) => new Promise(r => setTimeout(r, ms))
let last = 0
/** E-utilities allow three requests a second without a key; calls are spaced to stay under that. */
async function eutils(path: string): Promise<Response> {
  const wait = last + 400 - Date.now()
  if (wait > 0) await sleep(wait)
  last = Date.now()
  let res: Response
  try { res = await fetch(`${EUTILS}/${path}`) } catch { throw new Error('Could not reach NCBI E-utilities. Check the connection and try again.') }
  if (res.status === 429) { await sleep(1200); last = Date.now(); res = await fetch(`${EUTILS}/${path}`) }
  if (!res.ok) throw new Error(`NCBI returned ${res.status}.`)
  return res
}

/** Looks a variant up by HGVS and returns its ClinVar record, or null when ClinVar has no exact match. Only the HGVS description is sent. */
export async function fetchClinVar(v: { transcript?: string; gene?: string; hgvsC: string }): Promise<ClinVarRecord | null> {
  for (const term of clinvarTerms(v)) {
    const ids: string[] = (await (await eutils(`esearch.fcgi?db=clinvar&retmode=json&retmax=20&term=${encodeURIComponent(term)}`)).json())?.esearchresult?.idlist ?? []
    if (!ids.length) continue
    const sum = (await (await eutils(`esummary.fcgi?db=clinvar&retmode=json&id=${ids.join(',')}`)).json())?.result ?? {}
    const id = ids.find(i => titleMatches(String(sum[i]?.title ?? ''), v.hgvsC))
    if (!id) continue
    return parseVcvXml(await (await eutils(`efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=${id}`)).text())
  }
  return null
}
