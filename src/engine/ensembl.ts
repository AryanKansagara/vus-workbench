// Real gene models from the public Ensembl REST API (GRCh38): the MANE Select transcript and the reference sequence under it.
import type { RefSeq, Transcript } from './transcript.ts'

export const ENSEMBL = 'https://rest.ensembl.org'
/** Bases fetched either side of the transcript, enough for upstream variants and splice-site windows. */
export const FLANK = 300
/** Ensembl serves at most 10 Mb per region request; the largest human gene (DMD) is about 2.2 Mb. */
const MAX_REGION = 5_000_000

interface EnsExon { start: number; end: number }
interface EnsTranscript {
  id: string; version?: number; biotype?: string; is_canonical?: number; strand: number; seq_region_name: string
  Exon?: EnsExon[]; Translation?: { start: number; end: number; length?: number }
  MANE?: { type?: string; refseq_match?: string }[]
}
export interface EnsGene { display_name?: string; id?: string; assembly_name?: string; Transcript?: EnsTranscript[]; error?: string }

/**
 * Converts an Ensembl gene lookup (expand=1, mane=1) to a transcript model. MANE Select is used and named by its
 * RefSeq match, since the two are identical by definition and reports quote the NM_ accession; otherwise the Ensembl canonical transcript.
 */
export function transcriptFromEnsembl(gene: EnsGene, release?: number): Transcript {
  if (gene.error) throw new Error(gene.error)
  const coding = (gene.Transcript ?? []).filter(t => t.Translation && t.Exon?.length)
  if (!coding.length) throw new Error(`${gene.display_name ?? 'This gene'} has no protein-coding transcript in Ensembl.`)
  const mane = (t: EnsTranscript) => t.MANE?.find(m => m.type === 'MANE_Select')
  const t = coding.find(mane) ?? coding.find(t => t.is_canonical) ?? coding[0]
  const ens = `${t.id}${t.version ? '.' + t.version : ''}`
  const refseq = mane(t)?.refseq_match
  return {
    id: refseq ?? ens,
    gene: gene.display_name ?? gene.id ?? ens,
    chrom: t.seq_region_name,
    strand: t.strand === -1 ? -1 : 1,
    exons: t.Exon!.map(e => [e.start, e.end] as [number, number]).sort((a, b) => a[0] - b[0]),
    cdsStart: t.Translation!.start,
    cdsEnd: t.Translation!.end,
    source: `Ensembl${release ? ' ' + release : ''} ${refseq ? 'MANE Select' : 'canonical'} ${ens} (${gene.assembly_name ?? 'GRCh38'})`,
  }
}

async function get(path: string, type: 'json' | 'text') {
  let res: Response
  try { res = await fetch(`${ENSEMBL}${path}`, { headers: { Accept: type === 'json' ? 'application/json' : 'text/plain' } }) } catch { throw new Error('Could not reach rest.ensembl.org. Check the connection and try again.') }
  if (res.status === 429) throw new Error('Ensembl is rate-limiting requests. Wait a few seconds and try again.')
  if (type === 'text') { if (!res.ok) throw new Error(`Ensembl returned ${res.status} for the sequence request.`); return res.text() }
  const body = await res.json().catch(() => null)
  if (!res.ok) throw new Error(res.status === 400 ? 'not-found' : body?.error ?? `Ensembl returned ${res.status}.`)
  return body
}

/** Fetches the model and reference region for one gene symbol. Two requests; nothing about the worklist is sent, only the symbol. */
export async function fetchGeneModel(symbol: string): Promise<{ transcript: Transcript; region: RefSeq; key: string }> {
  const sym = symbol.trim()
  if (!/^[A-Za-z0-9-]{1,30}$/.test(sym)) throw new Error(`"${symbol}" is not a gene symbol.`)
  let gene: EnsGene
  try { gene = await get(`/lookup/symbol/homo_sapiens/${encodeURIComponent(sym)}?expand=1;mane=1`, 'json') } catch (e) {
    throw (e as Error).message === 'not-found' ? new Error(`Ensembl has no human gene called ${sym}. Use the current HGNC symbol.`) : e
  }
  const release = await get('/info/data', 'json').then(d => d?.releases?.[0] as number | undefined).catch(() => undefined)
  const transcript = transcriptFromEnsembl(gene, release)
  const start = Math.max(1, transcript.exons[0][0] - FLANK), end = transcript.exons[transcript.exons.length - 1][1] + FLANK
  if (end - start > MAX_REGION) throw new Error(`${sym} spans ${(end - start).toLocaleString()} bases, too large to load here.`)
  const seq = (await get(`/sequence/region/human/${transcript.chrom}:${start}..${end}:1?coord_system_version=GRCh38`, 'text') as string).replace(/\s+/g, '').toUpperCase()
  if (seq.length !== end - start + 1) throw new Error('Ensembl returned a sequence of unexpected length.')
  return { transcript, region: { chrom: transcript.chrom, start, seq }, key: `${transcript.chrom}:${start}-${end}` }
}
