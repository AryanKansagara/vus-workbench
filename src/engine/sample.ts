// Synthetic gene models and sample worklists. Every sequence, genotype, frequency and lab call here is synthetic or
// illustrative, not patient data. The real-gene sample uses public HGVS descriptions with made-up lab calls.
import { CODON, revcomp, rng } from './seq.ts'
import type { Genome, Transcript } from './transcript.ts'
import { cToG, hgvsCToG } from './transcript.ts'
import { parseC } from './hgvs.ts'
import type { GeneCfg } from './acmg.ts'

interface GeneSpec {
  id: string; gene: string; chrom: string; strand: 1 | -1; seed: number
  exons: number[]; introns: number[]; utr5: number; utr3: number
  codons: Record<number, string>
}

const STOPS = new Set(['TAA', 'TAG', 'TGA'])
const NON_STOP = Object.keys(CODON).filter(c => !STOPS.has(c))

/** Builds a gene in transcript orientation, then places it on its contig (reverse-complemented for minus strand). */
function build(s: GeneSpec): { tx: Transcript; ref: { chrom: string; start: number; seq: string } } {
  const r = rng(s.seed)
  const pick = (alpha: string, n: number) => Array.from({ length: n }, () => alpha[Math.floor(r() * alpha.length)]).join('')
  const cdsLen = s.exons.reduce((a, b) => a + b, 0) - s.utr5 - s.utr3
  if (cdsLen % 3) throw new Error(`${s.gene}: CDS ${cdsLen} not a multiple of 3`)
  const n = cdsLen / 3
  let cds = ''
  for (let i = 1; i <= n; i++) cds += i === 1 ? 'ATG' : i === n ? 'TAA' : s.codons[i] ?? NON_STOP[Math.floor(r() * NON_STOP.length)]
  const mrna = pick('ACGT', s.utr5) + cds + pick('ACGT', s.utr3)
  const FLANK = 300
  let seq = pick('ACGT', FLANK)
  const exTx: [number, number][] = []
  let off = 0
  s.exons.forEach((len, i) => {
    exTx.push([seq.length + 1, seq.length + len])
    seq += mrna.slice(off, off + len)
    off += len
    if (i < s.introns.length) {
      const L = s.introns[i]
      // Consensus donor GTAAGT, AT-rich body, branch-point-like A, polypyrimidine tract and AG acceptor.
      seq += 'GTAAGT' + pick('AATTTGCA', L - 6 - 22) + 'CTAAC' + 'TTTCTCTTCCCTTTTTCAG'.slice(-17)
    }
  })
  seq += pick('ACGT', FLANK)
  const START = 100001
  const L = seq.length
  const toG = (p: number) => s.strand === 1 ? START + p - 1 : START + L - p
  const exons = exTx.map(([a, b]) => { const x = toG(a), y = toG(b); return [Math.min(x, y), Math.max(x, y)] as [number, number] }).sort((a, b) => a[0] - b[0])
  const cdsTxStart = exTx[0][0] + s.utr5
  const lastEx = exTx[exTx.length - 1]
  const cdsTxEnd = lastEx[1] - s.utr3
  const cg = [toG(cdsTxStart), toG(cdsTxEnd)]
  const tx: Transcript = { id: s.id, gene: s.gene, chrom: s.chrom, strand: s.strand, exons, cdsStart: Math.min(...cg), cdsEnd: Math.max(...cg), source: 'Bundled synthetic model', synthetic: true }
  return { tx, ref: { chrom: s.chrom, start: START, seq: s.strand === 1 ? seq : revcomp(seq) } }
}

const SPECS: GeneSpec[] = [
  {
    id: 'SYN_000101.1', gene: 'SYNG1', chrom: 'chrS1', strand: 1, seed: 101,
    exons: [120, 96, 150, 87, 111, 300], introns: [400, 350, 500, 300, 450], utr5: 60, utr3: 201,
    codons: { 35: 'CAG', 44: 'GAC', 49: 'CTC', 50: 'AAA', 51: 'AAG', 52: 'CTG', 70: 'CGG', 90: 'GCC', 100: 'GCG', 101: 'GCA', 102: 'CAG', 120: 'GAT', 134: 'TCA', 135: 'CAA', 186: 'GAG', 187: 'CTG' },
  },
  {
    id: 'SYN_000202.1', gene: 'SYNG2', chrom: 'chrS2', strand: -1, seed: 202,
    exons: [90, 112, 135, 160], introns: [300, 420, 380], utr5: 40, utr3: 79,
    codons: { 34: 'GGC', 69: 'GAT', 70: 'CCA', 84: 'CGC' },
  },
]

const BUILT = SPECS.map(build)
export const BUNDLED_TRANSCRIPTS: Transcript[] = BUILT.map(b => b.tx)
export const BUNDLED_GENOME: Genome = Object.fromEntries(BUILT.map(b => [b.ref.chrom, b.ref]))

export const BUNDLED_GENES: Record<string, GeneCfg> = {
  SYNG1: { inheritance: 'AD', lofMechanism: true, pm2Max: 0.00002, bs1Min: 0.005, ba1Min: 0.05 },
  SYNG2: { inheritance: 'AR', lofMechanism: true, pm2Max: 0.0007, bs1Min: 0.01, ba1Min: 0.05 },
  // Protein lengths from the public RefSeq records (MANE Select), used to size truncations when no model is loaded.
  BRCA1: { inheritance: 'AD', lofMechanism: true, pm2Max: 0.00002, bs1Min: 0.001, ba1Min: 0.001, proteinLength: 1863 },
  BRCA2: { inheritance: 'AD', lofMechanism: true, pm2Max: 0.00002, bs1Min: 0.001, ba1Min: 0.001, proteinLength: 3418 },
  CFTR: { inheritance: 'AR', lofMechanism: true, pm2Max: 0.005, bs1Min: 0.02, ba1Min: 0.05, proteinLength: 1480 },
}

/** Genomic VCF line for a c. description on a bundled synthetic transcript (used to write the sample VCF). */
function vcfFor(c: string, extra: { info?: string; gt?: string; alt2?: string } = {}): string {
  const h = parseC(c)!
  const t = BUNDLED_TRANSCRIPTS.find(t => t.id === h.transcript)!
  // "N" as the reference base means "whatever the synthetic reference has here" (random intron/UTR sequence).
  if (h.kind === 'sub' && h.ref === 'N') {
    const g = cToG(t, h.start), r = BUNDLED_GENOME[t.chrom]
    const b = r.seq[g - r.start]
    h.ref = t.strand === 1 ? b : revcomp(b)
    if (h.alt === h.ref) h.alt = { A: 'G', C: 'T', G: 'A', T: 'C' }[h.ref]!
  }
  const v = hgvsCToG(t, BUNDLED_GENOME, h)
  if (!v) throw new Error(`sample: ${c} does not match the synthetic reference`)
  return [v.chrom, v.pos, '.', v.ref, v.alt + (extra.alt2 ? ',' + extra.alt2 : ''), 60, 'PASS', extra.info ?? '.', 'GT:PS', `${extra.gt ?? '0/1'}:${extra.gt?.includes('|') ? 1 : '.'}`].join('\t')
}

/** A deletion written one base to the right of where it normalises to, as some callers emit it. */
function unnormalisedDel(): string {
  const t = BUNDLED_TRANSCRIPTS[0]
  const g = cToG(t, { pos: 150, offset: 0 })
  const ref = BUNDLED_GENOME.chrS1.seq.slice(g - BUNDLED_GENOME.chrS1.start, g - BUNDLED_GENOME.chrS1.start + 2)
  return ['chrS1', g, '.', ref, ref[0], 55, 'PASS', 'gnomAD_AF_popmax=0;SpliceAI=A|SYNG1|0.00|0.00|0.01|0.00|0|0|0|0', 'GT:PS', '0/1:.'].join('\t')
}

const S1 = 'SYN_000101.1:'
const S2 = 'SYN_000202.1:'

export function sampleVcf(): string {
  const head = [
    '##fileformat=VCFv4.2',
    '##source=VUS Workbench synthetic sample (not patient data)',
    '##reference=synthetic contigs chrS1, chrS2',
    '##contig=<ID=chrS1,length=4500>',
    '##contig=<ID=chrS2,length=2400>',
    '##INFO=<ID=gnomAD_AF_popmax,Number=A,Type=Float,Description="Illustrative popmax allele frequency">',
    '##INFO=<ID=REVEL,Number=A,Type=Float,Description="Illustrative REVEL score">',
    '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYN-PROBAND-01',
  ]
  const rows = [
    vcfFor(`${S1}c.103C>T`, { info: 'gnomAD_AF_popmax=0' }),
    vcfFor(`${S1}c.208C>T`, { info: 'gnomAD_AF_popmax=0.000004,0;REVEL=0.81,0.52', alt2: 'G', gt: '1/2' }),
    vcfFor(`${S1}c.270C>T`, { info: 'gnomAD_AF_popmax=0.0003;SpliceAI=T|SYNG1|0.00|0.00|0.02|0.01|0|0|0|0' }),
    vcfFor(`${S1}c.306+1G>A`, { info: 'gnomAD_AF_popmax=0' }),
    vcfFor(`${S1}c.306+5G>A`, { info: 'gnomAD_AF_popmax=0.000008' }),
    vcfFor(`${S1}c.358G>A`, { info: 'gnomAD_AF_popmax=0.0071;REVEL=0.12' }),
    vcfFor(`${S1}c.130_132del`, { info: 'gnomAD_AF_popmax=0' }),
    unnormalisedDel(),
    vcfFor(`${S1}c.560dup`, { info: 'gnomAD_AF_popmax=0' }),
    vcfFor(`${S1}c.*35N>A`, { info: 'gnomAD_AF_popmax=0.118' }),
    vcfFor(`${S1}c.393+60N>C`, { info: 'gnomAD_AF_popmax=0.0021;SpliceAI=C|SYNG1|0.00|0.00|0.00|0.00|0|0|0|0' }),
    vcfFor(`${S2}c.100G>A`, { info: 'gnomAD_AF_popmax=0.00011;REVEL=0.70', gt: '0|1' }),
    vcfFor(`${S2}c.205del`, { info: 'gnomAD_AF_popmax=0.00002', gt: '1|0' }),
    vcfFor(`${S2}c.250C>T`, { info: 'gnomAD_AF_popmax=0.00005;REVEL=0.41', gt: '1|0' }),
  ]
  // VCF is sorted by contig and position.
  const sorted = rows.sort((a, b) => { const [ca, pa] = a.split('\t'), [cb, pb] = b.split('\t'); return ca.localeCompare(cb) || +pa - +pb })
  return [...head, ...sorted, ''].join('\n')
}

/** Some reference bases for the docs and tests. */
export const synthRef = (c: string) => { const h = parseC(c)!; const t = BUNDLED_TRANSCRIPTS.find(t => t.id === h.transcript)!; return hgvsCToG(t, BUNDLED_GENOME, h) }

export const REAL_VARIANTS_CSV = `gene,transcript,hgvs_c,hgvs_p,popmax_af,revel,spliceai,sample
BRCA2,NM_000059.4,c.68-7T>A,,0.021,,0.03,SYN-CASE-07
BRCA2,NM_000059.4,c.7522G>A,p.(Gly2508Ser),0.000012,0.68,0.01,SYN-CASE-07
BRCA2,NM_000059.4,c.9976A>T,p.(Lys3326Ter),0.0065,,0.00,SYN-CASE-11
BRCA1,NM_007294.4,c.5266dup,p.(Gln1756ProfsTer74),0.00009,,,SYN-CASE-03
BRCA1,NM_007294.4,c.5096G>A,p.(Arg1699Gln),0.00003,0.94,0.02,SYN-CASE-05
CFTR,NM_000492.4,c.1521_1523del,p.(Phe508del),0.012,,,SYN-CASE-09
`

export const LAB_CSV = `variant,gene,lab,classification,criteria,date_evaluated
NM_000059.4:c.68-7T>A,BRCA2,Synthetic Lab A,VUS,"PM2_Supporting, PP3",2023-04-12
NM_000059.4:c.68-7T>A,BRCA2,Synthetic Lab B,Benign,"BS1, BP4, BS3_Supporting",2025-02-03
NM_000059.4:c.68-7T>A,BRCA2,Synthetic Lab C,Likely benign,"BS1, BP4",2024-09-18
NM_000059.4:c.7522G>A,BRCA2,Synthetic Lab A,Likely pathogenic,"PM2_Supporting, PP3_Moderate, PS3_Moderate",2024-06-01
NM_000059.4:c.7522G>A,BRCA2,Synthetic Lab B,VUS,"PM2_Supporting, PP3",2022-11-20
NM_000059.4:c.9976A>T,BRCA2,Synthetic Lab A,Pathogenic,"PVS1, PM2",2019-08-30
NM_000059.4:c.9976A>T,BRCA2,Synthetic Lab C,Benign,"BS1, BS3, BP5",2025-05-14
NM_007294.4:c.5266dup,BRCA1,Synthetic Lab A,Pathogenic,"PVS1, PS4, PM2_Supporting",2024-01-09
NM_007294.4:c.5266dup,BRCA1,Synthetic Lab B,Pathogenic,"PVS1, PS3, PM2_Supporting",2023-07-22
NM_007294.4:c.5096G>A,BRCA1,Synthetic Lab B,Likely pathogenic,"PS3, PM2_Supporting, PP3_Strong",2024-10-02
NM_007294.4:c.5096G>A,BRCA1,Synthetic Lab C,VUS,"PM2_Supporting, PP3_Strong, BS3_Supporting",2025-03-27
NM_000492.4:c.1521_1523del,CFTR,Synthetic Lab A,Pathogenic,"PS3, PM3_VeryStrong, PM4",2021-12-15
NM_000492.4:c.1521_1523del,CFTR,Synthetic Lab B,Likely pathogenic,"PS3, PM4",2020-05-05
`

export const SAMPLE_GTF = BUNDLED_TRANSCRIPTS.flatMap(t => {
  const row = (type: string, a: number, b: number) => [t.chrom, 'synthetic', type, a, b, '.', t.strand === 1 ? '+' : '-', '.', `gene_id "${t.gene}"; transcript_id "${t.id}"; gene_name "${t.gene}";`].join('\t')
  return [
    ...t.exons.map(e => row('exon', e[0], e[1])),
    ...t.exons.filter(e => e[1] >= t.cdsStart && e[0] <= t.cdsEnd).map(e => row('CDS', Math.max(e[0], t.cdsStart), Math.min(e[1], t.cdsEnd))),
  ]
}).join('\n') + '\n'

export interface SampleDef { id: string; name: string; blurb: string }
export const SAMPLES: SampleDef[] = [
  { id: 'synthetic-panel', name: 'Synthetic panel: proband SYN-PROBAND-01', blurb: 'VCF on two synthetic genes: nonsense, canonical and near-splice, in-frame, unnormalised indel, multi-allelic site, and two VUS in one recessive gene with phase.' },
  { id: 'lab-discrepancy', name: 'BRCA1/2 and CFTR: lab discrepancies', blurb: 'Public HGVS descriptions with synthetic lab calls: P vs B on a last-exon stop, a hot VUS, and labs whose criteria do not add up to their call.' },
]
