import { CLASS_SHORT, vusTier, type ClassName } from './engine/acmg.ts'
import type { Annotation, Transcript } from './engine/transcript.ts'
import { cdsBounds, gToTx, txExons } from './engine/transcript.ts'

export function ClassBadge({ cls, points, small }: { cls: ClassName; points?: number; small?: boolean }) {
  const tier = cls === 'VUS' && points !== undefined ? vusTier(points) : null
  return (
    <span className={`cls cls-${CLASS_SHORT[cls].toLowerCase()}${small ? ' sm' : ''}`} title={`${cls}${points !== undefined ? `, ${points} points` : ''}${tier ? ` (${tier} VUS)` : ''}`}>
      {small ? CLASS_SHORT[cls] : cls}{tier && !small ? <em> · {tier}</em> : null}
    </span>
  )
}

// ---------- Bayesian points meter ----------

const LO = -12, HI = 14
const BANDS: { cls: ClassName; from: number; to: number }[] = [
  { cls: 'Benign', from: LO, to: -6.5 },
  { cls: 'Likely benign', from: -6.5, to: -0.5 },
  { cls: 'VUS', from: -0.5, to: 5.5 },
  { cls: 'Likely pathogenic', from: 5.5, to: 9.5 },
  { cls: 'Pathogenic', from: 9.5, to: HI },
]

export function Meter({ points, cls, compact }: { points: number; cls: ClassName; compact?: boolean }) {
  const x = (p: number) => ((Math.max(LO, Math.min(HI, p)) - LO) / (HI - LO)) * 100
  return (
    <div className={`meter${compact ? ' compact' : ''}`} role="meter" aria-valuemin={LO} aria-valuemax={HI} aria-valuenow={points} aria-label={`Bayesian points ${points}: ${cls}`}>
      <div className="meter-track">
        {BANDS.map(b => (
          <div key={b.cls} className={`meter-band band-${CLASS_SHORT[b.cls].toLowerCase()}${b.cls === cls ? ' on' : ''}`} style={{ left: `${x(b.from)}%`, width: `${x(b.to) - x(b.from)}%` }}>
            {!compact && <span>{CLASS_SHORT[b.cls]}</span>}
          </div>
        ))}
        <div className="meter-needle" style={{ left: `${x(points)}%` }}>
          <b>{points > 0 ? '+' : ''}{points}</b>
        </div>
      </div>
      {!compact && (
        <div className="meter-scale" aria-hidden="true">
          {[-10, -7, -1, 0, 6, 10].map(t => <span key={t} style={{ left: `${x(t)}%` }}>{t}</span>)}
        </div>
      )}
    </div>
  )
}

// ---------- Exon / intron track ----------

interface TrackProps { tx: Transcript | null; ann: Annotation | null; gpos?: number | null; height?: number; animate?: boolean }

/** Transcript drawn 5' to 3' with introns compressed; the variant is placed by its c. position. */
export function Track({ tx, ann, gpos, height = 132, animate }: TrackProps) {
  const W = 1000, PAD = 24, Y = height / 2
  if (!tx) return <SchematicTrack ann={ann} height={height} />
  const ex = txExons(tx)
  const lens = ex.map(e => e[1] - e[0] + 1)
  const INTRON = 46
  const total = lens.reduce((a, b) => a + b, 0)
  const scale = (W - PAD * 2 - INTRON * (ex.length - 1)) / total
  const cds = cdsBounds(tx)
  let x = PAD, acc = 0
  const boxes = ex.map((e, i) => {
    const w = lens[i] * scale
    const b = { i, x, w, txStart: acc + 1, txEnd: acc + lens[i], g: e }
    x += w + INTRON
    acc += lens[i]
    return b
  })
  const txX = (t: number) => { const b = boxes.find(b => t >= b.txStart && t <= b.txEnd) ?? boxes[boxes.length - 1]; return b.x + (t - b.txStart + 0.5) * scale }
  // Variant x: exonic by transcript position; intronic by offset within the compressed intron.
  let vx: number | null = null
  if (gpos != null) {
    const t = gToTx(tx, gpos)
    if (t !== null) vx = txX(t)
  }
  if (vx === null && ann?.intron !== undefined && ann.intronOffset !== undefined) {
    const b = boxes[ann.intron - 1]
    const frac = Math.min(0.45, Math.log10(Math.abs(ann.intronOffset) + 1) / 5 + 0.04)
    vx = ann.intronOffset > 0 ? b.x + b.w + INTRON * frac : boxes[ann.intron].x - INTRON * frac
  }
  const exonH = 30, utrH = 16
  const highlight = ann?.exon ?? (ann?.intron !== undefined ? (ann.intronOffset! > 0 ? ann.intron : ann.intron + 1) : null)
  return (
    <svg className={`track${animate ? ' animate' : ''}`} viewBox={`0 0 ${W} ${height}`} role="img" aria-label={`${tx.gene} ${tx.id}: ${ex.length} exons${ann ? `, variant ${ann.hgvsC}` : ''}`}>
      <line x1={PAD} x2={W - PAD} y1={Y} y2={Y} className="tr-intron" />
      {boxes.slice(0, -1).map(b => {
        const x1 = b.x + b.w, x2 = x1 + INTRON, mid = (x1 + x2) / 2
        return <path key={b.i} d={`M${x1} ${Y} L${mid} ${Y - 12} L${x2} ${Y}`} className="tr-splice" />
      })}
      {boxes.map(b => {
        // Split each exon into UTR and CDS parts.
        const a = Math.max(b.txStart, cds.start), z = Math.min(b.txEnd, cds.end)
        const parts: { x: number; w: number; coding: boolean }[] = []
        if (b.txStart < cds.start) parts.push({ x: b.x, w: (Math.min(b.txEnd, cds.start - 1) - b.txStart + 1) * scale, coding: false })
        if (a <= z) parts.push({ x: b.x + (a - b.txStart) * scale, w: (z - a + 1) * scale, coding: true })
        if (b.txEnd > cds.end) { const s = Math.max(b.txStart, cds.end + 1); parts.push({ x: b.x + (s - b.txStart) * scale, w: (b.txEnd - s + 1) * scale, coding: false }) }
        return (
          <g key={b.i} className={`tr-exon${highlight === b.i + 1 ? ' hi' : ''}`}>
            {parts.map((p, k) => <rect key={k} x={p.x} y={Y - (p.coding ? exonH : utrH) / 2} width={Math.max(p.w, 1)} height={p.coding ? exonH : utrH} rx={p.coding ? 5 : 3} className={p.coding ? 'cds' : 'utr'} />)}
            <text x={b.x + b.w / 2} y={Y + exonH / 2 + 16} className="tr-label">{b.i + 1}</text>
          </g>
        )
      })}
      <text x={PAD} y={16} className="tr-end">5′ {tx.gene} · {tx.id} · {tx.strand === 1 ? '+' : '−'} strand</text>
      <text x={W - PAD} y={16} className="tr-end" textAnchor="end">3′</text>
      {vx !== null && ann && (
        <g className="tr-var" style={{ transform: `translateX(${vx}px)` }}>
          <line x1={0} x2={0} y1={26} y2={height - 8} />
          <circle cx={0} cy={Y} r={7} />
          <text x={0} y={height - 2} textAnchor={vx > W * 0.8 ? 'end' : vx < W * 0.2 ? 'start' : 'middle'}>{ann.hgvsC.replace(/^.*:/, '')}</text>
        </g>
      )}
    </svg>
  )
}

/** Without a model: a generic exon-intron-exon schematic, placing the variant by the sign of its offset. */
function SchematicTrack({ ann, height }: { ann: Annotation | null; height: number }) {
  const Y = height / 2
  const off = ann?.intronOffset
  const vx = off === undefined ? 230 : off > 0 ? 400 + Math.min(180, Math.log10(Math.abs(off) + 1) * 90) : 600 - Math.min(180, Math.log10(Math.abs(off) + 1) * 90)
  return (
    <svg className="track" viewBox={`0 0 1000 ${height}`} role="img" aria-label="Schematic: no transcript model loaded">
      <line x1={60} x2={940} y1={Y} y2={Y} className="tr-intron" />
      <rect x={140} y={Y - 15} width={260} height={30} rx={5} className="cds schematic" />
      <rect x={600} y={Y - 15} width={260} height={30} rx={5} className="cds schematic" />
      <text x={500} y={18} textAnchor="middle" className="tr-end">Schematic only. Import a GTF subset to draw the real exon structure.</text>
      {ann && <g className="tr-var" style={{ transform: `translateX(${vx}px)` }}><line x1={0} x2={0} y1={26} y2={height - 8} /><circle cx={0} cy={Y} r={7} /><text x={0} y={height - 2} textAnchor="middle">{ann.hgvsC.replace(/^.*:/, '')}</text></g>}
    </svg>
  )
}

/** Sequence logo-ish view of a splice site: reference vs alternate with the changed base marked. */
export function SiteSeq({ refSeq, altSeq, site }: { refSeq: string; altSeq: string; site: 'donor' | 'acceptor' }) {
  const cut = site === 'donor' ? 3 : 15
  const row = (s: string, other: string) => [...s].map((b, i) => <span key={i} className={`nt nt-${b}${b !== other[i] ? ' diff' : ''}${i === cut ? ' cut' : ''}`}>{b}</span>)
  return (
    <div className="siteseq" aria-label={`${site} site, reference ${refSeq}, alternate ${altSeq}`}>
      <div><small>REF</small>{row(refSeq, altSeq)}</div>
      <div><small>ALT</small>{row(altSeq, refSeq)}</div>
    </div>
  )
}
