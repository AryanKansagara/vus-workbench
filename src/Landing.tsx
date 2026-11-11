import { useEffect, useMemo, useState } from 'react'
import { combine, critLabel, STRENGTH_SHORT, type Applied, type Strength } from './engine/acmg.ts'
import { parseVcf } from './engine/io.ts'
import { analyse } from './engine/pipeline.ts'
import { reconcile } from './engine/discrepancy.ts'
import { autoMap, LAB_FIELDS, parseCsv, validateLabCsv } from './engine/io.ts'
import { BUNDLED_GENOME, BUNDLED_TRANSCRIPTS, LAB_CSV, sampleVcf } from './engine/sample.ts'
import { ClassBadge, Meter, SiteSeq, Track } from './viz.tsx'
import { Link, SyntheticBadge } from './ui.tsx'

// The hero runs the real engine on the bundled synthetic SYNG1 c.306+5G>A.
let k = 0
const HERO = (() => {
  const v = parseVcf(sampleVcf(), 'hero.vcf', BUNDLED_GENOME, () => `h${k++}`).variants
  const all = v.map(r => analyse(r, BUNDLED_TRANSCRIPTS, BUNDLED_GENOME))
  return all.find(a => a.ann?.hgvsC.endsWith('c.306+5G>A'))!
})()

interface Chip { code: string; strength: Strength; why: string }
const CHIPS: Chip[] = [
  { code: 'PM2', strength: 'Supporting', why: 'Absent from the population database' },
  { code: 'PP3', strength: 'Supporting', why: `Donor ${HERO.splice?.canonical?.refScore} → ${HERO.splice?.canonical?.altScore} bits (local PWM)` },
  { code: 'PS3', strength: 'Strong', why: 'Minigene: exon 3 skipping, OddsPath 19.2' },
  { code: 'PP4', strength: 'Supporting', why: 'Phenotype specific for the gene' },
  { code: 'PP1', strength: 'Strong', why: 'Co-segregates across 7 meioses' },
  { code: 'BS3', strength: 'Strong', why: 'RNA from blood shows normal splicing' },
  { code: 'BP4', strength: 'Supporting', why: 'SpliceAI 0.04 (conflicts with PP3)' },
]
const SCRIPT = [['PM2', 'PP3'], ['PM2', 'PP3', 'PS3'], ['PM2', 'PP3', 'PS3', 'PP4'], ['PM2', 'PP3', 'PS3', 'PP4', 'PP1'], ['PM2', 'PP3', 'PS3', 'PP4', 'PP1', 'BS3'], ['PM2', 'BS3', 'BP4']]

function Hero() {
  const [on, setOn] = useState<string[]>(SCRIPT[0])
  const [auto, setAuto] = useState(() => !matchMedia('(prefers-reduced-motion: reduce)').matches)
  useEffect(() => {
    if (!auto) return
    let i = 0
    const t = setInterval(() => { i = (i + 1) % SCRIPT.length; setOn(SCRIPT[i]) }, 2600)
    return () => clearInterval(t)
  }, [auto])
  const applied: Applied[] = CHIPS.filter(c => on.includes(c.code)).map(c => ({ code: c.code, strength: c.strength }))
  const r = combine(applied)
  const c = HERO.splice?.canonical
  const toggle = (code: string) => { setAuto(false); setOn(o => o.includes(code) ? o.filter(x => x !== code) : [...o, code]) }
  return (
    <div className="hero-card" aria-label="Live curation of a synthetic splice-region variant">
      <div className="hero-card-head">
        <div>
          <div className="eyebrow light">Live engine · synthetic variant <SyntheticBadge /></div>
          <strong className="mono hero-var">{HERO.label}</strong>
        </div>
        <ClassBadge cls={r.cls} points={r.points} />
      </div>
      <Track tx={HERO.tx} ann={HERO.ann} gpos={HERO.g?.pos} animate />
      {c && (
        <div className="hero-splice">
          <SiteSeq refSeq={c.refSeq} altSeq={c.altSeq} site={c.site} />
          <div className="hero-splice-score"><span>Donor strength</span><b>{c.refScore}</b><i aria-hidden="true">→</i><b className="down">{c.altScore}</b><small>bits ({c.pct}%)</small></div>
        </div>
      )}
      <div className="chips" role="group" aria-label="Toggle evidence criteria">
        {CHIPS.map(ch => {
          const active = on.includes(ch.code)
          const excluded = active && r.excluded.some(e => e.code === ch.code)
          return (
            <button key={ch.code} type="button" aria-pressed={active} className={`chip ${ch.code.startsWith('B') ? 'benign' : 'path'}${excluded ? ' excluded' : ''}`} onClick={() => toggle(ch.code)} title={excluded ? r.excluded.find(e => e.code === ch.code)!.why : ch.why}>
              <b>{critLabel({ code: ch.code, strength: ch.strength })}</b>
              <span>{ch.why}</span>
              <em>{excluded ? 'not counted' : `${ch.code.startsWith('B') ? '−' : '+'}${{ VeryStrong: 8, Strong: 4, Moderate: 2, Supporting: 1, StandAlone: 8 }[ch.strength]}`}</em>
            </button>
          )
        })}
      </div>
      <Meter points={r.points} cls={r.cls} />
      <p className="hero-foot">{auto ? 'Playing a curation. Click any criterion to take over.' : `${applied.length} criteria · ${STRENGTH_SHORT.Supporting} = 1, Moderate = 2, Strong = 4, Very strong = 8 points (Tavtigian 2020).`}</p>
    </div>
  )
}

const QUOTES = [
  { src: 'r/ClinicalGenetics', title: 'Is AI helping in VUS classification yet?', text: 'Labs keep disagreeing on the same variant, and most of the work is still literature mining. Often the real answer is a wet-lab problem.' },
  { src: 'r/genetics', title: 'Two VUS in the same gene from exome sequencing', text: 'How do geneticists think about two uncertain variants in one gene? It depends on whether they sit on the same copy or opposite copies.' },
  { src: 'r/BRCA', title: 'Variant of uncertain significance in BRCA2', text: 'A BRCA2 VUS comes back and the family is left asking what it means, what to do now, and who will revisit it.' },
]

export function Landing() {
  const top = useMemo(() => {
    const csv = parseCsv(LAB_CSV)
    return reconcile(validateLabCsv(csv, autoMap(csv.rows[0], LAB_FIELDS), 'x', () => String(k++)).rows)[0]
  }, [])
  return (
    <>
      <section className="hero">
        <div className="hero-ring r1" aria-hidden="true" /><div className="hero-ring r2" aria-hidden="true" />
        <div className="hero-inner">
          <div className="hero-copy">
            <p className="eyebrow light">Variant interpretation workbench</p>
            <h1>Move a VUS with evidence you can audit, not a hunch.</h1>
            <p className="lede">Import a VCF or HGVS list. VUS Workbench maps every variant onto its transcript, checks splice sites, scores ACMG/AMP criteria with the Bayesian points system, and shows exactly why two labs disagree.</p>
            <div className="row gap wrap">
              <Link href="/app?sample=synthetic-panel" className="btn primary lg">Try the sample <span className="arrow" aria-hidden="true" /></Link>
              <Link href="/app?import=1" className="btn white-outline lg">Import your data</Link>
            </div>
            <p className="hero-note">Runs in your browser. Research use only, not a diagnostic. No account needed.</p>
          </div>
          <Hero />
        </div>
      </section>

      <section className="band">
        <div className="page">
          <p className="eyebrow">What curators keep running into</p>
          <h2>The hard part is not the call. It is defending it.</h2>
          <div className="quotes">
            {QUOTES.map(q => (
              <figure key={q.src} className="quote">
                <blockquote>{q.text}</blockquote>
                <figcaption><b>{q.src}</b> · {q.title} <span className="muted">(paraphrased)</span></figcaption>
              </figure>
            ))}
          </div>
        </div>
      </section>

      <section className="page">
        <p className="eyebrow">Workflow</p>
        <h2>From a VCF to a defensible classification</h2>
        <ol className="flow">
          <li><b>Import</b><span>VCF 4.x with multi-allelic splitting and left-normalisation, HGVS lists, lab tables, GTF subsets. A validation report names every bad line.</span></li>
          <li><b>Annotate</b><span>c. and p. from the transcript model with the HGVS 3′ rule, consequence, exon number, NMD, and a local splice-site strength delta.</span></li>
          <li><b>Score</b><span>PVS1 decision tree, PM2 from your frequency, calibrated PP3/BP4, PS3/BS3 from OddsPath, PM3/BP2 from phase. Every criterion logged with its evidence and editor.</span></li>
          <li><b>Reconcile</b><span>Lab calls side by side, criteria that differ, and the labs whose own criteria do not add up to the class they reported.</span></li>
          <li><b>Report</b><span>A printable curation report, a worklist CSV, curated calls in a shareable format, and a reclassification history.</span></li>
        </ol>
      </section>

      <section className="page split">
        <div>
          <p className="eyebrow">Discrepancy view</p>
          <h2>See the criterion behind the disagreement</h2>
          <p className="muted">In the bundled sample, two synthetic labs call <span className="mono">{top.variant}</span> Pathogenic and Benign. The workbench lines up their criteria, totals the points each lab's own evidence earns, and flags that a last-exon stop should not carry full PVS1.</p>
          <Link href="/app?sample=lab-discrepancy&tab=labs" className="btn outline">Open the discrepancy sample</Link>
        </div>
        <div className="card pad disc-preview">
          <div className="row between wrap gap-sm"><strong className="mono">{top.variant}</strong><span className="sev major">Major</span></div>
          {top.calls.map(c => (
            <div key={c.lab} className="disc-row">
              <span>{c.lab}</span>
              <span className="mono small">{c.criteria.map(critLabel).join(' ')}</span>
              <ClassBadge cls={c.cls} small />
              {c.mismatch && <span className="flag" title={`Criteria total ${c.points} points`}>criteria say {c.implied}</span>}
            </div>
          ))}
        </div>
      </section>

      <section className="page">
        <div className="features">
          {[
            ['Transcript-aware', 'Genomic to c. on either strand, intronic offsets, UTR numbering, dup vs ins, 3′ shifting.'],
            ['Splice checks', 'Canonical ±1/2, splice region, and a PWM donor/acceptor score for the reference and the variant, plus cryptic site scan.'],
            ['Points, not vibes', 'Supporting 1, Moderate 2, Strong 4, Very strong 8. P ≥ 10, LP 6 to 9, VUS 0 to 5, with VUS temperature.'],
            ['Two VUS, one gene', 'Phased genotypes decide trans vs cis, which drives PM3 in recessive genes and BP2 in cis.'],
            ['Local-first', 'Projects live in your browser. Sign in only if you want sync across devices.'],
            ['Audit trail', 'Every toggle, strength change and evidence entry is timestamped with the curator’s name.'],
          ].map(([h, t]) => <div key={h} className="feature"><h3>{h}</h3><p className="muted">{t}</p></div>)}
        </div>
      </section>

      <section className="cta">
        <div className="page row between wrap gap">
          <div><h2>Bring a VCF from your pipeline.</h2><p>Filter it to your panel, drop it in, and start curating. The docs show the exact bcftools commands.</p></div>
          <div className="row gap wrap"><Link href="/app" className="btn white lg">Open the workbench</Link><Link href="/docs" className="btn white-outline lg">Read the docs</Link></div>
        </div>
      </section>
    </>
  )
}
