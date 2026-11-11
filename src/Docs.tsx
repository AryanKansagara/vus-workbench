import { LAB_CSV, REAL_VARIANTS_CSV, SAMPLE_GTF, sampleVcf } from './engine/sample.ts'
import { Code, download, Link } from './ui.tsx'

const TOC = [['vcf', 'VCF from your pipeline'], ['csv', 'HGVS variant lists'], ['labs', 'Lab classification tables'], ['models', 'Transcript models and reference'], ['rules', 'Rules and dated assumptions'], ['limits', 'Limits']]

export default function Docs() {
  return (
    <section className="page docs">
      <aside className="docs-toc" aria-label="On this page">
        <p className="eyebrow">Docs</p>
        {TOC.map(([id, t]) => <a key={id} href={`#${id}`}>{t}</a>)}
      </aside>
      <div className="docs-body">
        <h1>Getting variants in</h1>
        <p className="lede-dark">Everything is parsed in your browser. Nothing is uploaded unless you sign in to sync, and then only the worklist (not raw files). Rules and thresholds below are as of November 2026.</p>

        <h2 id="vcf">VCF from your pipeline</h2>
        <p>Any VCF 4.x works: header lines, INFO and FORMAT are read, multi-allelic sites are split (genotypes are re-expressed per allele, so <span className="mono">1/2</span> becomes <span className="mono">1/0</span> and <span className="mono">0/1</span>), and indels are trimmed and left-aligned against loaded reference sequence. Symbolic alleles (<span className="mono">&lt;DEL&gt;</span>) are skipped with a warning. Keep the file to your panel or genes of interest: the importer stops at 20,000 variants.</p>
        <Code label="bash" code={`# 1. Restrict to panel regions, split multi-allelics, left-align (GRCh38 reference)
bcftools view -R panel.bed -Ou proband.vcf.gz \\
  | bcftools norm -m -any -f GRCh38.fa -Ov -o proband.panel.vcf

# 2. Optional: put popmax AF, REVEL and SpliceAI into INFO so they import with the variants
#    (VEP with the REVEL and SpliceAI plugins, or bcftools annotate from pre-scored VCFs)
bcftools annotate -a gnomad.v4.1.sites.vcf.gz -c INFO/AF_grpmax proband.panel.vcf -Ov -o proband.annot.vcf`} />
        <p>Recognised INFO keys: <span className="mono">AF_popmax, AF_grpmax, gnomAD_AF_popmax, gnomADe_AF, gnomADg_AF</span> (popmax frequency), <span className="mono">REVEL</span>, and <span className="mono">SpliceAI</span> in the standard <span className="mono">ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|…</span> format. Phased genotypes (<span className="mono">0|1</span>) with an optional <span className="mono">PS</span> phase set drive PM3 / BP2.</p>
        <button type="button" className="btn ghost sm" onClick={() => download('SYN-PROBAND-01.synthetic.vcf', sampleVcf(), 'text/plain')}>Download the synthetic sample VCF</button>

        <h2 id="csv">HGVS variant lists</h2>
        <p>A CSV or TSV with one variant per row. Map columns on import; common names are detected. HGVS c. is required unless chromosome, position, REF and ALT are all given. Prefixes like <span className="mono">NM_000059.4(BRCA2):</span> are accepted in the cell or via separate gene / transcript columns.</p>
        <Code label="csv" code={REAL_VARIANTS_CSV.trim()} />
        <p className="muted small">Columns: gene, transcript, hgvs_c, hgvs_p, popmax_af, revel, spliceai, genotype or zygosity (0/1, 0|1, het, hom), phase_set, sample. The values in this example are illustrative, not gnomAD lookups.</p>

        <h2 id="labs">Lab classification tables</h2>
        <p>One row per lab call. Classification accepts P, LP, VUS, LB, B or the full words. Criteria accept the usual shorthand: <span className="mono">PVS1, PM2_Supporting, PP3_Mod, PS3_P</span>. The discrepancy view totals each lab's own criteria with the points system and flags labs whose call does not match their evidence.</p>
        <Code label="csv" code={LAB_CSV.trim().split('\n').slice(0, 6).join('\n')} />
        <p>From ClinVar: copy the Submitted interpretations table from the variation page, or pull the rows for one VariationID from the weekly <span className="mono">submission_summary.txt.gz</span> (column 1 VariationID, 2 significance, 3 date evaluated, 10 submitter):</p>
        <Code label="bash" code={`curl -sO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz
zcat submission_summary.txt.gz \\
  | awk -F'\\t' -v id=51579 -v v='NM_000059.4:c.68-7T>A' \\
      'BEGIN{OFS=","; print "variant,lab,classification,date_evaluated"} $1==id {print v, "\\""$10"\\"", $2, $3}' \\
  > clinvar-calls.csv`} />
        <p className="muted small">ClinVar summaries rarely list criteria; add a criteria column by hand where a lab's report shows them.</p>

        <h2 id="models">Transcript models and reference</h2>
        <p>Two synthetic genes (SYNG1 on the plus strand, SYNG2 on the minus strand) ship with sequence so every feature works offline. For a real gene, import a GTF subset for its transcript, and a FASTA of the region to enable protein prediction, NMD and splice scoring.</p>
        <Code label="bash" code={`# GTF subset for one transcript (GENCODE / Ensembl, GRCh38)
zcat gencode.v47.annotation.gtf.gz | grep 'transcript_id "ENST00000380152' > brca2.gtf

# Reference region (samtools faidx output keeps the coordinates in the header)
samtools faidx GRCh38.fa chr13:32310000-32410000 > brca2.region.fa`} />
        <p>Contig names match with or without the <span className="mono">chr</span> prefix. Ensembl GTF CDS rows exclude the stop codon; include the <span className="mono">stop_codon</span> rows (grep keeps them) so the CDS length is a multiple of three.</p>
        <button type="button" className="btn ghost sm" onClick={() => download('synthetic-genes.gtf', SAMPLE_GTF, 'text/plain')}>Download the synthetic GTF</button>

        <h2 id="rules">Rules and dated assumptions</h2>
        <ul>
          <li><b>Points</b>: Supporting 1, Moderate 2, Strong 4, Very strong 8; benign negative. P ≥ 10, LP 6 to 9, VUS 0 to 5, LB −1 to −6, B ≤ −7 (Tavtigian et al., Hum Mutat 2020). BA1 is stand-alone benign.</li>
          <li><b>PVS1</b>: simplified Abou Tayoun et al. 2018 tree. NMD predicted unless the stop is in the last exon or within 50 nt of the last exon-exon junction; escapees get Strong if over 10% of the protein is removed or a critical region (gene setting) is lost, else Moderate. Canonical ±1/2 splice variants are modelled as skipping the adjacent exon. Initiation codon: Moderate.</li>
          <li><b>PM2</b>: applied at Supporting (ClinGen SVI, 2020) when popmax AF is at or below the gene threshold. BS1 and BA1 from gene thresholds (default 1% and 5%).</li>
          <li><b>PP3 / BP4</b>: REVEL calibration of Pejaver et al. 2022 for missense; SpliceAI ≥ 0.2 for PP3 and ≤ 0.1 for BP4 (ClinGen SVI splicing, Walker et al. 2023). PP3 is never combined with PVS1.</li>
          <li><b>Local splice model</b>: position weight matrices for the 9-nt donor and 16-nt acceptor from published consensus frequencies, scored in bits. A drop of 20% weakens and 60% abolishes; it supports PP3 only when no SpliceAI score is entered.</li>
          <li><b>PS3 / BS3</b>: OddsPath thresholds of Brnich et al. 2019 (over 18.7 Strong, 4.3 Moderate, 2.1 Supporting; under 0.053, 0.23, 0.48 for benign). Without OddsPath, Supporting.</li>
          <li><b>PM3 / BP2</b>: phase from phased genotypes. Recessive gene: in trans with P/LP is PM3 Moderate, phase unknown is PM3 Supporting. In cis with P/LP, or in trans in a dominant gene, is BP2.</li>
        </ul>

        <h2 id="limits">Limits</h2>
        <ul>
          <li>Research use only. Not a validated diagnostic system.</li>
          <li>No population database, ClinVar or literature lookups: frequencies and scores are entered or imported by you.</li>
          <li>Sequence variants only. CNVs and structural variants are skipped.</li>
          <li>The splice PWM is a MaxEnt-style approximation, not MaxEntScan or SpliceAI.</li>
          <li>Indels spanning an exon-intron boundary get a splice consequence but no protein prediction.</li>
        </ul>
        <p><Link href="/app" className="btn primary">Open the workbench</Link></p>
      </div>
    </section>
  )
}
