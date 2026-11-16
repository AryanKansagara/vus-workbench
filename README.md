# VUS Workbench

A variant interpretation workbench for genetic counsellors, clinical scientists and researchers curating variants of uncertain significance. Research use only, not a diagnostic.

**Live:** https://vus-workbench.vercel.app

## The pain
Labs disagree on the same variant, the evidence behind a classification lives in spreadsheets and heads, and two VUS in one gene raise phase questions nobody records. VUS Workbench turns a VCF or HGVS list into a worklist where every ACMG/AMP criterion is scored, logged with its evidence and editor, and combined with the Bayesian points system, and where lab discrepancies show the criterion behind the disagreement.

## Features
- Import VCF 4.x (multi-allelic splitting, genotype re-expression, left-normalisation, INFO popmax AF / REVEL / SpliceAI), including VEP `CSQ` and SnpEff `ANN` annotated files and ClinVar's own VCF; CSV/TSV HGVS lists with column mapping; ClinVar tabular downloads; lab classification tables; GTF/GFF3 transcript subsets; FASTA regions (samtools faidx). Validation report with line numbers.
- Real gene models on demand: the MANE Select transcript and GRCh38 sequence for any gene from the Ensembl REST API.
- Real lab calls on demand: every ClinVar submitter's classification for a worklist variant, via NCBI E-utilities.
- Engine (plain TypeScript, tested): HGVS c./p. parser and formatter, genomic to c. on either strand with intronic offsets and UTR numbering, HGVS 3' rule and dup detection, protein effect by translation, consequence and splice-region flags, NMD (50 nt rule), local PWM donor/acceptor strength delta with cryptic-site scan.
- ACMG/AMP: PVS1 decision tree, PM2 at Supporting, calibrated PP3/BP4 (REVEL, SpliceAI), PS3/BS3 from OddsPath, PM3/BP2 from phase, double-counting rules, Tavtigian points to class with posterior probability and VUS temperature.
- Discrepancy view: labs side by side, criteria spread, and labs whose own criteria do not add up to their call.
- Per-variant evidence log and reclassification history, dated snapshots with a compare view.
- Printable curation report, worklist CSV, curated calls CSV (re-importable as a lab table).
- Local-first (IndexedDB), optional Supabase sync, account export/delete, light and dark themes.

## Tested against
The VCF 4.3 specification example (samtools/hts-specs), an excerpt of ClinVar's GRCh38 VCF, ClinVar record VCV000038266, and Ensembl MANE Select models for BRCA1 and BRCA2 (c. to genomic mapping and protein changes checked against the public ClinVar descriptions). Point thresholds are checked against the Bayesian posteriors (prior 0.10, very strong odds 350) and REVEL intervals against Pejaver et al. 2022, Table 2.

## Known gaps
- Allele frequencies, REVEL and SpliceAI are not looked up; import them in the VCF or enter them.
- Ensembl models are GRCh38, MANE Select only. GRCh37 needs your own GTF and FASTA.
- ACMG codes read from ClinVar comments are free text and need checking.
- Sequence variants only; CNVs and structural variants are skipped.

Samples are synthetic and labelled: two synthetic genes with generated sequence, plus public HGVS descriptions for BRCA1, BRCA2 and CFTR with made-up lab calls and illustrative evidence values.

## Run locally
```
npm install
npm run dev     # http://localhost:4394
npm test        # engine tests
npm run build
```
Optional sync: put `VITE_SUPABASE_URL` and `VITE_SUPABASE_PUBLISHABLE_KEY` in `.env`.
