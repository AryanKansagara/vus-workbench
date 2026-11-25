import { choose, expect, importBtn, importFile, issues, MINI_FASTA, MINI_GFF3, MINI_GTF, newWorklist, pending, row, tab, test, toastText, v, vcf, withMini } from './fixtures'

test('VCF with CRLF: multi-allelic split, filters and wrong-build REF reported by line', async ({ page }) => {
  await withMini(page)
  const text = vcf([v(1026, 'G', 'A,T', 'AF_popmax=0,0.2', '1/2'), v(1086, 'G', 'T'), v(1131, 'G', 'A').replace('PASS', 'LowQual'), v(1189, 'G', '<DEL>'), v(1109, 'T', 'C', 'AF_popmax=0', '0/0')]).replace(/\n/g, '\r\n')
  await choose(page, 'a.vcf', text)
  await expect(issues(page, 'a.vcf')).toContainText('Multi-allelic site split into 2 records')
  await expect(issues(page, 'a.vcf').locator('li.error')).toContainText(['REF G does not match the reference (C) at chr9:1086; wrong genome build?'])
  await expect(issues(page, 'a.vcf')).toContainText('Symbolic/structural ALT <DEL> skipped')
  await expect(issues(page, 'a.vcf')).toContainText('1 allele skipped: P1 is homozygous reference or not called there')
  await importBtn(page, 'a.vcf').click()
  await expect(page.locator('table.worklist tbody tr')).toHaveCount(3)
  // TGG > TTG is Trp4Leu with AF 0.2 (second value): BA1.
  await expect(row(page, 'c.11G>T')).toContainText('p.(Trp4Leu)')
  await expect(row(page, 'c.11G>T')).toContainText('BA1')
  await expect(row(page, 'c.75+1G>A')).toContainText('FILTER LowQual')
})

test('VEP-annotated VCF is used before a model is loaded, and no model is said', async ({ page }) => {
  await newWorklist(page)
  const t = ['##fileformat=VCFv4.2', '##INFO=<ID=CSQ,Number=.,Type=String,Description="x. Format: Allele|SYMBOL|Feature|HGVSc|HGVSp|MANE_SELECT|REVEL">', '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
    '13\t100\t.\tA\tT\t.\tPASS\tCSQ=T|BRCA2|NM_000059.4|NM_000059.4:c.9976A>T|NP_000050.3:p.Lys3326Ter|NM_000059.4|', ''].join('\n')
  await importFile(page, 'vep.vcf', t)
  await expect(row(page, 'c.9976A>T')).toContainText('Nonsense')
  // 3418 - 3326 + 1 = 93 of 3418 residues = 2.7%: PVS1_Moderate.
  await expect(row(page, 'c.9976A>T').locator('.crit')).toHaveText(['PVS1_Moderate'])
})

test('ClinVar VCF adds the aggregate call to the lab comparison', async ({ page }) => {
  await newWorklist(page)
  await importFile(page, 'cv.vcf', ['##fileformat=VCFv4.1', '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', '13\t100\t1\tA\tT\t.\t.\tGENEINFO=BRCA2:675;CLNSIG=Benign/Likely_benign;CLNREVSTAT=criteria_provided,_multiple_submitters', ''].join('\n'))
  await tab(page, 'Lab discrepancies')
  await expect(page.locator('.disc .cls-row .cls').first()).toHaveText('LB')
})

test('HGVS list: tab-separated with CRLF, BOM and quoted cells; bad rows say why', async ({ page }) => {
  await withMini(page)
  const t = '﻿gene\thgvs_c\tpopmax_af\trevel\tzygosity\r\nMINI\tc.31C>T\t0\t0.95\thet\r\nMINI\t"c.11G>A"\t0\t\t\r\nMINI\tc.31C>\t\t\t\r\nMINI\tc.40G>A\t2\t\tmaybe\r\nMINI\tc.31C>T\t0\t\t\r\n'
  await choose(page, 'list.tsv', t)
  const i = issues(page, 'list.tsv')
  await expect(i.locator('li')).toHaveText([/line 4error.*"c\.31C>" is not a supported/, /line 5warning.*af "2" ignored \(expected 0 to 1\)/, /line 5warning.*Genotype "maybe" not recognised/, /line 6warning.*Duplicate variant c\.31C>T skipped/])
  await importBtn(page, 'list.tsv').click()
  await expect(page.locator('table.worklist tbody tr')).toHaveCount(3)
  await expect(row(page, 'c.31C>T')).toContainText('+5')
})

test('importing the same list again adds nothing', async ({ page }) => {
  await withMini(page)
  const t = 'gene,hgvs_c,popmax_af\nMINI,c.31C>T,0\n'
  await importFile(page, 'a.csv', t)
  await importFile(page, 'b.csv', t)
  await toastText(page, 'Added 0 variants, 1 already in the worklist')
  await expect(page.locator('table.worklist tbody tr')).toHaveCount(1)
})

test('ClinVar tabular download supplies gene, c., p. and the call', async ({ page }) => {
  await newWorklist(page)
  await importFile(page, 'clinvar_result.txt', 'Name\tGene(s)\tProtein change\tGermline classification\tGermline review status\nNM_000059.4(BRCA2):c.9976A>T (p.Lys3326Ter)\tBRCA2\tK3326*\tBenign\treviewed by expert panel\n')
  await expect(row(page, 'c.9976A>T')).toContainText('p.Lys3326Ter')
  await tab(page, 'Lab discrepancies')
  await expect(page.locator('.disc')).toContainText('ClinVar (reviewed by expert panel)')
})

test('lab table: criteria totals worked by hand and unknown classes refused', async ({ page }) => {
  await newWorklist(page)
  await choose(page, 'labs.csv', 'variant,lab,classification,criteria\nBRCA2:c.1A>G,Lab A,Pathogenic,"PS3, PM2"\nBRCA2:c.1A>G,Lab B,Likely benign,"BS1;BP4"\nBRCA2:c.1A>G,Lab C,maybe,\n')
  await expect(issues(page, 'labs.csv')).toContainText('Classification "maybe" not recognised')
  await importBtn(page, 'labs.csv').click()
  // 4 + 1 = 5 is VUS, not Pathogenic. -4 - 1 = -5 is likely benign.
  const card = page.locator('.disc')
  await expect(card).toContainText("Lab A's own criteria total 5 points (VUS), not Pathogenic.")
  await expect(card.locator('tr.sum td')).toHaveText(['Points', /^5/, '-5'])
})

for (const [name, text] of [['mini.gff3', MINI_GFF3], ['mini.gtf', MINI_GTF]] as const) {
  test(`${name} and FASTA build the same model`, async ({ page }) => {
    await newWorklist(page)
    await importFile(page, name, text)
    await importFile(page, 'mini.fa', MINI_FASTA)
    await tab(page, 'Genes')
    const r = page.locator('tr').filter({ hasText: 'NM_900001.1' })
    await expect(r.locator('td')).toHaveText(['NM_900001.1', 'MINI', 'chr9', '+', '3', '32 aa', 'loaded', name, 'Remove'])
  })
}

test('pasted text imports like a file', async ({ page }) => {
  await newWorklist(page)
  await page.getByText('Paste instead').click()
  await page.getByLabel('Pasted data').fill('gene,hgvs_c\nBRCA2,c.68-7T>A\n')
  await page.getByRole('button', { name: 'Read pasted text' }).click()
  await importBtn(page, 'pasted text').click()
  await expect(row(page, 'c.68-7T>A')).toContainText('Splice region')
})

const WRONG: [string, string, string][] = [
  ['notes.txt', 'hello there', 'could not tell the format'],
  ['empty.csv', '', 'could not tell the format'],
  ['data.json', '{"a":1}', 'could not tell the format'],
  ['a.vcf.gz', 'x', 'is compressed'],
]
for (const [name, text, message] of WRONG) test(`${name} is refused with a message`, async ({ page }) => {
  await newWorklist(page)
  await choose(page, name, text)
  await toastText(page, message)
  await expect(page.locator('.pending')).toHaveCount(0)
})

const NOTHING: [string, string, string][] = [
  ['empty.vcf', '', 'No #CHROM header line found'],
  ['cut.vcf', vcf([]) + 'chr9\t1026\t.\tG', 'Expected at least 8 tab-separated columns, found 4'],
  ['noheader.vcf', 'chr9\t1026\t.\tG\tA', 'Data before the #CHROM header line'],
  ['noheader.fa', 'ACGT\n', 'No FASTA records found'],
  ['exons.gtf', MINI_GTF.split('\n').filter((l) => l.includes('\texon\t')).join('\n') + '\n', 'No complete transcripts found'],
  ['cut.csv', 'gene,hgvs_c\nMINI,c.31C', 'is not a supported c. description'],
  ['labs.csv', 'variant,lab,classification\nBRCA2:c.1A>G,Lab A,Patho', 'Classification "Patho" not recognised'],
]
for (const [name, text, message] of NOTHING) test(`${name} gives an error and nothing to import`, async ({ page }) => {
  await newWorklist(page)
  await choose(page, name, text)
  await expect(issues(page, name)).toContainText(message)
  await expect(importBtn(page, name)).toBeDisabled()
  await pending(page, name).getByRole('button', { name: 'Discard' }).click()
  await expect(page.locator('.pending')).toHaveCount(0)
})

test('text from a file is shown as text', async ({ page }) => {
  await newWorklist(page)
  await importFile(page, 'x.csv', 'variant,lab,classification\nBRCA2:c.1A>G,<img src=x onerror=window.pwned=1>,VUS\n')
  await expect(page.locator('.disc')).toContainText('<img src=x')
  expect(await page.evaluate(() => 'pwned' in window)).toBe(false)
})
