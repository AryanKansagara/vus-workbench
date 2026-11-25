// Filed bugs, pinned with test.fail: each asserts the correct behaviour and starts passing when its issue is fixed.
import { choose, expect, importBtn, importFile, issues, MINI_FASTA, MINI_GFF3, MINI_GTF, miniHgvs, newWorklist, openVariant, row, saveAs, tab, test, v, vcf, withMini, cells } from './fixtures'

const crit = (page: Parameters<typeof row>[0], c: string) => row(page, c).locator('.crit')

/** A second model on the same bases, with a 6 nt middle exon that has no real splice sites. */
async function odd(page: Parameters<typeof row>[0], lines: string[]) {
  await newWorklist(page)
  await importFile(page, 'mini.fa', MINI_FASTA)
  const f = (t: string, a: number, b: number) => `chr9\ttest\t${t}\t${a}\t${b}\t.\t+\t.\tgene_id "ODD"; transcript_id "NM_900003.1"; gene_name "ODD";`
  await importFile(page, 'odd.gtf', [[1016, 1105], [1131, 1136], [1151, 1207]].flatMap(([a, b]) => [f('exon', a, b), f('CDS', a, b)]).join('\n') + '\n')
  await importFile(page, 'odd.csv', ['gene,hgvs_c', ...lines, ''].join('\n'))
}

test('#2 a deletion that removes a whole donor site is a null variant, not an in-frame deletion', async ({ page }) => {
  await miniHgvs(page, ['MINI,c.70_75+5del,0'])
  // Exon 2 is lost: 60 nt in frame, 20 of 32 residues. PVS1_Strong + PM2.
  await expect(crit(page, 'c.70_75+5del')).toHaveText(['PVS1_Strong', 'PM2'], { timeout: 5000 })
})
test('#2 a deletion of a whole exon is not intronic', async ({ page }) => {
  await miniHgvs(page, ['MINI,c.16-20_75+18del,0'])
  await expect(crit(page, 'c.16-20_75+18del')).toHaveText(['PVS1_Strong', 'PM2'], { timeout: 5000 })
})
test('#3 an in-frame insertion that brings a stop codon is a null variant', async ({ page }) => {
  // ATG GCT + GAT TAA: stop at codon 4, before the last 50 nt of exon 2. PVS1 8 + PM2 1.
  await miniHgvs(page, ['MINI,c.6_7insGATTAA,0'])
  await expect(crit(page, 'c.6_7ins')).toHaveText(['PVS1', 'PM2'], { timeout: 5000 })
})
test('#4 a decimal comma is not read as a frequency of 100%', async ({ page }) => {
  await withMini(page)
  await importFile(page, 'eu.csv', 'gene;hgvs_c;popmax_af\nMINI;c.31C>T;0,0001\n')
  await expect(crit(page, 'c.31C>T')).toHaveText(['PM2'], { timeout: 5000 })
})
test('#4 SpliceAI without a score is not benign evidence', async ({ page }) => {
  // SpliceAI writes "ALT|GENE|.|.|.|.|.|.|.|." for alleles it cannot score.
  await withMini(page)
  await importFile(page, 'a.vcf', vcf([v(1109, 'T', 'C', 'AF_popmax=0;SpliceAI=C|MINI|.|.|.|.|.|.|.|.')]))
  // BP7 rests on the local splice model here; BP4 needs a SpliceAI score.
  await expect(crit(page, 'c.54T>C')).toHaveText(['PM2', 'BP7'], { timeout: 5000 })
})
test('#4 scores outside 0 to 1 in a VCF are refused', async ({ page }) => {
  await withMini(page)
  await importFile(page, 'a.vcf', vcf([v(1086, 'C', 'T', 'AF_popmax=-1;REVEL=7')]))
  await expect(crit(page, 'c.31C>T')).toHaveCount(0, { timeout: 5000 })
})
test('#5 BP6 counts towards benign', async ({ page }) => {
  await newWorklist(page)
  await importFile(page, 'labs.csv', 'variant,lab,classification,criteria\nBRCA2:c.1A>G,Lab A,Likely benign,"BP4, BP6"\n')
  await expect(page.locator('tr.sum td').nth(1)).toHaveText('-2', { timeout: 5000 })
})
test('#5 a strength in brackets is kept, and a criterion marked not met is not counted', async ({ page }) => {
  await newWorklist(page)
  await importFile(page, 'labs.csv', 'variant,lab,classification,criteria\nBRCA2:c.1A>G,Lab A,VUS,"PS3 (Moderate), PM2 not met"\n')
  await expect(page.locator('tr.sum td').nth(1)).toHaveText('2', { timeout: 5000 })
})
test('#6 two spellings of one variant are compared', async ({ page }) => {
  await newWorklist(page)
  await importFile(page, 'labs.csv', 'variant,gene,lab,classification\nNM_007294.4:c.5266dupC,BRCA1,A,Pathogenic\nNM_007294.4:c.5266dup,BRCA1,B,Benign\n')
  await expect(page.locator('.disc')).toHaveCount(1, { timeout: 5000 })
  await expect(page.locator('.disc .sev')).toHaveText('Major (P/LP vs B/LB)')
})
test('#6 the latest call wins with ClinVar dates', async ({ page }) => {
  await newWorklist(page)
  await importFile(page, 'labs.csv', 'variant,lab,classification,date_evaluated\nBRCA2:c.1A>G,A,Benign,"Jan 09, 2024"\nBRCA2:c.1A>G,A,Pathogenic,"Sep 18, 2019"\n')
  await expect(page.locator('.disc .cls-row .cls')).toHaveText(['B'], { timeout: 5000 })
})
test('#7 a c. description on another transcript is not placed on the loaded one', async ({ page }) => {
  await miniHgvs(page, ['MINI,NM_900002.1,c.11G>A'], 'gene,transcript,hgvs_c')
  await expect(page.locator('table.worklist')).not.toContainText('p.(Trp4Ter)', { timeout: 5000 })
})
test.fail('#8 PM3 does not depend on row order', async ({ page }) => {
  await withMini(page)
  // Two likely pathogenic partners: one unphased, one in trans. SVI PM3 v1.0: in trans is 1.0 point, Moderate.
  await importFile(page, 'a.vcf', vcf([v(1022, 'G', 'T', 'AF_popmax=0', '0/1'), v(1026, 'G', 'A', 'AF_popmax=0', '0|1'), v(1086, 'C', 'T', 'AF_popmax=0', '1|0')]))
  await tab(page, 'Genes')
  await page.getByLabel('MINI inheritance').selectOption('AR')
  await tab(page, 'Worklist')
  await expect(crit(page, 'c.31C>T')).toHaveText(['PM2', 'PM3'], { timeout: 5000 })
})
test.fail('#8 an unphased likely pathogenic partner is 0.25 points, below Supporting', async ({ page }) => {
  await withMini(page)
  await importFile(page, 'a.vcf', vcf([v(1026, 'G', 'A'), v(1086, 'C', 'T')]))
  await tab(page, 'Genes')
  await page.getByLabel('MINI inheritance').selectOption('AR')
  await tab(page, 'Worklist')
  await expect(crit(page, 'c.31C>T')).toHaveText(['PM2'], { timeout: 5000 })
})
test.fail('#9 an OddsPath in the indeterminate range gives no PS3', async ({ page }) => {
  await miniHgvs(page, ['MINI,c.31C>T,'])
  await openVariant(page, 'c.31C>T')
  await page.getByRole('button', { name: 'Add functional result' }).click()
  await page.getByLabel('Assay').fill('Minigene')
  await page.getByLabel('OddsPath').fill('1')
  await page.getByRole('button', { name: 'Add', exact: true }).click()
  await expect(page.locator('.func-list')).toContainText('OddsPath 1')
  await expect(page.locator('.sticky')).toContainText('0 pts', { timeout: 5000 })
})
test.fail('#10 a percent sign in a GFF3 note does not take the page down', async ({ page }) => {
  await newWorklist(page)
  // Feeds the file without leaving an uncaught error behind: the crash itself removes the page.
  await page.evaluate(() => addEventListener('error', (e) => e.preventDefault()))
  await choose(page, 'a.gff3', MINI_GFF3.replace('ID=g1;', 'ID=g1;Note=100% match;'))
  await expect(importBtn(page, 'a.gff3')).toBeEnabled({ timeout: 5000 })
})
test.fail('#11 an unclosed quote does not pass another row off as criteria', async ({ page }) => {
  await newWorklist(page)
  await choose(page, 'labs.csv', 'variant,lab,classification,criteria\nBRCA2:c.1A>G,A,VUS,"PM2\nBRCA2:c.1A>G,B,Pathogenic,PVS1\n')
  await expect(issues(page, 'labs.csv').locator('li.error')).toHaveCount(1, { timeout: 5000 })
})
test.fail('#12 a FASTA shorter than its header says is reported', async ({ page }) => {
  await newWorklist(page)
  await choose(page, 'mini.fa', MINI_FASTA.trim().slice(0, -30))
  await expect(issues(page, 'mini.fa')).toBeVisible({ timeout: 5000 })
})
test.fail('#12 a VCF row with fewer columns than the header is reported', async ({ page }) => {
  await withMini(page)
  await choose(page, 'cut.vcf', vcf([]) + 'chr9\t1086\t.\tC\tT\t60\tPASS\tAF_popmax=0')
  await expect(issues(page, 'cut.vcf').locator('li.error')).toHaveCount(1, { timeout: 5000 })
})
test.fail('#13 a GTF without stop_codon rows does not turn an in-frame deletion into a stop', async ({ page }) => {
  await newWorklist(page)
  await importFile(page, 'mini.fa', MINI_FASTA)
  await importFile(page, 'mini.gtf', MINI_GTF.split('\n').filter((l) => !l.includes('stop_codon')).join('\n'))
  await importFile(page, 'list.csv', 'gene,hgvs_c\nMINI,c.5_7del\n')
  await expect(row(page, 'c.5_7del')).toContainText('p.(Ala2del)', { timeout: 5000 })
})
test.fail('#14 an insertion after the start codon is not start-lost', async ({ page }) => {
  await miniHgvs(page, ['MINI,c.3_4insA,'])
  await expect(row(page, 'c.3_4insA')).toContainText('Frameshift', { timeout: 5000 })
})
test.fail('#15 an intron offset from the middle of an exon is refused', async ({ page }) => {
  // c.30 is inside exon 2, so c.30+5 is not a position. The base there happens to be the stated one.
  await miniHgvs(page, ['MINI,c.30+5C>T,0'])
  await expect(row(page, 'c.30+5C>T').locator('.flag')).toHaveText('check', { timeout: 5000 })
})
test.fail('#15 a stated duplicated base that the reference does not have is refused', async ({ page }) => {
  await miniHgvs(page, ['MINI,c.11dupT,0'])
  await expect(row(page, 'c.11dup').locator('.flag')).toHaveText('check', { timeout: 5000 })
})
test.fail('#16 skipping an exon in a critical region is PVS1_Strong even under 10%', async ({ page }) => {
  // ODD: coding 90 + 6 + 57 = 153 nt = 50 residues + stop. Skipping exon 2 removes residues 31-32, 2 of 50 = 4%.
  await odd(page, ['ODD,c.96+1T>G'])
  await tab(page, 'Genes')
  await page.getByLabel('ODD critical regions').fill('30-35')
  await page.getByLabel('ODD critical regions').blur()
  await tab(page, 'Worklist')
  await expect(crit(page, 'c.96+1T>G')).toHaveText(['PVS1_Strong'], { timeout: 5000 })
})
test.fail('#17 a change that raises the score of a weak native site is not called abolishing', async ({ page }) => {
  // ODD's exon 3 acceptor has no AG, so its native score is negative. T>C at -3 moves towards the consensus (C 75%, T 19%).
  await odd(page, ['ODD,c.97-3T>C'])
  await expect(crit(page, 'c.97-3T>C')).toHaveCount(0, { timeout: 5000 })
})
test.fail('#18 BP7 starts at intron position +7', async ({ page }) => {
  await miniHgvs(page, ['MINI,c.75+7T>C,0.001,0'], 'gene,hgvs_c,popmax_af,spliceai')
  await expect(crit(page, 'c.75+7T>C')).toHaveText(['BP4', 'BP7'], { timeout: 5000 })
})
test.fail('#19 a duplicated residue is named dup', async ({ page }) => {
  await miniHgvs(page, ['MINI,c.31_33dup,0'])
  await expect(row(page, 'c.31_33dup')).toContainText('p.(Arg11dup)', { timeout: 5000 })
})
test.fail('#20 the same variant in two samples is compared twice', async ({ page }) => {
  await miniHgvs(page, ['MINI,c.31C>T,0,S1', 'MINI,c.31C>T,0.2,S2'], 'gene,hgvs_c,popmax_af,sample')
  await page.getByRole('button', { name: 'Save snapshot' }).click()
  await tab(page, 'History')
  await page.locator('.arrow-sep + label select').selectOption('current')
  await expect(page.getByText('2 unchanged')).toBeVisible({ timeout: 5000 })
})
test.fail('#21 a lab name cannot become a spreadsheet formula', async ({ page }) => {
  await newWorklist(page)
  await importFile(page, 'labs.csv', 'variant,lab,classification\nBRCA2:c.1A>G,"=HYPERLINK(""http://x.test"")",VUS\n')
  const [, csv] = await saveAs(page, page.getByRole('button', { name: 'Reconciliation CSV' }))
  expect(cells(csv)[1][2]).not.toMatch(/^=/)
})
test('#22 an evidence value survives a criterion clicked straight after it', async ({ page }) => {
  await miniHgvs(page, ['MINI,c.31C>T,'])
  await openVariant(page, 'c.31C>T')
  await page.getByRole('button', { name: 'Show all 26 criteria to add one' }).click()
  const box = page.locator('.crit-row').filter({ hasText: 'PP4' }).getByRole('checkbox')
  await page.getByLabel('Popmax AF').fill('0.2')
  // Same tick: the blur commits the value and the click runs before the page re-renders.
  await box.evaluate((el: HTMLInputElement) => { (document.activeElement as HTMLElement).blur(); el.click() })
  await expect(page.locator('.audit')).toContainText('PP4')
  await expect(page.getByLabel('Popmax AF')).toHaveValue('0.2', { timeout: 5000 })
})
test('#22 undo of one change keeps the later ones', async ({ page }) => {
  await miniHgvs(page, ['MINI,c.31C>T,'])
  await openVariant(page, 'c.31C>T')
  await page.getByRole('button', { name: 'Show all 26 criteria to add one' }).click()
  const box = page.locator('.crit-row').filter({ hasText: 'PP4' }).getByRole('checkbox')
  await box.click()
  await page.getByLabel('Popmax AF').fill('0.2')
  await page.getByLabel('Popmax AF').blur()
  await expect(page.locator('.audit')).toContainText('Popmax AF blank → 0.2')
  await page.locator('.toast').filter({ hasText: 'PP4' }).getByRole('button', { name: 'Undo' }).click()
  await expect(box).not.toBeChecked()
  await expect(page.locator('.audit')).toContainText('Undone')
  await expect(page.getByLabel('Popmax AF')).toHaveValue('0.2')
})
