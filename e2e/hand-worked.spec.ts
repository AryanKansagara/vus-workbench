// Every number here was derived on paper from the MINI gene in fixtures.ts and the published tables the Docs page cites
// (points: Tavtigian 2020; posterior: Tavtigian 2018, prior 0.1 and odds 350^(points/8); REVEL: Pejaver 2022; PVS1: Abou Tayoun 2018).
import { expect, importFile, openVariant, row, tab, test, v, vcf, withMini } from './fixtures'

const ROWS = [
  v(1026, 'G', 'A', 'AF_popmax=0'), // c.11G>A, TGG > TAG
  v(1086, 'C', 'T', 'AF_popmax=0;REVEL=0.95'), // c.31C>T, CGT > TGT
  v(1018, 'GGCT', 'G'), // deletes GCT, or equally CTG one base on
  v(1131, 'G', 'A'), // c.75+1G>A
  v(1189, 'G', 'T'), // c.94G>T, GAA > TAA in the last exon
  v(1109, 'T', 'C', 'AF_popmax=0.06;SpliceAI=C|MINI|0.00|0.00|0.00|0.00|0|0|0|0'), // c.54T>C, GCT > GCC
]

test.beforeEach(async ({ page }) => {
  await withMini(page)
  await importFile(page, 'mini.vcf', vcf(ROWS))
  await expect(page.locator('table.worklist tbody tr')).toHaveCount(6)
})

test('an early stop: codon 4 is before the last 50 nt of exon 2, so 8 + 1 points and 98.8%', async ({ page }) => {
  // Stop at transcript base 5 + 9 + 1 = 15. Last junction after base 20 + 60 = 80; 15 < 80 - 50. PVS1 (8) + PM2_Supporting (1) = 9.
  // 350^(9/8) = 727.9; 72.79 / (72.69 + 1) = 0.988.
  await openVariant(page, 'c.11G>A')
  await expect(page.locator('.var-title')).toContainText('p.(Trp4Ter)')
  const facts = page.locator('dl.facts')
  await expect(facts).toContainText('Nonsense')
  await expect(facts).toContainText('Exon 1 of 3')
  await expect(facts).toContainText('Codon 4')
  // Codons 4 to 32 are lost: 29 of 32 = 90.6%.
  await expect(facts).toContainText('90.6%')
  await expect(facts.locator('div').filter({ hasText: /^NMD/ })).toHaveText('NMDPredicted')
  const score = page.locator('.sticky')
  await expect(score).toContainText('+9 pts')
  await expect(score).toContainText('Likely pathogenic')
  await expect(score).toContainText('98.8%')
})

test('missense at REVEL 0.95: PP3_Strong (4) + PM2 (1) = 5, a hot VUS at 81.2%', async ({ page }) => {
  // 350^(5/8) = 38.91; 3.891 / (3.791 + 1) = 0.812.
  const r = row(page, 'c.31C>T')
  await expect(r).toContainText('p.(Arg11Cys)')
  await expect(r.locator('.crit')).toHaveText(['PM2', 'PP3_Strong'])
  await expect(r).toContainText('+5')
  await expect(r).toContainText('hot')
  await openVariant(page, 'c.31C>T')
  await expect(page.locator('.sticky')).toContainText('81.2%')
})

test('an in-frame deletion is named at its 3\' position: c.5_7del, p.(Ala2del), 3 points', async ({ page }) => {
  // ATG GCT GAA: deleting GCT (c.4_6) equals deleting CTG (c.5_7), which is the HGVS 3' form. Protein M A E > M E. PM4 (2) + PM2 (1).
  const r = row(page, 'c.5_7del')
  await expect(r).toContainText('p.(Ala2del)')
  await expect(r.locator('.crit')).toHaveText(['PM2', 'PM4'])
  await expect(r).toContainText('+3')
  await expect(r).toContainText('tepid')
})

test('a donor +1 change skips exon 2: 60 nt in frame, 20 of 32 residues, PVS1_Strong', async ({ page }) => {
  await openVariant(page, 'c.75+1G>A')
  const facts = page.locator('dl.facts')
  await expect(facts).toContainText('Splice donor')
  await expect(facts).toContainText('Intron 2, 1 nt from the donor')
  await expect(facts).toContainText('62.5%')
  await expect(page.locator('.note').filter({ hasText: 'skip exon 2' })).toContainText('60 nt, in frame')
  await expect(page.locator('.splice-verdict')).toContainText('abolished')
  await expect(page.locator('.sticky')).toContainText('+5 pts')
})

test('a last-exon stop escapes decay: 1 of 32 residues is 3.1%, PVS1_Moderate, 3 points', async ({ page }) => {
  await openVariant(page, 'c.94G>T')
  await expect(page.locator('.var-title')).toContainText('p.(Glu32Ter)')
  await expect(page.locator('dl.facts')).toContainText('Predicted to escape')
  await expect(page.locator('dl.facts')).toContainText('3.1%')
  await expect(page.locator('.sticky')).toContainText('+3 pts')
})

test('a common synonymous change is stand-alone benign: BA1, BP4 and BP7 total -10', async ({ page }) => {
  const r = row(page, 'c.54T>C')
  await expect(r).toContainText('p.(Ala18=)')
  await expect(r.locator('.crit')).toHaveText(['BA1', 'BP4', 'BP7'])
  await expect(r).toContainText('-10')
  await openVariant(page, 'c.54T>C')
  await expect(page.locator('.sticky')).toContainText('BA1 is stand-alone')
})

test('recessive and in trans: 5 points become 7 with PM3_Moderate', async ({ page }) => {
  // c.11G>A is 0|1 and likely pathogenic; c.31C>T is 1|0, so the two are on different haplotypes.
  await importFile(page, 'phased.vcf', vcf([v(1026, 'G', 'A', 'AF_popmax=0', '0|1'), v(1086, 'C', 'T', 'AF_popmax=0;REVEL=0.95', '1|0')], 'P2'))
  await tab(page, 'Genes')
  await page.getByLabel('MINI inheritance').selectOption('AR')
  await tab(page, 'Worklist')
  const r = row(page, 'c.31C>T').filter({ hasText: '+7' })
  await expect(r.locator('.crit')).toHaveText(['PM2', 'PM3', 'PP3_Strong'])
  await expect(page.locator('.pairs')).toContainText('P2: c.11G>A LP in trans with c.31C>T')
})

test('the stored genomic position matches the paper coordinates', async ({ page }) => {
  // c.31 is base 16 of exon 2: 1071 + 15 = 1086.
  await openVariant(page, 'c.31C>T')
  await expect(page.locator('.var-head')).toContainText('chr9:1086 C>T')
  await expect(page.locator('.var-head')).toContainText('from mini.vcf, line 5')
})
