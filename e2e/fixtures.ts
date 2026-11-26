import { test as base, expect, type Locator, type Page } from '@playwright/test'
import { readFile } from 'node:fs/promises'

const FONT_HOSTS = new Set(['fonts.googleapis.com', 'fonts.gstatic.com'])
/** The only services the README says the app contacts, and only when asked. Tests that need one add a mock on top. */
const ON_REQUEST = new Set(['rest.ensembl.org', 'eutils.ncbi.nlm.nih.gov'])

export const test = base.extend<{ external: string[]; onRequest: string[] }>({
  external: async ({}, use) => { await use([]) },
  /** Calls to the two documented services. A test that did not click a lookup expects none. */
  onRequest: async ({}, use) => { await use([]) },

  // Nothing leaves the machine. A request to any other host fails the test, and so does any uncaught page error.
  context: async ({ context, baseURL, external, onRequest }, use) => {
    const errors: string[] = []
    const origin = new URL(baseURL!).origin
    await context.route('**/*', (route) => {
      const url = new URL(route.request().url())
      if (url.origin === origin) return route.continue()
      if (FONT_HOSTS.has(url.host)) return route.fulfill({ status: 200, contentType: 'text/css', body: '' })
      ;(ON_REQUEST.has(url.host) ? onRequest : external).push(url.host)
      return route.abort('blockedbyclient')
    })
    context.on('weberror', (e) => errors.push(String(e.error())))
    await use(context)
    expect(errors, 'uncaught page errors').toEqual([])
    expect(external, 'requests to undocumented hosts').toEqual([])
  },
})

export { expect }

export const file = (name: string, data: string) => ({ name, mimeType: 'text/plain', buffer: Buffer.from(data) })
export const row = (page: Page, text: string | RegExp) => page.locator('table.worklist tbody tr').filter({ hasText: text })
export const pending = (page: Page, name: string) => page.locator('.pending').filter({ hasText: name })
export const issues = (page: Page, name: string) => pending(page, name).locator('.issues')
export const importBtn = (page: Page, name: string) => pending(page, name).getByRole('button', { name: /^Import \d/ })
export const tab = (page: Page, name: string) => page.getByRole('navigation', { name: 'Worklist sections' }).getByRole('link', { name }).click()
export const toastText = (page: Page, text: string | RegExp) => expect(page.locator('.toast').filter({ hasText: text })).toBeVisible()

/** Creates an empty worklist and stops on its Import tab. */
export async function newWorklist(page: Page, name = 'E2E case') {
  await page.goto('/app')
  await page.getByRole('button', { name: 'New worklist' }).click()
  await page.getByLabel('Case, panel or batch name').fill(name)
  await page.getByRole('button', { name: 'Create and import' }).click()
  await expect(page.getByRole('heading', { name: 'Drop files here' })).toBeVisible()
}

export const choose = (page: Page, name: string, data: string) => page.locator('input[type=file]').setInputFiles([file(name, data)])

/** Imports one file from wherever the worklist is open. */
export async function importFile(page: Page, name: string, data: string) {
  if (!(await page.getByRole('heading', { name: 'Drop files here' }).isVisible())) await tab(page, 'Import')
  await choose(page, name, data)
  await importBtn(page, name).click()
  await expect(pending(page, name)).toHaveCount(0)
}

/** An empty worklist with the hand-made MINI gene and its sequence loaded. */
export async function withMini(page: Page) {
  await newWorklist(page)
  await importFile(page, 'mini.fa', MINI_FASTA)
  await importFile(page, 'mini.gtf', MINI_GTF)
}

/** A worklist of MINI variants given by c. description. Popmax 0 unless stated. */
export async function miniHgvs(page: Page, lines: string[], header = 'gene,hgvs_c,popmax_af') {
  await withMini(page)
  await importFile(page, 'list.csv', [header, ...lines, ''].join('\n'))
  await expect(page.locator('table.worklist')).toBeVisible()
}

export async function openVariant(page: Page, c: string | RegExp) {
  await row(page, c).getByRole('link').first().click()
  await expect(page.locator('.var-title')).toBeVisible()
}

export async function saveAs(page: Page, trigger: Locator): Promise<[string, string]> {
  const [d] = await Promise.all([page.waitForEvent('download'), trigger.click()])
  return [d.suggestedFilename(), await readFile((await d.path())!, 'utf8')]
}

/** CSV text to rows of cells. */
export function cells(csv: string): string[][] {
  return csv.replace(/\r?\n$/, '').split(/\r?\n/).map((l) => [...l.matchAll(/(?:^|,)(?:"((?:[^"]|"")*)"|([^,]*))/g)].map((m) => m[1]?.replace(/""/g, '"') ?? m[2]))
}

// ---------- A gene small enough to work out on paper ----------
// Contig chr9 from 1001. Exon 1 = 1011-1030 (5 nt 5' UTR, then c.1-15), exon 2 = 1071-1130 (c.16-75),
// exon 3 = 1171-1200 (c.76-99 with the stop, then 6 nt 3' UTR). Both introns are 40 nt. 99 / 3 - 1 = 32 amino acids.
//   codons 1-5   ATG GCT GAA TGG CAA
//   codons 6-25  GAT CTG AAA GGT ACC CGT TCT GTT CAT AAC CCA ATG GCT TAT GAA CTG GAT AAA GGT TTC
//   codons 26-33 CGT GAT TCT CTG GCT CCA GAA TAA
export const MINI_SEQ = [
  'TTTTTTTTTT', 'CCACC', 'ATGGCTGAATGGCAA',
  'GTAAGTATTTATTTATTTACTAACTTTTCTTTCTTTTCAG',
  'GATCTGAAAGGTACCCGTTCTGTTCATAACCCAATGGCTTATGAACTGGATAAAGGTTTC',
  'GTAAGTTTTAATTTAATTTACTAACTTTTCTCTCTTTCAG',
  'CGTGATTCTCTGGCTCCAGAATAA', 'TTGCCT', 'AAAAAAAAAA',
].join('')
export const MINI_FASTA = `>chr9:1001-1210\n${MINI_SEQ}\n`
const gtf = (type: string, a: number, b: number) => `chr9\ttest\t${type}\t${a}\t${b}\t.\t+\t.\tgene_id "MINI"; transcript_id "NM_900001.1"; gene_name "MINI";`
const FEATURES: [string, number, number][] = [['exon', 1011, 1030], ['exon', 1071, 1130], ['exon', 1171, 1200], ['CDS', 1016, 1030], ['CDS', 1071, 1130], ['CDS', 1171, 1191], ['stop_codon', 1192, 1194]]
export const MINI_GTF = FEATURES.map((f) => gtf(...f)).join('\n') + '\n'
export const MINI_GFF3 = ['##gff-version 3', 'chr9\ttest\tgene\t1011\t1200\t.\t+\t.\tID=g1;Name=MINI', 'chr9\ttest\tmRNA\t1011\t1200\t.\t+\t.\tID=NM_900001.1;Parent=g1',
  ...FEATURES.filter(([t]) => t !== 'stop_codon').map(([t, a, b]) => `chr9\ttest\t${t}\t${a}\t${t === 'CDS' && b === 1191 ? 1194 : b}\t.\t+\t.\tParent=NM_900001.1`)].join('\n') + '\n'

export const vcf = (rows: string[], sample = 'P1') => ['##fileformat=VCFv4.2', '##INFO=<ID=AF_popmax,Number=A,Type=Float,Description="x">', `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sample}`, ...rows, ''].join('\n')
export const v = (pos: number, ref: string, alt: string, info = 'AF_popmax=0', gt = '0/1') => `chr9\t${pos}\t.\t${ref}\t${alt}\t60\tPASS\t${info}\tGT\t${gt}`
