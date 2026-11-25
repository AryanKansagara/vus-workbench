// Delimited-text parsing shared by the variant-list, lab-report and GTF importers.

export interface Csv { rows: string[][]; lines: number[]; /** Rows left out because their quotes do not close properly. */ bad?: { line: number; message: string }[] }

export function parseCsv(text: string): Csv {
  if (text.charCodeAt(0) === 0xfeff) text = text.slice(1)
  const delim = sniffDelimiter(text)
  const rows: string[][] = []
  const bad: NonNullable<Csv['bad']> = []
  // Line of the quote that opened the last quoted cell, the line to report when text follows its closing quote, and whether it is closed.
  let openLine = 0, stray = 0, closed = false
  const lines: number[] = []
  let row: string[] = []
  let field = ''
  let quoted = false
  let line = 1
  let rowLine = 1
  // Without a readable header row, the columns of the rows below it are unknown.
  let headerBad = false
  const end = () => {
    if (stray) { bad.push({ line: stray, message: `Text after the closing quote of a cell that starts on this line${line > rowLine ? ` (the row ends on line ${line})` : ''}; row skipped. Is a closing quote missing?` }); if (!rows.length) headerBad = true }
    else if (row.length > 1 || row[0] !== '') { rows.push(row); lines.push(rowLine) }
    row = []; stray = 0
  }
  for (let i = 0; i < text.length; i++) {
    const c = text[i]
    if (quoted) {
      if (c === '"') {
        if (text[i + 1] === '"') { field += '"'; i++ } else { quoted = false; closed = true }
      } else {
        if (c === '\n') line++
        field += c
      }
      continue
    }
    if (c === '"' && field === '' && !closed) { quoted = true; openLine = line }
    else if (c === delim) { row.push(field); field = ''; closed = false }
    else if (c === '\n' || c === '\r') {
      if (c === '\r' && text[i + 1] === '\n') i++
      row.push(field); field = ''; closed = false
      end()
      line++
      rowLine = line
    } else {
      // A cell such as "Foo" syndrome hides nothing; one that took in a delimiter or a line break may have swallowed others.
      if (closed && c !== ' ' && (field.includes(delim) || /[\r\n]/.test(field))) stray ||= openLine
      field += c
    }
  }
  if (quoted) {
    bad.push({ line: openLine, message: `A quote opened on this line is never closed; the row${line > rowLine ? ` and everything after it, up to line ${line},` : ''} was skipped` })
    if (!rows.length) headerBad = true
  } else if (field !== '' || row.length) { row.push(field); end() }
  return headerBad ? { rows: [], lines: [], bad } : { rows, lines, ...(bad.length ? { bad } : {}) }
}

function sniffDelimiter(text: string): string {
  const head = text.slice(0, text.indexOf('\n') === -1 ? text.length : text.indexOf('\n'))
  const counts = [',', '\t', ';', '|'].map(d => [d, head.split(d).length] as const)
  return counts.sort((a, b) => b[1] - a[1])[0][0]
}

/** A plain decimal or scientific number. A decimal comma is read when nothing else fits ("0,0001"); grouped or ambiguous forms ("1,000") and other notations are not numbers here. */
export function parseNum(raw: string): number {
  const s = raw.trim().replace(/^([+-]?)0,(\d+(?:e[+-]?\d+)?)$/i, '$10.$2')
  return /^[+-]?(\d+\.?\d*|\.\d+)(e[+-]?\d+)?$/i.test(s) ? Number(s) : NaN
}

/** Text a spreadsheet would run as a formula (OWASP CSV injection) is exported with a leading apostrophe. Plain numbers are not. */
const inert = (s: string) => /^[=+\-@\t\r]/.test(s) && Number.isNaN(parseNum(s)) ? `'${s}` : s

/** Takes the apostrophe off again, so an exported table reads back as it was. */
export const uninert = (s: string) => s[0] === "'" && inert(s.slice(1)) === s ? s.slice(1) : s

export function toCsv(rows: (string | number)[][]): string {
  return rows.map(r => r.map(v => {
    const s = typeof v === 'number' ? String(Math.round(v * 1e6) / 1e6) : inert(v)
    return /[",\r\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s
  }).join(',')).join('\n') + '\n'
}


export const norm = (h: string) => h.trim().toLowerCase().replace(/^"|"$/g, '').replace(/[\s\-.]+/g, '_')
