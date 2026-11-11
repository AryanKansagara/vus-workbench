import { useSyncExternalStore } from 'react'
import type { ClassName, GeneCfg } from './engine/acmg.ts'
import type { LabRow, VariantRec } from './engine/io.ts'
import type { Curation } from './engine/pipeline.ts'
import type { Genome, Transcript } from './engine/transcript.ts'

export interface ImportMeta { id: string; kind: string; fileName: string; rows: number; issues: number; at: number }

export interface SnapRow { key: string; label: string; cls: ClassName; points: number; criteria: string[] }
export interface Snapshot { id: string; label: string; createdAt: number; rows: SnapRow[] }

export interface Project {
  id: string
  name: string
  createdAt: number
  updatedAt: number
  variants: VariantRec[]
  labRows: LabRow[]
  /** Imported GTF/GFF3 models; the bundled synthetic models are always available. */
  transcripts: Transcript[]
  /** Imported reference regions (FASTA from samtools faidx). */
  genome: Genome
  genes: Record<string, GeneCfg>
  curation: Record<string, Curation>
  snapshots: Snapshot[]
  imports: ImportMeta[]
  sampleId?: string
  /** Set when the cloud copy was trimmed for size and this device lacks the full variant list. */
  partial?: boolean
}

// ---------- Tiny IndexedDB wrapper (falls back to memory in private modes) ----------

const DB = 'vus-workbench'
let dbp: Promise<IDBDatabase | null> | null = null
function db(): Promise<IDBDatabase | null> {
  return dbp ??= new Promise(res => {
    try {
      const req = indexedDB.open(DB, 1)
      req.onupgradeneeded = () => {
        req.result.createObjectStore('projects', { keyPath: 'id' })
        req.result.createObjectStore('data')
      }
      req.onsuccess = () => res(req.result)
      req.onerror = () => res(null)
    } catch { res(null) }
  })
}
const mem = { projects: new Map<string, unknown>(), data: new Map<string, unknown>() }

async function tx<T>(store: 'projects' | 'data', mode: IDBTransactionMode, fn: (s: IDBObjectStore) => IDBRequest<T>): Promise<T> {
  const d = await db()
  if (!d) throw new Error('no-idb')
  return new Promise((res, rej) => {
    const r = fn(d.transaction(store, mode).objectStore(store))
    r.onsuccess = () => res(r.result)
    r.onerror = () => rej(r.error)
  })
}

async function idbGetAll<T>(store: 'projects'): Promise<T[]> {
  try { return await tx(store, 'readonly', s => s.getAll()) as T[] } catch { return [...mem[store].values()] as T[] }
}
async function idbPut(store: 'projects' | 'data', value: unknown, key?: string) {
  try { await tx(store, 'readwrite', s => key === undefined ? s.put(value) : s.put(value, key)) } catch { mem[store].set(key ?? (value as { id: string }).id, value) }
}
async function idbDel(store: 'projects' | 'data', key: string) {
  try { await tx(store, 'readwrite', s => s.delete(key)) } catch { mem[store].delete(key) }
}

// ---------- In-memory project list with subscribers ----------

let projects: Project[] = []
let loaded = false
const subs = new Set<() => void>()
const emit = () => { projects = [...projects]; subs.forEach(f => f()) }
const changeHooks = new Set<(p: Project) => void>()

export const onProjectChange = (f: (p: Project) => void) => { changeHooks.add(f); return () => changeHooks.delete(f) }

export async function loadProjects() {
  projects = (await idbGetAll<Project>('projects')).sort((a, b) => b.updatedAt - a.updatedAt)
  loaded = true
  emit()
}

export function useProjects() {
  return useSyncExternalStore(f => { subs.add(f); return () => subs.delete(f) }, () => loaded ? projects : null)
}

export const getProjects = () => projects
export const getProject = (id: string) => projects.find(p => p.id === id)

export async function saveProject(p: Project, opts: { touch?: boolean; silent?: boolean } = {}) {
  const next = opts.touch === false ? p : { ...p, updatedAt: Date.now() }
  await idbPut('projects', next)
  const i = projects.findIndex(x => x.id === p.id)
  if (i >= 0) projects[i] = next
  else projects.unshift(next)
  projects.sort((a, b) => b.updatedAt - a.updatedAt)
  emit()
  if (!opts.silent) changeHooks.forEach(f => f(next))
  return next
}

/** Deletes locally; returns an undo function that restores the project. */
export async function deleteProject(id: string): Promise<() => Promise<void>> {
  const p = getProject(id)
  projects = projects.filter(x => x.id !== id)
  emit()
  await idbDel('projects', id)
  addTombstone(id)
  return async () => {
    if (!p) return
    removeTombstone(id)
    await saveProject(p)
  }
}

export const uuid = () => crypto.randomUUID()

// Tombstones keep a project deleted offline from coming back on the next sync.
const TOMB = 'vus-workbench.tombstones'
export const tombstones = (): string[] => { try { return JSON.parse(localStorage.getItem(TOMB) ?? '[]') } catch { return [] } }
const addTombstone = (id: string) => localStorage.setItem(TOMB, JSON.stringify([...new Set([...tombstones(), id])]))
const removeTombstone = (id: string) => localStorage.setItem(TOMB, JSON.stringify(tombstones().filter(x => x !== id)))
export const clearTombstones = (ids: string[]) => localStorage.setItem(TOMB, JSON.stringify(tombstones().filter(x => !ids.includes(x))))

// ---------- App defaults (settings page) ----------

export interface Defaults { editor: string; institution: string; inheritance: GeneCfg['inheritance']; pm2Max: number; bs1Min: number; ba1Min: number }
const DEF_KEY = 'vus-workbench.defaults'
const BASE_DEFAULTS: Defaults = { editor: '', institution: '', inheritance: 'AD', pm2Max: 0.0001, bs1Min: 0.01, ba1Min: 0.05 }
export function getDefaults(): Defaults {
  try { return { ...BASE_DEFAULTS, ...JSON.parse(localStorage.getItem(DEF_KEY) ?? '{}') } } catch { return BASE_DEFAULTS }
}
export const setDefaults = (d: Defaults) => localStorage.setItem(DEF_KEY, JSON.stringify(d))
export const resetDefaults = () => localStorage.removeItem(DEF_KEY)
export const editorName = () => getDefaults().editor.trim() || 'Unnamed curator'
export function defaultGene(): GeneCfg {
  const d = getDefaults()
  return { inheritance: d.inheritance, lofMechanism: true, pm2Max: d.pm2Max, bs1Min: d.bs1Min, ba1Min: d.ba1Min }
}
