import { lazy, Suspense, useEffect, useState } from 'react'
import { Landing } from './Landing.tsx'
import { Footer, Header, NotFound, Toasts } from './ui.tsx'

const Workspace = lazy(() => import('./Workspace.tsx'))
const ProjectView = lazy(() => import('./Project.tsx'))
const Report = lazy(() => import('./Report.tsx'))
const Docs = lazy(() => import('./Docs.tsx'))
const Settings = lazy(() => import('./Settings.tsx'))

export function usePath() {
  const [path, setPath] = useState(location.pathname + location.search)
  useEffect(() => {
    const onPop = () => setPath(location.pathname + location.search)
    addEventListener('popstate', onPop)
    return () => removeEventListener('popstate', onPop)
  }, [])
  return path
}

export type Route =
  | { page: 'landing' | 'docs' | 'settings' | 'workspace' | '404' }
  | { page: 'project'; id: string; tab: string; vid?: string }
  | { page: 'report'; id: string; vid?: string }

function route(full: string): Route {
  const path = full.split('?')[0].replace(/\/+$/, '') || '/'
  if (path === '/') return { page: 'landing' }
  if (path === '/docs') return { page: 'docs' }
  if (path === '/settings') return { page: 'settings' }
  if (path === '/app') return { page: 'workspace' }
  let m = path.match(/^\/app\/p\/([\w-]+)\/report(?:\/([\w-]+))?$/)
  if (m) return { page: 'report', id: m[1], vid: m[2] }
  m = path.match(/^\/app\/p\/([\w-]+)\/v\/([\w-]+)$/)
  if (m) return { page: 'project', id: m[1], tab: 'variant', vid: m[2] }
  m = path.match(/^\/app\/p\/([\w-]+)(?:\/(worklist|labs|import|history|genes))?$/)
  if (m) return { page: 'project', id: m[1], tab: m[2] ?? 'worklist' }
  return { page: '404' }
}

const Loading = () => <div className="boot" role="status"><span className="spinner" /> Loading</div>
const TITLES: Record<string, string> = { landing: 'VUS Workbench: variant interpretation you can audit', docs: 'Docs: getting variants in · VUS Workbench', settings: 'Settings · VUS Workbench', workspace: 'Worklists · VUS Workbench', '404': 'Not found · VUS Workbench' }

export function App() {
  const r = route(usePath())
  useEffect(() => { if (TITLES[r.page]) document.title = TITLES[r.page] }, [r.page])
  if (r.page === 'report') return <Suspense fallback={<Loading />}><Report id={r.id} vid={r.vid} /><Toasts /></Suspense>
  return (
    <>
      <a className="skip" href="#main">Skip to content</a>
      <Header />
      <main id="main" tabIndex={-1}>
        <Suspense fallback={<Loading />}>
          {r.page === 'landing' && <Landing />}
          {r.page === 'docs' && <Docs />}
          {r.page === 'settings' && <Settings />}
          {r.page === 'workspace' && <Workspace />}
          {r.page === 'project' && <ProjectView id={r.id} tab={r.tab} vid={r.vid} />}
          {r.page === '404' && <NotFound />}
        </Suspense>
      </main>
      <Footer />
      <Toasts />
    </>
  )
}
