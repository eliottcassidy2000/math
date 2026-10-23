"""Independent adjacency-array canonicalization and subset-DP path counts."""
import itertools as it, json, hashlib, os, subprocess
from pathlib import Path
import sympy as s

def need(x,msg):
 if not x: raise RuntimeError(msg)

def classify(n,fix=0):
 edges=list(it.combinations(range(n),2));perms=[p for p in it.permutations(range(n)) if all(p[i]==i for i in range(fix))]
 graphs={};tournaments={};pathcounts={}
 for x in range(1<<len(edges)):
  graph=[[0]*n for _ in range(n)];arcs=[[0]*n for _ in range(n)]
  for k,(i,j) in enumerate(edges):
   b=(x>>k)&1;graph[i][j]=graph[j][i]=b;arcs[i][j]=1-b;arcs[j][i]=b
  graphs[x]=min(sum(graph[p[i]][p[j]]<<k for k,(i,j) in enumerate(edges)) for p in perms)
  tournaments[x]=min(sum((1-arcs[p[i]][p[j]])<<k for k,(i,j) in enumerate(edges)) for p in perms)
  dp={ (1<<i,i):1 for i in range(n)}
  for mask in range(1,1<<n):
   for last in range(n):
    count=dp.get((mask,last),0)
    if count:
     for nxt in range(n):
      if not(mask>>nxt&1) and arcs[last][nxt]:
       key=(mask|(1<<nxt),nxt);dp[key]=dp.get(key,0)+count
  pathcounts[x]=sum(dp.get(((1<<n)-1,i),0) for i in range(n))
 return edges,graphs,tournaments,pathcounts

def certificate(words,g,t,paired=False):
 full=max(g)
 if paired:
  g={x:min(g[x],g[full^x]) for x in g};t={x:min(t[x],t[full^x]) for x in t}
 gs=sorted({g[x] for x in words});ts=sorted({t[x] for x in words})
 C=s.zeros(len(gs),len(ts));gi={x:i for i,x in enumerate(gs)};ti={x:i for i,x in enumerate(ts)}
 for x in words:C[gi[g[x]],ti[t[x]]]+=1
 r=[sum(C.row(i)) for i in range(C.rows)];c=[sum(C.col(j)) for j in range(C.cols)]
 B=s.diag(*[1/s.Integer(v) for v in c])*C.T*s.diag(*[1/s.Integer(v) for v in r])*C
 # AB and BA have the same nonzero spectrum; pad smaller-side polynomial.
 z=s.Symbol('lambda');poly=s.cancel(B.charpoly().as_expr()*z**(C.rows-C.cols))
 need(s.denom(poly)==1,'padded polynomial')
 return {'incidence':[[int(v) for v in C.row(i)] for i in range(C.rows)],'word_count':len(words),'charpoly':str(s.factor(poly)),'graph_states':len(gs),'tournament_states':len(ts)}

HERE=Path(__file__).resolve().parent
source=next((p for p in (HERE.parent/'coupling',HERE.parent/'session43-coupling') if (p/'MANIFEST.json').is_file()),None)
need(source is not None,'locate sibling coupling packet')
m=json.loads((source/'MANIFEST.json').read_text())
for f in m['files']:
 b=(source/f['path']).read_bytes();need(len(b)==f['bytes'] and hashlib.sha256(b).hexdigest()==f['sha256'],'manifest '+f['path'])
provenance=json.loads((source/'SOURCE_PROVENANCE.json').read_text())
repo_candidates=[]
if os.environ.get('MATH_REPOSITORY'):repo_candidates.append(Path(os.environ['MATH_REPOSITORY']))
for parent in HERE.parents:
 repo_candidates.extend([parent,parent/'session35-publication/clean-repo',parent/'math'])
repo_candidates=[p for p in repo_candidates if (p/'.git').exists()]
for f in provenance['sources']:
 # Pin the inherited source, rather than silently accepting a later working copy.
 suffix=f['path'].removeprefix('session42-')
 git_path='05-knowledge/results/session42_shared_lattices_square_colors_20261023/'+suffix
 data=None
 for repo in repo_candidates:
  result=subprocess.run(['git','-C',str(repo),'show',provenance['parent_publication']+':'+git_path],capture_output=True)
  if result.returncode==0:data=result.stdout;break
 if data is None:
  # Also permit a verified extracted source layout with no Git checkout.
  choices=[parent/f['path'] for parent in HERE.parents]
  choices += [parent/git_path for parent in HERE.parents]
  candidate=next((p for p in choices if p.is_file()),None)
  need(candidate is not None,'locate pinned provenance; set MATH_REPOSITORY to the repository')
  data=candidate.read_bytes()
 need(hashlib.sha256(data).hexdigest()==f['sha256'],'source provenance')
records=json.loads((source/'receipt.json').read_text())['reports'];out=[]
for n in range(2,6):
 E,G,T,H=classify(n);path=sum(1<<E.index((i,i+1)) for i in range(n-1))
 for r in records:
  if r['n']!=n or r['model']=='individually fixed vertices':continue
  model=r['model'];Tuse=T
  if model=='Hamiltonian observer':words=[x for x in G if x&path==0]
  elif model=='modulo-8 conditioning':words=[x for x in G if H[x]%8==r['residue']]
  elif model=='modulo-8 observed coarsening':words=list(G);Tuse={x:H[x]%8 for x in G}
  elif model.startswith('soft '):
   accept=(lambda x:x&path==0) if model=='soft observer' else (lambda x:H[x]%8==7)
   words=[x for x in G for _ in range(2 if accept(x) else 1)]
  else:raise RuntimeError('unhandled model')
  result=certificate(words,G,Tuse,r['paired'])
  for key in result:need(result[key]==r[key],f'independent {n} {model} {key}')
  out.append({'n':n,'model':model,'residue':r.get('residue'),'paired':r['paired'],'word_count':len(words),'exact_charpoly':result['charpoly']})
for r in records:
 if r['model']!='individually fixed vertices':continue
 E,G,T,H=classify(r['n'],r['fixed_vertices']);result=certificate(list(G),G,T,r['paired'])
 for key in result:need(result[key]==r[key],'independent marked '+key)
 out.append({'n':r['n'],'model':r['model'],'fixed_vertices':r['fixed_vertices'],'paired':r['paired'],'exact_charpoly':result['charpoly']})
print(json.dumps({'status':'PASS','independent_method':'adjacency-array permutation canonicalization, subset DP Hamilton counts, smaller Gram-side exact characteristic polynomial','reports_checked':len(out),'reports':out},sort_keys=True,indent=2))
