#!/usr/bin/env python3
"""Independent brute-force certificates, no third-party dependencies."""
import itertools,json
from pathlib import Path

def matrix(n,mask):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  a[i][j]=(mask>>k)&1;a[j][i]=1-a[i][j]
 return a

def modules(a):
 n=len(a);out=[]
 for bits in range(1,1<<n):
  s=frozenset(i for i in range(n) if bits>>i&1)
  if all(len({a[x][y] for y in s})==1 for x in range(n) if x not in s):out.append(s)
 return out

def hierarchy(a):
 mods=modules(a)
 strong=[s for s in mods if all(not(s&t) or s<=t or t<=s for t in mods)]
 child={s:[t for t in strong if t<s and not any(t<u<s for u in strong)] for s in strong}
 return strong,child

def hp(a):
 n=len(a); d={}
 for v in range(n):d[1<<v,v]=1
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1 and s!=(1<<v):d[s,v]=sum(d.get((s^(1<<v),u),0)*a[u][v] for u in range(n) if s>>u&1 and u!=v)
 return sum(d.get(((1<<n)-1,v),0) for v in range(n))

def canonical(a):
 return min(tuple(a[p[i]][p[j]] for i in range(len(a)) for j in range(i+1,len(a))) for p in itertools.permutations(range(len(a))))

def verify(a):
 strong,ch=hierarchy(a);n=len(a);root=frozenset(range(n))
 # Euler depth RMQ implementation separate from least-containing-set oracle.
 e=[];depth=[];first={}
 def tour(s,d):
  first.setdefault(s,len(e));e.append(s);depth.append(d)
  for t in sorted(ch[s],key=lambda z:min(z)):
   tour(t,d+1);e.append(s);depth.append(d)
 tour(root,0)
 for u,v in itertools.combinations(range(n),2):
  l,r=sorted((first[frozenset([u])],first[frozenset([v])]))
  w=e[min(range(l,r+1),key=lambda k:depth[k])]
  assert w==min((s for s in strong if u in s and v in s),key=len)
  cu=next(t for t in ch[w] if u in t);cv=next(t for t in ch[w] if v in t)
  # Representative quotient lookup, independent of queried endpoints.
  assert a[u][v]==a[min(cu)][min(cv)]
 for s,children in ch.items():
  if not children:continue
  q=[[int(i!=j)*a[min(x)][min(y)] for j,y in enumerate(children)] for i,x in enumerate(children)]
  transitive=not any(q[i][j] and q[j][k] and q[k][i] for i,j,k in itertools.permutations(range(len(q)),3))
  prime=all(len(t) in (1,len(q)) for t in modules(q))
  assert transitive or prime
 return strong,ch

report={'universe':'All labelled tournaments n=1..5; lexicographic pair bit=1 means lower vertex beats higher.','checks':[]}
prime_examples={}
for n in range(1,6):
 for m in range(1<<(n*(n-1)//2)):
  a=matrix(n,m);st,ch=verify(a)
  if n==5 and len(modules(a))==6:
   h=hp(a)
   if h not in prime_examples:prime_examples[h]={'mask':m,'H':h,'matrix':a,'strong_modules':[sorted(s) for s in st]}
 report['checks'].append({'n':n,'labelled_count':1<<(n*(n-1)//2),'passed':True})
report['prime_star_H_examples']=list(prime_examples.values())
for ex in report['prime_star_H_examples']:
 a=ex['matrix']
 assert ex['H']==sum(all(a[p[i]][p[i+1]] for i in range(len(a)-1)) for p in itertools.permutations(range(len(a))))
# Fixed Hamilton order 0 -> 1 -> ... -> n-1, exhaustive until first H collision.
report['cut_flux_search']=[]
for n in range(1,8):
 free=[(i,j) for i in range(n) for j in range(i+1,n) if j!=i+1];seen={};witness=None
 for m in range(1<<len(free)):
  a=[[int(i<j) for j in range(n)] for i in range(n)]
  for k,(i,j) in enumerate(free):a[i][j]=(m>>k)&1;a[j][i]=1-a[i][j]
  flux=tuple(sum(2*a[i][j]-1 for i in range(k) for j in range(k,n)) for k in range(1,n))
  h=hp(a)
  if flux in seen and seen[flux]['H']!=h:
   witness={'n':n,'flux':flux,'first':seen[flux],'second':{'free_mask':m,'H':h,'matrix':a}}
   break
  seen.setdefault(flux,{'free_mask':m,'H':h,'matrix':a})
 report['cut_flux_search'].append({'n':n,'possible_tournaments':1<<len(free),'completed':witness is None,'collision':witness})
 if witness:
  for key in ('first','second'):
   a=witness[key]['matrix']; scores=[sum(row) for row in a]
   witness[key]['ordered_scores']=scores
   assert flux==tuple(sum(2*scores[i]-(n-1) for i in range(k)) for k in range(1,n))
   assert witness[key]['H']==sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in itertools.permutations(range(n)))
  def cart(l,r):
   if l>=r:return None
   k=min(range(l,r),key=lambda j:(flux[j],j))
   return [k,cart(l,k),cart(k+1,r)]
  witness['leftmost_min_cartesian_tree']=cart(0,len(flux))
  witness['HP_independent_permutation_check']=True
  break
# A000568 input through n=7, independently recover strong counts by SEQ inverse.
t=[1,1,1,2,4,12,56,456];s=[0]*len(t)
for n in range(1,len(t)):s[n]=t[n]-sum(s[k]*t[n-k] for k in range(1,n))
assert s==[0,1,0,1,1,6,35,353]
report['sequence_inverse']={'tournament_counts':t,'strong_counts':s,'status':'Arithmetic check of stated known inputs, not independent tournament census past n=5'}
p=Path(__file__).with_name('certificate.json');p.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
