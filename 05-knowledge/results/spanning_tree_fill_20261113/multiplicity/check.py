"""Exact SCC-padding forests and permanent-mark counts, reusing the inherited atom atlas."""
from pathlib import Path
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from math import comb, factorial
import csv, json, subprocess
HERE=Path(__file__).resolve().parent
SOURCE_REF='dfaf7c7e88dc191ad5db42ccf6ff7b0c9c0346f9'
BANK='05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/extended_census.json'
OLD='05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/multiplicity/all_order_H_le375.json'

def need(p,msg):
 if not p:raise RuntimeError(msg)
def load(path):return json.loads(subprocess.check_output(['git','show',SOURCE_REF+':'+path],text=True))
def falling(n,k):
 out=1
 for j in range(k):out*=n-j
 return out

def auts(T):
 cells={}
 for v,row in enumerate(T):cells.setdefault(row.bit_count(),[]).append(v)
 order=list(cells.values());out=[]
 for blocks in product(*(list(permutations(c)) for c in order)):
  p=list(range(len(T)))
  for cell,im in zip(order,blocks):
   for v,w in zip(cell,im):p[v]=w
  if all(((T[i]>>j)&1)==((T[p[i]]>>p[j])&1) for i in range(len(T)) for j in range(len(T))):out.append(tuple(p))
 return out

def tables(T):
 n=len(T);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for m in range(1,1<<n):
  for v in range(n):
   if m>>v&1:
    rest=m^(1<<v)
    if rest:dp[m][v]=sum(dp[rest][u] for u in range(n) if rest>>u&1 and T[u]>>v&1)
 return [1]+[sum(z) for z in dp[1:]],dp

def opposite(T):
 full=(1<<len(T))-1
 return tuple(full^(1<<i)^row for i,row in enumerate(T))

def gaps(t,d):
 if d==1:yield (t,);return
 for a in range(t+1):
  for tail in gaps(t-a,d-1):yield (a,)+tail

def parent(g):
 j=next(i for i,x in enumerate(g) if x)
 out=list(g);out[j]-=1
 return tuple(out)

def children(g):
 j=next((i for i,x in enumerate(g) if x),len(g)-1)
 for i in range(j+1):
  h=list(g);h[i]+=1;yield tuple(h)

def main():
 raw=load(BANK);old=load(OLD);atoms=[]
 for rec in raw['records']:
  if not rec['strong'] or rec['H']>27:continue
  T=tuple(rec['adjacency_outmasks'])
  for flip in range(rec['ordinary_multiplicity']):
   U=T if not flip else opposite(T);A=auts(U)
   need(len(A)==rec['aut'],'inherited automorphism order')
   hp,_=tables(U);_,starts=tables(opposite(U));full=(1<<len(U))-1
   eligible=[r for r in range(len(U)) if sum(row>>r&1 for row in U)==1 and hp[full^(1<<r)]-starts[-1][r]==1]
   ready_orbits={min(a[r] for a in A) for r in eligible}
   atoms.append(dict(id=rec['node']+(':op' if flip else ''),H=rec['H'],n=len(U),rows=U,
     aut=len(A),fix=Counter(sum(i==v for i,v in enumerate(a)) for a in A),ready=len(ready_orbits)))
 words={1:[()]}
 def get(h):
  if h not in words:
   words[h]=[(i,)+w for i,a in enumerate(atoms) if h%a['H']==0 for w in get(h//a['H'])]
  return words[h]
 forest_checks=0
 for r in range(5):
  for t in range(7):
   layer=list(gaps(t,r+1));nxt=list(gaps(t+1,r+1));edges=[z for g in layer for z in children(g)]
   need(len(layer)==comb(t+r,r),'layer binomial')
   need(len(edges)==len(set(edges))==len(nxt),'canonical birth collision-free and complete')
   need(all(parent(z) in layer for z in nxt),'unique parent')
   if t:
    branch=sum(sum(1 for _ in children(g))>1 for g in layer)
    need(branch==(comb(t+r-1,r-1) if r else 0),'branching count')
   # unrestricted gap insertion has exactly #positive gaps incoming arcs
   indeg=Counter()
   for g in layer:
    for i in range(r+1):
     q=list(g);q[i]+=1;indeg[tuple(q)]+=1
   need(all(indeg[z]==sum(v>0 for v in z) for z in nxt),'unpruned collision multiplicity')
   forest_checks+=1
 rows=[];marked_checks=0
 for h in range(1,28,2):
  descriptors=[]
  for W in get(h):
   m=sum(atoms[i]['n'] for i in W);dist=Counter({0:1});a=1
   for i in W:
    a*=atoms[i]['aut'];nxt=Counter()
    for f,c in dist.items():
     for g,d in atoms[i]['fix'].items():nxt[f+g]+=c*d
    dist=nxt
   descriptors.append((W,m,a,dist))
  for n in range(1,15):
   N=0;mass=Fraction(0);marked=[Fraction(0)]*3
   for W,m,a,dist in descriptors:
    t=n-m;r=len(W)
    if t<0:continue
    B=comb(t+r,r);N+=B;mass+=Fraction(B,a)
    for k in range(3):
     direct=Fraction(sum(c*falling(t+f,k) for f,c in dist.items()),a)*B
     anchor=Fraction(0)
     for j in range(k+1):
      if t<j:continue
      core=Fraction(sum(c*falling(f,k-j) for f,c in dist.items()),a)
      anchor+=comb(k,j)*factorial(r+j)//factorial(r)*core*comb(t+r,r+j)
     need(direct==anchor,'permanent-anchor positive formula')
     marked[k]+=direct;marked_checks+=1
   inherited=old[str(h)]['core_words']
   Nold=sum(int(v['ordinary_coefficient'])*comb(n-v['total_size']+v['strong_components'],v['strong_components']) for v in inherited if n>=v['total_size'])
   need(N==Nold,'inherited complete fixed-H ordinary count')
   need(all(v.denominator==1 for v in marked),'marked class integrality')
   labelled=mass*factorial(n);need(labelled.denominator==1,'labelled class integrality')
   ready=sum(a['ready'] for a in atoms if a['H']==h and a['n']<=n)
   if h==1:ready=int(n>=2) # second vertex of the transitive tournament
   rows.append(dict(n=n,H=h,ordinary=N,labelled=int(labelled),rooted=int(marked[1]),two_ordered_marks=int(marked[2]),unit_ready_roots=ready))
 # Independent mark-orbit check on inherited representatives through5.
 # Converse partners are not regenerated for this small auxiliary control.
 orbit_checks=0
 for rec in raw['records']:
  if rec['n']>5:continue
  T=tuple(rec['adjacency_outmasks']);A=auts(T);n=len(T)
  for k in (1,2):
   reps={min(tuple(a[v] for v in mark) for a in A) for mark in permutations(range(n),k)}
   burn=sum(falling(sum(i==v for i,v in enumerate(a)),k) for a in A)//len(A)
   need(len(reps)==burn,'direct ordered-mark orbits');orbit_checks+=1
 (HERE/'growth_table_h_le27.csv').write_text('')
 with (HERE/'growth_table_h_le27.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
 out=dict(status='PASS',source_ref=SOURCE_REF,scope='Inherited globally complete strong atoms H<=27; ordinary not converse-merged; n1..14;k0,1,2 ordered distinct marks',
  atom_count=len(atoms),forest_controls=forest_checks,anchor_identity_controls=marked_checks,direct_mark_orbit_controls=orbit_checks,
  atoms=[{**a,'fix':dict(a['fix'])} for a in atoms],
  strong_ready_counts={h:sum(a['ready'] for a in atoms if a['H']==h) for h in range(3,28,2)},
  words_per_H={h:len(get(h)) for h in range(1,28,2)},table_rows=len(rows),
  h9_examples=[r for r in rows if r['H']==9 and r['n'] in(5,6,8,14)])
 (HERE/'receipt.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:out[k] for k in ('status','atom_count','forest_controls','anchor_identity_controls','direct_mark_orbit_controls','strong_ready_counts','h9_examples')},indent=2))

if __name__=='__main__':main()
