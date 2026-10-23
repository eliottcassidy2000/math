#!/usr/bin/env python3
from pathlib import Path
from math import isqrt
from fractions import Fraction
from itertools import product,combinations
from collections import Counter,defaultdict
import json
D=Path(__file__).resolve().parent

def ck(x,m):
 if not x:raise RuntimeError(m)
def parents(n):return [s*s-n for s in range(isqrt(n)+1,isqrt(2*n-1)+1)]
def edges(n):return [(a,b) for b in range(1,n+1) for a in parents(b)]
def canon(p):return p if p[0]<p[-1] else p[::-1]
def valid(p):return len(p)==len(set(p)) and set(p)==set(range(1,len(p)+1)) and all(isqrt(a+b)**2==a+b for a,b in zip(p,p[1:]))
def edgeword(p):return [tuple(sorted(e)) for e in zip(p,p[1:])]
def rank(p,q):
 ds={}
 def root(x):
  ds.setdefault(x,x)
  if ds[x]!=x:ds[x]=root(ds[x])
  return ds[x]
 r=0
 for a,b in zip(p,q):
  a,b=root(a),root(b)
  if a!=b:ds[a]=b;r+=1
 return r
F=[0,1,1]
while F[-1]<10000:F.append(sum(F[-2:]))
def zeck(n):
 ids=[]
 for i in range(len(F)-1,1,-1):
  if F[i]<=n:n-=F[i];ids.append(i)
 return ids
N=512
anc={};harm={};sel={};dcounts=[]
for n in range(1,N+1):
 pp=parents(n)
 if not pp:
  ck(n in (1,2,4),'unexpected founder');anc[n]={n};harm[n]=tuple(Fraction(int(n==j)) for j in (4,2,1))
 else:
  anc[n]=set.union(*(anc[p] for p in pp));harm[n]=tuple(sum(harm[p][j] for p in pp)/len(pp) for j in range(3))
 ck(sum(harm[n])==1,'probability');ck({(4,2,1)[j] for j,v in enumerate(harm[n]) if v}==anc[n],'support')
 for rule in ('min','max'):
  sel.setdefault(rule,{})[n]=n if not pp else sel[rule][min(pp) if rule=='min' else max(pp)]
for S,c in [({12,4,5,11},4),({2,18,7,9,16},2),({10,6,3,1,8,17},1)]:
 ck(all(anc[n]=={c} for n in S),'island ancestry')
# Exact finite certificate for the all-size ancestry tail.
ck(anc[70]=={1,4},'last-exception witness')
ck(all(anc[n]=={1,2,4} for n in range(71,100)),'induction base71..99')
ck((100-71)**2>4*(2*100-1),'largest-parent bound at100')
# Exhaust all selectors through25, preserving shared ancestry correlations.
colorings=Counter();forests=1
states=[{1:1,2:2,4:4}]
for n in range(1,26):
 pp=parents(n)
 if not pp:continue
 forests*=len(pp)
 states=[dict(s,**{})|{n:s[p]} for s in states for p in pp]
for s in states:colorings[tuple(s[i] for i in range(1,26))]+=1
ck(sum(colorings.values())==forests,'selector count')
for n in range(1,26):
 for j,c in enumerate((4,2,1)):
  ck(Fraction(sum(w for v,w in colorings.items() if v[n-1]==c),forests)==harm[n][j],'harmonic law')
ck(not any(s[12]==1 and s[22]==4 for s in colorings),'joint hostile')
# Load inherited complete finite path family, validate each path and totals.
raw=json.loads((D/'inherited_paths.json').read_text())
rows={r['n']:[canon(p) for p in r['paths']] for r in raw['rows'] if 25<=r['n']<=32}
for n,ps in rows.items():
 ck(all(valid(p) for p in ps),'path validation');ck(len(set(map(tuple,ps)))==len(ps),'path duplicate')
E=edges(80);eidx={e:i for i,e in enumerate(E)}
train=rows[32]
P=[([p[0],p[-1]],p,[eidx[e] for e in edgeword(p)]) for p in train]
extra=json.loads((D/'old_color_census.json').read_text())[-1]['first_collision']['paths']
host=json.loads((D/'hostile80.json').read_text())['paths']
ck(all(valid(p) for p in extra+host),'holdout valid')
# Each family is a prespecified finite linear-feature grammar over F2/F3.
def vf(n):
 z=zeck(n);return (n,n.bit_count(),isqrt(n),len(z),min(z),sum(z))
def ef(e):
 a,b=e;k=isqrt(a+b)
 return (a,b,k,a.bit_count()+b.bit_count(),isqrt(a)+isqrt(b),len(zeck(a))+len(zeck(b)),min(zeck(a))+min(zeck(b)))
features={'vertex':[vf(n) for n in range(1,81)],'edge':[ef(e) for e in E]}
summary={};best=[]
for kind,fs in features.items():
 dim=len(fs[0])
 for q in (2,3):
  seen=set();passed=[];maxlabels=0;hist=Counter();bestentry=None
  for coeff in product(range(q),repeat=dim):
   if not any(coeff):continue
   colors=tuple(sum(a*b for a,b in zip(coeff,f))%q for f in fs)
   # Canonicalize global color names for genuine distinct tested colorings.
   ren={};cc=tuple(ren.setdefault(x,len(ren)) for x in colors)
   if cc in seen:continue
   seen.add(cc)
   sig=[]
   for endpoints,p,es in P:
    word=tuple(colors[x-1] for x in p) if kind=='vertex' else tuple(colors[x] for x in es)
    sig.append((tuple(endpoints),word))
   count=len(set(sig));hist[count]+=1
   if count>maxlabels:maxlabels=count;bestentry=coeff
   if count==len(P):
    def distinguish(pair):
     if kind=='vertex':return tuple(colors[x-1] for x in pair[0])!=tuple(colors[x-1] for x in pair[1])
     return tuple(colors[eidx[e]] for e in edgeword(pair[0]))!=tuple(colors[eidx[e]] for e in edgeword(pair[1]))
    passed.append((coeff,distinguish(extra),distinguish(host)))
  summary[f'{kind}{q}']=dict(raw_nonzero=q**dim-1,distinct_up_to_color_names=len(seen),full_order32_decoders=len(passed),maximum_labels=maxlabels,best_coefficients=bestentry,holdout36_pair_pass=sum(x[1] for x in passed),holdout80_pair_pass=sum(x[2] for x in passed),both_holdout_pairs_pass=sum(x[1] and x[2] for x in passed),label_count_histogram=dict(sorted(hist.items())))
  best.extend(dict(kind=kind,q=q,coefficients=c,pass36=p36,pass80=p80) for c,p36,p80 in passed[:3])
# Exact collision probability ranks, conditioned on actual endpoint pair.
prob={}
for n,ps in rows.items():
 groups=defaultdict(list)
 for p in ps:groups[(p[0],p[-1])].append(p)
 for kind in ('vertex','edge'):
  ranks=Counter(rank(p if kind=='vertex' else edgeword(p),q if kind=='vertex' else edgeword(q)) for group in groups.values() for p,q in combinations(group,2))
  prob[f'{n}-{kind}']=dict(ranks=dict(sorted(ranks.items())),union_bounds={str(q):str(sum(Fraction(v,q**r) for r,v in ranks.items())) for q in (2,3)})
# First odd cycle via BFS bipartition, and exact first triangle.
first_triangle=next((n,[a,b,n]) for n in range(1,81) for a,b in combinations(parents(n),2) if isqrt(a+b)**2==a+b)
# First non-bipartite graph can precede first triangle.
nonbip=None
for n in range(1,81):
 G={i:[] for i in range(1,n+1)}
 for a,b in edges(n):G[a].append(b);G[b].append(a)
 col={};ok=True
 for s in G:
  if s in col:continue
  col[s]=0;todo=[s]
  for v in todo:
   for w in G[v]:
    if w not in col:col[w]=1-col[v];todo.append(w)
    elif col[w]==col[v]:ok=False
 if not ok:nonbip=n;break
ck(nonbip==15 and first_triangle[0]==30,'odd-cycle horizons')
cycle=[6,3,1,15,10]
ck(all(isqrt(a+b)**2==a+b for a,b in zip(cycle,cycle[1:]+cycle[:1])),'cycle witness')
out=dict(status='PASS',scope='Inherited complete path carrier25..32; pair-only heldouts36 and80; generated finite feature grammar, not all colorings.',birth=dict(founders=[1,2,4],first_two_ancestries=13,first_three_ancestries=22,all_three_tail_start=71,tail_induction_base=[71,99],last_exception70=sorted(anc[70]),parent_selector_forests25=forests,distinct_ancestry_colorings25=len(colorings),pure_counts512={str(k):sum(v=={k} for v in anc.values()) for k in (1,2,4)},mixed_counts512=dict(Counter(len(v) for v in anc.values())),harmonic_examples={str(n):list(map(str,harm[n])) for n in (13,14,15,22,23,25,32,64,128,512)},joint_impossible='color(13)=founder1 and color(23)=founder4'),first_triangle=first_triangle,first_nonbipartite=nonbip,odd_cycle15=[6,3,1,15,10],color_search=summary,example_successes=best,collision_probabilities=prob)
print(json.dumps(out,indent=2))
