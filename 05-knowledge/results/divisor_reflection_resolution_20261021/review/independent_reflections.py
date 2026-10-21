"""Independent set-edge audit; imports no producer implementation."""
from pathlib import Path
from itertools import combinations,permutations,product
from collections import Counter
import json,re,hashlib
BASE=Path(__file__).resolve().parent
PROOF=BASE.parent/'REFLECTION_RESOLUTION.md'
if not PROOF.exists():PROOF=BASE.parent/'session18-reflections'/'REFLECTION_RESOLUTION.md'
checks=Counter()
def require(p,k):
 if not p:raise RuntimeError(k)
 checks[k]+=1
def full_edges(m):return frozenset(combinations(range(m),2))
def turn_graph(g,m):return frozenset(tuple(sorted((m-1-a,m-1-b))) for a,b in g)
def graph(code,m):return frozenset(e for i,e in enumerate(combinations(range(m),2)) if code&(1<<i))
def C(x,k):
 a,b,c,e=x;return(k-a,1-b,1-c,e)
def Q(x):
 a,b,c,e=x;return(a,c,b,None if e is None else 1-e)
def src(k,resolved):
 out=set()
 for a,b,c in product(range(k+1),range(2),range(2)):
  if (a,b,c) in ((0,0,0),(k,1,1)):continue
  marks=(0,1) if resolved and k%2==0 and a*2==k and b+c==1 else (None,)
  out.update((a,b,c,e) for e in marks)
 return out
def groups(X,c,q):
 remain=set(X);o=[]
 while remain:
  x=next(iter(remain));s=frozenset((x,c(x),q(x),c(q(x))));require(s<=remain,'orbit_partition');remain-=s;o.append(s)
 return o
proof=PROOF.read_bytes();require(hashlib.sha256(proof).hexdigest()=='eecbe71f9cd36aba84ed3836d3f24c81b62ee50b7303eb430b4db6ac7f5da80c','frozen_proof')
rows=dict((a.strip(),int(b)) for a,b in re.findall(r'^\|([^|]+)\|\s*(\d+)\s*\|$',proof.decode(),re.M))
exps={'p':(1,0,0,None),'pqr':(1,1,1,None),'p^2':(2,0,0,None),'qr':(0,1,1,None),'q':(0,1,0,None),'r':(0,0,1,None),'p^2 r':(2,0,1,None),'p^2 q':(2,1,0,None),'(pq,0)':(1,1,0,0),'(pq,1)':(1,1,0,1),'(pr,0)':(1,0,1,0),'(pr,1)':(1,0,1,1)}
require(set(rows)==set(exps),'printed_source_labels');M={exps[s]:graph(code,4) for s,code in rows.items()};require(set(M)==src(2,True) and len(set(M.values()))==12,'printed_map_bijection')
E=full_edges(4)
def tournament(g,m):
 arcs=set()
 for i,j in combinations(range(m+1),2):
  arcs.add((j,i) if j>i+1 and (i,j-1) in g else (i,j))
 return frozenset(arcs)
blue=[]
for x,g in M.items():
 require(M[C(x,2)]==E-g,'printed_complement');require(M[Q(x)]==turn_graph(g,4),'printed_reflection');t=tournament(g,4)
 require(all((i,i+1) in t for i in range(4)),'directed_path')
 flipped=frozenset((b,a) if abs(a-b)>1 else (a,b) for a,b in t)
 require(tournament(E-g,4)==flipped,'actual_K')
 anti=frozenset((4-b,4-a) for a,b in t)
 require(tournament(turn_graph(g,4),4)==anti,'actual_R')
 require((t==anti)==(Q(x)==x),'actual_blue')
 if t==anti:blue.append(rows[next(n for n,v in exps.items() if v==x)])
require(len(blue)==4,'four_blue_states')
source_summary=[]
for k in range(2,51):
 X=src(k,False);chart={x:(2*x[0]-k*(x[1]+x[2]),x[1]-x[2]) for x in X};require(len(set(chart.values()))==len(X),'chart_injective')
 require(sum(C(x,k)==Q(x) for x in X)==(2 if k%2==0 else 0),'unresolved_CQ_fixed')
 require(sum(Q(x)==x for x in X)==2*k,'Q_fixed')
 Y=src(k,True)
 for x in Y:
  require(C(x,k) in Y and Q(x) in Y and C(Q(x),k)==Q(C(x,k)),'resolved_action')
  require(C(x,k)!=Q(x),'resolved_no_CQ_fixed')
  require(C(x,k)[:3]==C(x[:3]+(None,),k)[:3] and Q(x)[:3]==Q(x[:3]+(None,))[:3],'source_projection')
 o=groups(Y,lambda x:C(x,k),Q);require(Counter(map(len,o))=={2:k,4:k//2+1},'source_orbit_formula')
 require(len(Y)==4*k+2+(2 if k%2==0 else 0),'minimal_cardinality')
 source_summary.append({'k':k,'size':len(Y),'two_orbits':k,'free_orbits':k//2+1})
graph_summary=[];capacity_tests=[]
for m in range(2,7):
 E=full_edges(m);X={graph(code,m) for code in range(1<<len(E))};c=lambda g:E-g;q=lambda g:turn_graph(g,m)
 fixed_edges={e for e in E if tuple(sorted((m-1-e[0],m-1-e[1])))==e};require(len(fixed_edges)==m//2,'fixed_edge_count')
 for g in X:require(q(g)!=c(g),'no_complement_reflection_graph')
 orbs=groups(X,c,q);t=(len(E)+m//2)//2;actual=Counter(map(len,orbs));expect={2:2**(t-1),4:((1<<len(E))-(1<<t))//4};expect={s:n for s,n in expect.items() if n};require(actual==expect,'target_orbit_formula')
 for p in permutations(range(m)):
  if all(p[p[i]]==i for i in range(m)):
   require(any(tuple(sorted((p[a],p[b])))==(a,b) for a,b in E),'all_involutions_fixed_edge')
 if m<=5:
  pool={2:[o for o in orbs if len(o)==2],4:[o for o in orbs if len(o)==4]}
  for k in range(2,13):
   Y=src(k,True);ys=groups(Y,lambda x:C(x,k),Q);possible=k<=len(pool[2]) and k//2+1<=len(pool[4]);used=Counter();mapping={}
   if possible:
    for o in ys:
     x=next(iter(o));g=next(iter(pool[len(o)][used[len(o)]]));used[len(o)]+=1
     for xx,gg in [(x,g),(C(x,k),c(g)),(Q(x),q(g)),(C(Q(x),k),c(q(g)))]:
      require(xx not in mapping or mapping[xx]==gg,'capacity_well_defined');mapping[xx]=gg
    require(len(mapping)==len(Y) and len(set(mapping.values()))==len(Y),'capacity_constructed_injection')
    for x in Y:require(mapping[C(x,k)]==c(mapping[x]) and mapping[Q(x)]==q(mapping[x]),'capacity_equivariance')
   capacity_tests.append({'m':m,'k':k,'fits':possible})
 graph_summary.append({'m':m,'orbits':dict(sorted(actual.items())),'rho_fixed_graphs':sum(q(g)==g for g in X)})
print(json.dumps({'status':'PASS','method':'Independent set-edge, directed-arc and literal printed table implementation; no producer imports','proof_sha256':hashlib.sha256(proof).hexdigest(),'total_checks':sum(checks.values()),'checks':dict(sorted(checks.items())),'source_k_range':[2,50],'source_summary':source_summary,'graph_summary':graph_summary,'capacity_constructive_tests':capacity_tests,'printed_map':rows,'blue_codes':sorted(blue)},indent=2,sort_keys=True))
