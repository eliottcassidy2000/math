"""Independent exact source-cell, matching, and dual audit. No producer import."""
from pathlib import Path
from itertools import combinations,product,permutations
from collections import Counter
from functools import lru_cache
from fractions import Fraction as Q
import json,sys
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve();OUT=HERE.parent.parent/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
gates=0
def need(ok,msg):
 global gates
 gates+=1
 if not ok:raise RuntimeError(msg)
def boards(n):
 for rows in product(tuple(combinations(range(n),2)),repeat=n):
  deg=Counter(c for pair in rows for c in pair)
  if all(deg[c]==2 for c in range(n)):yield rows
def points(rows):return tuple((r,c) for r,pair in enumerate(rows) for c in pair)
def lines(pts,slopes):
 ans=[]
 for a in slopes:
  bank={}
  for i,(r,c) in enumerate(pts):bank[c-a*r]=bank.get(c-a*r,0)|(1<<i)
  ans.extend(bank.values())
 return ans
@lru_cache(None)
def tau(pts,slopes):
 active=[x for x in lines(pts,slopes) if x.bit_count()>2]
 for k in range(len(pts)+1):
  for subset in combinations(range(len(pts)),k):
   removed=sum(1<<i for i in subset)
   if all((mask&~removed).bit_count()<=2 for mask in active):return k
 raise RuntimeError('empty retained set must work')
def maximum_matching(n,edges):
 adj=[set() for _ in range(n)]
 for a,b in edges:adj[a].add(b);adj[b].add(a)
 @lru_cache(None)
 def visit(mask):
  if not mask:return 0
  a=(mask&-mask).bit_length()-1;rest=mask&~(1<<a)
  return max([visit(rest)]+[1+visit(rest&~(1<<b)) for b in adj[a] if rest>>b&1])
 return visit((1<<n)-1)
def independence(n,edges):
 return max((mask.bit_count() for mask in range(1<<n) if all(not(mask>>a&1 and mask>>b&1) for a,b in edges)),default=0)
def fractional_pair(n,edges):
 # Independent exact rational primal/dual optima: finite half-integer banks
 # are certified optimal by their equality, not by presumed integrality.
 @lru_cache(None)
 def primal(i,caps):
  if i==len(edges):return 0
  a,b=edges[i];best=0
  for amount in range(min(caps[a],caps[b])+1):
   after=list(caps);after[a]-=amount;after[b]-=amount
   best=max(best,amount+primal(i+1,tuple(after)))
  return best
 top=primal(0,(2,)*n)
 cover=min(sum(z) for z in product(range(3),repeat=n) if all(z[a]+z[b]>=2 for a,b in edges))
 need(top==cover,'exact fractional matching primal and vertex-cover dual agree')
 return Q(top,2)
def branch_value(pts,active):
 v=len(active)
 incidence=[tuple(j for j,mask in enumerate(active) if mask>>i&1) for i in range(len(pts))]
 H=[I for I in incidence if len(I)>=3];C=[I for I in incidence if 1<=len(I)<=2]
 @lru_cache(None)
 def keep(i,caps):
  if i==len(C):return 0
  best=keep(i+1,caps);I=C[i]
  if all(caps[j]>0 for j in I):
   nxt=list(caps)
   for j in I:nxt[j]-=1
   best=max(best,1+keep(i+1,tuple(nxt)))
  return best
 best=0
 for mask in range(1<<len(H)):
  caps=[2]*v
  for i,I in enumerate(H):
   if mask>>i&1:
    for j in I:caps[j]-=1
  if min(caps,default=0)>=0:best=max(best,mask.bit_count()+keep(0,tuple(caps)))
 return len(H)+len(C)-best,len(H),incidence
def graph_data(pts,slopes=(1,-1,2)):
 active=[m for m in lines(pts,slopes) if m.bit_count()>2]
 exact=tau(pts,slopes);branch,h,I=branch_value(pts,active)
 need(exact==branch,'all-cell deletion agrees with exceptional-cell branch')
 regime=all(m.bit_count()==3 for m in active) and all(len(j)<=2 for j in I)
 if not regime:return None,h
 edges=[j for j in I if len(j)==2];need(len(edges)==len(set(edges)),'distinct shared cells preserve actual line intersections')
 v=len(active);nu=maximum_matching(v,edges);alpha=independence(v,edges);nf=fractional_pair(v,tuple(edges))
 beta=0
 for mask in range(1<<v):
  union=0
  for j,line in enumerate(active):
   if mask>>j&1:union|=line
  beta=max(beta,union.bit_count()-2*mask.bit_count())
 need(exact==v-nu,'exact graph edge-cover completion')
 need(beta==alpha,'Boolean line selection equals graph independence')
 need(Q(exact)>=v-nf>=beta,'fractional and Boolean comparison paid by duality')
 return (v,tuple(edges),exact,beta,v-nf),h
census={}
for n in range(2,6):
 stats=Counter()
 for rows in boards(n):
  data,h=graph_data(points(rows));stats['boards']+=1;stats['exceptional']+=h>0
  if data is not None:
   stats['regime']+=1
   v,ed,t,beta,lp=data
   stats['fractional_gap']+=t>lp
   unseen=set(range(v));bad_cycle=False
   while unseen:
    seed=next(iter(unseen));component={seed};todo=[seed]
    while todo:
     a=todo.pop()
     for j,k in ed:
      if a in (j,k):
       other=k if a==j else j
       if other not in component:component.add(other);todo.append(other)
    unseen-=component
    if len(component)==5 and all(sum(a==j or b==j for a,b in ed)==2 for j in component):bad_cycle=True
   need(not bad_cycle,'no five-cycle component before six')
 need(stats['boards']=={2:1,3:6,4:90,5:2040}[n],'independent full row-product universe')
 need(stats['regime']=={2:1,3:6,4:73,5:1240}[n],'exact eligible regime counts')
 census[n]=dict(stats)
need(census[5]['exceptional']==16 and census[5]['fractional_gap']==6,'complete hostile incidence census')

named={
 'five_cycle':((0,3),(1,2),(0,4),(1,4),(3,5),(2,5)),
 'concurrent':((2,4),(0,1),(2,3),(3,4),(0,1)),
 'attached_triangle':((0,1),(0,2),(3,4),(1,3),(2,4)),
 'triangle':((0,3),(1,2),(0,4),(3,4),(1,2)),
 'safe_bonus':((0,1),(1,3),(0,5),(2,3),(2,4),(4,5))}
named_values={}
for name,rows in named.items():
 pts=points(rows);data,h=graph_data(pts);named_values[name]={'tau':tau(pts,(1,-1,2)),'exceptional':h}
 if data is not None:named_values[name].update(beta=data[3],lp=str(data[4]),vertices=data[0],edges=data[1])
 # Every singleton direction-block order on the actual named geometric controls.
 for order in permutations((1,-1,2,-2)):
  U=pts;old=[];charge=0
  for slope in order:
   charge+=tau(tuple(U),(slope,))
   old.append(slope)
   bad=0
   for mask in lines(pts,tuple(old)):
    if mask.bit_count()>2:bad|=mask
   U=tuple(p for i,p in enumerate(pts) if not bad>>i&1)
  need(tau(pts,order)>=charge,'all singleton block orders retain original-safe regions')
 # Three-plus-incidence branching remains exact when a fourth direction is added.
 graph_data(pts,(1,-1,2,-2))
need(named_values['five_cycle']['tau']==3 and named_values['five_cycle']['beta']==2 and named_values['five_cycle']['lp']=='5/2','actual C5 exact gap')
need(named_values['concurrent']['tau']==1 and named_values['concurrent']['exceptional']==1,'one concurrent hyperedge is one deletion')
need(named_values['attached_triangle']['tau']==2 and named_values['attached_triangle']['lp']=='2','odd cycle with attachment has no gap')

# Independent abstract subcubic graph bank, solved by literal set-cover DP.
abstract={}
for n in range(6):
 count=0;universe=list(combinations(range(n),2))
 for mask in range(1<<len(universe)):
  edges=[e for j,e in enumerate(universe) if mask>>j&1]
  degree=[sum(a==v or b==v for a,b in edges) for v in range(n)]
  if max(degree,default=0)>3:continue
  count+=1;nu=maximum_matching(n,edges);alpha=independence(n,edges);nf=fractional_pair(n,tuple(edges))
  cells=[(1<<a)|(1<<b) for a,b in edges]+[1<<v for v in range(n) for _ in range(3-degree[v])]
  dp={0:0}
  for cell in cells:
   nxt=dict(dp)
   for cover,cost in dp.items():nxt[cover|cell]=min(nxt.get(cover|cell,10**6),cost+1)
   dp=nxt
  need(dp[(1<<n)-1]==n-nu,'abstract literal minimum cell cover')
  need(alpha<=n-nf<=n-nu,'complete abstract fractionality inequalities')
 need(count==[1,1,2,8,64,768][n],'full labelled subcubic graph count')
 abstract[n]=count

cert={'status':'FINITE-EXACT independent audit controls; analytic proofs in companion referee',
 'universe':'all2137 two-regular boards n2..5; five named controls and their fourth-direction/block tests; all844 labelled subcubic graphs <=5vertices',
 'complete_census':census,'named':named_values,'abstract_graphs':abstract,'always_active_gates':gates}
(OUT/(HERE.stem+'_certificate.json')).write_text(json.dumps(cert,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
print('INDEPENDENT NO3: all2137 boards,1320 eligible; exact source deletion vs exceptional-cell branching')
print('GRAPH REGIME: tau=v-nu, LP=v-nu_f, Boolean=alpha;844 complete subcubic graph controls')
print('GEOMETRY: C5 gives2<5/2<3; concurrent triple needs1; attached odd cycle has no gap')
print('BLOCKS: all24 singleton direction orders on five actual boards preserve original-safe charges')
print('Always-active exact gates:',gates)
