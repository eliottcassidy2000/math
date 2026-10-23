import json,itertools,math,hashlib
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
checks=0
def ck(a,msg):
 global checks
 checks+=1
 if not a:raise RuntimeError(msg)
def subsets(xs):
 return (tuple(x for i,x in enumerate(xs) if b>>i&1) for b in range(1<<len(xs)))
def connected(vs,es):
 if not vs:return True
 reached={next(iter(vs))}
 while True:
  newer=reached|{v for u,v in es if u in reached}|{u for u,v in es if v in reached}
  if newer==reached:return set(vs)<=reached
  reached=newer
def rank(a,p=None):
 a=[[x%p if p else F(x) for x in r] for r in a];n=0
 for j in range(len(a[0]) if a else 0):
  k=next((k for k in range(n,len(a)) if a[k][j]),None)
  if k is None:continue
  a[n],a[k]=a[k],a[n]; inv=pow(a[n][j],-1,p) if p else 1/a[n][j]
  a[n]=[(x*inv)%p if p else x*inv for x in a[n]]
  for k in range(n+1,len(a)):
   b=a[k][j]
   if b:a[k]=[(x-b*y)%p if p else x-b*y for x,y in zip(a[k],a[n])]
  n+=1
 return n
# Direct assignment polytope's integer points, no augmenting-path or cut code.
profiles=0; demands=0
for bits in itertools.product(range(2),repeat=7):
 owners=[i+1 for i,v in enumerate(bits) if v]
 attainable={(0,0,0)}
 for mask in owners:
  attainable|={tuple(d[k]+(k==j) for k in range(3)) for d in list(attainable) for j in range(3) if mask>>j&1}
 r=[sum(bool(o&s) for o in owners) for s in range(8)]
 for d in itertools.product(range(4),repeat=3):
  served=max(sum(x) for x in attainable if all(a<=b for a,b in zip(x,d)))
  defect=max(sum(d[j] for j in range(3) if s>>j&1)-r[s] for s in range(8))
  ck(sum(d)-served==defect,'assignment vs Hall');demands+=1
 for s in range(8):
  ck(max(sum(d) for d in attainable if all(d[j]==0 for j in range(3) if not(s>>j&1)))==r[s],'rank recoverable')
 profiles+=1
# Exhaust every small integer set-function, not just known realizable profiles.
realizable=0
for vals in itertools.product(range(4),repeat=7):
 r=(0,)+vals;M=r[7];f=[M-r[7^t] for t in range(8)]
 c=[sum((-1)**((a^t).bit_count())*f[t] for t in range(8) if t&a==t) for a in range(8)]
 ck(all(sum(c[a] for a in range(1,8) if a&s)==r[s] for s in range(8)),'inverse identity arbitrary function')
 if min(c)>=0:realizable+=1
# Integer-only observer cannot distinguish real capacities .2 and .8 on one label.
ck([d for d in range(3) if d<=F(1,5)]==[d for d in range(3) if d<=F(4,5)],'real/integer guard')
ck(F(1,2)>F(1,5) and F(1,2)<=F(4,5),'real witness')
# Rational mask realization and inversion.
c=[F(0),F(1,2),F(2,3),F(3,7),F(5,11),F(0),F(1,13),F(2,17)]
r=[sum(c[a] for a in range(8) if a&s) for s in range(8)]
ck(all(sum((-1)**((a^t).bit_count())*(r[7]-r[7^t]) for t in range(8) if t&a==t)==c[a] for a in range(1,8)),'rational inversion')
cert=json.loads((ROOT/'ports/certificate.json').read_text());matching=[]
for item in cert['matching']:
 n=item['ports'];ms=item['matchings']
 matrix=[[int(connected(range(n),p+q)) for q in ms] for p in ms]
 ck(matrix==item['matrix'],'matching matrix')
 # Full rank mod a prime certifies full rational rank, independently of producer elimination.
 rq=rank(matrix,1000003);rf=rank(matrix,2)
 ck(rq==len(ms)==item['rank_Q'] and rf==item['rank_F2'],'matching ranks')
 ck(len(set(map(tuple,matrix)))==len(ms),'distinct matching responses')
 for i,j,ab,r in item['separators']:
  ck(connected(range(n),ms[j]+r) and not connected(range(n),ms[i]+r),'pair separator')
 matching.append([n,len(ms),rq,rf])
# Native G25 independent graph and DFS; unique leaf forces endpoint18.
es=[(u,v) for u in range(1,26) for v in range(u+1,26) if u+v in (4,9,16,25,36,49)]
adj={u:tuple(v if u==a else a for a,v in es if u in (a,v)) for u in range(1,26)}
ck([u for u in adj if len(adj[u])==1]==[18],'forced leaf')
paths=[]
def walk(path,left):
 if not left:paths.append(list(path));return
 for v in adj[path[-1]]:
  if v in left:walk(path+(v,),left-{v})
walk((18,),set(range(1,26))-{18})
ck(sorted(paths)==sorted(cert['G25']['paths']) and len(paths)==10,'all native paths')
blocks={'A':{2,7,9,16,18},'B':{4,5,11,12},'C':{1,3,6,8,10,17}}
core=set(range(1,26))-set.union(*blocks.values());ce=[e for e in es if set(e)<=core]
local={};ranks={}
for name,vs in blocks.items():
 b=cert['G25']['blocks'][name];ports=b['ports'];inner=[e for e in es if set(e)<=vs]
 found=set()
 for chosen in subsets(inner+[tuple(p) for p in ports]):
  deg={v:sum(v in e for e in chosen) for v in vs}
  if any(d not in (1,2) for d in deg.values()) or sum(d==1 for d in deg.values())>2:continue
  ie=[e for e in chosen if set(e)<=vs]; pi=tuple(i for i,p in enumerate(ports) if tuple(p) in chosen)
  # Every internal component must reach a selected boundary; internal forest via edge count per component.
  remaining=set(vs);ok=True
  while remaining:
   seed=next(iter(remaining));comp={seed}
   while True:
    new=comp|{v for e in ie if set(e)&comp for v in e}
    if new==comp:break
    comp=new
   remaining-=comp
   if sum(set(e)<=comp for e in ie)!=len(comp)-1 or not any(ports[i][0] in comp for i in pi):ok=False;break
  if ok:found.add((tuple(sorted(ie)),pi))
 keys=[(tuple(sorted(map(tuple,s['internal_edges']))),tuple(s['selected_ports'])) for s in b['states']]
 ck(set(keys)==found and len(keys)==len(found),'complete local profiles')
 local[name]=keys
sols=[]
for path in paths:
 pe={tuple(sorted(e)) for e in zip(path,path[1:])};ids=[]
 for name,vs in blocks.items():
  b=cert['G25']['blocks'][name]
  key=(tuple(sorted(e for e in pe if set(e)<=vs)),tuple(i for i,p in enumerate(b['ports']) if tuple(sorted(p)) in pe))
  ids.append(local[name].index(key))
 ids.append(list(subsets(ce)).index(tuple(e for e in ce if e in pe)));sols.append(ids)
ck(sorted(sols)==sorted(cert['G25']['solutions']),'solution indices')
for i,name in enumerate(blocks):
 rows=[set() for _ in local[name]]
 for s in sols:rows[s[i]].add(tuple(s[:i]+s[i+1:]))
 cols=sorted(set.union(*rows));a=[[int(c in r) for c in cols] for r in rows]
 classes=len(set(tuple(sorted(r)) for r in rows));rq=rank(a)
 b=cert['G25']['blocks'][name]
 ck(classes==b['class_count_including_dead'] and rq==b['rank_Q'],'context quotient')
 ranks[name]=[len(rows),classes,rq]
# Independent simultaneous quotient and endpoint-monomial contraction.
from collections import Counter
class_lists=[]
for name in blocks:
 b=cert['G25']['blocks'][name]
 ck(b['total_context_count']==16*math.prod(len(local[n]) for n in blocks if n!=name),'all contexts counted')
 class_lists.append([g['states'] for g in b['classes']])
endpoint_poly=Counter();weighted=0
for groups in itertools.product(*class_lists):
 for ci,core_edges in enumerate(subsets(ce)):
  rep=[g[0] for g in groups]+[ci]
  if rep not in sols:continue
  weighted+=math.prod(map(len,groups))
  for ids in itertools.product(*groups):
   edges=list(core_edges)
   for name,j in zip(blocks,ids):
    b=cert['G25']['blocks'][name];state=b['states'][j]
    edges+=list(map(tuple,state['internal_edges']))+[tuple(b['ports'][i]) for i in state['selected_ports']]
   degrees=Counter(v for e in edges for v in e)
   ck(connected(range(1,26),edges) and sorted(degrees.values())==[1,1]+[2]*23,'every quotient lift')
   endpoint_poly[tuple(sorted(v for v,d in degrees.items() if d==1))]+=1
ck(weighted==10,'weighted quotient count')
ck(endpoint_poly==Counter(tuple(sorted((p[0],p[-1]))) for p in paths),'endpoint quotient polynomial')
ck(endpoint_poly==Counter({tuple(x['endpoints']):x['coefficient'] for x in cert['G25']['endpoint_polynomial']}),'published endpoint polynomial')
# P4 and claw: all permutations test HP, exact Schur complement independent sympy.
import sympy as S
responses=[];hps=[]
for e in [[(0,1),(1,2),(2,3)],[(0,1),(1,2),(1,3)]]:
 L=S.zeros(4)
 for u,v in e:L[u,u]+=1;L[v,v]+=1;L[u,v]-=1;L[v,u]-=1
 B=[0,2];I=[1,3];ii=L.extract(I,I)
 responses.append((L.extract(B,B)-L.extract(B,I)*ii.inv()*L.extract(I,B),ii.det()))
 hps.append(sum(all(tuple(sorted(x)) in e for x in zip(p,p[1:])) for p in itertools.permutations(range(4)))//2)
ck(responses[0]==responses[1] and hps==[1,0],'electrical hostile')
out={'status':'PASS','explicit_checks':checks,'clock_profiles':profiles,'clock_demands':demands,'all_small_rank_tables':4**7,'nonnegative_inverse_tables':realizable,'matching_ports_states_Qrank_F2rank':matching,'G25_original_classes_Qrank':ranks,'paths':paths,'electrical_HP_counts':hps,'scope':'Independent exact audit. No producer function imports or Git writes.'}
Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
