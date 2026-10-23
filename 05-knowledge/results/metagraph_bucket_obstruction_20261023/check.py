"""Exact incidence-signature and port-capacity audit, no solver dependency.
Run python -B .../check.py; -O must return identical mathematical receipt.
"""
from pathlib import Path
from collections import Counter,defaultdict
from itertools import combinations,permutations
import json,hashlib,argparse,csv,math
from itertools import product
ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
checks=0
def need(b,msg):
 global checks
 checks+=1
 if not b:raise AssertionError(msg)
def tiles(n):return [(a,b) for b in range(1,n-1) for a in range(b+2,n+1)]
def matrix(n,m):
 A=[[0]*n for _ in range(n)]
 for i in range(1,n):A[i][i-1]=1
 for k,(a,b) in enumerate(tiles(n)):
  u,v=(b-1,a-1) if m>>k&1 else (a-1,b-1)
  A[u][v]=1
 return A
def code(A,p):return sum(A[p[i]][p[j]]<<k for k,(i,j) in enumerate(combinations(range(len(A)),2)))
def signatures(A):
 n=len(A);s=list(map(sum,A));c=[0]*n;edge=Counter()
 for i,j,k in combinations(range(n),3):
  if A[i][j]==A[j][k]==A[k][i]:
   for v in (i,j,k):c[v]+=1
   for e in combinations((i,j,k),2):edge[e]+=1
 old=(tuple(sorted(s)),tuple(sorted(c)),tuple(sorted(tuple(sorted(s[v] for v in range(n) if A[u][v])) for u in range(n))))
 joint=tuple(sorted((s[u],c[u],tuple(sorted((s[v],c[v]) for v in range(n) if A[u][v]))) for u in range(n)))
 moment=sum(t*(t-1)//2 for t in edge.values())
 v4=[0]*n
 for B in combinations(range(n),4):
  if sorted(sum(A[u][v] for v in B) for u in B)==[1,1,2,2]:
   for v in B:v4[v]+=1
 tags=[(s[v],c[v],v4[v]) for v in range(n)]
 rho1=[(tags[u],tuple(sorted(tags[v] for v in range(n) if A[u][v]))) for u in range(n)]
 joint4=tuple(sorted(rho1))
 joint42=tuple(sorted((rho1[u],tuple(sorted(rho1[v] for v in range(n) if A[u][v]))) for u in range(n)))
 return old,joint,(joint,moment),joint4,joint42
def reflect(n,m):
 ts=tiles(n);ix={e:i for i,e in enumerate(ts)}
 return sum((m>>i&1)<<ix[(n+1-b,n+1-a)] for i,(a,b) in enumerate(ts))
def peel(rowcap,colcap,allow):
 import heapq
 r=dict(rowcap);c=dict(colcap);rn=defaultdict(set);cn=defaultdict(set);fixed={};trace=[];queue=[]
 for a,b in allow:
  if r[a] and c[b]:rn[a].add(b);cn[b].add(a)
 for a in r:
  if len(rn[a])==1:heapq.heappush(queue,(0,a))
 for b in c:
  if len(cn[b])==1:heapq.heappush(queue,(1,b))
 def remove(a,b):
  rn[a].discard(b);cn[b].discard(a)
  if r[a] and len(rn[a])==1:heapq.heappush(queue,(0,a))
  if c[b] and len(cn[b])==1:heapq.heappush(queue,(1,b))
 while queue:
  side,v=heapq.heappop(queue)
  if side==0:
   if not r[v] or len(rn[v])!=1:continue
   a=v;b=next(iter(rn[a]));k=r[a];why='row'
  else:
   if not c[v] or len(cn[v])!=1:continue
   b=v;a=next(iter(cn[b]));k=c[b];why='column'
  need(k<=r[a] and k<=c[b],'capacity conflict');fixed[a,b]=fixed.get((a,b),0)+k;r[a]-=k;c[b]-=k;trace.append([a,b,k,why])
  if not r[a]:
   for bb in list(rn[a]):remove(a,bb)
  if not c[b]:
   for aa in list(cn[b]):remove(aa,b)
 need(all(not k or rn[a] for a,k in r.items()),'no residual source neighbor')
 need(all(not k or cn[b] for b,k in c.items()),'no residual target neighbor')
 return fixed,trace,{a:k for a,k in r.items() if k},{b:k for b,k in c.items() if k},sorted((a,b) for a in rn for b in rn[a])
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true',help='intentionally regenerate the receipt');args=parser.parse_args()
 out={'scope':'FINITE-EXACT, complete labelled tournaments and fixed-path tilings n=4..7','levels':{},'sources':{}}
 for n in (4,5,6,7):
  path=ROOT/f'04-computation/experiments/electric_metagraph_20261021/electric_n{n}.json';raw=path.read_bytes();data=json.loads(raw);out['sources'][str(path.relative_to(ROOT))]=hashlib.sha256(raw).hexdigest()
  classes=[]
  for v in data['nodes']:
   classes.append((v['id'],v['rep_mask']))
   if not v['SC']:classes.append((v['id'],v['conv_mask']))
  mats=[matrix(n,m) for u,m in classes];sg=[signatures(A) for A in mats];lookup={};auts=[]
  for ci,A in enumerate(mats):
   rep=code(A,range(n));aut=0
   for p in permutations(range(n)):
    x=code(A,p);aut+=x==rep
    if x in lookup:need(lookup[x]==ci,'two representative classes isomorphic')
    else:lookup[x]=ci
   auts.append(aut);need(aut==data['nodes'][classes[ci][0]]['Aut'],'automorphism mismatch')
  need(len(lookup)==1<<(n*(n-1)//2),'labelled universe not covered')
  classesigs=[]
  for level in range(4):
   by=defaultdict(list)
   for ci,s in enumerate(sg):by[s[level]].append(ci)
   classesigs.append(by)
  native=[];F=(1<<len(tiles(n)))-1;classfib=Counter();bluefib=Counter();openports=Counter();merged=[Counter(),Counter()]
  for m in range(F+1):
   A=matrix(n,m);ci=lookup[code(A,range(n))];cj=lookup[code(matrix(n,m^F),range(n))];s=list(map(sum,A));blue=m==reflect(n,m);closed=bool(A[0][n-1]);classfib[ci]+=1;bluefib[ci]+=blue
   if not closed:openports[ci,s[-1],s[0],blue]+=1
   if closed:native.append((ci,s[-1],s[0],blue,cj,m));merged[blue][classes[ci][0],classes[cj][0]]+=1
  for ci,(u,m) in enumerate(classes):
   need(classfib[ci]*auts[ci]==data['nodes'][u]['H'],'Hamilton fibre check')
   need(classfib[ci]%2==1,'ordinary fibre parity');need(bluefib[ci]==(data['nodes'][u]['s'] if data['nodes'][u]['SC'] else 0),'blue fibre check')
  for color,key in ((False,'black_closed_to_open'),(True,'blue_closed_to_open')):
   need(all(merged[color][i,j]==v for i,row in enumerate(data[key]) for j,v in enumerate(row)),'native metagraph mismatch')
  for A in mats:
   cyc=[B for B in combinations(range(n),3) if A[B[0]][B[1]]==A[B[1]][B[2]]==A[B[2]][B[0]]]
   cv=[sum(v in B for B in cyc) for v in range(n)];M=signatures(A)[2][1]
   strong4=sum(sorted(sum(A[u][v] for v in B) for u in B)==[1,1,2,2] for B in combinations(range(n),4))
   disjoint=sum(not(set(B)&set(C)) for B,C in combinations(cyc,2))
   need(M==strong4,'overlap is not strong-four count')
   need(disjoint==math.comb(len(cyc),2)-sum(math.comb(x,2) for x in cv)+M,'overlap/disjoint identity')
  ns={}
  for level,name in enumerate(('L3','neighbor_incidence','neighbor_incidence_plus_overlap','neighbor_incidence_four')):
   rows={};rowkeys=[];rowcap=Counter();actual=Counter();witness={}
   for ci,a,b,blue,cj,m in native:
    key=(ci,a,b,blue,sg[cj][level])
    if key not in rows:rows[key]=len(rows);rowkeys.append(key)
    ri=rows[key];rowcap[ri]+=1;actual[ri,(cj,n-a,n-2-b,blue)]+=1;witness.setdefault((ri,cj),m)
   cols=sorted(k for k,v in openports.items() if v);colidx={k:i for i,k in enumerate(cols)};colcap={colidx[k]:openports[k] for k in cols};allow=[];act=Counter()
   for ri,(ci,a,b,blue,sig) in enumerate(rowkeys):
    for cj in classesigs[level][sig]:
     col=(cj,n-a,n-2-b,blue)
     if col in colidx:allow.append((ri,colidx[col]))
   for (ri,col),v in actual.items():act[ri,colidx[col]]+=v
   allowset=set(allow);need(all(e in allowset for e in act),'actual port unsupported')
   fixed,trace,rr,cc,edges=peel(rowcap,colcap,allow)
   need(all(act[e]==v for e,v in fixed.items()),'forced entry differs from actual')
   collisions=[[classes[i][1] for i in vals] for vals in classesigs[level].values() if len(vals)>1]
   ns[name]={'distinct_class_signatures':len(classesigs[level]),'classes':len(classes),'collision_masks':collisions,'port_rows':len(rows),'candidate_entries':len(allow),'forced_entries':len(fixed),'unresolved_mass':sum(rr.values()),'unresolved_rows':len(rr),'unresolved_columns':len(cc),'peeling_steps':len(trace),'peeling_sha256':hashlib.sha256(json.dumps(trace,separators=(',',':')).encode()).hexdigest()}
   if n==6 and level==0:
    # Explicit reflected 2x2 switch: preserve source row and target column margins.
    rowfind={(k[0],k[1],k[2],k[3],k[4]):i for i,k in enumerate(rowkeys)}
    edits=[(31,4,2,22,-1),(31,4,2,33,1),(29,4,2,33,-1),(29,4,2,22,1),(30,3,1,21,-1),(30,3,1,32,1),(28,3,1,32,-1),(28,3,1,21,1)]
    delta=Counter();changed=[]
    for ci,a,b,cj,v in edits:
     ri=rowfind[ci,a,b,False,sg[cj][level]];co=colidx[cj,n-a,n-2-b,False];delta[ri,co]+=v;changed.append({'source_class':ci,'source_mask':classes[ci][1],'source_node':classes[ci][0],'end_scores':[a,b],'target_class':cj,'target_mask':classes[cj][1],'target_node':classes[cj][0],'delta':v,'existing_example_tiling':witness.get((ri,cj))})
    fake=act.copy();fake.update(delta);need(all(v>=0 for v in fake.values()),'negative fake entry');need(all(e in set(allow) for e in fake),'illegal fake support')
    need(all(sum(v for (a,b),v in delta.items() if a==r)==0 for r in rowcap),'fake changes source bucket');need(all(sum(v for (a,b),v in delta.items() if b==c)==0 for c in colcap),'fake changes target bucket')
    dm=Counter()
    for (ri,co),v in delta.items():dm[classes[rowkeys[ri][0]][0],classes[cols[co][0]][0]]+=v
    need(any(dm.values()),'fake quotient unchanged');need(all(v%2==0 for v in dm.values()),'fake breaks black twins')
    conv={i:next(j for j,(u,m) in enumerate(classes) if u==classes[i][0] and (i!=j or data['nodes'][u]['SC'])) for i in range(len(classes))}
    ag=Counter()
    for (ri,co),v in fake.items():
     if not rowkeys[ri][3]:ag[rowkeys[ri][0],cols[co][0]]+=v
    need(all(v==ag[conv[i],conv[j]] for (i,j),v in ag.items()),'fake breaks class converse reflection')
    out['counterfeit_n6']={'edits':changed,'merged_delta':[[list(k),v] for k,v in sorted(dm.items()) if v]}
  out['levels'][str(n)]=ns
  print(n,[(k,v['distinct_class_signatures'],v['unresolved_mass']) for k,v in ns.items()],flush=True)
 pair=[]
 for mask in (220,108):
  A=matrix(6,mask);s=signatures(A);pair.append({'mask':mask,'outgoing_rows':[sum(b<<i for i,b in enumerate(r)) for r in A],'overlap_moment':s[2][1]})
 need(signatures(matrix(6,220))[1]==signatures(matrix(6,108))[1],'hostile signatures differ');need(pair[0]['overlap_moment']==9 and pair[1]['overlap_moment']==12,'moment witness mismatch')
 out['one_moment_pair']=pair
 path=ROOT/'04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n8.tsv';out['sources'][str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
 groups=defaultdict(list);groups2=defaultdict(list);canons=set();mass=0;canon_checks=0
 for row in csv.DictReader(path.open(),delimiter='\t'):
  n=8;cc=int(row['canon']);A=[[0]*n for _ in range(n)]
  for k,(i,j) in enumerate(combinations(range(n),2)):A[i][j]=(cc>>(27-k))&1;A[j][i]=1-A[i][j]
  s=list(map(sum,A));blocks=[[i for i in range(n) if s[i]==v] for v in sorted(set(s))];best=None;aut=0
  for blocks_p in product(*(permutations(B) for B in blocks)):
   pp=sum(blocks_p,());value=0
   for i,j in combinations(range(n),2):value=2*value+A[pp[i]][pp[j]]
   canon_checks+=1
   if best is None or value<best:best=value;aut=1
   elif value==best:aut+=1
  need(best==cc,'n8 representative is not score canonical');need(aut==int(row['Aut']),'n8 aut mismatch');need(cc not in canons,'n8 duplicate representative');canons.add(cc);mass+=math.factorial(n)//aut
  sigs8=signatures(A);groups[sigs8[3]].append((row,A));groups2[sigs8[4]].append(cc)
 need(len(canons)==6880 and mass==1<<28,'n8 exhaustive orbit accounting')
 hostile=[]
 for g in groups.values():
  if len(g)<2:continue
  records=[]
  for row,A in g:
   dp=[[0]*8 for _ in range(256)]
   for u in range(8):dp[1<<u][u]=1
   for mask in range(1,256):
    for u in range(8):
     if not dp[mask][u]:continue
     for v in range(8):
      if not(mask>>v&1) and A[u][v]:dp[mask|1<<v][v]+=dp[mask][u]
   H=sum(dp[-1]);need(H==int(row['H']),'n8 hostile H mismatch');records.append({'canon':int(row['canon']),'H':H,'rows':[sum(b<<i for i,b in enumerate(r)) for r in A]})
  hostile.append(records)
 need(all(len({int(row['H'])//int(row['Aut']) for row,A in g})==len(g) for g in groups.values()),'n8 J4+h collision')
 need(len(groups2)==6880,'n8 second-round collision')
 out['n8_boundary']={'second_round_signatures':len(groups2),'with_h_signatures':6880,'scope':'all6880classes independently score-canonical, orbit sizes cover2^28 orientations','canonical_permutations_checked':canon_checks,'signatures':len(groups),'collision_groups':hostile}
 print('n8 J4',len(groups),'/',len(canons),'canonperm',canon_checks,flush=True)
 out['checks']=checks
 dest=HERE/('receipt.optimized.json' if not __debug__ else 'receipt.json')
 if args.write:dest.write_text(json.dumps(out,indent=2)+'\n')
 elif json.loads((HERE/'receipt.json').read_text())!=out:raise AssertionError('reproduction differs from frozen receipt')
 print('checks',checks,'receipt',('wrote ' if args.write else 'matched ')+dest.name)
if __name__=='__main__':main()
