#!/usr/bin/env python3
"""Exact consumer for actual endpoint repair and shared C3 chirality.
No removable assertions, no external package, no probabilistic signature hash.
"""
import argparse,collections,hashlib,itertools,json
from pathlib import Path
from analyze_collisions import matrix,refine,sign,path
HERE=Path(__file__).resolve().parent
checks=collections.Counter()
def need(x,msg):
 checks[msg]+=1
 if not x:raise ValueError(msg)
def literal(m,n=8):
 A=[[0]*n for _ in range(n)]
 for i in range(n-1):A[i+1][i]=1
 for k,(i,j) in enumerate(( (i,j) for i in range(n-2) for j in range(i+2,n))):A[i][j]=m>>k&1;A[j][i]=1-A[i][j]
 return A
def converse(A):return [[0 if i==j else 1-A[i][j] for j in range(len(A))] for i in range(len(A))]
def first_difference(A,B,marks):
 n=len(A);C=[[(tuple(int(v==u) for u in marks)) for v in range(n)] for _ in range(2)]
 for r in range(n+1):
  if sorted(C[0])!=sorted(C[1]):return r
  ids={x:i for i,x in enumerate(sorted(set(C[0]+C[1])))};C=[[ids[x] for x in cc] for cc in C]
  C=[[(cc[v],tuple(sorted(cc[w] for w in range(n) if M[v][w]))) for v in range(n)] for cc,M in zip(C,(A,B))]
 return None
def H(A):
 n=len(A);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for m in range(1<<n):
  for v in range(n):
   z=dp[m][v]
   if z:
    for w in range(n):
     if not(m>>w&1) and A[v][w]:dp[m|1<<w][w]+=z
 return sum(dp[-1])
def main():
 p=argparse.ArgumentParser();p.add_argument('--counterfeits',type=Path,default=HERE.parent/'session38-binary-counterfeits/receipt.json');p.add_argument('--write',action='store_true');a=p.parse_args()
 old=json.loads(a.counterfeits.read_text());pairs=sorted({(z['native_open_target_mask'],z['counterfeit_open_target_mask']) for row in old['results'] for z in row['literal_matching']});need(len(pairs)==244,'complete distinct comparison set')
 dist={};bad={}
 for marks in [(),(7,),(0,),(7,0)]:
  D=collections.Counter();bad[marks]=set()
  for x,y in pairs:
   r=first_difference(literal(x),literal(y),marks);D[str(r)]+=1
   if r is None:bad[marks].add((x,y))
  dist[str(marks)]=dict(sorted(D.items()))
 switches={tuple((z['native_open_target_mask'],z['counterfeit_open_target_mask']) for z in row['literal_matching']) for row in old['results']}
 survivors={str(marks):sum(all(pair in bad[marks] for pair in sw) for sw in switches) for marks in bad}
 need(len(switches)==61 and survivors=={'()':2,'(7,)':0,'(0,)':0,'(7, 0)':0},'actual whole switch survival under endpoint marking')
 need(not bad[(7,0)],'joint endpoint repair');need(not(bad[(7,)]&bad[(0,)]),'converse closed single endpoint repair');need(len(bad[()])==8,'unmarked hostile retained')
 need(dist['(7, 0)']=={'1':116,'2':116,'3':12},'exact joint endpoint round distribution')
 # Rebuild every surviving full-census collision from exact canonical matrices.
 groups=[]
 for line in (HERE/'endpoint_scan.txt').read_text().splitlines():
  if not line.startswith('collision '):continue
  group=[];hist=None;sigclass={}
  for item in line.split()[1:]:
   code,s,t=map(int,item.split(':'));A=matrix(code);h,cells=refine(A,s,t);eps=sign(A,cells);pp=path(A,s,t)
   need(pp is not None and len(set(pp))==8 and pp[0]==s and pp[-1]==t,'actual endpoint Hamiltonian witness')
   need(all(A[u][v] for u,v in zip(pp,pp[1:])),'literal witness arcs')
   need(sorted(map(len,cells))==[1,1,3,3],'every residual collision has two triples')
   if hist is None:hist=h
   need(h==hist,'unaugmented marked refinement collision')
   need(eps in (-1,1),'nonuniform shared selector exists')
   sigclass.setdefault(eps,set()).add(code)
   hs,cs=refine(converse(A),t,s);need(sign(converse(A),cs)==eps,'chirality converse even with endpoint swap')
   group.append({'code':code,'start':s,'end':t,'cells':cells,'sign':eps,'J_mod8':(3*eps)%8,'path':pp})
  need(set(sigclass)=={-1,1} and all(len(v)==1 for v in sigclass.values()),'one bit repairs entire collision group')
  groups.append(group)
 need(len(groups)==18,'all full census collisions consumed')
 # Unique regular-five tournament: all24 labelled regular matrices form one orbit.
 R=[]
 for code in range(1024):
  A=matrix(code,5)
  if all(sum(row)==2 for row in A):R.append(code)
 need(len(R)==24,'all labelled regular five tournaments')
 A=matrix(R[0],5);orbit=set()
 for perm in itertools.permutations(range(5)):
  x=0
  for i,j in itertools.combinations(range(5),2):x=2*x+A[perm[i]][perm[j]]
  orbit.add(x)
 need(orbit==set(R),'one regular five orbit')
 # Every orientation/polarity/permutation in the two-triple matching family.
 fourpoint=0
 for flipx,flipy in itertools.product((0,1),repeat=2):
  for perm in itertools.permutations(range(3)):
   for q in (1,2):
    B=[[0]*6 for _ in range(6)]
    for off,flip in [(0,flipx),(3,flipy)]:
     for v in range(3):B[off+v][off+(v+(2 if flip else 1))%3]=1
    for u,v in itertools.product(range(3),repeat=2):B[u][3+v]=int((perm[u]==v)==(q==1));B[3+v][u]=1-B[u][3+v]
    pi={u:3+perm[u] for u in range(3)}
    J=sum(2*B[pi[u]][pi[w]]-1 for u,w in itertools.permutations(range(3),2) if B[u][w]);need(J in(-3,3),'four vertex contraction range')
    Jc=sum(2*(1-B[pi[u]][pi[w]])-1 for u,w in itertools.permutations(range(3),2) if not B[u][w]);need(Jc==J,'shared rectangle converse even')
    need(J%8 in(3,5) and J%3==0 and J%2==1,'mod8 retains sign where mod2 mod3 lose it');fourpoint+=1
 # Five-cell boundary: all permutation matchings modulo C5 rotations.
 five_orbits={}
 for pp in itertools.permutations(range(5)):
  orb={tuple((pp[(x-a)%5]+b)%5 for x in range(5)) for a,b in itertools.product(range(5),repeat=2)}
  key=min(orb);par=sum(pp[i]>pp[j] for i in range(5) for j in range(i+1,5))%2
  five_orbits[key]=(len(orb),par)
 need(len(five_orbits)==8,'five cell matching has eight rotational double orbits')
 need(collections.Counter(z[0] for z in five_orbits.values())=={5:4,25:4},'five cell orbit sizes')
 need(collections.Counter(z[1] for z in five_orbits.values())=={0:4,1:4},'permutation parity leaves four five cell ambiguities')
 for aa in range(1,12,2):
  for bb in range(1,12,2):
   import math
   if math.gcd(aa,bb)==1:
    need([(q,p) for q in range(bb+1) for p in range(aa+1) if aa*q==bb*p]==[(0,0),(bb,aa)],'coprime quotient block must be uniform')
 # Literal target refinement transport on complete fixed-path universes through6.
 transport=0
 for n in range(2,7):
  L=(n-1)*(n-2)//2
  for m in range(1<<L):
   A=literal(m,n);B=literal(m^((1<<L)-1),n);c=[0 if v==n-1 else 1 if v==0 else 2 for v in range(n)]
   for roundno in range(n):
    keys=[]
    for v in range(n):
     z=[]
     for col in sorted(set(c)):
      C=[w for w in range(n) if c[w]==col]
      direct=sum(B[v][w] for w in C)
      transported=len(C)-int(v in C)-sum(A[v][w] for w in C)+int(v>0 and v-1 in C)-int(v<n-1 and v+1 in C)
      need(direct==transported,'actual path ownership target cell transport');z.append(direct)
     keys.append((c[v],tuple(z)))
    ids={x:i for i,x in enumerate(sorted(set(keys)))};d=[ids[x] for x in keys]
    if len(set(c))==len(set(d)):break
    c=d
   transport+=1
 # Sharp next-order boundary, actual common path and distinct H.
 hostile=[]
 for steps in [(1,2,3),(1,2,4)]:
  R=[[int(i!=j and (j-i)%7 in steps) for j in range(7)] for i in range(7)]
  B=[[0]*9 for _ in range(9)]
  for i in range(1,9):B[0][i]=1
  for i in range(1,8):
   B[i][8]=1
   for j in range(1,8):B[i][j]=R[i-1][j-1]
  need(all(B[i][i+1] for i in range(8)),'common explicit order nine Hamiltonian path')
  h,cells=refine(B,0,8);need(sorted(map(len,cells))==[1,1,7],'seven vertex invisible module')
  need(sign(B,cells) is None,'triple sign sidecar empty at boundary')
  need(all(R[i][j]==R[(-j)%7][(-i)%7] for i in range(7) for j in range(7)),'regular circulant self converse')
  hostile.append({'steps':steps,'H':H(R),'history':h})
 need(hostile[0]['history']==hostile[1]['history'],'same marked quotient at order nine');need(hostile[0]['H']!=hostile[1]['H'],'nonisomorphism by Hamiltonian count')
 out={'status':'PASS','scope':'Complete old literal switch repair; all18 full endpoint-census collision consumers; all fixed-path targets through6; sharp regular-seven boundary. Full n8 census reproduced separately by scan_endpoints.cpp.','old_target_comparisons':len(pairs),'whole_literal_switch_survivors':survivors,'refinement_rounds':dist,'single_start_survivors':sorted(bad[(7,)]),'single_end_survivors':sorted(bad[(0,)]),'full_census_collision_groups':groups,'regular5_labelled':24,'fourpoint_controls':fourpoint,'five_cell_double_orbits':len(five_orbits),'fixed_path_transport_universes':transport,'order9_hostile':[{'steps':z['steps'],'H':z['H']} for z in hostile],'checks':dict(sorted(checks.items())),'input_sha256':{'counterfeit_receipt':hashlib.sha256(a.counterfeits.read_bytes()).hexdigest(),'endpoint_scan.txt':hashlib.sha256((HERE/'endpoint_scan.txt').read_bytes()).hexdigest()}}
 dest=HERE/'receipt.json'
 if a.write:dest.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 else:need(json.loads(dest.read_text())==json.loads(json.dumps(out)),'frozen receipt exact reproduction')
 print(json.dumps({'status':'PASS','checks':sum(checks.values()),'old_comparisons':244,'full_census_collision_groups':18,'order9_hostile':out['order9_hostile']},indent=2))
if __name__=='__main__':main()
