#!/usr/bin/env python3
"""Independent small native forward-word metagraph and suspension/readout probe."""
import itertools,json,argparse
from pathlib import Path
CHECKS=0

def need(x,msg):
 global CHECKS
 CHECKS+=1
 if not x: raise RuntimeError(msg)
def matrix(n,w):
 B=[[0]*n for _ in range(n)]; k=0
 for i in range(n):
  for j in range(i+1,n):
   a=1 if j==i+1 else (w>>k)&1
   if j>i+1:k+=1
   B[i][j]=a;B[j][i]=1-a
 return B
def canonical(B):
 groups={}
 for i,r in enumerate(B):groups.setdefault(sum(r),[]).append(i)
 best=None
 for tup in itertools.product(*(itertools.permutations(groups[d]) for d in sorted(groups))):
  order=sum((tuple(t) for t in tup),());value=0
  for i in range(len(B)):
   for j in range(i+1,len(B)):value=(value<<1)|B[order[i]][order[j]]
  if best is None or value<best:best=value
 return best
def merged(B):return min(canonical(B),canonical(list(map(list,zip(*B)))))
def build(n):
 bits=(n-1)*(n-2)//2;full=(1<<bits)-1
 owners=[merged(matrix(n,w)) for w in range(full+1)]
 keys=sorted(set(owners));index={v:i for i,v in enumerate(keys)};owner=[index[v] for v in owners]
 L=[[0]*len(keys) for _ in keys];edges={}
 for w in range(full+1):
  if w>full^w:continue
  u,v=owner[w],owner[full^w];a,b=sorted((u,v));edges[a,b]=edges.get((a,b),0)+1
  if u!=v:L[u][u]+=1;L[v][v]+=1;L[u][v]-=1;L[v][u]-=1
 need(all(sum(r)==0 for r in L),'literal full metagraph row sums')
 return {'order':n,'keys':keys,'owner':owner,'L':L,'edges':[[u,v,w] for (u,v),w in sorted(edges.items())]}
def suspend(B):
 n=len(B);a=n;b=n+1;S=[[0]*(n+2) for _ in range(n+2)]
 for i in range(n):
  for j in range(n):S[i][j]=B[i][j]
  S[a][i]=S[i][b]=1
 S[a][b]=1;return S

def partner(B,p):
 edges={frozenset((a,b)) for a,b in zip(p,p[1:])}
 return [[0 if i==j else B[i][j] if frozenset((i,j)) in edges else B[j][i] for j in range(len(B))] for i in range(len(B))]
def mv(A,x,m):return tuple(sum(a*b for a,b in zip(r,x))%m for r in A)
def solve(A,b,m):
 """Exact local Smith elimination over Z/2^k, returning one actual solution."""
 original=[list(r) for r in A];original_b=tuple(b);A=[list(r) for r in A];b=list(b)
 rows=len(A);cols=len(A[0]);V=[[int(i==j) for j in range(cols)] for i in range(cols)];rank=0
 for p in range(min(rows,cols)):
  candidates=[((A[i][j]&-A[i][j]),i,j) for i in range(p,rows) for j in range(p,cols) if A[i][j]%m]
  if not candidates:break
  g,i,j=min(candidates)
  A[p],A[i]=A[i],A[p];b[p],b[i]=b[i],b[p]
  for r in A:r[p],r[j]=r[j],r[p]
  for r in V:r[p],r[j]=r[j],r[p]
  unit=A[p][p]//g;inv=pow(unit,-1,m)
  A[p]=[v*inv%m for v in A[p]];b[p]=b[p]*inv%m
  for i in range(p+1,rows):
   q=A[i][p]//g
   A[i]=[(a-q*c)%m for a,c in zip(A[i],A[p])];b[i]=(b[i]-q*b[p])%m
  for j in range(p+1,cols):
   q=A[p][j]//g
   for r in A:r[j]=(r[j]-q*r[p])%m
   for r in V:r[j]=(r[j]-q*r[p])%m
  rank+=1
 z=[0]*cols
 for i in range(rank):
  if b[i]%A[i][i]:return None
  z[i]=b[i]//A[i][i]
 if any(b[rank:]):return None
 x=mv(V,z,m);need(mv(original,x,m)==tuple(a%m for a in original_b),'solved vector rechecked literally')
 return x

def main():
 old=build(4);new=build(6);oldroot=0;idx={v:i for i,v in enumerate(new['keys'])};records=[];sourceimages={}
 for w,u in enumerate(old['owner']):
  B=matrix(4,w);S=suspend(B);p=(4,0,1,2,3,5);K=partner(S,p);v=old['owner'][7^w]
  ns=idx[merged(S)];nt=idx[merged(K)];sourceimages.setdefault(u,set()).add(ns)
  need([r[:4] for r in K[:4]]==partner(B,(0,1,2,3)),'literal inner deletion')
  records.extend([(ns,u),(nt,v)])
 newroot=min(sourceimages[oldroot]);coords=[i for i in range(len(new['keys'])) if i!=newroot]
 N=[[new['L'][i][j] for j in coords] for i in coords]
 A=[r[1:] for r in old['L'][1:]];readouts=sorted(set(records));constraints=[]
 for node,owner in readouts:constraints.append([int(j==node) for j in coords])
 M=N+constraints;tests=[]
 for m in (2,4,8):
  H=[x for x in itertools.product(range(m),repeat=len(A)) if not any(mv(A,x,m))];passed=0;failed=[];witnesses=[]
  for h in H:
   x=(0,)+h;b=[0]*len(N)+[x[owner] for node,owner in readouts];y=solve(M,b,m)
   if y is not None:passed+=1;witnesses.append({'old':h,'new':y})
   else:failed.append(h)
  tests.append({'modulus':m,'old_harmonics':len(H),'extended':passed,'failed_old_vectors':failed,'lifts':witnesses})
 # Build one exact integer readout lift, not merely independently selected residues.
 need(A==[[2,-2],[-2,3]],'old reduced matrix with Smith invariants 1,2')
 h=(1,0);b=[0]*len(N)+[int(owner==1) for node,owner in readouts]
 y=solve(M,b,2);need(y is not None,'old order-two generator extends')
 Z=[[a,0] for a in y]
 J=[[0,0] for _ in N]+[[int(owner==1),int(owner==2)] for node,owner in readouts]
 defect=[[sum(M[i][k]*Z[k][j] for k in range(len(coords)))-J[i][j] for j in range(2)] for i in range(len(M))]
 Q=[]
 for a,b in defect:
  need(a%2==0,'integer row-image criterion')
  Q.append([(3*a+2*b)//2,a+b])
 need([[sum(Q[i][k]*A[k][j] for k in range(2)) for j in range(2)] for i in range(len(Q))]==defect,'exact integer M Z - J = Q A')
 for m in range(2,13):
  for h in itertools.product(range(m),repeat=2):
   if not any(mv(A,h,m)):
    need(mv(M,mv(Z,h,m),m)==mv(J,h,m),'fixed integer lift at every small modulus')
 certificate={'identity':'M Z - J = Q A over integers','Z':Z,'Q':Q,'lift_support_new_nodes':[coords[i] for i,r in enumerate(Z) if r[0]],'all_integer_moduli':True}
 collisions={str(node):sorted({owner for n,owner in records if n==node}) for node in sorted({a for a,b in records}) if len({owner for n,owner in records if n==node})>1}
 return {'status':'FINITE-EXACT complete order4->6 forward-word suspension readout probe; not order8->10 transfer','old':old,'new':new,'old_root':oldroot,'new_root':newroot,'new_reduced_coordinates':coords,'readouts_new_old':readouts,'record_collisions':collisions,'tests':tests,'integer_certificate':certificate,'checks':CHECKS}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path);a=p.parse_args();text=json.dumps(main(),sort_keys=True,indent=2)+'\n'
 if a.out:a.out.write_text(text)
 else:print(text,end='')
