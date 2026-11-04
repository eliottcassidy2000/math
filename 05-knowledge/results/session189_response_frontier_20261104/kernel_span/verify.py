#!/usr/bin/env python3
"""Exact finite controls for the proved all-order kernel-span and matching lemmas."""
from pathlib import Path
import sys, json, itertools, math, hashlib
HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'session189-odd-transitive'
sys.path.insert(0,str(BASE))
from verify import exterior, plus, mul, formula, ThreeCover

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def graph(n,mask):
 rows=[0]*n
 for b,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if mask>>b&1:rows[i]|=1<<j
  else:rows[j]|=1<<i
 return rows

def maxmatch(rows,into,exclude=-1):
 O=[i for i in range(len(rows)) if i!=exclude and not(into>>i&1)]
 def visit(p,used):
  if p==len(O):return 0
  ans=visit(p+1,used);opts=rows[O[p]]&into&~used
  if exclude>=0:opts&=~(1<<exclude)
  while opts:
   b=opts&-opts;opts-=b;ans=max(ans,1+visit(p+1,used|b))
  return ans
 return visit(0,0)

def polynomial_controls():
 receipts=[]
 # Exact row reduction is over integers via fraction-free elimination.
 from fractions import Fraction
 for m in range(3,32,2):
  k=(m-1)//2;rs=list(range(2,k+1))+[m-1]
  polys=[formula(m,r) for r in rs]
  for p in polys:
   reflected=[sum(p[j]*(-1)**j*math.comb(j,h) for j in range(h,len(p))) for h in range(len(p))]
   need(reflected==p,f'reflection m={m}')
  rows=[[Fraction(x) for x in p] for p in polys];rank=0
  for c in range(m-2):
   pivot=next((i for i in range(rank,len(rows)) if rows[i][c]),None)
   if pivot is None:continue
   rows[rank],rows[pivot]=rows[pivot],rows[rank];a=rows[rank][c]
   rows[rank]=[v/a for v in rows[rank]]
   for i in range(len(rows)):
    if i!=rank and rows[i][c]:
     a=rows[i][c];rows[i]=[v-a*w for v,w in zip(rows[i],rows[rank])]
   rank+=1
  need(rank==k,f'full rank m={m}')
  receipts.append({'m':m,'selected_r':rs,'rank':rank,'invariant_dimension':k})
 return receipts

def controls(rows,into):
 n=len(rows);full=(1<<n)-1;d=(n+1)//2;q=maxmatch(rows,into)
 alpha,beta,tau,S,E,ell,rho=exterior(rows,into);ab=mul(alpha,beta)
 powers=[{0:1}]
 for j in range(d+2):powers.append(mul(powers[-1],tau))
 need(not powers[q+1],f'matching nilpotence {rows,into}')
 A=[plus(mul(alpha,S[i]),mul(beta,E[i]),-1) for i in range(n)]
 B=[mul(ab,plus(ell[i],rho[i],-1)) for i in range(n)]
 c=[[mul(A[i],powers[j]).get(full,0)+(j*mul(B[i],powers[j-1]).get(full,0) if j else 0) for j in range(d+2)] for i in range(n)]
 need(all(v[d+1]==0 for v in c),'raw cutoff')
 perfect=(2*q==n)
 if not perfect:
  j=q+1
  for i in range(n):
   signed=c[i][j]*(1 if into>>i&1 else -1)
   need(signed>=0,'nonperfect side sign')
   if maxmatch(rows,into,i)==q:need(signed>=math.factorial(q+1),'unmatched strict bound')
 elif q>=2:
  j=q;J={}
  for u in range(n):
   for v in range(n):
    if not(into>>u&1) and into>>v&1 and rows[u]>>v&1:J[(1<<u)|(1<<v)]=1
  Jp={0:1}
  for _ in range(q-1):Jp=mul(Jp,J)
  for i in range(n):
   residual=powers[q-1].get(full^(1<<i),0)
   if into>>i&1:
    residual+=sum(Jp.get(full^(1<<i)^(1<<u),0) for u in range(n) if not(into>>u&1) and rows[i]>>u&1)
    need(c[i][q]==q*residual,'perfect exact positive formula')
   else:
    residual+=sum(Jp.get(full^(1<<i)^(1<<v),0) for v in range(n) if into>>v&1 and rows[v]>>i&1)
    need(c[i][q]==-q*residual,'perfect exact negative formula')
  for side in [0,1]:
   marks=[i for i in range(n) if bool(into>>i&1)==bool(side)]
   need(any(c[i][q] for i in marks),'perfect strict side')
   for a,b in itertools.combinations(marks,2):need(c[a][q]+c[b][q]!=0,'perfect two-name strict group')
 else:
  j=0;need(sorted(v[0] for v in c)==[-2,2],'perfect q1 exceptional constants')
 sides=[[sum(c[i][j] for i in range(n) if bool(into>>i&1)==bool(side)) for j in range(d+1)] for side in [0,1]]
 need(any(x for s in sides for x in s),'nonzero homogeneous group')
 # All w-moments determine the raw jet, exactly.
 for row in c:
  w=[sum(math.comb(h,j-h)*row[j] for j in range(h,min(2*h,d)+1)) for h in range(d+1)]
  recovered=[0]*(d+1)
  for h in range(d,-1,-1):
   recovered[h]=w[h]-sum(math.comb(h,j-h)*recovered[j] for j in range(h+1,min(2*h,d)+1))
  need(recovered==row[:d+1],'unitriangular jet inverse')
 return {'n':n,'rows':rows,'into':into,'matching':q,'perfect':perfect,'raw_moments':[r[:d+1] for r in c],'side_moments':sides}

def physical_graph(m,erows,into,reverse=False):
 n=m+len(erows);rows=[0]*n
 for i in range(m):
  for j in range(i+1,m):
   if reverse:rows[j]|=1<<i
   else:rows[i]|=1<<j
 for u,row in enumerate(erows):rows[m+u]=row<<m
 for u in range(len(erows)):
  for i in range(m):
   if into>>u&1:rows[m+u]|=1<<i
   else:rows[i]|=1<<(m+u)
 return rows

def physical_controls():
 examples=[([0],0),([0],1),([2,0],2),([0,1],2),([2,4,1],2),([2,4,1],3),([14,12,8,0],12)]
 out=[]
 for erows,into in examples:
  n=len(erows);d=(n+1)//2;m=2*d+3;k=(m-1)//2;rs=list(range(2,k+1))+[m-1]
  ctrl=controls(erows,into)
  left=ThreeCover(physical_graph(m,erows,into));right=ThreeCover(physical_graph(m,erows,into,True))
  checks=[]
  for i in range(n):
   for r in rs:
    actual=left.quartic(0,1,r,m+i)-right.quartic(0,1,r,m+i)
    pred=sum(a*b for a,b in zip(formula(m,r),ctrl['raw_moments'][i]))
    need(actual==pred,'literal physical quartic agreement')
    checks.append({'i':m+i,'r':r,'quartic_difference':actual})
  need(any(c['quartic_difference'] for c in checks),'literal physical separation')
  out.append({'m':m,'exterior':ctrl,'checks':checks})
 return out

def main():
 polynomial=polynomial_controls();summary=[];examples=[]
 for n in range(1,6):
  count=0;perfect=0
  for mask in range(1<<(n*(n-1)//2)):
   rows=graph(n,mask)
   for into in range(1<<n):
    r=controls(rows,into);count+=1;perfect+=r['perfect']
    if mask==0 and into in [0,1,(1<<n)-1]:examples.append(r)
  summary.append({'n':n,'cases':count,'perfect_matching_cases':perfect})
  print(json.dumps(summary[-1]),flush=True)
 physical=physical_controls()
 # Unrestricted group counterexample: the two incidence sides in C3 cancel.
 cancellation=controls([2,0],2)
 need(all(sum(v[j] for v in cancellation['raw_moments'])==0 for j in range(2)),'mixed-incidence group counterexample')
 receipt={'novelty':'Restricted selected mixed-quartic battery, full kernel span/jet description, and constructive matching witness. Generic local module-converse detectability is prior incoming E4.','status':'PASS','kernel_span_controls':polynomial,'exhaustive_exterior_controls':summary,'total_exteriors':sum(x['cases'] for x in summary),'physical_checks':physical,'mixed_incidence_group_counterexample':cancellation,'scope':'All-order proofs are in PROOF.md. Finite controls validate full kernel rank, reflection, exact matching nilpotence, signed highest moments, strict homogeneous groups, finite-jet inversion, and literal physical quartics. No quotient-only or arbitrary-ambient source-root claim.','source_sha256':{str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [BASE/'PROOF.md',BASE/'verify.py',HERE/'PROOF.md',HERE/'verify.py']}}
 (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'status':'PASS','total_exteriors':receipt['total_exteriors'],'literal_physical_quartics':sum(len(x['checks']) for x in physical)}))
if __name__=='__main__':main()
