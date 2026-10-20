#!/usr/bin/env python3
"""Independent tuple-word census and integer multiplication-matrix norms."""
import argparse,json,math
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument("--all-shapes",action="store_true");args=ap.parse_args()
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise ArithmeticError(msg)
def edges(w):
 for i in range(len(w)):
  j=(i+1)%len(w)
  if w[i]!=w[j]:
   v=list(w);v[i],v[j]=v[j],v[i]
   yield i,tuple(v)
def necklace(w):return min(w[i:]+w[:i] for i in range(len(w)))
def B(w):
 positions=[i for i,b in enumerate(w) if b]
 return sum(2**i*3**j for j,i in enumerate(reversed(positions)))
def bareiss(A):
 A=[row[:] for row in A];n=len(A);sgn=1;prev=1
 for k in range(n-1):
  if A[k][k]==0:
   pivot=next((j for j in range(k+1,n) if A[j][k]),None)
   if pivot is None:return 0
   A[k],A[pivot]=A[pivot],A[k];sgn=-sgn
  pivot=A[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    value=A[i][j]*pivot-A[i][k]*A[k][j]
    need(value%prev==0,'nonintegral Bareiss division')
    A[i][j]=value//prev
   A[i][k]=0
  prev=pivot
 return sgn*A[-1][-1]
def norm(p,coeff):
 A=[[0]*p for unused in range(p)]
 for col in range(p):
  for a,c in coeff.items():
   quotient,row=divmod(col+a,p);A[row][col]+=c*2**quotient
 return bareiss(A)
rows=[];hits=[];paths=outputs=0
shapes=([(p,q) for p in range(1,11) for q in range(p+1,81)] if args.all_shapes else [(p,(3**p).bit_length()-negative) for p in range(1,41) for negative in [True,False]])
for p,q in shapes:
  if q<=p:continue
  w=tuple(((i+1)*p)//q-(i*p)//q for i in range(q));D=2**q-3**p
  n0=necklace(w);n1={necklace(v) for i,v in edges(w)}
  final={};npaths=0
  for i,u in edges(w):
   for j,v in edges(u):final.setdefault(v,(i,j));npaths+=1
  maxg=0;strict=0;nhits=0
  for v,path in final.items():
   b=B(v);g=math.gcd(abs(D),b);maxg=max(maxg,g);nv=necklace(v)
   distance=0 if nv==n0 else 1 if nv in n1 else 2;strict+=distance==2
   if g==abs(D):
    nhits+=1;hits.append(dict(p=p,q=q,word=''.join(map(str,v)),necklace_distance=distance,fixed_point=b//D))
  paths+=npaths;outputs+=len(final)
  rows.append(dict(p=p,q=q,paths=npaths,distinct_outputs=len(final),exact_necklace_distance_two_outputs=strict,max_content=maxg,integer_hits=nhits))
production=json.loads((Path(__file__).parent/('all-shapes.json' if args.all_shapes else 'probe100.json')).read_text())
prod={(r['p'],r['q']):r for r in production['rows']}
for row in rows:
 need(all(prod[row['p'],row['q']][k]==v for k,v in row.items()),'production census differs')
need(all(hit['necklace_distance']==0 for hit in hits),'strict two hit')
barriers=[]
for p,q,coeff,expected,path in [(4,7,{3:1,2:2,1:-3,0:1},343,[0,3]),(28,45,{27:-1,11:1,10:-1,0:3},23221006835657,None)]:
 R=norm(p,coeff);M=abs(R)
 for ell in [2,3]:
  while M%ell==0:M//=ell
 need(R==expected and M>3**p,'wrong exact norm barrier')
 w=tuple(((i+1)*p)//q-(i*p)//q for i in range(q));r=q-p
 # Reconstruct source one labels from the moving-one modular convention.
 labels={((j+(-(-j*r//p))-1)%q):(-j*r)%p for j in range(p)}
 if path is None:
  # Plus leader residue p-1-r, then follower residue p-1.
  leader=next(i for i,a in labels.items() if a==p-1-r)
  follower=next(i for i,a in labels.items() if a==p-1)
  path=[leader,follower]
 words=[w]
 for i in path:
  need(words[-1][i]!=words[-1][(i+1)%q],'illegal barrier path')
  v=list(words[-1]);j=(i+1)%q;v[i],v[j]=v[j],v[i];words.append(tuple(v))
 n0=necklace(w);n1={necklace(v) for i,v in edges(w)}
 need(necklace(words[-1]) not in n1|{n0},'barrier path not genuine necklace distance2')
 barriers.append(dict(p=p,q=q,coefficients=coeff,resultant=R,stripped_norm=M,three_to_p=3**p,path=path,words=[''.join(map(str,v)) for v in words],content=math.gcd(abs(2**q-3**p),B(words[-1]))))
print(json.dumps(dict(status='PASS',checks=checks,universe=('all shapes p<=10,q<=80' if args.all_shapes else 'critical shapes p<=40'),shape_count=len(rows),paths=paths,distinct_outputs=outputs,integer_hits=hits,exact_norm_barriers=barriers,rows=rows),indent=2))
