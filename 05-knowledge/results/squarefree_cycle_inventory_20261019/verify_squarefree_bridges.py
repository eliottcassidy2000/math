#!/usr/bin/env python3
"""Exact checks for square-zero prime-label and vertex-label carriers.
Run: python3 verify_squarefree_bridges.py [--output verification.json]
No network, floating point, or external mathematical packages are used.
"""
from fractions import Fraction
from itertools import product, permutations, combinations
from math import factorial, comb
from collections import defaultdict
import argparse, json

def require(condition, message):
 if not condition:raise RuntimeError(message)

def mul(p,q):
 out=defaultdict(Fraction)
 for a,u in p.items():
  for b,v in q.items():
   if not a&b:out[a|b]+=u*v
 return {a:v for a,v in out.items() if v}
def log1(p,n):
 q=dict(p);q.pop(0,None);out=defaultdict(Fraction);power={0:Fraction(1)}
 for k in range(1,n+1):
  power=mul(power,q)
  for m,c in power.items():out[m]+=Fraction((-1)**(k+1),k)*c
 return {m:c for m,c in out.items() if c}
def subsets(n):return range(1<<n)
def inds(mask,n):return [i for i in range(n) if mask>>i&1]
def conv(a,b,m):
 c=[Fraction(0)]*m
 for i,u in enumerate(a):
  for j,v in enumerate(b):c[(i+j)%m]+=u*v
 return c

def check_group(m,gs):
 n=len(gs);sums=[sum(gs[i] for i in inds(b,n))%m for b in subsets(n)]
 # Independent ordered, unused-label path DP in the regular Cayley representation.
 walk={(0,0):1};exp={0:1};trace_log={}
 for k in range(1,n+1):
  nxt=defaultdict(int)
  for (mask,g),c in walk.items():
   for i,a in enumerate(gs):
    if not(mask>>i&1):nxt[mask|1<<i,(g+a)%m]+=c
  walk=nxt
  for (b,g),c in walk.items():
   if g==0:exp[b]=Fraction(c,factorial(k));trace_log[b]=Fraction(c,k)
 want={b:1 for b in subsets(n) if sums[b]==0}
 require(exp==want, ("subset exponential",m,gs,exp,want))
 require(trace_log=={b:factorial(b.bit_count()-1) for b in want if b}, "cyclic-order trace logarithm")
 # Newton identities in the group algebra; repeated labels are corrected.
 E=[[Fraction(1)]+[Fraction(0)]*(m-1)]
 P=[]
 for j in range(1,n+1):
  row=[0]*m
  for g in gs:row[j*g%m]+=1
  P.append(row)
 for k in range(1,n+1):
  e=[Fraction(0)]*m
  for j in range(1,k+1):
   for g,c in enumerate(conv(E[k-j],P[j-1],m)):e[g]+=Fraction((-1)**(j-1),k)*c
  E.append(e)
 for k in range(n+1):
  direct=[sum(b.bit_count()==k and sums[b]==g for b in subsets(n)) for g in range(m)]
  require(E[k]==direct, ("Newton group-algebra recursion",m,gs,k))
 # Primitive product-one iff EVERY ordering gives distinct proper partial sums.
 checked=0
 for b in want:
  if not b:continue
  ids=inds(b,n)
  primitive=not any(c and c!=b and not(c&~b) and sums[c]==0 for c in subsets(n))
  every=True
  for order in permutations(ids):
   partial=0;seen={0};simple=True
   for i in order[:-1]:
    partial=(partial+gs[i])%m
    if partial in seen:simple=False;break
    seen.add(partial)
   if not simple:every=False;break
  require(primitive==every, ("minimality iff all-order simplicity",m,gs,b))
  checked+=1
 return checked

def sign(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
def det(a):
 n=len(a);ans=0
 for p in permutations(range(n)):
  q=sign(p)
  for i in range(n):q*=a[i][p[i]]
  ans+=q
 return ans

def cycles_on(a,vs):
 if not vs:return 0
 ans=0
 for tail in permutations(vs[1:]):
  w=[vs[0]]+list(tail)
  ans+=all(a[w[i]][w[(i+1)%len(w)]] for i in range(len(w)))
 return ans

def check_vertex(a):
 n=len(a);D={0:1};cycles={}
 for b in range(1,1<<n):
  vs=inds(b,n);minor=[[a[i][j] for j in vs] for i in vs]
  c=(-1)**len(vs)*det(minor)
  if c:D[b]=c
  c=cycles_on(a,vs)
  if c:cycles[b]=c
 require({b:-c for b,c in log1(D,n).items()}==cycles, "vertex-marked determinant cycle inventory")
 return D,cycles

def tournament(n,mask):
 a=[[0]*n for _ in range(n)]
 for j,(u,v) in enumerate(combinations(range(n),2)):
  if mask>>j&1:a[u][v]=1
  else:a[v][u]=1
 return a

def matrix_mul(a,b):
 return [[sum(a[i][k]*b[k][j] for k in range(len(a))) for j in range(len(a))] for i in range(len(a))]
def traces(a):
 p=a;out=[]
 for k in range(1,len(a)+1):out.append(sum(p[i][i] for i in range(len(a))));p=matrix_mul(p,a)
 return out

def binary_code_check(d,columns):
 n=len(columns);r=list(range(1<<d));image={0}
 for v in columns:image|={u^v for u in tuple(image)}
 C=[]
 for b in subsets(n):
  s=0
  for i in inds(b,n):s^=columns[i]
  if s==0:C.append(b)
 dual={sum((((u&v).bit_count()%2)<<i) for i,v in enumerate(columns)) for u in r}
 require(len(C)*len(dual)==1<<n, "kernel-rowspace dimension duality")
 # Exact MacWilliams identity at three rational Bernoulli biases, including 1/2.
 for rho in (Fraction(1,5),Fraction(1,2),Fraction(3,4)):
  direct=sum(rho**b.bit_count()*(1-rho)**(n-b.bit_count()) for b in C)
  rhs=sum((1-2*rho)**w.bit_count() for w in dual)/len(dual)
  require(direct==rhs, ("biased binary character identity",d,columns,rho))
 return len(C)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output');args=ap.parse_args()
 group_cases=0;zero_subsets=0
 for m in range(1,6):
  for n in range(5):
   for gs in product(range(m),repeat=n):
    zero_subsets+=check_group(m,gs);group_cases+=1
 graph_cases=0
 for n in range(1,5):
  for mask in range(1<<comb(n,2)):check_vertex(tournament(n,mask));graph_cases+=1
 witnesses=[]
 for mask in (73,83):
  a=tournament(6,mask);D,cs=check_vertex(a);coeff=[sum(v for b,v in D.items() if b.bit_count()==k) for k in range(7)]
  b=[row+[0] for row in a]+[[1]*6+[0]];D7,cs7=check_vertex(b)
  witnesses.append({'mask':mask,'adjacency':a,'characteristic_coefficients_descending':coeff,'traces':traces(a),'c6':cs.get(63,0),'extension_traces':traces(b),'extension_c6':sum(c for s,c in cs7.items() if s.bit_count()==6),'extension_c7':cs7.get(127,0)})
  graph_cases+=2
 require(witnesses[0]['characteristic_coefficients_descending']==witnesses[1]['characteristic_coefficients_descending'], 'tournament cospectrality')
 require([w['extension_c6'] for w in witnesses]==[0,1], 'sub-Hamiltonian six-cycle discrepancy')
 binary_cases=0
 for d in range(1,4):
  for n in range(5):
   for cols in product(range(1<<d),repeat=n):binary_code_check(d,cols);binary_cases+=1
 negative_log=[]
 for order in range(2,6):
  n=2*order;Z={b:1 for b in subsets(n) if b.bit_count()%order==0};c=log1(Z,n)[(1<<n)-1]
  require(c==1-Fraction(comb(2*order,order),2), "negative logarithm coefficient")
  negative_log.append({'element_order':order,'labels':n,'log_return_top_coefficient':int(c)})
 out={'status':'PASS','group_universe':'All sequences over C_m, 1<=m<=5, 0<=n<=4; repetitions and identity labels allowed','group_cases':group_cases,'zero_return_subsets_checked':zero_subsets,'graph_universe':'All labelled tournaments with 1<=n<=4, plus masks73/83 at n6 and their source extensions','graph_cases':graph_cases,'binary_universe':'All column sequences in F2^d, 1<=d<=3, 0<=n<=4','binary_cases':binary_cases,'negative_log_family':negative_log,'tournament_witnesses':witnesses,'minimality_hostile':{'group':'C5','labels':[1,2,3,4],'simple_order':[1,2,4,3],'proper_zero_subset':[1,4]}}
 text=json.dumps(out,indent=2)
 if args.output:open(args.output,'w').write(text+'\n')
 print(text)
if __name__=='__main__':main()
