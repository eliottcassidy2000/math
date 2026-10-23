#!/usr/bin/env python3
"""Literal truncated polynomial coefficient DP; independent of sorted-capacity recursions."""
import itertools,collections,math,json
from fractions import Fraction
from pathlib import Path

def require(b,msg):
 if not b:raise RuntimeError(msg)
def coefficient(n,k,symmetric=False):
 D={(0,)*n:1}
 if symmetric:
  factors=[[(i,)] for i in range(n)]+[[(i,j)] for i,j in itertools.combinations(range(n),2)]
 else:factors=[list(itertools.combinations(range(n),k)) for _ in range(n)]
 for terms in factors:
  E=collections.defaultdict(int)
  if symmetric:E.update(D)
  for a,v in D.items():
   for term in terms:
    b=list(a)
    for i in term:b[i]+=1
    if max(b)<=k:E[tuple(b)]+=v
  D=E
 return D.get((k,)*n,0)
def canon(B,p=5):
 def rotate(x,b):return ((x<<b)|(x>>(p-b)))&((1<<p)-1) if b else x
 return min(tuple(rotate(B[(i+a)%p],b) for i in range(p)) for a in range(p) for b in range(p))
def con(B,p=5):return tuple(sum(((B[(-j)%p]>>((-i)%p))&1)<<j for j in range(p)) for i in range(p))
rows=[]
for p in(3,5,7):
 layers=[]
 for k in range(1,(p+1)//2):
  A=coefficient(p,k);S=coefficient(p,k,True);N_num=A+(p-1)**2*math.comb(p,k);F_num=S+(p-1)*math.comb((p-1)//2,k//2)
  require(N_num%(p*p)==0 and F_num%p==0,'Burnside divisibility');N=N_num//(p*p);F=F_num//p;require((N+F)%2==0,'merged parity')
  layers.append({'k':k,'A':A,'S':S,'ordinary':N,'self_converse':F,'merged':(N+F)//2})
 rows.append({'p':p,'layers':layers,'all_layers_ordinary':1+sum(x['ordinary'] for x in layers),'all_layers_self':1+sum(x['self_converse'] for x in layers),'all_layers_merged':1+sum(x['merged'] for x in layers)})
# A different finite universe: unordered union of two disjoint perfect matchings.
P=list(itertools.permutations(range(5)));multiplicity=collections.Counter()
for a in P:
 for b in P:
  if all(x!=y for x,y in zip(a,b)):multiplicity[tuple((1<<x)|(1<<y) for x,y in zip(a,b))]+=1
require(len(multiplicity)==2040,'matrix count from decompositions');require(collections.Counter(multiplicity.values())=={2:1440,4:600},'two-factor decomposition multiplicities')
classes=collections.Counter(canon(B) for B in multiplicity);selfclasses=[B for B in classes if canon(con(B))==B]
require(len(classes)==88 and len(selfclasses)==24,'finite orbit counts');require(len({min(B,canon(con(B))) for B in classes})==56,'finite merged count')
require(collections.Counter(classes.values())=={5:8,25:80},'orbit sizes');require(sum(all(((B[i]>>j)&1)==((B[j]>>i)&1) for i in range(5) for j in range(5)) for B in multiplicity)==112,'symmetric count')
# Complement transpose swaps the two cells and maps degree k to5-k.
for B in multiplicity:
 C=tuple(sum((1-((B[j]>>i)&1))<<j for j in range(5)) for i in range(5))
 require(all(x.bit_count()==3 for x in C),'degree-three swap');require(tuple(sum((1-((C[j]>>i)&1))<<j for j in range(5)) for i in range(5))==B,'swap involution')
# Degree-two generating functions, evaluated through their component recurrences.
a=[Fraction(1),Fraction(0)];s=[1,0];egf=[]
for n in range(2,12):
 a.append(sum(a[:n-1])/(2*n))
 s.append(sum(math.comb(n-1,k-1)*(1 if k==2 else (k+1)*math.factorial(k-1)//2)*s[n-k] for k in range(2,n+1)))
 if n in(3,5,7,11):egf.append({'n':n,'A_n2':int(a[n]*math.factorial(n)**2),'S_n2':s[n]})
for row in rows:
 for layer in row['layers']:
  if layer['k']==2:
   val=next(x for x in egf if x['n']==row['p']);require((layer['A'],layer['S'])==(val['A_n2'],val['S_n2']),'EGF versus coefficient DP')
output={'status':'PASS','method':'literal truncated polynomial multiplication; independent matching-pair universe; component EGFs','rows':rows,'p5_degree2':{'matrices':2040,'cycle10':1440,'cycle4plus6':600,'ordinary':88,'self_converse':24,'merged':56,'orbit_size_distribution':{'5':8,'25':80}},'degree_two_egf':egf}
Path(__file__).with_name('receipt.json').write_text(json.dumps(output,indent=2)+'\n');print(output)
