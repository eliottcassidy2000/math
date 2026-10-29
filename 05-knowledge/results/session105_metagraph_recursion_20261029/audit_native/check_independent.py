#!/usr/bin/env python3
"""Independent small-order adjacency-matrix oracle; production imports are subjects only."""
import itertools as I
import importlib.util
import json
from pathlib import Path
from collections import Counter
from fractions import Fraction
import sys
import math
from hashlib import sha256
ROOT=Path(__file__).resolve().parents[1]
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok: raise RuntimeError(msg)
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def matrix(n,code):
 A=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(I.combinations(range(n),2)):
  A[i][j]=1-((code>>k)&1); A[j][i]=1-A[i][j]
 return tuple(map(tuple,A))
def code(A):return sum((1-A[i][j])<<k for k,(i,j) in enumerate(I.combinations(range(len(A)),2)))
def relabel(A,p):return tuple(tuple(A[i][j] for j in p) for i in p)
def aut(A):return tuple(p for p in I.permutations(range(len(A))) if relabel(A,p)==A)
C={}
def canon(A):
 if A not in C:
  orbit={relabel(A,p) for p in I.permutations(range(len(A)))}; k=min(map(code,orbit))
  for B in orbit:C[B]=k
 return C[A]
def mate(A):
 n=len(A)
 return tuple(tuple(A[i][j] if abs(i-j)<=1 else 1-A[i][j] for j in range(n)) for i in range(n))
def words(n):
 for x in range(1<<(n*(n-1)//2)):
  A=matrix(n,x)
  if all(A[i][i+1] for i in range(n-1)):yield A
def extension(A):
 n=len(A)
 for s in I.product((0,1),repeat=n-1):
  bits=s+(1,)
  yield tuple(A[i]+(bits[i],) for i in range(n))+(tuple(1-b for b in bits)+(0,),)
def folded(A):
 rev=tuple(tuple(A[j][i] for j in range(len(A))) for i in range(len(A)))
 return min(canon(A),canon(rev))
def signature(A,fold=False):
 return Counter(tuple(sorted((folded(B),folded(mate(B))))) if fold else (canon(B),canon(mate(B))) for B in extension(A))

def psd(A,label):
 A=[list(row) for row in A]
 for k in range(len(A)):
  need(A[k][k]>=0,label+' negative pivot')
  if A[k][k]==0:
   need(all(A[k][j]==0 for j in range(k+1,len(A))),label+' zero-pivot cross term')
   continue
  for i in range(k+1,len(A)):
   for j in range(i,len(A)):
    A[i][j]-=A[i][k]*A[k][j]/A[k][k];A[j][i]=A[i][j]

def proposal_mix(n,W,orders,H):
 classes=sorted(orders); m=len(classes);fac=math.factorial(n)
 weights=[Fraction(1,orders[c]) for c in classes]; mass=sum(weights); mu=[v/mass for v in weights]
 R=[[Fraction(W[a,b],fac)+(Fraction(1)-Fraction(H[a],fac) if a==b else 0) for b in classes] for a in classes]
 K=[]
 for a in classes:
  A=matrix(n,a); targets=Counter()
  for v in range(n):
   others=[j for j in range(n) if j!=v]
   for bits in I.product((0,1),repeat=n-1):
    B=[list(row) for row in A]
    for j,t in zip(others,bits):B[v][j]=t;B[j][v]=1-t
    targets[canon(tuple(map(tuple,B)))]+=1
  K.append([Fraction(targets[b],n*(1<<(n-1))) for b in classes])
 need(sum(mu[i]*Fraction(H[a],fac) for i,a in enumerate(classes))==Fraction(1,1<<(n-1)),'stationary native acceptance')
 for T,label in [(R,'proposal'),(K,'refresh')]:
  for i in range(m):
   need(sum(T[i])==1,label+' rows')
   for j in range(m):need(mu[i]*T[i][j]==mu[j]*T[j][i],label+' balance')
  psd([[mu[i]*T[i][j] for j in range(m)] for i in range(m)],label+' PSD')
 for i in range(m):need(R[i][i]>=Fraction(1,2),'proposal half-lazy')
 if n>=3:
  for alpha in (Fraction(1,4),Fraction(1,2),Fraction(1)):
   L=[[alpha*K[i][j]+(1-alpha)*R[i][j] for j in range(m)] for i in range(m)]
   bound=3*alpha/n
   # Weighted self-adjoint form of (1-bound)I-L+bound*stationary_projection.
   psd([[mu[i]*((1-bound)*(i==j)-L[i][j]+bound*mu[j]) for j in range(m)] for i in range(m)],'mixed gap')
 if n==3:need(R==[[Fraction(5,6),Fraction(1,6)],[Fraction(1,2),Fraction(1,2)]],'n3 proposal entries')
 return {'acceptance':str(Fraction(1,1<<(n-1))),'mixed_gap_alphas':['1/4','1/2','1'] if n>=3 else []}

def native():
 subject=load('native_subject',ROOT/'metagraph/check_weighted_transport.py')
 records=json.loads((ROOT/'metagraph/checks.json').read_text()); summary={}; prev=None
 for n in range(2,6):
  M=Counter();B=Counter(); fibres=Counter()
  for A in words(n):
   pair=(canon(A),canon(mate(A)));M[pair]+=1;fibres[pair[0]]+=1
   if all(A[i][j]==A[n-1-j][n-1-i] for i in range(n) for j in range(n)):B[pair]+=1
  need(M==Counter({tuple(r[:2]):r[2] for r in records[str(n)]['M']}),'native M mismatch')
  need(B==Counter({tuple(r[:2]):r[2] for r in records[str(n)]['B']}),'native B mismatch')
  W=Counter();orders={};H={}
  for c in sorted(fibres):
   A=matrix(n,c);orders[c]=len(aut(A));H[c]=0
   for p in I.permutations(range(n)):
    if all(A[p[i]][p[i+1]] for i in range(n-1)):
     H[c]+=1;W[c,canon(mate(relabel(A,p)))]+=1
   need(H[c]==orders[c]*fibres[c],'independent fibre mass')
  for a in fibres:
   for b in fibres:
    need(M[a,b]==M[b,a],'independent M symmetry')
    need(W[a,b]==orders[a]*M[a,b],'independent witness factor')
    need(Fraction(W[a,b],orders[a])==Fraction(W[b,a],orders[b]),'independent detailed balance')
  if prev is not None:
   total=Counter()
   for A in words(n-1):total.update(signature(A))
   need(total==M,'independent strip recurrence')
  prev=M;summary[n]={'classes':len(fibres),'mass':sum(M.values()),'blue_mass':sum(B.values()),'proposal_mix':proposal_mix(n,W,orders,H)}
 for n in (2,3):
  grouped={}
  for A in words(n):
   key=(canon(A),canon(mate(A))); sig=signature(A)
   need(key not in grouped or grouped[key]==sig,'earlier pair obstruction')
   grouped[key]=sig
 for x,y,f in [(4,18,False),(2,6,True)]:
  A,D=matrix(4,x),matrix(4,y)
  key=lambda T:tuple(sorted((folded(T),folded(mate(T))))) if f else (canon(T),canon(mate(T)))
  need(key(A)==key(D),'hostile source differs');need(signature(A,f)!=signature(D,f),'hostile does not obstruct')
 try:subject.require(False,'guard probe')
 except RuntimeError:pass
 else:raise RuntimeError('native guard disabled')
 return summary

def generator():
 sys.argv=[sys.argv[0],'--max-order','5','--repeats','1','--baseline-root',str(ROOT/'vendor')]
 g=load('generator_subject',ROOT/'generator/check_preparation.py');summary=[]
 # Every raw parent through four: independent complete permutation isomorphism sets.
 for n in range(5):
  parents=[tuple(sum(A[i][j]<<j for j in range(n)) for i in range(n)) for A in (matrix(n,x) for x in range(1<<(n*(n-1)//2)))]
  E,states,data,counts=g.prepare(parents,n)
  expected={}
  for c,raw,incoming,group in data:
   A=tuple(tuple(raw[i]>>j&1 for j in range(n)) for i in range(n)); D0=matrix(n,((1<<(n*(n-1)//2))-1)^c); normalized=tuple(sum(D0[i][j]<<j for j in range(n)) for i in range(n))
   D=tuple(tuple(normalized[i]>>j&1 for j in range(n)) for i in range(n))
   allmaps=tuple(p for p in I.permutations(range(n)) if relabel(A,p)==D)
   need(tuple(incoming)==allmaps,'incomplete canonicalizing torsor')
   need(tuple(group)==aut(D),'incomplete seeded group')
   need(E.cache[normalized]==(c,aut(D)),'incorrect cache seed')
   for s in I.product((0,1),repeat=n):
    orbit=[tuple(s[p[i]] for i in range(n)) for p in aut(D)]
    representative=min(sum(bit<<i for i,bit in enumerate(t)) for t in orbit)
    key=(c,representative)
    expected[key]=tuple(normalized[i]+((1-(representative>>i&1))<<n) for i in range(n))+(representative,)
  need(states==expected,'independent mask bank mismatch')
  summary.append({'parent_order':n,'raw_parents':len(parents),'pointed_states':len(states)})
 parents=[()]
 for n in range(5):
  result=g.build(parents,n,True); classkeys=[]
  for raw in result['bank']:
   A=tuple(tuple(raw[i]>>j&1 for j in range(n+1)) for i in range(n+1));classkeys.append(canon(A))
  allclasses={canon(matrix(n+1,x)) for x in range(1<<(n*(n+1)//2))}
  need(len(classkeys)==len(set(classkeys)) and set(classkeys)==allclasses,'independent full child bank')
  for key,raw in result['states'].items():
   A=tuple(tuple(raw[i]>>j&1 for j in range(n+1)) for i in range(n+1))
   need(result['automorphisms'][key]==aut(A),'complete child automorphisms')
   for vertex,(target,transports) in enumerate(result['root_transports'][key]):
    r=result['states'][target];D=tuple(tuple(r[i]>>j&1 for j in range(n+1)) for i in range(n+1))
    expected=tuple(p for p in I.permutations(range(n+1)) if p[-1]==vertex and relabel(A,p)==D)
    need(transports==expected,'complete aligned rooted transports')
  parents=list(result['bank'])
 for thunk in [lambda:g.old.ParentCanonicalizer(max_order=0).maps((0,)),lambda:g.base.check_row((('a',1),),'a',{'a':None},2)]:
  try:thunk()
  except RuntimeError:pass
  else:raise RuntimeError('inherited guard disabled')
 hostile=g.old.blowup_hostile(); rows=g.base.hostile_tests()
 return {'independent_smallbank':summary,'full_output_orders':[1,2,3,4,5],'inherited_81_27_replay':hostile,'inherited_row_guards':rows}
if __name__=='__main__':
 out={'native':native(),'generator':generator(),'checks':checks,'status':'PASS','python_optimized':not __debug__,'source_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [ROOT/'generator/check_preparation.py',ROOT/'generator/TRANSPORT_SEEDED_BANK.md',ROOT/'metagraph/check_weighted_transport.py',ROOT/'metagraph/WEIGHTED_NATIVE_RECURSION.md',ROOT/'updown/UPDOWN_HARMONIC_AUDIT.md',ROOT/'vendor/SOURCE_MANIFEST.json']}}
 print(json.dumps(out,sort_keys=True,indent=2))
