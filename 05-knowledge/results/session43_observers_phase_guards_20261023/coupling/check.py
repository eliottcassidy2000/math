"""Exact restricted couplings; Python 3 + SymPy, no floating-point decisions."""
import itertools,collections,json,hashlib
import sympy as s
from explore import classes,chain

def need(ok,message):
 if not ok:raise RuntimeError(message)
def psd(M):
 A=[list(M.row(i)) for i in range(M.rows)]
 for k in range(len(A)):
  p=A[k][k]
  if p<0:return False
  if not p:
   if any(A[k][j] for j in range(k+1,len(A))):return False
   continue
  for i in range(k+1,len(A)):
   for j in range(i,len(A)):
    A[i][j]=A[j][i]=A[i][j]-A[i][k]*A[k][j]/p
 return True

def verify(C,K):
 r=s.Matrix([sum(C.row(i)) for i in range(C.rows)]);c=[sum(C.col(j)) for j in range(C.cols)];D=s.diag(*r)
 need(K*s.ones(K.rows,1)==s.ones(K.rows,1),'stochasticity')
 need(D*K==(D*K).T,'reversibility');need(psd(D*K),'positive semidefinite kernel')
 seen=set();components=[]
 for a in range(K.rows):
  if a in seen:continue
  seen.add(a);todo=[a];block=[]
  while todo:
   i=todo.pop();block.append(i)
   for j in range(K.rows):
    if K[i,j] and j not in seen:seen.add(j);todo.append(j)
  components.append(sorted(block))
 poly=K.charpoly().as_poly(); z=poly.gens[0]; multiplicity=0
 while poly.eval(1)==0:
  poly=poly.exquo(s.Poly(z-1,z));multiplicity+=1
 need(len(components)==multiplicity,'components and fixed space')
 cert=D-D*K-s.Rational(8,9)*(D-r*r.T/sum(r))
 return {'graph_states':K.rows,'tournament_states':C.cols,'components':components,'word_count':int(sum(r)),'charpoly':str(s.factor(K.charpoly().as_expr())),'gap_at_least_eight_ninths':psd(cert),'incidence':[[int(x) for x in C.row(i)] for i in range(C.rows)]}

reports=[]
for n in range(2,6):
 E,G,T=classes(n);path=sum(1<<E.index((i,i+1)) for i in range(n-1));A=[x for x in G if x&path==0]
 def hp(x):
  return sum(all(((x>>E.index(tuple(sorted((p[i],p[i+1])))))&1)==(p[i]>p[i+1]) for i in range(n-1)) for p in itertools.permutations(range(n)))
 hs={t:hp(t) for t in set(T.values())}
 # Independently count automorphisms as number of maps fixing representative.
 for t in hs:
  orbit=sum(T[x]==t for x in T);aut=s.factorial(n)//orbit
  need(sum(T[x]==t for x in A)==hs[t]/aut,'Hamiltonian observer orbit weight')
 # Converse lift: index reflection, no bit complementation.
 for x in A:
  reflected=sum(1<<E.index(tuple(sorted((n-1-a,n-1-b)))) for k,(a,b) in enumerate(E) if x>>k&1)
  need(reflected in A and G[reflected]==G[x] and T[reflected]==T[max(T)^x],'normalized converse')
 for paired in [False,True]:
  gs,ts,C,K=chain(A,G,T,paired);r=verify(C,K);r.update(n=n,model='Hamiltonian observer',paired=paired)
  need(K[0,0]==1 and sum(C.row(0))==1 and sum(C.col(0))==1,'isolated transitive tile')
  if n>=3:need(not r['gap_at_least_eight_ninths'] and len(r['components'])>=2,'observer obstruction')
  if n==3:need(K==s.eye(2),'smallest counterexample')
  if n==4:
   expected=s.Matrix([[s.Rational(11,15),s.Rational(2,15),s.Rational(1,15),s.Rational(1,15)]]+[[s.Rational(1,5),s.Rational(2,5),s.Rational(1,5),s.Rational(1,5)]]*3)
   need(K[1:,1:]==expected,'order-four core matrix')
   need(expected.charpoly().as_expr().expand()==(s.Symbol('lambda')**2*(s.Symbol('lambda')-1)*(s.Symbol('lambda')-s.Rational(8,15))).expand(),'order-four core spectrum')
  reports.append(r)
  print('checked',r['model'],n,file=__import__('sys').stderr)
 for residue in [1,3,5,7]:
  B=[x for x in G if hs[T[x]]%8==residue]
  if not B:continue
  gs,ts,C,K=chain(B,G,T);r=verify(C,K);r.update(n=n,model='modulo-8 conditioning',residue=residue,paired=False)
  if n<=4:need(K*K==K and len(r['components'])==1,'small residue filter rank one')
  if n==5 and residue in [5,7]:
   lam=s.Rational(143903,1801800) if residue==5 else s.Rational(89,315)
   z=s.Symbol('lambda');need(K.charpoly().as_expr().expand()==(z**(K.rows-2)*(z-1)*(z-lam)).expand(),'exact residue eigenvalue')
   need(r['gap_at_least_eight_ninths']==(residue==5),'opposite-direction mod8 controls')
  reports.append(r)
 # Mod8 as observed summary, NOT conditioning: must retain the bound.
 Tc={x:hs[T[x]]%8 for x in T}
 gs,ts,C,K=chain(list(G),G,Tc);r=verify(C,K);r.update(n=n,model='modulo-8 observed coarsening',paired=False)
 need(r['gap_at_least_eight_ninths'],'coarsening bound');reports.append(r)
# Actual subgroups fixing named vertices, preserving the affine action.
for n,fixed,paired in [(4,2,False),(4,3,True)]:
 P=[p for p in itertools.permutations(range(n)) if all(p[i]==i for i in range(fixed))]
 E,G,T=classes(n,P);gs,ts,C,K=chain(list(G),G,T,paired);r=verify(C,K);r.update(n=n,model='individually fixed vertices',fixed_vertices=fixed,paired=paired)
 need(len(r['components'])>1,'marked common-information obstruction');reports.append(r)
# Soft guards keep every word, weights 2 on accepted and 1 elsewhere.
for n,guard in [(4,'observer'),(5,'mod8 residue7')]:
 E,G,T=classes(n);path=sum(1<<E.index((i,i+1)) for i in range(n-1))
 def hp_soft(x):
  return sum(all(((x>>E.index(tuple(sorted((p[i],p[i+1])))))&1)==(p[i]>p[i+1]) for i in range(n-1)) for p in itertools.permutations(range(n)))
 hs={t:hp_soft(t) for t in set(T.values())}
 accept=lambda x: x&path==0 if guard=='observer' else hs[T[x]]%8==7
 weighted_words=[x for x in G for _ in range(2 if accept(x) else 1)]
 gs,ts,C,K=chain(weighted_words,G,T);r=verify(C,K)
 rr=s.Matrix([sum(C.row(i)) for i in range(C.rows)]);D=s.diag(*rr)
 need(psd(D-D*K-s.Rational(4,9)*(D-rr*rr.T/sum(rr))),'soft guard comparison')
 need(len(r['components'])==1,'positive density connectivity')
 r.update(n=n,model='soft '+guard,paired=False,minimum_to_maximum_weight='1/2',certified_gap_lower_bound='4/9')
 reports.append(r)
print(json.dumps({'status':'PASS','arithmetic':'exact integers and SymPy rationals only','reports':reports},indent=2,sort_keys=True))
