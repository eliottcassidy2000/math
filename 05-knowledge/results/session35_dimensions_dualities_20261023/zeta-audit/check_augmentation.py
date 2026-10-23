from itertools import combinations, permutations, product
from collections import Counter,defaultdict
from math import gcd
from functools import reduce
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

def pf4(K,ids=(0,1,2,3)):
 a,b,c,d=ids
 return K[a][b]*K[c][d]-K[a][c]*K[b][d]+K[a][d]*K[b][c]
def pf_even(K):
 n=len(K)
 if n==0:return 1
 if n==2:return K[0][1]
 if n==4:return pf4(K)
 if n==6:return sum((-1)**(j+1)*K[0][j]*pf4(K,tuple(i for i in range(1,6) if i!=j)) for j in range(1,6))
 raise ValueError(n)
def form(M):
 return [[M[i][j]-M[i][-1]+M[j][-1] for j in range(len(M)-1)] for i in range(len(M)-1)]
def hp_signed(M):
 H=R=0;n=len(M)
 for p in permutations(range(n)):
  if all(M[p[i]][p[i+1]]==1 for i in range(n-1)):
   H+=1;R+=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
 return H,R

def verify(M):
 n=len(M);K=form(M);P=pf_even(K)
 cof=[];B=s.zeros(n,n-1)
 for i in range(n-1):B[i,i]=1;B[n-1,i]=-1
 for i in range(n):
  minor=[[M[r][c] for c in range(n) if c!=i] for r in range(n) if r!=i]
  cof.append((-1)**i*pf_even(minor))
  assert B.extract([r for r in range(n) if r!=i],range(n-1)).det()==(-1)**i
 assert P==sum(cof)
 mat=s.Matrix(M);w=s.Matrix(cof);km=s.Matrix(K)
 assert mat.adjugate()==w*w.T
 assert km.det()==P*P==(mat+s.ones(n)).det()==(s.ones(1,n)*(-mat).adjugate()*s.ones(n,1))[0]
 snf=smith_normal_form(km,domain=ZZ);diag=sorted(abs(int(v)) for v in snf.diagonal())
 assert all(diag[i]==diag[i+1] for i in range(0,n-1,2))
 return P,tuple(diag[::2]),hp_signed(M)
for out in [[18,4,9,19,6],[30,4,8,18,6]]:
 M=[[0 if i==j else 1 if out[i]>>j&1 else -1 for j in range(5)] for i in range(5)]
 assert all(M[i][j]==-M[j][i] for i in range(5) for j in range(5))
 print('EXACT witness',out,verify(M))
for n,R in [(5,{1,2}),(7,{1,2,4})]:
 M=[[0 if i==j else 1 if (j-i)%n in R else -1 for j in range(n)] for i in range(n)]
 print('EXACT cyclic/Paley',n,verify(M))
# Complete fixed-path small scans with independent SNF and literal Hamiltonian permutations.
for n in [3,5]:
 free=[(i,j) for i,j in combinations(range(n),2) if j!=i+1];prof=defaultdict(set)
 for signs in product([-1,1],repeat=len(free)):
  M=[[0 if i==j else 1 if j==i+1 else -1 if i==j+1 else 0 for j in range(n)] for i in range(n)]
  for (i,j),sign in zip(free,signs):M[i][j]=sign;M[j][i]=-sign
  P,ds,(H,R)=verify(M);prof[ds].add(H)
 print('EXACT fixed-path SNF/permutation scan',n,{k:sorted(v) for k,v in sorted(prof.items())})
# Independent closed Pf4/Pf6 full n7 profile census, no imported producer implementation.
n=7;free=[(i,j) for i,j in combinations(range(n),2) if j!=i+1];counts=Counter()
for signs in product([-1,1],repeat=len(free)):
 M=[[0]*n for _ in range(n)]
 for i in range(n-1):M[i][i+1]=1;M[i+1][i]=-1
 for (i,j),sign in zip(free,signs):M[i][j]=sign;M[j][i]=-sign
 K=form(M);p=abs(pf_even(K));g=reduce(gcd,(pf4(K,I) for I in combinations(range(6),4)))
 assert reduce(gcd,(abs(x) for row in K for x in row))==1
 counts[(1,g,p//g)]+=1
expected_d={1,3,5,7,9,11,13,15,17,19,21,23,25,27,31,33,35}
assert set(counts)=={(1,1,d) for d in expected_d}|{(1,7,7)}
assert counts[(1,7,7)]==9 and sum(counts.values())==32768
assert len({a*b*c for a,b,c in counts})==len(counts)
print('EXACT independent full n7 profiles',dict(sorted(counts.items())))
print('ALL AUGMENTATION AUDIT CHECKS PASS')
