#!/usr/bin/env python3
"""Exact univariate Fraction polynomial audit; no SymPy, numerical roots or asserts."""
from fractions import Fraction as F
from itertools import combinations, permutations
from functools import reduce
import json,sys
from math import comb

def norm(p):
 p=list(map(F,p))
 while len(p)>1 and p[-1]==0:p.pop()
 return tuple(p)
def add(p,q):return norm([(p[i] if i<len(p) else 0)+(q[i] if i<len(q) else 0) for i in range(max(len(p),len(q)))])
def neg(p):return tuple(-x for x in p)
def sub(p,q):return add(p,neg(q))
def mul(p,q):
 a=[F(0)]*(len(p)+len(q)-1)
 for i,x in enumerate(p):
  for j,y in enumerate(q):a[i+j]+=x*y
 return norm(a)
def scale(p,x):return norm([x*y for y in p])
def summ(ps):return reduce(add,ps,(F(0),))
def prod(ps):return reduce(mul,ps,(F(1),))
def const(x):return (F(x),)
def at(p,t):return sum(x*t**i for i,x in enumerate(p))
Z=const(0);O=const(1);T=(F(0),F(1));n=5
checks=0
def check(name,got,want=True):
 global checks
 checks+=1
 if got!=want:raise RuntimeError(f'{name}: got {got!r}, expected {want!r}')
def ident():return [[const(i==j) for j in range(n)] for i in range(n)]
def matmul(A,B):return [[summ(mul(A[i][k],B[k][j]) for k in range(n)) for j in range(n)] for i in range(n)]
def transpose(A):return [list(row) for row in zip(*A)]
def madd(A,B):return [[add(A[i][j],B[i][j]) for j in range(n)] for i in range(n)]
def det(A):
 r=len(A);ans=Z
 for perm in permutations(range(r)):
  sign=(-1)**sum(perm[i]>perm[j] for i in range(r) for j in range(i+1,r))
  ans=add(ans,scale(prod(A[i][perm[i]] for i in range(r)),sign))
 return ans

def coeffs(A):
 return [scale(summ(det([[A[i][j] for j in ss] for i in ss]) for ss in combinations(range(n),k)),(-1)**k) for k in range(n+1)]
def zmul(p,q):
 ans=[Z]*(len(p)+len(q)-1)
 for i,x in enumerate(p):
  for j,y in enumerate(q):ans[i+j]=add(ans[i+j],mul(x,y))
 return ans

def hp(A):return summ(prod(A[p[i]][p[i+1]] for i in range(n-1)) for p in permutations(range(n)))
def hdp(A):
 dp={}
 for v in range(n):dp[(1<<v,v)]=O
 for mask in range(1,1<<n):
  for v in range(n):
   cur=dp.get((mask,v),Z)
   for u in range(n):
    if not(mask>>u&1):dp[(mask|1<<u,u)]=add(dp.get((mask|1<<u,u),Z),mul(cur,A[v][u]))
 return summ(dp[((1<<n)-1,v)] for v in range(n))
def tri(A,ss):
 i,j,k=ss
 return add(prod([A[i][j],A[j][k],A[k][i]]),prod([A[i][k],A[k][j],A[j][i]]))
def invariants(A):
 V=summ(mul(A[i][j],A[j][i]) for i,j in combinations(range(n),2))
 tau=summ(tri(A,ss) for ss in combinations(range(n),3))
 M23=summ(prod([A[i][j],A[j][i],tri(A,[k for k in range(n) if k not in (i,j)])]) for i,j in combinations(range(n),2))
 H=hp(A);u=det(madd(ident(),A));v=neg(det(A));delta=add(sub(H,u),v)
 return dict(V=V,tau=tau,M23=M23,H=H,u=u,v=v,Delta=delta)

w=[1,1,-2,0,0];Q=[[const(F(i==j)-F(w[i]*w[j],3)) for j in range(n)] for i in range(n)]
K=[[Z for j in range(n)] for i in range(n)]
for i,j in [(0,3),(3,4),(4,0)]:K[i][j]=T;K[j][i]=neg(T)
L=matmul(matmul(Q,K),transpose(Q))
A=[[Z if i==j else scale(add(O,K[i][j]),F(1,2)) for j in range(n)] for i in range(n)]
B=[[Z if i==j else scale(add(O,L[i][j]),F(1,2)) for j in range(n)] for i in range(n)]
if '--mutate-conjugation' in sys.argv:B[0][3]=add(B[0][3],scale(T,F(1,100)))
check('orthogonal',matmul(Q,transpose(Q)),ident());check('fix ones',[summ(row) for row in Q],[O]*n)
check('conjugacy',matmul(matmul(Q,A),transpose(Q)),B)
check('K regular',[summ(row) for row in K],[Z]*n);check('L regular',[summ(row) for row in L],[Z]*n)
for name,M in [('K',K),('L',L)]:
 for i in range(n):
  for j in range(n):
   check(f'{name} skew {i,j}',add(M[i][j],M[j][i]),Z)
   check(f'{name} linear {i,j}',len(M[i][j])<=2)
   check(f'{name} domain bound {i,j}',sum(abs(c) for c in M[i][j])<=1)
for name,M in [('A',A),('B',B)]:
 check(f'{name} row sums',[summ(row) for row in M],[const(2)]*n)
 check(f'{name} path enumerator vs subset DP',hp(M),hdp(M))
check('cubic skew minimal polynomial',matmul(matmul(K,K),K),[[scale(mul((F(0),F(0),F(1)),K[i][j]),-3) for j in range(n)] for i in range(n)])
Pexpected=zmul(zmul([O,const(-2)],zmul([O,const(F(1,2))],[O,const(F(1,2))])),[O,O,(F(1,4),F(0),F(3,4))])
check('A det polynomial',coeffs(A),Pexpected);check('B det polynomial',coeffs(B),Pexpected)
reduced=zmul(zmul([O,const(F(1,2))],[O,const(F(1,2))]),[O,O,(F(1,4),F(0),F(3,4))])
functional=[Z]*5
for k,c in enumerate(reduced):
 sign=(-1)**k
 if '--mutate-functional-sign' in sys.argv:sign=1
 for j in range(5-k):functional[k+j]=add(functional[k+j],scale(c,sign*comb(4-k,j)))
check('reduced polynomial functional equation',functional,reduced)
# For lambda=-1/2+i*y, |lambda+1|^2=|lambda|^2: after reciprocity,
# every reduced-zeta pole lies on |z+1|=1. This is polynomial arithmetic.
lambda_real=const(F(-1,2));lambda_imag=T
lambda_norm=add(mul(lambda_real,lambda_real),mul(lambda_imag,lambda_imag))
shifted_real=add(lambda_real,O)
shifted_norm=add(mul(shifted_real,shifted_real),mul(lambda_imag,lambda_imag))
check('pole-circle norm difference',sub(shifted_norm,lambda_norm),Z)
a=invariants(A);b=invariants(B)
common=dict(V=(F(5,2),F(0),F(-3,4)),tau=(F(5,2),F(0),F(3,4)),u=(F(3,16),F(0),F(9,16)),v=(F(-1,8),F(0),F(-3,8)))
for key,value in common.items():check('A '+key,a[key],value);check('B '+key,b[key],value)
check('A H',a['H'],(F(15,2),F(0),F(9,4)))
check('A M23',a['M23'],const(F(5,8)))
check('M23 difference',sub(b['M23'],a['M23']),(F(0),F(0),F(0),F(0),F(-1,54)))
expected_diff=(F(0),F(0),F(0),F(0),F(-1,27))
if '--mutate-difference' in sys.argv:expected_diff=scale(expected_diff,-1)
check('H difference',sub(b['H'],a['H']),expected_diff)
check('Delta difference',sub(b['Delta'],a['Delta']),expected_diff)
for name,inv in [('A',a),('B',b)]:
 check(name+' spectral-incidence identity',inv['H'],sub(add(add(O,scale(inv['tau'],2)),scale(inv['M23'],2)),scale(Pexpected[5],2)))

# Exact circle parametrization and sharp product bound. The polynomial
# variable in these three identities is now the slope r, not the weight t.
D=(F(1),F(1),F(1));nums=[(F(0),F(-1)),(F(1),F(1)),(F(0),F(1),F(1))]
check('circle sum identity',summ(nums),D)
check('circle norm identity',summ(mul(x,x) for x in nums),mul(D,D))
N=prod(nums);square_base=prod([(F(-1),F(1)),(F(2),F(1)),(F(1),F(2))])
check('sharp product square identity',add(scale(prod([D,D,D]),4),scale(N,27)),mul(square_base,square_base))
for r in map(F,[-2,-1,0,1,2,3]):
 den=1+r+r*r;svec=[-r/den,(1+r)/den,r*(1+r)/den]
 X=[[Z for j in range(n)] for i in range(n)]
 for i,c in enumerate(svec):X[i][3]=scale(T,c);X[3][i]=scale(T,-c);X[i][4]=scale(T,-c);X[4][i]=scale(T,c)
 X[3][4]=T;X[4][3]=neg(T)
 C=[[Z if i==j else scale(add(O,X[i][j]),F(1,2)) for j in range(n)] for i in range(n)]
 iv=invariants(C);product_s=svec[0]*svec[1]*svec[2]
 check('circle P '+str(r),coeffs(C),Pexpected)
 check('circle M23 '+str(r),iv['M23'],norm([F(5,8),0,0,0,product_s/8]))
 check('circle H '+str(r),iv['H'],norm([F(15,2),0,F(9,4),0,product_s/4]))

samples=[]
for t in [F(-1),F(-1,4),F(0),F(1,4),F(1)]:
 row={'t':str(t),'A':{k:str(at(p,t)) for k,p in a.items()},'B':{k:str(at(p,t)) for k,p in b.items()}}
 check('v sign '+str(t),at(a['v'],t)<0)
 samples.append(row)
for m in [5,7,11,25,35,55,121,143,1001]:
 # Ring Z[1/6] -> Z/m; identity readout differs by a unit.
 d=(-pow(6912,-1,m))%m
 check('unit residue '+str(m),__import__('math').gcd(d,m),1)
print(json.dumps({'status':'PASS','checks':checks,'coefficient_field':'Q[t]','method':'Fraction coefficient arithmetic; permutation determinants; Hamilton permutation and subset-DP paths','H_difference':list(map(str,expected_diff)),'P_coefficients_t':[[str(x) for x in p] for p in Pexpected],'samples':samples},sort_keys=True,indent=2))
