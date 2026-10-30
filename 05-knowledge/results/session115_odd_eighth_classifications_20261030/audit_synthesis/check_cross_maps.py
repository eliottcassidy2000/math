"""Independent exact controls for the Session115 synthesis maps."""
from itertools import product
from fractions import Fraction as F
from math import comb
from collections import Counter
import json
import sympy as S
checks=0
def gate(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)
def compositions(n):
 if n==0:yield ();return
 for a in range(1,n+1):
  for tail in compositions(n-a):yield (a,)+tail
def address(w,q):
 mod=2**sum(w);z=0;prefix=0
 for i,a in enumerate(w,1):z-=2**prefix*pow(q,-i,mod);prefix+=a
 return z%mod
pi={0:0}
for A in range(1,4):
 for w in compositions(A):pi[address(w,3)*2**(3-A)]=address(w,5)*2**(3-A)
z=(1+S.I)/S.sqrt(2)
chi2={1:1,3:-1,5:-1,7:1};chim2={1:1,3:1,5:-1,7:-1}
for j in (1,3,5,7):
 raw=S.simplify((z**j+z**(7*j))/S.sqrt(2))
 image=S.simplify((z**(pi[1]*j)+z**(pi[7]*j))/S.sqrt(2))
 gate(S.simplify(raw-chi2[j])==0,('real character',j))
 gate(S.simplify(image+S.I*chim2[j])==0,('imaginary image',j))
 original=S.Rational(88,105)+S.Rational(32,105)*raw
 transported=S.Rational(88,105)+S.Rational(32,105)*image
 gate(bool(original>0),('positive source',j))
 gate(S.simplify(S.im(transported))!=0,('nonreal weight',j))
I=(1,2,4)
def act(A,x):
 r=0
 for i in range(3):
  if x>>i&1:r^=A[i]
 return r
def mul(A,B):return tuple(act(A,x) for x in B)
def order(A):
 B=I
 for n in range(1,169):
  B=mul(B,A)
  if B==I:return n
 raise RuntimeError('matrix order')
def affine_order(v,A):
 w=0;B=I
 for n in range(1,337):
  w^=act(B,v);B=mul(B,A)
  if w==0 and B==I:return n
 raise RuntimeError('affine order')
def cyc(A):
 C={I};B=A
 while B!=I:C.add(B);B=mul(B,A)
 return C
GL=[A for A in product(range(1,8),repeat=3) if len({act(A,x) for x in range(8)})==8]
gate(len(GL)==168,'GL order')
odd=[A for A in GL if order(A)%2]
profiles=Counter();normalized=Counter()
for A in odd:
 l=order(A);N=(0,0,0);B=I
 for _ in range(l):N=tuple(a^b for a,b in zip(N,B));B=mul(B,A)
 gate(mul(N,N)==N,('norm projection',A))
 f=[int(act(N,x)!=0) for x in range(8)]
 gate(sum(f)=={1:7,3:4,7:0}[l],('Boolean mass',A))
 gate(all(f[x]==int(affine_order(x,A)%2==0) for x in range(8)),('actual order predicate',A))
 r=sum(f[1<<i]<<i for i in range(3))
 s=[f[x]^((r&x).bit_count()%2) for x in range(8)]
 gate(all(s[1<<i]==0 for i in range(3)),('normalization',A))
 gate(sum(s[x]^((r&x).bit_count()%2) for x in range(8))==sum(f),('named character transport',A))
 if l==1:gate(sum(s)==3 and r==7, 'identity normalized gauge')
 else:gate(not any(s),('linear or zero normalized',A))
 profiles[(l,sum(f))]+=1;normalized[(l,r)]+=1
C7=cyc(next(A for A in GL if order(A)==7));a=next(A for A in C7 if A!=I)
N7={g for g in GL if any(mul(g,a)==mul(b,g) for b in C7)}
C3=cyc(next(A for A in N7 if order(A)==3));masses=[]
for H,expected in [({I},F(7,8)),(C3,F(5,8)),(C7,F(1,8)),(N7,F(3,8))]:
 observed=F(sum(affine_order(x,A)%2==0 for A in H for x in range(8)),8*len(H))
 gate(observed==expected,('group average',len(H)))
 masses.append({'complement_order':len(H),'even_mass':str(observed)})
for L in range(101):
 K,delta=divmod(L,2);count=0
 for j in range(L+1):
  a,eps=divmod(j,2);h=2;residue=True
  while h<=10*(L+1):
   residue &= j%h<=L%h;h*=5
  lucas=eps<=delta and comb(K,a)%5!=0
  gate(residue==lucas,('Lucas',L,j));count+=residue
 n=K;digits=1
 while n:digits*=n%5+1;n//=5
 gate(count==(delta+1)*digits,('phase count',L))
haar=[]
for E in range(1,6):
 M=5**E;vals=[]
 for c in (0,2):
  for j in range(5**(E-1)):
   a=4*j+c;v=(81-pow(16,a,M))%M;r=0
   while r<E and v%5==0:r+=1;v//=5
   vals.append(r)
 for k in range(1,E+1):gate(F(sum(r>=k for r in vals),len(vals))==F(1,5**(k-1)),('Haar finite layer',E,k))
 truncated=(F(sum(vals),len(vals))+F(1,4))/4
 gate(F(3,8)-truncated==F(1,16*5**(E-1)),('Haar tail',E))
 haar.append({'precision':E,'truncated_extended_mean':str(truncated)})
print(json.dumps({'status':'PASS','gates':checks,'odd_GL_matrices':len(odd),'Boolean_profiles':[{'matrix_order':k[0],'truth_weight':k[1],'matrices':v} for k,v in sorted(profiles.items())],'group_masses':masses,'Haar_controls':haar},indent=2))
