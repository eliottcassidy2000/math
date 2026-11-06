"""Independent physical-digit, language, finite-carry, and divisor audit.
No author-code imports; finite controls accompany separately reviewed proofs.
"""
from pathlib import Path
from math import isqrt
from itertools import product
from functools import lru_cache
import hashlib,json
ROOT=Path(__file__).resolve().parent

def need(x,m):
 if not x:raise RuntimeError(m)
def mm(a,b):return tuple(tuple(sum(a[i][k]*b[k][j]for k in range(len(b)))for j in range(len(b[0])))for i in range(len(a)))
def mv(a,b):return tuple(sum(x*y for x,y in zip(row,b))for row in a)
T=(1,0,2,1)
O=((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1))
Z=((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1))
J=((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0))
Q=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))
OPS={'O':O,'Z':Z,'J':J}
def word(w,s=T):
 for c in w:s=mv(OPS[c],s)
 return s
@lru_cache(None)
def g(u):
 need(u>=0,'nonnegative domain')
 if u<=1:return 1
 return word(''.join('O'if b=='1'else'Z'for b in bin(u-1)[3:]))[0]
def pair_sums(u,padding=0):
 m=u-1;bits=[(m>>i)&1 for i in range(max(1,m.bit_length())+padding)];L=len(bits)
 f1=1+sum((1-bits[i])*((m>>(i+1))+(bits[i+1]if i+1<L else 0))for i in range(L))
 f2=u-bits[0]-sum(bits[i]*bits[j]*2**(j-i-1)for i in range(L)for j in range(i+1,L))-sum(bits[i]*bits[i+1]for i in range(L-1))
 return f1,f2
T0=((1,0,0,0,0),(0,1,1,1,0),(0,0,2,0,0),(0,0,0,0,0),(0,0,0,0,1))
T1=((0,1,1,1,0),(0,1,0,0,0),(0,0,2,0,1),(0,0,0,0,1),(0,0,0,0,1))
def V(n):return g(n),g(n+1),n,n%2,1
for n in range(16385):
 need(mv(T0,V(n))==V(2*n)and mv(T1,V(n))==V(2*n+1),'five-coordinate transitions')
 s=V(0)
 for d in '000'+bin(n)[2:]:s=mv(T1 if d=='1'else T0,s)
 need(s==V(n),'MSB-first vector and common zero padding')
 if n:
  need(pair_sums(n)==(g(n),g(n))and pair_sums(n,3)==(g(n),g(n)),'two weighted formulas and high zero padding')
need(mv(T0,V(0))==V(0),'zero extension stable under zero digit')
# Independent kernel coefficient reconstruction, including n=0.
kernel_checks=0
for e in range(9):
 for r in range(1<<e):
  a=tuple(tuple(int(i==j)for j in range(5))for i in range(5))
  for digit in(format(r,'0'+str(e)+'b')if e else''):a=mm(T1 if digit=='1'else T0,a)
  for n in(0,1,2,7,13):
   need(mv(a,V(n))[0]==g((n<<e)+r),'integer kernel row identity');kernel_checks+=1
# Residue automata are finite and leading-zero invariant.
for mod in(1,2,3,5,8,17):
 for n in range(513):
  v=tuple(x%mod for x in V(0))
  for d in '00'+bin(n)[2:]:v=tuple(x%mod for x in mv(T1 if d=='1'else T0,v))
  need(v[0]==g(n)%mod,'finite residue automaton')
sparse=0
for a in range(2,81):
 for b in range(1,81):
  u=(1<<(a+b))+(1<<b)+1
  need(g(u)==u-(1<<(a-1)),'sparse formula');sparse+=1
for b in range(1,33):
 u=(1<<(b+1))+(1<<b)+1;need(u-g(u)==2,'adjacent bits excluded a=1')
for a in range(2,33):
 u=(1<<a)+2;need(u-g(u)==(1<<(a-1))+1,'b=0 boundary excluded')

ALPH={'A':(1,1),'Z':(0,0),'B':(1,0),'C':(0,1)}
def decode(w):
 u=y=0
 for c in w:a,b=ALPH[c];u=2*u+a;y=2*y+b
 return u,y
separations=0
for k in range(1,97):
 for m in range(1,97):
  text='A'+'Z'*k+'BC'+'Z'*(m-1)+'A';u,y=decode(text)
  need((u,y)==((1<<(k+m+2))+(1<<(m+1))+1,(1<<(k+m+2))+(1<<m)+1),'prefix/suffix decoding')
  need((g(u)==y)==(k==m),'pairwise prefix distinction')
  need(decode('ZZZ'+text)==(u,y),'arbitrary common leading zeros');separations+=1
language=0
for k in range(1,25):
 for r in range(1,25):
  for s in range(25):
   text='A'+'Z'*k+'B'+'C'*r+'Z'*s+'A';u,y=decode(text)
   need(u-y==1<<(s+1)and u-g(u)==1<<k,'exact regular-shape deficits')
   need((g(u)==y)==(s==k-1),'exact regular-language intersection');language+=1
# The K prefixes are distinguishable on the length-limited universe.
for K in range(1,33):
 for k in range(1,K+1):
  for ell in range(1,K+1):
   s='A'+'Z'*ell+'BC'+'Z'*(k-1)+'A'
   need(len(s)<=2*K+3 and (g(decode(s)[0])==decode(s)[1])==(k==ell),'finite-height separation bound')

# Fixed signed-coefficient relation with finite LSD-first carry.
# Initial carry includes delta and the one-time input parity contribution.
def carry_accept(v,y,t,coeff,padding):
 a,b,c,d=coeff;carry=d+c*(v&1);bound=max(abs(d)+abs(c),abs(a)+abs(b)+1)
 width=max(1,v.bit_length(),y.bit_length(),t.bit_length())+padding
 need(abs(carry)<=bound,'initial signed carry bound')
 for i in range(width):
  z=carry+a*((y>>i)&1)+b*((v>>i)&1)-((t>>i)&1)
  if z%2:return False
  carry=z//2;need(abs(carry)<=bound,'uniform signed carry bound')
 return carry==0
coefficients=[(1,0,0,0),(-1,2,0,0),(3,5,-2,-7),(-3,-5,2,7),(2,-1,4,-3),(0,1,-1,0),(7,0,-9,11)]
carry_controls=0
for coeff in coefficients:
 a,b,c,d=coeff
 for v,y,t in product(range(16),repeat=3):
  for pad in(0,3):
   need(carry_accept(v,y,t,coeff,pad)==(t==a*y+b*v+c*(v%2)+d),'signed carry recognizes exact linear relation');carry_controls+=1
# alpha=0 is deliberately allowed in this relation test but not in the recovery theorem.

# Exact positive coefficient of g(v) for every one of the twelve templates.
menu=json.loads((ROOT/'dependencies/MENUS.json').read_text())['5'];alpha_coeffs={}
for item in menu:
 w=item['word'];mat=tuple(tuple(int(i==j)for j in range(4))for i in range(4))
 for c in w:mat=mm(OPS[c],mat)
 A=mm(Q,mat);need(A==tuple(map(tuple,item['matrix'])),'reconstructed twelve matrix')
 col=tuple(A[i][0]for i in range(4));alpha_coeffs[w]=col
 # Express p via 1-p to expose a lower bound from g>=1,u>=2,0<=p<=1.
 cA,cp,cu,c1=col
 lower=cA+2*cu+c1+min(0,cp)
 need(cA>=0 and cu>=0 and lower>0,'all-order exact alpha positivity')
 for u in range(2,257):
  tau=(g(u),u%2,u,1);actual=sum(x*y for x,y in zip(col,tau))
  s=word(w[::-1].translate(str.maketrans('OZ','ZO')),tau)
  need(actual==s[3] and actual>=lower,'actual adjoint d coefficient')
need(alpha_coeffs['J']==(1,0,0,0),'J coefficient g(u)')
need(g(2)==1,'fixed u=2 J row')

# Complete factor-pair explanations of the two existing gaps.
def factor(n):
 out={};d=2
 while d*d<=n:
  while n%d==0:out[d]=out.get(d,0)+1;n//=d
  d+=1
 if n>1:out[n]=out.get(n,0)+1
 return out
def divisors(n):
 for d in range(1,isqrt(n)+1):
  if n%d==0:
   yield d,n//d
   if d*d!=n:yield n//d,d
@lru_cache(None)
def fiber(x):
 # Independent unrestricted range: g(u)>=1 implies u<x.
 return tuple(u for u in range(2,x)if u+g(u)==x)
factor_data={n:factor(n)for n in(2647,7939,13231,2134117,6402349,10670581)}
need(factor_data=={2647:{2647:1},7939:{17:1,467:1},13231:{101:1,131:1},2134117:{643:1,3319:1},6402349:{23:1,278363:1},10670581:{10670581:1}},'complete trial factorization')
small_fibers={x:fiber(x)for x in(6,8,156,643,664)}
need(small_fibers=={6:(3,),8:(),156:(),643:(),664:()},'exact source sum fibers')
# Select every residue-compatible divisor pair, then reject using a cheap empty fiber.
gates={}
for n in(2647,2134117):
 for c in(3,5):
  cases=[];witnesses=[]
  if n%c==0:
   for x,y in divisors(n//c):
    if x%2 and y%2:cases.append(('ee',x,y))
  for x,d in divisors(n):
   if x%2 and d%(2*c)==2*c-1:cases.append(('eo',x,(d+1)//c))
  for d,e in divisors(c*n-c+1):
   if d%(2*c)==e%(2*c)==2*c-1:cases.append(('oo',(d+1)//c,(e+1)//c))
  for parity,x,y in cases:
   # Evaluate the smaller fiber first to avoid an irrelevant huge empty-side scan.
   low,high=sorted((x,y));lo=fiber(low)
   if lo and fiber(high):witnesses.append((parity,x,y))
  need(not witnesses,'complete two-template exclusion at existing gap')
  gates[str(n)+':'+str(c)]=cases
print(json.dumps(dict(status='PASS',physical_g_and_state_controls_through=16384,kernel_row_controls=kernel_checks,
 sparse_bit_pairs=sparse,distinguishing_pairs=separations,regular_intersection_controls=language,
 finite_height_state_lower_bound='at least K states through length 2K+3',signed_carry_controls=carry_controls,
 alpha_coefficients_on_g_p_u_1=alpha_coeffs,gap_factorizations=factor_data,gap_source_sum_fibers=small_fibers,
 all_residue_compatible_gap_gates=gates,scope='Graph nonregularity only; value support and universal occupancy remain unproved.',
 verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
