"""Independent controls for refined SCC formulas; formal atoms are NOT tournaments."""
from collections import defaultdict
from functools import lru_cache
from math import prod,comb
from pathlib import Path
import hashlib,json
import sympy as sp
N=16
# name,order,h,b,opposite,e3,e5. Only u,v,r,s use certified tournament data.
L=[('u',1,1,1,0,0,0),('v',3,1,1,1,1,0),('r',5,3,3,2,0,1),
   ('s',4,5,3,3,0,0),('a',6,9,1,4,2,0),('b',6,9,3,5,1,1),
   ('c',6,9,0,7,1,0),('d',6,9,0,6,1,0),('e',8,27,9,8,0,2),
   ('f',9,27,27,9,2,0),('g',10,81,3,10,3,1)]
checks=0
def need(c,msg):
 global checks
 checks+=1
 if not c:raise RuntimeError(msg)
def zero():return [0]*(N+1)
def mon(n,c=1):
 v=zero()
 if n<=N:v[n]=c
 return v
def add(a,b):return [x+y for x,y in zip(a,b)]
def mul(a,b):return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(N+1)]
def sub2(a):return [a[n//2] if n%2==0 else 0 for n in range(N+1)]
def wt(i,s=1):return (2*3**L[i][5]*5**L[i][6])**s
P={};C={}
for scale in (1,2):
 P[scale]=defaultdict(zero);C[scale]=defaultdict(zero)
 for i,(_,n,h,b,op,e3,e5) in enumerate(L):
  P[scale][h][n]+=wt(i,scale)
  if i==op:C[scale][h][n]+=wt(i,scale)
 C[scale][1][0]+=1
@lru_cache(None)
def F(h,scale=1):
 if h==1:
  a=zero();a[0]=1
  for n in range(1,N+1):a[n]=wt(0,scale)*a[n-1]+(wt(1,scale)*a[n-3] if n>=3 else 0)
  return a
 out=zero()
 for j,v in P[scale].items():
  if j>1 and h%j==0:out=add(out,mul(v,F(h//j,scale)))
 return mul(F(1,scale),out)
@lru_cache(None)
def J(h):
 out=zero()
 for x in range(1,sp.integer_nthroot(h,2)[0]+1):
  if h%(x*x)==0:out=add(out,mul(sub2(F(x,2)),C[1][h//(x*x)]))
 return out
centers=[(0,1,1,1)]+[(n,h,b,wt(i)) for i,(_,n,h,b,op,_,_) in enumerate(L) if i==op]
@lru_cache(None)
def S(d,full=False):
 out=zero()
 for x in range(1,d+1,2):
  if d%x or x*(x-1)>d:continue
  e=zero()
  for n,h,b,w in centers:
   if x*(x*h-b)==d:e[n]+=w*(7**(x*x*h)*11**(x*b) if full else 1)
  out=add(out,mul(sub2(F(x,2)),e))
 return out
def B(k,full=False):return [s+(f-j)//2*(7**k if full else 1) if k%2 else s for s,f,j in zip(S(2*k,full),F(k),J(k))]
ordinary=defaultdict(zero);selfconv=defaultdict(zero);deficits=defaultdict(zero);degrees=defaultdict(zero)
full_deficits=defaultdict(zero);full_degrees=defaultdict(zero)
word_count=0
# Literal word enumerator uses raw multiplicative aut markers, not formula convolution.
def visit(word,n,h,c,e3,e5):
 global word_count
 word_count+=1;weight=2**c*3**e3*5**e5
 ordinary[h][n]+=weight
 opp=tuple(L[i][4] for i in word[::-1])
 if word==opp:
  x=prod(L[i][2] for i in word[:len(word)//2]);beta=1 if len(word)%2==0 else L[word[len(word)//2]][3]
  d=h-x*beta
  selfconv[h][n]+=weight;deficits[d][n]+=weight;degrees[d//2][n]+=weight
  full_deficits[d][n]+=weight*7**h*11**(x*beta);full_degrees[d//2][n]+=weight*7**h*11**(x*beta)
 elif word<opp:degrees[h][n]+=weight;full_degrees[h][n]+=weight*7**h
 for i,(_,size,j,_,_,a,b) in enumerate(L):
  if n+size<=N:visit(word+(i,),n+size,h*j,c+1,e3+a,e5+b)
visit((),0,1,0,0,0)
for h in range(1,201):
 need(F(h)==ordinary[h],('marked F',h));need(J(h)==selfconv[h],('marked J',h))
for d in range(2,201,2):need(S(d)==deficits[d],('marked S',d))
for k in range(1,101):need(B(k)==degrees[k],('marked B',k))
for d in range(2,201,2):need(S(d,True)==full_deficits[d],('fully marked S',d))
for k in range(1,101):need(B(k,True)==full_degrees[k],('fully marked B',k))
# Actual pure-blue boundary uses only the certified K1,C3,R5,S4 subalphabet.
real_B0=zero()
def real_visit(word,n,h,weight):
 opp=tuple(L[i][4] for i in word[::-1])
 if word==opp:
  x=prod(L[i][2] for i in word[:len(word)//2]);beta=1 if len(word)%2==0 else L[word[len(word)//2]][3]
  if h==x*beta:real_B0[n]+=weight*7**h*11**h
 for i in range(4):
  size=L[i][1]
  if n+size<=N:real_visit(word+(i,),n+size,h*L[i][2],weight*wt(i))
real_visit((),0,1,1)
center0=add(add(add(mon(0,7*11),mon(1,7*11*wt(0))),mon(3,7*11*wt(1))),mon(5,7**3*11**3*wt(2)))
need(mul(sub2(F(1,2)),center0)==real_B0,'fully marked actual B0 boundary')
# A deliberately incorrect half-substitution is caught at the two-triangle word.
wrong=sub2(F(1,1));correct=sub2(F(1,2))
need(wrong[6]!=correct[6],'hostile catches failure to square valuation and SCC markers')
# Exact rational algebra is performed independently of truncated series.
z=sp.symbols('z');p=1-z-z**3;q=1-z*z-z**6
need(sp.Poly(p,z).is_irreducible,'p irreducible Q')
need(sp.Poly(q,z).is_irreducible,'q irreducible Q')
need(sp.gcd(p,q)==1,'p q coprime')
need(sp.Poly(q-(1+z+z**3)**2,z,modulus=2).is_zero,'q mod2 is a square, not irreducible')
Ps=defaultdict(lambda:sp.Integer(0));Cs=defaultdict(lambda:sp.Integer(0))
centerss=[(0,1,1)]
for i,(_,n,h,b,op,_,_) in enumerate(L):
 Ps[h]+=z**n
 if i==op:Cs[h]+=z**n;centerss.append((n,h,b))
Cs[1]+=1
@lru_cache(None)
def fs(h):
 if h==1:return 1/p
 return sp.cancel(sum(Ps[j]*fs(h//j) for j in list(Ps) if j>1 and h%j==0)/p)
@lru_cache(None)
def rank(h):
 if h==1:return 0
 opts=[rank(h//j)+1 for j in list(Ps) if j>1 and h%j==0 and rank(h//j)>=0]
 return max(opts,default=-1)
def js(h):return sp.cancel(sum(fs(x).subs(z,z*z)*Cs[h//(x*x)] for x in range(1,sp.integer_nthroot(h,2)[0]+1) if h%(x*x)==0))
def ss(d):
 terms=[];ranks=[]
 for x in range(1,d+1,2):
  if d%x or x*(x-1)>d:continue
  E=sum(z**n for n,a,b in centerss if x*(x*a-b)==d)
  if E!=0 and rank(x)>=0:terms.append(fs(x).subs(z,z*z)*E);ranks.append(rank(x))
 return sp.cancel(sum(terms)),max(ranks,default=-1)
def denominator(expr):
 D=sp.fraction(sp.cancel(expr))[1]
 return sp.cancel(D/D.subs(z,0))
for h in range(1,151):
 R=rank(h)
 need(denominator(fs(h))==(sp.expand(p**(R+1)) if R>=0 else 1),('F exact denominator',h))
for d in range(2,101,2):
 expr,R=ss(d)
 need(denominator(expr)==(sp.expand(q**(R+1)) if R>=0 else 1),('S exact denominator',d))
orders=[]
for r in range(9):
 h=3**r;Sd,_=ss(2*h);expr=sp.cancel(Sd+(fs(h)-js(h))/2)
 e=2 if r==1 else r//2+1
 need(denominator(expr)==sp.expand(p**(r+1)*q**e),('B exact denominator',r))
 orders.append(3*(r+1)+6*e)
for r in range(101):
 m=r//2;Lr=z**(10*m)*(1+z+z**3) if r%2==0 else z**(10*m+5)
 need(not sp.rem(sp.Poly(Lr,z,modulus=2),sp.Poly(q,z,modulus=2)).is_zero,('top numerator mod2',r))
 allowed=[t for t in range(r+2) if 3**t*(3**t-1)<=2*3**r]
 need(max(allowed)==(1 if r==1 else m),('half bound',r))
# Raw-order coefficients, empty term and increasingly long support gaps.
for n in range(31):
 coeff={3**c:comb(n-2*c,c) for c in range(n//3+1)}
 need(all(v>0 for v in coeff.values()),('raw coefficients positive',n))
 need(max(coeff)==3**(n//3),('raw coefficient exponential degree',n))
need({3**c:comb(0-2*c,c) for c in range(1)}=={1:1},'empty raw aut contributes t')
receipt={'status':'PASS','checks':checks,'max_literal_order':N,'literal_words':word_count,'markers_evaluated':{'scc':2,'nu3_aut':3,'nu5_aut':5,'h':7,'b':11},'formal_scope':'Only u,v,r,s have certified tournament signatures; all remaining letters are synthetic algebra controls, not realized tournaments.','literal_F_J_h_through':200,'literal_S_even_deficit_through':200,'literal_B_degree_through':100,'exact_F_h_through':150,'exact_S_even_deficit_through':100,'exact_B_power3_r_through':8,'orders':orders,'mod2_and_half_bound_r_through':100,'proof_sha256':hashlib.sha256((Path(__file__).resolve().parents[1]/'session23-refined/PROOF.md').read_bytes()).hexdigest(),'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
print(json.dumps(receipt,indent=2))
