"""Independent exact affine subgroup and family-rate checks; no producer imports."""
from fractions import Fraction as F
from itertools import product
from math import comb
from collections import Counter
import json
INF=float('inf'); counts=Counter()
def gate(ok,label):
 counts[label]+=1
 if not ok:raise RuntimeError(label)
def valuation(n):
 if not n:return INF
 r=0
 while n%5==0:r+=1;n//=5
 return r
def coeff(w):
 d=len(w);clock=0;B=0
 for i,a in enumerate(w):B+=3**(d-1-i)*2**clock;clock+=a
 return 3**d,2**clock,B
def pair(w,mod):
 P,Q,B=coeff(w);z=pow(Q,-1,mod);return P*z%mod,B*z%mod
def compose(A,B,mod):
 a,b=A;c,d=B;return a*c%mod,(a*d+b)%mod
def subgroup(A,mod):
 S=set();x=(1,0)
 while x not in S:S.add(x);x=compose(A,x,mod)
 if x!=(1,0):raise RuntimeError('not invertible')
 return S
def words(n):
 yield ()
 for a in range(1,n+1):
  for v in words(n-a):yield (a,)+v
def order_params(w):
 P,Q,B=coeff(w);a=P*pow(Q,-1,5)%5
 t=next(i for i in (1,2,4) if pow(a,i,5)==1)
 r=min(valuation(P-Q),valuation(B)) if t==1 else valuation(P**t-Q**t)
 return t,r
def one_phase_depth(w,v):
 P,Q,B=coeff(w);R,S,C=coeff(v);t,r=order_params(w)
 comm=valuation((P-Q)*C-B*(R-S))
 if t==1:
  # Used only for constant-letter controllers; translation is a unit.
  gate(valuation(B)==0,'unit_translation')
  return comm
 if pow(R*pow(S,-1,5)%5,t,5)!=1:return 0
 beta=valuation(R**t-S**t)
 return min(comm,beta if beta<r else INF)
def geom(t,r,D):
 if D==INF:return F(4*r+1,4*t)
 D=int(D)
 return F(D,t) if D<=r else F(r,t)+F(1,4*t)*(1-F(1,5**(D-r)))
def exact_rate(w,u):
 t,r=order_params(w)
 return sum((geom(t,r,one_phase_depth(w,w[:j]+u)) for j in range(1,len(w)+1)),F(0))
def reset_max(L):
 r=valuation(9**(L+1)+16*4**L);z=F(4*r+1,16);h=2;e=1
 while h<=L:
  z+=F(L//h,4*5**max(0,e-r));e+=1;h*=5
 return z
def digits_le(a,b):
 while a or b:
  if a%5>b%5:return False
  a//=5;b//=5
 return True
def entry_predict(L,u):
 w=(1,)*L+(2,)
 for j in range(L+1):
  v=(1,)*j+u
  if len(v)%len(w)==0 and v==w*(len(v)//len(w)):
   return j%2<=L%2 and digits_le(j//2,L//2)
 return False
entries=list(words(6))
for L in range(21):
 w=(1,)*L+(2,);top=reset_max(L)
 gate(exact_rate(w,())==top,'reset_empty_mean')
 gate(top<F(L+2,2),'reset_clock_bound')
 for u in entries:
  R=exact_rate(w,u)
  gate(R<=top,'reset_entry_cap')
  gate((R==top)==entry_predict(L,u),'reset_entry_equality')
 for e in range(1,4):
  mod=5**e;S=subgroup(pair(w,mod),mod);A=subgroup(pair((1,),mod),mod)
  gate(S&A=={(1,0)},'literal_trivial_intersection')
  r=valuation(9**(L+1)+16*4**L)
  gate(len(S)==4*5**max(0,e-r),'literal_reset_order')
  h=2*5**(e-1)
  for u in entries:
   actual=sum(pair((1,)*j+u,mod) in S for j in range(L+1))
   expected=sum(one_phase_depth(w,(1,)*j+u)>=e for j in range(L+1))
   gate(actual==expected,'literal_reset_memberships')
   gate(actual<=1+L//h,'literal_layer_cap')
for L in range(201):
 K=L//2;delta=L%2;good=[]
 for j in range(L+1):
  h=2;ok=True
  while h<=2*(L+1):
   if j%h>L%h:ok=False
   h*=5
  gate(ok==(j%2<=delta and comb(K,j//2)%5!=0),'lucas_residue_equivalence')
  if ok:good.append(j)
 d=K;n=delta+1
 while d:n*=d%5+1;d//=5
 gate(len(good)==n,'lucas_phase_count')
for a in range(1,101):
 w=(a,);t,r=order_params(w);top=F(4*r+1,4*t)
 gate(top/F(a)<=F(5,8),'constant_clock_bound')
 gate((top/F(a)==F(5,8))==(a==1),'constant_clock_equality')
 for u in entries:
  R=exact_rate(w,u);power=all(x==a for x in u)
  gate(R<=top,'constant_entry_cap')
  gate((R==top)==power,'constant_entry_equality')
 if a<=20:
  for e in range(1,4):
   mod=5**e;S=subgroup(pair(w,mod),mod)
   gate(len(S)==t*5**max(0,e-r),'literal_constant_order')
   for u in entries:gate((pair(w+u,mod) in S)==(one_phase_depth(w,w+u)>=e),'literal_constant_membership')
# Complete finite exponent residue classes: all classes modulo4*5^3.
E=3
for residue in (0,1,2):
 rs=[]
 for k in range(5**E):
  a=residue+4*k
  if a==0:a+=4*5**E
  rs.append(order_params((a,))[1])
 for depth in range(1,E+2):gate(sum(r>=depth for r in rs)==5**(E-depth+1),'parameter_haar_layers')
print(json.dumps({'status':'PASS','total_gates':sum(counts.values()),'counts':dict(sorted(counts.items())),'universes':{'reset_L':'0..20','entries_clock_at_most':6,'lucas_L':'0..200','constant_a':'1..100','literal_depth':'1..3','finite_parameter_precision':3}},indent=2))
