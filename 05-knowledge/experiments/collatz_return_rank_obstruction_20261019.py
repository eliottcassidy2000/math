#!/usr/bin/env python3
"""Exact controls for dual local ranks and two-word rank recharge."""
from fractions import Fraction as Q
import json
CHECKS=0
def check(test,label):
 global CHECKS
 CHECKS+=1
 if not test:raise RuntimeError(label)
def vp(x,p):
 x=Q(x)
 if not x:return None
 n=abs(x.numerator);d=x.denominator;v=0
 while n%p==0:n//=p;v+=1
 while d%p==0:d//=p;v-=1
 return v

def T(x):
 x=Q(x)
 check(x.denominator%2==1,'2-adically integral source')
 return (3*x+1)/2 if x.numerator%2 else x/2

def path(x,n):
 x=Q(x);a=C=0;bits=[]
 for j in range(n):
  bit=x.numerator%2;bits.append(bit)
  if bit:a+=1;C=3*C+2**j
  x=T(x)
 return x,a,C,bits

def normalized_step(k,E,beta):
 if E%2==0:
  if beta==0:return k,E//2
  return k,(3*E+3**max(0,-k)-3**max(0,k))//2
 if beta==0:
  return k+1,(3*E+1)//2 if k>=0 else (E+3**(-k-1))//2
 return k-1,(E-3**(k-1))//2 if k>=1 else (3*E-1)//2

# Two actual oriented adjacent-pair return words.
F={'name':'F','n':2,'a':1,'C':1,'residue':1}
H={'name':'H','n':9,'a':4,'C':101,'residue':43}
for W in (F,H):
 n,a,C=W['n'],W['a'],W['C']
 W['center']=Q(C,2**n-3**a)
 _,aa,CC,bits=path(W['residue'],n)
 check((aa,CC)==(a,C),'actual return coefficients')
 W['bits']=bits
 for j in range(8):
  x=W['residue']+j*2**n
  v,b,Cv,bv=path(x,n);u,c,Cu,bu=path(x+1,n)
  check(u==v+1 and b==c==a and Cv==C,'lawful oriented adjacent return')
  check(v==Q(3**a*x+C,2**n),'actual affine return')
  k=0;E=1;U=Q(x+1);V=Q(x)
  for _ in range(n):
   k,E=normalized_step(k,E,V.numerator%2);U,V=T(U),T(V)
   check(E==3**max(0,-k)*(U-Q(3)**k*V),'integer normalized transition')
  check((k,E)==(0,1),'normalized state returns exactly')

# Exact forward and inverse valuation budgets for F and H.
budget_controls=[]
for W in (F,H):
 n,a,C,p=W['n'],W['a'],W['C'],W['center']
 for y in (1,7,9,19,43,57,115,127,169):
  rank=vp(Q(y)-p,3)
  if rank is None:
   check(Q(2**n*y-C,3**a)==y,'fixed inverse source')
   continue
  for j in range(rank//a+2):
   value=p+Q(2**n,3**a)**j*(Q(y)-p)
   check((value.denominator==1)==(a*j<=rank),'inverse integrality iff valuation budget')
  budget_controls.append({'word':W['name'],'seed':y,'three_adic_rank':rank,
                          'maximum_integer_inverse_blocks':rank//a})
 for x in (1,3,7,9,19,43,57,115,127,169):
  rank=vp(Q(x)-p,2)
  if rank is None:continue
  for j in range(rank//n+2):
   _,_,_,bits=path(x,n*j)
   check((bits==W['bits']*j)==(n*j<=rank),'forward legality iff valuation budget')

# Fixed-seed examples transport isolated atoms but exhaust the 3-adic rank.
seed_examples=[]
for W,y in ((F,7),(H,7),(F,43)):
 n,a,C,p=W['n'],W['a'],W['C'],W['center']
 z=Q(2**n*y-C,3**a)
 check(z.denominator==1 and z>0,'actual integer inverse source')
 lo,_,_,_=path(z,n);hi,_,_,_=path(z+1,n)
 check(lo==y and hi==y+1,'actual pair transport')
 u=Q(z+1);v=z;first=None;aa=bb=0
 for t in range(1,40):
  aa+=u.numerator%2;bb+=v.numerator%2;u,v=T(u),T(v)
  if u==v:first=t;break
 check(first is not None and aa!=bb,'transported isolated first collision')
 seed_examples.append({'word':W['name'],'seed':y,'lifted_lower_source':int(z),
                       'first_isolated_merge':first,'odd_difference':aa-bb})

# Arbitrarily large backward recharge in a fixed two-word bank.
backward_recharge=[]
p=H['center']
check(vp(p-1,3)==1,'centers have exact 3-adic separation one')
for M in (4,5,8,12,16,24,32,40):
 q=3**(M+1)
 r=(101*pow(431,-1,q)+3**M)%q
 z=r+q*((1-r)*pow(q,-1,4)%4)
 y=(3*z+1)//4
 check(z>1 and z%4==1 and 4*y==3*z+1,'positive lawful F step')
 check(vp(Q(z)-p,3)==M,'arbitrarily specified H rank')
 check(vp(y-1,3)==2 and vp(z-1,3)==1,'constant F ranks before and after switch')
 current=Q(z)
 for j in range(M//4):
  nxt=Q(512*current-101,81)
  check(nxt.denominator==1 and nxt>current>p,'positive integer H inverse lift')
  lo,_,_,_=path(nxt,9);hi,_,_,_=path(nxt+1,9)
  check(lo==current and hi==current+1,'lawful H pair lift')
  current=nxt
 check(Q(512*current-101,81).denominator!=1,'backward budget is exact')
 backward_recharge.append({'M':M,'starting_seed':y,'F_inverse_source':z,
                            'F_rank_before':2,'F_rank_after':1,
                            'H_rank_after_switch':M,'lawful_H_inverse_blocks':M//4,
                            'final_positive_integer_source':int(current)})

# Dual forward recharge: a fixed F budget can feed arbitrarily many H blocks.
forward_recharge=[]
check(vp(p-1,2)==1,'centers have exact 2-adic separation one')
for M in (9,10,18,27,36,45,54):
 q=2**(M+1)
 r=(101*pow(431,-1,q)+2**M)%q
 z=r+q*((1-r)*pow(q,-1,3)%3)
 x=(4*z-1)//3
 check(x>0 and 3*x==4*z-1,'positive integral F inverse for forward witness')
 check(vp(x-1,2)==3 and vp(z-1,2)==1,'constant F forward ranks')
 check(vp(Q(z)-p,2)==M,'arbitrary H forward rank')
 v,_,_,_=path(x,2)
 check(v==z,'actual F forward prefix')
 v,_,_,bits=path(z,9*(M//9))
 check(bits==H['bits']*(M//9) and v>0,'long legal positive H block sequence')
 _,_,_,bits=path(z,9*(M//9+1))
 check(bits!=H['bits']*(M//9+1),'forward budget is exact')
 forward_recharge.append({'M':M,'initial_source':x,'after_F':z,
                          'F_rank_before':3,'F_rank_after':1,
                          'H_rank_after_switch':M,'lawful_H_forward_blocks':M//9})

# Bounded normalized cycles distinguish finite positive prefixes, actual
# positive phases, and negative source centers.
negative_center_prefixes=[]
for L in range(1,21):
 v=Q(2**(L+1)-2);u=v+1;k=0;E=1
 for t in range(L+1):
  check(v>0 and u>0,'finite source prefix is positive')
  beta=v.numerator%2
  check(beta==(0 if t==0 else 1),'prefix zero then all ones')
  k,E=normalized_step(k,E,beta);u,v=T(u),T(v)
  check((k,E)==(1,2),'bounded nonidentity all-odd loop')
 negative_center_prefixes.append({'L':L,'positive_lower_source':2**(L+1)-2,
                                  'loop_state':[k,E]})
u=Q(2);v=Q(1);k=0;E=1
phase_states=[]
for t in range(4):
 phase_states.append([k,E]);check(u!=v,'positive phase loop never merges')
 k,E=normalized_step(k,E,v.numerator%2);u,v=T(u),T(v)
check(phase_states==[[0,1],[-1,1],[0,1],[-1,1]],'actual positive phase two-cycle')
check(Q(1,2-3)==-1,'all-odd loop source center is negative')

print(json.dumps({'status':'FINITE-EXACT controls for proved local-rank obstruction; no universal coverage',
                  'checks':CHECKS,
                  'return_words':[{**W,'center':str(W['center'])} for W in (F,H)],
                  'fixed_seed_budget_controls':budget_controls,
                  'transported_isolated_atoms':seed_examples,
                  'backward_rank_recharge':backward_recharge,
                  'forward_rank_recharge':forward_recharge,
                  'negative_center_positive_prefixes':negative_center_prefixes,
                  'positive_phase_cycle':phase_states},indent=2,sort_keys=True))
