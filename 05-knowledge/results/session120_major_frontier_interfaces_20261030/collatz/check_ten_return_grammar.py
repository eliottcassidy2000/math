#!/usr/bin/env python3
"""Arithmetic bounds and exact modulo-5^12 normalization of E12 candidates."""
from fractions import Fraction as F
import json

def need(b,s):
    if not b: raise RuntimeError(s)
def ceil(x): return -((-x.numerator)//x.denominator)
gamma=F(5878847,19531250)
eta=F(17587861,48828125)
need(1+2*gamma>F(9,8),'two nonone10 blocks excluded')
rows=[]
for j in range(3,10):
    amax=min(ceil(5*(F(j,8)-gamma))-1,(j-3)//2)
    need(amax<5,'no pure fifth return before nonone block')
    rows.append({'j11':j,'nonone10_blocks':1,'max_allone10_blocks':amax,
                 'possible_j10_core':[j+8-2*a for a in range(amax+1)]})
need([r['max_allone10_blocks'] for r in rows]==[0,0,1,1,2,2,3],'grammar bounds')
need(F(2,5)<gamma+F(1,5),'single narrow fee')
need(F(4,5)-eta<gamma+F(1,5),'double narrow fee')
# Direct modular affine composition independently verifies the two-digit
# commutative boundary law, using all625 coefficient pairs against U.
e=10; q=5**e; mod=25*q; T=2*5**9
m=pow(3,T,mod)*pow(pow(2,T,mod),-1,mod)%mod
t=(m-1)%mod
need((m-1)%q==0 and t%q==0,'allones return')
U=((m-1)//q,t//q)
need(tuple(x%5 for x in U)==(4,4),'allones first tangent')
need((4*F(11,8)+1+gamma)/5==1+eta,'core quantized fee')
need(F(1,5)+F(T,2*5**11)==F(6,25),'internal allone fee')
family_caps=[]
for a,b in [(0,1),(1,1),(2,1),(3,1),(0,2),(1,2)]:
 cap=ceil((1-F(a,5)-b*eta)/F(6,25))-1
 family_caps.append({'a12':a,'b11':b,'max_internal_allone10_blocks_total':cap})
need([z['max_internal_allone10_blocks_total'] for z in family_caps]==[2,1,0,0,1,0],'family internal caps')
checks=0
for alpha in range(25):
 for beta in range(25):
  C=(1+q*alpha,q*beta)
  for a in range(5):
   # c followed by u^a, versus u^a followed by c.
   ma=pow(m,a,mod); ta=(ma-1)%mod
   left=(ma*C[0]%mod,(ma*C[1]+ta)%mod)
   right=(C[0]*ma%mod,(C[0]*ta+C[1])%mod)
   predicted=((1+q*((alpha+a*U[0])%25))%mod,q*((beta+a*U[1])%25))
   need(left==right==predicted,'two digit addition')
   checks+=1
print(json.dumps({'status':'PASS','native10_nonone_fee':str(1+gamma),
 'candidate_grammar':rows,'family_internal_allone_caps':family_caps,'fee_per_internal_allone10':'6/25','allones10_clock':T,'allones10_two_digit_tangent_mod25':U,
 'modular_composition_checks':checks,'E12':'OPEN','universal_Collatz':'OPEN'},indent=2))
