#!/usr/bin/env python3
"""Exact controls for marked-observer lifts; no orbit-completeness claim."""
from fractions import Fraction
from math import gcd,lcm
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
def v2(n):return (n&-n).bit_length()-1
def step(n,sign):
 z=3*n+sign;a=v2(z)
 return z>>a,a
def order4(mod):
 x=4%mod;d=1
 while x!=1:x=4*x%mod;d+=1
 return d
def prefix(n,sign,r):
 states=[n];word=[]
 for _ in range(r):n,a=step(n,sign);states.append(n);word.append(a)
 return states,word
def setup(states,word,sign,observers):
 r=len(word);h=states[-1];A=sum(word);Ai=[0]
 for a in word:Ai.append(Ai[-1]+a)
 S=sum(3**(r-1-i)*2**Ai[i] for i in range(r))
 mods=[3**(r-i+1)*q//gcd(q,3*h+sign) for i,q in enumerate(observers)]
 M=lcm(*mods);d=order4(M);R=4**d
 K=sign*(R-1)*(2**A+3*S)//3**(r+1)
 return Ai,S,mods,M,d,R,K

def main():
 checks=0;minimality=0;largest_d=0;examples=[]
 for sign in [-1,1]:
  for start in [3,5,17]:
   for r in range(4):
    states,word=prefix(start,sign,r)
    for observers in [[1]*(r+1),[5]*(r+1),[[3,7,13][i%3] for i in range(r+1)]]:
     Ai,S,mods,M,d,R,K=setup(states,word,sign,observers)
     largest_d=max(largest_d,d);n=start
     for k in range(3):
      actual,actual_word=prefix(n,sign,r)
      assert actual_word==word
      assert all((a-b)%q==0 for a,b,q in zip(actual,states,observers))
      join,a_join=step(actual[-1],sign);old_join,old_a=step(states[-1],sign)
      assert join==old_join and a_join==old_a+2*k*d
      expected_h=4**(k*d)*states[-1]+sign*(4**(k*d)-1)//3
      assert actual[-1]==expected_h
      n=R*n+K;checks+=1
     # Independent rational backwards reconstruction tests the exact least period.
     for j in range(1,d+1):
      h=4**j*states[-1]+sign*(4**j-1)//3;back=[Fraction(h)]
      for a in reversed(word):back.append((2**a*back[-1]-sign)/3)
      back=back[::-1]
      legal=all(x.denominator==1 and x.numerator>0 and x.numerator%2==1 for x in back)
      compatible=legal and all((x.numerator-old)%q==0 for x,old,q in zip(back,states,observers))
      assert compatible==(j%d==0),(sign,states,observers,d,j)
      minimality+=1
 for sign,start in [(1,3),(-1,5)]:
  states,word=prefix(start,sign,1);observers=[4*n+sign for n in states]
  Ai,S,mods,M,d,R,K=setup(states,word,sign,observers)
  actual,actual_word=prefix(R*start+K,sign,1)
  tables=[4*n+sign for n in actual]
  assert actual_word==word
  assert all(new>old and new%old==0 for new,old in zip(tables,observers))
  examples.append({'sign':sign,'original_states':states,'word':word,'observers':observers,
                   'individual_clock_moduli':mods,'least_period':d,'R':R,'K':K,
                   'first_lifted_states':actual,'first_lifted_tables':tables,
                   'factor_cofactors':[new//old for new,old in zip(tables,observers)]})
 factor_checks=0
 for sign in [-1,1]:
  for n in range(1,32,2):
   m=4*n+sign
   for q in range(3,m+1,2):
    if m%q:continue
    d=order4(3*q)
    for j in range(2*d+1):
     mj=4**j*m+sign*(4**j-1)//3
     assert (mj%q==0)==(j%d==0)
     factor_checks+=1
 receipt={'status':'FINITE-EXACT controls of arbitrary-parameter theorems',
          'prefix_universe':'Both signs; starts 3,5,17; r=0,1,2,3; observer words all1, all5, or cycle3/7/13; k=0,1,2',
          'lifted_prefix_checks':checks,'least_period_rational_checks':minimality,'largest_period':largest_d,
          'factor_clock_universe':'Both signs; odd n<32; every odd divisor q>=3 of 4n+sign; 0<=j<=2ord_(3q)(4)',
          'factor_clock_checks':factor_checks,'composite_examples':examples,
          'hostile_controls':'Every smaller positive increment than the claimed least period is rejected by legal-head or observer congruence reconstruction. Minus examples retain their original cycle basin.',
          'universal_convergence_claim':False}
 (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
