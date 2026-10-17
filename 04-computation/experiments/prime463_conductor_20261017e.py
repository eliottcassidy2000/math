#!/usr/bin/env python3
"""Exact 463 norm clock, conductor-class order, and typed Paley control."""
import json,subprocess
from math import isqrt
from quadratic_context_20261017c import forms


def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def v2(x):
 need(x!=0,'nonzero valuation input');return (abs(x)&-abs(x)).bit_length()-1

def mul(x,y):
 a,b=x;c,d=y;return a*c-462*b*d,a*d+b*c

def power(x,n):
 r=(1,0)
 while n:
  if n&1:r=mul(r,x)
  x=mul(x,x);n//=2
 return r

def main():
 need(all(463%d for d in range(2,isqrt(463)+1)),'463 is prime')
 base=forms(-1848);need(len(base)==8 and all(b==0 or abs(b)==a or a==c for a,b,c in base),'maximal group C2 cubed')
 rows=[]
 for n in range(1,257):
  a,b=power((1,1),n)
  need(a%2==1 and v2(b)==v2(n),'exact conductor valuation of every power')
  need(a*a+462*b*b==463**n,'norm preserved by multiplication')
  if n&(n-1)==0:rows.append(dict(exponent=n,real=str(a),imaginary_coefficient=str(b),v2_imaginary=v2(b)))
 gp='for(k=0,8,D=-1848*4^k;Q=qfbprimeform(D,463);P=qfbpow(Q,2^k);if(k>0,if(qfbpow(Q,2^(k-1))==P,error("order too small")));print([k,Vec(qfbred(Q)),Vec(P)]))\n'
 p=subprocess.run(['gp','-fq'],input=gp,text=True,capture_output=True,check=True)
 need(not p.stderr.strip(),'PARI error-free exact form composition')
 classes=[json.loads(x) for x in p.stdout.splitlines()]
 need(len(classes)==9 and all(q[2]==[1,0,462*4**q[0]] for q in classes),'463-class has exact order2^k through conductor256')
 # 463 also occurs in the earlier Paley mu14 census. Reproduce its signed
 # first-moment invariant directly; no Collatz transfer is asserted.
 character=lambda a:0 if a%463==0 else (1 if pow(a%463,231,463)==1 else -1)
 forms463=forms(-463);h=len(forms463)
 first=sum(a*character(a) for a in range(1,463));half=sum(character(a) for a in range(1,232))
 need(h==7 and first==-463*h and half==h,'independent class-number / Paley character control')
 return dict(status='FINITE-EXACT controls; all-order conductor theorem in companion proof',base_forms=base,
  powers_checked=256,dyadic_power_rows=rows,pari_class_orders=classes,forms_minus463=forms463,
  paley463=dict(class_number=h,weighted_character_sum=first,half_interval_excess=half),
  warning='D=-463 and D=-1848 are distinct fields; common numeral is not a field identification')

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,separators=(',',':')))
