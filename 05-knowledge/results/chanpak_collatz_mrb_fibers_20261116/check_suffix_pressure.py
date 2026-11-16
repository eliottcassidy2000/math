"""Integer replay of one-block conditional pressure on the old suffix cone."""
from pathlib import Path
from fractions import Fraction
import json,sys
here=Path(__file__).resolve().parent
for p in ('eight_braid_root','decimal_rank_root','chanpak_eight_braid_20261116','chanpak_decimal_rank_four_20261116'):
 sys.path.insert(0,str(here.parent/p))
from mixed_controller import load,calibrate,need
from check_return_growth import op,h,mul

def upper(x,p,q,D):
 a,b=x.numerator,x.denominator
 C=a**p;R=D**q*b**p
 lo,hi=0,D
 need(hi**q*C>=R)
 while hi-lo>1:
  j=(hi+lo)//2
  if j**q*C>=R:hi=j
  else:lo=j
 need((hi-1)**q*C<R<=hi**q*C)
 return hi

def main():
 d=load();AA,_,_,_,U,_,_=calibrate(d)
 S=(0,3,5,2,4);P=op(S);f=h.reshape((1,-1));g=mul(f,P)[0]
 lam=[]
 for w in d['words']:
  a=mul(mul(f,op(tuple(w)+S)),P)[0]
  lam.append(min(Fraction(int(x),int(y)) for x,y in zip(a,g) if y>0))
 need(len(lam)==168 and min(lam)==Fraction(2439,5) and max(lam)==Fraction(39436,25))
 kinds=list({AA[i]:i for i in U}.values())
 need(len(kinds)==111)
 for name,ids,p,q in [('all',range(168),23,30),('mixed',kinds,25,36)]:
  D=10**12
  N=sum(upper(lam[i],p,q,D) for i in ids)
  need(N<D)
  print(name,len(ids),p,q,N,D)
if __name__=='__main__':main()
