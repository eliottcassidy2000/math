#!/usr/bin/env python3
"""All-r9-grid constructors for levels1,...,6; no assertion about s>=7."""
from realize9_lower import tail,realize as realize_lower
from module_tools import *
import argparse,json
HERE=Path(__file__).resolve().parent
SOURCES={
 2:((2,12,9,1),4,9),
 3:((2,28,25,17,1),1,9),
 4:((2,4,57,3,11,27),1,18),
 5:((2,12,25,49,35,7),1,36),
 6:((6,12,24,33,107,7,47),64,288),
}

def realize(r,s):
 if r<1 or s not in range(1,7):raise ValueError('require r>=1 and 1<=s<=6')
 if s==1:return realize_lower(r,1)
 B,O,U=SOURCES[s];j=(U&-U).bit_length()-1;t=r+s-j
 if t<1:raise RuntimeError('missing finite head')
 return tail(B,O,t)

def main():
 p=argparse.ArgumentParser();p.add_argument('--r',type=int);p.add_argument('--s',type=int);a=p.parse_args()
 if a.r is not None:
  if a.s is None:p.error('--s required with --r')
  T=realize(a.r,a.s);print(json.dumps({'r':a.r,'s':a.s,'rows':T}));return
 seeds=[]
 for s,(B,O,U) in SOURCES.items():
  h,w,u=m.data(B,O)
  if (h-w,u)!=(1+2**s,U):raise RuntimeError('offset source mismatch')
  seeds.append({'s':s,'rows':B,'O':O,'H':h,'W':w,'D':h-w,'U':u,'j':(u&-u).bit_length()-1})
 controls=[]
 for r in (1,2,3):
  for s in range(1,7):
   T=realize(r,s);h=sum(m.endpoints(T)[-1]);target=1+2**s+9*2**(r+s)
   if h!=target:raise RuntimeError(('grid mismatch',r,s,h,target))
   controls.append({'r':r,'s':s,'order':len(T),'H':h})
 fan=[]
 for p0 in range(4):
  for q in range(4):
   n=p0+q+1;T,_=m.inflate((2,4,1),0,0,n);O=1<<(2+p0);h,w,u=m.data(T,O)
   expected=2**(p0+1)*3**q+2**(q+1)-(1 if p0==0 else 0)
   if (h,w,u)!=(2**n+1,0,expected):raise RuntimeError('fan formula')
   if m.data(T,1)[2]!=2**n+1 or m.data(T,2)[2]!=2**(n+1)-1:raise RuntimeError('fan exterior formulas')
   fan.append({'p':p0,'q':q,'H':h,'U':u})
 out={'status':'PASS','scope':'All r>=1 at s1,...,6 by offset-tail proof; eighteen literal controls and sixteen fan controls; s>=7 remains open in this packet','sources':seeds,'literal_grid_controls':controls,'fan_controls':fan}
 HERE.joinpath('realize9_six_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
