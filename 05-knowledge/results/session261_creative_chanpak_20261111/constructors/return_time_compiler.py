#!/usr/bin/env python3
"""Effective all-length zero test for a fixed transitive-port return equation."""
from module_tools import *
from math import comb,factorial
import json
HERE=Path(__file__).resolve().parent

def run_coefficients(T,O,v):
 D=len(T)+1;C=[]
 for r in range(D+1):
  Z,P=m.inflate(T,O,v,r);rows=list(Z);start=len(T)-1;mask=((1<<r)-1)<<start
  for a in range(start,len(rows)):rows[a]&=~mask
  h,w,u=m.data(tuple(rows),P)
  R=tuple(sum(1<<i for i in range(len(rows)) if rows[i]>>j&1) for j in range(len(rows)))
  vv=m.data(R,P)[2];values=(h,h-w,u,vv)
  if any(x%factorial(r) for x in values):raise RuntimeError('skeleton divisibility')
  C.append(tuple(x//factorial(r) for x in values))
 return C

def power_coeffs(cs):
 out={}
 for r,c in enumerate(cs):
  for j in range(r+1):out[j]=out.get(j,0)+c*(-1)**(r-j)*comb(r,j)
 return {j:c for j,c in out.items() if c and j>0}

def all_positive_zeros(coeff):
 coeff={j:c for j,c in coeff.items() if c and j>0}
 if not coeff:return {'all_positive_integers':True}
 d=max(coeff);a=abs(coeff[d]);s=sum(abs(c) for j,c in coeff.items() if j<d)
 M=1
 while a*d**M<=s*(d-1)**M:M+=1
 roots=[n for n in range(1,M) if sum(c*j**n for j,c in coeff.items())==0]
 return {'all_positive_integers':False,'dominant_base':d,'first_uniform_nonzero_index':M,'all_positive_zeros':roots}

def main():
 T=(2,4,9,19,7);O=18;v=0;C=run_coefficients(T,O,v)
 eq=power_coeffs([r[2]-r[3] for r in C]);ans=all_positive_zeros(eq)
 if eq!={1:2,2:-7,3:4} or ans['all_positive_zeros']!=[1]:raise RuntimeError('balanced return classification')
 control=all_positive_zeros({2:0})
 if not control['all_positive_integers']:raise RuntimeError('zero identity control')
 out={'status':'PASS','rows':T,'selected_mask':O,'port':v,'run_coefficients_H_D_U_V':C,'U_minus_V_power_coefficients':eq,'exact_return_classification':ans,'scope':'Effective complete zero classification for each fixed finite integer positive-base exponential polynomial, not arbitrary changing-core return programs.'}
 HERE.joinpath('return_time_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
