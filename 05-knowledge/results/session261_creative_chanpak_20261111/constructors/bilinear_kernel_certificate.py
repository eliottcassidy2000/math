#!/usr/bin/env python3
"""Literal quotient and run coefficients for the two-port 329 norm kernel."""
from module_tools import *
from math import factorial
import json
HERE=Path(__file__).resolve().parent
Q=(198,196,136,211,167,207,20,64)

def substitute(A,B):
 blocks=[(0,)]*8;blocks[5]=A;blocks[7]=B;starts=[];n=0
 for Z in blocks:starts.append(n);n+=len(Z)
 out=[0]*n
 for i,Z in enumerate(blocks):
  for a,row in enumerate(Z):
   out[starts[i]+a]|=row<<starts[i]
   for j,Y in enumerate(blocks):
    if Q[i]>>j&1:out[starts[i]+a]|=((1<<len(Y))-1)<<starts[j]
 return tuple(out)

def profile(A):
 e=m.endpoints(A);h=sum(e[-1]);full=(1<<len(A))-1
 c2=sum(sum(e[S])*sum(e[full^S]) for S in range(1,full))
 return h,c2

def main():
 cs=[]
 for r in range(1,4):
  row=[]
  for s in range(1,4):
   value=sum(m.endpoints(substitute((0,)*r,(0,)*s))[-1]);den=factorial(r)*factorial(s)
   if value%den:raise RuntimeError('joint skeleton divisibility')
   row.append(value//den)
  cs.append(row)
 if cs!=[[149,72,0],[72,37,0],[0,0,0]]:raise RuntimeError('joint skeleton mismatch')
 mods=[(0,),(2,0),(6,4,0),(2,4,1)];controls=[]
 for A in mods:
  for B in mods:
   h,c=profile(A);g,d=profile(B);actual=sum(m.endpoints(substitute(A,B))[-1]);expected=149*h*g+72*(c*g+h*d)+37*c*d
   if actual!=expected:raise RuntimeError('arbitrary-module control')
   x=37*(c+2*h)-2*h;y=37*(d+2*g)-2*g
   if 37*actual!=x*y+329*h*g:raise RuntimeError('norm identity')
   controls.append({'A':A,'B':B,'H':actual})
 small=[]
 for A in mods:
  h,c=profile(A);offset=9*h-2*(c+2*h);small.append({'rows':A,'H':h,'Q2':c+2*h,'offset':offset})
 if [x['rows'] for x in small if x['offset']==1]!=[(2,0)]:raise RuntimeError('unit row classification')
 out={'status':'PASS','quotient_rows':Q,'ports':[5,7],'joint_run_coefficients_r_s_1_to_3':cs,'literal_controls':controls,'complete_nonnegative_offset_base':small,'scope':'Higher run coefficients vanish by the proved all-chain formula and positivity, not merely by the r,s<=3 table; the cut bound reduces all nonnegative-offset rows to this finite base.'}
 HERE.joinpath('bilinear_kernel_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
