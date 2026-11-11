#!/usr/bin/env python3
"""New finite controls for run-order filtration and closed ordinal response jets."""
from pathlib import Path
from functools import lru_cache
from math import comb
import importlib.util,json,sys
sys.dont_write_bytecode=True
import sympy as s
BASE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('kernel',BASE.parent/'core/response_kernel.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

def ordinal(A,B):
 n=len(A);m=len(B)
 return tuple([x|(((1<<m)-1)<<n) for x in A]+[x<<n for x in B])
def qjet(C,d):return [sum(comb(r,j)*C.get(j,0) for j in range(1,r+1)) for r in range(1,d+1)]
def triangle_module(m):
 # M -> a -> b -> M
 M=mod.tt(m)
 return tuple(x|(1<<m) for x in M)+(1<<(m+1),(1<<m)-1)
def clone(T,v,m):
 groups=[mod.tt(m) if i==v else (0,) for i in range(len(T))]
 starts=[];n=0
 for G in groups:starts.append(n);n+=len(G)
 R=[]
 for i,G in enumerate(groups):
  out=sum(((1<<len(groups[j]))-1)<<starts[j] for j in range(len(T)) if T[i]>>j&1)
  R.extend(out|(r<<starts[i]) for r in G)
 return tuple(R)

def main():
 # New tests of closure, on tiny independently chosen graph pairs.
 bank=[(0,),mod.tt(2),(2,4,1),mod.tt(3)]
 cases=[]
 for A in bank:
  for B in bank:
   qa,qb=qjet(mod.covers(A),4),qjet(mod.covers(B),4)
   actual=qjet(mod.covers(ordinal(A,B)),4)
   if actual!=[a*b for a,b in zip(qa,qb)]:raise RuntimeError('ordinal jet law')
   cases.append({'orders':[len(A),len(B)],'Q':actual})
 # The new sharp single-module bound p+1 for a two-exterior triangle.
 # Triangular exact solve in j!S(m,j), not fitting beyond an unproved cutoff.
 kernels=[]
 for p in [1,2,3]:
  d=p+1;coeff=[]
  for m in range(1,d+1):
   y=mod.covers(triangle_module(m)).get(p,0)
   F=mod.covers(mod.tt(m))
   residual=y-sum(coeff[j-1]*F.get(j,0) for j in range(1,m))
   den=F[m]
   if residual%den:raise RuntimeError('noninteger skeleton coefficient')
   coeff.append(residual//den)
  if coeff[-1]<=0:raise RuntimeError('sharpness missing')
  # One further physical module and a nontransitive module challenge the new kernel.
  for M in [mod.tt(d+1),(2,4,1)]:
   m=len(M);T=tuple(x|(1<<m) for x in M)+(1<<(m+1),(1<<m)-1)
   C=mod.covers(M)
   if mod.covers(T).get(p,0)!=sum(a*C.get(j,0) for j,a in enumerate(coeff,1)):raise RuntimeError('module kernel control')
  kernels.append({'observer_C':p,'max_run_order':d,'coefficients':coeff})
 h,q,r=s.symbols('h q r');c=q-2*h;f=r-3*q+3*h
 K=s.expand(149*h*h+144*h*c+37*c*c)
 U=s.expand(430*h*h+605*h*c+120*h*f+196*c*c+62*c*f)
 E=s.expand(14*h*h+7*h*c)
 family=[]
 for k in range(1,5):
  sub={h:17**k,q:224**k,r:1449**k}
  vals=[int(z.subs(sub)) for z in [K,U,E]]
  family.append({'k':k,'H_modules':sub[h],'Q2_modules':sub[q],'Q3_modules':sub[r],'K':vals[0],'U0':vals[1],'E0':vals[2],'bare_H':vals[0]+vals[1]-vals[2],'first_tail_H':vals[0]+2*vals[1],'independent_trade_ports':2*k})
 # Rooted operations are NOT included in the congruence. Measure their new guard debt.
 A=(2,4,9,19,39,15);B=(2,4,57,3,11,27)
 ca,cb=mod.covers(clone(A,1,2)),mod.covers(clone(B,0,2))
 hostile={'A_clone_C123':[ca[i] for i in [1,2,3]],'B_clone_C123':[cb[i] for i in [1,2,3]],'normalized_C2_cross_difference':ca[2]*cb[1]-cb[2]*ca[1],'normalized_C3_cross_difference':ca[3]*cb[1]-cb[3]*ca[1]}
 if not hostile['normalized_C2_cross_difference'] and not hostile['normalized_C3_cross_difference']:raise RuntimeError('expected clone guard loss absent')
 receipt={'status':'FINITE-EXACT new controls','ordinal_cases':cases,'sharp_two_exterior_kernels':kernels,'core_family':family,'rooted_clone_guard_hostile':hostile}
 (BASE/'response_filtration_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
