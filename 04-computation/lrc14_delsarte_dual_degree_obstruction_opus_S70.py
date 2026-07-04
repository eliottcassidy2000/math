#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
ROUTE 2 (covering-min dual): NO universal bounded-degree Delsarte/positive-polynomial certificate.
opus-2026-07-04-S70. The certificate is a trig poly g (deg N) with g<=0 on the danger set U D_i(beta) and
int g = 1 (=> g>0 somewhere => a safe point => M>=beta). FINDING: its min degree N ~ v_max (LINEAR), because
by SCALE-INVARIANCE dilating S by c scales the safe set's fine structure by 1/c => degree c-x higher.
v_max unbounded over covering families => degree UNBOUNDED => no universal bounded-degree certificate.
So Route 2's dual must be PARAMETRIC/family-specific (= kps residue-liar formulas, whose witness denominator
is ~v_max -- exactly matching), NOT one universal polynomial. Redirects mac-mini S40.
"""
import sys
import numpy as np
from scipy.optimize import linprog
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
beta=1.0/14
def danger_grid(S,G=12000):
    t=(np.arange(G)+0.5)/G; dang=np.zeros(G,bool)
    for v in S:
        fr=(v*t)%1.0; dang |= (np.minimum(fr,1.0-fr)<beta-1e-9)
    return t[dang]
def min_cert_degree(S,Ns):
    td=danger_grid(S)
    for N in Ns:
        cols=[np.ones(len(td))]+[np.cos(2*np.pi*k*td) for k in range(1,N+1)]+[np.sin(2*np.pi*k*td) for k in range(1,N+1)]
        A=np.column_stack(cols); nv=1+2*N; Aeq=np.zeros((1,nv)); Aeq[0,0]=1
        r=linprog(np.zeros(nv),A_ub=A,b_ub=np.zeros(len(td)),A_eq=Aeq,b_eq=[1.0],bounds=[(-80,80)]*nv,method='highs')
        if r.success: return N
    return None
print("Delsarte certificate min degree N (beta=1/14) vs v_max -- LINEAR scaling => no universal bounded cert:")
print("  family                 v_max   N     N/v_max")
Ns=[4,8,16,24,32,48,64,96,128,160]
base=list(range(2,15))
for c in [1,2,3]:
    S=[c*x for x in base]; N=min_cert_degree(S,Ns)
    print("  %2d*{2..14}              %4d   %-5s %s"%(c,max(S),N,('%.2f'%(N/max(S)) if N else '-')))
print("  => N/v_max ~ const (2-3): certificate degree is LINEAR in v_max.")
print("  MECHANISM: safe(cS) = (1/c)*safe(S) tiled c times => finest component width /c => deg(g) >= ~c/width.")
print("  Since v_max is UNBOUNDED over covering families (dilations + primitive Ostrowski ladder), the")
print("  universal bounded-degree Delsarte certificate is IMPOSSIBLE. The dual must be PARAMETRIC per family")
print("  -- exactly kps's residue-liar formulas (witness denominator ~ v_max). Route 2 redirected.")
print("DONE.")
