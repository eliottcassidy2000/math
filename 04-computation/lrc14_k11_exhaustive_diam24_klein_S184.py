#!/usr/bin/env python3
"""klein-2026-07-08-S184 -- EXHAUSTIVE compact check of the k=11 covering floor to prim-diam <= 24.
For every PRIMITIVE 11-set with primitive diameter D in [16,24], compute D3 (opus THM-661 covering
floor) and PZ (THM-660) via a batched float grid; report the min per diameter and the global min.
If min D3 >= bar for all D <= 24, the k=11 leg reduces to the prim-diam >= 25 spread tail (discrepancy)."""
import numpy as np
from math import gcd, comb
from functools import reduce
from itertools import combinations
import sys
TH=1/7; M=6/7; BAR=83549/252252
N=1200                      # grid (coarse scan; margins >> grid error ~1e-3)
x=(np.arange(N)+0.5)/N
def batch_floor(shapes):
    """shapes: (B,11) int array -> (D3[B], PZ[B])."""
    B=shapes.shape[0]
    ph=(shapes[:,:,None]*x[None,None,:])%1.0        # B x 11 x N
    ph.sort(axis=1)
    g=np.diff(ph,axis=1); wrap=ph[:,:1,:]+1.0-ph[:,-1:,:]
    gaps=np.concatenate([g,wrap],axis=1)             # B x 11 x N
    W=np.maximum(gaps-TH,0).sum(axis=1)              # B x N
    m1=W.mean(1); m2=(W*W).mean(1); m3=(W**3).mean(1)
    den=m2-m3/M
    D3=np.where(den>0, m1/M+(m1-m2/M)**2/np.where(den>0,den,1), m1/M)
    PZ=m1*m1/m2
    return D3,PZ
print(f"k=11 exhaustive covering floor, prim-diam 16..24; bar={BAR:.5f}; grid N={N}", flush=True)
gmin_d3=(9.9,None); gmin_pz=(9.9,None)
for D in range(16,25):
    mn_d3=(9.9,None); mn_pz=(9.9,None); cnt=0
    buf=[]
    def flush_buf():
        global mn_d3,mn_pz,gmin_d3,gmin_pz,cnt
        if not buf: return
        arr=np.array(buf,dtype=np.float64)
        d3,pz=batch_floor(arr)
        cnt+=len(buf)
        i=int(d3.argmin());
        if d3[i]<mn_d3[0]: mn_d3=(float(d3[i]),buf[i])
        j=int(pz.argmin())
        if pz[j]<mn_pz[0]: mn_pz=(float(pz[j]),buf[j])
        buf.clear()
    for mids in combinations(range(1,D),9):
        E=(0,)+mids+(D,)
        # primitive: gcd of all differences = gcd of the elements (since 0 in set) = gcd(E)
        if reduce(gcd,E)!=1: continue
        buf.append(E)
        if len(buf)>=3000: flush_buf()
    flush_buf()
    if mn_d3[0]<gmin_d3[0]: gmin_d3=mn_d3
    if mn_pz[0]<gmin_pz[0]: gmin_pz=mn_pz
    ok="OK" if mn_d3[0]>=BAR else "**BELOW BAR**"
    print(f"  prim-diam {D}: {cnt:>7} primitive shapes; min D3={mn_d3[0]:.5f} ({'+' if mn_d3[0]>=BAR else ''}{mn_d3[0]-BAR:+.5f}) [{ok}]; min PZ={mn_pz[0]:.5f}; D3-minimizer={mn_d3[1]}", flush=True)
print(f"\nGLOBAL min D3 over prim-diam 16..24 = {gmin_d3[0]:.5f} at {gmin_d3[1]} (bar {BAR:.5f}, margin {gmin_d3[0]-BAR:+.5f})")
print(f"GLOBAL min PZ = {gmin_pz[0]:.5f} at {gmin_pz[1]}")
print("VERDICT: "+("ALL prim-diam<=24 CLEAR -- k=11 leg reduces to prim-diam>=25 spread tail (discrepancy)" if gmin_d3[0]>=BAR else "SOME DIP BELOW -- investigate"))
