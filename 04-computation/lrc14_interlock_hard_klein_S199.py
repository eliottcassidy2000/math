#!/usr/bin/env python3
"""
klein-2026-07-09-S199 (corrected): interlock check on genuinely HARD sets.

Prior test erred: V=spread/0.85 gave spread<6V/7 (NON-hard, j=1 wraparound works).
The interlock is about HARD sets: spread>=6V/7 AND j=1 fails (all gaps<=V/7, i.e.
V/7-dense). A 2-block has a big middle gap => j=1 GOOD => NOT hard => excluded. So
hard sets are evenly-spread (low arc-count). CHECK: for HARD longest-AP=L sets
(L=2..6, k=13), is max c(L)=#arcs/spread < D3_inf^(L)? (route c airtight?)
"""
import numpy as np
from math import ceil
rng=np.random.default_rng(199199)
INV7=1/7
D3L={2:0.86,3:0.85,4:0.84,5:0.81,6:0.76}
def gaps1(E,V):
    p=np.sort(np.mod(np.array(E),V)); g=np.diff(p); return np.append(g,V-p[-1]+p[0])
def is_hard(E,V):  # j=1 fails: all gaps<=V/7
    return gaps1(E,V).max()<=V/7+1e-9
def arcs(E,V,Nx=200000):
    E=np.array(E); x=(np.arange(Nx)+0.5)/Nx*V
    ph=np.sort(np.mod(np.outer(x,E),V),axis=1)/V
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    gi=(g.max(axis=1)>INV7+1e-12).astype(int); ed=np.diff(np.concatenate([gi,gi[:1]]))
    nc=int((ed==1).sum())
    if gi.all():nc=1
    return nc
def longest_AP(E):
    Es=sorted(set(E)); Eset=set(Es); best=1
    for i in range(len(Es)):
        for jj in range(i+1,len(Es)):
            d=Es[jj]-Es[i]; L=2; x=Es[jj]+d
            while x in Eset: L+=1; x+=d
            best=max(best,L)
    return best

k=13
print(f"k={k}: max c(L)=#arcs/spread over HARD (V/7-dense, j=1 fails) longest-AP=L sets, vs D3_inf^L")
print(f"{'L':>3} {'#hard':>6} {'max c(L)':>9} {'D3_inf^L':>9} {'c<D3?':>7} {'margin':>7} {'worst (rel)':>22}")
for L in (2,3,4,5,6):
    maxc=0; worst=None; nh=0
    for _ in range(200000):
        # sample a set, find V making it hard with longest-AP=L
        rest=rng.choice(np.arange(1,140),k-1,replace=False)
        E=tuple(sorted([0]+[int(x) for x in rest]))
        sp=max(E)
        if sp<20: continue
        if longest_AP(E)!=L: continue
        # V in [sp+1, 7sp/6]; hard needs all gaps<=V/7 => V<=7*maxgap-at-j1... find V
        # try V just above sp (spread/V close to 1)
        for V in (sp+1, int(sp*1.05), int(sp*1.10), int(sp*7/6)):
            if V<=sp: continue
            if (max(E)>=6*V/7) and is_hard(E,V):
                nh+=1
                c=arcs(E,V)/sp
                if c>maxc: maxc=c; worst=(E,V)
                break
    d3=D3L[L]
    ws=str(list(worst[0][:7]))+f"...V={worst[1]}" if worst else "-"
    print(f"{L:>3} {nh:>6} {maxc:>9.3f} {d3:>9.3f} {str(maxc<d3):>7} {d3-maxc:>7.3f} {ws:>22}")
print("\n=> HARD L=5,6 sets are V/7-dense (no big gaps => not 2-blocks) => lower arc-count.")
print("   If max c(L)<D3_inf^L for HARD sets at all L<=6, the interlock IS airtight.")
