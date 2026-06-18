# Lean robustness check for M4, M6 forbidden classes at m=5. Fast windows.
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import importlib.util, os

spec = importlib.util.spec_from_file_location(
    "main", os.path.join(os.path.dirname(__file__), "lrc14_tourmap_danger-interval_kps-S2-wf.py"))
mn = importlib.util.module_from_spec(spec); spec.loader.exec_module(mn)

def gcd_list(xs):
    g=0
    for x in xs: g=gcd(g,x)
    return g

def all_classes(m):
    seen={}
    pairs=list(combinations(range(m),2))
    for bits in range(2**len(pairs)):
        adj=[[False]*m for _ in range(m)]
        for idx,(a,b) in enumerate(pairs):
            if (bits>>idx)&1: adj[a][b]=True
            else: adj[b][a]=True
        cn=mn.canon(adj,m)
        seen.setdefault(cn,(mn.ham_paths(adj,m),mn.num_3cycles(adj,m),mn.score_seq(adj,m)))
    return seen

m=5
allc=all_classes(m)
print(f"m={m}: full iso set = {len(allc)}",flush=True)

# M6 at OPTIMUM, growing windows
print("\n[M6 @ OPTIMUM tau0]",flush=True)
for vmax in [12,16,20,24]:
    real={}
    for S in (c for c in combinations(range(1,vmax+1),m) if gcd_list(c)==1):
        S=list(S); gap,tau0=mn.M(S)
        adj=mn.method6_adj(S,tau0); real.setdefault(mn.canon(adj,m),(tau0,tuple(S)))
    forb=[cn for cn in allc if cn not in real]
    print(f"  vmax={vmax}: realized {len(real)}/{len(allc)}, forbidden {len(forb)}",flush=True)
    if vmax==24:
        for cn in sorted(forb,key=lambda c:allc[c]):
            print("     FORBIDDEN",allc[cn],flush=True)

# M6 at ALL candidate times (small window only — tests if it's the optimum or the map itself)
print("\n[M6 @ ALL candidate tau, vmax=16] (does ANY tau realize the forbidden classes?)",flush=True)
real={}
for S in (c for c in combinations(range(1,17),m) if gcd_list(c)==1):
    S=list(S)
    for tau0 in mn.cand(S):
        adj=mn.method6_adj(S,tau0); real.setdefault(mn.canon(adj,m),tau0)
forb=[cn for cn in allc if cn not in real]
print(f"  realized {len(real)}/{len(allc)}, forbidden {len(forb)}",flush=True)
for cn in sorted(forb,key=lambda c:allc[c]):
    print("     STILL FORBIDDEN over all tau:",allc[cn],flush=True)

# M4 growing windows
print("\n[M4 section rotational @ grid]",flush=True)
for vmax in [13,18,24,30,40]:
    real={}
    for S in (c for c in combinations(range(1,vmax+1),m) if gcd_list(c)==1):
        S=list(S)
        for a in range(1,14):
            if gcd(a,14)!=1: continue
            adj=mn.method4_adj(S,a); real.setdefault(mn.canon(adj,m),(a,tuple(S)))
    forb=[cn for cn in allc if cn not in real]
    print(f"  vmax={vmax}: realized {len(real)}/{len(allc)}, forbidden {len(forb)}",flush=True)
    if vmax==40:
        for cn in sorted(forb,key=lambda c:allc[c]):
            print("     FORBIDDEN",allc[cn],flush=True)
