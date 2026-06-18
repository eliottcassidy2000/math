# Detailed forbidden-class analysis for the promising methods M3, M4, M6.
# Identify EXACTLY which iso classes are forbidden, and stress-test with larger families.

from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
import importlib.util, sys, os

# load the main module to reuse functions
spec = importlib.util.spec_from_file_location(
    "main", os.path.join(os.path.dirname(__file__), "lrc14_tourmap_danger-interval_kps-S2-wf.py"))
mn = importlib.util.module_from_spec(spec); spec.loader.exec_module(mn)

F_ = F
def gcd_list(xs):
    g=0
    for x in xs: g=gcd(g,x)
    return g

# ---- enumerate ALL iso classes on m vertices with canonical keys + invariants ----
def all_classes(m):
    seen = {}
    # iterate over all tournaments: for each unordered pair choose orientation
    pairs = list(combinations(range(m),2))
    for bits in range(2**len(pairs)):
        adj=[[False]*m for _ in range(m)]
        for idx,(a,b) in enumerate(pairs):
            if (bits>>idx)&1:
                adj[a][b]=True
            else:
                adj[b][a]=True
        cn = mn.canon(adj,m)
        if cn not in seen:
            h=mn.ham_paths(adj,m); c3=mn.num_3cycles(adj,m); sc=mn.score_seq(adj,m)
            seen[cn]=(h,c3,sc)
    return seen

def realized_set(adj_fn, m, vmax, uses_grid):
    sets=[c for c in combinations(range(1,vmax+1),m) if gcd_list(c)==1]
    real={}
    for S in sets:
        S=list(S)
        if uses_grid:
            for a in range(1,14):
                if gcd(a,14)!=1: continue
                adj=adj_fn(S,a); cn=mn.canon(adj,m)
                if cn not in real:
                    real[cn]=(mn.ham_paths(adj,m),mn.num_3cycles(adj,m),mn.score_seq(adj,m),(tuple(S),a))
        else:
            gap,tau0=mn.M(S)
            if tau0 is None: continue
            adj=adj_fn(S,tau0); cn=mn.canon(adj,m)
            if cn not in real:
                real[cn]=(mn.ham_paths(adj,m),mn.num_3cycles(adj,m),mn.score_seq(adj,m),(tuple(S),tau0))
    return real

methods=[("M3",mn.method3_adj,False),("M4",mn.method4_adj,True),("M6",mn.method6_adj,False)]
vmax_map={3:16,4:14,5:12}  # larger window than the first pass

for m in [4,5]:
    print("#"*72)
    print(f"VERTEX COUNT m={m}   (full iso set: {mn.FREE[m]})")
    print("#"*72)
    allc = all_classes(m)
    print(f"  enumerated {len(allc)} iso classes")
    for name,fn,grid in methods:
        real = realized_set(fn,m,vmax_map[m],grid)
        forb = [cn for cn in allc if cn not in real]
        print(f"\n  [{name}] realized {len(real)}/{len(allc)}  (vmax={vmax_map[m]})")
        print(f"    REALIZED invariants (H,c3,score):")
        for cn,(h,c3,sc,ex) in sorted(real.items(), key=lambda kv: (kv[1][0],kv[1][1])):
            print(f"       H={h} c3={c3} score={sc}   ex S={ex[0]} key={ex[1]}")
        if forb:
            print(f"    FORBIDDEN invariants (H,c3,score):")
            for cn in sorted(forb, key=lambda c:(allc[c][0],allc[c][1])):
                h,c3,sc=allc[cn]
                print(f"       H={h} c3={c3} score={sc}")
