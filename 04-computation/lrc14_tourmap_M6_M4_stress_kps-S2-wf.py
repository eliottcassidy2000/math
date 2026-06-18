# Stress test M4 and M6 forbidden classes with LARGE speed windows + clean analysis.
# Key question: is the forbidden-ness robust (structural) or a small-window artifact?
#
# M6 (difference-winding): vertex=runner; at optimum tau0, arc i->j iff
#   frac((v_i - v_j) * tau0) in (0,1/2). Forbids H=13 & H=15 score(1,2,2,2,3) at m=5.
# M4 (section rotational on grid): forbids 8 classes at m=5 (window vmax=13).
#
# We also add the FULL LRC constraint: instead of arbitrary speed sets, require the
# set to be primitive AND use its TRUE optimal tau0 (M6) / range a over grid (M4).

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
        if cn not in seen:
            seen[cn]=(mn.ham_paths(adj,m),mn.num_3cycles(adj,m),mn.score_seq(adj,m))
    return seen

def realized(adj_fn, m, vmax, uses_grid, also_all_grid_times=False):
    sets=[c for c in combinations(range(1,vmax+1),m) if gcd_list(c)==1]
    real={}
    for S in sets:
        S=list(S)
        if uses_grid:
            for a in range(1,14):
                if gcd(a,14)!=1: continue
                adj=adj_fn(S,a); cn=mn.canon(adj,m)
                real.setdefault(cn,(mn.ham_paths(adj,m),mn.num_3cycles(adj,m),mn.score_seq(adj,m),(tuple(S),a)))
        else:
            if also_all_grid_times:
                # use EVERY candidate tau (not just the optimum) to maximize realized set
                for tau0 in mn.cand(S):
                    adj=adj_fn(S,tau0); cn=mn.canon(adj,m)
                    real.setdefault(cn,(mn.ham_paths(adj,m),mn.num_3cycles(adj,m),mn.score_seq(adj,m),(tuple(S),tau0)))
            else:
                gap,tau0=mn.M(S)
                if tau0 is None: continue
                adj=adj_fn(S,tau0); cn=mn.canon(adj,m)
                real.setdefault(cn,(mn.ham_paths(adj,m),mn.num_3cycles(adj,m),mn.score_seq(adj,m),(tuple(S),tau0)))
    return real

m=5
allc=all_classes(m)
print(f"m={m}: full iso set = {len(allc)}")

# ---- M6 at OPTIMUM tau0, increasing window ----
print("\n[M6 at OPTIMUM tau0] window scan:")
for vmax in [12,16,20,26,32]:
    real=realized(mn.method6_adj,m,vmax,False)
    forb=[cn for cn in allc if cn not in real]
    print(f"  vmax={vmax}: realized {len(real)}/{len(allc)}, forbidden {len(forb)}")
    if vmax==32:
        print("   FORBIDDEN (H,c3,score):")
        for cn in sorted(forb,key=lambda c:(allc[c][0],allc[c][1])):
            print("     ",allc[cn])

# ---- M6 at ALL candidate times (does the optimum specifically forbid them?) ----
print("\n[M6 at ALL candidate tau] window scan (tests if it's the OPTIMUM or the map):")
for vmax in [12,20,32]:
    real=realized(mn.method6_adj,m,vmax,False,also_all_grid_times=True)
    forb=[cn for cn in allc if cn not in real]
    print(f"  vmax={vmax}: realized {len(real)}/{len(allc)}, forbidden {len(forb)}")
    if vmax==32 and forb:
        print("   FORBIDDEN even over ALL tau (H,c3,score):")
        for cn in sorted(forb,key=lambda c:(allc[c][0],allc[c][1])):
            print("     ",allc[cn])

# ---- M4 increasing window ----
print("\n[M4 section rotational on grid] window scan:")
for vmax in [13,18,24,30]:
    real=realized(mn.method4_adj,m,vmax,True)
    forb=[cn for cn in allc if cn not in real]
    print(f"  vmax={vmax}: realized {len(real)}/{len(allc)}, forbidden {len(forb)}")
    if vmax==30:
        print("   FORBIDDEN (H,c3,score):")
        for cn in sorted(forb,key=lambda c:(allc[c][0],allc[c][1])):
            print("     ",allc[cn])
