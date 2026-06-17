"""
tight_maxentry_probe_kps.py  (kind-pasteur, 2026-06-16)

FAST decisive probe of finiteness: among ALL primitive tight 13-sets reachable by
exhaustive small-window search, what is the MAX ENTRY and MAX LCM ever seen?

Strategy: the leave-one-out covering bound shows the 12 smallest speeds must
NEARLY cover [0,1] (a single extra speed covers <=1/7 of any gap). So tightness
forces small speeds. We exhaustively enumerate 13-subsets of {1..M} for modest M
restricted to those CONTAINING a long low run (so they have a chance to be tight),
and record every tight one. We also do full exhaustion of 13-subsets of {1..18}.
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
from itertools import combinations
import time, json

def danger_arcs(v):
    w = Fr(1, 14*v); A = []
    for k in range(v+1):
        c = Fr(k, v); lo = c - w; hi = c + w
        if lo < 0: A += [(Fr(0), hi), (1+lo, Fr(1))]
        elif hi > 1: A += [(lo, Fr(1)), (Fr(0), hi-1)]
        else: A.append((lo, hi))
    return A

# precompute arcs once
_ARC={}
def arcs(v):
    if v not in _ARC: _ARC[v]=danger_arcs(v)
    return _ARC[v]

def L_is_zero(S):
    A=[]
    for v in S: A.extend(arcs(v))
    A=sorted((a,b) for a,b in A if b>a)
    ch=Fr(0); ok=True
    cur=Fr(0)
    # check union covers [0,1]: sweep, track covered prefix
    cl=chh=None; covered=Fr(0)
    for a,b in A:
        if chh is None: cl,chh=a,b
        elif a<=chh: chh=max(chh,b)
        else: covered+=chh-cl; cl,chh=a,b
    if chh is not None: covered+=chh-cl
    return covered==1

def gcd_list(S): return reduce(gcd,S)
def lcm2(a,b): return a*b//gcd(a,b)
def lcm_list(S): return reduce(lcm2,S)

def main():
    t0=time.time()
    tight=[]
    # (A) full exhaustion of 13-subsets of {1..18}
    print("=== exhaustive 13-subsets of {1..18} (C(18,13)=8568) ===",flush=True)
    cnt=0
    for S in combinations(range(1,19),13):
        cnt+=1
        if gcd_list(S)!=1: continue
        if L_is_zero(S):
            tight.append({"S":list(S),"lcm":lcm_list(S),"max":max(S),"src":"exh<=18"})
    print(f"  scanned {cnt}, tight {len(tight)}, t={time.time()-t0:.1f}s",flush=True)
    # (B) exhaustive 13-subsets of {1..20} (C(20,13)=77520)
    print("=== exhaustive 13-subsets of {1..20} (C(20,13)=77520) ===",flush=True)
    cnt=0
    for S in combinations(range(1,21),13):
        cnt+=1
        if gcd_list(S)!=1: continue
        if L_is_zero(S):
            d={"S":list(S),"lcm":lcm_list(S),"max":max(S),"src":"exh<=20"}
            if d not in tight: tight.append(d)
    print(f"  scanned {cnt}, tight {len(tight)}, t={time.time()-t0:.1f}s",flush=True)
    # (C) exhaustive 13-subsets of {1..24} (C(24,13)=2496144) - bigger, decisive on max entry
    print("=== exhaustive 13-subsets of {1..24} (C(24,13)=2,496,144) ===",flush=True)
    cnt=0; seen={tuple(d["S"]) for d in tight}
    for S in combinations(range(1,25),13):
        cnt+=1
        if S in seen: continue
        if gcd_list(S)!=1: continue
        if L_is_zero(S):
            tight.append({"S":list(S),"lcm":lcm_list(S),"max":max(S),"src":"exh<=24"})
            seen.add(S)
    print(f"  scanned {cnt}, tight {len(tight)}, t={time.time()-t0:.1f}s",flush=True)

    ts=sorted(tight,key=lambda d:(d["max"],d["lcm"]))
    print("\n===== ALL TIGHT (exhaustive windows) =====",flush=True)
    for d in ts:
        print(f"  max={d['max']:>3} lcm={d['lcm']:>8} S={d['S']} [{d['src']}]",flush=True)
    print(f"\nTotal distinct tight: {len(ts)}",flush=True)
    if ts:
        print(f"MAX ENTRY over all tight: {max(d['max'] for d in ts)}",flush=True)
        print(f"MAX LCM over all tight:   {max(d['lcm'] for d in ts)}",flush=True)
    out={"tight":ts,"n":len(ts),
         "max_entry":max((d['max'] for d in ts),default=None),
         "max_lcm":max((d['lcm'] for d in ts),default=None),
         "elapsed":time.time()-t0}
    with open("05-knowledge/results/tight_maxentry_probe_kps.json","w") as f:
        json.dump(out,f,indent=2)
    print(f"Elapsed {time.time()-t0:.1f}s",flush=True)

if __name__=="__main__":
    main()
