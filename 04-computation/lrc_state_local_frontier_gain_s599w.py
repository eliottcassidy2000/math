"""Apply the 'state-local frontier-gain table' breakthrough to LRC. Global M(S)=max_t min_i||v_i t||
is O(m^3 B) per config (scan all arrangement vertices). STATE-LOCAL version: precompute, per speed
v, a KILL-MASK over Z/(2n-1) — the bit j set iff v·j mod (2n-1) ∈ {0,±1} (the discrete witness,
S599i). Loneliness = AND of survival masks (=~kill) over runners is NONZERO: ONE bitwise AND per
runner (frontier-gain), no continuous scan. Measure the edge-count-check reduction. opus-2026-06-03-S599w."""
from itertools import combinations
from fractions import Fraction as F
import random
def dist_f(x):
    x-=int(x)
    if x<0:x+=1
    return min(x,1-x)
def brute_M_ops(V):
    # exact-ish max-min by scanning breakpoints; return (M_float, op_count)
    ds=set(V)
    for a,b in combinations(V,2): ds.add(a+b); ds.add(abs(a-b))
    ds.discard(0); best=-1.0; ops=0
    for d in ds:
        for mden in range(d):
            t=mden/d
            v=min(dist_f(x*t) for x in V); ops+=len(V)
            if v>best: best=v
    return best,ops
def main():
    n=14; m=2*n-1; B=2*n; full=(1<<m)-1
    # FRONTIER-GAIN TABLE: survival mask per speed (state-local objective)
    survival={}
    for v in range(1,B+1):
        kill=0
        for j in range(m):
            r=(v*j)%m
            if r in (0,1,m-1): kill|=(1<<j)
        survival[v]=full & ~kill
    def state_local_lonely(V):
        s=full; ops=0
        for v in V:
            s&=survival[v]; ops+=1            # ONE bitwise AND per runner = the frontier-gain step
            if s==0: break
        return (s!=0), ops
    rng=random.Random(7); N=6000
    tot_brute_ops=0; tot_sl_ops=0; n_fast_resolved=0; n_need_cont=0; cont_ops=0
    for _ in range(N):
        V=tuple(sorted(rng.sample(range(1,B+1), n-1)))
        slL, slo = state_local_lonely(V); tot_sl_ops+=slo
        Mf, bo = brute_M_ops(V); tot_brute_ops+=bo
        if slL:                      # discrete witness fires => provably lonely (M>=2/(2n-1)) — DONE, skip continuous
            n_fast_resolved+=1
        else:                        # residual: must do the continuous scan (the 'Moser beam')
            n_need_cont+=1; cont_ops+=bo
    print(f"n={n}, modulus 2n-1={m}, {N} random configs:")
    print(f"  BRUTE total min-evals      = {tot_brute_ops:,}")
    print(f"  STATE-LOCAL frontier ANDs   = {tot_sl_ops:,}  (one AND per runner)")
    print(f"  discrete witness RESOLVED   = {n_fast_resolved}/{N} = {100*n_fast_resolved/N:.1f}% (provably lonely, no scan)")
    print(f"  residual needing scan       = {n_need_cont}/{N} = {100*n_need_cont/N:.1f}%")
    combined = tot_brute_ops / (tot_sl_ops + cont_ops)
    print(f"  EDGE-COUNT-CHECK REDUCTION  = brute / (frontier-ANDs + residual-scans) = {combined:.1f}x")
    print(f"  per-config fast-check speedup (avg brute ops / avg ANDs) = {(tot_brute_ops/N)/(tot_sl_ops/N):.0f}x")
if __name__=='__main__': main()
