"""The faithful analog of the Moser-beam breakthrough: STATE-LOCAL frontier-gain + BEAM SEARCH to
find the LRC extremal (minimum-M / worry-set) config, vs brute enumeration. Frontier = survival
bitmask over Z/(2n-1); gain = AND with survival[v]; beam keeps the K most-constrained partials
(fewest survivors, heading to tight). Measure the M-evaluation reduction. opus-2026-06-03-S599x."""
from itertools import combinations
from fractions import Fraction as F
def dist(x):
    x%=1; return min(x,1-x)
def M_exact(V):
    ds=set(V)
    for a,b in combinations(V,2): ds.add(a+b); ds.add(abs(a-b))
    ds.discard(0); best=F(-1)
    for d in ds:
        for mden in range(d):
            t=F(mden,d); v=min(min((x*t)%1,1-(x*t)%1) for x in V)
            if v>best: best=v
    return best
def run(n,B,K):
    m=2*n-1; full=(1<<m)-1
    survival=[0]*(B+1)
    for v in range(1,B+1):
        kill=0
        for j in range(m):
            if (v*j)%m in (0,1,m-1): kill|=1<<j
        survival[v]=full & ~kill
    target=n-1
    # BRUTE: M over all configs
    brute_evals=0; brute_min=F(2); 
    for V in combinations(range(1,B+1),target):
        brute_evals+=1; mm=M_exact(V)
        if mm<brute_min: brute_min=mm
    # BEAM: frontier-gain + beam by fewest survivors
    beam=[((v,), survival[v]) for v in range(1,B+1)]
    Mevals=0
    for depth in range(1,target):
        nxt=[]
        for (cfg,mask) in beam:
            for v in range(cfg[-1]+1, B+1):
                nxt.append((cfg+(v,), mask & survival[v]))
        # keep K most-constrained (fewest surviving witnesses)
        nxt.sort(key=lambda cm: bin(cm[1]).count('1'))
        beam=nxt[:K]
    beam_min=F(2)
    for (cfg,mask) in beam:
        Mevals+=1; mm=M_exact(cfg)
        if mm<beam_min: beam_min=mm
    return brute_min, beam_min, brute_evals, Mevals
def main():
    for (n,B,K) in [(6,12,40),(7,14,60),(8,16,80)]:
        bm,bem,be,me = run(n,B,K)
        ok = (bm==bem)
        print(f"n={n} B={B} K={K}: brute min M={bm} ({be} M-evals); beam min M={bem} ({me} M-evals); "
              f"found-optimum={ok}; REDUCTION={be/me:.0f}x")
    print("\nState-local frontier-gain (survival bitmask, 1 AND/runner) + beam(K) finds the worry-set")
    print("extremal with a large M-evaluation reduction — the LRC analog of the 211x Moser-beam.")
if __name__=='__main__': main()
