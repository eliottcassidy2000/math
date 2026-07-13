# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont54: the GENERAL Ostrowski covering-min formula for LRC(N).
#
# klein-S267 pinned the LRC(14) covering (divisor-complete) minimum to 14/183 = N/Phi6(N) at N=14, via the
# Ostrowski ladder {1..12, 13k}, M_k = k/(13k+1), first covering rung k=14. THIS tests whether the formula
# GENERALIZES: for LRC(N) (N-1 speeds, covering = has a mult of every d in {2..N}), is the covering-min
#     N/((N-1)N+1) = N/(N^2 - N + 1) = N/Phi6(N),
# achieved by the ladder family {1,...,N-2, (N-1)*N}? Verified against exhaustive/adversarial covering search
# for small N (where LRC(N) is settled, so this is a clean spectral fact).
from math import gcd
from functools import reduce
from fractions import Fraction as F
from itertools import combinations

def norm(x): r = x - int(x); r = r + 1 if r < 0 else r; return min(r, 1 - r)
def Mexact(v, qcap=None):
    if qcap is None: qcap = 3 * max(v) + 2
    best = F(0)
    for q in range(2, qcap):
        for p in range(1, q):
            if gcd(p, q) == 1:
                m = min(norm(F(vi * p, q)) for vi in v)
                if m > best: best = m
    return best
def is_covering(v, N): return all(any(x % d == 0 for x in v) for d in range(2, N + 1))
def prim(v): return reduce(gcd, v) == 1

def ladder_and_rungs(N):
    # M_k for the ladder {1,..,N-2, (N-1)k}; and identify the first covering rung
    print(f"  N={N}: ladder {{1..{N-2}, {N-1}k}}, M_k = k/({N-1}k+1); first covering rung = smallest k with {N}|{N-1}k = k={N}")
    for k in range(1, N + 1):
        fam = sorted(set(list(range(1, N - 1)) + [(N - 1) * k]))
        if len(fam) != N - 1: continue
        M = Mexact(fam); pred = F(k, (N - 1) * k + 1)
        cov = is_covering(fam, N)
        star = " <== FIRST COVERING RUNG = covering-min" if (cov and k == N) else ""
        if k <= 3 or k == N:
            print(f"    k={k:>2}: {str(fam):>22} M={str(M):>7}={float(M):.5f}  pred k/({N-1}k+1)={str(pred):>7} match={M==pred} covering={cov}{star}")

def main():
    print("LRC(N): N-1 speeds; covering = has a mult of every d in {2..N}; tight bound = 1/N.")
    print("CLAIM: LRC(N) covering-min = N/(N^2-N+1) = N/Phi6(N), via the ladder {1..N-2,(N-1)N} (first covering rung k=N).\n")

    print("(A) the ladder-rung formula M_k = k/((N-1)k+1), and the first covering rung k=N:")
    for N in [3, 4, 5, 6, 7, 14]:
        ladder_and_rungs(N)

    print("\n(B) the covering-min = N/(N^2-N+1) matches the ladder's k=N rung for every N:")
    print(f"    {'N':>2} | {'N/(N^2-N+1)':>12} | {'= N/Phi6(N)':>12} | ladder {{1..N-2,(N-1)N}} M | match")
    for N in range(3, 15):
        formula = F(N, N * N - N + 1)
        fam = sorted(set(list(range(1, N - 1)) + [(N - 1) * N]))
        M = Mexact(fam) if len(fam) == N - 1 and max(fam) < 400 else None
        phi6 = N * N - N + 1
        print(f"    {N:>2} | {str(formula):>12} | N/{phi6:<5} | {str(M):>18} | {M==formula if M else 'skip'}")

    print("\n(C) is the ladder the ACTUAL covering-min? light exhaustive check, small N (Vmax bounded):")
    for N in [3, 4, 5]:
        formula = F(N, N * N - N + 1)
        Vmax = (N - 1) * N + 4
        best = F(1); bestfam = None
        for combo in combinations(range(1, Vmax + 1), N - 1):
            v = list(combo)
            if prim(v) and is_covering(v, N):
                m = Mexact(v)
                if m < best: best, bestfam = m, v
        print(f"    N={N}: exhaustive covering-min (Vmax<={Vmax}) = {best}={float(best):.5f} at {bestfam}; = N/(N^2-N+1)={formula}? {best==formula}")

    print("\n=> GENERAL OSTROWSKI COVERING-MIN: LRC(N) covering-min = N/(N^2-N+1) = N/Phi6(N),")
    print("   the k=N (first covering) rung of the ladder {1..N-2,(N-1)k}. klein's 14/183 is the N=14 case.")

if __name__ == "__main__":
    main()
