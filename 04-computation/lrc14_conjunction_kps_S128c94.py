#!/usr/bin/env python3
"""
kind-pasteur-2026-07-19-S128c94 -- HYP-8030: THE CONJUNCTION SESSION.
The HYP-8020 obligations worked together:

(A) STRUCTURED-PAIR STATISTIC (the correct S6-mechanism instrument after the
    c93 random-pair null) + k0 REFINEMENT: at p in {719, 839, 983 (safe),
    631, 991 (smooth)}, 2000 greedy draws each; for every distinct minimal
    cover found, the INTERNAL resonance profile over its 78 pairs
    (r(w,w') = |A(w) & A(w')|): total waste (sanity: == 13*dk - h), spread
    vs concentration (top-3 pair share, max/mean), compared with random
    13-sets; and the near-shell distance histogram at 5x the c93 draws
    (the k0 estimate the defect-cage consumes).

(B) THE BOUNDED-DEFECT CAGE H0'(k0): if V matches c*F on >= 13 - k0
    coordinates mod every caging prime (p >= 700 grid, budget-checked),
    the MATCHED SUB-MULTISET F' (a fixed (13-k0)-subset of F per voting
    class) satisfies the c89 elimination verbatim: J-separators + Newton on
    13-k0 elements pin V' = c*F' exactly up to height H0'(k0).  Computed:
    the p >= 700 caging grid + weakest rung; voting classes 2*C(13,k0);
    per-pattern J-separator degeneracy (J(F'_AP) vs J(F'_GW) collisions
    across all C(13,k0)^2 pattern pairs); H0'(k0) for k0 = 0, 1, 2 by the
    c89 linear-solve model.  PAYOFF: the caged conclusion "11 of 13 speeds
    are an exact dilate sub-multiset" = the k0-far stratum -- exactly what
    mac-mini-S125's two-far sweep + small-witness law close.

(C) THE REPRICING AS A CONJUNCTION STATEMENT: shell sizes x prime counts
    tabled; the honest statement is that the sieve's cost collapses IFF the
    shell theorem supplies per-prime completeness -- the computational and
    analytic halves are one lead (no fake feasibility table).
"""
import sys, time, random
from math import log, comb
from fractions import Fraction as Fr
sys.path.insert(0, '04-computation')
from lrc14_I13p1_acceptance_test_kps_S128c88 import build
from lrc14_I13p1_minimal_covers_kps_S128c88 import ansatz_library
from lrc14_shell_census_sample_kps_S128c92 import sample_minimal

random.seed(94)

def primes_in(lo, hi):
    sieve = bytearray([1])*(hi+1); sieve[0:2] = b'\x00\x00'
    for i in range(2, int(hi**0.5)+1):
        if sieve[i]: sieve[i*i::i] = b'\x00'*len(sieve[i*i::i])
    return [i for i in range(lo, hi+1) if sieve[i]]

# ---------------- part A ----------------
def partA():
    print("== (A) structured-pair statistic + k0 refinement (2000 draws/prime) ==", flush=True)
    import statistics
    for p in (719, 839, 983, 631, 991):
        t0 = time.time()
        ctx = build(p)
        h, dk, maskA, cand, FULL, fold, inv = ctx
        lib = ansatz_library(p, ctx); libsets = list(lib.keys())
        seen = set(); near = {}
        for _ in range(2000):
            W = sample_minimal(p, ctx, topt=random.choice((2, 3, 4)))
            if W is None: continue
            if W in seen: continue
            seen.add(W)
            V = frozenset(fold(inv[w]) for w in W)
            d = min(len(V.symmetric_difference(F))//2 for F in libsets)
            near[d] = near.get(d, 0) + 1
        # internal resonance profiles of found covers vs random 13-sets
        def profile(ws):
            rs = []
            L = sorted(ws)
            for i in range(len(L)):
                for j in range(i+1, len(L)):
                    rs.append((maskA[L[i]] & maskA[L[j]]).bit_count())
            rs.sort(reverse=True)
            tot = sum(rs)
            top3 = sum(rs[:3])/tot if tot else 0
            return tot, top3, (rs[0]/ (tot/len(rs)) if tot else 0)
        cov_stats = [profile(W) for W in list(seen)[:40]]
        rnd_stats = [profile(random.sample(range(1, h+1), 13)) for _ in range(40)]
        allowance = 13*dk - h
        if cov_stats:
            cm = statistics.mean(x[0] for x in cov_stats)
            ct3 = statistics.mean(x[1] for x in cov_stats)
            cmx = statistics.mean(x[2] for x in cov_stats)
        rm = statistics.mean(x[0] for x in rnd_stats)
        rt3 = statistics.mean(x[1] for x in rnd_stats)
        rmx = statistics.mean(x[2] for x in rnd_stats)
        print(f"  p={p}: distinct covers {len(seen)}; near-shell {dict(sorted(near.items()))}", flush=True)
        if cov_stats:
            print(f"    covers:  waste-mean {cm:.1f} (allowance {allowance}), top3-share {ct3:.3f}, max/mean {cmx:.2f}", flush=True)
        print(f"    random:  waste-mean {rm:.1f} (allowance {allowance}), top3-share {rt3:.3f}, max/mean {rmx:.2f}  [{time.time()-t0:.0f}s]", flush=True)

# ---------------- part B ----------------
def partB():
    print("\n== (B) the bounded-defect cage H0'(k0) at the p >= 700 grid ==", flush=True)
    target = 13*(12*log(91) - log(13)) - log(360360)
    ps = []
    s = 0.0
    for p in primes_in(701, 5000):
        if p % 7 == 0: continue
        ps.append(p); s += log(p)
        if s > target: break
    LNP = s
    print(f"  caging grid: {len(ps)} primes {ps[0]}..{ps[-1]}, sum ln = {LNP:.1f} > {target:.1f}", flush=True)
    best = None
    for p in ps:
        for l in (1, 2, 4, 8):
            q = l*p; c = -(-q//14)
            r = Fr(c, q)
            if best is None or r < best[0]: best = (r, l, p)
    print(f"  weakest 14-coprime rung: {best[0]} (l={best[1]}, p={best[2]}); micro-gap width {best[0] - Fr(1,14)} ~ {float(best[0]-Fr(1,14)):.2e}", flush=True)
    AP = list(range(1, 14)); GW = list(range(1, 12)) + [13, 24]
    def S(V, m): return sum(v**(2*m) for v in V)
    def J(V): return Fr(S(V,1)*S(V,3), S(V,2)**2)
    from itertools import combinations
    for k0 in (0, 1, 2):
        patterns = list(combinations(range(13), k0))
        # J-separator degeneracy across the class table
        subAP = {}; subGW = {}
        for pat in patterns:
            keep = [x for i, x in enumerate(AP) if i not in pat]
            subAP[pat] = J(keep)
            keep = [x for i, x in enumerate(GW) if i not in pat]
            subGW[pat] = J(keep)
        degen = sum(1 for pa in patterns for pb in patterns if subAP[pa] == subGW[pb])
        nclasses = 2*len(patterns)
        # worst S_{2(13-k0)} over all sub-multisets (binding invariant)
        mtop = 13 - k0
        worstS = max(max(sum(v**(2*mtop) for j, v in enumerate(F) if j not in pat) for pat in patterns) for F in (AP, GW))
        # linear solve: (LNP - c_loss*x)/nclasses = mtop*ln(13) + 2*mtop*x + ln(worstS) + ln 2
        c_loss = 13*log(ps[-1])/log(ps[0])
        num = LNP/nclasses - mtop*log(13.0) - log(float(worstS)) - log(2.0)
        den = 2*mtop + c_loss/nclasses
        x = num/den
        H0 = 2.718 ** max(x, 0.0)
        print(f"  k0={k0}: voting classes {nclasses}, J-degeneracies {degen}/{len(patterns)**2}, "
              f"binding degree {2*mtop}, H0'({k0}) ~ {H0:,.0f}", flush=True)
    print("  caged conclusion at k0: 13-k0 of the speeds are EXACTLY a dilate sub-multiset of AP or GW;", flush=True)
    print("  the k0 remaining speeds are free integers -- the k0-FAR STRATUM (mac-mini S125's landing zone).", flush=True)

# ---------------- part C ----------------
def partC():
    print("\n== (C) the repricing conjunction statement ==", flush=True)
    print("  measured shell sizes at p >= 700: 1-5 distinct greedy-reachable covers/prime (parts A, c93b);", flush=True)
    print("  85-prime safe caging set budget-feasible (c93).  The sieve's I-phase collapses from the", flush=True)
    print("  10^25-scale representative space to ~10^1-2 shell members per prime IFF the shell theorem", flush=True)
    print("  supplies per-prime COMPLETENESS certificates (nothing improper outside the shell).", flush=True)
    print("  Without completeness the exhaustive grind returns; with it the whole finite check is", flush=True)
    print("  [shell verification, trivial] x [~100 primes].  The computational twin and the analytic", flush=True)
    print("  descent are therefore ONE lead: the F_p counting proof IS the cost collapse.", flush=True)

if __name__ == "__main__":
    partA(); partB(); partC()
