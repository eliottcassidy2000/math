#!/usr/bin/env python3
"""
kind-pasteur-2026-07-19-S128c88 -- HYP-7915 UPDATE: the I(13,p,1) ACCEPTANCE TEST.

Object (level-1 improper-tuple sieve, the S-T bottleneck; ledger
00-navigation/LRC14-FINITE-CHECK-FEASIBILITY-LEDGER-2026-07-19.md s3-s4):
a 13-tuple of nonzero residues mod prime p is IMPROPER MOD p if no time a/p is
a witness, i.e. for every a in Z_p* some speed v has |v*a mod p| < p/14.
Folding Z_p* by +-1 onto P = {1..(p-1)/2}: with dk = floor(p/14), speed v's
danger-time set is the folded 13-term AP A(w) = {fold(j*w) : j=1..dk},
w = fold(v^{-1}); a tuple is improper iff its distinct folded w-set COVERS P.
Dually, point x is covered exactly by the dk candidates w = fold(j^{-1} x).
WLOG 1 in W (global scaling). This is the natural level-1 discretization;
the (1,...,13) AP is improper mod every p coprime to 14 (its maximizers sit at
denominator 14), which is gate G2 below. Tuples with repeated/+- -collided
residues cover iff a distinct subset covers; we enumerate the distinct-13
stratum (every covering structure appears there by supersetting).

THE TEST (backlog acceptance test, one session): three searches on the SAME
enumeration semantics (each 13-subset of P containing 1 counted exactly once):
  A  naive-lex DFS + count cut (|U| <= dk * remaining)
  B  dual/MRV branching (branch on a least-candidates uncovered point, with
     include-exclude banning) + count cut  -- a fair stand-in for the
     published baseline's "branch-cutting via partial coverage analysis"
  C  B + two repo-flavored VALID cuts:
     capacity cut  -- coverage of U by any r remaining sets <= sum of the
        top-r values |A(w) & U| over allowed w  (comb-capacity; the discrete
        shadow of THM-1198's five-comb bound)
     run cut       -- same, restricted to the longest consecutive uncovered
        run R (the "slow gap"): covering U needs covering R, and
        |A(w) & R| is small for most w  (slow-gap geometry)
Both cuts only ever prune nodes from which no completion exists => A, B, C
enumerate identical cover sets (gate G4 checks equality where complete).

Metrics: nodes (comparable across variants), wall time, covers found
(counted with superset-completion C(avail, 13-d) when U empties early).
Completed primes give exact factors; p=191 runs under node budgets.
"""
import sys, time
from math import comb

def build(p, n_runners=14, k=13):
    h = (p-1)//2
    dk = p // n_runners                      # |r| < p/14  =>  |r| <= dk
    inv = [0]*(p)
    for a in range(1, p):
        inv[a] = pow(a, p-2, p)
    def fold(x):
        x %= p
        return x if x <= h else p - x
    maskA = [0]*(h+1)                        # maskA[w]: bit (x-1) for x in A(w)
    for w in range(1, h+1):
        m = 0
        for j in range(1, dk+1):
            m |= 1 << (fold(j*w) - 1)
        maskA[w] = m
    cand = [[] for _ in range(h+1)]          # cand[x] = {w : x in A(w)}
    for x in range(1, h+1):
        s = set()
        for j in range(1, dk+1):
            s.add(fold(inv[j]*x))
        cand[x] = sorted(s)
    FULL = (1 << h) - 1
    return h, dk, maskA, cand, FULL, fold, inv

# ---------------- search A: naive lex + count cut ----------------
def searchA(p, budget=None):
    h, dk, maskA, cand, FULL, fold, inv = build(p)
    nodes = 0; covers = 0; sols = []
    t0 = time.time()
    def rec(start, depth, U):
        nonlocal nodes, covers
        nodes += 1
        if budget and nodes > budget: raise TimeoutError
        if U == 0:
            covers += comb(h - start + 1, 13 - depth)  # any completion works
            if h <= 25: sols.append(("early", depth))
            return
        rem = 13 - depth
        if rem == 0: return
        if U.bit_count() > dk*rem: return
        for w in range(start, h+1):
            if h - w + 1 < rem: break
            rec(w+1, depth+1, U & ~maskA[w])
    done = True
    try:
        rec(2, 1, FULL & ~maskA[1])          # 1 is always chosen
    except TimeoutError:
        done = False
    return dict(nodes=nodes, covers=covers, time=time.time()-t0, done=done)

# ---------------- searches B / C ----------------
def searchBC(p, use_repo_cuts, budget=None, collect=None):
    h, dk, maskA, cand, FULL, fold, inv = build(p)
    nodes = 0; covers = 0
    t0 = time.time()
    used = [1]; usedset = {1}
    def rec(U, banned):
        nonlocal nodes, covers
        nodes += 1
        if budget and nodes > budget: raise TimeoutError
        depth = len(used)
        if U == 0:
            avail = h - depth - len(banned)
            covers += comb(avail, 13 - depth)
            if collect is not None and h <= 25:
                collect.append(tuple(sorted(usedset)))
            return
        rem = 13 - depth
        if rem == 0: return
        uc = U.bit_count()
        if uc > dk*rem: return
        if use_repo_cuts:
            # capacity cut (global) + run cut (longest consecutive run)
            caps = []
            for w in range(2, h+1):
                if w in usedset or w in banned: continue
                c = (maskA[w] & U).bit_count()
                if c: caps.append(c)
            caps.sort(reverse=True)
            if sum(caps[:rem]) < uc: return
            # longest run of uncovered points
            best_len = 0; best_lo = 0; run = 0; lo = 0
            V = U
            for i in range(h):
                if (V >> i) & 1:
                    if run == 0: lo = i
                    run += 1
                    if run > best_len: best_len, best_lo = run, lo
                else:
                    run = 0
            if best_len > 1:
                R = ((1 << best_len) - 1) << best_lo
                rcaps = []
                for w in range(2, h+1):
                    if w in usedset or w in banned: continue
                    c = (maskA[w] & R).bit_count()
                    if c: rcaps.append(c)
                rcaps.sort(reverse=True)
                if sum(rcaps[:rem]) < best_len: return
        # MRV pivot: among up to 24 lowest uncovered points, the fewest allowed candidates
        bestx = -1; bestcs = None
        V = U; seen = 0
        while V and seen < 24:
            i = (V & -V).bit_length() - 1
            x = i + 1
            cs = [w for w in cand[x] if w not in usedset and w not in banned]
            if not cs: return                 # x uncoverable
            if bestcs is None or len(cs) < len(bestcs):
                bestx, bestcs = x, cs
                if len(cs) <= 1: break
            V &= V - 1; seen += 1
        newban = set()
        for w in bestcs:
            usedset.add(w); used.append(w)
            rec(U & ~maskA[w], banned | newban)
            used.pop(); usedset.discard(w)
            newban.add(w)
    done = True
    try:
        rec(FULL & ~maskA[1], frozenset())
    except TimeoutError:
        done = False
    return dict(nodes=nodes, covers=covers, time=time.time()-t0, done=done)

def main():
    print("== gates ==")
    # G1: A(w)/cand consistency + G2: (1..13) improper + G3: GW improper
    for p in (29, 43, 61, 71, 101, 191):
        h, dk, maskA, cand, FULL, fold, inv = build(p)
        ok = all(w in cand[fold(j*w)] for w in (1, 2, h) for j in (1, dk))
        W0 = {fold(inv[v]) for v in range(1, 14)}
        m = 0
        for w in W0: m |= maskA[w]
        ap_ok = (m == FULL and len(W0) == 13)
        gw = [1,2,3,4,5,6,7,8,9,10,11,13,24]
        Wg = {fold(inv[v % p]) for v in gw}
        mg = 0
        for w in Wg: mg |= maskA[w]
        gw_ok = (mg == FULL)
        print(f"  p={p:3d}: h={h:2d} dk={dk:2d} |cand|={len(cand[1])} "
              f"dualG1={'OK' if ok else 'FAIL'} AP-improper={'OK' if ap_ok else 'FAIL'} "
              f"GW-improper={'OK' if gw_ok else 'FAIL'}")
    # G2b: k=7 flavor at p=31 (n=8): (1..7) improper with danger p/8
    h, dk, maskA, cand, FULL, fold, inv = build(31, n_runners=8, k=7)
    W0 = {fold(inv[v]) for v in range(1, 8)}
    m = 0
    for w in W0: m |= maskA[w]
    print(f"  k=7 p=31 (n=8): dk={dk} AP-improper={'OK' if m == FULL else 'FAIL'}")

    print("\n== completed primes (exact factors) ==")
    print(f"  {'p':>4} {'variant':>8} {'nodes':>12} {'covers':>14} {'time_s':>8} done")
    results = {}
    for p in (29, 43, 61, 71):
        sols_b, sols_c = [], []
        a = searchA(p)
        b = searchBC(p, use_repo_cuts=False, collect=sols_b)
        c = searchBC(p, use_repo_cuts=True, collect=sols_c)
        results[p] = (a, b, c)
        for name, r in (("A-lex", a), ("B-mrv", b), ("C-repo", c)):
            print(f"  {p:>4} {name:>8} {r['nodes']:>12,} {r['covers']:>14,} {r['time']:>8.2f} {r['done']}")
        agree = (a['covers'] == b['covers'] == c['covers'])
        print(f"       cover-count agreement A=B=C: {'OK' if agree else '!! FAIL'}"
              f"   factors: A/B={a['nodes']/max(b['nodes'],1):.1f}x"
              f"  B/C={b['nodes']/max(c['nodes'],1):.2f}x  A/C={a['nodes']/max(c['nodes'],1):.1f}x")

    print("\n== the target: p=191 under node budgets ==")
    BUD = 3_000_000
    for name, fn in (("A-lex", lambda: searchA(191, budget=BUD)),
                     ("B-mrv", lambda: searchBC(191, False, budget=BUD)),
                     ("C-repo", lambda: searchBC(191, True, budget=BUD))):
        r = fn()
        rate = r['nodes']/max(r['time'], 1e-9)
        print(f"  {name:>7}: nodes={r['nodes']:>11,}  covers={r['covers']:,}  "
              f"time={r['time']:.1f}s  rate={rate:,.0f}/s  completed={r['done']}")
    print("\n  (if C completes and A/B do not, the factor is >= budget/C-nodes;")
    print("   otherwise the completed-prime trend is the measured deliverable)")

if __name__ == "__main__":
    main()
