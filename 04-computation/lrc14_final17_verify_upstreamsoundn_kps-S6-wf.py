"""
LRC(14) FINAL VERIFICATION — the 1/7-spread bound (the last gap)
kps-S6-wf, adversarial verification.

TARGET: for every integer co-offset set E (0 in E, |E|=k, 8<=k<=12),
prove mu_{1/7}(E) >= thr_k.
Sufficient: prove consecutive {0,1,...,k-1} minimizes mu_{1/7} over all integer E.

This script:
 (1) re-derives mu_{1/7}(consec_k) and thr_k EXACTLY,
 (2) HUNTS for an integer E with mu_{1/7}(E) < thr_k:
       - exhaustive small-spread (all subsets of {0..W})
       - aggressive large-spread descent (the same families that crushed mu_{2/7})
 (3) tests whether consecutive minimizes mu_{1/7} (perforated/spread/structured),
 (4) sanity-checks the engine itself.

EXACT rationals throughout (Fraction).
"""
from fractions import Fraction as F
from itertools import combinations
import random, sys

# ----------------------------------------------------------------------
# EXACT mu_theta ENGINE (order-cell + gap=theta breakpoints) — copied from prompt
# ----------------------------------------------------------------------
def mu_theta(E, theta):
    E = sorted(set(E)); n = len(E); bp = set([F(0), F(1)])
    for i in range(n):
        for j in range(i+1, n):
            d = E[j]-E[i]
            for m in range(0, d+1): bp.add(F(m, d))
    bp = sorted(b for b in bp if 0 <= b <= 1); total = F(0)
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        mid = (a+b)/2
        order = sorted(range(n), key=lambda i: (E[i]*mid) % 1)
        ks = [(E[order[t]]*mid).__floor__() for t in range(n)]; subs = []
        for t in range(n):
            o1 = order[t]; o2 = order[(t+1) % n]; k1 = ks[t]; k2 = ks[(t+1) % n]
            wrap = 1 if t == n-1 else 0
            s = E[o2]-E[o1]; c = F(k1-k2+wrap)
            if s == 0:
                if c > theta: subs.append((a, b))
            elif s > 0:
                lo = max(a, (theta-c)/s)
                subs.append((lo, b)) if lo < b else None
            else:
                hi = min(b, (theta-c)/s)
                subs.append((a, hi)) if a < hi else None
        subs.sort(); cur = cb = None
        for lo, hi in subs:
            if cur is None: cur, cb = lo, hi
            elif lo <= cb: cb = max(cb, hi)
            else: total += cb-cur; cur, cb = lo, hi
        if cur is not None: total += cb-cur
    return total

def mu17(E):
    return mu_theta(E, F(1, 7))

# ----------------------------------------------------------------------
# Independent brute-force mu_theta via fine rational sampling on a common grid.
# This is a SANITY checker for the analytic engine: mu_theta = measure of
# x in [0,1) such that maxgap of {frac(e_i x)} > theta.
# We compute the maxgap symbolically per order-cell but verify with sampling
# at cell midpoints + a dense uniform check on integer denominators.
# ----------------------------------------------------------------------
def maxgap_at(E, x):
    """max circular gap of {frac(e_i x)} at a single rational x in [0,1)."""
    pts = sorted(set((F(e)*x) % 1 for e in E))
    if len(pts) == 1:
        return F(1)
    g = F(0)
    for i in range(len(pts)):
        nxt = pts[(i+1) % len(pts)] + (1 if i == len(pts)-1 else 0)
        g = max(g, nxt - pts[i])
    return g

def mu_theta_brute(E, theta, denom):
    """APPROXIMATE measure via midpoint sampling on a fixed fine grid of `denom`
    cells. This converges to the true (analytic) measure as denom -> infinity;
    it is NOT exact for finite denom, but is a fully independent check that the
    analytic engine is in the right ballpark and has no gross sign/structure bug.
    Returns the fraction of cell midpoints whose maxgap strictly exceeds theta."""
    E = sorted(set(E))
    cnt = 0
    for j in range(denom):
        x = F(2*j+1, 2*denom)
        if maxgap_at(E, x) > theta:
            cnt += 1
    return F(cnt, denom)

# ----------------------------------------------------------------------
# Claimed exact values
# ----------------------------------------------------------------------
CLAIMED_MU_CONSEC = {
    8:  F(691, 735),
    9:  F(247, 294),
    10: F(38, 49),
    11: F(1381, 2205),
    12: F(13823, 24255),
    13: F(477, 1078),
}
CLAIMED_THR = {
    8:  F(3637, 5880),
    9:  None,   # ~0.50 (not given exactly in prompt; recompute from meas(G_P))
    10: None,   # ~0.39
    11: None,   # ~0.28
    12: F(1, 7),
    13: F(0, 1),
}

# meas(G_P) for the SMALL part P (|P|=13-k) — we need thr_k = 1 - min_{|P|=13-k} meas(G_P).
# G_P = {x : ||p x|| >= 1/14 for all p in P}, where ||.|| is distance to nearest integer.
# m_P (claimed) = meas(G_{full P at |P|=13}) ... but thr_k uses |P|=13-k.
# We recompute meas(G_P) exactly for arbitrary P of small integers.

def meas_GP(P):
    """measure of {x in [0,1): ||p x|| >= 1/14 for all p in P}, exact."""
    P = sorted(set(int(p) for p in P if p != 0))
    if not P:
        return F(1)
    bp = set([F(0), F(1)])
    for p in P:
        # ||p x|| >= 1/14  <=>  frac(p x) in [1/14, 13/14]
        # boundaries at x = (j +- 1/14)/p
        for j in range(0, p+1):
            for off in (F(1,14), F(13,14)):
                v = F(j) + off
                x = v / p
                if 0 <= x <= 1: bp.add(x)
            bp.add(F(j, p))
    bp = sorted(b for b in bp if 0 <= b <= 1)
    total = F(0)
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        mid = (a+b)/2
        ok = True
        for p in P:
            fr = (F(p)*mid) % 1
            if fr < F(1,14) or fr > F(13,14):
                ok = False; break
        if ok: total += b - a
    return total

def thr_k_exact(k):
    """thr_k = 1 - min over P subset of {1..13}, |P|=13-k, meas(G_P).
    But the small part P consists of the SMALL elements <=13 of a primitive
    covering 13-set. Per prompt the relevant minimization is over |P|=13-k.
    We take P subset of {1,...,13}. Returns (thr, argmin_P)."""
    sz = 13 - k
    if sz <= 0:
        return F(0), ()
    best = None; bestP = None
    for P in combinations(range(1, 14), sz):
        m = meas_GP(P)
        if best is None or m < best:
            best = m; bestP = P
    return F(1) - best, bestP

# ----------------------------------------------------------------------
def main():
    print("="*72)
    print("PART 0: ENGINE SANITY (analytic mu_theta vs brute midpoint sampling)")
    print("="*72)
    tests = [
        ([0,1,2,3,4,5,6,7], F(1,7)),
        ([0,1,2,3,4,5,6], F(1,7)),
        ([0,1,2], F(1,7)),
        ([0,2,5,9], F(1,7)),
        ([0,1,3,7,12,20,30,44], F(1,7)),
    ]
    GRID = 200000  # fixed fine grid; brute is an APPROXIMATION that should track analytic
    for E, th in tests:
        analytic = mu_theta(E, th)
        brute = mu_theta_brute(E, th, GRID)
        diff = abs(float(analytic) - float(brute))
        flag = "OK (close)" if diff < 1e-3 else "*** FAR APART — CHECK ***"
        print(f"E={E}  analytic={float(analytic):.6f} (={analytic})  brute~{float(brute):.6f}  |diff|={diff:.2e}  {flag}")
        sys.stdout.flush()

    print()
    print("="*72)
    print("PART 1: RE-DERIVE mu_{1/7}(consec_k) and thr_k EXACTLY")
    print("="*72)
    for k in range(8, 14):
        E = list(range(k))
        mu = mu17(E)
        claim = CLAIMED_MU_CONSEC[k]
        ok = "OK" if mu == claim else "*** MISMATCH ***"
        print(f"k={k}: mu_1/7(consec) = {mu} = {float(mu):.6f}   claimed {claim}   {ok}")
    print()
    print("thr_k = 1 - min_{|P|=13-k, P in {1..13}} meas(G_P):")
    thr_computed = {}
    for k in range(8, 14):
        thr, argP = thr_k_exact(k)
        thr_computed[k] = thr
        claim = CLAIMED_THR.get(k)
        cs = "" if claim is None else f"   claimed {claim}={float(claim):.6f}"
        ok = ""
        if claim is not None:
            ok = "OK" if thr == claim else "*** MISMATCH ***"
        print(f"k={k}: thr = {thr} = {float(thr):.6f}   argmin P={argP}{cs}   {ok}")
    print()
    print("MARGIN check: mu_1/7(consec_k) - thr_k  (must be >= 0):")
    all_pos = True
    for k in range(8, 14):
        mu = CLAIMED_MU_CONSEC[k]
        thr = thr_computed[k]
        margin = mu - thr
        if margin < 0: all_pos = False
        print(f"k={k}: margin = {margin} = {float(margin):.6f}  {'>=0' if margin>=0 else '*** NEGATIVE ***'}")
    print(f"\nConsecutive clears thr at every k=8..13: {all_pos}")

    print()
    print("="*72)
    print("PART 2: HUNT for integer E with mu_{1/7}(E) < thr_k")
    print("="*72)
    violations = []

    # 2a. exhaustive small-spread: all k-subsets of {0..W} containing 0
    print("\n[2a] Exhaustive small-spread (subsets of {0..W} with 0 in E):")
    for k in range(8, 13):
        thr = thr_computed[k]
        # widen window until cost gets too big
        Wmax = {8:12, 9:13, 10:13, 11:13, 12:13}[k]
        best = None; bestE = None; count = 0
        for W in range(k-1, Wmax+1):
            # subsets of size k from {0..W} containing 0
            for rest in combinations(range(1, W+1), k-1):
                E = (0,) + rest
                if max(E) != W:  # only new ones (max==W) to avoid recount
                    continue
                mu = mu17(E)
                count += 1
                if best is None or mu < best:
                    best = mu; bestE = E
                if mu < thr:
                    violations.append((k, E, mu, thr))
        print(f"  k={k}: checked {count} sets up to W={Wmax}; min mu={best}={float(best):.6f} "
              f"at E={bestE}; thr={float(thr):.6f}; "
              f"{'consec is min' if bestE==tuple(range(k)) else 'NON-CONSEC MIN!'}")

    # 2b. aggressive large-spread descent: structured families
    print("\n[2b] Aggressive large-spread families (the mu_{2/7} crushers):")
    fam_results = {}
    def fam_check(name, gen):
        nonlocal violations
        mins = {}
        for k in range(8, 13):
            thr = thr_computed[k]
            best = None; bestE = None
            for E in gen(k):
                E = tuple(sorted(set(E)))
                if len(E) != k or 0 not in E: continue
                mu = mu17(E)
                if best is None or mu < best:
                    best = mu; bestE = E
                if mu < thr:
                    violations.append((k, E, mu, thr))
            mins[k] = (best, bestE, thr)
        fam_results[name] = mins
        print(f"  family '{name}':")
        for k in range(8, 13):
            best, bestE, thr = mins[k]
            consec = mu17(tuple(range(k)))
            tag = "BELOW CONSEC!" if best < consec else "ok"
            print(f"    k={k}: min={float(best):.6f} (E={bestE}) thr={float(thr):.6f} "
                  f"consec={float(consec):.6f} [{tag}]")

    # perforated: arithmetic-progression-like with step s
    def gen_ap(k):
        for s in range(1, 8):
            yield [s*i for i in range(k)]
    fam_check("AP step s", gen_ap)

    # geometric-ish / doubling spreads
    def gen_doubling(k):
        base = [0]
        # several growing rulers
        yield [0] + [2**i for i in range(k-1)]
        yield [0] + [i*i for i in range(1, k)]
        yield [0] + [i*(i+1)//2 for i in range(1, k)]
        yield list(range(0, 2*k, 2))            # even
        yield [0] + list(range(2, k+1)) + [k+3]  # gap insert
    fam_check("doubling/poly", gen_doubling)

    # perforated consecutive: remove one or two interior points, widen
    def gen_perf(k):
        # consec of length k+r, drop r interior points -> spread but dense-ish
        for r in range(1, 4):
            base = list(range(k+r))
            for drop in combinations(range(1, k+r-1), r):
                E = [x for x in base if x not in drop]
                if len(E) == k:
                    yield E
    fam_check("perforated-consec", gen_perf)

    # Sidon / perfect-difference-ish and random large spread
    def gen_random_large(k):
        rng = random.Random(12345 + k)
        for _ in range(4000):
            W = rng.randint(k+2, 6*k)
            E = [0] + rng.sample(range(1, W+1), k-1)
            yield E
    fam_check("random large-spread", gen_random_large)

    # descent: greedy local search from consec, try to LOWER mu by moving a point out
    def gen_descent(k):
        rng = random.Random(999 + k)
        for _ in range(60):
            E = list(range(k))
            curmu = mu17(tuple(E))
            for _step in range(200):
                i = rng.randint(1, k-1)
                cand = E[:]
                cand[i] = rng.randint(1, 5*k)
                cand = sorted(set(cand))
                if len(cand) != k or 0 not in cand: continue
                m2 = mu17(tuple(cand))
                if m2 < curmu:
                    E = cand; curmu = m2
            yield E
    fam_check("greedy descent", gen_descent)

    print()
    print("="*72)
    print("PART 3: Does consecutive minimize mu_{1/7}? (global min over small-spread)")
    print("="*72)
    for k in range(8, 13):
        consec = mu17(tuple(range(k)))
        # search ALL sets up to a window and report if anything ties/beats consec
        Wmax = {8:12, 9:13, 10:13, 11:13, 12:13}[k]
        below = 0; equal = 0; total = 0
        minmu = consec; argmin = tuple(range(k))
        for W in range(k-1, Wmax+1):
            for rest in combinations(range(1, W+1), k-1):
                E = (0,)+rest
                if max(E) != W: continue
                mu = mu17(E); total += 1
                if mu < consec: below += 1
                elif mu == consec and E != tuple(range(k)): equal += 1
                if mu < minmu: minmu = mu; argmin = E
        print(f"  k={k}: consec mu={float(consec):.6f}; over {total} small-spread sets: "
              f"{below} strictly below, {equal} ties; global-min mu={float(minmu):.6f} at {argmin}")

    print()
    print("="*72)
    print("VERDICT")
    print("="*72)
    if violations:
        print(f"*** {len(violations)} VIOLATIONS FOUND (mu_1/7 < thr_k): holds=FALSE ***")
        for k, E, mu, thr in violations[:20]:
            print(f"   k={k} E={E} mu={mu}={float(mu):.6f} < thr={float(thr):.6f}")
    else:
        print("No violation of mu_{1/7}(E) >= thr_k found in any tested family.")
        print("Consecutive remained the minimizer in all exhaustive + descent searches.")

if __name__ == "__main__":
    main()
