"""
LRC(14) — structural probe of "consecutive minimizes mu_{1/7}".
Three adversarial angles on the OPEN gap (consecutive = global minimizer):

 1. COMPRESSION MONOTONICITY: is mu_{1/7} monotone non-increasing under any
    single-point "compression" toward making E consecutive? If TRUE universally,
    that would PROVE consecutive is the minimizer. We look for a counterexample
    (a compression step that INCREASES mu) — if found, the easy monotonicity
    proof route is dead (but consecutive could still be the min).

 2. UNBOUNDED-SPREAD STRESS: does mu_{1/7}(E) -> 1 as spread -> infinity for
    fixed k? (If perforation only ever RAISES mu, the minimizer is the densest
    set = consecutive.) We track mu vs spread on random sets.

 3. NEAR-TIES: are there NON-consecutive sets achieving mu EQUAL to consec (which
    would still satisfy the bound but break strict uniqueness)? Report any.

EXACT rationals.
"""
from fractions import Fraction as F
from itertools import combinations
import random, sys

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
                lo = max(a, (theta-c)/s); subs.append((lo, b)) if lo < b else None
            else:
                hi = min(b, (theta-c)/s); subs.append((a, hi)) if a < hi else None
        subs.sort(); cur = cb = None
        for lo, hi in subs:
            if cur is None: cur, cb = lo, hi
            elif lo <= cb: cb = max(cb, hi)
            else: total += cb-cur; cur, cb = lo, hi
        if cur is not None: total += cb-cur
    return total

def mu17(E): return mu_theta(tuple(sorted(set(E))), F(1,7))

def normalize(E):
    E = sorted(set(int(x) for x in E)); g = 0
    import math
    m = E[0]
    E = [x-m for x in E]
    g = 0
    for x in E: g = math.gcd(g,x)
    if g>1: E=[x//g for x in E]
    return tuple(E)

def main():
    CONSEC = {k: mu17(range(k)) for k in range(8,13)}

    print("="*70)
    print("ANGLE 1: COMPRESSION MONOTONICITY (hunt for a compression that RAISES mu)")
    print("="*70)
    # A 'compression' = move one interior point one step toward its neighbor,
    # reducing total spread by 1 without creating a collision. If mu always
    # NON-INCREASES under compression, consecutive (max compression) is the min.
    rng = random.Random(101)
    raise_examples = []
    checked = 0
    for trial in range(8000):
        k = rng.randint(8,12)
        cap = rng.randint(k+1, 3*k)
        E = sorted(set([0]+rng.sample(range(1,cap+1),k-1)))
        if len(E)!=k: continue
        muE = mu17(E)
        # try every compression: pick a gap > 1 and shrink it by pulling the
        # upper part down by 1 (preserves order, reduces spread).
        for i in range(1,k):
            if E[i]-E[i-1] >= 2:
                E2 = E[:i] + [x-1 for x in E[i:]]
                E2 = tuple(E2)
                if len(set(E2))!=k: continue
                checked += 1
                mu2 = mu17(E2)
                if mu2 > muE:   # compression RAISED mu -> monotonicity FALSE
                    raise_examples.append((E, muE, E2, mu2))
                    if len(raise_examples) >= 8: break
        if len(raise_examples) >= 8: break
    print(f"compression steps checked: {checked}")
    if raise_examples:
        print(f"FOUND {len(raise_examples)} compression(s) that RAISE mu => "
              f"simple monotonicity is FALSE (consec-min proof needs more than compression):")
        for E,m,E2,m2 in raise_examples[:8]:
            print(f"   {E} mu={float(m):.6f}  --compress-->  {E2} mu={float(m2):.6f}  (rose by {float(m2-m):.6f})")
    else:
        print("NO compression raised mu in 8000 trials => mu is empirically MONOTONE "
              "non-increasing under compression. Strongly supports consec=min "
              "(but this is evidence, not a proof of the monotonicity lemma).")
    sys.stdout.flush()

    print()
    print("="*70)
    print("ANGLE 2: UNBOUNDED-SPREAD — does mu trend UP with spread? (perforation hurts)")
    print("="*70)
    for k in [8, 10, 12]:
        print(f"  k={k}: consec mu={float(CONSEC[k]):.6f}")
        rng = random.Random(7+k)
        buckets = {}
        for _ in range(3000):
            cap = rng.randint(k, 8*k)
            E = tuple(sorted(set([0]+rng.sample(range(1,cap+1),k-1))))
            if len(E)!=k: continue
            spread = max(E)
            band = spread // k  # roughly spread/k
            mu = mu17(E)
            buckets.setdefault(band, []).append(float(mu))
        for band in sorted(buckets)[:8]:
            vs = buckets[band]
            print(f"     spread/k~{band}: n={len(vs)} mean mu={sum(vs)/len(vs):.4f} "
                  f"min={min(vs):.4f} max={max(vs):.4f}")
    print("  (min mu should stay >= consec; mean should rise as spread grows.)")
    sys.stdout.flush()

    print()
    print("="*70)
    print("ANGLE 3: NEAR-TIES — non-consecutive sets with mu == consec? (uniqueness)")
    print("="*70)
    for k in range(8,13):
        ties = []
        # exhaustive over small window + scaled copies (scaled copies always tie)
        Wmax = 13
        cnt = 0
        for W in range(k-1, Wmax+1):
            for rest in combinations(range(1,W+1), k-1):
                E = (0,)+rest
                if max(E)!=W: continue
                cnt += 1
                if mu17(E) == CONSEC[k] and normalize(E)!=tuple(range(k)):
                    ties.append(E)
        print(f"  k={k}: {cnt} small-spread sets; non-consec ties with consec: "
              f"{len(ties)}" + (f"  e.g. {ties[:3]}" if ties else " (consec is STRICT min on this window)"))
    sys.stdout.flush()

if __name__=="__main__":
    main()
