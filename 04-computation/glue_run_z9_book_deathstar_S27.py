#!/usr/bin/env python3
"""death-star-2026-07-16-S27 (HYP-7082): (1) THE CAP-GLUE RUN — lattice-mass bound table +
explicit thresholds + pruned scan extension (second-witness closure of the residue-6 sign at
the stronger constant 535/7203; the sign itself already closed by codex-THM-910 at 32/343);
(2) THE Z9 CYCLIC 2-PAGE BOOK — page assignment by parallel classes vs Guy's Z(n).

GLUE DESIGN (cap A, the binding one; B and C analogous with more slack):
  value(V) = P_V(S_a)+P_V(S_b) <= 48/2401 + PAIR(V) + MASS(V) + REM(V), where
  - PAIR(V) = 24-assignment pair-layer, bounded per pair by the EXACT max deviation table
    dev*(p,q) = max over sector pairs |hit(p:q) - 1/49| for the reduced ratio (p:q)
    (computed exactly for pq <= PQ0; beyond: dev* <= 12/(7 p q), codex's D-law);
  - MASS(V) = sum over relations k (supp>=3, |k|<=K0) of the corner-max table M*(k)
    (THM-899/906/911 Bernoulli forms, maximized over sector assignments);
  - REM(V): the O(1/min-v) transient — handled by scanning the finite window exactly.
  The scan: primitive quadruples, maxspeed in (32, T], value evaluated EXACTLY by one
  breakpoint sweep per tuple (all 24+24 assignments in one pass); tuples with analytic
  bound U(V) < 1/12 are SKIPPED (pruning); report: any exact value > candidates?
"""
from fractions import Fraction as Fr
from math import gcd
from itertools import permutations, combinations, product as iproduct
import sys, time

SA = (1,3,5,6); SB = (2,3,4,6); SC = (1,2,4,5)

def value_sweep(v4, want_sets):
    """ONE pass: exact measures of {x: sector-tuple has set == S} for each S in want_sets
    (distinct sectors). Returns dict S -> measure."""
    vs = v4
    bps = sorted(set(Fr(k, 7*v) for v in vs for k in range(7*v+1)))
    out = {S: Fr(0) for S in want_sets}
    for i in range(len(bps)-1):
        a, b = bps[i], bps[i+1]
        mid = (a+b)/2
        secs = tuple(int((v*mid % 1)*7) for v in vs)
        if len(set(secs)) == len(vs):
            fs = frozenset(secs)
            for S in want_sets:
                if fs == frozenset(S):
                    out[S] += b - a
                    break
    return out

def capA_value(v4):
    m = value_sweep(v4, [SA, SB])
    return m[SA] + m[SB]

def capB_value(v4):
    m = value_sweep(v4, [SC])
    return m[SC]

# ---- pair deviation table (exact, reduced ratios) ----
_pair_cache = {}
def pair_dev_max(p, q):
    """max over sector pairs (s,t) of |meas{sec(px)=s, sec(qx)=t} - 1/49| for reduced p<q."""
    key = (p, q)
    if key in _pair_cache: return _pair_cache[key]
    best = Fr(0)
    bps = sorted(set(Fr(k, 7*p) for k in range(7*p+1)) | set(Fr(k, 7*q) for k in range(7*q+1)))
    meas = {}
    for i in range(len(bps)-1):
        a, b = bps[i], bps[i+1]
        mid = (a+b)/2
        st = (int((p*mid % 1)*7), int((q*mid % 1)*7))
        meas[st] = meas.get(st, Fr(0)) + b - a
    for s in range(7):
        for t in range(7):
            d = abs(meas.get((s,t), Fr(0)) - Fr(1,49))
            if d > best: best = d
    _pair_cache[key] = best
    return best

def B3(x): return x**3 - Fr(3,2)*x**2 + Fr(1,2)*x
def B4(x): return x**4 - 2*x**3 + x**2 - Fr(1,30)
def fr(x): return x - (x.numerator // x.denominator)

def corner_max(k):
    """max over singleton-sector assignments of |corner sum| for relation k (supp 3 or 4)."""
    n = len(k)
    prodk = 1
    for ki in k: prodk *= abs(ki)
    best = Fr(0)
    Bf, dc = (B3, 6) if n == 3 else (B4, 24)
    for secs in iproduct(range(7), repeat=n):
        tot = Fr(0)
        for eps in iproduct([0,1], repeat=n):
            arg = Fr(sum(ki*(ci+ei) for ki,ci,ei in zip(k,secs,eps)), 7)
            tot += (-1)**sum(eps) * Bf(fr(arg))
        val = abs(tot) / (dc * prodk)
        if val > best: best = val
    return best

def part1_tables():
    print("GLUE PART 1a: the corner-max table M*(k) (supp-3/4, heights <= 4)")
    ks = []
    for supp, rng in [(3, range(-3,4)), (4, range(-2,3))]:
        seen = set()
        for k in iproduct(rng, repeat=supp):
            if 0 in k: continue
            if k[0] < 0: continue
            g = 0
            for x in k: g = gcd(g, abs(x))
            if g != 1: continue
            kk = tuple(k)
            if kk in seen: continue
            seen.add(kk)
            ks.append(kk)
    tab = []
    for k in ks:
        tab.append((corner_max(k), k))
    tab.sort(reverse=True)
    tot3 = sum(t[0] for t in tab if len(t[1]) == 3)
    tot4 = sum(t[0] for t in tab if len(t[1]) == 4)
    print(f"  top-6: {[(str(t[0]), t[1]) for t in tab[:6]]}")
    print(f"  SUM over table: supp-3 total = {float(tot3):.5f}, supp-4 total = {float(tot4):.5f}")
    print(f"  crude all-relations mass bound (if EVERY tabled relation were present):"
          f" {float(tot3/7 + tot4):.5f} per assignment-set")
    # pair dev table extremes
    print("GLUE PART 1b: pair-deviation extremes dev*(p,q)")
    worst = []
    for q in range(2, 9):
        for p in range(1, q):
            if gcd(p,q) != 1: continue
            worst.append((pair_dev_max(p,q), (p,q)))
    worst.sort(reverse=True)
    print(f"  top-5: {[(f'{float(w[0]):.4f}', w[1]) for w in worst[:5]]}")
    return tab, worst

def part1_scan(T=44, budget=380):
    print(f"GLUE PART 1c: pruned exact scan of primitive quadruples, maxspeed in (32, {T}]")
    capA, capB = Fr(1,12), Fr(5,42)
    t0 = time.time()
    checked = skipped = 0
    worstA = (Fr(0), None); worstB = (Fr(0), None)
    for d in range(33, T+1):
        for abc in combinations(range(1, d), 3):
            v = list(abc) + [d]
            g = 0
            for x in v: g = gcd(g, x)
            if g != 1: continue
            # prune: analytic upper bound
            U = Fr(48, 2401)
            for i in range(4):
                for j in range(i+1, 4):
                    p, q = v[i], v[j]
                    g2 = gcd(p, q); p, q = p//g2, q//g2
                    if p*q <= 40:
                        U += pair_dev_max(min(p,q), max(p,q)) * Fr(2,49) * 2
                    else:
                        U += Fr(12, 7*p*q) * Fr(2,49) * 2
            U += Fr(1, 20)  # generous mass+transient allowance
            if U < capA:
                skipped += 1
                continue
            A = capA_value(v); B = capB_value(v)
            checked += 1
            if A > worstA[0]: worstA = (A, tuple(v))
            if B > worstB[0]: worstB = (B, tuple(v))
            if A > capA or B > capB:
                print(f"  *** EXCEEDANCE {v}: A={A} B={B}")
        if time.time() - t0 > budget:
            print(f"  [budget: reached d={d} of {T}]")
            break
    print(f"  checked={checked} skipped-by-prune={skipped}  "
          f"worstA={float(worstA[0]):.5f}@{worstA[1]}  worstB={float(worstB[0]):.5f}@{worstB[1]}  "
          f"caps A={float(capA):.5f} B={float(capB):.5f}  [{time.time()-t0:.0f}s]")

# ---------------- PART 2: the Z9 cyclic 2-page book ----------------

def crossing_pairs(n):
    """chords {a,b},{c,d} on a cyclic spine cross iff interleaved."""
    def cross(e, f):
        a, b = sorted(e); c, d = sorted(f)
        return (a < c < b < d) or (c < a < d < b)
    return cross

def z_book(n):
    """Parallel classes mod n (odd n): class s = chords {a,b}, a+b ≡ s. Class-crossing
    matrix + min crossings over 2-colorings of classes vs Guy Z(n)."""
    cross = crossing_pairs(n)
    classes = {s: [] for s in range(n)}
    for a in range(n):
        for b in range(a+1, n):
            classes[(a+b) % n].append((a, b))
    # within-class crossings (should be 0: parallel chords don't interleave on the cycle)
    within = {s: sum(1 for e, f in combinations(classes[s], 2) if cross(e, f)) for s in range(n)}
    X = [[0]*n for _ in range(n)]
    for s in range(n):
        for t in range(s+1, n):
            c = sum(1 for e in classes[s] for f in classes[t] if cross(e, f))
            X[s][t] = X[t][s] = c
    # circulant check
    circ = all(X[s][t] == X[0][(t-s) % n] for s in range(n) for t in range(n) if s != t)
    best, bestmask = None, None
    for mask in range(1 << (n-1)):   # fix class 0 on page 0
        pages = [(mask >> s) & 1 if s else 0 for s in range(n)]
        tot = sum(X[s][t] for s in range(n) for t in range(s+1, n) if pages[s] == pages[t])
        tot += sum(within.values())
        if best is None or tot < best:
            best, bestmask = tot, pages
    Zg = ((n//2)*((n-1)//2)*((n-2)//2)*((n-3)//2))//4
    print(f"  n={n}: within-class crossings = {sum(within.values())} (0 iff parallel⟹noncross); "
          f"circulant={circ}; row X[0] = {X[0]}; MIN over class-2-colorings = {best} "
          f"vs Z(n) = {Zg} {'== OPTIMAL' if best == Zg else ('(gap ' + str(best - Zg) + ')')}; "
          f"best pages = {bestmask}")

def part2():
    print("\nPART 2: THE CYCLIC 2-PAGE BOOK — page assignment by parallel classes (spine = Z_n)")
    for n in [5, 7, 9, 11, 13]:
        z_book(n)
        sys.stdout.flush()

if __name__ == "__main__":
    t0 = time.time()
    tab, worst = part1_tables()
    sys.stdout.flush()
    part2()
    sys.stdout.flush()
    part1_scan()
    print(f"[total {time.time()-t0:.1f}s]")
