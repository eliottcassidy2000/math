#!/usr/bin/env python3
"""collatz_procgen_20260922_endgame_hostile.py -- exact hostility prover for positive dyadic rationals in the
backward E-game, applied to the census, the generation-1 families and the Psi-candidates.

x in Bad_inf  <=>  no 3-adically legal reverse path from x ever has multiplier R = 2^K/3^s < 1.
Every node of x's reverse tree is dyadic; along a path z = R (x - beta), beta = sum_t 3^(t-1)/2^(K_t).

Proved ingredients (note, section 4):
  Lemma V   : a node z > 0 with z >= x has R > 1.
  Lemma S_b*: a positive integer node y with R > 3y/4 (S1) or 3R >= 2y+1 (S2) is SAFE (every continuation
              keeps R > 1).  (S1) holds at every positive integer node when x <= 4/3.
  (I)       : the integer children y_k of any node z (moves k >= f, z = N/2^f) satisfy (S2) as soon as
              2^k (3R - 2z) >= 1; for x < 3/2 we have 3R > 2z at every node.
The search below is therefore a complete decision procedure on dyadics x < 3/2 whose trees have no negative
nodes (negative nodes are searched for explicit descents); for x > 3/2 the root has infinitely many unsafe
integer children and the procedure can only find descents (reported 'undecided' otherwise).
Usage: python3 collatz_procgen_20260922_endgame_hostile.py [FAMILY_IMAX] [EMAX_CANDIDATES] [SECTIONS, default GHIJ]
"""
import sys, time, math
from fractions import Fraction as F
from collections import Counter

HALF = F(1, 2)
def K0(s): return (3 ** (s + 1)).bit_length() - 1
def Kstar(s): return 0 if s == 0 else (2 if s == 1 else K0(s))
def rho(half, k): return F(2 ** (Kstar(k - 1) + (1 if half else 0)), 3 ** k)
def c_of(s): return F(2 ** K0(s), 3 ** s)
def res3(x, mod):
    x = F(x); return (x.numerator * pow(x.denominator, -1, mod)) % mod
def v3(x):
    x = F(x); n, d = x.numerator, x.denominator; v = 0
    if n == 0: return 10 ** 9
    while n % 3 == 0: n //= 3; v += 1
    while d % 3 == 0: d //= 3; v -= 1
    return v
def psi(x):
    x = F(x); h = F(1) if res3(x, 3) == 1 else HALF
    if x == h: return None
    k = v3(x - h); r = rho(h == HALF, k)
    return (r * (x - h), h, k, r)
def itinerary(x, maxsteps=200):
    x = F(x); P = F(1); D = 0; steps = []
    for _ in range(maxsteps):
        r = psi(x)
        if r is None: return steps, ('1' if x == 1 else '1/2')
        y, h, k, rr = r; P *= rr; D += k
        steps.append((h, k, rr, P, D, y)); x = y
    return steps, None
def fmt(x):
    x = F(x); return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"

def split2(x):
    x = F(x); d = x.denominator; f = 0
    while d % 2 == 0: d //= 2; f += 1
    assert d == 1, "not dyadic"
    return x.numerator, f

def children(N, f, kmax):
    r9 = (N * pow(pow(2, f, 9), -1, 9)) % 9
    k = 0 if r9 % 3 == 1 else 1
    p2 = pow(2, k, 9)
    while k <= kmax:
        if (p2 * r9) % 9 in (4, 7):
            if k >= f:
                cn = (N * (1 << (k - f)) - 1) // 3; cf = 0
            else:
                cn = (N - (1 << (f - k))) // 3; cf = f - k
            yield k, cn, cf
        k += 2; p2 = (p2 * 4) % 9

def safe_int(y, K, s):
    R2, R3 = 1 << K, 3 ** s
    return (4 * R2 > 3 * R3 * y) or (3 * R2 >= (2 * y + 1) * R3)

def kthreshold(N, f, K, s):
    num = 3 * (1 << K) * (1 << f) - 2 * N * 3 ** s
    if num <= 0: return None
    den = 3 ** s * (1 << f); k = f
    while (1 << k) * num < den: k += 1
    return k

def neg_descent(N, f, K, s, limit):
    stack = [(N, f, K, s)]; cnt = 0
    while stack:
        N, f, K, s = stack.pop(); cnt += 1
        if cnt > limit: return None
        if (1 << K) < 3 ** s: return (s, K, f"{N}/2^{f}")
        for k, cn, cf in children(N, f, f + 8): stack.append((cn, cf, K + k, s + 1))
    return None

def prove(x, nodelimit=3_000_000, neg_limit=300_000, int_nodes=2_000_000, kcap_unsafe=None):
    """('hostile', info) complete proof | ('descends', cert) | ('undecided', why).
    kcap_unsafe: if set, nodes with 3R <= 2z are expanded with moves k <= kcap_unsafe (search only)."""
    x = F(x); assert x > 0
    N0, f0 = split2(x)
    if N0 % 3 == 0: return ('not a unit', None)
    if f0 == 0: return ('integer start', None)
    stack = [(N0, f0, 0, 0)]; nodes = 0; unsafe = 0; truncated = False
    while stack:
        N, f, K, s = stack.pop(); nodes += 1
        if nodes > nodelimit: return ('undecided', ('node limit', nodes))
        if s >= 1 and (1 << K) < 3 ** s: return ('descends', ('prefix', s, K, f"{N}/2^{f}"))
        if N < 0:
            r = neg_descent(N, f, K, s, neg_limit)
            if r is not None: return ('descends', ('negative', r))
            return ('undecided', ('negative node without certificate', N, f, K, s))
        if f == 0:
            if safe_int(N, K, s): continue
            unsafe += 1
            if unsafe > int_nodes: return ('undecided', ('unsafe integer nodes', unsafe))
        kt = kthreshold(N, f, K, s)
        if kt is None:
            if kcap_unsafe is None: return ('undecided', ('3R <= 2z', f"{N}/2^{f}", K, s))
            kt = kcap_unsafe; truncated = True
        for k, cn, cf in children(N, f, max(kt, f + 2)): stack.append((cn, cf, K + k, s + 1))
    if truncated: return ('undecided', ('no descent with capped moves', nodes))
    return ('hostile', ('nodes', nodes, 'unsafe integer nodes explored', unsafe))

def descent_search(x, maxdepth=20, nodelimit=400_000, kcap=40):
    x = F(x); N0, f0 = split2(x)
    for D in range(1, maxdepth + 1):
        stack = [(N0, f0, 0, 0)]; nodes = 0
        while stack:
            N, f, K, s = stack.pop(); nodes += 1
            if nodes > nodelimit: break
            if s >= 1 and (1 << K) < 3 ** s: return (s, K, f"{N}/2^{f}")
            if s == D or (1 << K) >= 3 ** D: continue
            for k, cn, cf in children(N, f, f + kcap):
                if (1 << (K + k)) >= 3 ** D: break
                stack.append((cn, cf, K + k, s + 1))
    return None

CENSUS_FILE_NOTE = "census data: collatz_procgen_20260922_endgame_chain.py (CENSUS)"

def main():
    IMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 22
    EMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    SECT = sys.argv[3] if len(sys.argv) > 3 else "GHIJ"
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from collatz_procgen_20260922_endgame_chain import census_points
    pts = census_points()
    if "G" in SECT: section_G(pts)
    if "H" in SECT: section_H(IMAX)
    if "I" in SECT: section_I(pts, EMAX)
    if "J" in SECT: section_J(pts)

def section_G(pts):
    print("=" * 100)
    print("G. exact prover on the 98-point census (97 non-integer points)")
    t = time.time(); res = {}
    for x in pts:
        if x == 1: continue
        res[x] = prove(x, nodelimit=12_000_000)
    c = Counter(r[0] for r in res.values())
    print(f"  results: {dict(c)}   ({time.time()-t:.0f}s)")
    for x, r in res.items():
        if r[0] != 'hostile': print(f"    {fmt(x)} = {float(x):.5f}: {r[0]} {r[1]}")
    hv = [float(x) for x, r in res.items() if r[0] == 'hostile']
    print(f"  PROVED hostile: {len(hv)} points, real values in [{min(hv):.4f}, {max(hv):.4f}]")
    sys.stdout.flush()

def section_H(IMAX):
    print("=" * 100)
    print(f"H. generation-1 families (Theorem F_b) and the T^-1(2) points, i = 3..{IMAX}")
    print("   x+(i) = 1/2 + 3^i/2^(K0(i-1)+2) = T_i^-1(1/2);  x-(i) = 1/2 + 3^i/2^(K0(i-1)+1) = T_i^-1(1);"
          "  t2(i) = 1/2 + 3^i/2^K0(i-1) = T_i^-1(2)")
    for i in range(3, IMAX + 1):
        K = K0(i - 1); line = f"  i={i:2d} c(i-1)={float(c_of(i-1)):.4f}:"
        for name, x in (("x+", HALF + F(3 ** i, 2 ** (K + 2))), ("x-", HALF + F(3 ** i, 2 ** (K + 1))),
                        ("t2", HALF + F(3 ** i, 2 ** K))):
            r = prove(x, nodelimit=6_000_000)
            if r[0] == 'undecided' and x > F(3, 2):
                d = descent_search(x, maxdepth=18)
                if d: r = ('descends', d)
            tag = {'hostile': 'H', 'descends': 'D', 'undecided': '?'}[r[0]]
            line += f"  {name}={float(x):.4f}:{tag}"
        print(line); sys.stdout.flush()

def section_I(pts, EMAX):
    print("=" * 100)
    print(f"I. Psi-candidates (Psi-preimages of 1, 1/2 with every prefix multiplier > 1), denominators <= 2^{EMAX}")
    def denexp(x):
        d = F(x).denominator; e = 0
        while d % 2 == 0: d //= 2; e += 1
        return e
    t = time.time()
    stack = [(F(1), F(1), F(0)), (HALF, F(1), F(0))]; valid = {}; nodes = 0
    while stack:
        y, Q, M = stack.pop(); nodes += 1; ey = denexp(y)
        for half in (False, True):
            h = HALF if half else F(1)
            for k in range(1, 400):
                a = Kstar(k - 1) + (1 if half else 0)
                if ey + a > EMAX: break
                r = rho(half, k); x = h + y / r
                Qx = r * Q; Mx = max(Q, M, F(1))
                if Qx * 2 ** ((EMAX - denexp(x)) // 5) <= Mx: continue    # can never become a valid start
                if Qx > Mx and F(3, 10) <= x <= 2: valid[x] = Qx
                stack.append((x, Qx, Mx))
    cand = sorted(valid)
    cset = set(pts)
    inc = [x for x in cand if x in cset]
    notin = [x for x in cand if x not in cset]
    missing = [x for x in pts if denexp(x) <= EMAX and x not in valid and x not in (1, HALF)]
    print(f"  {len(cand)} candidates ({nodes} backward nodes, {time.time()-t:.1f}s); census points with e <= {EMAX}"
          f" among them: {len(inc)}; census points with e <= {EMAX} missing: {len(missing)}")
    t = time.time(); cc = Counter(); und = []
    for x in notin:
        r = prove(x, nodelimit=300_000)
        if r[0] == 'undecided':
            d = descent_search(x, maxdepth=20)
            if d: r = ('descends', d)
        cc[r[0]] += 1
        if r[0] != 'descends': und.append((fmt(x), round(float(x), 4), r[0]))
    print(f"  candidates outside the census: {len(notin)}: {dict(cc)} ({time.time()-t:.0f}s); not descended: {und[:6]}")
    byv = lambda L: (sum(1 for x in L if x <= 1), sum(1 for x in L if 1 < x <= F(3, 2)), sum(1 for x in L if x > F(3, 2)))
    print(f"  census (e<={EMAX}) by value [<=1, (1,3/2], >3/2]: {byv(inc)};  non-census candidates: {byv(notin)}")
    sys.stdout.flush()

def capped_descent(x, depth, kcap, nodelimit):
    """depth-limited DFS for a descent from x, expanding every node (safe integer nodes pruned by S1/S2),
    moves k <= max(threshold, kcap) at nodes where 3R <= 2z.  Returns certificate, 'exhausted' or 'limit'."""
    N0, f0 = split2(x); stack = [(N0, f0, 0, 0)]; nodes = 0
    while stack:
        N, f, K, s = stack.pop(); nodes += 1
        if nodes > nodelimit: return 'limit', nodes
        if s >= 1 and (1 << K) < 3 ** s: return (s, K, f"{N}/2^{f}"), nodes
        if s >= depth or N < 0: continue
        if f == 0 and safe_int(N, K, s): continue
        kt = kthreshold(N, f, K, s)
        kk = kcap if kt is None else max(kt, f + 2)
        for k, cn, cf in children(N, f, kk): stack.append((cn, cf, K + k, s + 1))
    return 'exhausted', nodes

def prove_root_capped(x, kroot, nodelimit=2_000_000):
    """explore x's whole tree except the root's children with move k > kroot; nodes where 3R <= 2z are
    recorded as unresolved instead of expanded.  Returns (verdict, nodes, unresolved nodes)."""
    N0, f0 = split2(x); stack = [(N0, f0, 0, 0)]; nodes = 0; unres = []
    while stack:
        N, f, K, s = stack.pop(); nodes += 1
        if nodes > nodelimit: return ('node limit', nodes, unres)
        if s >= 1 and (1 << K) < 3 ** s: return ('DESCENT', (s, K, f"{N}/2^{f}"), unres)
        if N < 0: return ('negative node', (N, f, K, s), unres)
        if f == 0 and s >= 1 and safe_int(N, K, s): continue
        kt = kthreshold(N, f, K, s)
        if kt is None:
            if s == 0: kt = kroot
            else: unres.append((f"{N}/2^{f}", K, s)); continue
        for k, cn, cf in children(N, f, max(kt, f + 2)): stack.append((cn, cf, K + k, s + 1))
    return ('explored', nodes, unres)

def section_J(pts):
    print("=" * 100)
    print("J. census points above 3/2: depth-limited search for descents (the root has infinitely many unsafe children)")
    for x in [x for x in pts if x > F(3, 2)]:
        st, term = itinerary(x)
        out = []
        for depth, kcap in ((12, 60), (16, 44), (20, 34)):
            r, n = capped_descent(x, depth, kcap, 400_000)
            out.append(f"depth<={depth},k<={kcap}: {r if isinstance(r, str) else 'DESCENT '+str(r)} ({n} nodes)")
        print(f"  {fmt(x)} = {float(x):.5f}; itinerary {' '.join(('1/2' if s_[0]==HALF else '1')+'@'+str(s_[1]) for s_ in st)} -> {term}")
        for o in out: print("     " + o)
        v, n, unres = prove_root_capped(x, 60)
        depths = sorted(set(u[2] for u in unres))
        print(f"     rest of the tree (root children with k <= 60): {v}, {n} nodes; unresolved nodes (3R <= 2z): "
              f"{len(unres)}, all at depth(s) {depths}, e.g. {unres[:2]}")
        sys.stdout.flush()

if __name__ == '__main__':
    main()
