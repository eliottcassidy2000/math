#!/usr/bin/env python3
"""Independent verifier, part (A)1 of THM-5141 (draft).

Claim (A)1: no nonzero delta in {-2..2}^w has  sum delta_i 2^i = 0  and  sum delta_i 2^(w-1-i) = 0.

Four independent methods (all exact integer arithmetic):
  M1  brute force over all (2B+1)^w digit strings (odometer), B = 2 for w <= WB2, and B = 3, 4 for small w
      (B = 3, 4 are the NON-VACUOUS controls: nonzero solutions exist there, e.g. (2,-3,-3,2)).
  M2  meet in the middle: w = h1 + h2, enumerate (2B+1)^h half strings only.  Exhaustive for the full length w:
      D = D_lo + 2^h1 D_hi = 0 forces D_lo = 2^h1 A with |A| <= B-1 and D_hi = -A;  D* = 2^h2 D*_lo + D*_hi = 0
      forces D*_hi = 2^h2 A' with |A'| <= B-1 and D*_lo = -A'.  So the number of solutions is
      sum_{p,q} T1_h1[p,q] * T2_h2[-p,-q]  with
        T1_h[p,q] = #{e in [-B,B]^h : D(e) = p 2^h, D*(e) = q},   T2_h[p,q] = #{e : D(e) = p, D*(e) = q 2^h}.
  M3  finite automaton, valid for EVERY w: states (c, V), c = carry of D from the low end
      (c_0 = 0, c_(i+1) = (c_i + delta_i)/2, must be an integer, final 0), V = Horner value of D* from the top
      (V_0 = 0, V_(i+1) = 2 V_i + delta_i, must satisfy |V| <= B-1 to be able to return to 0, final 0).
      Number of solutions of length w = number of walks (0,0) -> (0,0) of length w.
  M4  the lattice of the theorem's proof: delta = (2z^2 - 5z + 2) p.  Depth-first search over integer p with
      p_0 = 1 (sign and shift normalised) keeping |delta_i| <= 2; asserts |u_i| = |2p_i - p_(i-1)| >= 2 at every node
      and that the closing state (p_(i-1), p_i) = (0, 0) is never reached; also checks the Gauss-lemma
      factorisation on every solution of the B = 3, 4 controls (exact division by 2z^2 - 5z + 2 in Z[z]).
Usage: python a1_injectivity.py [WB2 (default 12)] [HMAX (default 13)] [DFSDEPTH (default 22)]
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
import sys
import time
import numpy as np
from numba import njit, prange


# ---------------------------------------------------------------- M1 brute force
@njit(cache=False)
def brute_count(w, B):
    base = 2 * B + 1
    total = 1
    for _ in range(w):
        total *= base
    digits = np.empty(w, np.int64)
    D = 0
    Ds = 0
    for i in range(w):
        digits[i] = -B
        D += -B * (1 << i)
        Ds += -B * (1 << (w - 1 - i))
    cnt = 0
    for _ in range(total):
        if D == 0 and Ds == 0:
            cnt += 1
        i = 0
        while i < w:
            if digits[i] < B:
                digits[i] += 1
                D += (1 << i)
                Ds += (1 << (w - 1 - i))
                break
            else:
                digits[i] = -B
                D -= 2 * B * (1 << i)
                Ds -= 2 * B * (1 << (w - 1 - i))
                i += 1
    return cnt


def brute_solutions_py(w, B):
    """pure python, independent of the numba odometer: list of all solutions (small sizes only)."""
    import itertools
    sols = []
    for d in itertools.product(range(-B, B + 1), repeat=w):
        if sum(d[i] << i for i in range(w)) == 0 and sum(d[i] << (w - 1 - i) for i in range(w)) == 0:
            sols.append(d)
    return sols


# ---------------------------------------------------------------- M2 meet in the middle
@njit(cache=False)
def _half_chunk(h, B, top, T1, T2):
    """strings e of length h with e_(h-1) = top; accumulate T1, T2 (arrays (2B-1) x (2B-1))."""
    hl = h - 1
    mask = (1 << h) - 1
    digits = np.empty(max(hl, 1), np.int64)
    D = top * (1 << (h - 1))
    Ds = top
    total = 1
    for i in range(hl):
        digits[i] = -B
        D += -B * (1 << i)
        Ds += -B * (1 << (h - 1 - i))
        total *= (2 * B + 1)
    off = B - 1
    for _ in range(total):
        # T1: D = p 2^h, D* = q with |q| <= B-1
        if (D & mask) == 0 and Ds >= -off and Ds <= off:
            p = D >> h
            T1[p + off, Ds + off] += 1
        if (Ds & mask) == 0 and D >= -off and D <= off:
            q = Ds >> h
            T2[D + off, q + off] += 1
        i = 0
        while i < hl:
            if digits[i] < B:
                digits[i] += 1
                D += (1 << i)
                Ds += (1 << (h - 1 - i))
                break
            else:
                digits[i] = -B
                D -= 2 * B * (1 << i)
                Ds -= 2 * B * (1 << (h - 1 - i))
                i += 1


@njit(parallel=True, cache=False)
def half_tables(h, B):
    nb = 2 * B + 1
    sz = 2 * B - 1
    R1 = np.zeros((nb, sz, sz), np.int64)
    R2 = np.zeros((nb, sz, sz), np.int64)
    for t in prange(nb):
        _half_chunk(h, B, t - B, R1[t], R2[t])
    T1 = np.zeros((sz, sz), np.int64)
    T2 = np.zeros((sz, sz), np.int64)
    for t in range(nb):
        T1 += R1[t]
        T2 += R2[t]
    return T1, T2


def mitm_count(tables, h1, h2, B):
    T1 = tables[h1][0]
    T2 = tables[h2][1]
    sz = 2 * B - 1
    tot = 0
    for p in range(sz):
        for q in range(sz):
            tot += int(T1[p, q]) * int(T2[sz - 1 - p, sz - 1 - q])
    return tot


# ---------------------------------------------------------------- M3 automaton (exact, python ints)
def automaton_counts(B, wmax):
    """number of solutions for each w = 0..wmax, digits in [-B, B]."""
    counts = []
    cur = {(0, 0): 1}
    maxc = 0
    for w in range(0, wmax + 1):
        counts.append(cur.get((0, 0), 0))
        nxt = {}
        for (c, V), m in cur.items():
            for d in range(-B, B + 1):
                if (c + d) % 2:
                    continue
                c2 = (c + d) // 2
                V2 = 2 * V + d
                if abs(V2) > B - 1:
                    continue
                maxc = max(maxc, abs(c2))
                nxt[(c2, V2)] = nxt.get((c2, V2), 0) + m
        cur = nxt
    return counts, maxc


def automaton_edges(B):
    """all transitions among states reachable from (0,0)."""
    seen = {(0, 0)}
    stack = [(0, 0)]
    edges = []
    while stack:
        c, V = stack.pop()
        for d in range(-B, B + 1):
            if (c + d) % 2:
                continue
            c2 = (c + d) // 2
            V2 = 2 * V + d
            if abs(V2) > B - 1:
                continue
            edges.append(((c, V), d, (c2, V2)))
            if (c2, V2) not in seen:
                seen.add((c2, V2))
                stack.append((c2, V2))
    return seen, edges


# ---------------------------------------------------------------- M4 lattice / the proof's recursion
def poly_mul_q(p):
    """(2z^2 - 5z + 2) * p, coefficient lists (lowest first)."""
    n = len(p)
    out = [0] * (n + 2)
    for i, c in enumerate(p):
        out[i] += 2 * c
        out[i + 1] -= 5 * c
        out[i + 2] += 2 * c
    return out


def exact_div_q(delta):
    """divide delta (lowest first) by 2z^2 - 5z + 2 in Q[z]; return (integer quotient list or None, remainder zero?)."""
    from fractions import Fraction
    d = [Fraction(c) for c in delta]
    n = len(d)
    if n < 3:
        return ([] if all(c == 0 for c in d) else None), all(c == 0 for c in d)
    q = [Fraction(0)] * (n - 2)
    for i in range(n - 1, 1, -1):          # leading coefficient 2 at z^2
        q[i - 2] = d[i] / 2
        d[i] -= 2 * q[i - 2]
        d[i - 1] += 5 * q[i - 2]
        d[i - 2] -= 2 * q[i - 2]
    rem_zero = all(c == 0 for c in d)
    integral = all(c.denominator == 1 for c in q)
    return ([int(c) for c in q] if integral else None), rem_zero


@njit(cache=False)
def lattice_dfs(depth):
    """DFS over p_0 = 1, p_1, ..., with delta_i = 2p_i - 5p_(i-1) + 2p_(i-2) in [-2,2] for all i up to depth.
    returns (nodes, min |u_i| seen, closures) where a closure is a node with (p_(i-1), p_i) = (0,0)."""
    # explicit stack of (level, p_prev, p_cur)
    cap = 64 * (depth + 4)
    st_l = np.empty(cap, np.int64)
    st_a = np.empty(cap, np.int64)
    st_b = np.empty(cap, np.int64)
    top = 0
    # level 0: p_(-1) = 0, p_0 = 1 (delta_0 = 2)
    st_l[0] = 0
    st_a[0] = 0
    st_b[0] = 1
    top = 1
    nodes = 0
    minu = 1 << 60
    closures = 0
    maxp = 0
    while top > 0:
        top -= 1
        lev = st_l[top]
        a = st_a[top]
        b = st_b[top]
        nodes += 1
        u = 2 * b - a
        if u < 0:
            u = -u
        if u < minu:
            minu = u
        if a == 0 and b == 0:
            closures += 1
            continue
        if (b if b > 0 else -b) > maxp:
            maxp = b if b > 0 else -b
        if lev == depth:
            continue
        t = 5 * b - 2 * a            # need |2x - t| <= 2
        # x in [(t-2)/2, (t+2)/2]
        lo = (t - 2 + 1) // 2 if (t - 2) % 2 != 0 else (t - 2) // 2
        hi = (t + 2) // 2
        for x in range(lo, hi + 1):
            dlt = 2 * x - t
            if dlt >= -2 and dlt <= 2:
                st_l[top] = lev + 1
                st_a[top] = b
                st_b[top] = x
                top += 1
    return nodes, minu, closures, maxp


def main():
    WB2 = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    HMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 13
    DFSD = int(sys.argv[3]) if len(sys.argv) > 3 else 22
    ok = True
    print("=== M1 brute force: number of delta in [-B,B]^w with D = 0 and D* = 0 ===")
    brute = {}
    for B, wmax in ((2, WB2), (3, 9), (4, 8)):
        for w in range(1, wmax + 1):
            t0 = time.time()
            c = brute_count(w, B)
            brute[(B, w)] = c
            print(f"  B = {B}, w = {w:2d}: {(2*B+1)**w:>12d} strings, solutions = {c}   [{time.time()-t0:.1f}s]", flush=True)
            if B == 2 and c != 1:
                ok = False
                print("  *** (A)1 FAILS")
    # pure-python cross check of the odometer (small)
    for B, w in ((2, 7), (3, 6), (4, 5)):
        s = brute_solutions_py(w, B)
        same = len(s) == brute[(B, w)]
        ok &= same
        print(f"  pure python cross-check B = {B}, w = {w}: {len(s)} solutions (numba: {brute[(B, w)]}) match: {same}")
    print("  example nonzero solutions for B = 3, w = 4:", [s for s in brute_solutions_py(4, 3) if any(s)])

    print("=== M3 automaton (valid for every w) ===")
    for B in (2, 3, 4):
        counts, maxc = automaton_counts(B, 400 if B == 2 else 40)
        seen, edges = automaton_edges(B)
        agree = all(counts[w] == brute[(B, w)] for (b, w) in brute if b == B)
        ok &= agree
        print(f"  B = {B}: reachable states {len(seen)}, transitions {len(edges)}, max |carry| {maxc}; "
              f"counts w = 1..12: {counts[1:13]}; agrees with brute force on all {sum(1 for (b, w) in brute if b == B)} brute cases: {agree}")
        if B == 2:
            allone = all(c == 1 for c in counts)
            ok &= allone
            print(f"     B = 2: count = 1 for every w = 0..400: {allone}")
            print(f"     B = 2: transitions out of the start state (0,0): {[e for e in edges if e[0] == (0, 0)]}")
            print(f"     B = 2: reachable set from (0,0): {sorted(seen)}")

    print("=== M2 meet in the middle (exhaustive for the full length w = h1 + h2) ===")
    for B, hmax in ((3, 5), (4, 4), (2, HMAX)):
        tables = {}
        for h in range(1, hmax + 1):
            t0 = time.time()
            tables[h] = half_tables(h, B)
            if B == 2:
                print(f"  B = 2: half tables h = {h:2d} ({5**h} strings)  T1 = {tables[h][0].tolist()}  T2 = {tables[h][1].tolist()}  [{time.time()-t0:.1f}s]", flush=True)
        res = []
        for w in range(2, 2 * hmax + 1):
            h1 = (w + 1) // 2
            h2 = w // 2
            c = mitm_count(tables, h1, h2, B)
            res.append((w, c))
            if (B, w) in brute and brute[(B, w)] != c:
                ok = False
                print(f"  *** MITM disagrees with brute force at B = {B}, w = {w}: {c} vs {brute[(B, w)]}")
            if B == 2 and c != 1:
                ok = False
                print(f"  *** (A)1 FAILS at w = {w}: {c} solutions")
        # also unbalanced splits as a self-consistency test
        if B != 2:
            cnts, _ = automaton_counts(B, 2 * hmax)
            agree = all(c == cnts[w] for w, c in res)
            ok &= agree
            print(f"  B = {B} (control): MITM counts w = 2..{2*hmax}: {[c for _, c in res]}; equal to automaton counts: {agree}")
        else:
            print(f"  B = 2: MITM solution counts for w = 2..{2*hmax}: {[c for _, c in res]}")
            unb = []
            for w in range(2, 2 * hmax + 1):
                for h1 in range(max(1, w - hmax), min(hmax, w - 1) + 1):
                    unb.append(mitm_count(tables, h1, w - h1, B))
            a1 = all(c == 1 for c in unb)
            ok &= a1
            print(f"  B = 2: all {len(unb)} splits (h1, h2) with h1, h2 <= {hmax}: every count = 1: {a1}")

    print("=== M4 the proof's lattice: delta = (2z^2 - 5z + 2) p ===")
    # (a) Gauss-lemma factorisation on the control solutions
    for B, w in ((3, 7), (4, 6), (3, 8)):
        sols = [s for s in brute_solutions_py(w, B)]
        nz = [s for s in sols if any(s)]
        good = True
        for s in nz:
            q, rz = exact_div_q(list(s))
            if q is None or not rz or poly_mul_q(q) != list(s):
                good = False
        ok &= good
        print(f"  B = {B}, w = {w}: {len(nz)} nonzero solutions, every one = (2z^2-5z+2) * (integer polynomial): {good}")
    # (b) recursion identities on random integer p
    rng = np.random.default_rng(12345)
    good = True
    for _ in range(20000):
        n = int(rng.integers(1, 12))
        p = [int(v) for v in rng.integers(-6, 7, size=n)]
        if not any(p):
            continue
        m = next(i for i, c in enumerate(p) if c)
        d = max(i for i, c in enumerate(p) if c)
        delta = poly_mul_q(p)
        P = lambda i: p[i] if 0 <= i < n else 0
        Dl = lambda i: delta[i] if 0 <= i < len(delta) else 0
        u = lambda i: 2 * P(i) - P(i - 1)
        for i in range(-2, n + 5):
            if Dl(i) != 2 * P(i) - 5 * P(i - 1) + 2 * P(i - 2):
                good = False
            if u(i) != 2 * u(i - 1) + Dl(i):
                good = False
        if abs(u(m)) != 2 * abs(p[m]) or u(d + 2) != 0 or u(d + 1) != -p[d]:
            good = False
        # value checks
        from fractions import Fraction
        if sum(c * 2 ** i for i, c in enumerate(delta)) != 0 or sum(Fraction(c, 2 ** i) for i, c in enumerate(delta)) != 0:
            good = False
        if max(abs(c) for c in delta) < 3:
            good = False
            print("  *** nonzero p with max |delta_i| <= 2:", p, delta)
    ok &= good
    print(f"  20000 random p: delta_i = 2p_i - 5p_(i-1) + 2p_(i-2) = u_i - 2u_(i-1), |u_m| = 2|p_m|, u_(d+1) = -p_d, u_(d+2) = 0, "
          f"delta(2) = delta(1/2) = 0, max|delta_i| >= 3: {good}")
    # (c) exhaustive small p: all p in [-3,3]^n, n <= 7, p != 0: max |coefficient of (2z^2-5z+2)p| >= 3
    import itertools
    good = True
    cnt = 0
    minmax = 99
    for n in range(1, 8):
        for p in itertools.product(range(-3, 4), repeat=n):
            if p[0] == 0 or p[-1] == 0:
                continue
            cnt += 1
            mx = max(abs(c) for c in poly_mul_q(list(p)))
            minmax = min(minmax, mx)
            if mx < 3:
                good = False
    ok &= good
    print(f"  exhaustive: {cnt} polynomials p (coefficients in [-3,3], length <= 7, p_0 p_top != 0): min over p of max|delta_i| = {minmax} (>= 3: {good})")
    # (d) DFS over the lattice with |delta_i| <= 2
    t0 = time.time()
    for depth in (6, 10, 14, 18, DFSD):
        nodes, minu, closures, maxp = lattice_dfs(depth)
        g = (minu >= 2 and closures == 0)
        ok &= g
        print(f"  lattice DFS to depth {depth:2d} (delta of length up to {depth + 1}, p_0 = 1): nodes = {nodes}, min |u_i| = {minu}, "
              f"closures (p returns to 0 with |delta| <= 2) = {closures}, max |p_i| = {maxp}  -> {'OK' if g else 'FAIL'}  [{time.time()-t0:.1f}s]", flush=True)

    print("ALL CHECKS PASSED" if ok else "*** SOME CHECK FAILED ***")


if __name__ == "__main__":
    main()
