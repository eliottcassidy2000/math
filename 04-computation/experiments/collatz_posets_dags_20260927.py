#!/usr/bin/env python3
"""Posets and DAGs of a Collatz orbit: value-time poset, excursion forest, spine,
descent tree with its two-place (2-adic source / 3-adic landing) structure,
spine blocks = THM-4495's positive min-ending words, the Young-lattice form of
THM-4503's cells, and the Sturmian element of E_inf.

Session: opus, collatz-poset-dag-20260927.  Exact integer arithmetic throughout;
floats appear only in display columns.

Run:  python 04-computation/experiments/collatz_posets_dags_20260927.py
"""
from __future__ import annotations

import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction

# ----------------------------------------------------------------------------
# exact height comparison: h(o, j) = o*log2(3) - j
# ----------------------------------------------------------------------------

def height_lt(o1: int, j1: int, o2: int, j2: int) -> bool:
    """Exact test h(o1,j1) < h(o2,j2)  <=>  (o1-o2) log2 3 < j1 - j2."""
    do, dj = o1 - o2, j1 - j2
    if do == 0 and dj == 0:
        return False
    # 3^do < 2^dj   <=>   3^max(do,0) * 2^max(-dj,0) < 2^max(dj,0) * 3^max(-do,0)
    lhs = 3 ** max(do, 0) * 2 ** max(-dj, 0)
    rhs = 2 ** max(dj, 0) * 3 ** max(-do, 0)
    if lhs == rhs:
        raise ValueError("tie in heights: impossible for distinct (o,j) by irrationality")
    return lhs < rhs


def height_float(o: int, j: int) -> float:
    return o * math.log2(3) - j


# ----------------------------------------------------------------------------
# Collatz maps
# ----------------------------------------------------------------------------

def T(n: int) -> int:
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def t_orbit_to_one(n: int, cap: int = 10 ** 7) -> list[int]:
    orb = [n]
    while orb[-1] != 1 and len(orb) < cap:
        orb.append(T(orb[-1]))
    return orb


def word_heights(orb: list[int]) -> list[tuple[int, int]]:
    """(o_j, j) for the T-orbit prefix: o_j = number of odd entries among the first j."""
    out = [(0, 0)]
    o = 0
    for j in range(1, len(orb)):
        if orb[j - 1] % 2 == 1:
            o += 1
        out.append((o, j))
    return out


# ----------------------------------------------------------------------------
# value-time poset, records, future minima, excursion forest
# ----------------------------------------------------------------------------

def strict_lower_records(vals: list[int]) -> list[int]:
    rec, cur = [], None
    for i, v in enumerate(vals):
        if cur is None or v < cur:
            rec.append(i)
            cur = v
    return rec


def strict_upper_records(vals: list[int]) -> list[int]:
    rec, cur = [], None
    for i, v in enumerate(vals):
        if cur is None or v > cur:
            rec.append(i)
            cur = v
    return rec


def value_future_minima(vals: list[int]) -> list[int]:
    """j with vals[j] < vals[k] for all k > j (within the given window)."""
    out, cur = [], None
    for j in range(len(vals) - 1, -1, -1):
        if cur is None or vals[j] < cur:
            out.append(j)
            cur = vals[j]
    return out[::-1]


def coefficient_future_minima(hs: list[tuple[int, int]]) -> list[int]:
    """j with h_j < h_k for all k > j (window), exact comparison."""
    out, cur = [], None
    for j in range(len(hs) - 1, -1, -1):
        if cur is None or height_lt(hs[j][0], hs[j][1], cur[0], cur[1]):
            out.append(j)
            cur = hs[j]
    return out[::-1]


def excursion_parent(hs: list[tuple[int, int]]) -> list[int]:
    """Excursion order: i is the parent of j iff i < j, h_k > h_i for all k in (i, j],
    and i is the largest such.  parent[j] = -1 for roots.  Stack algorithm."""
    parent = [-1] * len(hs)
    stack: list[int] = []
    for j in range(len(hs)):
        while stack and not height_lt(hs[stack[-1]][0], hs[stack[-1]][1], hs[j][0], hs[j][1]):
            stack.pop()
        parent[j] = stack[-1] if stack else -1
        stack.append(j)
    return parent


def check_value_time_poset(vals: list[int]) -> dict:
    """Value-time poset: i < j and vals[i] > vals[j].  Verify minimal = strict upper
    records (leaders), maximal = strict value future minima, and that the
    cover-free relation has dimension <= 2 (it is the intersection of two orders,
    so this is by construction; we verify the extremal-element identification)."""
    n = len(vals)
    minimal = [i for i in range(n) if not any(vals[h] > vals[i] for h in range(i))]
    maximal = [j for j in range(n) if not any(vals[k] < vals[j] for k in range(j + 1, n))]
    assert minimal == strict_upper_records(vals)
    assert maximal == value_future_minima(vals)
    return {"n": n, "minimal": len(minimal), "maximal": len(maximal)}


# ----------------------------------------------------------------------------
# no-descent counts, spine blocks, generating functions
# ----------------------------------------------------------------------------

def no_descent_counts(K: int) -> list[int]:
    """W_k = # words of length k with 3^{o_j} > 2^j for all 1 <= j <= k (ballot DP)."""
    W = [1]
    layer = {0: 1}  # o -> count, at length 0
    for j in range(1, K + 1):
        new: dict[int, int] = defaultdict(int)
        for o, c in layer.items():
            for bit in (0, 1):
                o2 = o + bit
                if 3 ** o2 > 2 ** j:
                    new[o2] += c
        layer = dict(new)
        W.append(sum(layer.values()))
    return W


def binomial_tails(K: int) -> list[int]:
    B = [0]
    for n in range(1, K + 1):
        B.append(sum(math.comb(n, j) for j in range(n + 1) if 3 ** j > 2 ** n))
    return B


def series_inverse(a: list[int]) -> list[Fraction]:
    """coefficients of 1/A(t) for A = sum a_k t^k with a_0 = 1."""
    K = len(a) - 1
    inv = [Fraction(1)]
    for k in range(1, K + 1):
        s = sum(Fraction(a[i]) * inv[k - i] for i in range(1, k + 1))
        inv.append(-s)
    return inv


def series_exp(c: list[Fraction]) -> list[Fraction]:
    """exp of a series with c_0 = 0: e' = c' e."""
    K = len(c) - 1
    e = [Fraction(1)]
    for k in range(1, K + 1):
        s = sum(Fraction(i) * c[i] * e[k - i] for i in range(1, k + 1))
        e.append(s / k)
    return e


def enumerate_words(k: int):
    for x in range(1 << k):
        yield [(x >> (k - 1 - i)) & 1 for i in range(k)]


def prefix_heights(w: list[int]) -> list[tuple[int, int]]:
    hs = [(0, 0)]
    o = 0
    for j, b in enumerate(w, start=1):
        o += b
        hs.append((o, j))
    return hs


def is_no_descent(w: list[int]) -> bool:
    o = 0
    for j, b in enumerate(w, start=1):
        o += b
        if 3 ** o <= 2 ** j:
            return False
    return True


def is_spine_block(w: list[int]) -> bool:
    """all proper prefix heights (j >= 1) strictly above the final height, final height > 0."""
    hs = prefix_heights(w)
    L = len(w)
    o, j = hs[L]
    if not (3 ** o > 2 ** L):
        return False
    for i in range(1, L):
        if not height_lt(o, L, hs[i][0], hs[i][1]):
            return False
    return True


def spine_decompose(w: list[int]) -> list[list[int]]:
    hs = prefix_heights(w)
    fm = coefficient_future_minima(hs)  # includes 0 and len(w)
    assert fm[0] == 0 and fm[-1] == len(w)
    return [w[fm[i]:fm[i + 1]] for i in range(len(fm) - 1)]


# ----------------------------------------------------------------------------
# descent tree
# ----------------------------------------------------------------------------

def stopping_and_landing(n: int) -> tuple[int, int, int, int]:
    """(sigma, D(n), kappa, o) for the T-map: sigma = first k with T^k n < n,
    D = T^sigma n, kappa = first k with 3^{o_k} < 2^k (coefficient stopping time),
    o = number of odd steps among the first sigma."""
    if n == 1:
        return (0, 1, 0, 0)  # convention: 1 is the root
    x, k, o = n, 0, 0
    kappa = None
    while True:
        if x % 2 == 1:
            x = (3 * x + 1) // 2
            o += 1
        else:
            x //= 2
        k += 1
        if kappa is None and 3 ** o < 2 ** k:
            kappa = k
        if x < n:
            return (k, x, kappa, o)


def first_descent_words(K: int):
    """Yield (word, o, C) for words w of length k <= K with first coefficient descent
    exactly at k: prefix of length k-1 no-descent, 3^{o} < 2^k.  DFS."""
    stack = [([], 0, 0)]  # (word, o, C) with C the carry of the word so far
    while stack:
        w, o, C = stack.pop()
        k = len(w)
        if k >= K:
            continue
        for bit in (0, 1):
            w2 = w + [bit]
            o2 = o + bit
            # carry: C_w = sum_{j: w_j = 1} 2^j 3^{#ones after j}; appending a 1 at
            # position k multiplies earlier contributions by 3 and adds 2^k.
            C2 = 3 * C + 2 ** k if bit else C
            if 3 ** o2 > 2 ** (k + 1):
                stack.append((w2, o2, C2))
            elif bit == 0:
                yield (w2, o2, C2)
            # a 1 never causes the first descent: 3^{o+1} > 3^{o} > 2^{k} and 3 > 2


def source_residue(o: int, C: int, k: int) -> int:
    return (-C * pow(3, -o, 2 ** k)) % (2 ** k)


# ----------------------------------------------------------------------------
# THM-4503 cells as Young-lattice intervals
# ----------------------------------------------------------------------------

def f_threshold(j: int) -> int:
    o = 0
    while not (3 ** o > 2 ** (o + j)):
        o += 1
    return o


def cell_words_prefix(a: int, b: int) -> list[tuple[int, ...]]:
    """Words with a ones and b zeros, every nonempty prefix with i ones, j zeros
    having 3^i > 2^{i+j}."""
    out = []
    for w in enumerate_words(a + b):
        if sum(w) != a:
            continue
        i = j = 0
        ok = True
        for bit in w:
            if bit:
                i += 1
            else:
                j += 1
            if not 3 ** i > 2 ** (i + j):
                ok = False
                break
        if ok:
            out.append(tuple(w))
    return out


def word_to_partition(w: tuple[int, ...]) -> tuple[int, ...]:
    """e_i = number of zeros before the i-th one."""
    e, zeros = [], 0
    for bit in w:
        if bit:
            e.append(zeros)
        else:
            zeros += 1
    return tuple(e)


def carry(w) -> int:
    ones_after = [0] * (len(w) + 1)
    for j in range(len(w) - 1, -1, -1):
        ones_after[j] = ones_after[j + 1] + (1 if w[j] else 0)
    return sum(2 ** j * 3 ** ones_after[j + 1] for j in range(len(w)) if w[j])


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def part1_orbits():
    print("== P1: value-time poset, records, excursion forest on famous orbits ==")
    for n in [27, 703, 6171, 77031, 837799, 8400511, 63728127, 670617279]:
        orb = t_orbit_to_one(n)
        hs = word_heights(orb)
        info = check_value_time_poset(orb)
        lows = strict_lower_records(orb)
        highs = strict_upper_records(orb)
        vfm = value_future_minima(orb)
        cfm = coefficient_future_minima(hs)
        par = excursion_parent(hs)
        roots = [j for j in range(len(hs)) if par[j] == -1]
        # roots of the excursion forest (coefficient heights) = strict lower records
        # of the height walk; compare with strict lower records of the values.
        hlows = []
        cur = None
        for j in range(len(hs)):
            if cur is None or height_lt(hs[j][0], hs[j][1], cur[0], cur[1]):
                hlows.append(j)
                cur = hs[j]
        assert roots == hlows
        # value future minima within the window: only the final 1 (proposition)
        assert vfm == [len(orb) - 1]
        # coefficient future minima within the window: subset of value future minima
        # except possibly at the end where carry matters; check inclusion for j < len-1
        assert all(j in vfm for j in cfm if j < len(orb) - 1)
        # landing constraint along the chain of lower records
        for a_, b_ in zip(lows, lows[1:]):
            assert orb[a_] // 2 <= orb[b_] < orb[a_] or (orb[a_] % 2 == 1 and (orb[a_] + 1) // 2 <= orb[b_] < orb[a_])
        print(f"n={n:>10}: steps={len(orb)-1:>5} peak={max(orb):>14} leaders={len(highs):>3} "
              f"lower-records={len(lows):>3} floor(log2 n)+1={int(math.log2(n))+1:>3} "
              f"coef-future-minima(window)={len(cfm):>2} forest-roots={len(roots):>3} "
              f"max-depth={max_depth(par):>3}")


def max_depth(par: list[int]) -> int:
    depth = [0] * len(par)
    for j in range(len(par)):
        depth[j] = 0 if par[j] < 0 else depth[par[j]] + 1
    return max(depth) if depth else 0


def part2_blocks(KENUM: int = 16, KSER: int = 60):
    print("\n== P2: spine blocks = positive min-ending words; W = 1/(1-b); b = 1 - exp(-sum B_n t^n/n) ==")
    W = no_descent_counts(KSER)
    B = binomial_tails(KSER)
    inv = series_inverse(W)
    b_from_W = [Fraction(0)] + [-inv[k] for k in range(1, KSER + 1)]  # 1 - 1/W
    ser = [Fraction(0)] + [Fraction(-B[n], n) for n in range(1, KSER + 1)]
    e = series_exp(ser)  # exp(-sum B_n t^n / n)
    b_from_B = [Fraction(0)] + [-e[k] for k in range(1, KSER + 1)]
    assert b_from_W == b_from_B, "the two block generating functions differ"
    assert all(x.denominator == 1 and x >= 0 for x in b_from_W)
    # direct enumeration for small lengths
    for k in range(1, KENUM + 1):
        nd = [w for w in enumerate_words(k) if is_no_descent(w)]
        assert len(nd) == W[k]
        blocks = [w for w in nd if is_spine_block(w)]
        assert len(blocks) == b_from_W[k], (k, len(blocks), b_from_W[k])
        # unique decomposition: every no-descent word is a concatenation of blocks
        # and the concatenation of any blocks is no-descent with that decomposition
        for w in nd:
            dec = spine_decompose(w)
            assert all(is_spine_block(x) for x in dec)
            assert sum(dec, []) == w
    # concatenation test: random pairs of blocks (all pairs for k <= 8)
    small_blocks = [w for k in range(1, 9) for w in enumerate_words(k) if is_spine_block(w)]
    for x in small_blocks:
        for y in small_blocks:
            w = x + y
            assert is_no_descent(w)
            assert spine_decompose(w) == [x, y]
    print("   W_k, b_k (spine blocks), B_k for k <= 24:")
    for k in range(1, 25):
        print(f"   k={k:>2} W={W[k]:>8} b={int(b_from_W[k]):>8} B={B[k]:>8} b/W={float(b_from_W[k]/W[k]):.4f}")
    # block lengths: l >= 2 is a block length iff {l log_3 2} > log_3 2, i.e. the
    # largest a with 3^a < 2^l has 3^a < 2^(l-1); then every block has a+1 ones.
    for l in range(2, KSER + 1):
        a_max = 0
        while 3 ** (a_max + 1) < 2 ** l:
            a_max += 1
        exists = 3 ** a_max < 2 ** (l - 1)
        assert (b_from_W[l] > 0) == exists, (l, b_from_W[l], exists)
    for l in range(2, KENUM + 1):
        a_max = 0
        while 3 ** (a_max + 1) < 2 ** l:
            a_max += 1
        for w in enumerate_words(l):
            if is_spine_block(w):
                assert sum(w) == a_max + 1 and w[0] == 1
                assert height_lt(sum(w), l, 1, 1)  # H(w) < log2 3 - 1
    lengths = [l for l in range(1, KSER + 1) if b_from_W[l] > 0]
    beta = math.log2(3) / (math.log2(3) - 1)
    assert lengths[1:] == [1 + int(b * beta) for b in range(1, len(lengths))]
    print(f"   block lengths l <= {KSER}: {lengths}")
    print(f"   = {{1}} U {{1 + floor(b beta)}}, beta = log2(3)/(log2(3)-1) = {beta:.6f}; density 1 - log_3 2 = {1 - math.log2(3)**-1:.4f}")
    ratio = [float(b_from_W[k] / W[k]) for k in range(40, KSER + 1) if b_from_W[k] > 0]
    print(f"   b_k/W_k over block lengths in 40..{KSER}: min {min(ratio):.4f} max {max(ratio):.4f}")
    print(f"   b_k 2^(-h k) k^(3/2) at k = 40, 50, 60: " +
          ", ".join(f"{float(b_from_W[k]) * 2 ** (-0.9499555 * k) * k ** 1.5:.3f}" for k in (40, 50, 60)))
    return W, [int(x) for x in b_from_W]


def part3_descent_tree(NMAX: int = 2_000_000, KWORDS: int = 26, MCHECK: int = 100_000):
    print("\n== P3: the descent tree D(m) = first T-iterate below m ==")
    sigma = [0] * (NMAX + 1)
    D = [0] * (NMAX + 1)
    kappa = [0] * (NMAX + 1)
    odd_steps = [0] * (NMAX + 1)
    terras_fail = []
    for m in range(2, NMAX + 1):
        s, d, k, o = stopping_and_landing(m)
        sigma[m], D[m], kappa[m], odd_steps[m] = s, d, k, o
        assert m // 2 <= d < m, (m, d)  # landing in [m/2, m)
        if k != s:
            terras_fail.append(m)
    print(f"   all 2 <= m <= {NMAX}: D(m) in [m/2, m)  OK; max sigma = {max(sigma)} at m = {sigma.index(max(sigma))}")
    print(f"   Terras equality kappa = sigma failures in [2, {NMAX}]: {terras_fail[:10]} (count {len(terras_fail)})")
    # lower-record count = length of the D-chain; bound floor(log2 n) + 1, equality iff power of 2
    chain_len = [0] * (NMAX + 1)
    chain_len[1] = 1
    eq = []
    viol = 0
    for n in range(2, NMAX + 1):
        chain_len[n] = 1 + chain_len[D[n]]
        bound = int(math.log2(n)) + 1
        if chain_len[n] < bound:
            viol += 1
        if chain_len[n] == bound:
            eq.append(n)
    powers = [1 << j for j in range(1, int(math.log2(NMAX)) + 1)]
    print(f"   #lower records >= floor(log2 n) + 1 violated {viol} times; equality set = powers of two: {eq == powers} ({len(eq)} values)")
    # in-degrees for landings m <= NMAX/2 (complete: sources are <= 2m <= NMAX)
    indeg = Counter(D[m] for m in range(2, NMAX + 1))
    half = NMAX // 2
    dist = Counter(indeg.get(m, 0) for m in range(1, half + 1))
    print(f"   in-degree distribution over landings 1..{half}: " +
          ", ".join(f"{d}:{c}" for d, c in sorted(dist.items())))
    mean = sum(indeg.get(m, 0) for m in range(1, half + 1)) / half
    mx = max(indeg.get(m, 0) for m in range(1, half + 1))
    argmx = [m for m in range(1, half + 1) if indeg.get(m, 0) == mx][:5]
    print(f"   mean in-degree {mean:.4f}, max {mx} at {argmx}")
    # two-place structure: first-descent words up to length KWORDS
    words = list(first_descent_words(KWORDS))
    Wc = no_descent_counts(KWORDS)
    Fk = Counter(len(w) for w, _, _ in words)
    for k in range(1, KWORDS + 1):
        assert Fk[k] == 2 * Wc[k - 1] - Wc[k], (k, Fk[k], 2 * Wc[k - 1] - Wc[k])
    print(f"   first-descent words with length <= {KWORDS}: {len(words)} (F_k = 2 W_(k-1) - W_k checked)")
    weight = sum(Fraction(1, 3 ** o) for _, o, _ in words)
    weight2 = sum(Fraction(1, 2 ** len(w)) for w, _, _ in words)
    print(f"   sum_w 3^(-o(w)) = {float(weight):.6f}   sum_w 2^(-|w|) = {float(weight2):.6f} (both over |w| <= {KWORDS})")
    # predicted in-degree from the 3-adic landing classes
    pred = Counter()
    for w, o, C in words:
        k = len(w)
        r = source_residue(o, C, k)
        mod2 = 2 ** k
        assert (3 ** o * r + C) % mod2 == 0
        m0 = (3 ** o * r + C) // mod2
        N_thr = Fraction(C, mod2 - 3 ** o)  # source must exceed this
        mod3 = 3 ** o
        # landings m = m0 + t*mod3 with source (2^k m - C)/3^o = r + t*2^k > N_thr
        m = m0
        t = 0
        while m <= MCHECK:
            src = r + t * mod2
            if src > N_thr:
                pred[m] += 1
            m += mod3
            t += 1
    direct_short = Counter(D[m] for m in range(2, NMAX + 1) if sigma[m] <= KWORDS)
    long_src = [m for m in range(2, NMAX + 1) if sigma[m] > KWORDS and D[m] <= MCHECK]
    bad = [m for m in range(1, MCHECK + 1) if pred.get(m, 0) != direct_short.get(m, 0)]
    print(f"   landings m <= {MCHECK}: predicted in-degree from 3-adic classes (sigma <= {KWORDS}) "
          f"matches direct count for all m: {not bad} (mismatches {bad[:5]}); "
          f"sources with sigma > {KWORDS} landing <= {MCHECK}: {len(long_src)}")
    # example: the class of a first-descent word and its landing progression
    for w, o, C in words:
        if len(w) == 5 and o == 3:
            k = len(w)
            r = source_residue(o, C, k)
            m0 = (3 ** o * r + C) // 2 ** k
            N_thr = Fraction(C, 2 ** k - 3 ** o)
            ts = [t for t in range(6) if r + t * 2 ** k > N_thr]
            lands = [D[r + t * 2 ** k] for t in ts]
            print(f"   example word {''.join(map(str, w))}: o={o} C={C} r={r} (mod {2**k}) N={float(N_thr):.2f} -> landings "
                  f"{lands} = {m0} + t*{3**o} for t in {ts}: {lands == [m0 + t * 3 ** o for t in ts]}")
            break
    return sigma, D


def part4_cells(AB: int = 12):
    print("\n== P4: THM-4503 cells as intervals of Young's lattice; carry strictly monotone ==")
    total_checked = 0
    for a in range(1, AB):
        for b in range(1, AB - a + 1):
            ws = cell_words_prefix(a, b)
            # boundary partition g: e_i <= g_i where g_i = min{j-1 : f(j) >= i} (j <= b), else b
            g = []
            for i in range(1, a + 1):
                gi = b
                for j in range(1, b + 1):
                    if f_threshold(j) >= i:
                        gi = j - 1
                        break
                g.append(gi)
            parts = set()
            # all partitions e_1 <= ... <= e_a with e_i <= g_i; the j-th zero needs
            # f(j) ones before it, impossible when f(b) > a (THM-4503: a < m_b, no words)
            def rec(i, lo, cur):
                if i == a:
                    parts.add(tuple(cur))
                    return
                for v in range(lo, g[i] + 1):
                    rec(i + 1, v, cur + [v])
            if f_threshold(b) <= a:
                rec(0, 0, [])
            wp = set(word_to_partition(w) for w in ws)
            assert wp == parts, (a, b)
            # carry monotone: adding one cell (e_i -> e_i + 1, staying a partition) increases C
            by_part = {word_to_partition(w): carry(w) for w in ws}
            for e, c in by_part.items():
                for i in range(a):
                    e2 = list(e)
                    e2[i] += 1
                    e2 = tuple(e2)
                    if e2 in by_part:
                        assert by_part[e2] > c
                        total_checked += 1
            if (a, b) == (5, 1):
                res = [(''.join(map(str, w)), (-carry(w) * pow(3, -a, 64)) % 64) for w in ws]
                print(f"   cell (5,1): {res}  (THM-4503: 27, 39, 47, 31)")
    print(f"   all cells a+b <= {AB}: words = partitions inside g(a,b); {total_checked} covering pairs with strictly larger carry")


def part5_sturmian(J: int = 20000):
    print("\n== P5: the Sturmian element of E_inf (upper mechanical word of slope log_3 2) ==")
    # o_j = ceil(j * alpha), alpha = log_3 2; exact via 3^o > 2^j >= 3^(o-1)... o_j = min{o : 3^o > 2^j}? that is ceil(j log_3 2) when not integer
    hs = [(0, 0)]
    o = 0
    word = []
    for j in range(1, J + 1):
        while not (3 ** o > 2 ** j):
            o += 1
            # o increases by at most 1 per step
        word.append(1 if o > hs[-1][0] else 0)
        hs.append((o, j))
    assert all(3 ** o_ > 2 ** j_ for o_, j_ in hs[1:])
    hf = [height_float(o_, j_) for o_, j_ in hs]
    print(f"   heights in [0, log2 3): min {min(hf[1:]):.6f} max {max(hf):.6f} (log2 3 = {math.log2(3):.6f})")
    # coefficient future minima within [0, J]: their gaps
    fm = coefficient_future_minima(hs)
    gaps = [b - a for a, b in zip(fm, fm[1:])]
    print(f"   spine (future minima) times up to J={J}: {fm[:12]} ... ; block lengths {gaps[:12]} ...")
    print(f"   distinct block lengths: {sorted(set(gaps))}")
    # least positive residue of the cylinder of the first k letters (T-map coding: n mod 2^k)
    # r_k = -C_w 3^{-o} mod 2^k
    rs = []
    for k in (8, 16, 24, 32, 40, 48, 56, 64):
        w = word[:k]
        C = carry(w)
        o_ = sum(w)
        r = (-C * pow(3, -o_, 2 ** k)) % 2 ** k
        rs.append((k, r, r / 2 ** k))
    print("   least residues r_k of the Sturmian cylinders: " + ", ".join(f"k={k}: {r} ({fr:.3f} 2^k)" for k, r, fr in rs))
    # real series diverges: partial sums of 2^{-h_k}
    S = sum(2.0 ** (-h) for h in hf[:J + 1])
    print(f"   sum_(k<={J}) 2^(-h_k) = {S:.1f} (diverges linearly; every term >= 1/3)")


def part6_drift_control(steps: int = 3000):
    print("\n== P6: DRIFT control: 5x+1 orbit of 7, spine and blocks within the window ==")
    q = 5
    m = 7
    vals = [m]
    hs = [(0, 0)]
    o = 0
    for j in range(1, steps + 1):
        if m % 2 == 1:
            m = (q * m + 1) // 2
            o += 1
        else:
            m //= 2
        vals.append(m)
        hs.append((o, j))

    def lt5(o1, j1, o2, j2):
        do, dj = o1 - o2, j1 - j2
        if do == 0 and dj == 0:
            return False
        lhs = 5 ** max(do, 0) * 2 ** max(-dj, 0)
        rhs = 2 ** max(dj, 0) * 5 ** max(-do, 0)
        return lhs < rhs

    cfm, cur = [], None
    for j in range(len(hs) - 1, -1, -1):
        if cur is None or lt5(hs[j][0], hs[j][1], cur[0], cur[1]):
            cfm.append(j)
            cur = hs[j]
    cfm = cfm[::-1]
    vfm = value_future_minima(vals)
    # coefficient future minima are value future minima (positive carry); check within window
    inner = [j for j in cfm if j < steps - 200]
    assert all(j in vfm for j in inner)
    # step out of a value future minimum is a rising step
    for j in vfm:
        if j + 1 < len(vals):
            assert vals[j + 1] > vals[j]
    gaps = [b - a for a, b in zip(cfm, cfm[1:])]
    print(f"   window {steps} steps: coefficient future minima {len(cfm)}, value future minima {len(vfm)}, inclusion OK;")
    print(f"   block lengths: max {max(gaps)}, mean {sum(gaps)/len(gaps):.2f}, first 20: {gaps[:20]}")
    print(f"   log2 of the value at the end: {math.log2(vals[-1]):.1f}")


if __name__ == "__main__":
    part1_orbits()
    part2_blocks()
    part3_descent_tree()
    part4_cells()
    part5_sturmian()
    part6_drift_control()
    print("\nALL CHECKS PASSED")
