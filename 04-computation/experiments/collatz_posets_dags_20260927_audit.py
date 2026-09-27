#!/usr/bin/env python3
"""Independent audit of collatz_posets_dags_20260927 (candidate THM-4514).

Written WITHOUT importing or copying 04-computation/experiments/collatz_posets_dags_20260927.py;
every routine below was re-derived from the definitions in section 2 of the note
(T-map, parity word, heights h_j = o_j log_2 3 - j, carry, E_inf, sigma, kappa, D).
All assertions use exact integer arithmetic; floats are display only.

Run:  python3 04-computation/experiments/collatz_posets_dags_20260927_audit.py
Output: 05-knowledge/results/collatz_posets_dags_20260927_audit.out
"""
from __future__ import annotations

import math
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations

RESULTS: list[tuple[str, bool]] = []


def check(name: str, ok, detail: str = "") -> bool:
    ok = bool(ok)
    RESULTS.append((name, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    return ok


# --------------------------------------------------------------------------------------
# exact heights: h(o, j) = o log_2 q - j, compared through q^o 2^j
# --------------------------------------------------------------------------------------
_POW: dict[int, list[int]] = {3: [1], 5: [1]}


def qpow(q: int, n: int) -> int:
    L = _POW[q]
    while len(L) <= n:
        L.append(L[-1] * q)
    return L[n]


def h_lt(o1: int, j1: int, o2: int, j2: int, q: int = 3) -> bool:
    """h(o1, j1) < h(o2, j2)  <=>  q^o1 2^j2 < q^o2 2^j1 (exact; no ties for (o,j) distinct)."""
    return (qpow(q, o1) << j2) < (qpow(q, o2) << j1)


def h_float(o: int, j: int, q: int = 3) -> float:
    return o * math.log2(q) - j


def a_of(l: int, q: int = 3) -> int:
    """least a with q^a > 2^l (= ceil(l log_q 2) for l >= 1)."""
    a = 0
    while qpow(q, a) <= (1 << l):
        a += 1
    return a


# --------------------------------------------------------------------------------------
# maps, words, carries
# --------------------------------------------------------------------------------------
def Tmap(x: int, q: int = 3) -> int:
    return x >> 1 if (x & 1) == 0 else (q * x + 1) >> 1


def parity_word(n: int, k: int, q: int = 3) -> list[int]:
    w, x = [], n
    for _ in range(k):
        w.append(x & 1)
        x = Tmap(x, q)
    return w


def carry(w, q: int = 3) -> int:
    """C_w = sum over ones at 0-indexed position p of 2^p q^(#ones after p)."""
    C, after = 0, 0
    for p in range(len(w) - 1, -1, -1):
        if w[p]:
            C += (1 << p) * qpow(q, after)
            after += 1
    return C


def prefix_oj(w):
    out, o = [(0, 0)], 0
    for j, b in enumerate(w, 1):
        o += b
        out.append((o, j))
    return out


def is_no_descent(w, q: int = 3) -> bool:
    o = 0
    for j, b in enumerate(w, 1):
        o += b
        if qpow(q, o) <= (1 << j):
            return False
    return True


def is_block(w, q: int = 3) -> bool:
    """spine block: H = h_l > 0 and h_j > H for 1 <= j < l."""
    l = len(w)
    pre = prefix_oj(w)
    o = pre[l][0]
    if qpow(q, o) <= (1 << l):
        return False
    for j in range(1, l):
        if not h_lt(o, l, pre[j][0], j, q):
            return False
    return True


def future_minima(pre, q: int = 3) -> list[int]:
    """indices j in [0, K] with h_j < h_k for all k in (j, K] (window version), exact."""
    out, cur = [], None
    for j in range(len(pre) - 1, -1, -1):
        if cur is None or h_lt(pre[j][0], pre[j][1], cur[0], cur[1], q):
            out.append(j)
            cur = pre[j]
    return out[::-1]


def value_future_minima(vals) -> list[int]:
    out, cur = [], None
    for j in range(len(vals) - 1, -1, -1):
        if cur is None or vals[j] < cur:
            out.append(j)
            cur = vals[j]
    return out[::-1]


def strict_lower_records(vals) -> list[int]:
    out, cur = [], None
    for i, v in enumerate(vals):
        if cur is None or v < cur:
            out.append(i)
            cur = v
    return out


def strict_upper_records(vals) -> list[int]:
    out, cur = [], None
    for i, v in enumerate(vals):
        if cur is None or v > cur:
            out.append(i)
            cur = v
    return out


def excursion_parents(pre, q: int = 3) -> list[int]:
    """parent(j) = nearest i < j with h_i < h_j (proved equal to the largest ancestor); -1 for roots."""
    par, stack = [-1] * len(pre), []
    for j in range(len(pre)):
        while stack and not h_lt(pre[stack[-1]][0], pre[stack[-1]][1], pre[j][0], pre[j][1], q):
            stack.pop()
        par[j] = stack[-1] if stack else -1
        stack.append(j)
    return par


# --------------------------------------------------------------------------------------
# A. sanity of the two-place identity
# --------------------------------------------------------------------------------------
def part_A():
    print("== A. two-place identity T^k n = (3^o n + C_w)/2^k (sanity of the carry convention) ==")
    ok = True
    for n in [1, 3, 7, 27, 97, 703, 12345, 999999, 2 ** 40 + 1]:
        for k in (1, 5, 17, 40):
            w = parity_word(n, k)
            o = sum(w)
            x = n
            for _ in range(k):
                x = Tmap(x)
            ok &= (qpow(3, o) * n + carry(w)) == (x << k)
    check("T^k n = (3^o n + C_w)/2^k for sample n, k", ok)


# --------------------------------------------------------------------------------------
# B. W_k (two routes), B_n, the block generating function, block enumeration, tightness
# --------------------------------------------------------------------------------------
def W_by_dp(K: int, q: int = 3) -> list[int]:
    W, layer = [1], {0: 1}
    for j in range(1, K + 1):
        nxt: dict[int, int] = defaultdict(int)
        for o, c in layer.items():
            for b in (0, 1):
                if qpow(q, o + b) > (1 << j):
                    nxt[o + b] += c
        layer = nxt
        W.append(sum(layer.values()))
    return W


def B_tails(K: int, q: int = 3) -> list[int]:
    return [0] + [sum(math.comb(n, j) for j in range(n + 1) if qpow(q, j) > (1 << n)) for n in range(1, K + 1)]


def W_by_spitzer(B: list[int], K: int) -> list[int]:
    W = [1]
    for k in range(1, K + 1):
        s = sum(B[n] * W[k - n] for n in range(1, k + 1))
        if s % k:
            return []
        W.append(s // k)
    return W


def blocks_from_W(W: list[int]) -> list[Fraction]:
    K = len(W) - 1
    inv = [Fraction(1)]
    for k in range(1, K + 1):
        inv.append(-sum(Fraction(W[i]) * inv[k - i] for i in range(1, k + 1)))
    return [Fraction(0)] + [-inv[k] for k in range(1, K + 1)]


def blocks_from_B(B: list[int]) -> list[Fraction]:
    K = len(B) - 1
    e = [Fraction(1)]
    for k in range(1, K + 1):
        e.append(-sum(Fraction(B[n]) * e[k - n] for n in range(1, k + 1)) / k)
    return [Fraction(0)] + [-e[k] for k in range(1, K + 1)]


def all_words(l: int):
    for x in range(1 << l):
        yield [(x >> (l - 1 - i)) & 1 for i in range(l)]


def no_descent_words_upto(L: int, q: int = 3) -> dict[int, list[list[int]]]:
    out: dict[int, list[list[int]]] = defaultdict(list)
    stack = [([], 0)]
    while stack:
        w, o = stack.pop()
        k = len(w)
        if k == L:
            continue
        for b in (0, 1):
            if qpow(q, o + b) > (1 << (k + 1)):
                w2 = w + [b]
                out[k + 1].append(w2)
                stack.append((w2, o + b))
    return out


def floor_b_beta(b: int) -> int:
    """floor(b beta), beta = log_2 3/(log_2 3 - 1), exactly: largest n with 3^(n-b) < 2^n."""
    n = b
    while qpow(3, n + 1 - b) < (1 << (n + 1)):
        n += 1
    return n


def part_B(KSER: int = 60, LENUM: int = 20, LBRUTE: int = 16):
    print("\n== B. no-descent counts, block generating function, block enumeration, tightness ==")
    W = W_by_dp(KSER)
    B = B_tails(KSER)
    W2 = W_by_spitzer(B, KSER)
    check("W_k by ballot DP == W_k by Spitzer recurrence k W_k = sum B_n W_(k-n), k <= 60", W == W2)
    check("W_1..W_24 as printed in the note's table", W[1:25] == [1, 1, 2, 3, 4, 8, 13, 19, 38, 64, 128, 226, 367, 734, 1295, 2114, 4228, 7495, 14990, 27328, 46611, 93222, 168807, 286581])
    bW = blocks_from_W(W)
    bB = blocks_from_B(B)
    check("1 - 1/W(t) == 1 - exp(-sum B_n t^n/n) coefficientwise to t^60", bW == bB)
    check("all b_l are non-negative integers (l <= 60)", all(x.denominator == 1 and x >= 0 for x in bW))
    b = [int(x) for x in bW]
    note_b = [1, 0, 1, 0, 0, 2, 0, 0, 7, 0, 30, 0, 0, 113, 0, 0, 525, 0, 2652, 0, 0, 11433, 0, 0]
    check("b_1..b_24 == the note's list", b[1:25] == note_b, str(b[1:25]))

    # direct enumeration: brute force over ALL 2^l words for l <= LBRUTE, DFS over no-descent words to LENUM
    brute_nd, brute_blk = {}, {}
    for l in range(1, LBRUTE + 1):
        nd = blk = 0
        for w in all_words(l):
            if is_no_descent(w):
                nd += 1
                if is_block(w):
                    blk += 1
        brute_nd[l], brute_blk[l] = nd, blk
    check(f"brute force over all 2^l words: #no-descent == W_l for l <= {LBRUTE}", all(brute_nd[l] == W[l] for l in range(1, LBRUTE + 1)))
    check(f"brute force over all 2^l words: #blocks == b_l (from 1 - 1/W) for l <= {LBRUTE}", all(brute_blk[l] == b[l] for l in range(1, LBRUTE + 1)))
    ndw = no_descent_words_upto(LENUM)
    check(f"DFS: #no-descent words == W_l for l <= {LENUM}", all(len(ndw[l]) == W[l] for l in range(1, LENUM + 1)))
    blocks = {l: [w for w in ndw[l] if is_block(w)] for l in range(1, LENUM + 1)}
    counts = [len(blocks[l]) for l in range(1, LENUM + 1)]
    check(f"DFS: #blocks == b_l (from 1 - 1/W) for l <= {LENUM}", counts == b[1:LENUM + 1], f"b_1..b_{LENUM} = {counts}")
    # every block is a no-descent word (a block has all prefix heights > H > 0): consistency of using the DFS
    check("every block of length <= 16 found by brute force is no-descent (so the DFS misses none)",
          all(brute_blk[l] == sum(1 for w in all_words(l) if is_block(w)) for l in range(1, LBRUTE + 1)))

    # tightness (Theorem 3 (ii)): first letter 1, a(l) ones, H < log2 3 - 1 for l >= 2; carry bound (iii)
    first_ok = ones_ok = height_ok = carry_ok = True
    for l in range(2, LENUM + 1):
        a = a_of(l)
        for w in blocks[l]:
            first_ok &= w[0] == 1
            ones_ok &= sum(w) == a
            height_ok &= qpow(3, a - 1) < (1 << (l - 1))           # a log2 3 - l < log2 3 - 1
            carry_ok &= 2 * carry(w) < a * qpow(3, a)                # C_w / 3^a < a/2
    check(f"every block of length 2..{LENUM} starts with 1", first_ok)
    check(f"every block of length 2..{LENUM} has exactly a(l) = ceil(l log_3 2) ones", ones_ok)
    check(f"every block of length 2..{LENUM} has 0 < H < log_2 3 - 1", height_ok)
    check(f"carry bound C_w < a 3^a / 2 on every block of length <= {LENUM}", carry_ok)
    # the one-letter block
    check("l = 1: the unique block is '1' with H = log_2 3 - 1 exactly", blocks[1] == [[1]])

    # block length law (Theorem 3 (ii)), exact, l = 2..60
    beatty = {qpow(3, a).bit_length() - 1 for a in range(1, 80)}       # floor(a log_2 3)
    rayleigh = {floor_b_beta(bb) for bb in range(1, 80)}                 # floor(b beta)
    law1 = law2 = law3 = law4 = True
    for l in range(2, KSER + 1):
        exists = b[l] > 0
        a = a_of(l)
        crit_frac = (1 << (l - 1)) > qpow(3, a - 1)      # {l log_3 2} > log_3 2  <=>  (a(l)-1) log_2 3 < l - 1
        law1 &= exists == crit_frac
        law2 &= exists == ((l - 1) not in beatty)
        law3 &= exists == ((l - 1) in rayleigh)
        if exists:
            law4 &= is_block([1] * a + [0] * (l - a))
    check("block length l (2..60) iff {l log_3 2} > log_3 2 (exact form 2^(l-1) > 3^(a(l)-1))", law1)
    check("block length l (2..60) iff l - 1 not a Beatty number floor(a log_2 3)", law2)
    check("block length l (2..60) iff l = 1 + floor(b beta) (Rayleigh complement), exact", law3)
    check("1^a 0^(l-a) is a block at every block length <= 60", law4)
    b100 = {x for x in beatty if x <= 100}
    r100 = {x for x in rayleigh if x <= 100}
    check("Beatty(log_2 3) and Beatty(beta) partition 1..100 (Rayleigh)", (b100 | r100) == set(range(1, 101)) and not (b100 & r100))
    lengths = [l for l in range(1, KSER + 1) if b[l] > 0]
    check("block lengths <= 60 == note's list", lengths == [1, 3, 6, 9, 11, 14, 17, 19, 22, 25, 28, 30, 33, 36, 38, 41, 44, 47, 49, 52, 55, 57, 60], str(lengths))
    dens = (len(lengths) - 1) / 59
    print(f"  block-length density on 2..60: {len(lengths) - 1}/59 = {dens:.4f}; 1 - log_3 2 = {1 - 1 / math.log2(3):.5f}; beta = {math.log2(3) / (math.log2(3) - 1):.6f}")
    ratios = [b[l] / W[l] for l in range(40, 61) if b[l] > 0]
    print(f"  b_l/W_l on block lengths 40..60: min {min(ratios):.6f} max {max(ratios):.6f} (note, 4 dp: [0.0705, 0.1124])")
    check("b_l/W_l on block lengths in 40..60: min rounds to 0.0705, max to 0.1124 (4 dp)", f"{min(ratios):.4f}" == "0.0705" and f"{max(ratios):.4f}" == "0.1124")

    # Theorem 3 (i): unique decomposition of no-descent words (<= 14) and free concatenation (blocks <= 8)
    dec_ok = True
    for l in range(1, 15):
        for w in ndw[l]:
            fm = future_minima(prefix_oj(w))
            dec_ok &= fm[0] == 0 and fm[-1] == l
            parts = [w[fm[i]:fm[i + 1]] for i in range(len(fm) - 1)]
            dec_ok &= all(is_block(p) for p in parts)
    check("every no-descent word of length <= 14 cuts at its window future minima into blocks", dec_ok)
    small = [w for l in range(1, 9) for w in blocks[l]]
    cat_ok = True
    for x in small:
        for y in small:
            w = x + y
            fm = future_minima(prefix_oj(w))
            cat_ok &= is_no_descent(w) and fm == [0, len(x), len(w)]
    check("concatenation of any two blocks of length <= 8 is no-descent with exactly that decomposition", cat_ok)
    # dominance: no-descent iff o_j >= ceil(j log_3 2) for all j (principal filter above the mechanical word), l <= 16
    dom_ok = True
    for l in range(1, 17):
        for w in all_words(l):
            pre = prefix_oj(w)
            dom_ok &= is_no_descent(w) == all(pre[j][0] >= a_of(j) for j in range(1, l + 1))
    check("no-descent words of length <= 16 == words dominating the upper mechanical word o_j = ceil(j log_3 2)", dom_ok)
    return W, b


# --------------------------------------------------------------------------------------
# C. descent tree scan
# --------------------------------------------------------------------------------------
def scan_descent(NMAX: int):
    P3 = [3 ** i for i in range(700)]
    sigma = [0] * (NMAX + 1)
    D = [0] * (NMAX + 1)
    kappa = [0] * (NMAX + 1)
    defects = []
    for m in range(2, NMAX + 1):
        x, k, o, kap = m, 0, 0, 0
        while x >= m:
            if x & 1:
                x = (3 * x + 1) >> 1
                o += 1
            else:
                x >>= 1
            k += 1
            if not kap and P3[o] < (1 << k):
                kap = k
        sigma[m], D[m], kappa[m] = k, x, kap
        if kap != k:
            defects.append(m)
    return sigma, D, kappa, defects


def part_C(NMAX: int = 2_000_000):
    print(f"\n== C. descent tree D(m) = T^sigma(m) m for 2 <= m <= {NMAX} ==")
    t0 = time.time()
    sigma, D, kappa, defects = scan_descent(NMAX)
    print(f"  scan time {time.time() - t0:.1f} s")
    rng_ok = all((m + 1) // 2 <= D[m] <= m - 1 for m in range(2, NMAX + 1))
    check("D(m) in [ceil(m/2), m-1] for all 2 <= m <= NMAX", rng_ok)
    ms = max(sigma)
    check("max sigma = 224 at m = 1126015", ms == 224 and sigma.index(ms) == 1126015, f"max sigma {ms} at {sigma.index(ms)}")
    check("Terras equality kappa = sigma for all 2 <= m <= NMAX (no defect)", not defects, f"defects {defects[:5]}")
    check("kappa <= sigma everywhere", all(kappa[m] <= sigma[m] for m in range(2, NMAX + 1)))
    # record bound: #strict lower records = chain length to 1 >= floor(log2 m) + 1 = bit_length(m); equality iff power of two
    chain = [0] * (NMAX + 1)
    chain[1] = 1
    viol, eq = 0, []
    for n in range(2, NMAX + 1):
        chain[n] = 1 + chain[D[n]]
        bl = n.bit_length()
        if chain[n] < bl:
            viol += 1
        elif chain[n] == bl:
            eq.append(n)
    check("#lower records >= floor(log_2 m) + 1 for all m <= NMAX", viol == 0)
    powers = [1 << j for j in range(1, NMAX.bit_length()) if (1 << j) <= NMAX]
    check("equality set == powers of two (20 values)", eq == powers and len(eq) == 20, f"{len(eq)} values")
    # every dyadic shell below m contains a record (direct check, m <= 300000)
    shell_ok = True
    for m in range(2, 300_001):
        seen, x = set(), m
        while True:
            seen.add(x.bit_length() - 1)
            if x == 1:
                break
            x = D[x]
        shell_ok &= seen == set(range(m.bit_length()))
    check("every shell [2^j, 2^(j+1)), j <= floor(log2 m), contains a lower record (m <= 300000)", shell_ok)
    # in-degrees over landings 1..NMAX/2
    indeg = Counter(D[m] for m in range(2, NMAX + 1))
    half = NMAX // 2
    dist = Counter(indeg.get(m, 0) for m in range(1, half + 1))
    tot = sum(indeg.get(m, 0) for m in range(1, half + 1))
    mx = max(indeg.get(m, 0) for m in range(1, half + 1))
    argmx = [m for m in range(1, half + 1) if indeg.get(m, 0) == mx]
    print(f"  in-degree distribution over landings 1..{half}: " + ", ".join(f"{d}:{c}" for d, c in sorted(dist.items())))
    print(f"  total sources landing <= {half}: {tot}; mean in-degree {tot / half:.6f}; max {mx} at {argmx}")
    note_dist = {1: 531332, 2: 337608, 3: 78502, 4: 28890, 5: 15341, 6: 4769, 7: 2102, 8: 874, 9: 348, 10: 129, 11: 61, 12: 26, 13: 13, 15: 3, 16: 1, 17: 1}
    check("in-degree distribution over landings <= 10^6 == note's table", dict(dist) == note_dist)
    check("mean in-degree 1.6903 (exact 1690291/10^6), max 17 at 293501", tot == 1690291 and mx == 17 and argmx == [293501])
    check("no landing m <= 10^6 has in-degree 0 (2m -> m always) and none has 14", 0 not in dist and 14 not in dist)
    return sigma, D, kappa


# --------------------------------------------------------------------------------------
# D. first-descent words, the 2-adic -> 3-adic bijection, in-degree formula, c_D
# --------------------------------------------------------------------------------------
def first_descent_words(K: int):
    """(k, o, C, wordint) for all words of length k <= K whose FIRST coefficient descent is at k.
    Own derivation: extend only no-descent prefixes; a 1 keeps 3^o > 2^k (3 > 2); a 0 descends iff 3^o < 2^(k+1)."""
    out = []
    stack = [(0, 0, 0, 0)]       # (k, o, C, wordint) with the k-prefix no-descent (k = 0: empty)
    while stack:
        k, o, C, wi = stack.pop()
        if k == K:
            continue
        # append 1: position k, all earlier ones get one more one after them
        stack.append((k + 1, o + 1, 3 * C + (1 << k), (wi << 1) | 1))
        # append 0
        if qpow(3, o) > (1 << (k + 1)):
            stack.append((k + 1, o, C, wi << 1))
        else:
            out.append((k + 1, o, C, wi << 1))
    return out


def word_bits(wi: int, k: int) -> list[int]:
    return [(wi >> (k - 1 - i)) & 1 for i in range(k)]


def part_D(sigma, D, W, KW: int = 26, MCHECK: int = 100_000, MSMALL: int = 20_000):
    print(f"\n== D. first-descent words |w| <= {KW}: 2-adic source class -> 3-adic landing class; in-degree formula; c_D ==")
    t0 = time.time()
    words = first_descent_words(KW)
    print(f"  {len(words)} first-descent words of length <= {KW} ({time.time() - t0:.1f} s)")
    check("190069 first-descent words of length <= 26", len(words) == 190069)
    Fk = Counter(k for k, _, _, _ in words)
    check("F_k = 2 W_(k-1) - W_k for k <= 26", all(Fk[k] == 2 * W[k - 1] - W[k] for k in range(1, KW + 1)))
    # carry convention: incremental C equals carry() from the letters (sample) and the word is first-descent
    samp_ok = True
    for k, o, C, wi in words[::997]:
        w = word_bits(wi, k)
        samp_ok &= carry(w) == C and sum(w) == o and is_no_descent(w[:-1]) and qpow(3, o) < (1 << k) and w[-1] == 0
    check("sampled words: carry, ones, no-descent prefix, descent at k, last letter 0", samp_ok)
    # overshoot: 3^o > 2^(k-1) for every first-descent word EXCEPT w = '0' (k = 1, o = 0: equality)
    viol = [(k, o) for k, o, _, _ in words if qpow(3, o) <= (1 << (k - 1))]
    check("3^o > 2^(k-1) fails only for w = '0' (equality 1 = 1)", viol == [(1, 0)], str(viol[:5]))
    # partial sums, exact
    OM = 40
    S3 = S2 = 0
    per_k3 = defaultdict(int)
    per_k2 = defaultdict(int)
    for k, o, _, _ in words:
        per_k3[k] += qpow(3, OM - o)
        per_k2[k] += 1 << (KW - k)
    print("  partial sums (exact -> 9 decimals):  K   sum_{|w|<=K} 3^(-o)   sum_{|w|<=K} 2^(-|w|)   1 - W_K/2^K")
    acc3 = acc2 = 0
    sums2_ok = True
    for K in range(1, KW + 1):
        acc3 += per_k3[K]
        acc2 += per_k2[K]
        f3 = Fraction(acc3, qpow(3, OM))
        f2 = Fraction(acc2, 1 << KW)
        sums2_ok &= f2 == 1 - Fraction(W[K], 1 << K)
        if K <= 22 or K == KW:
            print(f"    {K:>2}   {float(f3):.9f}          {float(f2):.9f}          {float(1 - Fraction(W[K], 1 << K)):.9f}")
    check("sum_{|w|<=K} 2^(-|w|) == 1 - W_K/2^K exactly for every K <= 26", sums2_ok)
    c26 = Fraction(acc3, qpow(3, OM))
    p26 = Fraction(acc2, 1 << KW)
    check("sum_{|w|<=26} 3^(-o(w)) = 1.669582 (6 dp) and sum 2^(-|w|) = 0.984542 (6 dp)",
          f"{float(c26):.6f}" == "1.669582" and f"{float(p26):.6f}" == "0.984542", f"{float(c26):.9f}, {float(p26):.9f}")
    tail = 2 * (1 - p26)                         # tail <= sum_{|w|>26} 2^(1-|w|)
    lo, hi = c26, c26 + tail
    print(f"  exact tail bound 2 W_26/2^26 = {float(tail):.9f}; from |w| <= 26 alone: c_D in [{float(lo):.6f}, {float(hi):.6f}]; W_26 = {W[26]}")
    check("upper end: 1.669582 + 0.030916 = 1.700498 <= 1.7005", float(hi) <= 1.7005)
    check("lower end as derived in the note: partial sum 1.669582 >= 1.6696 ? (NO: 1.6696 is the partial sum rounded UP)", float(lo) >= 1.6696)
    # rigorous bracket from a longer enumeration (accumulate only, no storage)
    K2 = 28
    S3b, S2b, nw = 0, 0, 0
    stack = [(0, 0, 0)]
    while stack:
        k, o, C = stack.pop()
        if k == K2:
            continue
        stack.append((k + 1, o + 1, 3 * C + (1 << k)))
        if qpow(3, o) > (1 << (k + 1)):
            stack.append((k + 1, o, C))
        else:
            S3b += qpow(3, OM - o)
            S2b += 1 << (K2 - (k + 1))
            nw += 1
    c28 = Fraction(S3b, qpow(3, OM))
    p28 = Fraction(S2b, 1 << K2)
    check(f"|w| <= {K2}: sum 2^(-|w|) == 1 - W_{K2}/2^{K2}", p28 == 1 - Fraction(W[K2], 1 << K2))
    hi28 = c28 + 2 * (1 - p28)
    print(f"  |w| <= {K2}: {nw} words; sum 3^(-o) = {float(c28):.9f}; sum 2^(-|w|) = {float(p28):.9f}; rigorous c_D in [{float(c28):.6f}, {float(hi28):.6f}]")
    check("rigorous bracket from |w| <= 28 confirms 1.6696 <= c_D <= 1.7005 (the note's interval is true; its lower end needs |w| = 27)", float(c28) >= 1.6696 and float(hi28) <= 1.7005)
    check("1 < partial sum <= c_D and upper bound < 2", 1 < c26 and hi < 2)
    # sharper tail from the actual overshoot is not needed; check 3^(-o) < 2^(1-k) for k >= 2 (used in the tail bound)
    check("2^(-k) < 3^(-o) < 2^(1-k) for every first-descent word with k >= 2", all((1 << k) > qpow(3, o) > (1 << (k - 1)) for k, o, _, _ in words if k >= 2))

    # affine bijection data: r, m0, N; thresholds
    data = []
    maxN = Fraction(0)
    below = []
    for k, o, C, wi in words:
        mod2 = 1 << k
        r = (-C * pow(3, -o, mod2)) % mod2 if k >= 1 else 0
        num = qpow(3, o) * r + C
        if num % mod2:
            check("m0 integrality", False, f"word {wi:b} k={k}")
        m0 = num // mod2
        N = Fraction(C, mod2 - qpow(3, o))
        if N > maxN:
            maxN = N
        if 1 <= r <= N:                      # r = 0 only for w = '0' (source 0 is not a positive integer)
            below.append((k, o, C, r, N))
        data.append((k, o, C, r, m0, mod2 - qpow(3, o)))
    print(f"  max N(w) over |w| <= 26: {float(maxN):.3f}; words whose least positive residue r <= N(w): {[(k, r, str(N)) for k, o, C, r, N in below]}")
    check("N(w) < 2^|w| for every first-descent word |w| <= 26 (at most the least residue is uncertified)", all(Fraction(C, gap) < (1 << k) for k, o, C, r, m0, gap in data))
    check("the only least residue r <= N(w) is r = 1 for w = 10 (the fixed point 1 = N(w)), so no Terras-defect source exists for |w| <= 26",
          below == [(2, 1, 1, 1, Fraction(1))])

    # in-degree prediction from the 3-adic landing classes, landings m <= MCHECK, sources with sigma <= KW
    pred = Counter()
    for k, o, C, r, m0, gap in data:
        mod3 = qpow(3, o)
        mod2 = 1 << k
        m, t = m0, 0
        while m <= MCHECK:
            if (r + t * mod2) * gap > C:          # source > N(w)
                pred[m] += 1
            m += mod3
            t += 1
    direct = Counter(D[m] for m in range(2, 2 * MCHECK + 1) if sigma[m] <= KW)
    bad = [m for m in range(1, MCHECK + 1) if pred.get(m, 0) != direct.get(m, 0)]
    check(f"3-adic class prediction == direct in-degree restricted to sigma <= {KW}, all landings m <= {MCHECK}", not bad, f"mismatches {bad[:5]}")
    long_src = sum(1 for m in range(2, 2 * MCHECK + 1) if sigma[m] > KW and D[m] <= MCHECK)
    check(f"sources with sigma > {KW} landing <= {MCHECK}: 2071", long_src == 2071, str(long_src))
    check("landing m0(w) + t 3^o is never counted for t < 0 (all sources > N >= 0 are >= r)", all(d[4] >= 0 for d in data))

    # FULL formula for landings m <= MSMALL, every source m' <= 2 MSMALL, using each source's OWN word (any sigma)
    full_ok = True
    by_land = Counter()
    for mp in range(2, 2 * MSMALL + 1):
        k = sigma[mp]
        w = parity_word(mp, k)
        o, C = sum(w), carry(w)
        mod2, mod3 = 1 << k, qpow(3, o)
        r = (-C * pow(3, -o, mod2)) % mod2
        m0 = (mod3 * r + C) // mod2
        t = (mp - r) // mod2
        cond = (is_no_descent(w[:-1]) and mod3 < mod2                     # w is a first-descent word (kappa = sigma)
                and (mod3 * r + C) % mod2 == 0
                and mp % mod2 == r and (mp - r) % mod2 == 0 and t >= 0
                and mp * (mod2 - mod3) > C                                # m' > N(w)
                and D[mp] == m0 + t * mod3                                # landing in the 3-adic class
                and (mod3 * mp + C) == D[mp] * mod2)
        full_ok &= cond
        if D[mp] <= MSMALL:
            by_land[D[mp]] += 1
    check(f"every source m' <= {2 * MSMALL}: its own first-descent word puts it in class r mod 2^k above N(w) and D(m') = m0 + t 3^o (full formula, any sigma)", full_ok)
    direct_all = Counter(D[mp] for mp in range(2, 2 * MSMALL + 1))
    check(f"in-degree of every landing m <= {MSMALL} == #sources m' <= 2m with D(m') = m (sources of m lie in (m, 2m])", all(by_land[m] == direct_all[m] for m in range(1, MSMALL + 1)))
    print(f"  max sigma among sources <= {2 * MSMALL}: {max(sigma[2:2 * MSMALL + 1])} (words far longer than {KW} enter the full formula)")
    # the example word 11100
    w = [1, 1, 1, 0, 0]
    o, C = sum(w), carry(w)
    r = (-C * pow(3, -o, 32)) % 32
    m0 = (27 * r + C) // 32
    lands = [D[r + 32 * t] for t in range(6)]
    check("example w = 11100: o=3, C=19, r=23, N=3.8, landings 20+27t, t=0..5", (o, C, r, m0) == (3, 19, 23, 20) and Fraction(C, 32 - 27) == Fraction(19, 5) and lands == [20 + 27 * t for t in range(6)], str(lands))


# --------------------------------------------------------------------------------------
# E. THM-4503 cells as intervals of Young's lattice; carry monotone; (5,1) residues
# --------------------------------------------------------------------------------------
def f_thr(j: int) -> int:
    o = 0
    while qpow(3, o) <= (1 << (o + j)):
        o += 1
    return o


def part_E(AB: int = 12):
    print("\n== E. cells (a, b) with a + b <= 12: words == partitions e <= g; carry strictly increasing under 'add a cell' ==")
    covers, cells_ok, empty_ok = 0, True, True
    for a in range(1, AB):
        for bz in range(1, AB - a + 1):
            l = a + bz
            words = []
            for pos in combinations(range(l), a):
                w = [0] * l
                for p in pos:
                    w[p] = 1
                if is_no_descent(w):
                    words.append(tuple(w))
            parts_w = set()
            for w in words:
                e, z = [], 0
                for bit in w:
                    if bit:
                        e.append(z)
                    else:
                        z += 1
                parts_w.add(tuple(e))
            # boundary g
            g = []
            for i in range(1, a + 1):
                gi = bz
                for j in range(1, bz + 1):
                    if f_thr(j) >= i:
                        gi = j - 1
                        break
                g.append(gi)
            parts_g = set()
            if f_thr(bz) <= a:
                def rec(i, lo, cur):
                    if i == a:
                        parts_g.add(tuple(cur))
                        return
                    for v in range(lo, g[i] + 1):
                        rec(i + 1, v, cur + [v])
                rec(0, 0, [])
            else:
                empty_ok &= not words
            cells_ok &= parts_w == parts_g
            C_of = {}
            for w in words:
                e, z = [], 0
                for bit in w:
                    if bit:
                        e.append(z)
                    else:
                        z += 1
                C_of[tuple(e)] = carry(w)
            for e, c in C_of.items():
                for i in range(a):
                    e2 = list(e)
                    e2[i] += 1
                    e2 = tuple(e2)
                    if e2 in C_of:
                        covers += 1
                        cells_ok &= C_of[e2] > c
                        # the increment is exactly 3^(a-i-1) 2^i 2^(e_i)  (i is 0-indexed here)
                        cells_ok &= C_of[e2] - c == qpow(3, a - i - 1) * (1 << i) * (1 << e[i])
    check("cell (a,b) == {partitions e <= g} for all a+b <= 12; empty iff f(b) > a", cells_ok and empty_ok)
    check("870 covering pairs, carry strictly larger by exactly 3^(a-i) 2^(i-1) 2^(e_i)", covers == 870, str(covers))
    res = {}
    for w in ["110111", "111011", "111101", "111110"]:
        wl = [int(c) for c in w]
        res[w] = (-carry(wl) * pow(3, -5, 64)) % 64
    print(f"  (5,1) residues: {res}")
    check("(5,1): residues 27, 39, 47, 31 for 110111, 111011, 111101, 111110 (THM-4503)", [res[w] for w in ["110111", "111011", "111101", "111110"]] == [27, 39, 47, 31])
    chain = ["111110", "111101", "111011", "110111"]
    check("chain 111110 < 111101 < 111011 < 110111 has residues 31, 47, 39, 27 and strictly increasing carries",
          [res[w] for w in chain] == [31, 47, 39, 27] and all(carry([int(c) for c in chain[i]]) < carry([int(c) for c in chain[i + 1]]) for i in range(3)))
    check("f(1) = 2 (the zero of a (5,1) word needs two ones before it)", f_thr(1) == 2)


# --------------------------------------------------------------------------------------
# F. the Sturmian element of E_inf
# --------------------------------------------------------------------------------------
def part_F(J: int = 20000):
    print(f"\n== F. upper mechanical word of slope log_3 2 (o_j = ceil(j log_3 2)) to j = {J} ==")
    pre, word, o = [(0, 0)], [], 0
    for j in range(1, J + 1):
        while qpow(3, o) <= (1 << j):
            o += 1
        word.append(1 if o > pre[-1][0] else 0)
        pre.append((o, j))
    check("o_j increases by 0 or 1 (a valid word) and all heights h_j > 0 (in E_inf)", all(pre[j][0] - pre[j - 1][0] in (0, 1) for j in range(1, J + 1)) and all(qpow(3, o_) > (1 << j_) for o_, j_ in pre[1:]))
    hf = [h_float(o_, j_) for o_, j_ in pre]
    print(f"  heights j>=1: min {min(hf[1:]):.6f} at j={hf.index(min(hf[1:]))}, max {max(hf):.6f}; log2 3 = {math.log2(3):.6f}")
    check("min height 0.000063, max 1.584621 (6 dp)", f"{min(hf[1:]):.6f}" == "0.000063" and f"{max(hf):.6f}" == "1.584621")
    check("all heights < log_2 3 (so 2^(-h_j) > 1/3: the real series diverges linearly)", all(h_lt(o_, j_, 1, 0) for o_, j_ in pre[1:]))
    fm = future_minima(pre)
    gaps = [fm[i + 1] - fm[i] for i in range(len(fm) - 1)]
    print(f"  window future minima: {fm[:6]} ... {fm[-6:]}; gaps (distinct) {sorted(set(gaps))}; first gaps {gaps[:5]}")
    check("time 0 IS a strict future minimum (h_0 = 0 < all later heights): the spine is {0}, not empty", fm[0] == 0)
    check("window future minima after 0 sit at multiples of 1054, then gaps 569/84/19/1 near the window end; distinct gaps {1,19,84,569,1054}", sorted(set(gaps)) == [1, 19, 84, 569, 1054] and gaps[0] == 1054)
    # continued fraction of log_3 2: denominators of the convergents above alpha and the semiconvergent 359/569
    fr = [(1, 1), (2, 3), (12, 19), (53, 84), (359, 569), (665, 1054)]
    above = all(qpow(3, p) > (1 << q) for p, q in fr)    # p/q > log_3 2  <=>  3^p > 2^q
    check("1/1, 2/3, 12/19, 53/84, 359/569, 665/1054 all lie above log_3 2 (one-sided approximations from above)", above)
    # children of the root 0: running minima of the walk on (0, inf) restricted to heights > 0 = all running minima
    runmin, cur = [], None
    for j in range(1, J + 1):
        if cur is None or h_lt(pre[j][0], j, cur[0], cur[1]):
            runmin.append(j)
            cur = pre[j]
    print(f"  children of the root 0 within the window (running height minima): {runmin}")
    check("root 0 has more than 5 children in the window (the tree at 0 is infinite but not locally finite)", len(runmin) >= 6)
    rs = []
    for k in (8, 16, 24, 32, 40, 48, 56, 64):
        w = word[:k]
        C, o_ = carry(w), sum(w)
        r = (-C * pow(3, -o_, 1 << k)) % (1 << k)
        rs.append((k, r, r / 2 ** k))
    print("  least residues of the cylinders: " + ", ".join(f"k={k}: {r} ({fr_:.3f} 2^k)" for k, r, fr_ in rs))
    check("cylinder residues k=8..64 at 0.980, 0.359, 0.744, 0.788, 0.671, 0.253, 0.005, 0.676 of 2^k",
          [f"{fr_:.3f}" for _, _, fr_ in rs] == ["0.980", "0.359", "0.744", "0.788", "0.671", "0.253", "0.005", "0.676"])
    S = sum(2.0 ** (-h) for h in hf)
    print(f"  sum_(j<=J) 2^(-h_j) = {S:.1f} (> J/3 = {J / 3:.1f})")


# --------------------------------------------------------------------------------------
# G. controls: 5x+1 orbit of 7 (drift), 5x+1 orbit of 5 (the Prop 1(c) counterexample), record orbits
# --------------------------------------------------------------------------------------
def embeds(s: str, t: str) -> bool:
    i = 0
    for ch in s:
        i = t.find(ch, i)
        if i < 0:
            return False
        i += 1
    return True


def part_G(steps: int = 3000):
    print(f"\n== G. controls ==")
    # --- 5x+1 orbit of 5: eventually periodic with a strict value future minimum at time 0
    vals, x = [5], 5
    for _ in range(40):
        x = Tmap(x, 5)
        vals.append(x)
    cyc = vals.index(13, 2)
    check("5x+1: orbit of 5 enters the cycle through 13 (5 -> 13 -> ... -> 13)", vals[1] == 13 and 13 in vals[2:])
    check("5x+1: orbit of 5 is eventually periodic, yet time 0 is a strict value future minimum (5 < every later value)", all(v > 5 for v in vals[1:]) and vals[cyc] == vals[1])
    print(f"  5x+1 orbit of 5: {vals[:10]} ... (cycle min 13 > 5): Prop 1(c) 'iff diverges' fails in the note's own T_b setting")
    # --- 5x+1 orbit of 7
    q = 5
    vals, pre, x, o = [7], [(0, 0)], 7, 0
    for j in range(1, steps + 1):
        if x & 1:
            o += 1
        x = Tmap(x, q)
        vals.append(x)
        pre.append((o, j))
    cfm = future_minima(pre, q)
    vfm = value_future_minima(vals)
    print(f"  5x+1 orbit of 7, {steps} steps: coefficient future minima {len(cfm)} (incl. 0 and {steps}? {0 in cfm}, {steps in cfm}), value future minima {len(vfm)}")
    check("5x+1: 512 coefficient future minima == 512 value future minima in the window", len(cfm) == 512 and len(vfm) == 512)
    check("5x+1: coefficient future minima (window) are value future minima (window)", all(j in set(vfm) for j in cfm))
    check("5x+1: the step out of every value future minimum is odd (value rises)", all(vals[j + 1] > vals[j] for j in vfm if j + 1 < len(vals)))
    gaps = [cfm[i + 1] - cfm[i] for i in range(len(cfm) - 1)]
    print(f"  block lengths: max {max(gaps)}, mean {sum(gaps) / len(gaps):.2f}, first 20 {gaps[:20]}")
    check("5x+1: block lengths max 174, mean 5.87", max(gaps) == 174 and f"{sum(gaps) / len(gaps):.2f}" == "5.87")
    bits = vals[-1].bit_length()
    print(f"  value at step {steps}: {bits} bits (log2 = {math.log2(vals[-1]):.3f})")
    check("value at step 3000 has 451 bits (log2 = 450.97; the .out's '451.0' is log2 rounded)", bits == 451, str(bits))
    # Theorem 3 (ii)/(iii) for q = 5 on blocks ending strictly inside the window (end < steps): forced ones, height, growth
    ones_ok = ht_ok = grow_ok = first_ok = True
    nblk = 0
    for i in range(len(cfm) - 1):
        a_, b_ = cfm[i], cfm[i + 1]
        if b_ >= steps:
            continue
        l = b_ - a_
        oo = pre[b_][0] - pre[a_][0]
        s_i, s_n = vals[a_], vals[b_]
        nblk += 1
        # 2^H s_i < s_(i+1) < 2^H (s_i + l/2),  2^H = 5^oo / 2^l  (exact: 5^oo s_i < s_n 2^l and 2 s_n 2^l < 5^oo (2 s_i + l))
        grow_ok &= qpow(5, oo) * s_i < s_n * (1 << l)
        grow_ok &= 2 * s_n * (1 << l) < qpow(5, oo) * (2 * s_i + l)
        if l >= 2:
            first_ok &= vals[a_ + 1] > vals[a_]
            ones_ok &= oo == a_of(l, 5)
            ht_ok &= qpow(5, oo - 1) < (1 << (l - 1))       # H < log2 5 - 1
    check(f"5x+1: {nblk} interior blocks: l >= 2 starts with an odd step, has ceil(l log_5 2) odd steps, H < log2 5 - 1", first_ok and ones_ok and ht_ok)
    check("5x+1: spine growth 2^H s_i < s_(i+1) < 2^H (s_i + l/2) on every interior block (Theorem 3(iii) shape)", grow_ok)
    # Higman probe on the first 300 spine points
    sp = [vals[j] for j in cfm[:300]]
    waits = []
    for i in range(len(sp)):
        bi = bin(sp[i])[2:]
        w = None
        for j in range(i + 1, len(sp)):
            if embeds(bi, bin(sp[j])[2:]):
                w = j - i
                break
        waits.append(w)
    found = [w for w in waits if w is not None]
    print(f"  Higman-good partners for {len(found)}/300 spine points; wait max {max(found)}, mean {sum(found) / len(found):.2f}, #wait=1: {sum(1 for w in found if w == 1)}")
    check("Higman: 165/300, max wait 141, mean 61.0, never 1", len(found) == 165 and max(found) == 141 and f"{sum(found) / len(found):.1f}" == "61.0" and all(w != 1 for w in found))

    # --- record orbits (T-map, 3x+1)
    print("  record orbits (T-map): n, steps, leaders, lower records, bound, forest roots == records, max depth, value future minima in window")
    exp_rec = {27: (8, 5, 16, 17), 703: (14, 10, 17, 15), 6171: (20, 13, 16, 15), 77031: (27, 17, 15, 14), 837799: (33, 20, 23, 23), 8400511: (35, 24, 27, 25), 63728127: (38, 26, 25, 26), 670617279: (44, 30, 19, 25)}
    all_ok = True
    for n, (nrec, bound, nlead, depth) in exp_rec.items():
        orb, x = [n], n
        while x != 1:
            x = Tmap(x)
            orb.append(x)
        pre = prefix_oj([v & 1 for v in orb[:-1]])
        lows = strict_lower_records(orb)
        highs = strict_upper_records(orb)
        par = excursion_parents(pre)
        roots = [j for j in range(len(par)) if par[j] == -1]
        dep = [0] * len(par)
        for j in range(len(par)):
            dep[j] = 0 if par[j] < 0 else dep[par[j]] + 1
        vfm = value_future_minima(orb)
        cfm = future_minima(pre)
        # minimal/maximal elements of the value-time poset by definition
        minimal = [i for i in range(len(orb)) if not any(orb[h] > orb[i] for h in range(i))]
        maximal = [j for j in range(len(orb)) if not any(orb[k] < orb[j] for k in range(j + 1, len(orb)))]
        ok = (len(lows), n.bit_length(), len(highs), max(dep)) == (nrec, bound, nlead, depth) and roots == lows and vfm == [len(orb) - 1] and minimal == highs and maximal == vfm and cfm == [len(orb) - 1]
        ok &= {orb[i].bit_length() - 1 for i in lows} == set(range(n.bit_length()))     # records hit every shell
        all_ok &= ok
        print(f"    {n:>10}: steps {len(orb) - 1:>4}, leaders {len(highs):>2}, records {len(lows):>2}, bound {n.bit_length():>2}, roots==records {roots == lows}, depth {max(dep):>2}, vfm {vfm == [len(orb) - 1]}, ok={ok}")
    check("record orbits: all eight lines of section 7 reproduced (records, bound, leaders, roots = records, depth, minimal = leaders, maximal = future minima)", all_ok)
    # --- excursion order from its bare definition (Prop 1(b)) on the orbits of 27 and 703 and the 5x+1 orbit of 7 (400 steps)
    forest_ok = True
    for (n, q, K) in [(27, 3, None), (703, 3, None), (7, 5, 250)]:
        orb, x = [n], n
        if K is None:
            while x != 1:
                x = Tmap(x, q)
                orb.append(x)
        else:
            for _ in range(K):
                x = Tmap(x, q)
                orb.append(x)
        pre = prefix_oj([v & 1 for v in orb[:-1]])
        M = len(pre)
        # rel[i][j]: i <= j and h_k > h_i for all k in (i, j]   (the bare definition, built incrementally in j)
        rel = [[False] * M for _ in range(M)]
        for i in range(M):
            rel[i][i] = True
            for j in range(i + 1, M):
                if h_lt(pre[i][0], pre[i][1], pre[j][0], pre[j][1], q):
                    rel[i][j] = True
                else:
                    break
        par = excursion_parents(pre, q)
        for j in range(M):
            anc = [i for i in range(j) if rel[i][j]]
            chain, p = [], par[j]
            while p >= 0:
                chain.append(p)
                p = par[p]
            forest_ok &= sorted(anc) == sorted(chain)                                   # ancestors = parent chain
            forest_ok &= all(rel[anc[s]][anc[s + 1]] for s in range(len(anc) - 1))    # ancestors form a chain
            kids = [k for k in range(j + 1, M) if par[k] == j]
            runmin, cur = [], None
            for k in range(j + 1, M):
                if cur is None or h_lt(pre[k][0], pre[k][1], cur[0], cur[1], q):
                    cur = pre[k]
                    if h_lt(pre[j][0], pre[j][1], pre[k][0], pre[k][1], q):
                        runmin.append(k)
            forest_ok &= kids == runmin                                                 # children = running minima above h_j
        roots = [j for j in range(M) if par[j] == -1]
        hl, cur = [], None
        for j in range(M):
            if cur is None or h_lt(pre[j][0], pre[j][1], cur[0], cur[1], q):
                hl.append(j)
                cur = pre[j]
        forest_ok &= roots == hl                                                        # roots = strict height lower records
        forest_ok &= all((not rel[i][j]) or (not rel[j][k]) or rel[i][k] for i in range(M) for j in range(i, M) if rel[i][j] for k in range(j, M))  # transitive
    check("excursion order from the bare definition: transitive; ancestors of every node form a chain = parent chain; roots = strict height lower records; children = running minima above h_j (orbits 27, 703; 5x+1 orbit of 7, 250 steps)", forest_ok)
    # --- section 6.5's modulus: at T-time f, x_f = 2^(-f) C_(w,f) mod 3^(o_f), not "mod 3^f"; orbit of 7 (coprime to 3;
    #     the literal formula holds iff 3^(f - o_f) | n, i.e. here only while the word is all ones, f <= 3)
    orb, x = [7], 7
    for _ in range(40):
        x = Tmap(x)
        orb.append(x)
    w = [v & 1 for v in orb[:-1]]
    mod_ok, wrong_f, n_small = True, [], 0
    for f in range(1, 41):
        o = sum(w[:f])
        C = carry(w[:f])
        mod_ok &= (orb[f] * (1 << f)) % qpow(3, o) == C % qpow(3, o)          # x_f = 2^(-f) C mod 3^(o_f)
        least_right = (C * pow(2, -f, qpow(3, o))) % qpow(3, o) if o > 0 else 0
        least_wrong = (C * pow(2, -f, qpow(3, f))) % qpow(3, f)
        if orb[f] < qpow(3, o):
            n_small += 1
            mod_ok &= least_right == orb[f]
        if least_wrong != orb[f]:
            wrong_f.append(f)
    check(f"orbit of 7: x_f = 2^(-f) C_(w,f) (mod 3^(o_f)) at every f <= 40 and equals the least residue at the {n_small} times with x_f < 3^(o_f); the literal 'mod 3^f' of section 6.5 gives x_f at no time f >= 4",
          mod_ok and n_small >= 30 and all(f in wrong_f for f in range(4, 41)), f"literal formula correct only at f in {sorted(set(range(1, 41)) - set(wrong_f))}")
    # value lower records == height lower records for n <= 10^4 (Terras along the chain)
    same = True
    for n in range(2, 10001):
        orb, x = [n], n
        while x != 1:
            x = Tmap(x)
            orb.append(x)
        pre = prefix_oj([v & 1 for v in orb[:-1]])
        hl, cur = [], None
        for j in range(len(pre)):
            if cur is None or h_lt(pre[j][0], pre[j][1], cur[0], cur[1]):
                hl.append(j)
                cur = pre[j]
        same &= hl == strict_lower_records(orb)
    check("value lower records == height lower records (forest roots) for all 2 <= n <= 10^4", same)
    # 3x+1: orbits reaching 1 have exactly one window future minimum (the final 1); a hypothetical cycle min c = 2 mod 3 would have the odd preimage (2c-1)/3 < c
    check("if a cycle minimum c is odd with c = 2 mod 3 then (2c-1)/3 is an odd integer below c mapping to c (so 'iff diverges' is unproven for 3x+1)",
          all(((2 * c - 1) % 3 == 0) and (((2 * c - 1) // 3) % 2 == 1) and Tmap((2 * c - 1) // 3) == c for c in range(5, 2000, 6)))


# --------------------------------------------------------------------------------------
def main():
    t0 = time.time()
    part_A()
    W, b = part_B()
    sigma, D, kappa = part_C()
    part_D(sigma, D, W)
    part_E()
    part_F()
    part_G()
    npass = sum(1 for _, ok in RESULTS if ok)
    nfail = sum(1 for _, ok in RESULTS if not ok)
    print(f"\n{npass} checks passed, {nfail} failed ({time.time() - t0:.1f} s)")
    for name, ok in RESULTS:
        if not ok:
            print(f"  FAILED: {name}")
    print("AUDIT SCRIPT DONE")


if __name__ == "__main__":
    main()
