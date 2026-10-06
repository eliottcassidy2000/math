#!/usr/bin/env python3
"""Owed audits (shadow theorem, segment ledger), the Beatty cone construction,
extended cone counts, the descent census to 2^28, basin-minima mass, rising
prefix counts versus the first-descent census, and the kernel-readout
independence check (bank-based readouts never certify outside the bank).

Companion note:
  05-knowledge/results/collatz_measurement_independence_20261005.md

Conventions as in collatz_refinement_floor_shadows_20261005.py: U(n)=oddpart(3n+1);
word w=(a_1..a_l) rising iff 2^A<3^l; carry c_w; cycle point x_w=c_w/(2^A-3^l);
class(w) = c_w * 2^(-A) mod 3^l; 2-adic shadow = x_w mod 2^(A+1).
W(n)=2K!(L+1)!/(L+K+2)! on rooted n; nu(m)=8/(3*4^bitlength(m)); p_j=lambda(6j+3).

Usage: python3 <this file> [--cone-depth 16] [--sieve-bits 28] [--head-bits 22]
"""
from fractions import Fraction as F
from math import factorial, lgamma, exp, log
import argparse
import json
import time

import numpy as np

try:
    from numba import njit
except Exception:  # pragma: no cover
    def njit(*a, **k):
        def wrap(f):
            return f
        return wrap if not a or not callable(a[0]) else a[0]

CHECKS = 0


def require(cond, witness=None):
    global CHECKS
    CHECKS += 1
    if not cond:
        raise RuntimeError(witness)


def v2(n):
    return (n & -n).bit_length() - 1


def U(n):
    t = 3 * n + 1
    return t >> v2(t)


def root_word(n):
    w = []
    while n != 1:
        t = 3 * n + 1
        a = v2(t)
        w.append(a)
        n = t >> a
    return w


def counters(word):
    if not word:
        return 0, 0
    return len(word) - 1, sum((a - 1) // 2 for a in word)


def wW(L, K):
    return F(2 * factorial(K) * factorial(L + 1), factorial(L + K + 2))


def carry(word):
    l = len(word)
    c = 0
    Apre = 0
    for i, a in enumerate(word):
        c += 3 ** (l - 1 - i) * (1 << Apre)
        Apre += a
    return c, Apre


def rising_words(max_len):
    out = []

    def rec(prefix, A, l):
        if l > 0 and (1 << A) < 3 ** l:
            out.append(tuple(prefix))
        if l == max_len:
            return
        for a in range(1, 64):
            if (1 << (A + a + (max_len - l - 1))) >= 3 ** max_len and (1 << (A + a)) >= 3 ** (l + 1):
                break
            rec(prefix + [a], A + a, l + 1)

    rec([], 0, 0)
    return out


# ---------------------------------------------------------------------------
# S1. Exhaustive-residue audit of the shadow theorem (independent of the formulas)
# ---------------------------------------------------------------------------
def follows(m, word):
    v = m
    for a in word:
        t = 3 * v + 1
        if v2(t) != a:
            return None
        v = t >> a
    return v


def backward(n, word):
    u = n
    for a in reversed(word):
        num = (1 << a) * u - 1
        if num % 3 != 0 or num <= 0:
            return None
        u = num // 3
    return u


def section_shadow_audit(max_len=7):
    print("== S1. Shadow theorem: exhaustive residue audit, words of length <= %d ==" % max_len)
    words = rising_words(max_len)
    fwd = bwd = 0
    for w in words:
        l = len(w)
        c, A = carry(w)
        d = (1 << A) - 3 ** l
        r2 = (c * pow(d % (1 << (A + 1)), -1, 1 << (A + 1))) % (1 << (A + 1))
        r3 = (c * pow(1 << A, -1, 3 ** l)) % 3 ** l       # the simpler class formula c*2^(-A)
        require(r3 == (c * pow(d % 3 ** l, -1, 3 ** l)) % 3 ** l, ("class formulas agree", w))
        # forward: every odd m below 4*2^(A+1) follows w exactly iff m = r2 mod 2^(A+1)
        for m in range(1, 4 * (1 << (A + 1)), 2):
            n = follows(m, w)
            require((n is not None) == (m % (1 << (A + 1)) == r2), ("forward", w, m))
            if n is not None:
                require(n == (3 ** l * m + c) >> A and n > m, ("forward value", w, m))
            fwd += 1
        # backward: every odd n below 6*3^l has an integral positive chain iff n = r3 mod 3^l,
        # except finitely many small n where the chain would hit a nonpositive integer
        for n in range(1, 6 * 3 ** l, 2):
            u = backward(n, w)
            inclass = (n % 3 ** l == r3)
            # integrality is exactly the class condition; positivity and descent are automatic
            require((u is not None) == inclass, ("backward", w, n))
            if u is not None:
                require(follows(u, w) == n and u < n, ("round trip", w, n))
            bwd += 1
    print("  %d rising words; %d forward residue checks; %d backward residue checks; no exceptions"
          % (len(words), fwd, bwd))


# ---------------------------------------------------------------------------
# S2. Direct audit of the segment-closure ledger (predecessor enumeration)
# ---------------------------------------------------------------------------
def nu(m):
    return F(8, 3 * 4 ** m.bit_length())


def section_ledger_audit(M=1 << 16, targets=1 << 12):
    print("== S2. Segment ledger: direct predecessor audit, starts below 2^%d ==" % (M.bit_length() - 1))
    from collections import defaultdict
    for name in ('stopping', 'single-rise'):
        Wsig = defaultdict(F)
        exits = defaultdict(F)
        through = defaultdict(F)
        for m in range(3, M, 2):
            v = m
            seg = [m]
            while True:
                t = 3 * v + 1
                a = v2(t)
                u = t >> a
                stop = (u < m or u == 1) if name == 'stopping' else (a >= 2 or u == 1)
                if stop:
                    exits[u] += nu(m)
                    break
                seg.append(u)
                v = u
            for p in seg:
                Wsig[p] += nu(m)
            for y in seg[1:]:
                through[y] += nu(m)
        # incoming sum at y by enumerating the keys p of Wsig with U(p) = y
        incoming = defaultdict(F)
        for p, val in Wsig.items():
            incoming[U(p)] += val
        bad = []
        for y in range(3, targets, 2):
            require(incoming[y] == through[y] + exits[y], (name, y))
            if Wsig[y] - incoming[y] != nu(y) - exits[y]:
                bad.append(y)
        require(not bad, (name, bad[:5]))
        print("  rule %-12s: K W_sigma(y) = through(y) + exits(y) and defect = nu(y) - exits(y)"
              " at %d targets (direct enumeration of %d recorded predecessors)" % (name, (targets - 3) // 2, len(Wsig)))


# ---------------------------------------------------------------------------
# S3. Beatty construction of primitive rising words (sufficiency)
# ---------------------------------------------------------------------------
def ceil_log23(j):
    """ceil(j log_2 3) = bit_length(3^j) (3^j is never a power of two)."""
    return (3 ** j).bit_length() if j > 0 else 0


def section_beatty(max_len=60):
    print("== S3. Beatty cones: primitive rising word exists at depth l iff {l log2 3} < log2(3/2) ==")
    exists = []
    for l in range(1, max_len + 1):
        cond = 1 + ceil_log23(l - 1) <= (3 ** l).bit_length() - 1   # 1 + ceil((l-1)a) <= floor(l a)
        if cond:
            # construct: a_1 = 1, suffix sums exactly ceil(j a)
            word = [1] + [0] * (l - 1)
            for j in range(1, l):
                word[l - j] = ceil_log23(j) - ceil_log23(j - 1)
            require(all(a in (1, 2) for a in word[1:]), (l, word))
            A = sum(word)
            require((1 << A) < 3 ** l, ("rising", l))
            for j in range(1, l):
                require((1 << sum(word[l - j:])) >= 3 ** j, ("suffix", l, j))
            exists.append(l)
        else:
            # necessity: any primitive word needs a_1 = 1 and 1 + ceil((l-1)a) <= floor(l a)
            pass
    print("  depths with primitive cones (constructive), l<=%d:" % max_len, exists)
    frac_ok = [l for l in range(1, max_len + 1) if 1 + ceil_log23(l - 1) <= (3 ** l).bit_length() - 1]
    require(frac_ok == exists)
    return exists


# ---------------------------------------------------------------------------
# S4. Extended cone counts (numba DFS) and the descent census to 2^28
# ---------------------------------------------------------------------------
@njit(cache=True)
def _powmod(b, e, m):
    r = 1
    b %= m
    while e > 0:
        if e & 1:
            r = (r * b) % m
        b = (b * b) % m
        e >>= 1
    return r


@njit(cache=True)
def _mark_cones(L, pow3, own):
    """DFS over words of length <= L that can still rise by length L; mark own[l][class]."""
    # explicit stack: arrays of (depth, carry, A, next letter)
    maxd = L + 1
    st_c = np.zeros(maxd + 1, dtype=np.int64)
    st_A = np.zeros(maxd + 1, dtype=np.int64)
    st_a = np.zeros(maxd + 1, dtype=np.int64)
    depth = 0
    st_c[0] = 0
    st_A[0] = 0
    st_a[0] = 1
    count = 0
    while depth >= 0:
        if depth == L:
            depth -= 1
            continue
        a = st_a[depth]
        st_a[depth] += 1
        A = st_A[depth] + a
        # prune: minimal all-ones extension must rise by length L: 2^(A + L - (depth+1)) < 3^L
        if (A + L - depth - 1) >= 64 or (1 << (A + L - depth - 1)) >= pow3[L]:
            depth -= 1
            continue
        c = 3 * st_c[depth] + (1 << st_A[depth])
        l = depth + 1
        if (1 << A) < pow3[l]:
            inv2 = (pow3[l] + 1) // 2
            cls = (c % pow3[l]) * _powmod(inv2, A, pow3[l]) % pow3[l]
            own[l][cls] = 1
            count += 1
        depth += 1
        st_c[depth] = c
        st_A[depth] = A
        st_a[depth] = 1
    return count


@njit(cache=True)
def _sieve_D(limit):
    m_odd = (limit + 1) // 2
    mark = np.zeros(m_odd, dtype=np.uint8)
    for i in range(m_odd):
        m = 2 * i + 1
        v = m
        while True:
            t = 3 * v + 1
            while (t & 1) == 0:
                t >>= 1
            v = t
            if v < m or v == 1:
                break
            if v < limit:
                mark[(v - 1) // 2] = 1
    return mark


def section_cones_and_census(L, sieve_bits):
    print("== S4. Cone classes to depth %d (numba) and the descent census to 2^%d ==" % (L, sieve_bits))
    t0 = time.time()
    pow3 = np.array([3 ** i for i in range(L + 1)], dtype=np.int64)
    own = [np.zeros(1, dtype=np.uint8)] + [np.zeros(3 ** l, dtype=np.uint8) for l in range(1, L + 1)]
    from numba.typed import List
    own_t = List()
    for arr in own:
        own_t.append(arr)
    nwords = _mark_cones(L, pow3, own_t)
    cover_prev = np.zeros(1, dtype=np.uint8)
    series = []
    new_counts = []
    cum = F(0)
    for l in range(1, L + 1):
        cover = np.array(own_t[l]) | np.tile(cover_prev, 3) if l > 1 else np.array(own_t[1])
        new = int(cover.sum()) - 3 * int(cover_prev.sum()) if l > 1 else int(cover.sum())
        new_counts.append(new)
        cum += F(new, 3 ** l)
        series.append(float(cum))
        cover_prev = cover
    print("  rising words visited: %d (%.1fs); new classes per depth:" % (nwords, time.time() - t0), new_counts)
    print("  cone density series:", " ".join("%d:%.5f" % (l + 1, s) for l, s in enumerate(series)))
    # depth-12 values must reproduce the previous note (1,1,0,1,0,2,8,0,28,0,124,602)
    ref = [1, 1, 0, 1, 0, 2, 8, 0, 28, 0, 124, 602]
    require(new_counts[:12] == ref[:len(new_counts)], new_counts[:12])
    # zeros exactly where {l log2 3} >= log2(3/2)
    for l in range(1, L + 1):
        cond = 1 + ceil_log23(l - 1) <= (3 ** l).bit_length() - 1
        require((new_counts[l - 1] > 0) == cond, ("beatty", l, new_counts[l - 1]))
    t0 = time.time()
    mark = _sieve_D(1 << sieve_bits)
    dens = float(mark.mean())
    n = 2 * np.arange(len(mark)) + 1
    dens_units = float(mark[n % 3 != 0].mean())
    blocks = [(b, float(mark[(n >= (1 << b)) & (n < (1 << (b + 1)))].mean())) for b in range(sieve_bits - 4, sieve_bits)]
    print("  descent set D below 2^%d (%.1fs): density %.6f (units %.6f); last blocks:"
          % (sieve_bits, time.time() - t0, dens, dens_units),
          " ".join("[2^%d): %.5f" % (b, d) for b, d in blocks))
    require(dens >= series[-1] - 1e-9)
    return dict(new_counts=new_counts, series=series, density=dens, density_units=dens_units), mark


# ---------------------------------------------------------------------------
# S5. W-mass on basin minima (units with no smaller ancestor)
# ---------------------------------------------------------------------------
@njit(cache=True)
def _counters_head(limit):
    m = (limit + 1) // 2
    Ls = np.empty(m, dtype=np.int64)
    Ks = np.empty(m, dtype=np.int64)
    for i in range(m):
        n = 2 * i + 1
        if n == 1:
            Ls[i] = 0
            Ks[i] = 0
            continue
        L = -1
        K = 0
        while n != 1:
            t = 3 * n + 1
            a = 0
            while (t & 1) == 0:
                t >>= 1
                a += 1
            K += (a - 1) // 2
            L += 1
            n = t
        Ls[i] = L
        Ks[i] = K
    return Ls, Ks


def section_minima_mass(head_bits, markD):
    print("== S5. W-mass by class below 2^%d: leaves / descent set / basin minima ==" % head_bits)
    limit = 1 << head_bits
    Ls, Ks = _counters_head(limit)
    n = 2 * np.arange(len(Ls)) + 1
    L = Ls.astype(np.float64)
    K = Ks.astype(np.float64)
    lg = np.vectorize(lgamma)
    Wv = np.exp(log(2) + lg(K + 1) + lg(L + 2) - lg(L + K + 3))
    D = markD[:len(Ls)].astype(bool)
    leaf = (n % 3 == 0)
    unit = ~leaf
    minima = unit & ~D & (n > 1)
    tot = Wv.sum()
    parts = {
        'root': float(Wv[n == 1].sum()),
        'leaves (lambda head)': float(Wv[leaf].sum()),
        'descent set D': float(Wv[D].sum()),
        'basin minima (units, no smaller ancestor)': float(Wv[minima].sum()),
    }
    for k, v in parts.items():
        print("  %-44s %.6f  (%.2f%% of head %.6f)" % (k, v, 100 * v / tot, tot))
    require(abs(sum(parts.values()) - tot) < 1e-9)
    cnt = dict(units=int(unit.sum()), minima=int(minima.sum()), D=int(D.sum()))
    print("  counts: units %d, basin minima %d (%.4f of units), D %d" % (cnt['units'], cnt['minima'],
                                                                      cnt['minima'] / cnt['units'], cnt['D']))
    return parts, cnt


# ---------------------------------------------------------------------------
# S6. Rising parity prefixes versus the first-descent census
# ---------------------------------------------------------------------------
def rising_prefix_counts(Amax):
    """N(A) = number of parity words b_1..b_A (b_1=1) with 3^{k_j} > 2^j for all prefixes j."""
    # dp over (j, k): number of words of length j with k ones, all prefixes rising
    counts = []
    dp = {1: 1}  # after j=1: k=1 (3 > 2)
    counts.append(sum(dp.values()))
    for j in range(2, Amax + 1):
        new = {}
        for k, cnt in dp.items():
            for bit in (0, 1):
                k2 = k + bit
                if 3 ** k2 > (1 << j):
                    new[k2] = new.get(k2, 0) + cnt
        dp = new
        counts.append(sum(dp.values()))
    return counts


@njit(cache=True)
def _stopping_T(limit):
    """first-descent time in T-steps for odd n < limit (0 if never within cap)."""
    m = (limit + 1) // 2
    out = np.zeros(m, dtype=np.int32)
    for i in range(m):
        n = 2 * i + 1
        if n == 1:
            continue
        v = n
        j = 0
        while True:
            if v & 1:
                v = (3 * v + 1) >> 1
            else:
                v >>= 1
            j += 1
            if v < n:
                break
        out[i] = j
    return out


def section_prefix_profile(head_bits=24, Amax=120):
    print("== S6. Rising parity prefixes N(A) versus the exact first-descent census below 2^%d ==" % head_bits)
    N = rising_prefix_counts(Amax)
    hstar = -(log(2) / log(3)) * log(log(2) / log(3)) / log(2) - (1 - log(2) / log(3)) * log(1 - log(2) / log(3)) / log(2)
    print("  h* = H(log_3 2) = %.6f; log2 N(A)/A at A=30,60,120: %.4f %.4f %.4f"
          % (hstar, log(N[29], 2) / 30, log(N[59], 2) / 60, log(N[119], 2) / 120))
    sig = _stopping_T(1 << head_bits)
    odd_count = len(sig) - 1
    rows = []
    for A in (4, 8, 12, 16, 20, 23, 30, 40, 60, 80, 100):
        actual = int((sig > A).sum())
        predicted = N[A - 1] / 2 ** (A - 1) * odd_count if A - 1 < len(N) else float('nan')
        rows.append((A, actual, predicted))
    print("  A : #{n<2^%d odd, sigma_T(n) > A} vs N(A) 2^-(A-1) * (odd count)" % head_bits)
    for A, actual, pred in rows:
        print("   %3d : %10d  vs %14.1f   ratio %.4f" % (A, actual, pred, actual / pred if pred else float('nan')))
    # for A <= head_bits-1 the two agree up to least representatives
    for A, actual, pred in rows:
        if A <= head_bits - 1:
            require(abs(actual - pred) <= 0.02 * pred + 50, (A, actual, pred))
    print("  max first-descent time below 2^%d: %d T-steps" % (head_bits, int(sig.max())))
    return rows


# ---------------------------------------------------------------------------
# S7. Kernel readouts from a certified bank never certify outside the bank
# ---------------------------------------------------------------------------
def section_readout(bank_bits=18):
    print("== S7. Codex kernel readout A_(m,d) = 9H_(d+1) - 8H_d from a certified bank ==")
    limit = 1 << bank_bits
    Ls, Ks = _counters_head(limit)
    n = 2 * np.arange(len(Ls)) + 1
    leaf = (n % 3 == 0)
    L = Ls.astype(np.float64)
    K = Ks.astype(np.float64)
    lg = np.vectorize(lgamma)
    Wv = np.exp(log(2) + lg(K + 1) + lg(L + 2) - lg(L + K + 3))
    idx = ((n[leaf] - 3) // 6).astype(np.float64)   # indices j of bank leaves
    pj = Wv[leaf]
    certified_mass = pj.sum()
    u = 1.0 - certified_mass   # exact total mass is 1 (Codex P5), so u bounds the uncertified mass
    print("  bank: rooted leaves below 2^%d, certified mass %.6f, uncertified/unenumerated mass u = %.6f"
          % (bank_bits, certified_mass, u))

    bank_set = set(idx.astype(int).tolist())
    j_max = int(idx.max())

    def readouts(m, d):
        # h = 4t/(1+t)^2 with t = 2^(m-j) equals 1/cosh^2((m-j) ln2 / 2); cosh overflow gives h = 0 exactly
        with np.errstate(over='ignore'):
            h = 1.0 / np.cosh((m - idx) * (log(2) / 2)) ** 2
        c = float((pj * h ** (d + 1)).sum())            # H_(d+1) >= bank part (omitted atoms are >= 0)
        if m in bank_set:
            with np.errstate(over='ignore'):
                sup_out = float(1.0 / np.cosh((m - (j_max + 1)) * (log(2) / 2)) ** 2)  # h decreases away from m
        else:
            sup_out = 1.0                                # the target itself is outside the bank
        b = float((pj * h ** d).sum()) + u * sup_out ** d   # H_d <= bank part + unknown mass * sup
        return 9 * c - 8 * b

    # targets outside the bank: the first few leaves above the bank limit
    targets_out = [((limit + 3) // 6) + k for k in range(0, 40, 13)]
    for m in targets_out:
        vals = [readouts(m, d) for d in (0, 5, 10, 20, 40)]
        require(all(v < 0 for v in vals), (m, vals))
    print("  targets outside the bank (indices %s): readout 9c-8b < 0 for d in {0,5,10,20,40} (never certified)" % targets_out)
    # a target inside the bank: readout turns positive once d is large (Codex (4)), using the bank only
    m_in = int(((27 - 3) // 6))   # the leaf 27, inside the bank
    vals = [readouts(m_in, d) for d in (0, 10, 20, 40, 80, 120)]
    print("  target 27 inside the bank: readouts at d=0,10,20,40,80,120:", ["%.2e" % v for v in vals])
    print("  (lambda(27) = %.3e; the readout needs 8(16/25)^d below the atom: d >~ %d)"
          % (float(wW(*counters(root_word(27)))), int(log(8 / float(wW(*counters(root_word(27))))) / log(25 / 16)) + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cone-depth', type=int, default=16)
    ap.add_argument('--sieve-bits', type=int, default=28)
    ap.add_argument('--head-bits', type=int, default=22)
    ap.add_argument('--audit-len', type=int, default=7)
    ap.add_argument('--json', type=str, default='')
    args = ap.parse_args()
    t0 = time.time()
    section_shadow_audit(args.audit_len)
    section_ledger_audit()
    exists = section_beatty(60)
    cones, markD = section_cones_and_census(args.cone_depth, args.sieve_bits)
    parts, cnt = section_minima_mass(args.head_bits, markD)
    rows = section_prefix_profile(24, 120)
    section_readout(18)
    print("== Summary ==")
    print("  checks: %d, total time %.1fs" % (CHECKS, time.time() - t0))
    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(dict(checks=CHECKS, beatty_depths=exists, cones=cones, mass_parts=parts,
                           counts=cnt, prefix_rows=rows,
                           status="PROVED scoped; FINITE-EXACT audits and census; VERIFIED readouts"),
                      fh, indent=1, default=str)
        print("  json written:", args.json)


if __name__ == '__main__':
    main()
