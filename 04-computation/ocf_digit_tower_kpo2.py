#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
ocf_digit_tower_kpo2.py
=======================
THE 2-ADIC OCF DIGIT TOWER (HYP-2379, THM-466, tangent T007 upgrade).
Session: kind-pasteur-2026-06-10-S2, Thread B lab (subagent).

THEOREM CONTEXT (PROVED canon, THM-002 / CONJ-001; Grinberg-Stanley
arXiv:2307.05569 + Irving-Omar arXiv:2412.10572):
    H(T) = I(Omega(T), 2) = sum_k alpha_k 2^k
where Omega(T) is the conflict graph on DIRECTED ODD CYCLES of T (two cycles
adjacent iff they share a vertex) and alpha_k = number of collections of k
pairwise vertex-disjoint directed odd cycles (alpha_0 = 1).

DIGIT LEMMA (THM-466(i)):   H == sum_{k<m} alpha_k 2^k  (mod 2^m), all m.
FINITE IDENTITY (THM-466(ii)): a k-collection needs >= 3k vertices, so
    n <= 8  =>  H = 1 + 2*alpha_1 + 4*alpha_2   exactly,
    n <= 5  =>  H = 1 + 2*alpha_1               exactly.
REVERSAL INVARIANCE (THM-466(iii)): Omega(T) ~= Omega(T^op), so every
    alpha_k -- hence every 2-adic digit of H -- is reversal-invariant.

WHAT THIS SCRIPT DOES
  P0  LUT_3, LUT_5: directed-Hamiltonian-cycle counts of all 3- and
      5-tournaments by explicit rooted-permutation enumeration.
  P1  LUT_7 over all 2^21 7-tournaments by vectorized Held-Karp-type DP;
      cross-checked against brute force on random codes + exact total.
  P2  FULL census n = 3..7 (all 2^C(n,2) labeled tournaments): c3, c5, c7,
      alpha_1 = c3+c5+c7, alpha_2 (disjoint (3,3) and (3,5) pairs), and H by
      an INDEPENDENT Held-Karp path DP; ASSERT H == 1 + 2*alpha_1 + 4*alpha_2
      everywhere. Exact total checks for H, c3, c5, c7, alpha_2.
  P3  n = 8: 200,000 random tournaments, same assertion ((3,5) pairs live).
  P4  pure-Python brute-force cross-validation (explicit directed-odd-cycle
      lists + n!-permutation path scans) on random codes at every n.
  P5  T007 payload: distribution of v2(H-1) over all labeled tournaments
      n = 4..7; score-class analysis: is alpha_1 mod 2 (equivalently
      H mod 4) constant on every sorted-score class? Same test for c5, c7,
      c5+c7, alpha_2 parity and for the full valuation v2(H-1).
      Counterexamples printed explicitly if any test fails.

MISTAKE GUARDS HONORED
  - MISTAKE-001: nothing reused from old scripts; fresh code from definitions.
  - MISTAKE-023 (k>=5 trap): cycle counts are DIRECTED-CYCLE counts, never
    vertex-set counts. Every odd-subset contributes the number of directed
    Hamiltonian cycles of its induced subtournament (rooted enumeration/DP),
    so a 5-set carrying 3 distinct directed 5-cycles contributes 3.
  - MISTAKE-054: orientation handled by ONE audited arc_bit() used
    everywhere; LUT path cross-checked by brute force; the global
    H == 1+2*a1+4*a2 assertion is itself an end-to-end orientation check
    (H computed by an independent DP that never looks at cycles).
  - MISTAKE-028/036/055: each claim is labeled with its verified range.

BIT CONVENTION
  A tournament on [n] is an integer `code` of C(n,2) bits. Pairs (i<j) are
  indexed lexicographically: (0,1),(0,2),...,(0,n-1),(1,2),...
  Bit 1 at pair (i,j), i<j, means i -> j; bit 0 means j -> i.
  Induced subtournaments inherit the convention through the order-preserving
  relabeling, so sub-codes are consistent with the LUTs.

All arithmetic is exact integer arithmetic (numpy int32/int64, far below
overflow: max H at n=8 < 2^12, cycle-count products < 2^10).
"""

import sys
import time
from itertools import combinations, permutations

import numpy as np

RNG_SEED = 20260610
CHUNK = 1 << 17
N8_SAMPLE = 200_000

T0 = time.time()


def log(msg=""):
    print(f"[{time.time()-T0:7.1f}s] {msg}", flush=True)


# ---------------------------------------------------------------- bit layer
def pair_index_map(n):
    pid = {}
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            pid[(i, j)] = c
            c += 1
    return pid, c


def arc_bit(code, i, j, pid):
    """1 iff i -> j in tournament `code` (pure-python ints)."""
    if i < j:
        return (code >> pid[(i, j)]) & 1
    return 1 - ((code >> pid[(j, i)]) & 1)


def arc_arrays(n, codes, pid):
    """dict (i,j) -> int32 array, entry 1 iff i -> j. codes: int64 array."""
    a = {}
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if i < j:
                a[(i, j)] = ((codes >> pid[(i, j)]) & 1).astype(np.int32)
            else:
                a[(i, j)] = (1 - ((codes >> pid[(j, i)]) & 1)).astype(np.int32)
    return a


def subcode_vec(codes, S, pid):
    """Induced-subtournament codes on sorted tuple S, vectorized (int64)."""
    sub = np.zeros(codes.shape[0], dtype=np.int64)
    b = 0
    k = len(S)
    for x in range(k):
        for y in range(x + 1, k):
            sub |= ((codes >> pid[(S[x], S[y])]) & 1) << b
            b += 1
    return sub


# ------------------------------------------------------------- cycle layer
def build_lut_bruteforce(k):
    """LUT[code] = number of directed Hamiltonian cycles of k-tournament.

    Explicit DIRECTED-cycle enumeration (MISTAKE-023 guard): root at vertex 0,
    permute the remaining k-1 vertices, check all k arcs cyclically. Each
    directed Hamiltonian cycle is counted exactly once (unique rooted form).
    """
    pid, m = pair_index_map(k)
    lut = np.zeros(1 << m, dtype=np.int32)
    perms = list(permutations(range(1, k)))
    for code in range(1 << m):
        cnt = 0
        for p in perms:
            seq = (0,) + p
            ok = True
            for t in range(k):
                if not arc_bit(code, seq[t], seq[(t + 1) % k], pid):
                    ok = False
                    break
            if ok:
                cnt += 1
        lut[code] = cnt
    return lut


def ham_cycles_dp_vec(k, codes):
    """Directed Hamiltonian cycle counts of k-tournaments, vectorized DP.

    Rooted at vertex 0: dp[(S,last)] = # directed paths 0 -> ... -> last
    visiting {0} u S exactly; close with arc last -> 0.
    """
    pid, _ = pair_index_map(k)
    a = arc_arrays(k, codes, pid)
    cur = {}
    for v in range(1, k):
        cur[(1 << (v - 1), v)] = a[(0, v)].copy()
    for _ in range(k - 2):
        nxt = {}
        for (S, last), vec in cur.items():
            for w in range(1, k):
                wb = 1 << (w - 1)
                if S & wb:
                    continue
                key = (S | wb, w)
                contrib = vec * a[(last, w)]
                if key in nxt:
                    nxt[key] += contrib
                else:
                    nxt[key] = contrib
        cur = nxt
    out = np.zeros(codes.shape[0], dtype=np.int32)
    for (S, last), vec in cur.items():
        out += vec * a[(last, 0)]
    return out


def ham_paths_dp_vec(n, codes):
    """Directed Hamiltonian PATH counts (Held-Karp DP over (subset,last)).

    Completely independent of the cycle machinery: this is the H side of the
    identity. Each directed Hamiltonian path is counted once (as a sequence).
    """
    pid, _ = pair_index_map(n)
    a = arc_arrays(n, codes, pid)
    N = codes.shape[0]
    cur = {}
    for v in range(n):
        cur[(1 << v, v)] = np.ones(N, dtype=np.int32)
    for _ in range(n - 1):
        nxt = {}
        for (S, last), vec in cur.items():
            for w in range(n):
                wb = 1 << w
                if S & wb:
                    continue
                key = (S | wb, w)
                contrib = vec * a[(last, w)]
                if key in nxt:
                    nxt[key] += contrib
                else:
                    nxt[key] = contrib
        cur = nxt
    out = np.zeros(N, dtype=np.int32)
    for vec in cur.values():
        out += vec
    return out


# ------------------------------------------------------------ census layer
def census(n, codes, luts, chunk=CHUNK):
    """c3,c5,c7,a1,a2,H,score-key for every tournament in `codes` (int64)."""
    pid, _ = pair_index_map(n)
    N = codes.shape[0]
    out = {q: np.zeros(N, dtype=np.int32) for q in ("c3", "c5", "c7", "a2", "H")}
    keys = np.zeros(N, dtype=np.int64)
    subs3 = list(combinations(range(n), 3))
    subs5 = list(combinations(range(n), 5)) if n >= 5 else []
    subs7 = list(combinations(range(n), 7)) if n >= 7 else []
    pairs33 = [(A, B) for A, B in combinations(subs3, 2) if not set(A) & set(B)]
    pairs35 = [(A, B) for A in subs3 for B in subs5 if not set(A) & set(B)]
    for lo in range(0, N, chunk):
        hi = min(N, lo + chunk)
        ch = codes[lo:hi]
        cyc3 = {S: luts[3][subcode_vec(ch, S, pid)] for S in subs3}
        cyc5 = {S: luts[5][subcode_vec(ch, S, pid)] for S in subs5}
        cyc7 = {S: luts[7][subcode_vec(ch, S, pid)] for S in subs7}
        out["c3"][lo:hi] = sum(cyc3.values())
        if subs5:
            out["c5"][lo:hi] = sum(cyc5.values())
        if subs7:
            out["c7"][lo:hi] = sum(cyc7.values())
        a2 = np.zeros(hi - lo, dtype=np.int32)
        for A, B in pairs33:
            a2 += cyc3[A] * cyc3[B]
        for A, B in pairs35:
            a2 += cyc3[A] * cyc5[B]
        out["a2"][lo:hi] = a2
        out["H"][lo:hi] = ham_paths_dp_vec(n, ch)
        # sorted score key (base 16; scores < 16 for n <= 8)
        sc = np.zeros((hi - lo, n), dtype=np.int8)
        for i in range(n):
            s = np.zeros(hi - lo, dtype=np.int32)
            for j in range(n):
                if j == i:
                    continue
                if i < j:
                    s += ((ch >> pid[(i, j)]) & 1).astype(np.int32)
                else:
                    s += (1 - ((ch >> pid[(j, i)]) & 1)).astype(np.int32)
            sc[:, i] = s.astype(np.int8)
        sc.sort(axis=1)
        kk = np.zeros(hi - lo, dtype=np.int64)
        for t in range(n):
            kk = kk * 16 + sc[:, t]
        keys[lo:hi] = kk
    out["a1"] = out["c3"] + out["c5"] + out["c7"]
    out["key"] = keys
    out["n_pairs33"] = len(pairs33)
    out["n_pairs35"] = len(pairs35)
    return out


# --------------------------------------------------- brute-force validators
def brute_H(n, code, pid):
    cnt = 0
    for p in permutations(range(n)):
        ok = True
        for t in range(n - 1):
            if not arc_bit(code, p[t], p[t + 1], pid):
                ok = False
                break
        if ok:
            cnt += 1
    return cnt


def brute_odd_cycles(n, code, pid):
    """Explicit list of ALL directed odd cycles (one entry per directed cycle,
    stored as its vertex frozenset; multiplicity preserved). MISTAKE-023:
    distinct directed cycles on the same vertex set appear separately."""
    cycles = []
    for k in (3, 5, 7):
        if k > n:
            continue
        for S in combinations(range(n), k):
            v0 = S[0]
            for p in permutations(S[1:]):
                seq = (v0,) + p
                ok = True
                for t in range(k):
                    if not arc_bit(code, seq[t], seq[(t + 1) % k], pid):
                        ok = False
                        break
                if ok:
                    cycles.append(frozenset(S))
    return cycles


def brute_profile(n, code):
    pid, _ = pair_index_map(n)
    cycles = brute_odd_cycles(n, code, pid)
    a1 = len(cycles)
    a2 = 0
    for x in range(a1):
        for y in range(x + 1, a1):
            if not (cycles[x] & cycles[y]):
                a2 += 1
    by_len = {}
    for c in cycles:
        by_len[len(c)] = by_len.get(len(c), 0) + 1
    scores = sorted(
        sum(arc_bit(code, i, j, pid) for j in range(n) if j != i) for i in range(n)
    )
    return dict(
        H=brute_H(n, code, pid),
        a1=a1,
        a2=a2,
        c3=by_len.get(3, 0),
        c5=by_len.get(5, 0),
        c7=by_len.get(7, 0),
        scores=scores,
    )


def fmt_tournament(n, code):
    pid, _ = pair_index_map(n)
    lines = []
    for i in range(n):
        beats = [j for j in range(n) if j != i and arc_bit(code, i, j, pid)]
        lines.append(f"    {i} -> {beats}")
    return "\n".join(lines)


# --------------------------------------------------------------- v2 layer
def v2_code(H):
    """v2(H-1) clipped to 14; 15 encodes infinity (H == 1, transitive)."""
    x = H.astype(np.int64) - 1
    out = np.full(x.shape, 15, dtype=np.int64)
    nz = x != 0
    y = x[nz].copy()
    v = np.zeros(y.shape, dtype=np.int64)
    while np.any(y % 2 == 0):
        even = y % 2 == 0
        v[even] += 1
        y[even] >>= 1
    out[nz] = np.minimum(v, 14)
    return out


def v2_hist(H):
    vc = v2_code(H)
    vals, cts = np.unique(vc, return_counts=True)
    return {("inf" if int(v) == 15 else int(v)): int(c) for v, c in zip(vals, cts)}


# ----------------------------------------------------- score-class analysis
def score_class_nonconstant(key, parity):
    """Indices (into np.unique(key)) of score classes containing BOTH parities."""
    uniq, inv = np.unique(key, return_inverse=True)
    nc = uniq.shape[0]
    cnt = np.bincount(inv * 2 + parity.astype(np.int64), minlength=2 * nc)
    both = np.where((cnt[0::2] > 0) & (cnt[1::2] > 0))[0]
    return uniq, inv, both


def decode_key(k, n):
    digs = []
    for _ in range(n):
        digs.append(int(k % 16))
        k //= 16
    return tuple(reversed(digs))


def cyc3_py(code, T3, pid, lut3):
    i, j, k = T3
    sub = (
        arc_bit(code, i, j, pid)
        | (arc_bit(code, i, k, pid) << 1)
        | (arc_bit(code, j, k, pid) << 2)
    )
    return int(lut3[sub])


def find_triangle_reversal_pair(n, members, parities, lut3):
    """Within one score class: a pair differing by reversing ONE directed
    3-cycle (score-preserving move) with opposite alpha_1 parity."""
    pid, _ = pair_index_map(n)
    pdict = {int(c): int(p) for c, p in zip(members, parities)}
    for code in members[:5000]:
        code = int(code)
        for T3 in combinations(range(n), 3):
            if cyc3_py(code, T3, pid, lut3) == 1:
                i, j, k = T3
                mask = (
                    (1 << pid[(i, j)]) | (1 << pid[(j, k)]) | (1 << pid[(i, k)])
                )
                fl = code ^ mask
                if fl in pdict and pdict[fl] != pdict[code]:
                    return code, fl, T3
    return None


def exhibit_pair(n, code0, code1, label):
    p0 = brute_profile(n, code0)
    p1 = brute_profile(n, code1)
    print(f"  COUNTEREXAMPLE ({label}), n={n}:", flush=True)
    for tag, code, pr in (("A", code0, p0), ("B", code1, p1)):
        print(f"  tournament {tag}: code={code}  scores={pr['scores']}")
        print(fmt_tournament(n, code))
        x = pr["H"] - 1
        v2 = "inf" if x == 0 else (x & -x).bit_length() - 1
        print(
            f"    c3={pr['c3']} c5={pr['c5']} c7={pr['c7']} "
            f"alpha1={pr['a1']} alpha2={pr['a2']} H={pr['H']} "
            f"H mod 4 = {pr['H'] % 4}  v2(H-1) = {v2}"
        )
    assert p0["scores"] == p1["scores"], "score sequences differ!"
    return p0, p1


# ------------------------------------------------------------------- main
def main():
    summary = []
    rng = np.random.default_rng(RNG_SEED)

    # ---------------- P0: LUT_3, LUT_5
    log("P0: building LUT_3, LUT_5 by rooted-permutation brute force")
    lut3 = build_lut_bruteforce(3)
    lut5 = build_lut_bruteforce(5)
    assert int(lut3.sum()) == 2 and int(lut3.max()) == 1, "LUT_3 totals wrong"
    assert int(lut5.sum()) == 768, "LUT_5 total wrong (expect 24*2^5=768)"
    log(f"P0 ok: sum LUT_3 = {int(lut3.sum())} (=2), sum LUT_5 = {int(lut5.sum())}"
        f" (=768), max directed 5-cycles on one 5-set = {int(lut5.max())}")

    # ---------------- P1: LUT_7 (vectorized DP) + brute cross-check
    log("P1: building LUT_7 over all 2^21 7-tournaments (vectorized DP)")
    N7 = 1 << 21
    lut7 = np.zeros(N7, dtype=np.int32)
    for lo in range(0, N7, CHUNK):
        hi = min(N7, lo + CHUNK)
        lut7[lo:hi] = ham_cycles_dp_vec(7, np.arange(lo, hi, dtype=np.int64))
    tot7 = int(lut7.astype(np.int64).sum())
    assert tot7 == 720 * (1 << 14), f"LUT_7 total {tot7} != 11796480"
    pid7, _ = pair_index_map(7)
    perms7 = list(permutations(range(1, 7)))
    for code in rng.integers(0, N7, size=30):
        code = int(code)
        cnt = 0
        for p in perms7:
            seq = (0,) + tuple(p)
            ok = True
            for t in range(7):
                if not arc_bit(code, seq[t], seq[(t + 1) % 7], pid7):
                    ok = False
                    break
            if ok:
                cnt += 1
        assert cnt == int(lut7[code]), f"LUT_7 brute mismatch at code {code}"
    log(f"P1 ok: sum LUT_7 = {tot7} (= 720*2^14), 30/30 brute-force matches, "
        f"max directed 7-cycles on 7 vertices = {int(lut7.max())}")
    luts = {3: lut3, 5: lut5, 7: lut7}

    # ---------------- P2: full census n = 3..7
    results = {}
    log("P2: FULL census n=3..7, identity H == 1 + 2*a1 + 4*a2")
    for n in (3, 4, 5, 6, 7):
        pid, m = pair_index_map(n)
        codes = np.arange(1 << m, dtype=np.int64)
        r = census(n, codes, luts)
        results[n] = r
        H, a1, a2 = r["H"], r["a1"], r["a2"]
        # Redei machine-check
        assert np.all(H % 2 == 1), f"n={n}: H even somewhere (Redei violated?!)"
        bad = np.where(H != 1 + 2 * a1 + 4 * a2)[0]
        nfail = int(bad.size)
        # exact totals (independent of the DP/LUT plumbing)
        fact = 1
        for t in range(2, n + 1):
            fact *= t
        totH_exp = fact * (1 << (m - n + 1))
        c = lambda a, b: (
            0 if b > a else __import__("math").comb(a, b)
        )
        totc3_exp = c(n, 3) * (1 << m) // 4
        totc5_exp = c(n, 5) * 24 * (1 << (m - 5)) if n >= 5 else 0
        totc7_exp = c(n, 7) * 720 * (1 << (m - 7)) if n >= 7 else 0
        tota2_exp = r["n_pairs33"] * (1 << (m - 4)) if n >= 6 else 0
        totH = int(H.astype(np.int64).sum())
        totc3 = int(r["c3"].astype(np.int64).sum())
        totc5 = int(r["c5"].astype(np.int64).sum())
        totc7 = int(r["c7"].astype(np.int64).sum())
        tota2 = int(a2.astype(np.int64).sum())
        assert totH == totH_exp, f"n={n}: sum H {totH} != {totH_exp}"
        assert totc3 == totc3_exp, f"n={n}: sum c3 {totc3} != {totc3_exp}"
        assert totc5 == totc5_exp, f"n={n}: sum c5 {totc5} != {totc5_exp}"
        assert totc7 == totc7_exp, f"n={n}: sum c7 {totc7} != {totc7_exp}"
        if n <= 7:
            assert tota2 == tota2_exp, f"n={n}: sum a2 {tota2} != {tota2_exp}"
        ntrans = int(np.sum(H == 1))
        assert ntrans == fact, f"n={n}: #(H==1)={ntrans} != n!={fact}"
        if n <= 5:
            assert int(a2.max()) == 0, f"n={n}: alpha_2 nonzero below n=6!"
        log(
            f"  n={n}: {1<<m} tournaments, identity failures = {nfail}; "
            f"sum-checks H/c3/c5/c7/a2 ALL EXACT; #transitive = {ntrans} = n!; "
            f"H in [{int(H.min())},{int(H.max())}], a1 max={int(a1.max())}, "
            f"a2 max={int(a2.max())}"
        )
        summary.append(
            f"n={n}: H = 1+2a1+4a2 on ALL {1<<m} labeled tournaments, "
            f"{nfail} failures"
        )
        if nfail:
            for b in bad[:5]:
                print(f"    FAIL code={int(b)}: {brute_profile(n, int(b))}")

    # ---------------- P3: n = 8 random sample
    log(f"P3: n=8, {N8_SAMPLE} random tournaments (28-bit codes)")
    codes8 = rng.integers(0, 1 << 28, size=N8_SAMPLE, dtype=np.int64)
    r8 = census(8, codes8, luts)
    results[8] = r8
    H8, a18, a28 = r8["H"], r8["a1"], r8["a2"]
    assert np.all(H8 % 2 == 1)
    bad8 = np.where(H8 != 1 + 2 * a18 + 4 * a28)[0]
    log(
        f"  n=8: identity failures = {int(bad8.size)} / {N8_SAMPLE}; "
        f"sample means c3={r8['c3'].mean():.3f} (theory 14), "
        f"c5={r8['c5'].mean():.3f} (theory 42), "
        f"c7={r8['c7'].mean():.3f} (theory 45), "
        f"H={H8.mean():.2f} (theory 315); "
        f"#tournaments with (3,5)-pairs counted: a2 max={int(a28.max())}"
    )
    # how many sampled tournaments actually HAVE a (3,5) pair (recompute split)
    pid8, _ = pair_index_map(8)
    subs3_8 = list(combinations(range(8), 3))
    a2_35 = np.zeros(N8_SAMPLE, dtype=np.int32)
    for A in subs3_8:
        B = tuple(sorted(set(range(8)) - set(A)))
        a2_35 += luts[3][subcode_vec(codes8, A, pid8)] * luts[5][
            subcode_vec(codes8, B, pid8)
        ]
    log(
        f"  n=8: (3,5)-pair contribution: present in "
        f"{int(np.sum(a2_35 > 0))}/{N8_SAMPLE} samples, max {int(a2_35.max())} "
        f"-> the (3,5) digit channel is exercised"
    )
    summary.append(
        f"n=8: H = 1+2a1+4a2 on {N8_SAMPLE} random tournaments, "
        f"{int(bad8.size)} failures ((3,5) pairs active in "
        f"{int(np.sum(a2_35 > 0))} samples)"
    )

    # ---------------- P4: brute-force cross-validation
    log("P4: pure-python brute-force cross-validation on random codes")
    rngv = np.random.default_rng(RNG_SEED + 4)
    for n in (4, 5, 6, 7):
        pid, m = pair_index_map(n)
        ncheck = 12 if n < 7 else 8
        for code in rngv.integers(0, 1 << m, size=ncheck):
            code = int(code)
            pr = brute_profile(n, code)
            r = results[n]
            assert pr["H"] == int(r["H"][code]), (n, code, "H")
            assert pr["c3"] == int(r["c3"][code]), (n, code, "c3")
            assert pr["c5"] == int(r["c5"][code]), (n, code, "c5")
            assert pr["c7"] == int(r["c7"][code]), (n, code, "c7")
            assert pr["a1"] == int(r["a1"][code]), (n, code, "a1")
            assert pr["a2"] == int(r["a2"][code]), (n, code, "a2")
        log(f"  n={n}: {ncheck}/{ncheck} brute-force profiles match census")
    for t in range(4):
        code = int(codes8[t])
        pr = brute_profile(8, code)
        assert pr["H"] == int(H8[t]) and pr["a1"] == int(a18[t]) and pr["a2"] == int(
            a28[t]
        ), (8, code)
    log("  n=8: 4/4 brute-force profiles match (incl. explicit odd-cycle lists)")
    summary.append("brute-force validation: all profiles match (n=4..8)")

    # ---------------- P5: T007 payload
    log("P5a: distribution of v2(H-1) over ALL labeled tournaments")
    print()
    print("  v2(H-1) distribution (full labeled censuses; 'inf' = transitive):")
    for n in (4, 5, 6, 7):
        h = v2_hist(results[n]["H"])
        tot = sum(h.values())
        keys_sorted = sorted([k for k in h if k != "inf"]) + (
            ["inf"] if "inf" in h else []
        )
        row = ", ".join(f"v2={k}: {h[k]} ({h[k]/tot:.4f})" for k in keys_sorted)
        print(f"    n={n}: {row}", flush=True)
    print()

    log("P5b: score-class determinism tests (full labeled enumeration)")
    A000571 = {4: 4, 5: 9, 6: 22, 7: 59}
    verdicts = {}
    for n in (4, 5, 6, 7):
        r = results[n]
        key = r["key"]
        tests = {
            "a1 mod 2 (= H mod 4)": (r["a1"] % 2).astype(np.int64),
            "c3 mod 2": (r["c3"] % 2).astype(np.int64),
            "c5 mod 2": (r["c5"] % 2).astype(np.int64),
            "c7 mod 2": (r["c7"] % 2).astype(np.int64),
            "(c5+c7) mod 2": ((r["c5"] + r["c7"]) % 2).astype(np.int64),
            "a2 mod 2": (r["a2"] % 2).astype(np.int64),
        }
        uniq = np.unique(key)
        nclass = uniq.shape[0]
        assert nclass == A000571[n], f"n={n}: {nclass} score classes != A000571"
        print(f"    n={n} ({nclass} score classes):")
        for name, par in tests.items():
            u, inv, both = score_class_nonconstant(key, par)
            const = "CONSTANT on every class" if both.size == 0 else (
                f"NON-CONSTANT on {both.size}/{nclass} classes"
            )
            nodd = int(par.sum())
            print(f"      {name:24s}: {const}   (#odd = {nodd})", flush=True)
            verdicts[(n, name)] = (both.size, inv, u)
        # full valuation v2(H-1) score-determined?
        vc = v2_code(r["H"])
        u, inv = np.unique(key, return_inverse=True)
        cnt = np.bincount(inv * 16 + vc, minlength=16 * u.shape[0]).reshape(-1, 16)
        ncv2 = int(np.sum((cnt > 0).sum(axis=1) > 1))
        print(
            f"      {'v2(H-1) (full valuation)':24s}: "
            f"{'CONSTANT on every class' if ncv2 == 0 else f'NON-CONSTANT on {ncv2}/{u.shape[0]} classes'}",
            flush=True,
        )
        verdicts[(n, "v2")] = (ncv2, inv, u)

    # smallest counterexample for the headline question
    print()
    headline = "a1 mod 2 (= H mod 4)"
    smallest_fail = None
    for n in (4, 5, 6, 7):
        nbad, inv, uniq = verdicts[(n, headline)]
        if nbad > 0:
            smallest_fail = n
            break
    if smallest_fail is None:
        print("  VERDICT: alpha_1 mod 2 IS score-determined at n = 4,5,6,7 (full).")
        summary.append("alpha_1 mod 2 (= H mod 4) IS a score function at n<=7")
    else:
        n = smallest_fail
        r = results[n]
        par = (r["a1"] % 2).astype(np.int64)
        uniq, inv, both = score_class_nonconstant(r["key"], par)
        sizes = np.bincount(inv)
        k = int(both[np.argmin(sizes[both])])
        members = np.where(inv == k)[0].astype(np.int64)
        pmem = par[members]
        scoreseq = decode_key(int(uniq[k]), n)
        print(
            f"  VERDICT: alpha_1 mod 2 is NOT score-determined; smallest n = {n}."
        )
        print(
            f"  smallest offending class: scores {scoreseq}, size {members.size}, "
            f"parity split {int(np.sum(pmem == 0))}/{int(np.sum(pmem == 1))}"
        )
        code0 = int(members[pmem == 0].min())
        code1 = int(members[pmem == 1].min())
        exhibit_pair(n, code0, code1, "lexicographically smallest of each parity")
        tri = find_triangle_reversal_pair(n, members, pmem, lut3)
        if tri is not None:
            cA, cB, T3 = tri
            print(
                f"  MINIMAL MOVE: reversing the single directed 3-cycle on "
                f"vertices {T3} flips alpha_1 parity (hence H mod 4):"
            )
            exhibit_pair(n, cA, cB, f"one triangle reversal on {T3}")
            summary.append(
                f"alpha_1 mod 2 is NOT score-determined: smallest n = {n}, "
                f"single-triangle-reversal counterexample on scores {scoreseq}"
            )
        else:
            summary.append(
                f"alpha_1 mod 2 is NOT score-determined: smallest n = {n}, "
                f"class scores {scoreseq}"
            )

    # n=8 sampled score evidence (existence-only)
    par8 = (a18 % 2).astype(np.int64)
    u8, inv8, both8 = score_class_nonconstant(r8["key"], par8)
    print(
        f"  n=8 (sampled, existence-only): alpha_1 parity non-constant on "
        f"{both8.size}/{u8.shape[0]} sampled score classes", flush=True,
    )

    # is c5+c7 globally even? (the YES-branch formula hunt, tested regardless)
    print()
    for n in (5, 6, 7):
        r = results[n]
        nodd5 = int(np.sum(r["c5"] % 2 == 1))
        nodd7 = int(np.sum(r["c7"] % 2 == 1))
        nodd57 = int(np.sum((r["c5"] + r["c7"]) % 2 == 1))
        tot = r["H"].shape[0]
        print(
            f"    n={n}: #odd c5 = {nodd5}/{tot}, #odd c7 = {nodd7}/{tot}, "
            f"#odd (c5+c7) = {nodd57}/{tot}"
            + ("  -> c5+c7 ALWAYS EVEN" if nodd57 == 0 else
               "  -> c5+c7 NOT always even"),
            flush=True,
        )
    nodd57_8 = int(np.sum((r8["c5"] + r8["c7"]) % 2 == 1))
    print(f"    n=8 (sampled): #odd (c5+c7) = {nodd57_8}/{N8_SAMPLE}")

    print()
    log("==================== SUMMARY VERDICTS ====================")
    for s in summary:
        print(f"  * {s}", flush=True)
    log("done.")


if __name__ == "__main__":
    main()
