#!/usr/bin/env python3
"""
walsh_arc_usage.py -- Connecting Walsh coefficients to arc usage in H-paths

QUESTION: Why does the degree-4 Walsh contribution flip sign at p mod 4?

APPROACH: Each Walsh coefficient h_hat[S] measures the effect of simultaneously
flipping chords in S. Flipping chord i replaces gap (i+1) with gap (p-i-1).

For a single chord flip, the H difference is:
  delta_H_i = H(sigma with sigma_i=+1) - H(sigma with sigma_i=-1)
            = h_hat[{i}] = 0 (by odd-degree vanishing)

Wait -- odd degree vanishes! So single chord flips have NO net effect on H
when averaged over all other chord choices. This is remarkable.

For a pair flip:
  Corr(sigma_i, sigma_j; H) = h_hat[{i,j}]
This measures whether chords i and j "cooperate" (+) or "compete" (-).

KEY IDEA: The Walsh coefficient h_hat[{i,j}] should be related to the
number of Hamiltonian paths that use BOTH arcs from chord i AND chord j.

If H-paths using both chords are more numerous than those using neither,
then flipping both simultaneously hurts H (negative coefficient).

This script computes:
1. Arc usage counts in all Hamiltonian paths for each chord
2. Pairwise arc co-usage counts
3. Comparison with Walsh coefficients

Author: kind-pasteur-2026-03-12-S59c
"""

import time
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths_with_arc_usage(A, n):
    """Count Hamiltonian paths AND track which arcs each path uses.

    Returns:
      total: total H-path count
      arc_count: arc_count[(u,v)] = number of H-paths using arc u->v
      pair_count: pair_count[((u1,v1),(u2,v2))] = paths using both arcs
    """
    # For small n, enumerate paths explicitly
    # For larger n, we need a different approach
    # Let's use Held-Karp with arc tracking at p=7

    if n > 9:
        # Too expensive for explicit enumeration
        return None, None, None

    # Enumerate all Hamiltonian paths by backtracking
    total = 0
    arc_count = defaultdict(int)
    # pair_count too expensive for large n, skip

    def backtrack(path, visited):
        nonlocal total
        if len(path) == n:
            total += 1
            # Record arcs used
            for idx in range(len(path)-1):
                arc_count[(path[idx], path[idx+1])] += 1
            return
        for w in range(n):
            if w in visited:
                continue
            if not path or A[path[-1]][w]:
                if path:  # check arc exists
                    visited.add(w)
                    path.append(w)
                    backtrack(path, visited)
                    path.pop()
                    visited.discard(w)
                else:
                    visited.add(w)
                    path.append(w)
                    backtrack(path, visited)
                    path.pop()
                    visited.discard(w)

    for start in range(n):
        backtrack([start], {start})

    return total, dict(arc_count), None


def chord_arc_usage(arc_count, p):
    """Group arc usage by chord (gap type).

    Chord i corresponds to gap (i+1). Each arc with gap (i+1) is
    of the form v -> v+(i+1) mod p.

    For a circulant tournament, arc_count[(v, (v+g) mod p)] should be
    the same for all v (by circulant symmetry). So we can define:

    chord_usage[g] = (1/p) * sum_v arc_count[(v, (v+g) mod p)]
    """
    chord_usage = {}
    for g in range(1, p):
        total = sum(arc_count.get((v, (v+g) % p), 0) for v in range(p))
        chord_usage[g] = total / p
    return chord_usage


def chord_pair_co_usage(A, p, n_sample=None):
    """Compute co-usage of chord pairs in Hamiltonian paths.

    For chords g1, g2: how many H-paths use an arc with gap g1 AND an arc with gap g2?

    For circulant: co_usage[g1, g2] is the average over all vertex pairs (v1, v2) of
    the number of H-paths using both v1->v1+g1 and v2->v2+g2.

    Actually, we want: for each H-path, define its "chord signature" as the multiset
    of gaps used. Then co_usage[g1,g2] = E[count(g1) * count(g2)] per H-path.
    """
    # Enumerate H-paths and track gap signatures
    if p > 9:
        return None

    co_usage = defaultdict(int)
    single_usage = defaultdict(int)
    total = 0

    def backtrack(path, visited, gap_count):
        nonlocal total
        if len(path) == p:
            total += 1
            # Record all gap pairs
            for g1 in gap_count:
                single_usage[g1] += gap_count[g1]
                for g2 in gap_count:
                    co_usage[(g1, g2)] += gap_count[g1] * gap_count[g2]
            return
        for w in range(p):
            if w in visited:
                continue
            if not path or A[path[-1]][w]:
                if path:
                    g = (w - path[-1]) % p
                    visited.add(w)
                    path.append(w)
                    gap_count[g] = gap_count.get(g, 0) + 1
                    backtrack(path, visited, gap_count)
                    gap_count[g] -= 1
                    if gap_count[g] == 0:
                        del gap_count[g]
                    path.pop()
                    visited.discard(w)
                else:
                    visited.add(w)
                    path.append(w)
                    backtrack(path, visited, gap_count)
                    path.pop()
                    visited.discard(w)

    for start in range(p):
        backtrack([start], {start}, {})

    # Normalize by H and p (circulant symmetry)
    # co_usage[g1,g2] is summed over all starting vertices and all paths
    # For circulant: normalized = co_usage[g1,g2] / total

    return total, dict(single_usage), dict(co_usage)


def walsh_from_H(p):
    """Compute Walsh decomposition from H values."""
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]

    H_dict = {}
    for bits in range(1 << m):
        sigma = []
        S = []
        for i, (a, b) in enumerate(pairs):
            if bits & (1 << i):
                S.append(a)
                sigma.append(+1)
            else:
                S.append(b)
                sigma.append(-1)
        S = sorted(S)
        A = build_adj(p, S)
        H = count_ham_paths_hk(A, p)
        H_dict[tuple(sigma)] = H

    # Walsh transform
    n = 1 << m
    h_hat = {}
    for bits in range(n):
        S_idx = [i for i in range(m) if bits & (1 << i)]
        coeff = 0
        for sigma_bits in range(n):
            sigma = [(1 if sigma_bits & (1 << i) else -1) for i in range(m)]
            H = H_dict[tuple(sigma)]
            prod = 1
            for i in S_idx:
                prod *= sigma[i]
            coeff += H * prod
        h_hat[tuple(S_idx)] = coeff / n

    return h_hat, H_dict


def count_ham_paths_hk(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def main():
    print("=" * 70)
    print("WALSH COEFFICIENTS AND ARC USAGE IN H-PATHS")
    print("=" * 70)

    for p in [7]:
        m = (p - 1) // 2
        S_int = list(range(1, m + 1))
        A = build_adj(p, S_int)

        print(f"\np={p}, m={m}, S_int={S_int}")

        # 1. Arc usage in Interval H-paths
        print(f"\n  1. ARC USAGE (Interval)")
        t0 = time.time()
        H, arc_count, _ = count_ham_paths_with_arc_usage(A, p)
        t1 = time.time()
        print(f"    H = {H} (enumerated in {t1-t0:.1f}s)")

        chord_usage = chord_arc_usage(arc_count, p)
        print(f"    Chord usage (avg per H-path per vertex):")
        for g in sorted(chord_usage):
            avg = chord_usage[g] / (H / p)  # normalize by H-paths per starting vertex
            print(f"      gap {g}: {chord_usage[g]:.0f} total, {avg:.4f} per path")

        # 2. Gap co-usage
        print(f"\n  2. GAP CO-USAGE (Interval)")
        t2 = time.time()
        H2, single, co_use = chord_pair_co_usage(A, p)
        t3 = time.time()
        print(f"    Verified H = {H2} ({t3-t2:.1f}s)")

        print(f"\n    Single gap usage (total across all paths & vertices):")
        for g in sorted(single):
            print(f"      gap {g}: {single[g]}")

        # Co-usage matrix
        print(f"\n    Co-usage matrix (normalized by H):")
        gaps_used = sorted(set(g for g in range(1, p) if single.get(g, 0) > 0))
        header = "     " + " ".join(f"g={g:>2}" for g in gaps_used)
        print(f"    {header}")
        for g1 in gaps_used:
            row = f"    g={g1:>2}"
            for g2 in gaps_used:
                val = co_use.get((g1, g2), 0) / H2 if H2 > 0 else 0
                row += f" {val:>5.2f}"
            print(row)

        # 3. Walsh coefficients
        print(f"\n  3. WALSH COEFFICIENTS")
        h_hat, H_dict = walsh_from_H(p)

        print(f"    Degree-2 coefficients:")
        for S, c in sorted(h_hat.items()):
            if len(S) == 2 and abs(c) > 0.01:
                i, j = S
                g_i = i + 1
                g_j = j + 1
                print(f"      h[{{{i},{j}}}] = {c:>+8.2f}  "
                      f"(chords with gaps {g_i},{g_j} or {p-g_i},{p-g_j})")

        # 4. Connection: Walsh coeff vs co-usage
        print(f"\n  4. WALSH-USAGE CONNECTION")
        print(f"    For chord pair (i,j), define:")
        print(f"      U_{i,j} = co-usage when both chords 'positive' (gap=i+1, j+1)")
        print(f"      V_{i,j} = co-usage when both chords 'negative' (gap=p-i-1, p-j-1)")
        print(f"      W_{i,j} = co-usage when mixed (one + one -)")
        print(f"    Walsh coeff h[{{i,j}}] should be related to (U+V-2W)/4")
        print()

        for S, c in sorted(h_hat.items()):
            if len(S) != 2 or abs(c) < 0.01:
                continue
            i, j = S
            g_pos_i = i + 1
            g_pos_j = j + 1
            g_neg_i = p - i - 1
            g_neg_j = p - j - 1

            # U = co-usage of (g_pos_i, g_pos_j) = paths with both positive gaps
            U = co_use.get((g_pos_i, g_pos_j), 0) / H2

            # V = co-usage of (g_neg_i, g_neg_j)
            V = co_use.get((g_neg_i, g_neg_j), 0) / H2

            # W = average of mixed: (g_pos_i, g_neg_j) and (g_neg_i, g_pos_j)
            W1 = co_use.get((g_pos_i, g_neg_j), 0) / H2
            W2 = co_use.get((g_neg_i, g_pos_j), 0) / H2
            W = (W1 + W2) / 2

            print(f"    chord pair ({i},{j}): h_hat = {c:>+6.2f}")
            print(f"      gaps: pos=({g_pos_i},{g_pos_j}), neg=({g_neg_i},{g_neg_j})")
            print(f"      U={U:.4f}, V={V:.4f}, W={W:.4f}")
            print(f"      (U+V-2W)/4 = {(U+V-2*W)/4:.4f}")
            print()

    # ====== DEEPER: CHORD PARITY AT p=11 ======
    print(f"\n{'='*70}")
    print("CHORD INTERACTION STRUCTURE AT p=11")
    print("=" * 70)

    p = 11
    m = (p - 1) // 2

    # Walsh coefficients
    h_hat, H_dict = walsh_from_H(p)

    # Analyze the structure of degree-2 coefficients
    print(f"\n  Degree-2 Walsh coefficients at p={p}:")
    print(f"    {'pair':>10} {'h_hat':>10} {'gap_pos':>10} {'gap_neg':>10}")
    for S, c in sorted(h_hat.items()):
        if len(S) != 2:
            continue
        if abs(c) < 0.01:
            continue
        i, j = S
        gp = f"({i+1},{j+1})"
        gn = f"({p-i-1},{p-j-1})"
        label = f"{{{i},{j}}}"
        print(f"    {label:>10} {c:>+10.2f} {gp:>10} {gn:>10}")

    # What determines the sign?
    # Look at the GAP DIFFERENCE |g_i - g_j| vs |g_i + g_j|
    print(f"\n  Gap difference analysis:")
    for S, c in sorted(h_hat.items()):
        if len(S) != 2 or abs(c) < 0.01:
            continue
        i, j = S
        g_i = i + 1
        g_j = j + 1
        diff = abs(g_i - g_j)
        summ = (g_i + g_j) % p
        print(f"    {{{i},{j}}}: gap ({g_i},{g_j}), |diff|={diff}, sum mod p={summ}, h_hat={c:>+.2f}")

    # What about the QR structure?
    # QR at p=11: {1,3,4,5,9}
    QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    NQR = set(range(1, p)) - QR
    print(f"\n  QR({p}) = {sorted(QR)}")
    print(f"  NQR({p}) = {sorted(NQR)}")

    print(f"\n  QR/NQR classification of chord pairs:")
    for S, c in sorted(h_hat.items()):
        if len(S) != 2 or abs(c) < 0.01:
            continue
        i, j = S
        g_i = i + 1
        g_j = j + 1
        qi = "QR" if g_i in QR else "NQR"
        qj = "QR" if g_j in QR else "NQR"
        # Also check (g_i * g_j) mod p
        prod = (g_i * g_j) % p
        qp = "QR" if prod in QR else "NQR"
        print(f"    {{{i},{j}}}: ({g_i}:{qi}, {g_j}:{qj}), product={prod}:{qp}, h_hat={c:>+.2f}")

    # ====== THE KEY: MULTIPLICATIVE CHARACTER CONNECTION ======
    print(f"\n{'='*70}")
    print("MULTIPLICATIVE CHARACTER CONNECTION")
    print("=" * 70)
    print("""
  The Walsh degree-2 coefficient h_hat[{i,j}] depends on whether
  the gap product g_i * g_j is a QR or NQR mod p.

  CONJECTURE: |h_hat[{i,j}]| has exactly 2 values:
    Large: when g_i * g_j mod p is in one residue class
    Small: when in the other

  For p=11:
    Large |h_hat| = 272.25 pairs: let's check their gap products
    Small |h_hat| = 8.25 pairs: let's check their gap products
""")

    # Check the conjecture
    large_prods = []
    small_prods = []
    for S, c in h_hat.items():
        if len(S) != 2 or abs(c) < 0.01:
            continue
        i, j = S
        g_i = i + 1
        g_j = j + 1
        prod = (g_i * g_j) % p
        if abs(abs(c) - 272.25) < 1:
            large_prods.append((g_i, g_j, prod, prod in QR))
        elif abs(abs(c) - 8.25) < 1:
            small_prods.append((g_i, g_j, prod, prod in QR))

    print(f"  Large |h| pairs: {[(g1,g2,p2,'QR' if qr else 'NQR') for g1,g2,p2,qr in large_prods]}")
    print(f"  Small |h| pairs: {[(g1,g2,p2,'QR' if qr else 'NQR') for g1,g2,p2,qr in small_prods]}")

    large_qr_count = sum(1 for _,_,_,qr in large_prods if qr)
    small_qr_count = sum(1 for _,_,_,qr in small_prods if qr)
    print(f"\n  Large: {large_qr_count}/{len(large_prods)} have QR product")
    print(f"  Small: {small_qr_count}/{len(small_prods)} have QR product")

    # ====== VERIFY AT p=13 ======
    print(f"\n{'='*70}")
    print("VERIFICATION AT p=13")
    print("=" * 70)

    p = 13
    m = (p - 1) // 2
    h_hat13, H_dict13 = walsh_from_H(p)

    QR13 = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    print(f"  QR({p}) = {sorted(QR13)}")

    print(f"\n  Degree-2 Walsh coefficients:")
    mags_13 = defaultdict(list)
    for S, c in sorted(h_hat13.items()):
        if len(S) != 2 or abs(c) < 0.01:
            continue
        i, j = S
        g_i = i + 1
        g_j = j + 1
        prod = (g_i * g_j) % p
        qr_label = "QR" if prod in QR13 else "NQR"
        mag = round(abs(c), 2)
        mags_13[mag].append((g_i, g_j, prod, qr_label, c > 0))
        print(f"    {{{i},{j}}}: gap ({g_i},{g_j}), prod={prod}:{qr_label}, h_hat={c:>+.2f}")

    print(f"\n  Magnitude grouping:")
    for mag in sorted(mags_13, reverse=True):
        items = mags_13[mag]
        qr_count = sum(1 for _,_,_,ql,_ in items if ql == "QR")
        pos_count = sum(1 for _,_,_,_,s in items if s)
        print(f"    |h|={mag}: {len(items)} pairs, {qr_count} QR, {pos_count} positive")
        for g1, g2, prod, ql, sign in items:
            print(f"      ({g1},{g2}): prod={prod}:{ql}, {'+'if sign else '-'}")


if __name__ == '__main__':
    main()
