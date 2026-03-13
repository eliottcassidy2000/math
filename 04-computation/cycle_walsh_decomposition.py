#!/usr/bin/env python3
"""
cycle_walsh_decomposition.py -- Walsh decomposition of individual cycle counts

KEY QUESTION: The degree-2 Walsh of c_3 is ZERO. Does this hold for c_5, c_7?
If yes, then ALL degree-2 Walsh energy in H comes from cycle INTERACTIONS
(alpha_j for j >= 2 in the OCF), not from individual cycle counts.

This would mean: the overlap/conflict structure of Omega(T) is ENTIRELY
responsible for the degree-2 Walsh content of H.

PLAN:
1. Compute c_k(sigma) for each orientation sigma and each odd k
2. Compute Walsh transform of c_k at degree 2
3. Compute Walsh transform of alpha_j at degree 2
4. Verify Parseval: sum |h_hat_alpha_j|^2 * 4^j = |h_hat_H|^2 (degree 2)

Author: kind-pasteur-2026-03-12-S60
"""

from collections import defaultdict
from itertools import combinations
import math


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def classify_resonance(a, b, p):
    resonances = []
    for k in range(1, p):
        q = 2*k - 1
        if q >= p:
            break
        if (q*a - b) % p == 0:
            resonances.append((q, f"{q}a=b"))
        if (q*a + b) % p == 0:
            resonances.append((q, f"{q}a=-b"))
        if (a - q*b) % p == 0 and q != 1:
            resonances.append((q, f"a={q}b"))
        if (a + q*b) % p == 0 and q != 1:
            resonances.append((q, f"a=-{q}b"))
    return resonances


def count_kcycles_circulant(S, p, k):
    """Count directed k-cycles in circulant tournament with connection set S.
    Uses gap sequence enumeration: sum_{t=0}^{p-1} S_hat(t)^k / k
    where S_hat(t) = sum_{s in S} omega^{ts}.
    """
    import cmath
    omega = cmath.exp(2j * cmath.pi / p)
    S_hat = [sum(omega ** (t * s) for s in S) for t in range(p)]
    total = sum(s**k for s in S_hat)
    return round(total.real / k)


def count_ham_paths_fast(adj_list, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            cnt = dp[mask][v]
            if cnt == 0:
                continue
            for w in adj_list[v]:
                if not (mask & (1 << w)):
                    dp[mask | (1 << w)][w] += cnt
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def enumerate_odd_cycles(A, p, max_k):
    """Enumerate all directed odd cycles up to length max_k.
    Returns list of frozensets (vertex sets)."""
    cycles_by_len = defaultdict(list)
    for k in range(3, max_k + 1, 2):
        for subset in combinations(range(p), k):
            verts = list(subset)
            n_cyc = count_directed_ham_cycles(A, verts)
            for _ in range(n_cyc):
                cycles_by_len[k].append(frozenset(subset))
    return cycles_by_len


def count_directed_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        fwd = A[a][b] * A[b][c] * A[c][a]
        bwd = A[a][c] * A[c][b] * A[b][a]
        return fwd + bwd

    start = 0
    dp = {}
    dp[(1 << start, start)] = 1
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


def compute_alpha_decomp(cycles_all, p):
    """Compute alpha_j from list of cycle frozensets using backtracking."""
    n = len(cycles_all)
    if n > 25:
        return None  # too large

    # Build conflict graph adjacency
    nbr = [0] * n
    for i in range(n):
        for j in range(i+1, n):
            if cycles_all[i] & cycles_all[j]:
                nbr[i] |= (1 << j)
                nbr[j] |= (1 << i)

    alpha = [0] * (n + 1)
    def backtrack(v, mask, size):
        alpha[size] += 1
        for w in range(v + 1, n):
            if not (mask & (1 << w)):
                backtrack(w, mask | nbr[w], size + 1)
    backtrack(-1, 0, 0)
    return alpha


def main():
    print("=" * 70)
    print("CYCLE-COUNT WALSH DECOMPOSITION")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        half = 1 << m

        print(f"\n{'='*60}")
        print(f"p = {p}, m = {m}")
        print(f"{'='*60}")

        # For each orientation, compute c_k for k = 3, 5, 7, ..., p
        # Use Fourier method for speed
        ck_vals = {}  # ck_vals[k][bits] = c_k value
        H_vals = {}
        alpha_vals = {}  # alpha_vals[bits] = [alpha_0, alpha_1, ...]

        for bits in range(half):
            S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                       for i in range(m))

            # Compute cycle counts via Fourier
            for k in range(3, p + 1, 2):
                if k not in ck_vals:
                    ck_vals[k] = {}
                ck_vals[k][bits] = count_kcycles_circulant(S, p, k)

            # Compute H via bitmask DP
            adj_list = [[] for _ in range(p)]
            for v in range(p):
                for s in S:
                    adj_list[v].append((v + s) % p)
            H_vals[bits] = count_ham_paths_fast(adj_list, p)

            # Compute alpha decomposition (only for p=7 where it's feasible)
            if p <= 7:
                A = [[0]*p for _ in range(p)]
                for v in range(p):
                    for s in S:
                        A[v][(v + s) % p] = 1
                cycles_by_len = {}
                for k in range(3, p + 1, 2):
                    cycles_by_len[k] = []
                    for subset in combinations(range(p), k):
                        verts = list(subset)
                        nc = count_directed_ham_cycles(A, verts)
                        for _ in range(nc):
                            cycles_by_len[k].append(frozenset(subset))

                all_cycles = []
                for k in range(3, p + 1, 2):
                    all_cycles.extend(cycles_by_len[k])

                alpha = compute_alpha_decomp(all_cycles, p)
                if alpha:
                    alpha_vals[bits] = alpha

        # Show cycle count ranges
        print(f"\n  Cycle count ranges across orientations:")
        for k in range(3, p + 1, 2):
            vals = list(ck_vals[k].values())
            print(f"    c_{k}: min={min(vals)}, max={max(vals)}, "
                  f"mean={sum(vals)/len(vals):.2f}")

        # Compute degree-2 Walsh of each c_k
        print(f"\n  --- Degree-2 Walsh of c_k ---")
        for k in range(3, p + 1, 2):
            print(f"\n  c_{k} degree-2 Walsh:")
            any_nonzero = False
            q_groups = defaultdict(list)
            for a in range(m):
                for b in range(a + 1, m):
                    total = 0
                    for bits in range(half):
                        sa = 1 if bits & (1 << a) else -1
                        sb = 1 if bits & (1 << b) else -1
                        total += ck_vals[k][bits] * sa * sb
                    h_ck = total / half
                    if abs(h_ck) > 0.001:
                        any_nonzero = True
                        gap_a, gap_b = a + 1, b + 1
                        res = classify_resonance(gap_a, gap_b, p)
                        min_q = min(qq for qq, t in res) if res else 'inf'
                        chi_ab = legendre(gap_a * gap_b, p)
                        q_groups[min_q].append(((a, b), h_ck, chi_ab))

            if not any_nonzero:
                print(f"    ALL ZERO (c_{k} has no degree-2 Walsh content)")
            else:
                for q in sorted(q_groups.keys()):
                    items = q_groups[q]
                    mags = sorted(set(abs(v) for _, v, _ in items))
                    signs_ok = all((1 if v > 0 else -1) == c
                                  for _, v, c in items)
                    print(f"    q={q}: |h_hat_c{k}| = {mags}, "
                          f"sign=chi(ab): {'ALL OK' if signs_ok else 'FAILS'}")
                    for (a, b), v, c in items:
                        print(f"      ({a},{b}): h_hat_c{k}={v:>10.4f}, "
                              f"chi(ab)={c:+d}")

        # Compute degree-2 Walsh of H
        print(f"\n  --- Degree-2 Walsh of H ---")
        q_groups_H = defaultdict(list)
        for a in range(m):
            for b in range(a + 1, m):
                total = 0
                for bits in range(half):
                    sa = 1 if bits & (1 << a) else -1
                    sb = 1 if bits & (1 << b) else -1
                    total += H_vals[bits] * sa * sb
                h_H = total / half
                gap_a, gap_b = a + 1, b + 1
                res = classify_resonance(gap_a, gap_b, p)
                min_q = min(qq for qq, t in res) if res else 'inf'
                chi_ab = legendre(gap_a * gap_b, p)
                q_groups_H[min_q].append(((a, b), h_H, chi_ab))

        for q in sorted(q_groups_H.keys()):
            items = q_groups_H[q]
            mags = sorted(set(abs(v) for _, v, _ in items))
            print(f"    q={q}: |h_hat_H| = {mags}")

        # At p=7: verify alpha decomposition Walsh
        if p <= 7 and alpha_vals:
            print(f"\n  --- Alpha-decomposition Walsh (p={p}) ---")
            max_j = max(max(j for j in range(len(alpha))
                           if alpha[j] > 0)
                        for alpha in alpha_vals.values())

            for j in range(max_j + 1):
                print(f"\n    alpha_{j} degree-2 Walsh:")
                q_groups_a = defaultdict(list)
                for a in range(m):
                    for b in range(a + 1, m):
                        total = 0
                        for bits in range(half):
                            sa = 1 if bits & (1 << a) else -1
                            sb = 1 if bits & (1 << b) else -1
                            total += alpha_vals[bits][j] * sa * sb
                        h_aj = total / half
                        if abs(h_aj) > 0.001:
                            gap_a, gap_b = a + 1, b + 1
                            res = classify_resonance(gap_a, gap_b, p)
                            min_q = min(qq for qq, t in res) if res else 'inf'
                            chi_ab = legendre(gap_a * gap_b, p)
                            q_groups_a[min_q].append(((a, b), h_aj, chi_ab))

                if not q_groups_a:
                    print(f"      ALL ZERO")
                else:
                    for q in sorted(q_groups_a.keys()):
                        items = q_groups_a[q]
                        mags = sorted(set(abs(v) for _, v, _ in items))
                        print(f"      q={q}: |h_hat_alpha{j}| = {mags}")
                        for (a, b), v, c in items:
                            print(f"        ({a},{b}): val={v:>10.4f}, "
                                  f"chi(ab)={c:+d}, "
                                  f"sign={'OK' if (v>0)==(c>0) else 'FAIL'}")

            # Verify: H = sum 2^j * alpha_j
            print(f"\n    Verification: H = sum 2^j * alpha_j")
            for bits in range(half):
                H_check = sum(alpha_vals[bits][j] * (2**j)
                             for j in range(len(alpha_vals[bits])))
                if H_check != H_vals[bits]:
                    print(f"      MISMATCH at bits={bits}: H={H_vals[bits]}, "
                          f"check={H_check}")
                    break
            else:
                print(f"      ALL MATCH")

            # Degree-2 Walsh consistency: h_hat_H = sum 2^j * h_hat_alpha_j
            print(f"\n    Degree-2 consistency: h_hat_H = sum 2^j * h_hat_alpha_j")
            for a in range(m):
                for b in range(a + 1, m):
                    # h_hat_H
                    total_H = 0
                    for bits in range(half):
                        sa = 1 if bits & (1 << a) else -1
                        sb = 1 if bits & (1 << b) else -1
                        total_H += H_vals[bits] * sa * sb
                    h_H = total_H / half

                    # sum 2^j * h_hat_alpha_j
                    h_sum = 0
                    for j in range(max_j + 1):
                        total_aj = 0
                        for bits in range(half):
                            sa = 1 if bits & (1 << a) else -1
                            sb = 1 if bits & (1 << b) else -1
                            total_aj += alpha_vals[bits][j] * sa * sb
                        h_sum += (2**j) * total_aj / half

                    gap_a, gap_b = a + 1, b + 1
                    if abs(h_H - h_sum) > 0.001:
                        print(f"      MISMATCH ({a},{b}): h_H={h_H:.4f}, "
                              f"sum={h_sum:.4f}")
                    else:
                        print(f"      ({a},{b}): h_H={h_H:>8.4f} = "
                              f"sum 2^j * h_alpha_j = {h_sum:>8.4f} OK")

        # NEW: Degree-2 Walsh of total cycle count alpha_1
        print(f"\n  --- Degree-2 Walsh of alpha_1 = total cycles ---")
        total_cycles = {}
        for bits in range(half):
            total_cycles[bits] = sum(ck_vals[k][bits]
                                    for k in range(3, p + 1, 2))

        q_groups_tc = defaultdict(list)
        for a in range(m):
            for b in range(a + 1, m):
                total = 0
                for bits in range(half):
                    sa = 1 if bits & (1 << a) else -1
                    sb = 1 if bits & (1 << b) else -1
                    total += total_cycles[bits] * sa * sb
                h_tc = total / half
                if abs(h_tc) > 0.001:
                    gap_a, gap_b = a + 1, b + 1
                    res = classify_resonance(gap_a, gap_b, p)
                    min_q = min(qq for qq, t in res) if res else 'inf'
                    chi_ab = legendre(gap_a * gap_b, p)
                    q_groups_tc[min_q].append(((a, b), h_tc, chi_ab))

        if not q_groups_tc:
            print(f"    ALL ZERO (total cycle count has no deg-2 Walsh)")
        else:
            for q in sorted(q_groups_tc.keys()):
                items = q_groups_tc[q]
                mags = sorted(set(abs(v) for _, v, _ in items))
                signs_ok = all((1 if v > 0 else -1) == c for _, v, c in items)
                print(f"    q={q}: |h_hat_total| = {mags}, "
                      f"sign=chi(ab): {'ALL OK' if signs_ok else 'FAILS'}")
                for (a, b), v, c in items:
                    print(f"      ({a},{b}): val={v:>10.4f}, chi(ab)={c:+d}")


if __name__ == '__main__':
    main()
