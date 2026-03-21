#!/usr/bin/env python3
"""
overnight_cartan_s14.py -- kind-pasteur-2026-03-21-S14

OVERNIGHT SESSION: Unifying the Cartan Bridge

Three independently discovered principles must be reconciled:
1. SPECTRAL FLATNESS (S13): min tr(S^4) => max H among regular
2. CARTAN DECOUPLING (S94/95): [A,S]=0 for regular => max H
3. GRAND TRICHOTOMY (S94): gl(n,R) = INERT(2) + RAMIFIED(3) + SPLIT(7)

CRITICAL TEST: n=9
- d = (9-3)/4 = 1.5 NOT integer => NO DRT exists
- Regular tournaments have [A,S]=0 (Cartan decoupled) but NOT spectral flat
- This is where the two principles DIVERGE
- Does spectral flatness or Cartan decoupling better predict max H?

Also: Lie group structure, BCH composition, Trichotomy verification.

Author: kind-pasteur-2026-03-21-S14
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
import sys


def signed_adjacency(A):
    n = A.shape[0]
    S = 2.0 * A - 1.0 + np.eye(n)
    np.fill_diagonal(S, 0)
    return S


def count_ham_paths(A):
    n = A.shape[0]
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp.get((mask, v), 0) == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[(mask | (1 << w), w)] = dp.get((mask | (1 << w), w), 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def circulant_tournament(n, conn_set):
    """Build circulant tournament on Z_n with given connection set."""
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for d in conn_set:
            j = (i + d) % n
            A[i][j] = 1
    return A


def random_regular_tournament(n, num_switches=1000, seed=None):
    """Generate a random regular tournament via arc-switching from a circulant."""
    rng = np.random.RandomState(seed)
    # Start with a circulant regular tournament
    half = (n - 1) // 2
    conn = list(range(1, half + 1))
    A = circulant_tournament(n, conn)

    # Random switching: find two arcs i->j and k->l where i->l and k->j don't exist,
    # and swapping maintains regularity
    for _ in range(num_switches):
        i, j = rng.randint(0, n, 2)
        if i == j or not A[i][j]:
            continue
        k, l = rng.randint(0, n, 2)
        if k == l or not A[k][l]:
            continue
        if i == k or j == l or i == l or k == j:
            continue
        # Check if swap preserves scores
        # Swap: i->j becomes j->i, k->l becomes l->k
        # Score changes: i loses 1, j gains 1, k loses 1, l gains 1
        # Only works if we also swap another pair to compensate
        # Actually, a SCORE-PRESERVING switch: swap i->j and k->l for i->l and k->j
        # IF i->l doesn't exist (l->i exists) and k->j doesn't exist (j->k exists)
        if A[i][l] == 0 and A[k][j] == 0:
            # Swap: remove i->j, k->l; add i->l, k->j
            A[i][j] = 0; A[j][i] = 1
            A[k][l] = 0; A[l][k] = 1
            A[i][l] = 1; A[l][i] = 0
            A[k][j] = 1; A[j][k] = 0

    # Verify regularity
    scores = A.sum(axis=1)
    if not all(s == half for s in scores):
        # If switching broke regularity, start over with circulant
        return circulant_tournament(n, list(range(1, half + 1)))
    return A


# ============================================================
# PART A: The n=9 Critical Test
# ============================================================

def n9_critical_test():
    """
    At n=9, d=(9-3)/4=1.5 is NOT integer.
    No DRT exists. But regular tournaments still have [A,S]=0 (S_2=0).

    Question: Among regular n=9 tournaments, does min tr(S^4) => max H?
    This tests if the Spectral Flatness Principle survives beyond DRT.
    """
    print("=" * 70)
    print("PART A: THE n=9 CRITICAL TEST")
    print("=" * 70)
    print("d = (9-3)/4 = 1.5 => NO DRT exists at n=9")
    print("All regular tournaments have [A,S]=0 (Cartan decoupled)")
    print("But spectral flatness VARIES\n")

    n = 9

    # Generate all 16 circulant regular tournaments
    print("Circulant regular tournaments on Z_9:")
    pairs = [(1, 8), (2, 7), (3, 6), (4, 5)]
    circulants = []
    for mask in range(16):
        conn = []
        for bit, (a, b) in enumerate(pairs):
            if mask & (1 << bit):
                conn.append(a)
            else:
                conn.append(b)
        A = circulant_tournament(n, conn)
        # Verify regular
        scores = A.sum(axis=1)
        if not all(s == 4 for s in scores):
            continue

        H = count_ham_paths(A)
        S = signed_adjacency(A)
        tr_S4 = float(np.real(np.trace(np.linalg.matrix_power(S, 4))))
        tr_S6 = float(np.real(np.trace(np.linalg.matrix_power(S, 6))))

        # Check DRT: S^2 = -9I + J?
        S2 = S @ S
        target = -9 * np.eye(n) + np.ones((n, n))
        is_drt = np.max(np.abs(S2 - target)) < 0.01

        # Eigenvalues
        eigvals = np.linalg.eigvals(S)
        mags = sorted(set(round(abs(e), 3) for e in eigvals))

        # Commutator norm
        A_float = A.astype(float)
        A_anti = (A_float - A_float.T) / 2
        A_sym = (A_float + A_float.T) / 2
        scalar = np.trace(A_float) / n
        A_sym_tl = A_sym - scalar * np.eye(n)
        comm = A_anti @ A_sym_tl - A_sym_tl @ A_anti
        comm_norm = np.linalg.norm(comm, 'fro')

        circulants.append({
            'conn': tuple(conn), 'H': H, 'tr_S4': tr_S4,
            'tr_S6': tr_S6, 'is_drt': is_drt, 'mags': mags,
            'comm_norm': comm_norm
        })

    circulants.sort(key=lambda x: -x['H'])
    print(f"  {'conn':>15s} {'H':>6s} {'tr(S^4)':>10s} {'tr(S^6)':>12s} "
          f"{'DRT':>4s} {'|comm|':>8s} {'|eig| mags':>20s}")
    for c in circulants:
        print(f"  {str(c['conn']):>15s} {c['H']:6d} {c['tr_S4']:10.0f} "
              f"{c['tr_S6']:12.0f} {'Y' if c['is_drt'] else 'N':>4s} "
              f"{c['comm_norm']:8.4f} {str(c['mags']):>20s}")

    # Test spectral flatness principle
    max_H = max(c['H'] for c in circulants)
    min_tr4 = min(c['tr_S4'] for c in circulants)
    max_H_tours = [c for c in circulants if c['H'] == max_H]
    min_tr4_tours = [c for c in circulants if c['tr_S4'] == min_tr4]

    print(f"\n  Max H = {max_H} at conn = {[c['conn'] for c in max_H_tours]}")
    print(f"  Min tr(S^4) = {min_tr4} at conn = {[c['conn'] for c in min_tr4_tours]}")
    print(f"  SPECTRAL FLATNESS PRINCIPLE: "
          f"{'HOLDS' if max_H_tours[0]['tr_S4'] == min_tr4 else 'FAILS'}")

    # Also generate random regular tournaments
    print(f"\n  Random regular tournaments (50 samples):")
    random_results = []
    for seed in range(50):
        A = random_regular_tournament(n, num_switches=5000, seed=seed)
        scores = A.sum(axis=1)
        if not all(s == 4 for s in scores):
            continue
        H = count_ham_paths(A)
        S = signed_adjacency(A)
        tr_S4 = float(np.real(np.trace(np.linalg.matrix_power(S, 4))))
        random_results.append({'H': H, 'tr_S4': tr_S4})

    if random_results:
        by_H = defaultdict(list)
        for r in random_results:
            by_H[r['H']].append(r['tr_S4'])
        print(f"  {'H':>6s} {'count':>6s} {'mean tr(S^4)':>14s} {'range':>20s}")
        for H in sorted(by_H.keys(), reverse=True):
            vals = by_H[H]
            print(f"  {H:6d} {len(vals):6d} {np.mean(vals):14.1f} "
                  f"[{min(vals):.0f}, {max(vals):.0f}]")


# ============================================================
# PART B: Lie Group Structure — exp(t * S_T) periodicity
# ============================================================

def lie_group_analysis():
    """
    The exponential map exp(t * S_T) gives a curve in SO(n).
    For Paley: periodic with period 2*pi/sqrt(p).
    For non-Paley: quasi-periodic (incommensurate frequencies).

    Question: Does the periodicity/quasi-periodicity of exp(t*S_T)
    predict H or tournament structure?

    Also: the ORBIT under exp(t*S_T) action on R^n decomposes into
    circles of different radii. The radii are the eigenvalue magnitudes.
    """
    print("\n" + "=" * 70)
    print("PART B: LIE GROUP STRUCTURE — exp(t * S_T)")
    print("=" * 70)

    for n in [7, 9]:
        print(f"\nn = {n}:")

        # Collect regular tournaments
        if n == 7:
            # Exhaustive at n=7
            tournaments = []
            m = n * (n - 1) // 2
            for bits in range(1 << m):
                A = np.zeros((n, n), dtype=int)
                pos = 0
                for i in range(n):
                    for j in range(i + 1, n):
                        if bits & (1 << pos):
                            A[i][j] = 1
                        else:
                            A[j][i] = 1
                        pos += 1
                if all(A.sum(axis=1) == 3):
                    H = count_ham_paths(A)
                    tournaments.append((A, H))
                    if len(tournaments) >= 100:
                        break
            # Get one per H class
            by_H = {}
            for A, H in tournaments:
                if H not in by_H:
                    by_H[H] = A
        else:
            # Circulant at n=9
            by_H = {}
            pairs = [(1, 8), (2, 7), (3, 6), (4, 5)]
            for mask in range(16):
                conn = []
                for bit, (a, b) in enumerate(pairs):
                    conn.append(a if mask & (1 << bit) else b)
                A = circulant_tournament(n, conn)
                if all(A.sum(axis=1) == (n - 1) // 2):
                    H = count_ham_paths(A)
                    if H not in by_H:
                        by_H[H] = A

        for H in sorted(by_H.keys(), reverse=True):
            A = by_H[H]
            S = signed_adjacency(A)

            # Eigenvalues
            eigvals = np.linalg.eigvals(S)
            # Frequencies: omega_k = Im(lambda_k)
            freqs = sorted(set(round(abs(e.imag), 4) for e in eigvals if abs(e.imag) > 0.001))

            # Check periodicity: all frequencies commensurate?
            if len(freqs) == 1:
                period = 2 * np.pi / freqs[0]
                print(f"  H={H}: PERIODIC, omega={freqs[0]:.4f}, "
                      f"period={period:.4f}")
            elif len(freqs) > 1:
                # Check if ratios are rational
                ratios = [freqs[i] / freqs[0] for i in range(1, len(freqs))]
                # A ratio is "rational" if close to a ratio of small integers
                rational = True
                for r in ratios:
                    # Check if r = p/q with p,q < 20
                    found = False
                    for q in range(1, 20):
                        p = round(r * q)
                        if abs(r - p / q) < 0.001:
                            found = True
                            break
                    if not found:
                        rational = False

                print(f"  H={H}: {'QUASI-PERIODIC' if not rational else 'PERIODIC'}, "
                      f"freqs={freqs}, ratios={[round(r, 4) for r in ratios]}")
            else:
                print(f"  H={H}: TRIVIAL (no oscillation)")

            # Compute exp(t*S) at t = pi/sqrt(max_freq^2)
            if freqs:
                t = np.pi / freqs[-1]
                expS = np.real(np.linalg.matrix_power(
                    np.eye(n) + S * (t / 100), 100))  # crude approximation
                # Better: use scipy or eigendecomposition
                # exp(t*S) via eigendecomposition
                eigvals_full, V = np.linalg.eig(S)
                D = np.diag(np.exp(t * eigvals_full))
                expS = np.real(V @ D @ np.linalg.inv(V))
                det_exp = np.real(np.linalg.det(expS))
                print(f"         det(exp(pi*S/omega_max)) = {det_exp:.6f} "
                      f"(=1 for SO(n))")


# ============================================================
# PART C: BCH Tournament Composition
# ============================================================

def bch_composition():
    """
    When two tournaments T_1, T_2 are composed via BCH:
    exp(S_1) * exp(S_2) = exp(S_1 + S_2 + [S_1,S_2]/2 + ...)

    The commutator [S_1, S_2] is SYMMETRIC (Ghost 13).
    So BCH generates cooperation from tournament composition.

    Question: How much cooperation is generated? Does it scale
    with the "distance" between T_1 and T_2?
    """
    print("\n" + "=" * 70)
    print("PART C: BCH TOURNAMENT COMPOSITION")
    print("=" * 70)

    n = 5
    print(f"\nn={n}: Composing pairs of regular tournaments")

    # All regular n=5 tournaments
    regulars = []
    m = n * (n - 1) // 2
    for bits in range(1 << m):
        A = np.zeros((n, n), dtype=int)
        pos = 0
        for i in range(n):
            for j in range(i + 1, n):
                if bits & (1 << pos):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                pos += 1
        if all(A.sum(axis=1) == 2):
            regulars.append(A)

    print(f"  {len(regulars)} regular tournaments found")

    # For each pair, compute BCH corrections
    results = []
    for i in range(min(len(regulars), 10)):
        for j in range(i + 1, min(len(regulars), 10)):
            S1 = signed_adjacency(regulars[i])
            S2 = signed_adjacency(regulars[j])

            # First BCH correction: [S1, S2] / 2
            comm = (S1 @ S2 - S2 @ S1) / 2

            # Verify Ghost 13: [S1, S2] is symmetric
            asym_err = np.linalg.norm(comm - comm.T, 'fro')
            is_sym = asym_err < 0.001

            # How much cooperation is generated?
            coop_norm = np.linalg.norm(comm, 'fro')

            # Hamming distance between tournaments
            ham_dist = np.sum(regulars[i] != regulars[j]) // 2

            # H values
            H1 = count_ham_paths(regulars[i])
            H2 = count_ham_paths(regulars[j])

            results.append({
                'H1': H1, 'H2': H2, 'ham_dist': ham_dist,
                'coop_norm': coop_norm, 'is_sym': is_sym
            })

    # Ghost 13 verification
    all_sym = all(r['is_sym'] for r in results)
    print(f"\n  Ghost 13 ([S1,S2] is symmetric): "
          f"{'VERIFIED' if all_sym else 'FAILS'} "
          f"({sum(r['is_sym'] for r in results)}/{len(results)})")

    # Correlation between distance and cooperation
    if results:
        dists = [r['ham_dist'] for r in results]
        coops = [r['coop_norm'] for r in results]
        if len(set(dists)) > 1 and len(set(coops)) > 1:
            corr = np.corrcoef(dists, coops)[0, 1]
            print(f"  Corr(Hamming_distance, ||[S1,S2]||) = {corr:.4f}")
        print(f"\n  Sample pairs:")
        print(f"  {'H1':>4s} {'H2':>4s} {'dist':>5s} {'||[S1,S2]||':>12s}")
        for r in results[:10]:
            print(f"  {r['H1']:4d} {r['H2']:4d} {r['ham_dist']:5d} "
                  f"{r['coop_norm']:12.4f}")


# ============================================================
# PART D: Grand Trichotomy Verification
# ============================================================

def grand_trichotomy():
    """
    The Grand Cartan-Trichotomy Map:
    gl(n,R) = R*I (INERT/2) + so(n) (RAMIFIED/3) + p (SPLIT/7)

    Verify this at the COMPUTATIONAL level:
    1. The scalar sector controls parity (p=2): H mod 2 = 1 always
    2. The tournament sector controls cycles (p=3): c_3 determines 3-cycle structure
    3. The cooperation sector controls self-knowledge (p=7): H=7 is forbidden

    Also: compute mod-p reductions of each sector.
    """
    print("\n" + "=" * 70)
    print("PART D: GRAND TRICHOTOMY VERIFICATION")
    print("=" * 70)

    for n in [5, 7]:
        print(f"\nn={n}:")
        m = n * (n - 1) // 2

        # For each tournament, compute:
        # - H mod 2 (parity from scalar sector)
        # - c_3 mod 3 (cycles from tournament sector)
        # - H mod 7 (forbidden value from cooperation sector)

        parity_ok = 0
        c3_data = defaultdict(list)
        h_mod7 = defaultdict(int)
        total = 0

        limit = min(1 << m, 50000)
        np.random.seed(42)

        for trial in range(limit):
            if limit == 1 << m:
                bits = trial
            else:
                bits = np.random.randint(0, 1 << m)

            A = np.zeros((n, n), dtype=int)
            pos = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if bits & (1 << pos):
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    pos += 1

            H = count_ham_paths(A)
            total += 1

            # Parity
            if H % 2 == 1:
                parity_ok += 1

            # 3-cycles
            c3 = 0
            for triple in combinations(range(n), 3):
                a, b, c = triple
                if A[a][b] and A[b][c] and A[c][a]:
                    c3 += 1
                if A[a][c] and A[c][b] and A[b][a]:
                    c3 += 1
            c3_data[c3 % 3].append(H)

            # H mod 7
            h_mod7[H % 7] += 1

        print(f"  Parity (INERT/p=2): H odd for {parity_ok}/{total} "
              f"= {100*parity_ok/total:.1f}%"
              f"  {'REDEI CONFIRMED' if parity_ok == total else 'VIOLATION!'}")

        print(f"\n  Cycle structure (RAMIFIED/p=3):")
        print(f"  c_3 mod 3 distribution:")
        for r in sorted(c3_data.keys()):
            hs = c3_data[r]
            print(f"    c3 = {r} mod 3: {len(hs)} tournaments, "
                  f"mean H = {np.mean(hs):.1f}")

        print(f"\n  Self-knowledge (SPLIT/p=7):")
        print(f"  H mod 7 distribution:")
        for r in sorted(h_mod7.keys()):
            print(f"    H = {r} mod 7: {h_mod7[r]} tournaments "
                  f"({100*h_mod7[r]/total:.1f}%)")
        # Check: is H=7 actually forbidden?
        h7_count = sum(1 for trial in range(limit)
                       if count_ham_paths(np.zeros((1, 1))) == 7)
        # Actually just check from the data
        print(f"  H=7 is {'FORBIDDEN' if 7 not in [H for r in c3_data.values() for H in r] else 'ACHIEVABLE'} at n={n}")


# ============================================================
# PART E: Commutator-Spectral-H Triangle
# ============================================================

def commutator_spectral_triangle(n=7):
    """
    Three quantities form a triangle of invariants:
    1. ||[A_anti, A_sym_tl]||^2 = n*S_2/2 (commutator norm)
    2. tr(S_T^4) (spectral kurtosis)
    3. H (Hamiltonian path count)

    For regular tournaments: (1) = 0 always.
    So the triangle collapses to a LINE: tr(S^4) -- H.

    For general tournaments: all three are independent.
    Question: What does the triangle look like? Is H determined
    by the pair (commutator, kurtosis)?
    """
    print("\n" + "=" * 70)
    print(f"PART E: COMMUTATOR-SPECTRAL-H TRIANGLE (n={n})")
    print("=" * 70)

    m = n * (n - 1) // 2
    total = min(1 << m, 30000)

    results = []
    np.random.seed(42)

    for trial in range(total):
        if total == 1 << m:
            bits = trial
        else:
            bits = np.random.randint(0, 1 << m)

        A = np.zeros((n, n), dtype=int)
        pos = 0
        for i in range(n):
            for j in range(i + 1, n):
                if bits & (1 << pos):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                pos += 1

        H = count_ham_paths(A)
        S = signed_adjacency(A)
        tr_S4 = float(np.real(np.trace(np.linalg.matrix_power(S, 4))))

        # Commutator norm
        A_float = A.astype(float)
        scores = A_float.sum(axis=1)
        S2 = np.sum((scores - (n - 1) / 2) ** 2)
        comm_norm_sq = n * S2 / 2

        results.append({
            'H': H, 'tr_S4': tr_S4, 'comm_sq': comm_norm_sq, 'S2': S2
        })

    # Correlations
    Hs = [r['H'] for r in results]
    tr4s = [r['tr_S4'] for r in results]
    comms = [r['comm_sq'] for r in results]

    print(f"\n  Pairwise correlations ({total} tournaments):")
    print(f"    Corr(H, tr(S^4))     = {np.corrcoef(Hs, tr4s)[0,1]:.4f}")
    print(f"    Corr(H, ||[A,S]||^2) = {np.corrcoef(Hs, comms)[0,1]:.4f}")
    print(f"    Corr(tr(S^4), ||[A,S]||^2) = {np.corrcoef(tr4s, comms)[0,1]:.4f}")

    # Conditional analysis: within each score class, how does tr(S^4) predict H?
    by_S2 = defaultdict(list)
    for r in results:
        by_S2[round(r['S2'], 1)].append(r)

    print(f"\n  Within score classes:")
    print(f"  {'S2':>6s} {'count':>6s} {'H range':>12s} {'tr(S^4) range':>20s} "
          f"{'corr(H,tr4)':>12s}")
    for S2 in sorted(by_S2.keys()):
        group = by_S2[S2]
        if len(group) < 10:
            continue
        hs = [r['H'] for r in group]
        t4s = [r['tr_S4'] for r in group]
        if len(set(hs)) > 1 and len(set(t4s)) > 1:
            c = np.corrcoef(hs, t4s)[0, 1]
        else:
            c = float('nan')
        print(f"  {S2:6.1f} {len(group):6d} [{min(hs):3d}, {max(hs):3d}] "
              f"[{min(t4s):7.0f}, {max(t4s):7.0f}] {c:12.4f}")


# ============================================================
# PART F: The Spectral Flatness Index
# ============================================================

def spectral_flatness_index(n=7):
    """
    Define: F(T) = tr(S^2)^2 / (n * tr(S^4))

    F = 1 means maximally flat (DRT/Paley).
    F < 1 means spectrally peaked.

    For n=7 regular:
    - Paley: F = 42^2 / (7*294) = 1764/2058 = 0.857
    - H=175: F = 1764 / (7*742) = 0.340
    - H=171: F = 1764 / (7*486) = 0.518

    Does F predict H across ALL tournaments, not just regular?
    """
    print("\n" + "=" * 70)
    print(f"PART F: SPECTRAL FLATNESS INDEX (n={n})")
    print("=" * 70)

    m = n * (n - 1) // 2
    total = min(1 << m, 30000)

    results = []
    np.random.seed(42)

    for trial in range(total):
        if total == 1 << m:
            bits = trial
        else:
            bits = np.random.randint(0, 1 << m)

        A = np.zeros((n, n), dtype=int)
        pos = 0
        for i in range(n):
            for j in range(i + 1, n):
                if bits & (1 << pos):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                pos += 1

        H = count_ham_paths(A)
        S = signed_adjacency(A)
        tr2 = float(np.real(np.trace(S @ S)))
        tr4 = float(np.real(np.trace(np.linalg.matrix_power(S, 4))))

        F = tr2 ** 2 / (n * tr4) if abs(tr4) > 0.001 else 0

        results.append({'H': H, 'F': F, 'tr4': tr4})

    corr = np.corrcoef([r['H'] for r in results], [r['F'] for r in results])[0, 1]
    print(f"\n  Correlation(H, F) = {corr:.4f}")
    print(f"  (Positive means: flatter spectrum => more H)")

    # Bin by F
    bins = np.linspace(0, 1, 11)
    for i in range(len(bins) - 1):
        group = [r for r in results if bins[i] <= r['F'] < bins[i + 1]]
        if group:
            hs = [r['H'] for r in group]
            print(f"  F in [{bins[i]:.1f}, {bins[i+1]:.1f}): "
                  f"n={len(group):5d}, mean H={np.mean(hs):7.1f}, "
                  f"max H={max(hs):5d}")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    n9_critical_test()
    lie_group_analysis()
    bch_composition()

    if len(sys.argv) > 1 and sys.argv[1] == '--full':
        grand_trichotomy()
        commutator_spectral_triangle(n=5)
        spectral_flatness_index(n=5)

    commutator_spectral_triangle(n=7)
    spectral_flatness_index(n=7)
