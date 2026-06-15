#!/usr/bin/env python3
"""Moment-shadow versus compatibility wall.

monad-explorer-2026-06-15-S10.  Connects the new baby-Hodge moment-problem
reframe (T820 / OPEN-Q-099) to the local Faulhaber odd-moment packets
(HYP-2457/HYP-2458/HYP-2454) and the repunit atom packet (HYP-2529).

Two exact local computations:

1. Exhaustive n=6 baby-Hodge holes in the (c3,c5) plane admit explicit convex
   certificates inside a fixed-c3 fiber.  In particular, the flagship hole
   (c3,c5)=(8,10) sits between realized spectral points (8,8) and (8,11), so
   the simplest moment/flag relaxation permits it.  The obstruction therefore
   lives in the compatibility layer, not the one-particle moment cone.
2. The odd Faulhaber moments S_{2r+1}(n) form a Stieltjes moment sequence
   S_{2r+1}(n)=sum_{i=1..n} i*(i^2)^r, so their Hankel matrices are positive
   semidefinite.  Exact tower failure for p>=3 is therefore not a positivity
   failure either; it is again an integrality/compatibility wall.

Tournament Analysis is meaningful here because the script compares proof
carriers rather than raw tournaments.
"""

from __future__ import annotations

import itertools
from collections import defaultdict
from fractions import Fraction
from math import gcd


def all_tournaments(n: int):
    pairs = list(itertools.combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield bits, A


def matmul(A: list[list[int]], B: list[list[int]]) -> list[list[int]]:
    n = len(A)
    C = [[0] * n for _ in range(n)]
    for i in range(n):
        Ai = A[i]
        Ci = C[i]
        for k in range(n):
            aik = Ai[k]
            if aik == 0:
                continue
            Bk = B[k]
            for j in range(n):
                Ci[j] += aik * Bk[j]
    return C


def trace(A: list[list[int]]) -> int:
    return sum(A[i][i] for i in range(len(A)))


def cycle_counts(A: list[list[int]]) -> tuple[int, int]:
    A2 = matmul(A, A)
    A3 = matmul(A2, A)
    A4 = matmul(A2, A2)
    A5 = matmul(A4, A)
    return trace(A3) // 3, trace(A5) // 5


def score_sequence(A: list[list[int]]) -> tuple[int, ...]:
    return tuple(sorted(sum(row) for row in A))


def hamiltonian_paths(A: list[list[int]]) -> int:
    n = len(A)
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            cur = row[v]
            if cur == 0:
                continue
            av = A[v]
            for w in range(n):
                if (mask >> w) & 1 == 0 and av[w]:
                    dp[mask | (1 << w)][w] += cur
    return sum(dp[-1])


def hole_certificates(n: int):
    fibers: dict[int, dict[int, list[dict[str, object]]]] = defaultdict(lambda: defaultdict(list))
    for bits, A in all_tournaments(n):
        c3, c5 = cycle_counts(A)
        H = hamiltonian_paths(A)
        record = {
            "bits": bits,
            "scores": score_sequence(A),
            "H": H,
            "D": (H - 1 - 2 * (c3 + c5)) // 4,
        }
        fibers[c3][c5].append(record)

    certificates = []
    for c3 in sorted(fibers):
        realized = sorted(fibers[c3])
        for lower, upper in zip(realized, realized[1:]):
            gap = upper - lower
            if gap <= 1:
                continue
            for hole in range(lower + 1, upper):
                lam = Fraction(upper - hole, gap)
                mu = Fraction(hole - lower, gap)
                certificates.append(
                    {
                        "c3": c3,
                        "hole": hole,
                        "lower": lower,
                        "upper": upper,
                        "lambda": lam,
                        "mu": mu,
                        "lower_rep": fibers[c3][lower][0],
                        "upper_rep": fibers[c3][upper][0],
                    }
                )
    return certificates


def bareiss_det(M: list[list[int]]) -> int:
    n = len(M)
    if n == 0:
        return 1
    A = [row[:] for row in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        pivot = A[k][k]
        if pivot == 0:
            swap = None
            for i in range(k + 1, n):
                if A[i][k] != 0:
                    swap = i
                    break
            if swap is None:
                return 0
            A[k], A[swap] = A[swap], A[k]
            sign *= -1
            pivot = A[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * pivot - A[i][k] * A[k][j]) // prev
        prev = pivot
        for i in range(k + 1, n):
            A[i][k] = 0
    return sign * A[-1][-1]


def odd_power_sum(n: int, r: int) -> int:
    power = 2 * r + 1
    return sum(i**power for i in range(1, n + 1))


def hankel_matrix_odd_moments(n: int, size: int) -> list[list[int]]:
    return [[odd_power_sum(n, i + j) for j in range(size)] for i in range(size)]


def leading_hankel_dets(n: int, max_size: int) -> list[int]:
    return [bareiss_det(hankel_matrix_odd_moments(n, size)) for size in range(1, max_size + 1)]


def proof_carrier_tournament():
    carriers = [
        "spectral_convex_witness",
        "alpha2_compatibility_wall",
        "faulhaber_hankel_positivity",
        "repunit_prime_atom_supply",
    ]
    features = {
        "spectral_convex_witness": (3, 2, 2, 2),
        "alpha2_compatibility_wall": (3, 3, 3, 2),
        "faulhaber_hankel_positivity": (2, 3, 2, 3),
        "repunit_prime_atom_supply": (2, 2, 2, 3),
    }
    # Observable = majority over:
    # (exact local certificate, cross-thread reach, sharpness of obstruction,
    #  scalable proof leverage).
    edges = {u: set() for u in carriers}
    for i, u in enumerate(carriers):
        for v in carriers[i + 1 :]:
            wins_u = sum(a > b for a, b in zip(features[u], features[v]))
            wins_v = sum(a < b for a, b in zip(features[u], features[v]))
            if wins_u > wins_v:
                edges[u].add(v)
            elif wins_v > wins_u:
                edges[v].add(u)
            else:
                if carriers.index(u) < carriers.index(v):
                    edges[u].add(v)
                else:
                    edges[v].add(u)
    return carriers, edges


def score_histogram(carriers: list[str], edges: dict[str, set[str]]) -> dict[int, int]:
    hist: dict[int, int] = defaultdict(int)
    for u in carriers:
        hist[len(edges[u])] += 1
    return dict(sorted(hist.items()))


def directed_3cycles(carriers: list[str], edges: dict[str, set[str]]) -> int:
    count = 0
    for a, b, c in itertools.combinations(carriers, 3):
        if b in edges[a] and c in edges[b] and a in edges[c]:
            count += 1
        if c in edges[a] and a in edges[b] and b in edges[c]:
            count += 1
    return count


def scc_sizes(carriers: list[str], edges: dict[str, set[str]]) -> list[int]:
    reverse = {u: set() for u in carriers}
    for u in carriers:
        for v in edges[u]:
            reverse[v].add(u)

    seen = set()
    order = []

    def dfs1(u: str):
        seen.add(u)
        for v in edges[u]:
            if v not in seen:
                dfs1(v)
        order.append(u)

    for u in carriers:
        if u not in seen:
            dfs1(u)

    seen.clear()
    sizes = []

    def dfs2(u: str) -> int:
        seen.add(u)
        total = 1
        for v in reverse[u]:
            if v not in seen:
                total += dfs2(v)
        return total

    for u in reversed(order):
        if u not in seen:
            sizes.append(dfs2(u))
    return sorted(sizes, reverse=True)


def carrier_hamiltonian_paths(carriers: list[str], edges: dict[str, set[str]]) -> int:
    count = 0
    for perm in itertools.permutations(carriers):
        if all(perm[i + 1] in edges[perm[i]] for i in range(len(perm) - 1)):
            count += 1
    return count


def print_hole_certificate(cert: dict[str, object]) -> None:
    lower = cert["lower_rep"]
    upper = cert["upper_rep"]
    print(
        f"  hole (c3,c5)=({cert['c3']},{cert['hole']})"
        f" = {cert['lambda']}*({cert['c3']},{cert['lower']})"
        f" + {cert['mu']}*({cert['c3']},{cert['upper']})"
    )
    print(
        "    lower witness:"
        f" score={lower['scores']} H={lower['H']} D={lower['D']} bits={lower['bits']}"
    )
    print(
        "    upper witness:"
        f" score={upper['scores']} H={upper['H']} D={upper['D']} bits={upper['bits']}"
    )


def main():
    print("Moment shadow versus compatibility wall")
    print("HYP-2530 / T822 / OPEN-Q-101 candidate")
    print()

    print("[1] Exhaustive n=6 baby-Hodge holes and convex spectral certificates")
    certs = hole_certificates(6)
    for cert in certs:
        print_hole_certificate(cert)
    flagship = [c for c in certs if c["c3"] == 8 and c["hole"] == 10][0]
    print()
    print("  Flagship hole:")
    print("    (c3,c5)=(8,10) is spectrally moment-interior: it sits in the same")
    print("    c3-fiber between realized points (8,8) and (8,11).")
    print("    Therefore the one-particle spectral/moment shadow permits it; the")
    print("    obstruction must live in the compatibility/conflict layer (THM-499 D).")
    print()

    print("[2] Odd Faulhaber moments as a Stieltjes cone")
    print("  Exact identity: S_{2r+1}(n) = sum_{i=1..n} i*(i^2)^r")
    print("  So the odd-moment Hankel matrices are PSD; leading determinants are:")
    for n in range(1, 7):
        dets = leading_hankel_dets(n, 4)
        print(f"    n={n}: leading Hankel determinants {dets}")
    print("  Interpretation: the odd-moment lane already lives in a positive moment cone;")
    print("  the failure of exact p>=3 towers is again beyond positivity, at the")
    print("  Bernoulli/integrality/compatibility layer.")
    print()

    print("[3] Cross-thread synthesis")
    print("  Tournament side: c5=10 is a baby-Hodge hole that survives the simplest")
    print("  moment/convex relaxation.")
    print("  Faulhaber side: odd moments satisfy exact Hankel positivity automatically.")
    print("  Shared lesson: moment positivity is only the one-particle shadow.")
    print("  Real obstructions are packet/compatibility/integrality constraints.")
    print()

    print("[4] Proof-carrier Tournament Analysis")
    carriers, edges = proof_carrier_tournament()
    print("  vertices =", carriers)
    print("  observable = majority(exact_local_certificate, cross_thread_reach,")
    print("                       obstruction_sharpness, scalable_proof_leverage)")
    print("  tie path =", carriers)
    print("  score_hist =", score_histogram(carriers, edges))
    print("  directed_3cycles =", directed_3cycles(carriers, edges))
    print("  scc_sizes =", scc_sizes(carriers, edges))
    print("  hamiltonian_paths =", carrier_hamiltonian_paths(carriers, edges))
    print("  edges:")
    for u in carriers:
        print(f"    {u} -> {sorted(edges[u])}")
    print()

    print("[5] Next questions")
    print("  1. Upgrade the convex baby-Hodge certificates to actual flag-algebra or")
    print("     compatibility inequalities cutting the hole c5=10.")
    print("  2. Make the HYP-2458 packet idea explicit: which compatibility variable")
    print("     plays the role of alpha_2 / D on the Faulhaber side?")
    print("  3. Compare HYP-2529 prime-length atoms with moment-feasible but")
    print("     atom-forbidden overlap lengths.")


if __name__ == "__main__":
    main()
