#!/usr/bin/env python3
"""S645: perfect numbers as fixed points of the aliquot/divisor carrier.

Carrier:

    divisor lattice D(n) -> sigma(n) -> s(n)=sigma(n)-n.

Perfect numbers are fixed points of the aliquot map s(n)=n, equivalently
sigma(n)=2n or abundancy sigma(n)/n=2.  The useful side channel is the divisor
factorization product:

    sigma(n)/n = prod_{p^a || n} (1 + 1/p + ... + 1/p^a).

The script scans a small finite window, records fixed points and sociable
cycles, prices the near-fixed defect, and runs Tournament Analysis on carrier
lenses rather than on raw integers.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt


LIMIT = 100_000
TIE_PATH = [
    "divisor_lattice_fixed_point",
    "sigma_product_abundancy",
    "aliquot_function_graph",
    "euclid_euler_mersenne_section",
    "local_prime_power_ledger",
    "odd_perfect_side_channel",
    "near_fixed_defect_scout",
    "raw_sequence_list",
]


def sigma_sieve(limit: int) -> list[int]:
    sigma = [0] * (limit + 1)
    for d in range(1, limit + 1):
        for m in range(d, limit + 1, d):
            sigma[m] += d
    return sigma


def factor(n: int) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            out.append((d, e))
        d += 1 if d == 2 else 2
    if n > 1:
        out.append((n, 1))
    return out


def fmt_factor(factors: list[tuple[int, int]]) -> str:
    pieces = []
    for p, e in factors:
        pieces.append(str(p) if e == 1 else f"{p}^{e}")
    return " * ".join(pieces) if pieces else "1"


def sigma_from_factor(factors: list[tuple[int, int]]) -> int:
    total = 1
    for p, e in factors:
        total *= (p ** (e + 1) - 1) // (p - 1)
    return total


def tau_from_factor(factors: list[tuple[int, int]]) -> int:
    total = 1
    for _, e in factors:
        total *= e + 1
    return total


def abundancy_from_factor(factors: list[tuple[int, int]]) -> Fraction:
    out = Fraction(1, 1)
    for p, e in factors:
        out *= Fraction(p ** (e + 1) - 1, p**e * (p - 1))
    return out


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def canonical_cycle(cycle: list[int]) -> tuple[int, ...]:
    if not cycle:
        return ()
    rotations = [tuple(cycle[i:] + cycle[:i]) for i in range(len(cycle))]
    return min(rotations)


def aliquot_cycles(s: list[int], limit: int, max_steps: int = 256) -> set[tuple[int, ...]]:
    cycles: set[tuple[int, ...]] = set()
    for start in range(1, limit + 1):
        seen: dict[int, int] = {}
        path: list[int] = []
        n = start
        for _ in range(max_steps):
            if n <= 0 or n > limit:
                break
            if n in seen:
                cycles.add(canonical_cycle(path[seen[n] :]))
                break
            seen[n] = len(path)
            path.append(n)
            n = s[n]
    return cycles


@dataclass(frozen=True)
class Lens:
    name: str
    fixed_predicate: int
    divisor_side_channel: int
    dynamics: int
    transfer: int
    scalar_risk: int


LENSES = [
    Lens("divisor_lattice_fixed_point", 5, 5, 4, 5, 1),
    Lens("sigma_product_abundancy", 5, 4, 3, 4, 1),
    Lens("aliquot_function_graph", 5, 3, 5, 4, 1),
    Lens("euclid_euler_mersenne_section", 4, 5, 2, 4, 1),
    Lens("local_prime_power_ledger", 4, 5, 2, 4, 2),
    Lens("odd_perfect_side_channel", 4, 4, 3, 5, 2),
    Lens("near_fixed_defect_scout", 3, 2, 4, 3, 2),
    Lens("raw_sequence_list", 2, 1, 1, 1, 4),
]


GAUGES = {
    "fixed": (4, 3, 1, 2, -2),
    "dynamics": (3, 1, 4, 2, -2),
    "transfer": (2, 3, 2, 4, -2),
}


def gauge_score(lens: Lens, gauge: str) -> int:
    w = GAUGES[gauge]
    values = (
        lens.fixed_predicate,
        lens.divisor_side_channel,
        lens.dynamics,
        lens.transfer,
        lens.scalar_risk,
    )
    return sum(a * b for a, b in zip(w, values))


def beats_in_gauge(a: Lens, b: Lens, gauge: str) -> bool:
    sa = gauge_score(a, gauge)
    sb = gauge_score(b, gauge)
    if sa != sb:
        return sa > sb
    return TIE_PATH.index(a.name) < TIE_PATH.index(b.name)


def beats(a: Lens, b: Lens) -> bool:
    votes = sum(1 for gauge in GAUGES if beats_in_gauge(a, b, gauge))
    if votes != len(GAUGES) / 2:
        return votes > len(GAUGES) / 2
    return TIE_PATH.index(a.name) < TIE_PATH.index(b.name)


def score_hist(lenses: list[Lens]) -> Counter[int]:
    hist: Counter[int] = Counter()
    for a in lenses:
        out = sum(1 for b in lenses if a != b and beats(a, b))
        hist[out] += 1
    return hist


def directed_triangles(lenses: list[Lens]) -> int:
    total = 0
    for i in range(len(lenses)):
        for j in range(i + 1, len(lenses)):
            for k in range(j + 1, len(lenses)):
                wins = []
                for x, y in ((i, j), (i, k), (j, k)):
                    wins.append((x, y) if beats(lenses[x], lenses[y]) else (y, x))
                out = Counter(x for x, _ in wins)
                if sorted(out.values()) == [1, 1, 1]:
                    total += 1
    return total


def strongly_connected_components(lenses: list[Lens]) -> list[list[str]]:
    n = len(lenses)
    graph = {i: [j for j in range(n) if i != j and beats(lenses[i], lenses[j])] for i in range(n)}
    rev = {i: [] for i in range(n)}
    for i, outs in graph.items():
        for j in outs:
            rev[j].append(i)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    comps: list[list[str]] = []
    seen.clear()

    def rdfs(v: int, comp: list[int]) -> None:
        seen.add(v)
        comp.append(v)
        for w in rev[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[int] = []
            rdfs(v, comp)
            comps.append([lenses[i].name for i in comp])
    return comps


def hamiltonian_paths(lenses: list[Lens]) -> int:
    n = len(lenses)
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp.get((mask, last), 0)
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if beats(lenses[last], lenses[nxt]):
                    key = (mask | (1 << nxt), nxt)
                    dp[key] = dp.get(key, 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, i), 0) for i in range(n))


def edge_flips(lenses: list[Lens], g1: str, g2: str) -> int:
    flips = 0
    for i, a in enumerate(lenses):
        for b in lenses[i + 1 :]:
            if beats_in_gauge(a, b, g1) != beats_in_gauge(a, b, g2):
                flips += 1
    return flips


def print_perfect_packet(sigma: list[int], s: list[int]) -> None:
    perfect = [n for n in range(2, LIMIT + 1) if s[n] == n]
    print("Fixed points")
    print("------------")
    print(f"aliquot fixed points n<= {LIMIT}: {perfect}")
    for n in perfect:
        factors = factor(n)
        abundance = abundancy_from_factor(factors)
        print(
            f"  n={n:<5} factors={fmt_factor(factors):<16} "
            f"tau={tau_from_factor(factors):<2} sigma={sigma[n]:<6} "
            f"s(n)={s[n]:<5} A={abundance}"
        )
    print()


def print_euclid_section() -> None:
    print("Euclid-Euler even section")
    print("-------------------------")
    rows = []
    for p in range(2, 32):
        mersenne = 2**p - 1
        if is_prime(mersenne):
            n = 2 ** (p - 1) * mersenne
            rows.append((p, mersenne, n, n <= LIMIT))
    for p, mersenne, n, in_window in rows:
        mark = "in window" if in_window else "beyond window"
        print(f"  p={p:<2} M=2^p-1={mersenne:<8} n=2^(p-1)M={n:<12} {mark}")
    print("  identity: sigma(2^(p-1))*sigma(M) = (2^p-1)*2^p = 2n")
    print()


def print_cycle_packet(cycles: set[tuple[int, ...]]) -> None:
    print("Aliquot cycles in the finite window")
    print("-----------------------------------")
    by_len: dict[int, list[tuple[int, ...]]] = defaultdict(list)
    for cyc in cycles:
        by_len[len(cyc)].append(cyc)
    for length in sorted(by_len):
        sample = sorted(by_len[length])[:8]
        print(f"  length {length}: count={len(by_len[length])}, sample={sample}")
    print("  fixed points are length-1 cycles; amicable pairs are length-2 cycles.")
    print()


def print_defect_packet(sigma: list[int], s: list[int]) -> None:
    print("Defect and near-fixed rows")
    print("--------------------------")
    counts = Counter()
    for n in range(1, LIMIT + 1):
        defect = sigma[n] - 2 * n
        counts["perfect" if defect == 0 else "abundant" if defect > 0 else "deficient"] += 1
    print(f"status counts n<= {LIMIT}: {dict(counts)}")
    almost = [n for n in range(2, LIMIT + 1) if sigma[n] - 2 * n == -1]
    quasi = [n for n in range(2, LIMIT + 1) if sigma[n] - 2 * n == 1]
    print(f"  defect -1 rows (almost-perfect in this window): {almost[:18]} ...")
    print(f"  defect +1 rows (quasi-perfect candidates in this window): {quasi}")

    odd_near = []
    for n in range(3, LIMIT + 1, 2):
        defect = sigma[n] - 2 * n
        if defect == 0:
            continue
        odd_near.append((Fraction(abs(defect), n), n, defect, sigma[n], s[n], fmt_factor(factor(n))))
    odd_near.sort()
    print("  closest odd rows to A=2 by |sigma-2n|/n:")
    for ratio, n, defect, sig, aliquot, factors in odd_near[:10]:
        print(f"    n={n:<5} defect={defect:<5} sigma={sig:<6} s(n)={aliquot:<6} ratio={ratio} factors={factors}")
    print()


def print_tournament() -> None:
    print("Tournament Analysis")
    print("-------------------")
    print("Pairwise observable: carrier lens quality for proving or exploiting s(n)=n.")
    print("Switch/gauge: fixed-point proof, aliquot dynamics, and transfer-to-repo side channels.")
    print(f"Tie Hamiltonian path: {' > '.join(TIE_PATH)}")
    ranking = sorted(LENSES, key=lambda lens: sum(beats(lens, other) for other in LENSES if lens != other), reverse=True)
    print(f"majority ranking: {' > '.join(l.name for l in ranking)}")
    print(f"score_hist={dict(sorted(score_hist(LENSES).items()))}")
    print(f"directed_3_cycles={directed_triangles(LENSES)}")
    print(f"sccs={strongly_connected_components(LENSES)}")
    print(f"hamiltonian_paths={hamiltonian_paths(LENSES)}")
    pairs = list(GAUGES)
    for i, g1 in enumerate(pairs):
        for g2 in pairs[i + 1 :]:
            print(f"edge_flips[{g1},{g2}]={edge_flips(LENSES, g1, g2)}")
    print()


def main() -> None:
    sigma = sigma_sieve(LIMIT)
    s = [0] * (LIMIT + 1)
    for n in range(1, LIMIT + 1):
        s[n] = sigma[n] - n

    print("S645 Perfect Numbers as Aliquot Fixed Points")
    print("============================================")
    print()
    print("Carrier")
    print("-------")
    print("divisor lattice D(n) -> sigma(n) -> aliquot map s(n)=sigma(n)-n")
    print("perfect means s(n)=n, equivalently sigma(n)=2n or abundancy A(n)=2")
    print("the side channel is the prime-power product for sigma(n)/n")
    print()
    print_perfect_packet(sigma, s)
    print_euclid_section()
    print_cycle_packet(aliquot_cycles(s, LIMIT))
    print_defect_packet(sigma, s)
    print_tournament()
    print("Synthesis")
    print("---------")
    print("Perfect numbers are not merely special values of sigma; they are loops in")
    print("the aliquot dynamical graph.  The scalar fixed equation sigma(n)=2n is")
    print("usable only with divisor side channels: prime powers, parity, Mersenne")
    print("primality for the even section, and the still-missing odd-perfect packet.")


if __name__ == "__main__":
    main()
