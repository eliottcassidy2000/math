#!/usr/bin/env python3
"""Scout the Fibonacci/Farey/additive-basis carrier bridge.

This is a synthesis computation, not a proof.  It connects four existing repo
threads:

* Fermat polygonal numbers as bounded-arity additive bases.
* Binary/ternary Goldbach as prime representation hypergraphs.
* Zeckendorf/Fibonacci as sparse-carry normal forms.
* Farey fraction addresses with separate sum/product/power payloads.

Tournament Analysis is over proof carriers rather than runners.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb, gcd, log10


N = 500


def fibs(count: int) -> list[int]:
    out = [0, 1]
    while len(out) <= count:
        out.append(out[-1] + out[-2])
    return out


def sparse_fibonacci_rows(rows: int = 12) -> list[tuple[int, list[int], int, int]]:
    fs = fibs(rows + 3)
    table = []
    for n in range(rows + 1):
        row = [comb(n - k, k) for k in range(n // 2 + 1)]
        table.append((n, row, sum(row), fs[n + 1]))
    return table


def pascal_slope_row(d: int, n: int) -> list[int]:
    """Slope-d Pascal carry row, normalized so d=2 gives F_{n+1}."""
    return [comb(n - (d - 1) * k, k) for k in range(n // d + 1)]


def pascal_slope_family(max_d: int = 6, rows: int = 12) -> list[dict[str, object]]:
    names = {
        1: "full Pascal / binary subsets / 2^n",
        2: "Fibonacci / Zeckendorf no-adjacent carry",
        3: "Narayana cows / gap-3 carry",
        4: "gap-4 carry gas",
        5: "plastic-number slope",
        6: "post-Pisot warning slope",
    }
    out: list[dict[str, object]] = []
    for d in range(1, max_d + 1):
        seq = [sum(pascal_slope_row(d, n)) for n in range(rows)]
        recurrence_ok = all(seq[n] == seq[n - 1] + seq[n - d] for n in range(d, rows))
        out.append(
            {
                "d": d,
                "name": names.get(d, f"gap-{d} carry"),
                "seq": seq,
                "sample_rows": [pascal_slope_row(d, n) for n in range(min(7, rows))],
                "recurrence_ok": recurrence_ok,
            }
        )
    return out


def farey(order: int) -> list[tuple[int, int]]:
    fracs: list[tuple[int, int]] = []
    for q in range(1, order + 1):
        for p in range(q + 1):
            if gcd(p, q) == 1:
                fracs.append((p, q))
    return sorted(fracs, key=lambda x: x[0] / x[1])


def farey_adjacent_count(order: int) -> tuple[int, int]:
    fracs = farey(order)
    total = 0
    det_one = 0
    for (a, b), (c, d) in zip(fracs, fracs[1:]):
        total += 1
        if b * c - a * d == 1:
            det_one += 1
    return total, det_one


def golden_spine(rows: int = 9) -> list[dict[str, object]]:
    fs = fibs(rows + 5)
    out: list[dict[str, object]] = []
    for i in range(2, rows + 2):
        p, q = fs[i], fs[i + 1]
        sum_payload = p + q
        product_payload = p * q
        denpow_log = p * log10(q) if q > 0 else 0.0
        numpow_log = q * log10(p) if p > 0 else 0.0
        if abs(denpow_log - numpow_log) < 1e-12:
            power_winner = "tie"
        else:
            power_winner = "q^p" if denpow_log > numpow_log else "p^q"
        out.append(
            {
                "i": i,
                "fraction": f"{p}/{q}",
                "sum": sum_payload,
                "next_fib": fs[i + 2],
                "product": product_payload,
                "log10_q^p": denpow_log,
                "log10_p^q": numpow_log,
                "power_winner": power_winner,
            }
        )
    return out


def unit_excess_lrc_chain(limit: int = 8, n: int = 14) -> list[dict[str, int]]:
    out: list[dict[str, int]] = []
    for p in range(1, limit + 1):
        q = n * p - 1
        out.append(
            {
                "p": p,
                "q": q,
                "excess": n * p - q,
                "sum": p + q,
                "product": p * q,
            }
        )
    return out


def sieve(limit: int) -> list[int]:
    is_prime = [True] * (limit + 1)
    if limit >= 0:
        is_prime[0] = False
    if limit >= 1:
        is_prime[1] = False
    p = 2
    while p * p <= limit:
        if is_prime[p]:
            for q in range(p * p, limit + 1, p):
                is_prime[q] = False
        p += 1
    return [i for i, ok in enumerate(is_prime) if ok]


def goldbach_counts(limit: int) -> tuple[list[int], list[int]]:
    primes = sieve(limit)
    prime_set = set(primes)
    missing_binary = []
    missing_ternary = []
    for n in range(4, limit + 1, 2):
        if not any((n - p) in prime_set for p in primes if p <= n):
            missing_binary.append(n)
    for n in range(7, limit + 1, 2):
        found = False
        for p in primes:
            if p > n:
                break
            for q in primes:
                if p + q > n:
                    break
                if (n - p - q) in prime_set:
                    found = True
                    break
            if found:
                break
        if not found:
            missing_ternary.append(n)
    return missing_binary, missing_ternary


def polygonal(sides: int, k: int) -> int:
    return ((sides - 2) * k * k - (sides - 4) * k) // 2


def polygonal_terms(sides: int, limit: int) -> list[int]:
    terms = []
    k = 0
    while True:
        v = polygonal(sides, k)
        if v > limit:
            break
        terms.append(v)
        k += 1
    return terms


def min_terms(vals: list[int], limit: int) -> int:
    inf = 10**9
    dp = [inf] * (limit + 1)
    dp[0] = 0
    atoms = [v for v in vals if v > 0]
    for n in range(1, limit + 1):
        dp[n] = 1 + min(dp[n - v] for v in atoms if v <= n)
    return max(dp[1:])


def zeckendorf(n: int, atoms: list[int]) -> list[int]:
    rem = n
    out = []
    for f in reversed(atoms):
        if f <= rem:
            out.append(f)
            rem -= f
    assert rem == 0
    return list(reversed(out))


@dataclass(frozen=True)
class Carrier:
    name: str
    retained: int
    proof_power: int
    scalar_safety: int
    lrc_transfer: int

    @property
    def score(self) -> tuple[int, int, int, int]:
        return (self.retained, self.proof_power, self.scalar_safety, self.lrc_transfer)


@dataclass(frozen=True)
class Economy:
    name: str
    carrier: str
    preserves: str
    forgets: str
    lrc_pull: str


def tournament_fingerprint(carriers: list[Carrier]) -> dict[str, object]:
    n = len(carriers)
    outdeg = [0] * n
    c3 = 0
    edges: dict[tuple[int, int], int] = {}
    for i in range(n):
        for j in range(i + 1, n):
            if carriers[i].score >= carriers[j].score:
                winner = i
            else:
                winner = j
            loser = j if winner == i else i
            edges[(winner, loser)] = 1
            outdeg[winner] += 1
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                wins = {
                    (a, b)
                    for a, b in edges
                    if a in (i, j, k) and b in (i, j, k)
                }
                cyclic = (
                    ((i, j) in wins and (j, k) in wins and (k, i) in wins)
                    or ((i, k) in wins and (k, j) in wins and (j, i) in wins)
                )
                if cyclic:
                    c3 += 1
    order = sorted(range(n), key=lambda idx: carriers[idx].score, reverse=True)
    return {
        "score_hist": sorted((d, outdeg.count(d)) for d in set(outdeg)),
        "directed_3cycles": c3,
        "hamiltonian_path": " > ".join(carriers[i].name for i in order),
    }


def main() -> None:
    print("=== Farey/Fibonacci/additive-basis synthesis S169 ===")
    print()

    print("Pascal-slope carry family: a_d(n)=sum_k binom(n-(d-1)k,k)")
    for record in pascal_slope_family():
        seq = ",".join(map(str, record["seq"]))
        print(
            f"  d={record['d']}: {record['name']:<42s}"
            f" recurrence a(n)=a(n-1)+a(n-d): {record['recurrence_ok']}"
        )
        print(f"       seq n=0..11: {seq}")
        if record["d"] in (1, 2, 3):
            rows = ["+".join(map(str, row)) for row in record["sample_rows"]]
            print(f"       first rows: {'; '.join(rows)}")
    print("  readout: the named Fibonacci row is d=2; d=1 is full Pascal, d=3 is Narayana.")
    print("           The parameter d is the minimum gap forced between carry sites.")
    print()

    print("Fibonacci sparse-carry rows: F_{n+1} = sum_k binom(n-k,k)")
    for n, row, total, fib in sparse_fibonacci_rows(11):
        print(f"  n={n:2d}: {'+'.join(map(str, row)):<18s} = {total:3d}  F={fib:3d}")
    print()

    print("Golden Farey/Stern-Brocot spine: p/q = F_i/F_{i+1}")
    print("  i  p/q      p+q nextF  p*q    power winner  log10(q^p) log10(p^q)")
    for row in golden_spine(10):
        print(
            f" {row['i']:2d}  {row['fraction']:<7s}"
            f" {row['sum']:4d} {row['next_fib']:5d}"
            f" {row['product']:6d} {row['power_winner']:^12s}"
            f" {row['log10_q^p']:11.3f} {row['log10_p^q']:11.3f}"
        )
    print("  readout: additive payload is literally the next Fibonacci on the all-1 continued-fraction spine.")
    print("           product is the K_{p,q} area; powers are ordered magnitude-stress channels.")
    print()

    print("Farey adjacency sanity")
    for order in [5, 8, 13, 21]:
        total, det_one = farey_adjacent_count(order)
        print(f"  F_{order}: adjacent pairs={total:3d}, determinant-one pairs={det_one:3d}")
    print("  readout: adjacent Farey edges are unimodular; mediants add numerator and denominator vectors.")
    print()

    print("Unit-excess LRC/Farey scheduler: q=14p-1")
    print("  p   q  e=14p-q  p+q  p*q")
    for row in unit_excess_lrc_chain():
        print(
            f" {row['p']:2d} {row['q']:3d} {row['excess']:8d}"
            f" {row['sum']:4d} {row['product']:5d}"
        )
    print("  readout: q is the binding scale, p+q is a linear recursion clock,")
    print("           and p*q is a quadratic K_{p,q} incidence clock.")
    print()

    missing_binary, missing_ternary = goldbach_counts(N)
    print("Additive-basis finite audit")
    print(f"  Goldbach binary missing evens <= {N}: {missing_binary[:8]} count={len(missing_binary)}")
    print(f"  Goldbach ternary missing odds <= {N}: {missing_ternary[:8]} count={len(missing_ternary)}")
    for sides in range(3, 9):
        mt = min_terms(polygonal_terms(sides, N), N)
        print(f"  {sides}-gonal max summands <= {N}: {mt}")
    fs = [1, 2]
    while fs[-1] + fs[-2] <= N:
        fs.append(fs[-1] + fs[-2])
    max_z = max(len(zeckendorf(n, fs)) for n in range(1, N + 1))
    print(f"  Zeckendorf atoms <= {N}: {len(fs)}, max digits={max_z}, atoms={fs}")
    print()

    economies = [
        Economy(
            "Goldbach",
            "prime pair graph",
            "local residue and singular-series side data",
            "which prime-pair branch caused the coverage",
            "smoothing only after residue clocks are retained",
        ),
        Economy(
            "Ternary Goldbach",
            "prime triple hypergraph",
            "extra-summand smoothing and major/minor arc labels",
            "binary branch identity",
            "model for admissible smoothing dispatchers",
        ),
        Economy(
            "Fermat polygonal",
            "s-gonal atom multiset",
            "bounded summand budget and local obstruction invoice",
            "representation multiplicity beyond the arity cap",
            "bounded residual debt rather than raw counts",
        ),
        Economy(
            "Zeckendorf",
            "Fibonacci carry automaton",
            "canonical no-adjacent normal form",
            "redundant representations killed by confluence",
            "normal-form quotient if the carry rule survives",
        ),
        Economy(
            "Farey address",
            "fraction vector (p,q)",
            "q, p+q, p*q, and ordered power clocks separately",
            "order and incidence data after scalarization",
            "packet classifier field for sequence shadows",
        ),
    ]
    print("Representation-economy ledger")
    for e in economies:
        print(f"  {e.name}:")
        print(f"    carrier   = {e.carrier}")
        print(f"    preserves = {e.preserves}")
        print(f"    forgets   = {e.forgets}")
        print(f"    LRC pull  = {e.lrc_pull}")
    print()

    carriers = [
        Carrier("zeckendorf_sparse_normal_form", 5, 5, 5, 4),
        Carrier("farey_address_vector", 5, 4, 4, 5),
        Carrier("fermat_polygonal_bounded_arity", 4, 4, 4, 3),
        Carrier("ternary_goldbach_smoothing", 3, 5, 3, 3),
        Carrier("binary_goldbach_pair_graph", 3, 3, 2, 3),
        Carrier("farey_product_Kpq_area", 4, 3, 3, 5),
        Carrier("farey_power_stress_channel", 2, 2, 5, 2),
        Carrier("raw_scalar_rep_count", 1, 1, 1, 1),
    ]
    fp = tournament_fingerprint(carriers)
    print("Tournament Analysis over proof carriers")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  hamiltonian_path={fp['hamiltonian_path']}")
    print()

    print("Assumption challenge / alternate tournament vertices")
    print("  considered: integers, summands, primes, polygonal atoms, Fibonacci carry sites,")
    print("              Farey fractions, K_{p,q} incidences, power-order tests, proof obligations.")
    print("  chosen: proof carriers, because LRC needs to know which quotient preserves the predicate.")
    print("  destroyed by raw scalarization: residue clocks, arity budgets, carry adjacency,")
    print("                                 Farey order, Kpq incidence, and cocycle debt.")
    print()

    print("Synthesis")
    print("  Goldbach: many prime branches plus local residue/singular-series side data.")
    print("  Ternary Goldbach: one extra summand turns residue noise into enough smoothing.")
    print("  Fermat polygonal: bounded arity pays a finite local-obstruction invoice.")
    print("  Zeckendorf: no-adjacent sparse carries give a unique normal form.")
    print("  Farey: numerator/denominator vector addition is mediant recursion; product is")
    print("         incidence area; powers are noncommutative stress tests for lost order.")


if __name__ == "__main__":
    main()
