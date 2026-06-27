#!/usr/bin/env python3
"""Perfect-number packet merge for the LRC14 automaton stack.

This is a synthesis computation, not a proof.  It imports the older perfect
number / aliquot fixed-point carrier into the current LRC14 automatic-gap
packet stack.  The useful comparison is not a scalar product value: it is the
exact Euclid-Euler unit-excess control at n=2 versus the LRC14 n=14 shadow,
with divisor-lattice and prime-factor labels retained.

Tournament Analysis is over proof carriers and packet side channels, not
runners, arcs, or raw sequence entries.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
from math import isqrt


MAX_MERSENNE_EXPONENT = 17
MAX_POWER_K = 12


def digits(n: int, base: int) -> tuple[int, ...]:
    if n == 0:
        return (0,)
    out: list[int] = []
    while n:
        out.append(n % base)
        n //= base
    return tuple(reversed(out))


def digit_string(n: int, base: int) -> str:
    alphabet = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    return "".join(alphabet[d] for d in digits(n, base))


def is_fibbinary(n: int) -> bool:
    return "11" not in digit_string(n, 2)


def is_moser_de_bruijn(n: int) -> bool:
    return all(d in (0, 1) for d in digits(n, 4))


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    step = 2
    d = 5
    while d * d <= n:
        if n % d == 0:
            return False
        d += step
        step = 6 - step
    return True


def factorization(n: int) -> tuple[tuple[int, int], ...]:
    out: list[tuple[int, int]] = []
    d = 2
    while d <= isqrt(n):
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            out.append((d, e))
        d = 3 if d == 2 else d + 2
    if n > 1:
        out.append((n, 1))
    return tuple(out)


def factor_string(n: int) -> str:
    parts = []
    for p, e in factorization(n):
        parts.append(str(p) if e == 1 else f"{p}^{e}")
    return " * ".join(parts) if parts else "1"


def sigma_from_factorization(factors: tuple[tuple[int, int], ...]) -> int:
    total = 1
    for p, e in factors:
        total *= (p ** (e + 1) - 1) // (p - 1)
    return total


def sigma(n: int) -> int:
    return sigma_from_factorization(factorization(n))


@dataclass(frozen=True)
class PerfectControlRow:
    r: int
    a: int
    q2: int
    n: int
    abundancy: Fraction
    defect: Fraction
    a_fibbinary: bool
    a_moser: bool


@dataclass(frozen=True)
class LrcShadowRow:
    k: int
    a: int
    q14: int
    n: int
    q14_prime: bool
    q14_factorization: str
    abundancy: Fraction
    defect: Fraction
    prime_formula_defect: Fraction | None
    a_fibbinary: bool
    a_moser: bool
    q14_fibbinary: bool
    q14_moser: bool


def perfect_control_rows() -> list[PerfectControlRow]:
    rows: list[PerfectControlRow] = []
    for r in range(2, MAX_MERSENNE_EXPONENT + 1):
        q2 = 2**r - 1
        if not is_prime(q2):
            continue
        a = 2 ** (r - 1)
        n = a * q2
        abundancy = Fraction(sigma(n), n)
        rows.append(
            PerfectControlRow(
                r=r,
                a=a,
                q2=q2,
                n=n,
                abundancy=abundancy,
                defect=Fraction(2, 1) - abundancy,
                a_fibbinary=is_fibbinary(a),
                a_moser=is_moser_de_bruijn(a),
            )
        )
    return rows


def lrc_shadow_rows() -> list[LrcShadowRow]:
    rows: list[LrcShadowRow] = []
    for k in range(MAX_POWER_K + 1):
        a = 2**k
        q14 = 14 * a - 1
        n = a * q14
        abundancy = Fraction(sigma(n), n)
        q14_prime = is_prime(q14)
        rows.append(
            LrcShadowRow(
                k=k,
                a=a,
                q14=q14,
                n=n,
                q14_prime=q14_prime,
                q14_factorization=factor_string(q14),
                abundancy=abundancy,
                defect=Fraction(2, 1) - abundancy,
                prime_formula_defect=Fraction(12, q14) if q14_prime else None,
                a_fibbinary=is_fibbinary(a),
                a_moser=is_moser_de_bruijn(a),
                q14_fibbinary=is_fibbinary(q14),
                q14_moser=is_moser_de_bruijn(q14),
            )
        )
    return rows


@dataclass(frozen=True)
class Carrier:
    name: str
    lrc_predicate: int
    exact_scale: int
    divisor_lattice: int
    graph_minor: int
    automaton_state: int
    power_exception: int
    anti_scalar: int
    proof_maturity: int
    destroys: str

    def vector(self) -> dict[str, int]:
        return {
            "lrc_predicate": self.lrc_predicate,
            "exact_scale": self.exact_scale,
            "divisor_lattice": self.divisor_lattice,
            "graph_minor": self.graph_minor,
            "automaton_state": self.automaton_state,
            "power_exception": self.power_exception,
            "anti_scalar": self.anti_scalar,
            "proof_maturity": self.proof_maturity,
        }


FIELDS = (
    "lrc_predicate",
    "exact_scale",
    "divisor_lattice",
    "graph_minor",
    "automaton_state",
    "power_exception",
    "anti_scalar",
    "proof_maturity",
)

TIE_PATH = (
    "labelled_lrc_packet_sheaf",
    "divisor_lattice_abundancy_packet",
    "exact_farey_unit_excess",
    "lrc14_n14_deficient_shadow",
    "perfect_n2_fixed_control",
    "k33_product_incidence",
    "automatic_power_of_two_state",
    "fermat_catalan_power_budget",
    "raw_product_scalar",
)


def carrier_bank() -> list[Carrier]:
    return [
        Carrier(
            "labelled_lrc_packet_sheaf",
            lrc_predicate=4,
            exact_scale=4,
            divisor_lattice=4,
            graph_minor=4,
            automaton_state=4,
            power_exception=4,
            anti_scalar=4,
            proof_maturity=3,
            destroys="nothing if exact M, divisor factorization, automaton state, and route labels are retained",
        ),
        Carrier(
            "divisor_lattice_abundancy_packet",
            lrc_predicate=2,
            exact_scale=3,
            divisor_lattice=4,
            graph_minor=1,
            automaton_state=1,
            power_exception=3,
            anti_scalar=4,
            proof_maturity=3,
            destroys="endpoint owners and K33 route unless reattached",
        ),
        Carrier(
            "exact_farey_unit_excess",
            lrc_predicate=4,
            exact_scale=4,
            divisor_lattice=2,
            graph_minor=3,
            automaton_state=1,
            power_exception=2,
            anti_scalar=3,
            proof_maturity=3,
            destroys="divisor-chain differences between prime and composite q",
        ),
        Carrier(
            "lrc14_n14_deficient_shadow",
            lrc_predicate=4,
            exact_scale=3,
            divisor_lattice=3,
            graph_minor=2,
            automaton_state=1,
            power_exception=3,
            anti_scalar=3,
            proof_maturity=2,
            destroys="perfect fixed-point equality but keeps the n=14 defect label",
        ),
        Carrier(
            "perfect_n2_fixed_control",
            lrc_predicate=1,
            exact_scale=4,
            divisor_lattice=4,
            graph_minor=1,
            automaton_state=1,
            power_exception=4,
            anti_scalar=3,
            proof_maturity=4,
            destroys="the LRC14 threshold and endpoint geometry",
        ),
        Carrier(
            "k33_product_incidence",
            lrc_predicate=3,
            exact_scale=2,
            divisor_lattice=1,
            graph_minor=4,
            automaton_state=1,
            power_exception=1,
            anti_scalar=3,
            proof_maturity=3,
            destroys="aliquot/divisor data and automatic carries",
        ),
        Carrier(
            "automatic_power_of_two_state",
            lrc_predicate=2,
            exact_scale=1,
            divisor_lattice=1,
            graph_minor=1,
            automaton_state=4,
            power_exception=3,
            anti_scalar=3,
            proof_maturity=2,
            destroys="prime-q status and abundancy defect",
        ),
        Carrier(
            "fermat_catalan_power_budget",
            lrc_predicate=1,
            exact_scale=1,
            divisor_lattice=2,
            graph_minor=1,
            automaton_state=1,
            power_exception=4,
            anti_scalar=2,
            proof_maturity=2,
            destroys="exact Farey endpoint data and Kpq incidence",
        ),
        Carrier(
            "raw_product_scalar",
            lrc_predicate=0,
            exact_scale=0,
            divisor_lattice=0,
            graph_minor=1,
            automaton_state=0,
            power_exception=0,
            anti_scalar=0,
            proof_maturity=0,
            destroys="all packet coordinates except a magnitude shadow",
        ),
    ]


def compare(a: Carrier, b: Carrier) -> str:
    va = a.vector()
    vb = b.vector()
    wins_a = sum(1 for field in FIELDS if va[field] > vb[field])
    wins_b = sum(1 for field in FIELDS if vb[field] > va[field])
    if wins_a > wins_b:
        return a.name
    if wins_b > wins_a:
        return b.name
    order = {name: i for i, name in enumerate(TIE_PATH)}
    return a.name if order[a.name] < order[b.name] else b.name


def edge_map(carriers: list[Carrier]) -> dict[tuple[str, str], str]:
    edges: dict[tuple[str, str], str] = {}
    for a, b in combinations(carriers, 2):
        winner = compare(a, b)
        loser = b.name if winner == a.name else a.name
        edges[(winner, loser)] = winner
    return edges


def beats(edges: dict[tuple[str, str], str], a: str, b: str) -> bool:
    return (a, b) in edges


def score_histogram(carriers: list[Carrier], edges: dict[tuple[str, str], str]) -> dict[int, int]:
    scores = {carrier.name: 0 for carrier in carriers}
    for winner, _loser in edges:
        scores[winner] += 1
    return dict(sorted(Counter(scores.values()).items()))


def directed_3cycles(names: list[str], edges: dict[tuple[str, str], str]) -> int:
    count = 0
    for a, b, c in combinations(names, 3):
        if (
            beats(edges, a, b)
            and beats(edges, b, c)
            and beats(edges, c, a)
            or beats(edges, a, c)
            and beats(edges, c, b)
            and beats(edges, b, a)
        ):
            count += 1
    return count


def scc_sizes(names: list[str], edges: dict[tuple[str, str], str]) -> list[int]:
    graph: dict[str, list[str]] = defaultdict(list)
    reverse: dict[str, list[str]] = defaultdict(list)
    for a, b in edges:
        graph[a].append(b)
        reverse[b].append(a)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for name in names:
        if name not in seen:
            dfs(name)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: str) -> int:
        seen.add(v)
        size = 1
        for w in reverse[v]:
            if w not in seen:
                size += rdfs(w)
        return size

    for name in reversed(order):
        if name not in seen:
            sizes.append(rdfs(name))
    return sorted(sizes, reverse=True)


def hamiltonian_paths(names: list[str], edges: dict[tuple[str, str], str]) -> tuple[int, tuple[str, ...] | None]:
    count = 0
    first: tuple[str, ...] | None = None
    for perm in permutations(names):
        if all(beats(edges, perm[i], perm[i + 1]) for i in range(len(perm) - 1)):
            count += 1
            if first is None:
                first = perm
    return count, first


def format_fraction(frac: Fraction) -> str:
    return str(frac.numerator) if frac.denominator == 1 else f"{frac.numerator}/{frac.denominator}"


def print_perfect_controls(rows: list[PerfectControlRow]) -> None:
    print("Even perfect controls: n=2 unit-excess fixed points")
    print("r a q2=2a-1 N sigma(N)/N defect a_fibbinary a_moser")
    for row in rows:
        print(
            row.r,
            row.a,
            row.q2,
            row.n,
            format_fraction(row.abundancy),
            format_fraction(row.defect),
            int(row.a_fibbinary),
            int(row.a_moser),
        )
    print()


def print_lrc_shadows(rows: list[LrcShadowRow]) -> None:
    print("LRC14 power-of-two shadows: q=14a-1")
    print(
        "k a q14 prime factorization sigma(aq)/aq defect "
        "prime_formula_defect a_fibbinary a_moser q_fibbinary q_moser"
    )
    for row in rows:
        print(
            row.k,
            row.a,
            row.q14,
            int(row.q14_prime),
            row.q14_factorization,
            format_fraction(row.abundancy),
            format_fraction(row.defect),
            format_fraction(row.prime_formula_defect) if row.prime_formula_defect else "NA",
            int(row.a_fibbinary),
            int(row.a_moser),
            int(row.q14_fibbinary),
            int(row.q14_moser),
        )
    print()


def print_summary(control_rows: list[PerfectControlRow], shadow_rows: list[LrcShadowRow]) -> None:
    prime_shadow_rows = [row for row in shadow_rows if row.q14_prime]
    composite_shadow_rows = [row for row in shadow_rows if not row.q14_prime]
    formula_ok = all(row.defect == row.prime_formula_defect for row in prime_shadow_rows)
    composite_abundant = [row for row in composite_shadow_rows if row.defect < 0]
    print("Merge readout")
    print(f"perfect_control_rows={len(control_rows)}")
    print(f"lrc14_shadow_rows={len(shadow_rows)}")
    print(f"prime_q14_shadow_rows={len(prime_shadow_rows)}")
    print(f"composite_q14_shadow_rows={len(composite_shadow_rows)}")
    print(f"prime_q14_defect_formula_ok={int(formula_ok)}")
    print(f"composite_q14_abundant_rows={[(row.k, row.a, row.q14, format_fraction(row.defect)) for row in composite_abundant]}")
    print(
        "packet_fields=unit_excess_apex, perfect_control_status, abundancy_defect, "
        "divisor_lattice_factorization, prime_q_flag, product_incidence_rank, "
        "automaton_transition_state"
    )
    print(
        "guardrail=the n=2 perfect chain is exact, while prime q=14a-1 shadows "
        "are deficient by 12/(14a-1); composite q rows prove the factorization "
        "side channel is load-bearing."
    )
    print()


def print_tournament() -> None:
    carriers = carrier_bank()
    names = [carrier.name for carrier in carriers]
    edges = edge_map(carriers)
    hp_count, hp_first = hamiltonian_paths(names, edges)
    print("Tournament Analysis")
    print("vertices=proof carriers and side channels, not runners/arcs/raw sequence entries")
    print(f"pairwise_observable={','.join(FIELDS)}")
    print(f"tie_hamiltonian_path={'>'.join(TIE_PATH)}")
    print(f"score_hist={score_histogram(carriers, edges)}")
    print(f"directed_3cycles={directed_3cycles(names, edges)}")
    print(f"scc_sizes={scc_sizes(names, edges)}")
    print(f"hamiltonian_path_count={hp_count}")
    print(f"first_hamiltonian_path={'>'.join(hp_first) if hp_first else 'NONE'}")
    print()
    print("Assumption challenge")
    print(
        "Alternate vertices considered: runners, gaps, fixed circle sections, "
        "section boundaries, wall-crossing events, residues, cover arcs, "
        "Fourier modes, matroid circuits, divisor atoms, and proof obligations."
    )
    print(
        "Chosen quotient preserves the LRC predicate only through labelled packet "
        "fields; it destroys raw runner identity and refuses to collapse divisor "
        "or automaton information to a product scalar."
    )


def main() -> None:
    control_rows = perfect_control_rows()
    shadow_rows = lrc_shadow_rows()
    print("LRC14 perfect-number packet merge (codex-2026-06-26-S174)")
    print("source_threads=HYP-2221,HYP-2941,HYP-2945,HYP-3008,HYP-3009,HYP-3011,HYP-3012")
    print()
    print_perfect_controls(control_rows)
    print_lrc_shadows(shadow_rows)
    print_summary(control_rows, shadow_rows)
    print_tournament()


if __name__ == "__main__":
    main()
