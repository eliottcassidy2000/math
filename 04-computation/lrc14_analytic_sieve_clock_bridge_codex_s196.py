#!/usr/bin/env python3
"""LRC14 analytic sieve-clock bridge.

This is a proof-interface computation, not a proof of LRC14.

HYP-2982/HYP-2983 imported Mobius/totient sums, large-sieve/circle-method
language, smoothing, exponential sums, and the Kaczynski boundary viewpoint.
HYP-3024/HYP-3027 then made the local lesson sharper: clocks such as
Erdos-Turan or exact M are useful only when their quotient is guarded.

This script asks the same question for analytic sieve quantities:

    what may mu(n)/n, mu(n)^2/phi(n), large-sieve clocks,
    smoothing choices, and Kaczynski boundary labels forget?

The finite audit uses the HYP-3026 named carrier bank plus the two residual
mixed pairs isolated by HYP-3027.  Tournament Analysis is over analytic proof
clocks and repair packets, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import log
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
SWITCHBOARD_PATH = REPO / "04-computation" / "lrc14_carrier_fusion_switchboard_codex_s189.py"


def load_switchboard():
    spec = spec_from_file_location("lrc14_carrier_fusion_switchboard_s189", SWITCHBOARD_PATH)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules["lrc14_carrier_fusion_switchboard_s189"] = module
    spec.loader.exec_module(module)
    return module


sw = load_switchboard()


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def ffloat(x: float) -> str:
    return f"{x:.6f}"


def factor_word(factors: dict[int, int]) -> str:
    if not factors:
        return "1"
    return "*".join(
        str(p) if e == 1 else f"{p}^{e}" for p, e in sorted(factors.items())
    )


def linear_mu_phi(limit: int) -> tuple[list[int], list[int], list[int], list[bool]]:
    primes: list[int] = []
    is_comp = [False] * (limit + 1)
    is_prime = [False] * (limit + 1)
    mu = [0] * (limit + 1)
    phi = [0] * (limit + 1)
    mu[1] = 1
    phi[1] = 1
    for i in range(2, limit + 1):
        if not is_comp[i]:
            primes.append(i)
            is_prime[i] = True
            mu[i] = -1
            phi[i] = i - 1
        for p in primes:
            v = i * p
            if v > limit:
                break
            is_comp[v] = True
            if i % p == 0:
                mu[v] = 0
                phi[v] = phi[i] * p
                break
            mu[v] = -mu[i]
            phi[v] = phi[i] * (p - 1)
    return primes, mu, phi, is_prime


def factorization(n: int, primes: list[int]) -> dict[int, int]:
    out: dict[int, int] = {}
    m = n
    for p in primes:
        if p * p > m:
            break
        while m % p == 0:
            out[p] = out.get(p, 0) + 1
            m //= p
    if m > 1:
        out[m] = out.get(m, 0) + 1
    return out


@dataclass(frozen=True)
class Prefix:
    mertens: list[int]
    mu_over_n: list[float]
    mu2_over_phi: list[float]
    phi_sum: list[int]
    theta: list[float]
    prime_recip: list[float]


def prefix_tables(limit: int, mu: list[int], phi: list[int], is_prime: list[bool]) -> Prefix:
    mertens = [0] * (limit + 1)
    mu_over_n = [0.0] * (limit + 1)
    mu2_over_phi = [0.0] * (limit + 1)
    phi_sum = [0] * (limit + 1)
    theta = [0.0] * (limit + 1)
    prime_recip = [0.0] * (limit + 1)
    for n in range(1, limit + 1):
        mertens[n] = mertens[n - 1] + mu[n]
        mu_over_n[n] = mu_over_n[n - 1] + mu[n] / n
        mu2_over_phi[n] = mu2_over_phi[n - 1]
        if mu[n] != 0:
            mu2_over_phi[n] += 1.0 / phi[n]
        phi_sum[n] = phi_sum[n - 1] + phi[n]
        theta[n] = theta[n - 1]
        prime_recip[n] = prime_recip[n - 1]
        if is_prime[n]:
            theta[n] += log(n)
            prime_recip[n] += 1.0 / n
    return Prefix(mertens, mu_over_n, mu2_over_phi, phi_sum, theta, prime_recip)


def augmented_rows() -> list[sw.Row]:
    rows = list(sw.named_rows())
    extra = [
        sw.Row(
            "residual_petal_drop10_13_add20_26",
            tuple(sorted((set(range(1, 14)) - {10, 13}) | {20, 26})),
            "BOUNDARY-PETAL-SPORADIC",
        ),
        sw.Row(
            "residual_cover_drop8_12_add16_24",
            tuple(sorted((set(range(1, 14)) - {8, 12}) | {16, 24})),
            "COVERING-MOMENT",
        ),
        sw.Row(
            "residual_k33_drop12_13_add26_36",
            tuple(sorted((set(range(1, 14)) - {12, 13}) | {26, 36})),
            "K33-STATE-LIFT",
        ),
        sw.Row(
            "residual_cover_single_12_to_72",
            tuple(list(range(1, 12)) + [13, 72]),
            "COVERING-MOMENT",
        ),
    ]
    seen: set[tuple[int, ...]] = set()
    out: list[sw.Row] = []
    for row in rows + extra:
        if row.speeds in seen:
            continue
        seen.add(row.speeds)
        out.append(row)
    return out


def sieve_clock(q: int, factors: dict[int, int], mu_q: int) -> str:
    if q == 14:
        return "apex_14_squarefree_unit"
    if q == 27:
        return "c27_primepower_mu2phi_blind"
    if q == 41:
        return "k33_prime_unit"
    if mu_q == 0:
        return "repeated_prime_mu2phi_blind"
    if len(factors) == 1:
        return "prime_squarefree_unit"
    return "composite_squarefree_unit"


def smoothing_label(route: str, status: str, q: int) -> str:
    if route == "BOUNDARY-AP-GW":
        return "endpoint_kaczynski_boundary_atom"
    if route == "BOUNDARY-PETAL-SPORADIC":
        return "fejer_ramanujan_c27_packet"
    if route == "K33-STATE-LIFT":
        return "k33_exponential_sum_state_lift"
    if route == "COVERING-MOMENT":
        return "boundary_moment_large_sieve_precondition"
    if route == "AUTOMATIC-OPEN-CONTROL":
        return "automatic_open_control_clock"
    if status == "open" and q < 14:
        return "direct_q_witness_major_arc"
    return "direct_q_witness_minor_arc"


def boundary_label(route: str, status: str) -> str:
    if status == "boundary":
        return "true_wide_kaczynski_ambiguous_approach"
    if route == "COVERING-MOMENT":
        return "finite_grid_boundary_moment_debt"
    if route == "K33-STATE-LIFT":
        return "nonunit_state_lift_boundary_debt"
    if route == "BOUNDARY-PETAL-SPORADIC":
        return "c27_endpoint_transfer_debt"
    return "ordinary_open_fatou_approach"


@dataclass(frozen=True)
class AnalyticRow:
    name: str
    route: str
    status: str
    exact_m: Fraction
    exact_t: Fraction
    safe_mu: Fraction
    q: int
    p: int
    factors: dict[int, int]
    mu_q: int
    phi_q: int
    mu2_phi_value: Fraction
    squarefree: bool
    clock: str
    smoothing: str
    boundary: str
    automatic_word: str
    bar_count: int
    first_chart_den: int | None
    boundary_count: int
    zero_sum_pairs: int


def build_rows() -> tuple[list[AnalyticRow], Prefix, int]:
    raw_features = [sw.row_features(row) for row in augmented_rows()]
    max_q = max(item.exact_m.denominator for item in raw_features)
    limit = max(5000, 4312, max_q)
    primes, mu, phi, is_prime = linear_mu_phi(limit)
    prefix = prefix_tables(limit, mu, phi, is_prime)

    rows: list[AnalyticRow] = []
    for item in raw_features:
        q = item.exact_m.denominator
        factors = factorization(q, primes)
        route = item.row.route
        status = item.status
        mu2_phi_value = Fraction(0) if mu[q] == 0 else Fraction(1, phi[q])
        rows.append(
            AnalyticRow(
                name=item.row.name,
                route=route,
                status=status,
                exact_m=item.exact_m,
                exact_t=item.exact_t,
                safe_mu=item.safe_mu,
                q=q,
                p=item.exact_m.numerator,
                factors=factors,
                mu_q=mu[q],
                phi_q=phi[q],
                mu2_phi_value=mu2_phi_value,
                squarefree=mu[q] != 0,
                clock=sieve_clock(q, factors, mu[q]),
                smoothing=smoothing_label(route, status, q),
                boundary=boundary_label(route, status),
                automatic_word=item.automatic_word,
                bar_count=item.bar_count,
                first_chart_den=item.first_chart_den,
                boundary_count=item.boundary_count,
                zero_sum_pairs=item.zero_sum_pairs,
            )
        )
    return rows, prefix, limit


def fiber_report(rows: list[AnalyticRow]) -> None:
    key_bank = {
        "mu2_phi_blind_flag": lambda r: "live" if r.squarefree else "blind",
        "mu2_phi_value_only": lambda r: r.mu2_phi_value,
        "sieve_clock_only": lambda r: r.clock,
        "sieve_clock_plus_exact_denominator": lambda r: (r.clock, r.q),
        "sieve_clock_plus_smoothing": lambda r: (r.clock, r.smoothing),
        "nonroute_clock_packet": lambda r: (
            r.clock,
            r.q,
            r.status,
            r.bar_count,
            r.boundary_count,
            r.zero_sum_pairs,
            r.first_chart_den,
        ),
        "declared_analytic_exit": lambda r: (r.clock, r.smoothing, r.boundary),
    }
    print("[2] Quotient stress")
    print("key fibers mixed_route mixed_status largest_fiber largest_mixed")
    for key_name, key_func in key_bank.items():
        groups: dict[object, list[AnalyticRow]] = defaultdict(list)
        for row in rows:
            groups[key_func(row)].append(row)
        mixed_route = 0
        mixed_status = 0
        largest = 0
        largest_mixed = 0
        examples: list[str] = []
        for key, group in groups.items():
            largest = max(largest, len(group))
            routes = {row.route for row in group}
            statuses = {row.status for row in group}
            if len(routes) > 1:
                mixed_route += 1
                largest_mixed = max(largest_mixed, len(group))
                if len(examples) < 3:
                    sample = ",".join(row.name for row in group[:5])
                    examples.append(f"{repr(key)}=>{sample}")
            if len(statuses) > 1:
                mixed_status += 1
        print(key_name, len(groups), mixed_route, mixed_status, largest, largest_mixed, sep=" | ")
        if examples:
            print(f"  mixed_route_examples: {' ; '.join(examples)}")
    print()


def arithmetic_checkpoint_readout(rows: list[AnalyticRow], prefix: Prefix) -> None:
    checkpoints = sorted({14, 27, 41, 84, 168, 4312} | {row.q for row in rows})
    print("[1] Analytic arithmetic checkpoints")
    print("q M(q) sum_mu/n G=sum_mu2/phi Phi(q) theta(q) sum_1/p")
    for q in checkpoints:
        if q >= len(prefix.mertens):
            continue
        print(
            q,
            prefix.mertens[q],
            ffloat(prefix.mu_over_n[q]),
            ffloat(prefix.mu2_over_phi[q]),
            prefix.phi_sum[q],
            ffloat(prefix.theta[q]),
            ffloat(prefix.prime_recip[q]),
            sep=" | ",
        )
    print()


def row_readout(rows: list[AnalyticRow], prefix: Prefix) -> None:
    print("[0] Row-level analytic packet readout")
    print(
        "name route status M@t q factors mu phi mu2/phi A(q) G(q) clock smoothing boundary"
    )
    for row in rows:
        print(
            row.name,
            row.route,
            row.status,
            f"{fmt(row.exact_m)}@{fmt(row.exact_t)}",
            row.q,
            factor_word(row.factors),
            row.mu_q,
            row.phi_q,
            fmt(row.mu2_phi_value),
            ffloat(prefix.mu_over_n[row.q]),
            ffloat(prefix.mu2_over_phi[row.q]),
            row.clock,
            row.smoothing,
            row.boundary,
            sep=" | ",
        )
    blind = [row.name for row in rows if not row.squarefree]
    print()
    print(f"squarefree_blind_rows={blind}")
    print(
        "interpretation=mu2/phi is a capacity clock with an attached blindness "
        "certificate; it is not an admissible LRC quotient by itself."
    )
    print()


@dataclass(frozen=True)
class Carrier:
    name: str
    vector: tuple[int, ...]
    destroys: str


TOURNAMENT_FIELDS = (
    "status_retention",
    "route_retention",
    "denominator_retention",
    "squarefree_blindness",
    "smoothing_retention",
    "exponential_sum_declared",
    "boundary_approach",
    "finite_checkability",
    "noncircularity",
)

TIE_PATH = (
    "labelled_repair_ladder_packet",
    "analytic_sieve_clock_bridge",
    "kaczynski_boundary_approach",
    "smoothing_explicit_formula_packet",
    "exponential_sum_checksum",
    "circle_method_major_minor_split",
    "large_sieve_minor_arc_gate",
    "mu2_phi_inverse_unit_clock",
    "mobius_mu_over_n_tail",
    "raw_prime_count",
)


def carrier_bank() -> list[Carrier]:
    return [
        Carrier("raw_prime_count", (0, 0, 0, 0, 0, 0, 0, 3, 1), "all residues, packets, and smoothing debt"),
        Carrier("mobius_mu_over_n_tail", (0, 0, 1, 1, 0, 0, 0, 3, 2), "positive density and endpoint geometry"),
        Carrier("mu2_phi_inverse_unit_clock", (1, 0, 1, 2, 0, 0, 0, 3, 2), "prime powers, repeated-prime packets, and route labels"),
        Carrier("large_sieve_minor_arc_gate", (2, 1, 2, 2, 1, 1, 0, 2, 3), "major-arc local packet and boundary approach class"),
        Carrier("circle_method_major_minor_split", (2, 2, 2, 2, 2, 1, 1, 2, 3), "smoothing defect and exact exceptional packets"),
        Carrier("exponential_sum_checksum", (2, 2, 2, 2, 2, 3, 1, 2, 3), "which named packet owns a resonant checksum"),
        Carrier("smoothing_explicit_formula_packet", (3, 2, 2, 2, 4, 3, 2, 2, 4), "Kaczynski approach ambiguity unless labelled"),
        Carrier("kaczynski_boundary_approach", (3, 2, 2, 3, 4, 3, 4, 2, 4), "finite family label without repair ladder"),
        Carrier("analytic_sieve_clock_bridge", (3, 3, 3, 4, 4, 4, 4, 3, 4), "full packet labels and topology if used alone"),
        Carrier("labelled_repair_ladder_packet", (4, 4, 4, 4, 4, 4, 4, 4, 4), "nothing relevant in the audited bridge"),
    ]


def compare(a: Carrier, b: Carrier) -> str:
    aw = sum(1 for x, y in zip(a.vector, b.vector) if x > y)
    bw = sum(1 for x, y in zip(a.vector, b.vector) if y > x)
    if aw > bw:
        return a.name
    if bw > aw:
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


def score_histogram(names: list[str], edges: dict[tuple[str, str], str]) -> dict[int, int]:
    scores = Counter({name: 0 for name in names})
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
    for winner, loser in edges:
        graph[winner].append(loser)
        reverse[loser].append(winner)
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
        total = 1
        for w in reverse[v]:
            if w not in seen:
                total += rdfs(w)
        return total

    for name in reversed(order):
        if name not in seen:
            sizes.append(rdfs(name))
    return sorted(sizes, reverse=True)


def hamiltonian_count(names: list[str], edges: dict[tuple[str, str], str]) -> int:
    index = {name: i for i, name in enumerate(names)}
    n = len(names)
    dp: dict[tuple[int, int], int] = {(1 << i, i): 1 for i in range(n)}
    for mask_size in range(2, n + 1):
        for mask in range(1 << n):
            if mask.bit_count() != mask_size:
                continue
            for last in range(n):
                if not ((mask >> last) & 1):
                    continue
                prev_mask = mask ^ (1 << last)
                total = 0
                for prev in range(n):
                    if ((prev_mask >> prev) & 1) and beats(edges, names[prev], names[last]):
                        total += dp.get((prev_mask, prev), 0)
                if total:
                    dp[(mask, last)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, index[name]), 0) for name in names)


def tournament_readout() -> None:
    carriers = carrier_bank()
    names = [carrier.name for carrier in carriers]
    edges = edge_map(carriers)
    print("[3] Tournament Analysis")
    print("vertices=analytic proof clocks and repair packets, not runners/arcs")
    print(f"pairwise_observable={','.join(TOURNAMENT_FIELDS)}")
    print(f"tie_hamiltonian_path={'>'.join(TIE_PATH)}")
    print(f"score_hist={score_histogram(names, edges)}")
    print(f"directed_3cycles={directed_3cycles(names, edges)}")
    print(f"scc_sizes={scc_sizes(names, edges)}")
    print(f"hamiltonian_path_count={hamiltonian_count(names, edges)}")
    for carrier in carriers:
        print(f"  {carrier.name}: destroys={carrier.destroys}")
    print()


def proof_readout() -> None:
    print("[4] Proof-angle synthesis")
    print(
        "mu_over_n_role=cancellation/tail guard; it can remove a smoothed kernel "
        "term but cannot certify positive lonely density alone."
    )
    print(
        "mu2_phi_role=inverse primitive-unit capacity; it is valuable exactly "
        "when squarefree blindness is carried as a side-channel."
    )
    print(
        "large_sieve_role=minor-arc family majorant; admissible only after "
        "major-arc packet labels and exceptional fibers are retained."
    )
    print(
        "exponential_sum_role=checksum for resonance packets; an unexplained "
        "large sum must route to AP/GW, C27, K33, covering, or F7/THM-572 debt."
    )
    print(
        "smoothing_role=kernel homotopy ledger; changing the smoothing function "
        "is legal only when boundary defects are named."
    )
    print(
        "kaczynski_role=approach-class label for true-wide boundary ambiguity; "
        "ordinary approaches should open, ambiguous approaches must be classified."
    )
    print(
        "candidate_lemma=inside a fixed automatic/residue/fusion fiber, the first "
        "nonzero analytic clock among mu/n tail, mu2/phi capacity, large-sieve "
        "minor-arc budget, exponential-sum checksum, smoothing defect, and "
        "Kaczynski approach class either opens a strict component, descends to "
        "AP/GW/C27/K33/covering, is dual-annihilated by Fejer/Ramanujan/Haar, "
        "or emits F7/THM-572 residual debt."
    )
    print()


def assumption_challenge() -> None:
    print("[5] Assumption challenge")
    print(
        "Alternate vertices considered: runners, primes, denominators, residues, "
        "Fourier modes, smoothing kernels, exponential phases, boundary approach "
        "classes, local obstructions, and proof obligations."
    )
    print(
        "Chosen vertices are analytic proof clocks and repair packets.  The "
        "quotient preserves open/boundary status and theorem-route purity only "
        "when squarefree blindness, exact denominator, smoothing choice, "
        "exponential checksum, Kaczynski approach class, and packet labels remain "
        "attached or are discharged."
    )
    print(
        "Challenged assumption: a number-theoretic scalar such as sum mu(n)/n or "
        "sum mu(n)^2/phi(n) is not an LRC proof object.  It becomes useful only "
        "as a clock inside the HYP-3027 first-nonzero repair ladder."
    )


def main() -> None:
    print("LRC14 analytic sieve-clock bridge (codex-2026-06-26-S196)")
    print("source_threads=HYP-2982,HYP-2983,HYP-3024,HYP-3026,HYP-3027")
    print("external_prompt=Mobius/totient sums, large sieve/circle method, smoothing, Kaczynski")
    print()
    rows, prefix, limit = build_rows()
    print(f"rows={len(rows)} arithmetic_limit={limit}")
    row_readout(rows, prefix)
    arithmetic_checkpoint_readout(rows, prefix)
    fiber_report(rows)
    tournament_readout()
    proof_readout()
    assumption_challenge()


if __name__ == "__main__":
    main()
