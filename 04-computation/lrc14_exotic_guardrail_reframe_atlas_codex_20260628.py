#!/usr/bin/env python3
"""HYP-3408: exotic guardrail reframe atlas for LRC14.

This is a proof-facing router, not a proof.  It absorbs the prompt's
Ramanujan-Soldner, Sophie Germain, Hermite-Lindemann-Weierstrass, Krasner, and
Meissel-Mertens suggestions into the current HYP-3404/HYP-3406 program:

    residue word -> first failure -> owner/height repair -> exact proof exit.

The point is to test which named ideas become exact LRC14 sidecars and which
are only guardrails against illegal scalar shortcuts.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, gcd, log

import lrc14_three_coordinate_obstruction_exactness_codex_20260628 as collar


EULER_GAMMA = 0.57721566490153286060
SOLDNER_KNOWN_WINDOW = (1.45136923488338, 1.45136923488339)
AP = tuple(range(1, 14))


def frac_word(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def v_p(n: int, p: int) -> int:
    n = abs(n)
    if n == 0:
        return 99
    out = 0
    while n % p == 0:
        out += 1
        n //= p
    return out


def odd_part(n: int) -> int:
    n = abs(n)
    while n and n % 2 == 0:
        n //= 2
    return n


def replace_speed(old: int, new: int) -> tuple[int, ...]:
    return tuple(sorted(new if speed == old else speed for speed in AP))


def safe_mass(speeds: tuple[int, ...]) -> Fraction:
    return collar.interval_measure(collar.safe_open_components(speeds))


def li_series(x: float, terms: int = 80) -> float:
    """Principal logarithmic integral via li(x)=gamma+log|log x|+sum."""
    z = log(x)
    total = EULER_GAMMA + log(abs(z))
    for k in range(1, terms + 1):
        total += z**k / (k * factorial(k))
    return total


def soldner_root_from_series() -> float:
    lo, hi = 1.05, 2.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if li_series(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def sophie_germain_residual(a: int, b: int) -> int:
    left = a**4 + 4 * b**4
    plus = a * a + 2 * b * b + 2 * a * b
    minus = a * a + 2 * b * b - 2 * a * b
    return left - plus * minus


def sophie_germain_audit(limit: int = 9) -> dict[str, object]:
    residuals = []
    positive_factor_examples = []
    for a in range(-limit, limit + 1):
        for b in range(-limit, limit + 1):
            residuals.append(sophie_germain_residual(a, b))
            if a > 0 and b > 0 and len(positive_factor_examples) < 6:
                plus = a * a + 2 * b * b + 2 * a * b
                minus = a * a + 2 * b * b - 2 * a * b
                positive_factor_examples.append((a, b, plus, minus, plus * minus))
    live_height_examples = []
    for a, b in [(13, 1), (12, 1), (12, 2), (13, 2), (7, 2)]:
        plus = a * a + 2 * b * b + 2 * a * b
        minus = a * a + 2 * b * b - 2 * a * b
        live_height_examples.append((a, b, plus, minus, plus * minus))
    return {
        "checked_pairs": (2 * limit + 1) ** 2,
        "max_abs_residual": max(abs(x) for x in residuals),
        "factor_examples": positive_factor_examples,
        "live_height_examples": live_height_examples,
    }


def primes_upto(n: int) -> list[int]:
    sieve = bytearray(b"\x01") * (n + 1)
    if n >= 0:
        sieve[0:2] = b"\x00\x00"
    p = 2
    while p * p <= n:
        if sieve[p]:
            start = p * p
            sieve[start : n + 1 : p] = b"\x00" * (((n - start) // p) + 1)
        p += 1
    return [i for i in range(n + 1) if sieve[i]]


def meissel_mertens_table(limit: int = 200000) -> list[tuple[int, int, float, float]]:
    primes = primes_upto(limit)
    rows = []
    cursor = 0
    reciprocal_sum = 0.0
    for x in [14, 41, 4312, 200000]:
        while cursor < len(primes) and primes[cursor] <= x:
            reciprocal_sum += 1.0 / primes[cursor]
            cursor += 1
        rows.append((x, cursor, reciprocal_sum, reciprocal_sum - log(log(x))))
    return rows


@dataclass(frozen=True)
class MoveAudit:
    name: str
    old: int
    new: int
    status_hint: str

    @property
    def row(self) -> tuple[int, ...]:
        return replace_speed(self.old, self.new)

    @property
    def mass(self) -> Fraction:
        return safe_mass(self.row)

    @property
    def delta(self) -> int:
        return self.new - self.old

    @property
    def same_mod14(self) -> bool:
        return self.old % 14 == self.new % 14

    @property
    def same_mod7(self) -> bool:
        return self.old % 7 == self.new % 7

    @property
    def exit_status(self) -> str:
        return "boundary-tight" if self.mass == 0 else "strict-open"


def move_audits() -> list[MoveAudit]:
    return [
        MoveAudit("Goddyn-Wong hinge 12->24", 12, 24, "known AP/GW equality"),
        MoveAudit("same residue 12->26", 12, 26, "same mod14 but strict"),
        MoveAudit("same residue 2->16", 2, 16, "same mod14 but strict"),
        MoveAudit("unit height 13->27", 13, 27, "HYP-3405 first unit-height leak"),
        MoveAudit("near hinge 12->36", 12, 36, "minimum AP-collar strict floor"),
        MoveAudit("petal owner row 13->26", 13, 26, "HYP-3406 owner-support side"),
        MoveAudit("single owner row 1->26", 1, 26, "HYP-3406 owner-support side"),
        MoveAudit("single owner row 3->26", 3, 26, "HYP-3406 owner-support side"),
        MoveAudit("single owner row 5->26", 5, 26, "HYP-3406 owner-support side"),
    ]


@dataclass(frozen=True)
class Carrier:
    cid: str
    title: str
    proof_use: str
    preserved: tuple[str, ...]
    destroyed_if_naive: tuple[str, ...]
    first_test: str
    scores: dict[str, int]

    @property
    def total(self) -> int:
        weights = {
            "owner_repair": 4,
            "exact_finite_test": 4,
            "formal_theorem_hook": 3,
            "algebraic_exactness": 3,
            "global_tail": 2,
            "guardrail_value": 2,
            "risk_penalty": -3,
        }
        return sum(weights[key] * value for key, value in self.scores.items())


def carriers() -> list[Carrier]:
    return [
        Carrier(
            "X00",
            "Owner-residue-height first-failure theorem",
            "Use the exotic ideas only as gates around the HYP-3406 theorem route: residue failures must be repaired by owner support, height/flex, or named debt.",
            ("exit_status", "residue_word", "endpoint_owner_support", "height_flex"),
            ("raw named constants",),
            "Extend HYP-3406 beyond (72,20) until residue+owner_support first fails; classify the next leak.",
            {
                "owner_repair": 3,
                "exact_finite_test": 3,
                "formal_theorem_hook": 3,
                "algebraic_exactness": 2,
                "global_tail": 1,
                "guardrail_value": 2,
                "risk_penalty": 0,
            },
        ),
        Carrier(
            "X01",
            "Krasner owner-lift stability gate",
            "Treat p-adic closeness as useful only when the owner/contact root packet is stable; instability emits owner or height debt.",
            ("2-adic height", "7-adic residue", "endpoint_owner_support", "contact_roots"),
            ("theorem exit if owner support is forgotten",),
            "For same-residue lifts such as 12->26, 2->16, and 13->27, test whether owner/contact roots change before accepting a local p-adic lift.",
            {
                "owner_repair": 3,
                "exact_finite_test": 3,
                "formal_theorem_hook": 3,
                "algebraic_exactness": 2,
                "global_tail": 0,
                "guardrail_value": 2,
                "risk_penalty": 0,
            },
        ),
        Carrier(
            "X02",
            "Sophie-Germain quartic factor split",
            "Split quartic height/flex obstructions into two quadratic sign channels before feeding the covering-flex Hessian route.",
            ("quartic_height_obstruction", "two_quadratic_channels", "sign_payload"),
            ("endpoint owner", "off-grid floor"),
            "Apply the identity to the first unit-height leak 13->27 and to Worpitzky/Faulhaber fourth-difference packets.",
            {
                "owner_repair": 1,
                "exact_finite_test": 3,
                "formal_theorem_hook": 3,
                "algebraic_exactness": 3,
                "global_tail": 0,
                "guardrail_value": 2,
                "risk_penalty": 0,
            },
        ),
        Carrier(
            "X03",
            "Meissel-Mertens denominator entropy normalizer",
            "Use prime reciprocal/loglog scale only as a labelled tail normalizer for denominator sieves, never as a finite AP-collar certificate.",
            ("denominator_tail_scale", "smoothing_label", "exceptional_set_label"),
            ("prime powers", "endpoint owner", "exact period"),
            "Attach HYP-2982 Phi/G labels before using any Mertens or large-sieve tail estimate.",
            {
                "owner_repair": 0,
                "exact_finite_test": 1,
                "formal_theorem_hook": 2,
                "algebraic_exactness": 1,
                "global_tail": 3,
                "guardrail_value": 3,
                "risk_penalty": 1,
            },
        ),
        Carrier(
            "X04",
            "HLW no-scalar-shadow guardrail",
            "Use Hermite-Lindemann-Weierstrass as a warning: a transcendental analytic shadow cannot certify a finite rational LRC packet until translated into exact algebraic inequalities.",
            ("algebraic_shadow_audit", "transverse_payload_warning"),
            ("finite packet labels", "owner support"),
            "Before using exp/log/pi-like shadows, demand an exact rational sidecar or a named discarded-coordinate debt.",
            {
                "owner_repair": 0,
                "exact_finite_test": 1,
                "formal_theorem_hook": 2,
                "algebraic_exactness": 2,
                "global_tail": 1,
                "guardrail_value": 3,
                "risk_penalty": 1,
            },
        ),
        Carrier(
            "X05",
            "Ramanujan-Soldner zero-renormalization anchor",
            "Use the li zero as a metaphor for choosing a declared zero level in analytic boundary integrals; it supplies no LRC14 inequality by itself.",
            ("renormalization_basepoint", "sign_change_warning"),
            ("owner support", "height flex", "exact period"),
            "If a boundary integral appears, record its subtraction point and prove the sign with rational/interval data.",
            {
                "owner_repair": 0,
                "exact_finite_test": 1,
                "formal_theorem_hook": 1,
                "algebraic_exactness": 0,
                "global_tail": 2,
                "guardrail_value": 3,
                "risk_penalty": 1,
            },
        ),
        Carrier(
            "X06",
            "p-adic HLW claim firewall",
            "Do not accept a p-adic transcendence slogan as a proof step; require the actual convergence domain, theorem statement, and finite packet translation.",
            ("hypothesis_hygiene", "p_adic_convergence_domain"),
            ("LRC predicate if used as terminal proof"),
            "Reject SHA placeholders and route any p-adic exponential claim to an exact sidecar ledger.",
            {
                "owner_repair": 0,
                "exact_finite_test": 0,
                "formal_theorem_hook": 1,
                "algebraic_exactness": 1,
                "global_tail": 1,
                "guardrail_value": 3,
                "risk_penalty": 2,
            },
        ),
        Carrier(
            "X07",
            "Raw exotic-constant scalar",
            "Negative control: a bag of named constants with no preserved LRC predicate.",
            (),
            ("exit_status", "owner support", "height", "period", "finite witness"),
            "This should lose every pairwise comparison.",
            {
                "owner_repair": 0,
                "exact_finite_test": 0,
                "formal_theorem_hook": 0,
                "algebraic_exactness": 0,
                "global_tail": 0,
                "guardrail_value": 0,
                "risk_penalty": 3,
            },
        ),
    ]


def beats(left: Carrier, right: Carrier) -> bool:
    return (left.total, left.cid) > (right.total, right.cid)


def tournament_fingerprint(items: list[Carrier]) -> dict[str, object]:
    graph = {item.cid: set() for item in items}
    by_id = {item.cid: item for item in items}
    for left, right in combinations(items, 2):
        if beats(left, right):
            graph[left.cid].add(right.cid)
        else:
            graph[right.cid].add(left.cid)

    cycles = 0
    for a, b, c in combinations(graph, 3):
        if b in graph[a] and c in graph[b] and a in graph[c]:
            cycles += 1
        if c in graph[a] and b in graph[c] and a in graph[b]:
            cycles += 1

    hpaths = 0
    for order in permutations(items):
        if all(order[i + 1].cid in graph[order[i].cid] for i in range(len(order) - 1)):
            hpaths += 1

    priority = sorted(items, key=lambda item: (item.total, item.cid), reverse=True)
    return {
        "score_hist": dict(sorted(Counter(item.total for item in items).items())),
        "directed_3cycles": cycles,
        "hamiltonian_path_count": hpaths,
        "priority_path": [f"{item.cid}:{by_id[item.cid].title}" for item in priority],
    }


def print_soldner() -> None:
    root = soldner_root_from_series()
    lo, hi = SOLDNER_KNOWN_WINDOW
    print("RAMANUJAN-SOLDNER ZERO CHECK")
    print(f"  computed_root={root:.15f}")
    print(f"  li(root)={li_series(root):+.3e}")
    print(f"  in_known_window={lo < root < hi}")
    print("  LRC14 use: zero-level/renormalization guardrail only, not a packet proof.")
    print()


def print_sophie() -> None:
    audit = sophie_germain_audit()
    print("SOPHIE-GERMAIN QUARTIC FACTOR CHECK")
    print(f"  checked_pairs={audit['checked_pairs']}")
    print(f"  max_abs_identity_residual={audit['max_abs_residual']}")
    print("  live height/flex examples: a,b -> plus_factor, minus_factor, product")
    for a, b, plus, minus, product in audit["live_height_examples"]:
        print(f"    ({a:2d},{b:1d}) -> {plus:5d}, {minus:5d}, {product:8d}")
    print("  LRC14 use: exact split of quartic height debt into two quadratic channels.")
    print()


def print_mertens() -> None:
    print("MEISSEL-MERTENS PRIME-RECIPROCAL CHECK")
    print("  x       pi(x)   sum_p<=x 1/p     sum-loglog(x)")
    for x, pix, reciprocal_sum, constant_estimate in meissel_mertens_table():
        print(f"  {x:6d}  {pix:6d}   {reciprocal_sum:12.9f}   {constant_estimate:12.9f}")
    print("  LRC14 use: denominator tail entropy; HYP-2982 labels prime-power loss.")
    print()


def print_krasner_move_table() -> None:
    print("KRASNER-STYLE LOCAL-LIFT STABILITY TABLE")
    print("  move                         exit            mass          mod14 mod7  v2 old->new  v7 old->new  v2(delta) v7(delta)")
    for item in move_audits():
        print(
            f"  {item.name:28s} {item.exit_status:14s} {frac_word(item.mass):12s} "
            f"{str(item.same_mod14):5s} {str(item.same_mod7):5s} "
            f"{v_p(item.old, 2):2d}->{v_p(item.new, 2):2d}        "
            f"{v_p(item.old, 7):2d}->{v_p(item.new, 7):2d}        "
            f"{v_p(item.delta, 2):2d}        {v_p(item.delta, 7):2d}"
        )
    print("  Readout: same mod14 lifts can be strict, while 12->24 is tight but")
    print("  changes residue.  Krasner-style stability must include contact/owner")
    print("  roots, not raw p-adic closeness alone.")
    print()


def print_carriers() -> None:
    items = carriers()
    print("RANKED EXOTIC GUARDRAIL CARRIERS")
    for item in sorted(items, key=lambda carrier: (carrier.total, carrier.cid), reverse=True):
        print(f"  {item.cid} score={item.total:3d}  {item.title}")
        print(f"      use={item.proof_use}")
        print(f"      first_test={item.first_test}")
    print()

    fp = tournament_fingerprint(items)
    print("TOURNAMENT ANALYSIS")
    print("  vertices=proof carriers/guardrails, not constants, runners, or arcs")
    print("  pairwise_observable=weighted proof leverage with exact sidecar retention")
    print("  tie_hamiltonian_path=carrier id")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print("  priority_path=")
    for step in fp["priority_path"]:
        print(f"    {step}")
    print()


def print_assumption_challenge() -> None:
    print("ASSUMPTION CHALLENGE")
    print("  Rejected vertex set: named constants as proof vertices.")
    print("  Adopted vertex set: sidecar obligations and theorem exits.")
    print("  Preserved predicate: boundary-tight, strict-open, positive-Haar-open,")
    print("    unit-petal-named, AP/GW equality, or named debt.")
    print("  Destroyed by naive scalar quotient: endpoint owner support, height/flex,")
    print("    exact period, p-adic contact-root stability, and finite witness mass.")
    print("  Concrete next theorem target: residue+owner_support exactness on enlarged")
    print("    HYP-2963 banks, with Krasner-style lift instability and Sophie-Germain")
    print("    quartic factors used only after they name the lost coordinate.")


def main() -> None:
    print("HYP-3408 EXOTIC GUARDRAIL REFRAME ATLAS")
    print("=" * 78)
    print("status=synthesis/evidence; not an LRC14 proof")
    print("source=HYP-3404 creative atlas + HYP-3406 owner-support repair")
    print()
    print_soldner()
    print_sophie()
    print_krasner_move_table()
    print_mertens()
    print_carriers()
    print_assumption_challenge()


if __name__ == "__main__":
    main()
