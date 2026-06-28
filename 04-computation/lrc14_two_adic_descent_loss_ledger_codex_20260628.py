#!/usr/bin/env python3
"""HYP-3428: 2-adic descent loss ledger for the LRC14 floor route.

S259/HYP-3418 corrected the critical path: the remaining covering-floor
inequality is controlled by even speeds and 2-adic descent, not by the apex-7
census/equioscillation story.  This scout makes that correction operational on
the exact HYP-3410/HYP-3417 AP-collar mixed-fiber substrate.  Those rows are
q-witness/non-covering proxies, not arbitrary covering rows; the point is to
name what the HYP-3422/HYP-3425 covering-floor halving route is not allowed to
forget.

It reconstructs the AP-collar row names when possible, computes the exact LRC
maximin value for those finite speed sets, and records:

* failure of the naive odd/coprime witness at t=1/2;
* whether the exact optimum is carried by even binders;
* whether the even child under u=2t carries the same bottleneck;
* whether v2 profiles alone separate theorem exits inside HYP-3410 fibers;
* how the HYP-3417 owner-current cuts expose an even-hinge label.

This is not an LRC14 proof.  It is a proof-angle scout: a quotient using
2-adic descent must carry a loss ledger for odd blockers, even children, and
owner-current hinge labels before it can feed the decorrelation floor theorem.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
import math
import re

import lrc14_bring_sc_bdh_menger_charal_recursion_codex_20260628 as h3410
import lrc14_owner_cut_dual_certificate_synthesis_codex_20260628 as h3417


BASE_AP = tuple(range(1, 14))
THRESHOLD = Fraction(1, 14)


def v2(n: int) -> int:
    out = 0
    while n and n % 2 == 0:
        out += 1
        n //= 2
    return out


def frac_part(x: Fraction) -> Fraction:
    return x - math.floor(x)


def circle_dist(s: int, t: Fraction) -> Fraction:
    y = frac_part(Fraction(s) * t)
    return min(y, 1 - y)


def reconstruct_speed_set(name: str) -> tuple[int, ...] | None:
    """Reconstruct AP-collar speed sets from HYP-3410 row names.

    Rows with opaque names such as P10+GW are deliberately left unreconstructed;
    the ledger keeps them as named debt instead of guessing.
    """

    base = set(BASE_AP)
    m2 = re.search(r"two drop\((\d+),\s*(\d+)\)->add\((\d+),\s*(\d+)\)", name)
    if m2:
        old_a, old_b, new_a, new_b = map(int, m2.groups())
        return tuple(sorted((base - {old_a, old_b}) | {new_a, new_b}))

    m1 = re.search(r"(\d+)->(\d+)", name)
    if m1:
        old, new = map(int, m1.groups())
        return tuple(sorted((base - {old}) | {new}))

    return None


def q_witness_status(speeds: tuple[int, ...]) -> tuple[str, int | None]:
    """Return the first q<=14 missed by all speeds, or covering if none."""

    for q in range(2, 15):
        if all(s % q != 0 for s in speeds):
            return ("q-witness", q)
    return ("covering", None)


def affine_piece(s: int, mid: Fraction) -> tuple[int, int]:
    x = Fraction(s) * mid
    k = math.floor(x)
    y = x - k
    if y <= Fraction(1, 2):
        return s, -k
    return -s, k + 1


def exact_lrc_maximin(speeds: tuple[int, ...]) -> tuple[Fraction, Fraction, tuple[int, ...]]:
    """Exact max of min_s ||s t|| over t in [0,1].

    The distance functions are piecewise affine with breakpoints m/(2s).  On
    each common cell, the lower envelope maximum is attained at a cell endpoint
    or at an intersection of two affine pieces.
    """

    breaks = {Fraction(0), Fraction(1)}
    for s in speeds:
        for m in range(0, 2 * s + 1):
            breaks.add(Fraction(m, 2 * s))
    points = sorted(breaks)
    candidates: set[Fraction] = set(points)

    for a, b in zip(points, points[1:]):
        if a == b:
            continue
        mid = (a + b) / 2
        lines = [affine_piece(s, mid) for s in speeds]
        for (m1, c1), (m2, c2) in combinations(lines, 2):
            if m1 == m2:
                continue
            t = Fraction(c2 - c1, m1 - m2)
            if a <= t <= b:
                candidates.add(t)

    best_t = Fraction(0)
    best_val = Fraction(-1)
    for t in candidates:
        val = min(circle_dist(s, t) for s in speeds)
        if val > best_val or (val == best_val and t < best_t):
            best_val = val
            best_t = t

    binders = tuple(s for s in speeds if circle_dist(s, best_t) == best_val)
    return best_val, best_t, binders


def hist_v2(speeds: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(v2(s) for s in speeds).items()))


def layer_survival(speeds: tuple[int, ...]) -> tuple[int, ...]:
    max_depth = max(v2(s) for s in speeds)
    return tuple(sum(1 for s in speeds if v2(s) >= depth) for depth in range(1, max_depth + 1))


def even_child(speeds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(s // 2 for s in speeds if s % 2 == 0))


@dataclass(frozen=True)
class RowLedger:
    fiber: str
    name: str
    exit: str
    speeds: tuple[int, ...] | None
    q_status: tuple[str, int | None] | None
    value: Fraction | None
    t_star: Fraction | None
    binders: tuple[int, ...]
    binder_parity: tuple[int, int]
    v2_hist: tuple[tuple[int, int], ...]
    survival: tuple[int, ...]
    half_min: Fraction | None
    child_at_u: Fraction | None
    odd_at_t: Fraction | None
    loss_class: str


def audit_row(fiber: str, row: h3410.Row) -> RowLedger:
    speeds = reconstruct_speed_set(row.name)
    if speeds is None:
        return RowLedger(
            fiber=fiber,
            name=row.name,
            exit=row.exit,
            speeds=None,
            q_status=None,
            value=None,
            t_star=None,
            binders=(),
            binder_parity=(0, 0),
            v2_hist=(),
            survival=(),
            half_min=None,
            child_at_u=None,
            odd_at_t=None,
            loss_class="opaque_named_row_debt",
        )

    value, t_star, binders = exact_lrc_maximin(speeds)
    even_binders = sum(1 for s in binders if s % 2 == 0)
    odd_binders = len(binders) - even_binders
    child = even_child(speeds)
    u = frac_part(2 * t_star)
    child_at_u = min((circle_dist(s, u) for s in child), default=Fraction(1, 2))
    odd_speeds = tuple(s for s in speeds if s % 2)
    odd_at_t = min((circle_dist(s, t_star) for s in odd_speeds), default=Fraction(1, 2))
    half_min = min(circle_dist(s, Fraction(1, 2)) for s in speeds)

    if even_binders and child_at_u == value and odd_at_t > value:
        loss_class = "pure_even_child_bottleneck"
    elif even_binders and child_at_u == value:
        loss_class = "even_child_with_odd_coupling"
    elif even_binders:
        loss_class = "even_binder_but_child_not_exact"
    elif odd_binders:
        loss_class = "odd_binder_after_even_shift"
    else:
        loss_class = "no_binder"

    return RowLedger(
        fiber=fiber,
        name=row.name,
        exit=row.exit,
        speeds=speeds,
        q_status=q_witness_status(speeds),
        value=value,
        t_star=t_star,
        binders=binders,
        binder_parity=(odd_binders, even_binders),
        v2_hist=hist_v2(speeds),
        survival=layer_survival(speeds),
        half_min=half_min,
        child_at_u=child_at_u,
        odd_at_t=odd_at_t,
        loss_class=loss_class,
    )


def fmt_fraction(x: Fraction | None) -> str:
    if x is None:
        return "NA"
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def fmt_q_status(status: tuple[str, int | None] | None) -> str:
    if status is None:
        return "opaque"
    kind, q = status
    return kind if q is None else f"{kind}:{q}"


def mixed_key_count(ledgers: list[RowLedger], key_name: str) -> int:
    fibers: dict[object, set[str]] = defaultdict(set)
    for ledger in ledgers:
        key = getattr(ledger, key_name)
        fibers[key].add(ledger.exit)
    return sum(1 for exits in fibers.values() if len(exits) > 1)


FEATURE_WEIGHTS = {
    "critical_path_floor": 14,
    "two_adic_descent": 13,
    "exact_row_maximin": 11,
    "loss_ledger": 11,
    "owner_current_bridge": 9,
    "inductive_child": 8,
    "q_witness_exit": 7,
    "off_path_guard": 6,
    "raw_scalar_shadow": 1,
}


@dataclass(frozen=True)
class Carrier:
    name: str
    keeps: frozenset[str]
    loses: frozenset[str]
    priority: int

    def score(self) -> int:
        kept = sum(FEATURE_WEIGHTS[x] for x in self.keeps)
        lost = sum(max(1, FEATURE_WEIGHTS[x] // 3) for x in self.loses)
        return kept - lost + self.priority


def carriers() -> list[Carrier]:
    return [
        Carrier(
            "two_adic_descent_loss_ledger",
            frozenset({"critical_path_floor", "two_adic_descent", "exact_row_maximin", "loss_ledger"}),
            frozenset(),
            6,
        ),
        Carrier(
            "even_child_induction_packet",
            frozenset({"two_adic_descent", "inductive_child", "loss_ledger"}),
            frozenset({"owner_current_bridge"}),
            5,
        ),
        Carrier(
            "owner_current_even_hinge",
            frozenset({"owner_current_bridge", "loss_ledger", "critical_path_floor"}),
            frozenset({"inductive_child"}),
            4,
        ),
        Carrier(
            "odd_half_witness_failure_gate",
            frozenset({"critical_path_floor", "loss_ledger"}),
            frozenset({"two_adic_descent"}),
            3,
        ),
        Carrier(
            "q_witness_noncovering_exit",
            frozenset({"q_witness_exit", "exact_row_maximin"}),
            frozenset({"critical_path_floor"}),
            2,
        ),
        Carrier(
            "apex7_census_offpath_guard",
            frozenset({"off_path_guard", "q_witness_exit"}),
            frozenset({"critical_path_floor", "two_adic_descent"}),
            1,
        ),
        Carrier(
            "raw_coprime_to_14_reduction",
            frozenset({"raw_scalar_shadow"}),
            frozenset({"critical_path_floor", "two_adic_descent", "loss_ledger", "inductive_child"}),
            0,
        ),
    ]


def tournament_summary(nodes: list[Carrier]) -> tuple[Counter[int], int, int, list[Carrier]]:
    ordered = sorted(nodes, key=lambda x: (-x.score(), x.name))
    score_hist = Counter(node.score() for node in nodes)
    # The score order orients every pair, so cycles are impossible.
    directed_3cycles = 0
    hamiltonian_path_count = 1 if len({n.score() for n in nodes}) == len(nodes) else 0
    return score_hist, directed_3cycles, hamiltonian_path_count, ordered


def main() -> None:
    ledgers: list[RowLedger] = []
    for fiber, rows in h3410.FIBERS.items():
        for row in rows:
            ledgers.append(audit_row(fiber, row))

    reconstructed = [x for x in ledgers if x.speeds is not None]
    opaque = [x for x in ledgers if x.speeds is None]
    print("HYP-3428 2-adic descent loss ledger scout")
    print("=" * 78)
    print("substrate=HYP-3410 AP-collar mixed fibers + HYP-3417 owner-current cuts")
    print("upstream=HYP-3418 says the covering floor is 2-adic / even-speed controlled")
    print("scope=q-witness AP-collar proxy; pair with HYP-3422/HYP-3425 covering lift")
    print("goal=name the loss ledger needed before a halving descent becomes a theorem")
    print()

    print("## Exact reconstructed AP-collar rows")
    print(f"rows_total={len(ledgers)} reconstructed={len(reconstructed)} opaque_named_debt={len(opaque)}")
    print(f"threshold=1/14")
    print()
    for ledger in reconstructed:
        odd_binders, even_binders = ledger.binder_parity
        print(f"### {ledger.fiber} :: {ledger.name}")
        print(f"exit={ledger.exit}")
        print(f"speeds={ledger.speeds}")
        print(f"q_status={fmt_q_status(ledger.q_status)}")
        print(
            "M="
            f"{fmt_fraction(ledger.value)} t*={fmt_fraction(ledger.t_star)} "
            f"above_1_14={ledger.value is not None and ledger.value > THRESHOLD}"
        )
        print(
            f"binders={ledger.binders} binder_parity=(odd:{odd_binders}, even:{even_binders})"
        )
        print(f"v2_hist={ledger.v2_hist} layer_survival={ledger.survival}")
        print(
            "half_witness_min="
            f"{fmt_fraction(ledger.half_min)} child_at_u=2t*:{fmt_fraction(ledger.child_at_u)} "
            f"odd_at_t*={fmt_fraction(ledger.odd_at_t)}"
        )
        print(f"loss_class={ledger.loss_class}")
        print()

    print("## Fiber-level loss ledger")
    for fiber in h3410.FIBERS:
        sub = [x for x in ledgers if x.fiber == fiber]
        sub_recon = [x for x in sub if x.speeds is not None]
        selected = h3417.audit_fiber(fiber, h3410.FIBERS[fiber]).selected
        cut_group_hist = Counter(h3410.owner_group(label) for label in selected.cut)
        exit_hist = Counter(x.exit for x in sub)
        loss_hist = Counter(x.loss_class for x in sub)
        binder_even_rows = sum(1 for x in sub_recon if x.binder_parity[1] > 0)
        half_fail_rows = sum(1 for x in sub_recon if x.half_min == 0)
        print(f"### {fiber}")
        print(f"exit_hist={dict(exit_hist)}")
        print(f"selected_owner_current={selected.mode} cut={selected.cut} groups={dict(cut_group_hist)}")
        print(f"loss_hist={dict(loss_hist)}")
        print(f"reconstructed_rows={len(sub_recon)} half_witness_failure_rows={half_fail_rows}")
        print(f"rows_with_even_binders={binder_even_rows}")
        print(
            f"v2_hist_mixed_fibers={mixed_key_count(sub_recon, 'v2_hist')} "
            f"survival_mixed_fibers={mixed_key_count(sub_recon, 'survival')}"
        )
        print()

    all_half_fail = sum(1 for x in reconstructed if x.half_min == 0)
    all_even_binders = sum(1 for x in reconstructed if x.binder_parity[1] > 0)
    all_child_exact = sum(
        1 for x in reconstructed if x.value is not None and x.child_at_u == x.value
    )
    q_hist = Counter(fmt_q_status(x.q_status) for x in reconstructed)
    print("## Synthesis")
    print(f"q_status_hist={dict(q_hist)}")
    print(f"half_witness_failure_rows={all_half_fail}/{len(reconstructed)}")
    print(f"rows_with_even_binders={all_even_binders}/{len(reconstructed)}")
    print(f"rows_where_even_child_carries_M_at_u=2t*={all_child_exact}/{len(reconstructed)}")
    print(f"opaque_named_rows={[x.name for x in opaque]}")
    print(
        "main_readout=the t=1/2 odd/coprime witness is exactly the wrong witness "
        "whenever any even speed is present; the useful descent object is the "
        "even-child packet plus a ledger of odd blockers and owner-current hinge labels."
    )
    print(
        "frontier_bridge=HYP-3417's hardest current uses cut {2:g2,11:g1,13:g1}, "
        "so the first size-3 owner debt already contains one even-cover label."
    )
    print(
        "covering_proof_interface=HYP-3425 Helly/odd-bad-core target plus this "
        "loss ledger for owner-current/even-hinge/off-grid/state-lift debt."
    )
    print()

    score_hist, cycles, hp_count, path = tournament_summary(carriers())
    print("## Tournament Analysis")
    print("vertices=proof carriers and descent ledgers, not runners")
    print("pairwise_observable=retained critical-path floor data minus forgotten descent debt")
    print("switch_gauge=higher weighted retained payload; ties by name")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={cycles}")
    print(f"hamiltonian_path_count={hp_count}")
    print("priority_hamiltonian_path=")
    for node in path:
        print(
            f"  {node.name}: score={node.score()} "
            f"keeps={sorted(node.keeps)} loses={sorted(node.loses)}"
        )
    print()

    print("## Assumption Challenge")
    print(
        "considered_vertices=runners, odd speeds, even speeds, v2 layers, owner labels, "
        "halved child packets, q-witness exits, off-grid witnesses, and proof obligations"
    )
    print(
        "chosen_vertices=proof carriers/descent ledgers; this preserves the covering-floor "
        "predicate and destroys raw row order"
    )
    print(
        "challenged_assumption=the coprime-to-14 or odd-speed subproblem can supply a witness "
        "without recording how even speeds move the optimum away from t=1/2"
    )


if __name__ == "__main__":
    main()
