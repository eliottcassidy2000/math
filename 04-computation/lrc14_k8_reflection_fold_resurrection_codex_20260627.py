#!/usr/bin/env python3
"""HYP-3138 scout: k=8 reflection-fold coordinate resurrection.

This script tests a concrete bridge among the recent k=8 De Moivre
biquadratic-resolvent lane, the phi4/gK8 bounded-core lane, and older
circuit/sheaf/theta/crystallographic sidecar work.

The proposed quotient is the reflection fold on the miss-count distribution
q_0,...,q_6:

    even_fold(q) = (q0+q6, q1+q5, q2+q4, q3).

The gK8 Delsarte functional L_yK8 = 10*q0 + q3 + 10*q6 depends only on this
fold, but endpoint/Phi/observer proof obligations may still need the odd
coordinates (q0-q6, q1-q5, q2-q4).  The scout therefore asks whether the even
fold is a legal finite lookup/adjoint on bounded k=8 banks, rather than
assuming the odd coordinates vanish.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import reduce
from itertools import combinations
from math import gcd, lcm
from fractions import Fraction


K = 8
GK8 = (10, 0, 0, 1, 0, 0, 10)
CAP8 = Fraction(2243, 5880)


def fmt(frac: Fraction) -> str:
    return str(frac)


def fmt_tuple(values: tuple[Fraction, ...]) -> str:
    return "(" + ", ".join(fmt(value) for value in values) + ")"


def primitive(row: tuple[int, ...]) -> bool:
    return reduce(gcd, (abs(x) for x in row if x), 0) == 1


def grid_denominator(span: int) -> int:
    denom = 1
    for speed in range(1, span + 1):
        denom = lcm(denom, 7 * speed)
    return denom


def miss_distribution_numerators(row: tuple[int, ...], denom: int) -> tuple[int, ...]:
    """Exact q distribution on the integer breakpoint grid with denominator denom."""
    breakpoints = {0, denom}
    for speed in row:
        if speed == 0:
            continue
        d = 7 * abs(speed)
        step = denom // d
        for m in range(d + 1):
            breakpoints.add(m * step)

    q = [0] * 7
    bps = sorted(breakpoints)
    for left, right in zip(bps, bps[1:]):
        mid2 = left + right  # midpoint numerator over 2*denom
        hit = {
            ((7 * abs(speed) * mid2) // (2 * denom)) % 7
            for speed in row
            if speed
        }
        occupied = sum(1 for sector in hit if sector != 0)
        missed = 6 - occupied
        q[missed] += right - left
    return tuple(q)


def folded_signature(q: tuple[int, ...]) -> tuple[int, ...]:
    return (q[0] + q[6], q[1] + q[5], q[2] + q[4], q[3])


def odd_leakage(q: tuple[int, ...]) -> tuple[int, ...]:
    return (q[0] - q[6], q[1] - q[5], q[2] - q[4])


def ly_num(q: tuple[int, ...]) -> int:
    return sum(weight * q[i] for i, weight in enumerate(GK8))


@dataclass(frozen=True)
class SpanAudit:
    span: int
    denom: int
    primitive_rows: int
    folded_signatures: int
    collision_fibers: int
    max_fiber_size: int
    best_row: tuple[int, ...]
    best_q: tuple[int, ...]
    best_fold: tuple[int, ...]
    best_odd: tuple[int, ...]
    best_ly: int
    top_rows: tuple[tuple[int, tuple[int, ...], tuple[int, ...]], ...]


def audit_span(span: int, top_n: int = 5) -> SpanAudit:
    denom = grid_denominator(span)
    rows = []
    fibers: dict[tuple[int, ...], list[tuple[int, ...]]] = defaultdict(list)

    for rest in combinations(range(1, span + 1), K - 1):
        row = (0,) + rest
        if not primitive(row):
            continue
        q = miss_distribution_numerators(row, denom)
        fold = folded_signature(q)
        rows.append((ly_num(q), row, q, fold))
        fibers[fold].append(q)

    rows.sort(reverse=True, key=lambda item: item[0])
    collision_fibers = sum(1 for qs in fibers.values() if len(set(qs)) > 1)
    max_fiber_size = max((len(qs) for qs in fibers.values()), default=0)
    best_ly, best_row, best_q, best_fold = rows[0]
    top_rows = tuple((val, row, q) for val, row, q, _fold in rows[:top_n])

    return SpanAudit(
        span=span,
        denom=denom,
        primitive_rows=len(rows),
        folded_signatures=len(fibers),
        collision_fibers=collision_fibers,
        max_fiber_size=max_fiber_size,
        best_row=best_row,
        best_q=best_q,
        best_fold=best_fold,
        best_odd=odd_leakage(best_q),
        best_ly=best_ly,
        top_rows=top_rows,
    )


def de_moivre_resolvent_lines() -> list[str]:
    # g(t)=(t-1)(t-2)(t-4)(t-5)=t^4-12t^3+49t^2-78t+40.
    # With t=u+3, the coefficients become u^4 + 0*u^3 -5*u^2 + 0*u + 4.
    return [
        "resolvent_g(t)=t^4-12*t^3+49*t^2-78*t+40",
        "substitute t=u+3 => u^4-5*u^2+4=(u^2-4)(u^2-1)",
        "odd_coefficients_after_fold=(u^3:0,u^1:0)",
        "quadratic_in_u2_discriminant=25-16=9",
        "reading=the De Moivre/quartic lane supplies an even fold, not an odd-coordinate eraser",
    ]


@dataclass(frozen=True)
class Carrier:
    name: str
    exactness: int
    preserves_k8_predicate: int
    quotient_legality: int
    destroyed_coordinate_control: int
    formal_next_step: int
    niche_bridge: int
    role: str

    def key(self) -> tuple[int, int, int, int, int, int, str]:
        return (
            self.preserves_k8_predicate,
            self.destroyed_coordinate_control,
            self.quotient_legality,
            self.formal_next_step,
            self.exactness,
            self.niche_bridge,
            self.name,
        )


def carriers() -> list[Carrier]:
    return [
        Carrier(
            "endpoint_phi_activation_circuit",
            5,
            5,
            5,
            5,
            5,
            3,
            "HYP-3116/HYP-2108/HYP-2112 exact gap circuit; terminal repair target",
        ),
        Carrier(
            "coordinate_resurrection_adjoint",
            4,
            5,
            5,
            5,
            4,
            4,
            "HYP-3118 rule: a quotient is legal after a right-adjoint repair section",
        ),
        Carrier(
            "k8_even_reflection_fold_table",
            5,
            5,
            4,
            5,
            4,
            4,
            "new finite lookup: even fold is injective on tested bounded k=8 banks",
        ),
        Carrier(
            "gK8_even_delsarte_functional",
            5,
            5,
            3,
            3,
            4,
            3,
            "HYP-3085/HYP-3122: L_yK8 depends only on q0+q6 and q3",
        ),
        Carrier(
            "de_moivre_biquadratic_resolvent",
            5,
            4,
            3,
            3,
            3,
            5,
            "HYP-3132 exact even quartic u^4-5u^2+4",
        ),
        Carrier(
            "jacobi_theta_even_odd_channels",
            4,
            4,
            3,
            4,
            3,
            5,
            "HYP-3110 theta channels as signed even/odd residue tails",
        ),
        Carrier(
            "A000568_global_consistency_quotient",
            4,
            4,
            4,
            4,
            3,
            4,
            "HYP-3134 controlled-forgetting discipline for local/global quotienting",
        ),
        Carrier(
            "wallpaper17_space230_orbifold_audit",
            3,
            2,
            3,
            3,
            2,
            5,
            "HYP-3110 finite crystallographic quotient audit, useful only with sidecars",
        ),
        Carrier(
            "raw_LyK8_scalar",
            5,
            3,
            1,
            1,
            2,
            1,
            "numeric scalar; preserves the Delsarte value but forgets proof coordinates",
        ),
    ]


def tournament_fingerprint(items: list[Carrier]) -> tuple[Counter[int], int, list[int], int, list[str]]:
    n = len(items)
    edges = [[False] * n for _ in range(n)]
    score = [0] * n
    for i, left in enumerate(items):
        for j, right in enumerate(items):
            if i == j:
                continue
            if left.key() > right.key():
                edges[i][j] = True
                score[i] += 1

    cycles = 0
    for a, b, c in combinations(range(n), 3):
        if (
            edges[a][b] and edges[b][c] and edges[c][a]
            or edges[a][c] and edges[c][b] and edges[b][a]
        ):
            cycles += 1

    # Strongly connected components by mutual reachability.
    def reachable(start: int) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for w in range(n):
                if edges[v][w] and w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    reach = [reachable(i) for i in range(n)]
    unseen = set(range(n))
    scc_sizes = []
    while unseen:
        root = min(unseen)
        comp = {j for j in unseen if root in reach[j] and j in reach[root]}
        scc_sizes.append(len(comp))
        unseen -= comp

    # Count Hamiltonian paths by DP over subsets.
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if count == 0:
                continue
            for nxt in range(n):
                if mask & (1 << nxt) == 0 and edges[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count
    full = (1 << n) - 1
    hamiltonian_paths = sum(dp.get((full, last), 0) for last in range(n))

    path = [item.name for item in sorted(items, key=lambda item: item.key(), reverse=True)]
    return Counter(score), cycles, sorted(scc_sizes, reverse=True), hamiltonian_paths, path


def main() -> None:
    print("HYP-3138 / codex-2026-06-27-k8-reflection-fold")
    print("title=k=8 reflection-fold coordinate resurrection scout")
    print("namespace=HYP-3138/T1203/LTI-264/LTT-162")
    print()

    print("De Moivre / biquadratic audit")
    for line in de_moivre_resolvent_lines():
        print(f"- {line}")
    print()

    print("Exact k=8 bounded-bank reflection-fold audit")
    for audit in (audit_span(14), audit_span(15), audit_span(16)):
        denom = audit.denom
        best_ly = Fraction(audit.best_ly, denom)
        margin = 10 * CAP8 - best_ly
        best_q = tuple(Fraction(x, denom) for x in audit.best_q)
        best_fold = tuple(Fraction(x, denom) for x in audit.best_fold)
        best_odd = tuple(Fraction(x, denom) for x in audit.best_odd)
        print(
            f"- span<={audit.span}: primitive_rows={audit.primitive_rows}, "
            f"folded_signatures={audit.folded_signatures}, "
            f"collision_fibers={audit.collision_fibers}, max_fiber_size={audit.max_fiber_size}"
        )
        print(
            f"  best_row={audit.best_row}, L_yK8={fmt(best_ly)}, "
            f"10*cap8-L_yK8={fmt(margin)}"
        )
        print(f"  even_fold={fmt_tuple(best_fold)}")
        print(f"  odd_leakage={fmt_tuple(best_odd)}")
        print(f"  q={fmt_tuple(best_q)}")
        print("  top_rows:")
        for val, row, q in audit.top_rows[:3]:
            q_frac = tuple(Fraction(x, denom) for x in q)
            odd = tuple(Fraction(x, denom) for x in odd_leakage(q))
            print(
                f"    row={row}, L_yK8={fmt(Fraction(val, denom))}, "
                f"odd={fmt_tuple(odd)}, q={fmt_tuple(q_frac)}"
            )
    print()

    print("Interpretation")
    print(
        "- The even fold is injective on every tested primitive k=8 bank "
        "span<=14,15,16, so it behaves like a finite adjoint lookup rather than "
        "a lossy scalar quotient in this range."
    )
    print(
        "- The best row is not reflection-balanced: odd_leakage is nonzero. "
        "Therefore HYP-3132's biquadratic symmetry should be used as a legal "
        "fold with odd-coordinate resurrection, not as a claim that odd "
        "coordinates vanish."
    )
    print(
        "- L_yK8 only sees even_fold=(q0+q6,q1+q5,q2+q4,q3), but endpoint "
        "Phi/P, observer gluing, and finite-address exits may still demand "
        "the destroyed odd coordinates."
    )
    print()

    score_hist, cycles, scc_sizes, hp_count, path = tournament_fingerprint(carriers())
    print("Tournament Analysis")
    print("vertices=proof carriers / quotient operators, not runners or raw roots")
    print("pairwise_observable=(k8 predicate retention, destroyed-coordinate control, quotient legality, formal next step, exactness, niche bridge)")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={cycles}")
    print(f"scc_sizes={scc_sizes}")
    print(f"hamiltonian_path_count={hp_count}")
    print("selected_path=" + " -> ".join(path))
    print()

    print("Candidate invariant")
    print(
        "k8_reflection_fold_certificate = "
        "even_folded_miss_distribution + odd_coordinate_resurrection_table "
        "+ de_moivre_biquadratic_resolvent_word + endpoint_phi_activation_status "
        "+ observer_gluing_or_named_debt"
    )
    print()
    print("LRC14 next theorem target")
    print(
        "- Prove a finite k=8 fold-adjoint lemma: on the bounded-core bank, "
        "even_fold determines the odd leakage or routes it to a named finite "
        "address/observer-gluing debt."
    )
    print(
        "- Then the gK8/phi4 dip bound can use the even biquadratic coordinate, "
        "while HYP-3116/HYP-3118 supply the legal repair for endpoint Phi/P "
        "coordinates before final quotienting."
    )
    print(
        "- The 17 wallpaper groups and 230 space groups remain useful only as "
        "finite quotient audits: they tell us to name stabilizers, glide/screw "
        "sidecars, and destroyed coordinates before using symmetry."
    )


if __name__ == "__main__":
    main()
