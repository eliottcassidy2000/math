#!/usr/bin/env python3
"""HYP-3421 scout: off-grid resonance transparency plus Rprime closure.

This is an exact-rational proof-route scout, not an LRC14 proof.

User prompt distilled:
  the useful optima sit off the 14-grid, with visible floors 1/12, 1/8,
  and 1/9; resonant speeds (multiples of 2 or 7, including 14Q tips)
  are safe at those optima; hence the feared resonant survivor should be
  treated as a full-optimum transparency ledger feeding the signed-SPEC
  Rprime >= 0.64178 certificate of HYP-3129.  HYP-3418 corrects the tempting
  coprime-to-14 reduction: the floor is still 2-adic because even speeds bind.

Tournament Analysis declaration:
  vertices: proof obligations / theorem carriers, not runners, arcs, raw
            residues, or scalar floor values;
  pairwise observable: majority over retained LRC predicate, finite exactness,
            resonance transparency, Rprime interface, formalization readiness,
            quotient legality, and failure guardrails;
  switch/gauge: orient A -> B when A wins more axes; ties use the declared
            Hamiltonian path;
  tie path: finite off-grid transparency, canonical 84m binding formula,
            2-adic even-speed descent, signed-SPEC Rprime chase, fiber PGF,
            edge witness, core rigidity, owner cuts, raw resonance worry.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Dict, Iterable, List, Sequence, Tuple


F = Fraction
THRESHOLD = F(1, 14)
UNIT_GRID = tuple(a for a in range(1, 14) if F(a, 14).denominator == 14 and __import__("math").gcd(a, 14) == 1)


def norm(x: F) -> F:
    r = x % 1
    return min(r, 1 - r)


def min_distance(speeds: Sequence[int], t: F) -> F:
    return min(norm(F(s) * t) for s in speeds)


def candidates(speeds: Sequence[int]) -> List[F]:
    """All pair-crossing and half-period candidates in [0,1/2]."""
    speeds = sorted(set(speeds))
    out = {F(1, 2)}

    for speed in speeds:
        k = 0
        while F(2 * k + 1, 2 * speed) <= F(1, 2):
            out.add(F(2 * k + 1, 2 * speed))
            k += 1

    for i, a in enumerate(speeds):
        for b in speeds[i + 1 :]:
            for denom in (a + b, b - a):
                if denom <= 0:
                    continue
                k = 1
                while F(k, denom) <= F(1, 2):
                    out.add(F(k, denom))
                    k += 1

    return sorted(out)


def exact_m(speeds: Sequence[int]) -> Tuple[F, F, Tuple[int, ...]]:
    best_m = F(-1)
    best_t = F(0)
    for t in candidates(speeds):
        value = min_distance(speeds, t)
        if value > best_m:
            best_m = value
            best_t = t
    active = tuple(s for s in sorted(set(speeds)) if norm(F(s) * best_t) == best_m)
    return best_m, best_t, active


@dataclass(frozen=True)
class TransparencyCase:
    name: str
    speeds: Tuple[int, ...]
    expected_m: F | None
    note: str


CASES: Tuple[TransparencyCase, ...] = (
    TransparencyCase(
        name="q12_seed_core",
        speeds=tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13]),
        expected_m=F(1, 12),
        note="The zeta-12 core C={1..11,13}; q=12 witness survives.",
    ),
    TransparencyCase(
        name="even_cover_floor_1_8",
        speeds=tuple(range(2, 15)),
        expected_m=F(1, 8),
        note="A covering-heavy bounded row; optimum is off the 14-grid.",
    ),
    TransparencyCase(
        name="multi_apex_floor_1_9",
        speeds=(1, 2, 3, 4, 5, 6, 7, 8, 14, 28, 42, 56, 70),
        expected_m=F(1, 9),
        note="A many-14Q row from the reduction probes; all 14Q tips are safe at optimum.",
    ),
    TransparencyCase(
        name="multi_apex_floor_1_8",
        speeds=(1, 2, 3, 4, 5, 6, 14, 28, 42, 56, 70, 84, 98),
        expected_m=F(1, 8),
        note="A second many-14Q row; resonant distances run from 1/8 to 1/2.",
    ),
    TransparencyCase(
        name="canonical_84m_m1",
        speeds=tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 84]),
        expected_m=F(7, 89),
        note="Canonical one-tail covering tower at m=1; binding pair (5,84).",
    ),
)


def resonant_speeds(speeds: Sequence[int]) -> Tuple[int, ...]:
    return tuple(s for s in sorted(set(speeds)) if s % 2 == 0 or s % 7 == 0)


def fourteen_q_speeds(speeds: Sequence[int]) -> Tuple[int, ...]:
    return tuple(s for s in sorted(set(speeds)) if s % 14 == 0)


def grid_killed_by_q14(q14: Sequence[int]) -> bool | None:
    if not q14:
        return None
    return all(norm(F(s) * F(a, 14)) == 0 for s in q14 for a in UNIT_GRID)


def audit_case(case: TransparencyCase) -> Dict[str, object]:
    m_value, t_star, active = exact_m(case.speeds)
    resonant = resonant_speeds(case.speeds)
    q14 = fourteen_q_speeds(case.speeds)
    distances = tuple((s, norm(F(s) * t_star)) for s in resonant)
    q14_distances = tuple((s, norm(F(s) * t_star)) for s in q14)
    min_res = min((d for _, d in distances), default=None)
    max_res = max((d for _, d in distances), default=None)
    min_q14 = min((d for _, d in q14_distances), default=None)
    off_grid = (14 * t_star).denominator != 1

    return {
        "case": case,
        "M": m_value,
        "t": t_star,
        "active": active,
        "off_grid": off_grid,
        "resonant": resonant,
        "q14": q14,
        "min_res": min_res,
        "max_res": max_res,
        "min_q14": min_q14,
        "unsafe_resonant": tuple((s, d) for s, d in distances if d < THRESHOLD),
        "grid_killed_by_q14": grid_killed_by_q14(q14),
        "expected_ok": case.expected_m is None or m_value == case.expected_m,
        "distances": distances,
        "q14_distances": q14_distances,
    }


def canonical_tower_rows(limit: int = 8) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    base = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13)
    for m in range(1, limit + 1):
        speed = 84 * m
        speeds = base + (speed,)
        predicted_t = F(35 * m + 2, 84 * m + 5)
        predicted_m = F(7 * m, 84 * m + 5)
        exact_value, exact_t, active = exact_m(speeds)
        resonant = resonant_speeds(speeds)
        min_res = min(norm(F(s) * predicted_t) for s in resonant)
        rows.append(
            {
                "m": m,
                "speed": speed,
                "predicted_t": predicted_t,
                "predicted_m": predicted_m,
                "exact_t": exact_t,
                "exact_m": exact_value,
                "active": active,
                "formula_ok": exact_t == predicted_t and exact_value == predicted_m,
                "min_resonant_at_predicted": min_res,
                "margin": predicted_m - THRESHOLD,
            }
        )
    return rows


AXES: Tuple[str, ...] = (
    "lrc_predicate_retention",
    "finite_exactness",
    "resonance_transparency",
    "rprime_interface",
    "core_rigidity_interface",
    "formalization_readiness",
    "quotient_legality",
    "failure_guard",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: Dict[str, int]
    preserves: str
    destroys: str
    next_hook: str


def carrier(
    name: str,
    scores: Sequence[int],
    preserves: str,
    destroys: str,
    next_hook: str,
) -> Carrier:
    if len(scores) != len(AXES):
        raise ValueError(name)
    return Carrier(name, dict(zip(AXES, scores)), preserves, destroys, next_hook)


CARRIERS: Tuple[Carrier, ...] = (
    carrier(
        "finite_offgrid_transparency_lemma",
        (10, 10, 10, 7, 6, 9, 8, 9),
        "exact witness predicate M(S)>1/14 for off-grid denominator families and resonant packets",
        "global all-packet coverage unless linked to HYP-2963/HYP-3266 residual classes",
        "classify every O12 residual into q-witness, 84m binding, 2-adic descent, positive-open owner packet, or Rprime row",
    ),
    carrier(
        "canonical_84m_binding_formula",
        (9, 10, 10, 5, 5, 9, 7, 9),
        "infinite resonant one-tail tower C union {84m}, with exact M=7m/(84m+5)",
        "multi-tail interactions and owner-side collisions",
        "use the (5,84m) binding flank as the model for resonant tips becoming guard rails",
    ),
    carrier(
        "two_adic_even_speed_descent",
        (9, 8, 8, 8, 6, 8, 8, 9),
        "HYP-3418 correction that even speeds are the binding floor obstruction and halve under u=2t",
        "odd-speed coupling and terminal all-packet descent proof are not yet closed",
        "turn full-optimum transparency into a 2-adic descent theorem for even-R union 14Q",
    ),
    carrier(
        "signed_SPEC_Rprime_constant_chase",
        (9, 8, 7, 10, 5, 8, 8, 8),
        "the HYP-3129 conditional decorrelation floor Rprime >= 0.64178 on representative rows",
        "closed-form all-row low-frequency bound is still a finite constant chase",
        "prove uniform SPEC_low and Parseval-tail constants over all legal (R,Q)",
    ),
    carrier(
        "fiber_PGF_conditional_moment",
        (8, 9, 6, 9, 5, 8, 9, 7),
        "Rprime as E[N_R | Q]/E[N_R] over 14 sheets",
        "individual sheet identity if collapsed to one scalar Rprime",
        "attach fiber_pgf_word and Q_masked_fiber_pgf_word to each edge-floor packet",
    ),
    carrier(
        "edge_witness_tail_tip_packet",
        (8, 8, 7, 8, 5, 8, 10, 8),
        "R-safe tail, Q-safe tip, four-sector deck, and deletion children",
        "analytic constant values unless joined to SPEC/PGF",
        "route each transparency row through HYP-3124/HYP-3125 endpoint-child schema",
    ),
    carrier(
        "unit_core_equioscillation_rigidity",
        (8, 8, 5, 4, 10, 7, 8, 8),
        "six AP/GW unit contacts in three antipodal binding pairs",
        "covering-layer off-grid floors and nonunit magnitude flex",
        "finish the AP/GW-only finite equioscillation census with blind sidecars retained",
    ),
    carrier(
        "owner_cut_terminal_router",
        (7, 8, 6, 5, 7, 8, 9, 10),
        "endpoint-owner cut separation for mixed theorem-exit fibers",
        "raw analytic Rprime constants unless the row reaches the floor branch",
        "use HYP-3414 owner transversals as terminal routing after transparency/floor classification",
    ),
    carrier(
        "raw_resonant_survivor_worry",
        (2, 2, 2, 1, 2, 2, 1, 3),
        "a useful alarm that resonant speeds can kill the 14-grid",
        "the off-grid witness, exact distances, and 2-adic/Rprime floor mechanism",
        "retire as a theorem carrier; keep only as a negative-control alarm",
    ),
)

TIE_PATH = tuple(c.name for c in CARRIERS)


def compare(a: Carrier, b: Carrier) -> int:
    aw = bw = 0
    for axis in AXES:
        if a.scores[axis] > b.scores[axis]:
            aw += 1
        elif b.scores[axis] > a.scores[axis]:
            bw += 1
    if aw != bw:
        return 1 if aw > bw else -1
    return 1 if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else -1


def tournament_edges(carriers: Sequence[Carrier]) -> Dict[str, set[str]]:
    edges = {c.name: set() for c in carriers}
    for a, b in combinations(carriers, 2):
        if compare(a, b) > 0:
            edges[a.name].add(b.name)
        else:
            edges[b.name].add(a.name)
    return edges


def count_directed_3cycles(names: Sequence[str], edges: Dict[str, set[str]]) -> int:
    total = 0
    for a, b, c in combinations(names, 3):
        if b in edges[a] and c in edges[b] and a in edges[c]:
            total += 1
        if c in edges[a] and b in edges[c] and a in edges[b]:
            total += 1
    return total


def hamiltonian_path_count(names: Sequence[str], edges: Dict[str, set[str]]) -> int:
    n = len(names)
    index = {name: i for i, name in enumerate(names)}
    dp: Dict[Tuple[int, int], int] = {}
    for name in names:
        dp[(1 << index[name], index[name])] = 1
    for mask in range(1 << n):
        for last in range(n):
            current = dp.get((mask, last), 0)
            if not current:
                continue
            last_name = names[last]
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                nxt_name = names[nxt]
                if nxt_name in edges[last_name]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + current
    return sum(dp.get(((1 << n) - 1, last), 0) for last in range(n))


def selected_hamiltonian_path(names: Sequence[str], edges: Dict[str, set[str]]) -> List[str]:
    remaining = set(names)
    path: List[str] = []
    current = max(remaining, key=lambda name: len(edges[name]))
    path.append(current)
    remaining.remove(current)
    while remaining:
        candidates = [name for name in remaining if name in edges[current]]
        if not candidates:
            candidates = list(remaining)
        current = max(candidates, key=lambda name: (len(edges[name]), -TIE_PATH.index(name)))
        path.append(current)
        remaining.remove(current)
    return path


def print_fraction(value: F | None) -> str:
    if value is None:
        return "n/a"
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    print("HYP-3421 / codex-2026-06-28")
    print("LRC14 off-grid resonance transparency + Rprime closure scout")
    print("status: exact-rational synthesis; not a proof")
    print()

    print("1. ASSUMPTION CHALLENGE")
    print(
        "alternate vertex sets considered=runners,gaps,fixed circle sections,"
        "section boundaries,wall-crossing events,residues,cover arcs,Fourier modes,"
        "witness times,off-grid cells,14Q tips,PGF coefficients,proof obligations"
    )
    print("chosen tournament vertices=proof obligations / theorem carriers")
    print("preserved predicate=existence of a time t with all speeds at distance >= 1/14")
    print(
        "destroyed by raw resonance quotient=off-grid witness location, endpoint owners,"
        " tail/tip child packets, sheet counts, and SPEC low/tail constants"
    )
    print("challenged assumption=resonant speeds are the obstruction; on the checked packets they are transparent at the off-grid optimum")
    print()

    print("2. EXACT OFF-GRID TRANSPARENCY CASES")
    print(
        "case | |S| | M | t_star | M_den | off_14_grid | active | resonant_min..max | q14_min | grid_killed_by_q14 | expected_ok"
    )
    audits = [audit_case(case) for case in CASES]
    for audit in audits:
        case = audit["case"]
        print(
            f"{case.name} | {len(case.speeds)} | {print_fraction(audit['M'])} | "
            f"{print_fraction(audit['t'])} | {audit['M'].denominator} | "
            f"{audit['off_grid']} | {audit['active']} | "
            f"{print_fraction(audit['min_res'])}..{print_fraction(audit['max_res'])} | "
            f"{print_fraction(audit['min_q14'])} | {audit['grid_killed_by_q14']} | {audit['expected_ok']}"
        )
        if audit["unsafe_resonant"]:
            print(f"  UNSAFE_RESONANT={audit['unsafe_resonant']}")
        print(f"  note={case.note}")
    print()
    global_min = min(audit["min_res"] for audit in audits if audit["min_res"] is not None)
    global_max = max(audit["max_res"] for audit in audits if audit["max_res"] is not None)
    print(
        "readout: resonant speeds in these off-grid witnesses have exact distances "
        f"{print_fraction(global_min)}..{print_fraction(global_max)}, all >= 1/14; "
        "the visible floor denominators are 12, 8, and 9."
    )
    print()

    print("3. CANONICAL 84m RESONANT TOWER")
    print("m | speed | predicted_t | exact_t | predicted_M | exact_M | active | min_resonant | margin_over_1_14 | formula_ok")
    for row in canonical_tower_rows(8):
        print(
            f"{row['m']} | {row['speed']} | {print_fraction(row['predicted_t'])} | "
            f"{print_fraction(row['exact_t'])} | {print_fraction(row['predicted_m'])} | "
            f"{print_fraction(row['exact_m'])} | {row['active']} | "
            f"{print_fraction(row['min_resonant_at_predicted'])} | "
            f"{print_fraction(row['margin'])} | {row['formula_ok']}"
        )
    print(
        "tower_readout: the most resonant one-tail family C union {84m} is not an obstruction; "
        "the 14Q runner is one binding flank and the exact margin is 7m/(84m+5)-1/14."
    )
    print()

    print("4. PROOF-FACING OBLIGATIONS")
    obligations = (
        "O12a finite_offgrid_transparency: prove every O12 residual reaches one of the checked denominator floors, the 84m formula, or a named positive-open packet.",
        "O12b divisor_grid_localization: 14Q tips kill the 14-grid but are safe on the certified off-grid cells.",
        "O12c signed_SPEC_constant_chase: turn HYP-3129's per-row exact low part plus Parseval tail into a closed-form all-packet Rprime >= c theorem.",
        "O12d fiber_PGF_translation: state Rprime as E[N_R | Q]/E[N_R] before scalarizing.",
        "O12e two_adic_even_descent: incorporate HYP-3418's correction that coprime-to-14 reduction fails and even speeds drive the covering floor.",
        "O15 finite_equioscillation_core: finish AP/GW rigidity from the three antipodal unit binding pairs plus blind sidecars.",
        "GLUE terminal_factorization: combine finite transparency or L=Rprime*R_safe*Q_lonely with terminal owner/state-lift routers.",
    )
    for item in obligations:
        print(f"- {item}")
    print()

    print("5. TOURNAMENT ANALYSIS OVER PROOF CARRIERS")
    names = [c.name for c in CARRIERS]
    edges = tournament_edges(CARRIERS)
    score_hist = Counter(len(edges[name]) for name in names)
    print(f"vertices={len(names)} edges={sum(len(v) for v in edges.values())}")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={count_directed_3cycles(names, edges)}")
    print(f"hamiltonian_path_count={hamiltonian_path_count(names, edges)}")
    print("selected_hamiltonian_path=" + " -> ".join(selected_hamiltonian_path(names, edges)))
    print()
    for carrier_obj in CARRIERS:
        print(
            f"- {carrier_obj.name}: preserves={carrier_obj.preserves}; "
            f"destroys_if_scalarized={carrier_obj.destroys}; next={carrier_obj.next_hook}"
        )
    print()

    print("6. SYNTHESIS")
    print(
        "The resonant-survivor worry should be downgraded from a separate analytic obstruction "
        "to a finite transparency ledger plus a named constant chase.  Resonance attacks the "
        "14-grid core; the tested optima move off-grid, where every resonant speed in the packet "
        "is already at distance at least 1/12, 1/9, 1/8, or 7m/(84m+5).  HYP-3418 then "
        "blocks the naive coprime-to-14 reduction: odd runners alone prefer t=1/2, where even "
        "runners die.  After finite transparency, the remaining floor is the 2-adic even-speed "
        "descent plus signed SPEC / fiber PGF / edge-witness Rprime gluing.  This is not a "
        "completed proof; it identifies the strongest completion obligation as the all-packet "
        "transparency classification joined to the 2-adic descent and HYP-3129's closed-form "
        "Rprime constant chase."
    )


if __name__ == "__main__":
    main()
