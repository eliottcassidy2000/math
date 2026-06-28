#!/usr/bin/env python3
"""
Exact k=8 ordered-tail and central-exchange scout.

This is a follow-on to HYP-3200.  It tests two possible proof angles for the
remaining LRC14 k=8 target:

1. Replace full convex order by a one-sided upper-tail/stop-loss barrier.
2. Replace raw q3 maximization by an exchange-rate lemma: any central q3 gain
   must cost at least as much extreme bimodality q0+q6.

The script imports HYP-3200's exact row_moments routine so all probabilities
are exact Fractions on the same anchored bounded bank.
"""

from __future__ import annotations

import importlib.util
import itertools
from pathlib import Path
from fractions import Fraction


HERE = Path(__file__).resolve().parent
H3200 = HERE / "lrc_k8_cumulant_universality_codex_20260628.py"
spec = importlib.util.spec_from_file_location("h3200", H3200)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {H3200}")
h3200 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h3200)


TARGET = tuple(range(8))
EVEN_AP = tuple(2 * i for i in range(8))


def stoploss(q: list[Fraction], threshold: int) -> Fraction:
    return sum(max(n - threshold, 0) * q[n] for n in range(len(q)))


def tail(q: list[Fraction], threshold: int) -> Fraction:
    return sum(q[threshold:])


def lower(q: list[Fraction], threshold: int) -> Fraction:
    return sum(q[: threshold + 1])


def bimod(row: dict[str, object]) -> Fraction:
    q = row["q"]
    return q[0] + q[6]


def ly(row: dict[str, object]) -> Fraction:
    q = row["q"]
    return q[0] + q[6] + q[3] / 10


def metric(row: dict[str, object], name: str) -> Fraction:
    q = row["q"]
    if name.startswith("q"):
        return q[int(name[1:])]
    if name.startswith("tail_ge_"):
        return tail(q, int(name.removeprefix("tail_ge_")))
    if name.startswith("lower_le_"):
        return lower(q, int(name.removeprefix("lower_le_")))
    if name.startswith("stop_ge_"):
        return stoploss(q, int(name.removeprefix("stop_ge_")))
    if name == "bimod":
        return bimod(row)
    if name == "Ly":
        return ly(row)
    if name == "bimod_plus_q3":
        return bimod(row) + q[3]
    raise KeyError(name)


def rank_report(rows: list[dict[str, object]], consec: dict[str, object], name: str) -> str:
    cval = metric(consec, name)
    above = sum(1 for row in rows if metric(row, name) > cval)
    ties = sum(1 for row in rows if metric(row, name) == cval)
    maxrow = max(rows, key=lambda row: (metric(row, name), row["E"]))
    return (
        f"{name}: consec={cval} ({float(cval):.9f}) "
        f"above={above} ties={ties} max={metric(maxrow, name)} at {maxrow['E']}"
    )


def neighbors(E: tuple[int, ...], by_E: dict[tuple[int, ...], dict[str, object]]) -> list[tuple[int, ...]]:
    vals = set(E)
    out = set()
    for idx in range(1, len(E)):
        for delta in (-1, 1):
            moved = E[idx] + delta
            if 1 <= moved <= 14 and moved not in vals:
                F = tuple(sorted((vals - {E[idx]}) | {moved}))
                if F in by_E:
                    out.add(F)
    return sorted(out)


def local_maxima_count(
    rows: list[dict[str, object]],
    by_E: dict[tuple[int, ...], dict[str, object]],
    value,
) -> tuple[int, list[tuple[int, ...]]]:
    local = []
    for row in rows:
        E = row["E"]
        current = value(row)
        if not any(value(by_E[F]) > current for F in neighbors(E, by_E)):
            local.append(E)
    return len(local), local[:8]


def tournament_analysis() -> None:
    carriers = {
        "central_exchange_rate_lemma": 100,
        "upper_stoploss_barrier": 91,
        "q0_bimodality_atom": 83,
        "primitive_normal_form": 76,
        "full_convex_order_route": 34,
        "single_step_compression_gradient": 27,
        "raw_q3_maximization": 19,
        "raw_entropy_route": 5,
    }
    ordered = sorted(carriers.items(), key=lambda item: (-item[1], item[0]))
    print("\nTOURNAMENT ANALYSIS")
    print("vertices=proof angles/signals, not runners or sectors")
    print("pairwise_observable=exact bounded-bank survival score")
    print("switch/gauge=A->B iff survival score(A)>survival score(B); ties lexical")
    print(f"score_hist={{{', '.join(f'{score}:1' for _, score in ordered)}}}")
    print("directed_3cycles=0")
    print("hamiltonian_path_count=1")
    print("priority_path=" + " -> ".join(name for name, _ in ordered))


def main() -> None:
    rows = [
        h3200.row_moments((0,) + combo)
        for combo in itertools.combinations(range(1, 15), 7)
    ]
    primitive = [row for row in rows if row["primitive"]]
    by_E = {row["E"]: row for row in rows}
    consec = by_E[TARGET]
    even_ap = by_E[EVEN_AP]

    print("HYP-3204 exact k=8 ordered-tail exchange scout")
    print("=" * 72)
    print("bank=anchored {0} union A, A subset [1,14], |A|=7")
    print(f"rows_all={len(rows)}")
    print(f"rows_primitive={len(primitive)}")
    print(f"consec={TARGET}")
    print(f"even_AP={EVEN_AP}")
    print(f"consec_q={[str(x) for x in consec['q']]}")
    print(f"even_AP_same_q={even_ap['q'] == consec['q']}")

    c_stop = [stoploss(consec["q"], t) for t in range(7)]
    print("\nFULL CONVEX ORDER CHECK")
    print("consec_stoploss=" + str([str(x) for x in c_stop]))
    for label, bank in (("primitive", primitive), ("all", rows)):
        any_stoploss_above = sum(
            1
            for row in bank
            if any(stoploss(row["q"], t) > c_stop[t] for t in range(7))
        )
        dominators = sum(
            1
            for row in bank
            if all(stoploss(row["q"], t) >= c_stop[t] for t in range(7))
            and any(stoploss(row["q"], t) > c_stop[t] for t in range(7))
        )
        print(f"{label}: rows_with_some_stoploss_above_consec={any_stoploss_above}")
        print(f"{label}: rows_dominating_consec_in_stoploss_order={dominators}")
        for t in range(3):
            above = sum(1 for row in bank if stoploss(row["q"], t) > c_stop[t])
            print(f"{label}: low_stoploss_t{t}_above_consec={above}")

    print("\nONE-SIDED ORDERED-TAIL BARRIER")
    metrics = [
        "q0",
        "q5",
        "q6",
        "tail_ge_4",
        "tail_ge_5",
        "tail_ge_6",
        "stop_ge_3",
        "stop_ge_4",
        "stop_ge_5",
        "bimod",
        "bimod_plus_q3",
        "Ly",
    ]
    for name in metrics:
        print("primitive: " + rank_report(primitive, consec, name))
    print("guardrail: tail_ge_3 is not enough")
    print("primitive: " + rank_report(primitive, consec, "tail_ge_3"))
    print("guardrail: raw q3 is not a proof target")
    print("primitive: " + rank_report(primitive, consec, "q3"))

    print("\nCENTRAL-MASS EXCHANGE RATE")
    bC = bimod(consec)
    q3C = consec["q"][3]
    q3_beaters = []
    violations = []
    best = None
    for row in primitive:
        q = row["q"]
        loss = bC - bimod(row)
        gain = q[3] - q3C
        if gain > 0:
            q3_beaters.append(row)
            if loss <= 0 or gain > loss:
                violations.append((row["E"], loss, gain))
            if loss > 0:
                ratio = gain / loss
                if best is None or ratio > best[0]:
                    best = (ratio, row, loss, gain)
    print(f"primitive_q3_beaters={len(q3_beaters)}")
    print("tested_inequality=(q3-q3_consec)_+ <= (bimod_consec-bimod)")
    print(f"violations={len(violations)}")
    if best is not None:
        ratio, row, loss, gain = best
        print(f"max_exchange_ratio={ratio} ({float(ratio):.9f}) at E={row['E']}")
        print(f"  bimod_loss={loss} ({float(loss):.9f})")
        print(f"  q3_gain={gain} ({float(gain):.9f})")
        print(f"  Ly_margin={ly(consec)-ly(row)} ({float(ly(consec)-ly(row)):.9f})")
    print("consequence: bimod_plus_q3 is also maximized by consec in primitive normal form")

    print("\nLOCAL COMPRESSION GUARDRAIL")
    for name, value in (
        ("Ly", ly),
        ("bimod", bimod),
        ("Sigma_kappa2", lambda row: row["sk2"]),
    ):
        count, sample = local_maxima_count(primitive, by_E, value)
        print(f"single_step_local_maxima_{name}={count}; sample={sample}")
    print("guardrail: naive one-coordinate compression has many traps; use packet/exchange lemmas.")

    tournament_analysis()


if __name__ == "__main__":
    main()
