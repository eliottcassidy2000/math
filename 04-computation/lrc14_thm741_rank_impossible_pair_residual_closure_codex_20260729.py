#!/usr/bin/env python3
r"""Exact pair-residual closure of 38 rank-impossible THM-741 roots.

For a nine-speed root ``E`` in ``{1,...,14}``, put

    G_E = its lonely-time set,       m_E=|G_E|,
    c_E(w)=|G_E intersect D_w|,      w>=15.

The 38 roots below are exactly the non-flood residual roots on which the
scalar two-head/two-tail inequality

    q_1(E)+q_2(E)+2m_E/7 < m_E

fails after the global top coverages are sealed.  Their canonical body digest
is inherited as a hostile membership control.  This verifier closes every one
of them by retaining the exact union of a pair.

Let

    U_2(E)=sup_{15<=a<b}|G_E intersect (D_a union D_b)|.

The globally largest single coverage q_1 is sealed by the exact W=600 atlas.
The strict inequality q_1<4m_E/7 makes the pair tail finite: if

    w > (99/70)r_E / (7(4m_E/7-q_1)),

then q_1+c_E(w)<5m_E/7.  Every pair in 15..2500 is branch-and-bound
maximized exactly; each paid carrier is cross-checked by direct eleven-comb
reconstruction.  Combining the head maximum with q_1+u_E(2501), where

    u_E(w)=m_E/7+(99/70)r_E/(7w),

gives a global rational cap Ubar_2(E)<5m_E/7.

Choose W_E so Ubar_2(E)+2u_E(W_E+1)<m_E.  If a quadruple has at least two
speeds above W_E, take the other two speeds as one pair and union-bound only
the two tails: it is closed.  Thus a possible obstruction has three smallest
speeds in 15..W_E.

Exactly one tail is handled output-sensitively.  Order the three head speeds
by descending c_E, compute the exact union of the top pair, and construct the
three-comb residual only if the pair union plus the third coverage and the
root tail cap can reach m_E.  THM-732's one-speed estimate on that literal
residual closes every t>W_E.

If all four speeds lie in the head, order them by descending c_E.  The top
pair is exact; the top-three residual is constructed only when its fourth
single coverage can reach it; every remaining quadruple is rebuilt both by
nested subtraction and direct thirteen-comb union.  All survivors are
strictly positive.

If an added speed is at most 14, the completed family has ten speeds in
{1,...,14} and is already lonely by THM-738.  Hence the pure-tail computation
closes all 38 whole THM-741 roots.  It does not prove global THM-741 or
LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
FIRST_EXTERNAL = 15
RANK_HORIZON = 600
PAIR_HORIZON = 2500
BODIES = (
    (1, 2, 3, 4, 5, 7, 8, 11, 13),
    (1, 2, 3, 4, 5, 7, 9, 11, 13),
    (1, 2, 3, 4, 7, 8, 10, 11, 13),
    (1, 2, 3, 4, 7, 9, 10, 11, 13),
    (1, 2, 3, 4, 8, 10, 11, 13, 14),
    (1, 2, 3, 5, 7, 8, 9, 12, 13),
    (1, 2, 3, 5, 7, 8, 10, 11, 13),
    (1, 2, 3, 5, 7, 8, 11, 12, 13),
    (1, 2, 3, 5, 7, 9, 10, 11, 13),
    (1, 2, 3, 6, 7, 8, 9, 10, 13),
    (1, 2, 3, 6, 7, 9, 10, 11, 13),
    (1, 2, 3, 7, 8, 9, 10, 12, 13),
    (1, 2, 3, 7, 8, 10, 11, 12, 13),
    (1, 2, 3, 7, 8, 10, 11, 13, 14),
    (1, 2, 4, 5, 6, 7, 8, 9, 13),
    (1, 2, 4, 6, 7, 8, 9, 10, 13),
    (1, 2, 4, 6, 8, 9, 10, 13, 14),
    (1, 2, 4, 6, 9, 10, 12, 13, 14),
    (1, 2, 4, 8, 9, 10, 12, 13, 14),
    (1, 2, 5, 6, 7, 8, 9, 10, 13),
    (1, 2, 6, 7, 8, 9, 10, 12, 13),
    (1, 2, 6, 8, 9, 10, 12, 13, 14),
    (1, 3, 4, 5, 7, 8, 10, 11, 13),
    (1, 3, 4, 5, 8, 10, 11, 13, 14),
    (1, 3, 4, 6, 7, 8, 9, 10, 13),
    (1, 3, 4, 7, 8, 9, 10, 12, 13),
    (1, 3, 4, 7, 8, 10, 11, 13, 14),
    (1, 3, 4, 8, 9, 10, 12, 13, 14),
    (1, 4, 5, 6, 7, 8, 9, 10, 13),
    (1, 4, 5, 6, 7, 8, 9, 11, 13),
    (1, 4, 5, 8, 9, 10, 12, 13, 14),
    (1, 4, 6, 7, 8, 9, 10, 11, 13),
    (1, 4, 6, 7, 8, 9, 10, 13, 14),
    (1, 4, 6, 8, 9, 10, 12, 13, 14),
    (1, 5, 6, 7, 8, 9, 10, 11, 13),
    (2, 3, 4, 5, 7, 8, 10, 11, 13),
    (2, 3, 4, 7, 8, 10, 11, 13, 14),
    (2, 4, 6, 8, 9, 10, 12, 13, 14),
)
EXPECTED_BODY_DIGEST = "81e0807f389f0b10ac0ff4bfd9e3c84f86b003edc9b5469d1c581ca6e25b0c4d"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(sha256(CORE_PATH) == CORE_SHA256, "THM-741 core hash changed")
    spec = importlib.util.spec_from_file_location("thm741_rank2_pair_core", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CORE = load_core()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def coverage(good: list[tuple[F, F]], root_m: F, speed: int) -> F:
    return root_m - CORE.subtract_sparse(good, speed)


def direct_carrier(body: tuple[int, ...], extras: tuple[int, ...]):
    family = tuple(sorted((*body, *extras)))
    require(len(family) == len(set(family)), f"duplicate family {family}")
    return family, CORE.good_norm(family)


def probe_body(body: tuple[int, ...]) -> dict[str, object]:
    good, root_r, root_m = CORE.good_norm(body)
    require(root_r == len(good) and root_m > 0, f"bad root carrier {body}")

    covers = {
        speed: coverage(good, root_m, speed)
        for speed in range(FIRST_EXTERNAL, RANK_HORIZON + 1)
    }
    ranked600 = sorted(covers, key=lambda speed: (-covers[speed], speed))
    q1 = covers[ranked600[0]]
    q2 = covers[ranked600[1]]
    q3 = covers[ranked600[2]]
    q4 = covers[ranked600[3]]
    require(q4 > root_m / 7, f"fourth rank below limiting cap {body}")
    rank_threshold = CORE.S2 * root_r / (7 * (q4 - root_m / 7))
    require(rank_threshold < RANK_HORIZON + 1, f"global ranks unsealed {body}")
    require(
        root_m - q1 - q2 - q3 <= root_m / 7,
        f"rank-impossible type changed {body}",
    )
    require(q1 + q2 >= F(5, 7) * root_m, f"scalar A2 type changed {body}")

    q1_margin = F(4, 7) * root_m - q1
    require(q1_margin > 0, f"q1 reaches 4m/7 {body}")
    pair_entry_threshold = CORE.S2 * root_r / (7 * q1_margin)
    require(pair_entry_threshold < PAIR_HORIZON + 1, f"pair horizon too short {body}")

    for speed in range(RANK_HORIZON + 1, PAIR_HORIZON + 1):
        covers[speed] = coverage(good, root_m, speed)
    ranked_pair = sorted(covers, key=lambda speed: (-covers[speed], speed))

    head_pair_cap = F(0)
    exact_pair_checks = 0
    pair_rows = []
    for i, a in enumerate(ranked_pair[:-1]):
        if covers[a] + covers[ranked_pair[i + 1]] <= head_pair_cap:
            break
        _, _, after_a = CORE.subtract(good, a)
        for b in ranked_pair[i + 1 :]:
            if covers[a] + covers[b] <= head_pair_cap:
                break
            pair_r, pair_m, pair_good = CORE.subtract(after_a, b)
            family, (direct_good, direct_r, direct_m) = direct_carrier(body, (a, b))
            require(
                direct_good == pair_good
                and direct_r == pair_r
                and direct_m == pair_m,
                f"pair paths disagree {family}",
            )
            pair_union = root_m - pair_m
            head_pair_cap = max(head_pair_cap, pair_union)
            exact_pair_checks += 1
            pair_rows.append(
                f"E={','.join(map(str,body))};P={a},{b};"
                f"U={ftext(pair_union)};h={ftext(pair_m)};r={pair_r}\n"
            )

    root_pair_tail = (
        q1
        + root_m / 7
        + CORE.S2 * root_r / (7 * (PAIR_HORIZON + 1))
    )
    pair_cap = max(head_pair_cap, root_pair_tail)
    pair_margin = F(5, 7) * root_m - pair_cap
    require(pair_margin > 0, f"global pair cap reaches 5m/7 {body}")

    W = max(14, floor_fraction(2 * CORE.S2 * root_r / (7 * pair_margin)))
    root_tail_cap = root_m / 7 + CORE.S2 * root_r / (7 * (W + 1))
    multi_tail_margin = root_m - pair_cap - 2 * root_tail_cap
    require(multi_tail_margin > 0, f"two-tail partition failed {body}")

    for speed in range(PAIR_HORIZON + 1, W + 1):
        covers[speed] = coverage(good, root_m, speed)
    ranked = sorted(
        (speed for speed in covers if speed <= W),
        key=lambda speed: (-covers[speed], speed),
    )

    scalar_triples = 0
    one_tail_pairs = 0
    one_tail_triples = 0
    tail_minimum = None
    tail_rows = []
    target = root_m - root_tail_cap

    for i, a in enumerate(ranked[:-2]):
        if covers[a] + covers[ranked[i + 1]] + covers[ranked[i + 2]] < target:
            break
        _, _, after_a = CORE.subtract(good, a)
        for j in range(i + 1, len(ranked) - 1):
            b = ranked[j]
            if covers[a] + covers[b] + covers[ranked[j + 1]] < target:
                break
            k0 = j + 1
            while k0 < len(ranked) and covers[a] + covers[b] + covers[ranked[k0]] >= target:
                scalar_triples += 1
                k0 += 1

            pair_r, pair_m, pair_good = CORE.subtract(after_a, b)
            one_tail_pairs += 1
            pair_union = root_m - pair_m
            if pair_union + covers[ranked[j + 1]] < target:
                continue
            for c in ranked[j + 1 :]:
                if pair_union + covers[c] < target:
                    break
                triple_r, triple_m, triple_good = CORE.subtract(pair_good, c)
                family, (direct_good, direct_r, direct_m) = direct_carrier(body, (a, b, c))
                require(
                    direct_good == triple_good
                    and direct_r == triple_r
                    and direct_m == triple_m,
                    f"triple paths disagree {family}",
                )
                tail_survivor = (
                    F(6, 7) * triple_m
                    - CORE.S2 * triple_r / (7 * (W + 1))
                )
                require(tail_survivor > 0, f"one-tail bound failed {family}")
                one_tail_triples += 1
                record = (tail_survivor, family, triple_m, triple_r)
                if tail_minimum is None or record < tail_minimum:
                    tail_minimum = record
                tail_rows.append(
                    f"E={','.join(map(str,body))};T={a},{b},{c};"
                    f"h={ftext(triple_m)};r={triple_r};"
                    f"tail={ftext(tail_survivor)}\n"
                )

    all_head_pairs = 0
    all_head_triples = 0
    dangerous_quadruples = 0
    head_minimum = None
    head_rows = []
    for i, a in enumerate(ranked[:-3]):
        if (
            covers[a]
            + covers[ranked[i + 1]]
            + covers[ranked[i + 2]]
            + covers[ranked[i + 3]]
            < root_m
        ):
            break
        _, _, after_a = CORE.subtract(good, a)
        for j in range(i + 1, len(ranked) - 2):
            b = ranked[j]
            if (
                covers[a]
                + covers[b]
                + covers[ranked[j + 1]]
                + covers[ranked[j + 2]]
                < root_m
            ):
                break
            pair_r, pair_m, pair_good = CORE.subtract(after_a, b)
            all_head_pairs += 1
            pair_union = root_m - pair_m
            if (
                pair_union
                + covers[ranked[j + 1]]
                + covers[ranked[j + 2]]
                < root_m
            ):
                continue
            for k in range(j + 1, len(ranked) - 1):
                c = ranked[k]
                if pair_union + covers[c] + covers[ranked[k + 1]] < root_m:
                    break
                triple_r, triple_m, triple_good = CORE.subtract(pair_good, c)
                all_head_triples += 1
                if covers[ranked[k + 1]] < triple_m:
                    continue
                for d in ranked[k + 1 :]:
                    if covers[d] < triple_m:
                        break
                    nested_r, survivor, nested_good = CORE.subtract(triple_good, d)
                    family, (direct_good, direct_r, direct_m) = direct_carrier(
                        body, (a, b, c, d)
                    )
                    require(
                        len(family) == 13
                        and direct_good == nested_good
                        and direct_r == nested_r
                        and direct_m == survivor,
                        f"quadruple paths disagree {family}",
                    )
                    require(survivor > 0, f"nonpositive survivor {family}")
                    dangerous_quadruples += 1
                    record = (survivor, family, direct_r)
                    if head_minimum is None or record < head_minimum:
                        head_minimum = record
                    head_rows.append(
                        f"E={','.join(map(str,body))};"
                        f"Q={a},{b},{c},{d};L={ftext(survivor)};r={direct_r}\n"
                    )

    return {
        "body": body,
        "r": root_r,
        "m": root_m,
        "q1_margin": q1_margin,
        "pair_entry_threshold": pair_entry_threshold,
        "head_pair_cap": head_pair_cap,
        "root_pair_tail": root_pair_tail,
        "pair_cap": pair_cap,
        "pair_margin": pair_margin,
        "W": W,
        "multi_tail_margin": multi_tail_margin,
        "exact_pair_checks": exact_pair_checks,
        "scalar_triples": scalar_triples,
        "one_tail_pairs": one_tail_pairs,
        "one_tail_triples": one_tail_triples,
        "tail_minimum": tail_minimum,
        "all_head_pairs": all_head_pairs,
        "all_head_triples": all_head_triples,
        "dangerous_quadruples": dangerous_quadruples,
        "head_minimum": head_minimum,
        "pair_rows": tuple(pair_rows),
        "tail_rows": tuple(tail_rows),
        "head_rows": tuple(head_rows),
    }


def minimum_row(rows: list[dict[str, object]], key: str):
    candidates = [row for row in rows if row[key] is not None]
    require(candidates, f"no hostile controls for {key}")
    return min(candidates, key=lambda row: (row[key], row["body"]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(6, os.cpu_count() or 1))
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    require(CORE.S2 == F(99, 70) and CORE.S2**2 > 2, "sqrt(2) majorant changed")
    body_text = "\n".join(",".join(map(str, body)) for body in BODIES)
    body_digest = hashlib.sha256(body_text.encode()).hexdigest()
    require(len(BODIES) == 38 and body_digest == EXPECTED_BODY_DIGEST, "body atlas changed")

    if args.workers == 1:
        rows = list(map(probe_body, BODIES))
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = pool.map(probe_body, BODIES, chunksize=1)
    require(tuple(row["body"] for row in rows) == BODIES, "parallel order changed")

    min_q1 = min(rows, key=lambda row: (row["q1_margin"], row["body"]))
    max_entry = max(rows, key=lambda row: (row["pair_entry_threshold"], row["body"]))
    min_pair = min(rows, key=lambda row: (row["pair_margin"], row["body"]))
    max_W = max(rows, key=lambda row: (row["W"], row["body"]))
    min_multi = min(rows, key=lambda row: (row["multi_tail_margin"], row["body"]))
    tail_controls = [row for row in rows if row["tail_minimum"] is not None]
    head_controls = [row for row in rows if row["head_minimum"] is not None]
    require(tail_controls and head_controls, "hostile controls disappeared")
    min_tail = min(tail_controls, key=lambda row: (row["tail_minimum"], row["body"]))
    min_head = min(head_controls, key=lambda row: (row["head_minimum"], row["body"]))

    pair_ledger = "THM741/rank2-pair/v1\n" + "".join(
        line for row in rows for line in row["pair_rows"]
    )
    tail_ledger = "THM741/rank2-one-tail/v1\n" + "".join(
        line for row in rows for line in row["tail_rows"]
    )
    head_ledger = "THM741/rank2-all-head/v1\n" + "".join(
        line for row in rows for line in row["head_rows"]
    )
    consequence_ledger = "THM741/rank2-consequence/v1\n" + "".join(
        "E="
        + ",".join(map(str, row["body"]))
        + f";U2={ftext(row['pair_cap'])};W={row['W']};"
        + f"tailN={row['one_tail_triples']};headN={row['dangerous_quadruples']}\n"
        for row in rows
    )

    print("THM-741 RANK-IMPOSSIBLE PAIR-RESIDUAL CLOSURE")
    print("status=FINITE-EXACT+GLOBAL-TAIL-SEALED")
    print("root_universe=38 scalar-A2-failing non-flood rank-impossible roots")
    print("body_digest_sha256=" + body_digest)
    print(
        "minimum_q1_margin_below_4m7="
        + ftext(min_q1["q1_margin"])
        + ";body="
        + ",".join(map(str, min_q1["body"]))
    )
    print(
        "maximum_pair_entry_threshold="
        + ftext(max_entry["pair_entry_threshold"])
        + ";body="
        + ",".join(map(str, max_entry["body"]))
        + f";horizon={PAIR_HORIZON}"
    )
    print(
        f"exact_pair_head_checks={sum(row['exact_pair_checks'] for row in rows)};"
        "pair_tail_cap=q1+u(2501)"
    )
    print(
        "minimum_global_pair_margin_below_5m7="
        + ftext(min_pair["pair_margin"])
        + ";body="
        + ",".join(map(str, min_pair["body"]))
    )
    print(
        f"maximum_root_cutoff={max_W['W']};body="
        + ",".join(map(str, max_W["body"]))
    )
    print(
        "minimum_multi_tail_margin="
        + ftext(min_multi["multi_tail_margin"])
        + ";body="
        + ",".join(map(str, min_multi["body"]))
        + f";W={min_multi['W']}"
    )
    print(
        f"one_tail_scalar_triples={sum(row['scalar_triples'] for row in rows)};"
        f"exact_pair_carriers={sum(row['one_tail_pairs'] for row in rows)};"
        f"exact_triple_carriers={sum(row['one_tail_triples'] for row in rows)};"
        "failures=0"
    )
    print(
        "minimum_one_tail_survivor="
        + ftext(min_tail["tail_minimum"][0])
        + ";family="
        + ",".join(map(str, min_tail["tail_minimum"][1]))
        + ";h="
        + ftext(min_tail["tail_minimum"][2])
        + f";r={min_tail['tail_minimum'][3]}"
    )
    print(
        f"all_head_pair_carriers={sum(row['all_head_pairs'] for row in rows)};"
        f"triple_carriers={sum(row['all_head_triples'] for row in rows)};"
        f"dangerous_quadruples={sum(row['dangerous_quadruples'] for row in rows)};"
        "positive=all;tight=0"
    )
    print(
        "minimum_head_survivor="
        + ftext(min_head["head_minimum"][0])
        + ";family="
        + ",".join(map(str, min_head["head_minimum"][1]))
        + f";r={min_head['head_minimum'][2]}"
    )
    print("pair_ledger_sha256=" + hashlib.sha256(pair_ledger.encode()).hexdigest())
    print("one_tail_ledger_sha256=" + hashlib.sha256(tail_ledger.encode()).hexdigest())
    print("all_head_ledger_sha256=" + hashlib.sha256(head_ledger.encode()).hexdigest())
    print(
        "consequence_ledger_sha256="
        + hashlib.sha256(consequence_ledger.encode()).hexdigest()
    )
    print("whole_root_composition=38 pure tails + THM-738 small-speed chamber")
    print("scope=38 new whole THM-741 roots; global THM-741 and LRC14 remain open")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
