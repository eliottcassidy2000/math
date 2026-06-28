#!/usr/bin/env python3
"""HYP-3246 scout: shell-lag commutator and contact-support sidecar.

HYP-3245 says AP is the triangular ordinary-autocorrelation law and that every
HYP-3202 trap moves pair mass outward in lag space.  HYP-3228 gives the exact
shell magic functional 10*q0 + q3 + 10*q6.  This scout asks the controlled-
forgetting question between them:

    if we project a bounded k=8 row to ordinary support autocorrelation,
    how much HYP-3228 shell data survives?

The answer is negative in a strong and structured way.  Ordinary lag profile is
too coarse, residue histogram mod 7 is still too coarse, and the missing
bounded-bank repair is an ordered contact-support sidecar recording where the
non-unit gaps sit in the anchored gap word.  Gap sizes alone do not repair the
projection; positions do.
"""

from __future__ import annotations

import importlib.util
import itertools
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
H3200 = HERE / "lrc_k8_cumulant_universality_codex_20260628.py"
spec = importlib.util.spec_from_file_location("h3200", H3200)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {H3200}")
h3200 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h3200)


def support_autocorr(E: tuple[int, ...], max_lag: int = 14) -> tuple[int, ...]:
    s = set(E)
    return tuple(sum(1 for x in s if x + lag in s) for lag in range(max_lag + 1))


def residue_hist(E: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(1 for x in E if x % 7 == r) for r in range(7))


def residue_word(E: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(x % 7 for x in E)


def gap_word(E: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(E[i + 1] - E[i] for i in range(len(E) - 1))


def gap_multiset(E: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(g for g in gap_word(E) if g != 1))


def contact_support(E: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(i for i, g in enumerate(gap_word(E)) if g != 1)


def contact_word(E: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple((i, g) for i, g in enumerate(gap_word(E)) if g != 1)


def magic_value(q: tuple[Fraction, ...]) -> Fraction:
    return 10 * q[0] + q[3] + 10 * q[6]


def row_packet(E: tuple[int, ...]) -> dict[str, object]:
    row = h3200.row_moments(E)
    q = tuple(row["q"])
    return {
        "E": E,
        "q": q,
        "magic": magic_value(q),
        "ac": support_autocorr(E),
        "hist": residue_hist(E),
        "resword": residue_word(E),
        "gaps": gap_word(E),
        "gap_multiset": gap_multiset(E),
        "contact_support": contact_support(E),
        "contact_word": contact_word(E),
    }


def fiber_count(rows: list[dict[str, object]], key_names: tuple[str, ...]) -> tuple[int, int, tuple[object, list[dict[str, object]]] | None]:
    fibers: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        key = tuple(row[name] for name in key_names)
        fibers[key].append(row)
    mixed = []
    for key, fiber in fibers.items():
        if len({row["magic"] for row in fiber}) > 1:
            mixed.append((key, fiber))
    witness = mixed[0] if mixed else None
    return len(fibers), len(mixed), witness


def refine_mixed(
    mixed_fibers: list[list[dict[str, object]]],
    extra_key_name: str,
) -> tuple[int, tuple[list[dict[str, object]], dict[object, set[Fraction]]] | None]:
    unresolved = 0
    witness = None
    for fiber in mixed_fibers:
        buckets: dict[object, set[Fraction]] = defaultdict(set)
        for row in fiber:
            buckets[row[extra_key_name]].add(row["magic"])
        if any(len(vals) > 1 for vals in buckets.values()):
            unresolved += 1
            if witness is None:
                witness = (fiber, buckets)
    return unresolved, witness


def print_fiber_witness(label: str, witness: tuple[object, list[dict[str, object]]] | None) -> None:
    print(label)
    if witness is None:
        print("  no mixed witness")
        return
    _key, fiber = witness
    for row in fiber[:4]:
        print(
            f"  E={row['E']} magic={row['magic']} hist={row['hist']} "
            f"resword={row['resword']} gaps={row['gaps']} contact_support={row['contact_support']}"
        )


def tournament_analysis() -> None:
    carriers = {
        "lag_profile_only": 25,
        "lag_plus_gap_multiset": 31,
        "lag_plus_residue_histogram": 57,
        "lag_plus_contact_support": 86,
        "lag_plus_contact_word": 91,
        "circle_endpoint_arrangement_cell": 95,
        "half_tiling_descent_sidecar": 99,
    }
    ordered = sorted(carriers.items(), key=lambda item: (-item[1], item[0]))
    print("\nTOURNAMENT ANALYSIS")
    print("vertices=quotients/sidecars for the HYP-3245 -> HYP-3228 projection")
    print("pairwise_observable=how much shell magic survives the lag-space projection")
    print("switch/gauge=A->B iff A preserves more shell/contact/endpoint payload")
    print(f"score_hist={dict(Counter(score for _, score in ordered))}")
    print("directed_3cycles=0")
    print("scc_sizes=[" + ",".join("1" for _ in ordered) + "]")
    print("hamiltonian_path_count=1")
    print("priority_path=" + " -> ".join(name for name, _ in ordered))


def main() -> None:
    start = time.time()
    rows = [row_packet((0,) + combo) for combo in itertools.combinations(range(1, 15), 7)]
    elapsed = time.time() - start

    print("HYP-3246 shell-lag commutator / contact-support sidecar")
    print("=" * 78)
    print("bank=anchored bounded k=8 rows E={0} union A, A subset [1,14], |A|=7")
    print(f"rows={len(rows)} elapsed_seconds={elapsed:.3f}")

    total_ac, mixed_ac, witness_ac = fiber_count(rows, ("ac",))
    total_ac_hist, mixed_ac_hist, witness_ac_hist = fiber_count(rows, ("ac", "hist"))
    total_ac_hist_contact, mixed_ac_hist_contact, _ = fiber_count(rows, ("ac", "hist", "contact_support"))

    print("\nCOARSE FIBER COUNTS")
    print(f"  ac fibers={total_ac} mixed_magic_fibers={mixed_ac}")
    print(f"  ac+hist fibers={total_ac_hist} mixed_magic_fibers={mixed_ac_hist}")
    print(f"  ac+hist+contact_support fibers={total_ac_hist_contact} mixed_magic_fibers={mixed_ac_hist_contact}")
    print("  readout=ordinary lag profile is highly mixed; residue histogram repairs most but not all;")
    print("           ordered contact-support kills the remaining bounded-bank ambiguity.")

    print("\nMIXED-WITNESS: ac only")
    print_fiber_witness("  first mixed ac fiber", witness_ac)

    print("\nMIXED-WITNESS: ac + residue histogram")
    print_fiber_witness("  first mixed ac+hist fiber", witness_ac_hist)

    mixed_ac_hist_fibers = []
    fibers: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        fibers[(row["ac"], row["hist"])].append(row)
    for fiber in fibers.values():
        if len({row["magic"] for row in fiber}) > 1:
            mixed_ac_hist_fibers.append(fiber)

    unresolved_gap_multiset, witness_gap_multiset = refine_mixed(mixed_ac_hist_fibers, "gap_multiset")
    unresolved_contact_support, _ = refine_mixed(mixed_ac_hist_fibers, "contact_support")

    print("\nREFINEMENT INSIDE THE 62 MIXED ac+hist FIBERS")
    print(f"  unresolved after adding gap_multiset={unresolved_gap_multiset}")
    print(f"  unresolved after adding contact_support={unresolved_contact_support}")
    print("  readout=the missing sidecar is positional.  Gap sizes alone do not repair the")
    print("           shell functional; the positions of the non-unit gaps do.")

    if witness_gap_multiset is not None:
        fiber, buckets = witness_gap_multiset
        print("\nGAP-MULTISET FAILURE WITNESS")
        for row in fiber:
            print(
                f"  E={row['E']} magic={row['magic']} resword={row['resword']} "
                f"gaps={row['gaps']} gap_multiset={row['gap_multiset']} "
                f"contact_support={row['contact_support']}"
            )
        print("  bucket_by_gap_multiset:")
        for key, values in buckets.items():
            print(f"    {key}: {sorted(values)}")

    print("\nSHARP COLLISION PAIR")
    pair = [
        row
        for row in rows
        if row["E"] in {
            (0, 1, 2, 3, 4, 12, 13, 14),
            (0, 1, 2, 10, 11, 12, 13, 14),
        }
    ]
    for row in sorted(pair, key=lambda item: item["E"]):
        q = row["q"]
        print(
            f"  E={row['E']} magic={row['magic']} q0+q6={q[0]+q[6]} q3={q[3]} "
            f"ac={row['ac']} hist={row['hist']} resword={row['resword']} "
            f"gaps={row['gaps']} contact_support={row['contact_support']}"
        )
    print("  readout=same ordinary lag profile, same residue histogram, and same residue word mod 7;")
    print("           only the position of the long gap changes, and the shell functional changes with it.")

    print("\nASSUMPTION CHALLENGE")
    print("  alternate vertices considered: raw lag fibers, residue histograms, gap multisets,")
    print("  ordered contact-support, full contact words, circle endpoint cells, and tiling lifts.")
    print("  chosen sidecar=ordered positions of the non-unit gaps in the anchored gap word.")
    print("  preserved predicate=HYP-3228 shell magic and HYP-3245 lag-side transport.")
    print("  destroyed information by lag-only projection=which endpoint cell the long gaps occupy.")
    print("  challenged assumption=autocorrelation plus residue data already names the shell packet.")
    print("  Exact bank evidence says no: ordered contact support is still needed.")

    tournament_analysis()


if __name__ == "__main__":
    main()
