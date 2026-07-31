#!/usr/bin/env python3
"""Saturate the K25 proof path under THM-2897 q5+2B2 activation.

Every inactive apex at every actual scalar-closed state is tested.  A cheap
exact lower witness U({x,y})<=B2 refutes the strict activation whenever
q5+2U({x,y})>=h.  Only candidates surviving that hostile test pay for a
globally exact B2 cap.  At a joint scalar/rank-pair fixed point, the next
apex from the known scalar seed bank is reserved for parity.
"""

from __future__ import annotations

import importlib.util
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (
    ROOT
    / ".scratch/lrc_j6_full_pipeline_design_20260729/"
    "lrc14_j6_k25_closed_state_parity_closure_codex_20260729.py"
)
SEED_BANK = (23, 27, 19, 46, 18, 17)


def load(path: Path):
    spec = importlib.util.spec_from_file_location("k25_rankpair_source", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load(SOURCE)


def main() -> None:
    base = S.A.profile_body(S.BODY)
    good, _, _ = S.A.T.CORE.good_norm(S.BODY)
    root = {**base, "good": good, "K": base["adaptive_k"]}
    gate = tuple(speed for _, speed in root["top"][: root["K"]])
    profiles = {
        apex: S.M.profile_apex(root, gate, apex) for apex in gate
    }
    index = {speed: bit for bit, speed in enumerate(gate)}
    carriers = {
        apex: S.M.M.T.CORE.good_norm(
            tuple(sorted((*S.BODY, apex)))
        )[0]
        for apex in gate
    }

    mask, rounds, _ = S.scalar_closure(gate, profiles, 0)
    if rounds:
        raise RuntimeError("empty prefix unexpectedly scalar-activates")
    prefix: list[int] = []
    proof_steps: list[tuple[object, ...]] = []
    audit_rows: list[tuple[object, ...]] = []
    paid_index = 0
    full = (1 << len(gate)) - 1
    while mask != full:
        passing: list[int] = []
        local: list[tuple[object, ...]] = []
        for apex in gate:
            if mask & (1 << index[apex]):
                continue
            excluded = set(prefix) | {apex}
            top5 = S.available_top(profiles[apex], excluded, 5)
            q5 = top5[4][0]
            first = top5[0][1]
            second = top5[1][1]
            residual = S.Q.H.R.subtract_local_multi(
                carriers[apex],
                (first, second),
            )
            pair_union = profiles[apex]["m"] - S.interval_mass(residual)
            hostile_margin = profiles[apex]["m"] - q5 - 2 * pair_union
            if hostile_margin <= 0:
                local.append(
                    (
                        apex,
                        "LOWER_REFUTED",
                        first,
                        second,
                        S.ftext(q5),
                        S.ftext(pair_union),
                        S.ftext(hostile_margin),
                    )
                )
                continue
            exact = S.Q.H.pair_cap(carriers[apex], excluded)
            exact_margin = profiles[apex]["m"] - q5 - 2 * exact["cap"]
            local.append(
                (
                    apex,
                    "EXACT_PASS" if exact_margin > 0 else "EXACT_FAIL",
                    first,
                    second,
                    S.ftext(q5),
                    S.ftext(pair_union),
                    S.ftext(hostile_margin),
                    S.ftext(exact["cap"]),
                    S.ftext(exact_margin),
                    exact["paid"],
                )
            )
            if exact_margin > 0:
                passing.append(apex)
        audit_rows.extend(
            (tuple(prefix), *row) for row in local
        )
        print(
            f"P={tuple(prefix)};inactive={len(local)};"
            f"lower_refuted={sum(row[1]=='LOWER_REFUTED' for row in local)};"
            f"exact={sum(row[1]!='LOWER_REFUTED' for row in local)};"
            f"passing={tuple(passing)}"
        )
        if passing:
            for apex in passing:
                prefix.append(apex)
                mask |= 1 << index[apex]
            target, scalar_rounds, _ = S.scalar_closure(gate, profiles, mask)
            for row in scalar_rounds:
                prefix.extend(row)
            proof_steps.append(
                ("RANKPAIR", tuple(passing), scalar_rounds)
            )
            mask = target
            continue

        while (
            paid_index < len(SEED_BANK)
            and mask & (1 << index[SEED_BANK[paid_index]])
        ):
            paid_index += 1
        if paid_index == len(SEED_BANK):
            raise RuntimeError("seed bank exhausted before full closure")
        apex = SEED_BANK[paid_index]
        paid_index += 1
        prefix.append(apex)
        mask |= 1 << index[apex]
        target, scalar_rounds, _ = S.scalar_closure(gate, profiles, mask)
        for row in scalar_rounds:
            prefix.extend(row)
        proof_steps.append(("PARITY", apex, scalar_rounds))
        mask = target

    print(f"proof_steps={tuple(proof_steps)}")
    print(f"parity_count={sum(row[0]=='PARITY' for row in proof_steps)}")
    print(f"rankpair_count={sum(row[0]=='RANKPAIR' for row in proof_steps)}")
    print(f"audit_rows={len(audit_rows)}")
    print(f"ordered_prefix={tuple(prefix)}")


if __name__ == "__main__":
    main()
