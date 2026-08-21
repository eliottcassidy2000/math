#!/usr/bin/env python3
"""Exact dyadic Rule-A transition and adjoint atlas for THM-3597."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT / "04-computation/amm12592_r512_truncated_adjoint_pascal_horizons_thm3588.py"
EXPECTED_PARENT_SHA256 = "692ac6411f9a99c2734fe85995915b144a9752528816bae6df73cd0943c53985"
EXPECTED = {
    256: {"failed": {0: (61, 25)}, "survivor": 1},
    512: {
        "failed": {0: (107, 36), 1: (110, 39), 2: (113, 40),
                   3: (116, 43), 4: (121, 45)},
        "survivor": 5,
    },
    1024: {
        "failed": {
            0: (207, 64), 1: (209, 66), 2: (211, 67),
            3: (213, 69), 4: (216, 70), 5: (219, 72),
            6: (221, 74), 7: (223, 76), 8: (227, 78),
            9: (229, 81), 10: (233, 82), 11: (236, 84),
            12: (240, 88), 13: (244, 90), 14: (250, 95),
        },
        "survivor": 15,
    },
}
EXPECTED_SEMANTIC_SHA256 = "4169c97f724ae4a510dbc2b258f2f0039387d7758a552ee783d30d6b24b29fbe"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


require(lf_sha256(PARENT) == EXPECTED_PARENT_SHA256, "THM-3588 parent hash")
source = PARENT.read_text(encoding="utf-8")
parts = source.split("\nR = 512\n", 1)
require(len(parts) == 2, "THM-3588 definition split")
namespace = {"__name__": "thm3597_parent_definitions"}
exec(compile(parts[0], str(PARENT), "exec"), namespace)

comb = namespace["comb"]
floor_gamma_star = namespace["floor_gamma_star"]
two_g_coeffs = namespace["two_g_coeffs"]
t4_row_load = namespace["t4_row_load"]
kernel = namespace["kernel"]
baseline_before_death = namespace["baseline_before_death"]
certificate = namespace["certificate"]


def clamp(load: list[int], degree: int):
    state = []
    for t, value in enumerate(load):
        lower = -comb(degree - 1, t) if t < degree else 0
        upper = comb(degree - 1, t - 1) if t else 0
        correction = min(upper, max(lower, value))
        state.append(value - correction)
    return tuple(state[:degree]), state[degree]


def rule_a_first_death(r: int, offset: int):
    profile = tuple(floor_gamma_star(r + i) + offset for i in range(r))
    require(all(profile[i] - profile[i - 1] in (0, 1) for i in range(1, r)),
            (r, offset, "degree word"))
    initial_full = t4_row_load(r, profile[0])
    require(all(value % 2 == 0 for value in initial_full), (r, offset, "initial parity"))
    state, top = clamp([value // 2 for value in initial_full], profile[0])
    if top:
        return profile, 0
    rows = [state]
    feed = [value // 2 for value in two_g_coeffs(r)]
    for i in range(1, r):
        degree, previous = profile[i], profile[i - 1]
        transition = kernel(degree, previous)
        load = [0] * (degree + 1)
        for t, value in enumerate(rows[-1]):
            for q, coefficient in enumerate(transition):
                load[t + q] += coefficient * value
        if degree + i <= r - 1:
            load[0] += feed[degree + i]
        if degree > previous and degree - 1 + i <= r - 1:
            value = feed[degree - 1 + i]
            load[0] += value
            load[1] += value
        state, top = clamp(load, degree)
        if top:
            return profile, i
        rows.append(state)
    return profile, None


def certificate_sweep(profile, rows, loads, death):
    """Compute every cut value in one exact backward-adjoint sweep."""
    d = profile[death]
    transition = kernel(d, profile[death - 1])
    level = {}
    for q, value in enumerate(transition):
        t = d - q
        if 0 <= t < profile[death - 1]:
            level[t] = level.get(t, 0) + value

    rhs = loads[death][d]
    distinct_rows = 0
    multiplier_mass = 0
    max_multiplier = 0
    ledger = {}
    for i in range(death - 1, 0, -1):
        d_i = profile[i]
        distinct_rows += len(level)
        multiplier_mass += sum(level.values())
        max_multiplier = max(max_multiplier, *level.values())
        for t, multiplier in level.items():
            base_load = loads[i][t]
            base_state = rows[i][t]
            base_correction = base_load - base_state
            lower = -comb(d_i - 1, t)
            rhs += multiplier * (base_correction - lower)
        ledger[i] = (rhs, distinct_rows, multiplier_mass, max_multiplier)

        transition = kernel(d_i, profile[i - 1])
        previous = {}
        for t, multiplier in level.items():
            for q, value in enumerate(transition):
                u = t - q
                if 0 <= u < profile[i - 1]:
                    previous[u] = previous.get(u, 0) + multiplier * value
        level = previous

    require(tuple(ledger) == tuple(range(death - 1, 0, -1)),
            "complete adjoint sweep")
    return ledger


def boundary_certificate(r: int, offset: int, death: int, first_negative: int):
    profile = [floor_gamma_star(r + i) + offset for i in range(r)]
    rows, loads = baseline_before_death(r, profile, death)
    fatal = loads[death][profile[death]]
    require(fatal < 0, (r, offset, "fatal sign", fatal))

    sweep = certificate_sweep(profile, rows, loads, death)
    values = [sweep[cut][0] for cut in range(1, first_negative + 1)]
    metadata = [sweep[cut][1:] for cut in range(1, first_negative + 1)]
    # Two independently assembled legacy ledgers at the hardest offset of
    # each epoch guard the sweep indexing and cumulative metadata exactly at
    # the strict sign wall without restoring the old quadratic-in-cuts cost.
    if (r, offset) in {(256, 0), (512, 4), (1024, 14)}:
        for cut in (first_negative - 1, first_negative):
            require(sweep[cut] == certificate(profile, rows, loads, death, cut),
                    (r, offset, "sweep/legacy certificate", cut))
    require(all(value >= 0 for value in values[:-1]),
            (r, offset, "pre-wall sign"))
    require(values[-1] < 0, (r, offset, "wall sign"))
    require(first_negative >= 2 and values[-2] > 0,
            (r, offset, "strict boundary"))

    # B_s is the fatal load plus nonnegative row charges from i=s onward.
    # Therefore B_s-B_(s+1)>=0.  This exact monotonicity turns the first
    # negative cut into a wall extending through the death row.
    require(all(values[i] >= values[i + 1] for i in range(len(values) - 1)),
            (r, offset, "adjoint monotonicity"))
    boundary = (values[-2], values[-1])
    return {
        "fatal": fatal,
        "boundary": boundary,
        "boundary_digest": digest_json(boundary),
        "boundary_bits": tuple(abs(value).bit_length() for value in boundary),
        "active_metadata": metadata[-1],
        "profile_head_tail": (profile[0], profile[death], profile[-1]),
    }


def main() -> None:
    ledger = []
    transitions = []
    for r, expected in EXPECTED.items():
        failed = expected["failed"]
        survivor = expected["survivor"]
        require(tuple(failed) == tuple(range(survivor)), (r, "offset interval"))
        deaths = []
        must_differ = []
        for offset, (expected_death, first_negative) in failed.items():
            profile, death = rule_a_first_death(r, offset)
            require(death == expected_death, (r, offset, "death", death))
            data = boundary_certificate(r, offset, death, first_negative)
            deaths.append(death)
            must_differ.append(first_negative - 1)
            ledger.append((r, offset, death, first_negative, data))
        survivor_profile, survivor_death = rule_a_first_death(r, survivor)
        require(survivor_death is None, (r, survivor, "survival"))
        transitions.append((
            r, tuple(deaths), survivor, tuple(must_differ),
            (survivor_profile[0], survivor_profile[-1]),
        ))

    require(transitions[1][3] == (35, 38, 39, 42, 44),
            "THM-3588 inherited vector")
    semantic = digest_json((EXPECTED_PARENT_SHA256, transitions, ledger))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== THM-3597 AMM dyadic Rule-A transition and adjoint atlas ==")
    print(f"parent_sha256_lf={EXPECTED_PARENT_SHA256}")
    for r, deaths, survivor, must_differ, survivor_degrees in transitions:
        print(
            f"R={r}|failed_offsets=0..{survivor-1}|death_rows={deaths}|"
            f"first_surviving_offset={survivor}|must_differ_by={must_differ}|"
            f"survivor_degree_endpoints={survivor_degrees}"
        )
    for r, offset, death, first_negative, data in ledger:
        print(
            f"row=R{r}/D{offset}|death={death}|first_negative={first_negative}|"
            f"must_differ_by={first_negative-1}|boundary_bits={data['boundary_bits']}|"
            f"boundary_sha256={data['boundary_digest']}|"
            f"active=(cells={data['active_metadata'][0]},mass={data['active_metadata'][1]},"
            f"max={data['active_metadata'][2]})|profile={data['profile_head_tail']}"
        )
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("status=FINITE-EXACT+VERIFIED-EXACT;positive integer Farkas walls")
    print("scope=fixed dyadic epochs/fixed Rule A;no alternative feasibility/uniform extractor/Cstar")
    print("PASS")


if __name__ == "__main__":
    main()
