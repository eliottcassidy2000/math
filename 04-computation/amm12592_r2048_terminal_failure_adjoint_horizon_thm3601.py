#!/usr/bin/env python3
"""Exact R=2048 terminal Rule-A failure and adjoint wall for THM-3601."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT / "04-computation/amm12592_dyadic_rule_a_adjoint_atlas_thm3597.py"
EXPECTED_PARENT_SHA256 = "bf6a429fb9785a02f3433d1bd2cdca8f3b8f2872bbbdc61ea8f6bcf514d5cc1e"
EXPECTED = {
    "death": 508,
    "first_negative": 196,
    "failure_profile": (1261, 1565, 2485),
    "survivor_profile": (1262, 2486),
    "boundary_bits": (1211, 1212),
    "boundary_digest": "f38b3838b74bc6ab557fd6919efcb42901350073745abd16d95d9f4fe7d2ab95",
    "active": (
        48828,
        383906924668156805798067812893197399113080992965104527379505084660104758411637611723600286539507103880742384734627445584706149259400431166642522540376,
        10248578385150114262857028294077844120966011618655540198805391799146850862676148277125424815376100473925059468454633373051883752561439644748608883448,
    ),
}
EXPECTED_SEMANTIC_SHA256 = "aac4e238138a4acde9618854eee894d08131830e4f6745dd9351432c434c1247"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


require(lf_sha256(PARENT) == EXPECTED_PARENT_SHA256, "THM-3597 parent hash")
namespace = {"__name__": "thm3601_parent", "__file__": str(PARENT)}
exec(compile(PARENT.read_text(encoding="utf-8"), str(PARENT), "exec"), namespace)

rule_a_first_death = namespace["rule_a_first_death"]
baseline_before_death = namespace["baseline_before_death"]
certificate_sweep = namespace["certificate_sweep"]
certificate = namespace["certificate"]


def main() -> None:
    r, failed_offset, survivor_offset = 2048, 37, 38
    profile, death = rule_a_first_death(r, failed_offset)
    require(death == EXPECTED["death"], ("death", death))
    rows, loads = baseline_before_death(r, profile, death)
    sweep = certificate_sweep(profile, rows, loads, death)
    values = tuple(sweep[cut][0] for cut in range(1, death))
    negative = tuple(cut for cut in range(1, death) if sweep[cut][0] < 0)
    first_negative = EXPECTED["first_negative"]
    require(negative == tuple(range(first_negative, death)), "single sign wall")
    require(all(values[index] >= values[index + 1] for index in range(len(values) - 1)),
            "adjoint cut monotonicity")

    for cut in (first_negative - 1, first_negative):
        require(sweep[cut] == certificate(profile, rows, loads, death, cut),
                ("sweep/legacy boundary", cut))

    boundary = (sweep[first_negative - 1][0], sweep[first_negative][0])
    boundary_bits = tuple(abs(value).bit_length() for value in boundary)
    require(boundary[0] > 0 > boundary[1], "strict boundary signs")
    require(boundary_bits == EXPECTED["boundary_bits"], "boundary bits")
    require(digest_json(boundary) == EXPECTED["boundary_digest"], "boundary digest")
    active = sweep[first_negative][1:]
    require(active == EXPECTED["active"], "active multiplier metadata")
    failure_profile = (profile[0], profile[death], profile[-1])
    require(failure_profile == EXPECTED["failure_profile"], "failure profile")

    survivor_profile, survivor_death = rule_a_first_death(r, survivor_offset)
    require(survivor_death is None, "D0=38 survives full epoch")
    survivor_endpoints = (survivor_profile[0], survivor_profile[-1])
    require(survivor_endpoints == EXPECTED["survivor_profile"], "survivor profile")

    semantic = digest_json((
        EXPECTED_PARENT_SHA256, r, failed_offset, survivor_offset, death,
        first_negative, boundary, active, failure_profile, survivor_endpoints,
    ))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== THM-3601 AMM R=2048 terminal-failure adjoint horizon ==")
    print(f"parent_sha256_lf={EXPECTED_PARENT_SHA256}")
    print(
        f"R={r}|failed_offset={failed_offset}|death={death}|"
        f"first_negative={first_negative}|must_differ_by={first_negative-1}|"
        f"failure_profile={failure_profile}"
    )
    print(
        f"boundary_bits={boundary_bits}|boundary_sha256={EXPECTED['boundary_digest']}|"
        f"active=(cells={active[0]},mass={active[1]},max={active[2]})"
    )
    print(
        f"surviving_offset={survivor_offset}|survivor_profile={survivor_endpoints}|"
        "sweep=legacy at both boundary cuts"
    )
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("status=FINITE-EXACT+VERIFIED-EXACT;positive integer Farkas wall")
    print("scope=R2048/D0=37,38/fixed Rule A;no alternative feasibility/uniform extractor/Cstar")
    print("PASS")


if __name__ == "__main__":
    main()
