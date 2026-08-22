#!/usr/bin/env python3
"""Exact compressed AMM adjoint horizons at R=2048 and R=4096.

THM-3616 companion.  This standard-library implementation hash-pins the
THM-3601 theorem, script, and output, but does not import or execute either
THM-3601 or its parent.  It independently rebuilds the exact Fibonacci--Lucas
degree words and uses the top-distance recurrences

  delta=0: M_prev=(1+z)M,
  delta=1: M_prev=Pi_{>=0}(z^-1(1+z)^2 M)

for one complete integer adjoint sweep.
"""

from __future__ import annotations

from functools import lru_cache
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import sys


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


ROOT = Path(__file__).resolve().parents[1]
PINNED_THM3601 = (
    (
        "theorem",
        ROOT / "01-canon/theorems/THM-3601-amm-r2048-terminal-failure-adjoint-horizon.md",
        "e7574bd0a2fbd75d71c9fbb4dcdfc9eef9023a0cc5dcbae3c2c1e2fc0be0447d",
    ),
    (
        "script",
        ROOT / "04-computation/amm12592_r2048_terminal_failure_adjoint_horizon_thm3601.py",
        "071db60cd42ce1d1fff7c693e2b2baeeb782046cd93169401636022da06b0205",
    ),
    (
        "output",
        ROOT / "05-knowledge/results/amm12592_r2048_terminal_failure_adjoint_horizon_thm3601.out",
        "51caae743b5fa55a4c30840baf83232ddd62626c16607179c2e1b9f129d79dc8",
    ),
)

EXPECTED_2048 = {
    "death": 508,
    "first_negative": 196,
    "q": 195,
    "profile": (1261, 1565, 2485),
    "boundary_bits": (1211, 1212),
    "boundary_digest": "f38b3838b74bc6ab557fd6919efcb42901350073745abd16d95d9f4fe7d2ab95",
    "active": (
        48828,
        383906924668156805798067812893197399113080992965104527379505084660104758411637611723600286539507103880742384734627445584706149259400431166642522540376,
        10248578385150114262857028294077844120966011618655540198805391799146850862676148277125424815376100473925059468454633373051883752561439644748608883448,
    ),
}

EXPECTED_4096 = {
    "death": 1014,
    "first_negative": 383,
    "q": 382,
    "profile": (2537, 3143, 4986),
    "boundary_bits": (2452, 2448),
    "boundary_digest": "b92b6be9ea381a4a1789f02765a561bd93b115235351aa60432ad9e50d6cc728",
    "active_cells": 199396,
    "active_digest": "fb837984d57ac02681619ae75c9c00a53918f85a9b7bf0517b21dec5682c2d64",
}

EXPECTED_SEMANTIC_SHA256 = "4cd0b160559ea7b6820b90f7ccb204d76206d3271517b6e6afa06e5dce6a5d9a"


def require(condition: bool, payload: object) -> None:
    """Optimization-safe unconditional gate."""
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n")
    return sha256(data).hexdigest()


def digest_json(value: object) -> str:
    payload = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return sha256(payload).hexdigest()


@lru_cache(maxsize=None)
def fib_pair(index: int) -> tuple[int, int]:
    """Return (F_index,F_(index+1)) by exact fast doubling."""
    if index == 0:
        return 0, 1
    first, second = fib_pair(index >> 1)
    doubled = first * (2 * second - first)
    successor = first * first + second * second
    if index & 1:
        return successor, doubled + successor
    return doubled, successor


def five_power_below_phi_2m(power: int, m_value: int) -> bool:
    """Decide 5^power <= phi^(2m) without floating point."""
    if power < 0:
        return True
    fibonacci, fibonacci_next = fib_pair(2 * m_value)
    lucas = 2 * fibonacci_next - fibonacci
    gap = 2 * 5**power - lucas
    return gap <= 0 or gap * gap < 5 * fibonacci * fibonacci


def golden_floor(m_value: int) -> int:
    power = 3 * m_value // 5
    while five_power_below_phi_2m(power + 1, m_value):
        power += 1
    while not five_power_below_phi_2m(power, m_value):
        power -= 1
    require(five_power_below_phi_2m(power, m_value), (m_value, power, "lower"))
    require(
        not five_power_below_phi_2m(power + 1, m_value),
        (m_value, power, "upper"),
    )
    return power


def degree_profile(r_value: int, offset: int) -> tuple[int, ...]:
    profile = tuple(golden_floor(r_value + row) + offset for row in range(r_value))
    require(
        all(profile[row] - profile[row - 1] in (0, 1) for row in range(1, r_value)),
        (r_value, offset, "binary degree word"),
    )
    return profile


def initial_halved_load(r_value: int, degree: int, index: int) -> int:
    raw = (
        (-1) ** (degree - index) * comb(r_value - 2 - index, degree - index)
        - comb(degree + 1, index + 1)
        + 2 * comb(degree, index)
    )
    require(raw % 2 == 0, (r_value, degree, index, "initial parity"))
    return raw // 2


def top_clamp(load: int, degree: int, distance: int) -> tuple[int, int, int]:
    """Return residual, correction, and nonnegative charge at top distance r."""
    lower_magnitude = comb(degree - 1, distance)
    upper = comb(degree - 1, distance + 1)
    correction = min(upper, max(-lower_magnitude, load))
    residual = load - correction
    charge = correction + lower_magnitude
    require(charge >= 0, (degree, distance, "charge"))
    return residual, correction, charge


def initial_top_state(r_value: int, degree: int, width: int) -> tuple[int, ...]:
    require(width <= degree - 2, (r_value, degree, width, "feed separation"))
    state = []
    for distance in range(width):
        index = degree - 1 - distance
        residual, _, _ = top_clamp(
            initial_halved_load(r_value, degree, index), degree, distance
        )
        state.append(residual)
    return tuple(state)


def first_death_from_top(
    r_value: int, profile: tuple[int, ...], scan_limit: int
) -> tuple[int | None, int, int]:
    """Find the first fatal row through scan_limit by a shrinking exact tail.

    The initial tail has scan_limit+1 cells.  Finite propagation consumes one
    remote cell per row, so every candidate fatal top output through the scan
    horizon is exact despite discarding all lower cells.
    """
    require(0 < scan_limit < r_value, (r_value, scan_limit, "scan horizon"))
    degree = profile[0]
    initial_fatal_load = initial_halved_load(r_value, degree, degree)
    initial_fatal = initial_fatal_load - min(1, max(0, initial_fatal_load))
    if initial_fatal:
        return 0, initial_fatal, scan_limit + 1

    width = scan_limit + 1
    state = initial_top_state(r_value, degree, width)
    for row in range(1, scan_limit + 1):
        fatal = state[0]
        if fatal:
            return row, fatal, width

        degree = profile[row]
        delta = degree - profile[row - 1]
        require(delta in (0, 1), (row, delta, "degree step"))
        new_width = len(state) - 1
        require(new_width <= degree - 2, (row, new_width, degree, "feed separation"))
        next_state = []
        choose_r = 1
        for distance in range(new_width):
            if delta == 0:
                load = state[distance] + state[distance + 1]
            else:
                load = (
                    (state[distance - 1] if distance else 0)
                    + 2 * state[distance]
                    + state[distance + 1]
                )
            choose_next = choose_r * (degree - 1 - distance) // (distance + 1)
            correction = min(choose_next, max(-choose_r, load))
            next_state.append(load - correction)
            choose_r = choose_next
        state = tuple(next_state)
    return None, 0, width


def previous_multiplier_level(
    level: tuple[int, ...], delta: int
) -> tuple[int, ...]:
    """Apply one compressed M_i -> M_(i-1) top-distance step."""
    result = [0] * (len(level) + 1)
    if delta == 0:
        for distance, multiplier in enumerate(level):
            result[distance] += multiplier
            result[distance + 1] += multiplier
    else:
        require(delta == 1, (delta, "adjoint degree step"))
        for distance, multiplier in enumerate(level):
            if distance:
                result[distance - 1] += multiplier
            result[distance] += 2 * multiplier
            result[distance + 1] += multiplier
    require(
        result[0] > 0 and result[-1] > 0 and all(value > 0 for value in result),
        "positive contiguous adjoint level",
    )
    return tuple(result)


def multiplier_levels(
    profile: tuple[int, ...], death: int
) -> tuple[tuple[int, ...], ...]:
    levels: list[tuple[int, ...] | None] = [None] * death
    level = (1,)
    levels[death - 1] = level
    for row in range(death - 1, 0, -1):
        delta = profile[row] - profile[row - 1]
        level = previous_multiplier_level(level, delta)
        levels[row - 1] = level
    result = tuple(level for level in levels if level is not None)
    require(len(result) == death, "complete multiplier walk")
    require(
        all(len(result[row]) == death - row for row in range(death)),
        "multiplier support widths",
    )
    return result


def adjoint_horizon(
    r_value: int, profile: tuple[int, ...], death: int
) -> dict[str, object]:
    levels = multiplier_levels(profile, death)
    state = initial_top_state(r_value, profile[0], len(levels[0]))
    contributions = [0] * death

    for row in range(1, death):
        degree = profile[row]
        delta = degree - profile[row - 1]
        width = len(levels[row])
        require(len(state) == width + 1, (row, "primal/dual widths"))
        require(width <= degree - 2, (row, width, degree, "feed separation"))
        next_state = []
        pairing = 0
        choose_r = 1
        for distance in range(width):
            if delta == 0:
                load = state[distance] + state[distance + 1]
            else:
                load = (
                    (state[distance - 1] if distance else 0)
                    + 2 * state[distance]
                    + state[distance + 1]
                )
            choose_next = choose_r * (degree - 1 - distance) // (distance + 1)
            correction = min(choose_next, max(-choose_r, load))
            next_state.append(load - correction)
            charge = correction + choose_r
            require(charge >= 0, (row, distance, "charge"))
            pairing += levels[row][distance] * charge
            choose_r = choose_next
        require(pairing >= 0, (row, "row charge pairing"))
        contributions[row] = pairing
        state = tuple(next_state)

    require(len(state) == 1, "one fatal predecessor cell")
    fatal = state[0]
    require(fatal < 0, (r_value, "fatal sign", fatal))

    ledger: dict[int, tuple[int, int, int, int]] = {}
    rhs = fatal
    active_cells = 0
    multiplier_mass = 0
    maximum_multiplier = 0
    for row in range(death - 1, 0, -1):
        rhs += contributions[row]
        level = levels[row]
        active_cells += len(level)
        multiplier_mass += sum(level)
        maximum_multiplier = max(maximum_multiplier, max(level))
        ledger[row] = (rhs, active_cells, multiplier_mass, maximum_multiplier)

    values = tuple(ledger[cut][0] for cut in range(1, death))
    negative = tuple(cut for cut in range(1, death) if ledger[cut][0] < 0)
    require(negative, "nonempty negative wall")
    first_negative = negative[0]
    require(
        negative == tuple(range(first_negative, death)),
        "single contiguous negative wall",
    )
    require(
        all(values[index] >= values[index + 1] for index in range(len(values) - 1)),
        "cut monotonicity",
    )
    boundary = (ledger[first_negative - 1][0], ledger[first_negative][0])
    require(boundary[0] > 0 > boundary[1], "strict boundary signs")
    return {
        "death": death,
        "fatal": fatal,
        "first_negative": first_negative,
        "q": first_negative - 1,
        "boundary": boundary,
        "boundary_bits": tuple(abs(value).bit_length() for value in boundary),
        "boundary_digest": digest_json(boundary),
        "active": ledger[first_negative][1:],
    }


def certify_epoch(
    r_value: int, offset: int, expected: dict[str, object]
) -> dict[str, object]:
    profile = degree_profile(r_value, offset)
    scan_limit = r_value // 3
    death, scan_fatal, initial_width = first_death_from_top(
        r_value, profile, scan_limit
    )
    require(death == expected["death"], (r_value, "death", death))
    result = adjoint_horizon(r_value, profile, death)
    require(result["fatal"] == scan_fatal, (r_value, "fatal agreement"))
    require(result["first_negative"] == expected["first_negative"], (r_value, "wall"))
    require(result["q"] == expected["q"], (r_value, "q"))
    require(result["boundary_bits"] == expected["boundary_bits"], (r_value, "bits"))
    require(
        result["boundary_digest"] == expected["boundary_digest"],
        (r_value, "boundary digest"),
    )
    profile_control = (profile[0], profile[death], profile[-1])
    require(profile_control == expected["profile"], (r_value, "profile"))
    result["profile"] = profile_control
    result["scan_limit"] = scan_limit
    result["initial_width"] = initial_width
    return result


def main() -> None:
    pinned = []
    for label, path, expected_hash in PINNED_THM3601:
        actual_hash = lf_sha256(path)
        require(actual_hash == expected_hash, (label, actual_hash, expected_hash))
        pinned.append((label, expected_hash))

    control = certify_epoch(2048, 37, EXPECTED_2048)
    require(control["active"] == EXPECTED_2048["active"], "R2048 active metadata")

    result = certify_epoch(4096, 88, EXPECTED_4096)
    require(result["active"][0] == EXPECTED_4096["active_cells"], "R4096 active cells")
    active_digest = digest_json(result["active"])

    # q_2048/R_2048=195/2048 and q_4096/R_4096=191/2048.
    # Both lie below (3-sqrt(5))/8.  The latter exact discrepancy is
    # (-577+256sqrt(5))/2048; squaring positive sides proves its sign.
    require(573 * 573 - 5 * 256 * 256 == 649, "R2048 target-side square gap")
    require(577 * 577 - 5 * 256 * 256 == 5249, "R4096 target-side square gap")
    require(2 * 195 - 382 == 8, "normalized horizon drop")
    require(4096 == 8 * 512, "absolute-error growth denominator")
    # Therefore |e_4096|-|e_2048|=1/512>0: finite-scale monotone
    # sharpening toward the golden target fails at the preregistered test.

    semantic_payload = (
        tuple(pinned),
        (
            2048,
            37,
            control["death"],
            control["first_negative"],
            control["q"],
            control["profile"],
            control["boundary"],
            control["active"],
        ),
        (
            4096,
            88,
            result["death"],
            result["first_negative"],
            result["q"],
            result["profile"],
            result["boundary"],
            result["active"],
        ),
        ("golden discrepancy", -577, 256, 2048, 5249),
        ("absolute error growth", 1, 512),
    )
    semantic = digest_json(semantic_payload)

    # Keep these two gates last during construction so a failed initial pin
    # reports both exact deterministic digests in a single payload.
    require(
        active_digest == EXPECTED_4096["active_digest"]
        and semantic == EXPECTED_SEMANTIC_SHA256,
        (
            "pin exact digests",
            active_digest,
            semantic,
            EXPECTED_4096["active_digest"],
            EXPECTED_SEMANTIC_SHA256,
        ),
    )

    print("== THM-3616 AMM binary-Sturmian R=4096 adjoint horizon ==")
    print(
        "THM3601_sha256_lf="
        f"theorem:{pinned[0][1]}|script:{pinned[1][1]}|output:{pinned[2][1]}"
    )
    print(
        "R2048_control="
        f"D0:37|death:{control['death']}|first_negative:{control['first_negative']}|"
        f"q:{control['q']}|profile:{control['profile']}"
    )
    print(
        f"R2048_boundary=signs:(+,-)|bits:{control['boundary_bits']}|"
        f"sha256:{control['boundary_digest']}|active_cells:{control['active'][0]}"
    )
    print(
        "R4096_result="
        f"D0:88|death:{result['death']}|first_negative:{result['first_negative']}|"
        f"q:{result['q']}|profile:{result['profile']}"
    )
    print(
        f"R4096_boundary=signs:(+,-)|bits:{result['boundary_bits']}|"
        f"sha256:{result['boundary_digest']}|active_cells:{result['active'][0]}|"
        f"active_sha256:{active_digest}"
    )
    print(
        "top_distance_controls="
        f"R2048_scan:{control['scan_limit']}/width:{control['initial_width']}|"
        f"R4096_scan:{result['scan_limit']}/width:{result['initial_width']}"
    )
    print(
        "golden_discrepancy="
        "q/R-phi^-2/4=(-577+256*sqrt(5))/2048<0|square_gap:5249"
    )
    print(
        "monotone_finite_scale_golden_sharpening=FAILS|"
        "abs_error_R4096-abs_error_R2048=1/512>0"
    )
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("status=FINITE-EXACT+VERIFIED-EXACT;compressed positive-integer adjoint wall")
    print(
        "scope=R2048/D0=37 control and R4096/D0=88 fixed Rule A;"
        "no alternative feasibility/asymptotic limit/uniform extractor/Cstar"
    )
    print("PASS")


if __name__ == "__main__":
    main()
