#!/usr/bin/env python3
"""Exact pre-word ghost-midpoint atlas for THM-2706.

The old delayed sectors and pair-prefix clocks are empty at the forced
midpoint phases.  This referee therefore does not reuse their aggregate
THM-2640 unit rows.  It omits that old word, checks every remaining
pointwise source-one packet field, and carries every survivor through both
degree-one affine carry/root laws.  A new transit grammar and its recomputed
private row remain deliberately undefined.
"""

from collections import Counter
from fractions import Fraction as F
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation" / "lrc14_central_half_odometer_full_local_cycle_thm2698.py"
SPEC = importlib.util.spec_from_file_location("half", BASE)
half = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(half)

P = 13
R = P ** 6
SAFE_SPEEDS = (1, 14, 27, 40, 53, 66, 13, P ** 3, 2 * P ** 5)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def centered_distance(value):
    value = half.frac(value)
    return min(value, 1 - value)


def safe_slack(value):
    return min(
        (centered_distance(speed * value)
         - (F(1, 7) if index == 0 else F(1, 14))) / speed
        for index, speed in enumerate(SAFE_SPEEDS)
    )


def old_word_profile(y, pair_prefixes, sector_words, sector_starts):
    h = half.floor_fraction(P * y)
    kappa = half.floor_fraction(2 * P * y) - 2 * h
    sectors = tuple(
        sector for sector in range(2)
        if half.strict_interval_member(
            y * half.T, sector_words[sector], sector_starts[sector]
        )
    )
    clocks = tuple(
        (sector, tuple(
            ell for ell in range(7)
            if half.strict_interval_member(
                y, half.prefix_intervals(
                    pair_prefixes[sector][ell][h][kappa]
                )
            )
        ))
        for sector in range(2)
    )
    return h, kappa, sectors, clocks


def local_cells(y, desired_clock, rail_js, module, rails, present,
                present_starts):
    h = half.floor_fraction(P * y)
    kappa = half.floor_fraction(2 * P * y) - 2 * h
    first_n = half.floor_fraction(F(6 * R, P) - y) + 1
    last_n = half.floor_fraction(F(7 * R, P) - y)
    rail_starts = [[left for left, _, _ in rail[3]] for rail in rails[:14]]
    counts = Counter()
    cells = []
    for n in range(first_n, last_n + 1):
        x = F(n, R) + F(y, R)
        counts["source_n"] += 1
        if half.shallow(x) != desired_clock:
            continue
        counts["shallow"] += 1
        if half.owner(x) != desired_clock:
            continue
        counts["owner"] += 1
        coordinate = x * half.T
        if not half.strict_interval_member(
                coordinate, present[desired_clock, (-h) % P],
                present_starts[desired_clock, (-h) % P]):
            continue
        counts["present"] += 1
        carry = n % P
        for j in rail_js:
            if not half.strict_weighted_member(
                    coordinate, rails[j][3], rail_starts[j]):
                continue
            counts["rail"] += 1
            require(rails[j][0] == 1 and rails[j][1] == desired_clock,
                    "source-one rail/owner typing changed")
            for edge in (0, 1):
                root = (2 * carry + (2 * h + kappa) // P
                        + (edge == 0)) % P
                if root == 0 or not half.half_membership(module, x, edge, root):
                    continue
                counts["local_root"] += 1
                cells.append({
                    "x": x, "N": n, "y": y, "j": j,
                    "carry": carry, "h": h, "kappa": kappa,
                    "edge": edge, "root": root,
                })
    return cells, counts


def split_audit(source, target, macro_lift, clock, rail_js, module, rails,
                present, present_starts):
    y = half.frac(P * source["y"])
    cells, counts = local_cells(
        y, clock, rail_js, module, rails, present, present_starts
    )
    compatible = []
    source_base_carry = half.floor_fraction(P * source["y"])
    middle_base_carry = half.floor_fraction(P * y)
    for middle in cells:
        a = (middle["N"] - P * source["N"] - source_base_carry) % R
        b = (target["N"] - P * middle["N"] - middle_base_carry) % R
        require((P * a + b) % R == macro_lift,
                "macro split congruence changed")
        require(half.frac(P * source["x"] + F(a, R)) == middle["x"]
                and half.frac(P * middle["x"] + F(b, R)) == target["x"],
                "degree-one affine point maps changed")
        require(middle["carry"] == (a + source_base_carry) % P
                and target["carry"] == (b + middle_base_carry) % P,
                "degree-one carry law changed")
        source_roots = half.half_roots(
            module, half.frac(P * source["x"]), middle["edge"]
        )
        middle_roots = half.half_roots(
            module, half.frac(P * middle["x"]), target["edge"]
        )
        if (any((root + 2 * a) % P == middle["root"]
                for root in source_roots)
                and any((root + 2 * b) % P == target["root"]
                        for root in middle_roots)):
            compatible.append((middle, a, b, source_roots, middle_roots))
    require(len(compatible) == len(cells),
            "some pre-word local root failed an affine split root law")
    safe = tuple((item, safe_slack(item[0]["x"])) for item in compatible
                 if safe_slack(item[0]["x"]) > 0)
    return y, tuple(counts.items()), len(compatible), compatible[0], safe


def main():
    module, _, _, _, rails, present, present_starts = half.core.build_carrier_data()
    pair_prefixes = half.private.build_pair_prefixes(module)
    sector_words = half.prior.sector_words(module)
    sector_starts = [[left for left, _ in word] for word in sector_words]
    x0 = {
        "x": F(39123022, 82055753), "N": 2301354, "y": F(4, 17),
        "carry": 3, "edge": 0, "root": 7,
    }
    x1 = {
        "x": F(41305372, 82055753), "N": 2429727, "y": F(13, 17),
        "carry": 1, "edge": 0, "root": 4,
    }
    forward = split_audit(
        x0, x1, 4472391, 4, (8, 9), module, rails, present, present_starts
    )
    reverse = split_audit(
        x1, x0, 1956127, 1, (2, 3), module, rails, present, present_starts
    )
    forward_profile = old_word_profile(
        forward[0], pair_prefixes, sector_words, sector_starts
    )
    reverse_profile = old_word_profile(
        reverse[0], pair_prefixes, sector_words, sector_starts
    )
    require(forward_profile == (0, 1, (), ((0, ()), (1, ())))
            and reverse_profile == (12, 0, (), ((0, ()), (1, ()))),
            "forced midpoint unexpectedly entered an inherited word")

    print("LRC14 RELATIVE-SEGAL GHOST PRE-WORD TRANSIT EXACT REFEREE")
    print(f"forward_profile={forward_profile}; cascade={forward[1]}; split_compatible={forward[2]}; first={forward[3]}; safe_count={len(forward[4])}; first_safe={forward[4][0] if forward[4] else None}")
    print(f"reverse_profile={reverse_profile}; cascade={reverse[1]}; split_compatible={reverse[2]}; first={reverse[3]}; safe_count={len(reverse[4])}; first_safe={reverse[4][0] if reverse[4] else None}")
    print("scope=old sectors/prefixes are empty; no inherited or recomputed transit-word unit row is asserted")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
