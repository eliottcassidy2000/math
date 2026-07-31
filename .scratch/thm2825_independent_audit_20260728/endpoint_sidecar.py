#!/usr/bin/env python3
"""Selected exact endpoint-origin hostile for the THM-2825 collar.

The geometry is reconstructed from the lower physical module.  The endpoint
bank is the independent 169-address E3 quotient.  Profiles retain exact local
intersection intervals, not only Boolean hit counts.
"""

from bisect import bisect_right
from hashlib import sha256
from pathlib import Path
import ast
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    "lrc14_root_zero_full_target_semantic_clutch_20260728.py":
        "208f71020efa19fa47f66d2da061ab03fa7bc87beeb077b4008c069f499736d8",
    "lrc14_extended_carrier_endpoint_lib.py":
        "4b3f9f195b1634e1e84a1bc8bccb878a1c8c44aec13f24d197f92547c9e36c57",
}
for name, expected in PINNED.items():
    payload = (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    require(
        sha256(payload).hexdigest() == expected,
        f"endpoint dependency changed: {name}",
    )


import lrc14_extended_carrier_endpoint_lib as endpoint
import lrc14_root_zero_full_target_semantic_clutch_20260728 as physical


P = 13
T = physical.T
SHIFT = physical.SHIFT
H = T // (2 * P**5)
CELL = (1, 0, 3)
LENGTH = 26_444_880


def weighted_intersection(pieces, intervals):
    return tuple(
        physical.relative.private.old.intersect_weighted_union(
            pieces, intervals
        )
    )


def build_cell():
    (
        module,
        _prefixes,
        _whole,
        _masses,
        rails,
        present,
        _starts,
    ) = physical.relative.lift.m.core.build_carrier_data()
    pair_prefixes = physical.relative.private.build_pair_prefixes(module)
    _sv, _tv, _pairs, details = physical.overlap.overlap_vectors(
        module, pair_prefixes, rails, present, rail_index=8
    )
    full = physical.target.load_present_module()
    e3 = physical.target.exclusive_source(full, 3)
    clocks = tuple(
        full.make_comb(
            full.C1, 182, 26 * clock - 13, 26 * clock + 13
        )
        for clock in range(7)
    )
    clock, s, t = CELL
    section = physical.target.source_present_section(
        full, e3, clock, s, t, clocks
    )
    source = weighted_intersection(details[clock][0], section)
    target = weighted_intersection(details[clock][1], section)
    pullback = physical.overlap.shift_weighted(target, -SHIFT)
    aligned = physical.overlap.intersect_weighted_profiles(
        source, pullback
    )
    common = tuple(
        (left, right, a)
        for left, right, a, b in aligned
        if a == b
    )
    right = physical.subtract_weighted(pullback, common)
    require(
        len(common) == 239
        and len(right) == 2
        and all(piece[1] - piece[0] == LENGTH
                for piece in common + right),
        "selected weighted cell changed",
    )
    return common, right


def index_present(bank):
    return {
        address: (
            intervals,
            tuple(right for _left, right in intervals),
        )
        for address, intervals in bank.items()
    }


def local_profile(interval, indexed):
    intervals, rights = indexed
    left, right = interval
    index = bisect_right(rights, left)
    rows = []
    while index < len(intervals) and intervals[index][0] < right:
        a = max(left, intervals[index][0])
        b = min(right, intervals[index][1])
        if a < b:
            rows.append((a - left, b - left))
        index += 1
    return tuple(rows)


def profile(piece, bank, target=False):
    left, right = piece[:2]
    if target:
        left += SHIFT
        right += SHIFT
    return tuple(
        local_profile((left, right), bank[address])
        for address in endpoint.KEYS
    )


def translated_address_deltas(before, after):
    before_bank = dict(zip(endpoint.KEYS, before))
    after_bank = dict(zip(endpoint.KEYS, after))
    return tuple(
        delta
        for delta in endpoint.KEYS
        if all(
            after_bank[address]
            == before_bank[
                (
                    (address[0] + delta[0]) % P,
                    (address[1] + delta[1]) % P,
                )
            ]
            for address in endpoint.KEYS
        )
    )


def profile_type(row):
    return (
        sum(not value for value in row),
        sum(value == ((0, LENGTH),) for value in row),
        sum(
            bool(value) and value != ((0, LENGTH),)
            for value in row
        ),
    )


def main():
    common, right = build_cell()
    common_by_left = {piece[0]: piece for piece in common}
    present = index_present(endpoint.present_cache())
    rows = []
    for r_index, r_piece in enumerate(right):
        for step in (1, 2, 14):
            target_piece = common_by_left.get(
                (r_piece[0] + step * H) % T
            )
            require(
                target_piece is not None,
                f"selected {step}h collar left common support",
            )
            source_before = profile(r_piece, present)
            source_after = profile(target_piece, present)
            target_before = profile(r_piece, present, target=True)
            target_after = profile(target_piece, present, target=True)
            row = (
                r_index,
                step,
                r_piece,
                target_piece,
                profile_type(source_before),
                profile_type(source_after),
                sum(a != b for a, b in zip(
                    source_before, source_after
                )),
                translated_address_deltas(
                    source_before, source_after
                ),
                profile_type(target_before),
                profile_type(target_after),
                sum(a != b for a, b in zip(
                    target_before, target_after
                )),
                translated_address_deltas(
                    target_before, target_after
                ),
            )
            rows.append(row)
    require(
        len(rows) == 6
        and all(row[6] > 0 and row[10] > 0 for row in rows)
        and all(not row[7] and not row[11] for row in rows),
        "endpoint hostile acquired a fixed or translated identification",
    )
    tree = ast.parse(Path(__file__).read_text())
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "endpoint audit contains an assert node",
    )
    print("THM-2825 SELECTED ENDPOINT-ORIGIN HOSTILE")
    print(f"dependency_sha256={tuple(PINNED.items())}")
    print(f"cell={CELL};common={len(common)};right={len(right)}")
    print(
        "columns=(right_index,step,piece,image,"
        "source_types_before_after,source_hamming,source_address_deltas,"
        "target_types_before_after,target_hamming,target_address_deltas)"
    )
    for row in rows:
        print(f"endpoint_row={row}")
    print(
        "verdict=no +h,+2h,or surviving +14h collar on the selected "
        "cell preserves the fixed endpoint profile or becomes one under "
        "any F13^2 address translation"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
