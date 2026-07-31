#!/usr/bin/env python3
"""Exact affine absolute-tooth reroute probe for THM-2635.

The clock-two intersection isolates two neighboring half labels at h=3,

    (epsilon,r)=(0,11),(1,12),

which have the same absolute tooth coordinate a=r-epsilon=11=1-h.  This
companion tests the all-h section

    a=1-h,                 r=1-h+epsilon                 (mod 13)

on the complete THM-2635 half-carrier.  It uses the already proved global
primitive content 26, reports positivity and unit sections for the canonical
THM-2600 rail and for either rail, and keeps r=0 as a genuine puncture rather
than inventing an inverse.

This is a coefficient/common-x section test.  The THM-2630 predecessor digit
on it is j=7(a-kappa), so at h=3 it is j=12, not the adjacent j=4 of the
opposite-graph chronological closure.
"""

from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_old_wall_successor_sector_thm2630 as wall


P = 13
Q7 = 7
CONTENT = 26
SHARDS = cross.SHARDS
EXPECTED_SUMMARIES = {
    0: (
        (3, 4, 5, 6, 7, 8, 9, 10), (9,),
        (3, 4, 5, 6, 7, 8, 9, 10), (3, 5, 6, 7, 9),
        ((8, 19), (9, 65)),
        ((4, 1), (7, 9), (8, 22), (9, 52)),
        ((8, 2), (9, 82)),
        ((7, 2), (8, 6), (9, 76)),
    ),
    1: (
        (3, 4, 5, 6, 7, 8, 9, 10), (9,),
        (3, 4, 5, 6, 7, 8, 9, 10), (3, 5, 6, 7, 9),
        ((8, 19), (9, 65)),
        ((4, 1), (7, 9), (8, 21), (9, 53)),
        ((8, 2), (9, 82)),
        ((7, 2), (8, 6), (9, 76)),
    ),
}
EXPECTED_DIGEST = "fa251e5427a6cd90ab807c3d976afd59e5ff86fada71bae9c8f2f12111127a03"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compute_shard(bounds):
    start, stop = bounds
    (module, prefixes, _whole_prefixes, _digit_masses, rails,
     present, present_starts) = cross.build_carrier_data()
    rows = {}
    for s, ell4, root, pieces in rails[start:stop]:
        theta = (root - 12) % P
        require(theta in (0, 1), "rail theta chart changed")
        by_half = [[None] * P for _ in range(2)]
        for h in range(1, 12):
            shift = (-h) % P
            absolute_tooth = (1 - h) % P
            vectors = [[], []]
            for ell5 in range(Q7):
                overlap = cross.old.intersect_weighted_union(
                    pieces, present[ell5, shift], present_starts[ell5, shift]
                )
                for epsilon in range(2):
                    probe = (absolute_tooth + epsilon) % P
                    if probe == 0:
                        vectors[epsilon].append(0)
                        continue
                    if epsilon == 0:
                        low, high = 14 * probe, 14 * probe + 13
                    else:
                        low, high = 14 * probe - 13, 14 * probe
                    half = cross.old.intersect_weighted_comb(
                        overlap, module.C3, 182, low, high
                    )
                    value = cross.old.delayed_weighted_numerator(
                        half, prefixes[ell5][h]
                    )
                    require(value % CONTENT == 0,
                            "THM-2635 global content does not divide reroute")
                    vectors[epsilon].append(value)
            for epsilon in range(2):
                by_half[epsilon][h] = tuple(vectors[epsilon])
        key = (s, ell4, theta)
        require(key not in rows, "duplicate reroute rail key")
        rows[key] = tuple(tuple(values) for values in by_half)
    return bounds, rows


def unit(values, probe):
    if probe == 0:
        return False
    normalized = tuple(
        (value // CONTENT) * pow(probe, -1, P) % P for value in values
    )
    reduced = tuple((normalized[i] - normalized[-1]) % P
                    for i in range(Q7 - 1))
    return cross.old.sat.multiplication_determinant_7(reduced) != 0


def main():
    with ProcessPoolExecutor(max_workers=len(SHARDS)) as pool:
        results = list(pool.map(compute_shard, SHARDS))
    require(tuple(bounds for bounds, _ in results) == SHARDS,
            "reroute shards returned out of order")
    rows = {}
    for _bounds, shard in results:
        require(not set(rows).intersection(shard), "duplicate assembled reroute rail")
        rows.update(shard)
    require(len(rows) == 162, "reroute rail census changed")

    positive = defaultdict(set)
    units = defaultdict(set)
    digest = sha256()
    for key in sorted(rows):
        s, ell4, theta = key
        digest.update(bytes(key))
        for epsilon in range(2):
            for h in range(1, 12):
                values = rows[key][epsilon][h]
                require(values is not None and len(values) == Q7,
                        "missing reroute vector")
                absolute_tooth = (1 - h) % P
                probe = (absolute_tooth + epsilon) % P
                if any(values):
                    positive[s, ell4, theta, epsilon].add(h)
                if unit(values, probe):
                    units[s, ell4, theta, epsilon].add(h)
                for value in values:
                    payload = str(value).encode("ascii")
                    digest.update(len(payload).to_bytes(4, "big"))
                    digest.update(payload)

    cells = [(s, ell4) for s in range(1, P) for ell4 in range(Q7)]
    summaries = {}
    for epsilon in range(2):
        canonical_positive = [
            positive[s, ell4, wall.selected_theta(s, ell4), epsilon]
            for s, ell4 in cells
        ]
        canonical_units = [
            units[s, ell4, wall.selected_theta(s, ell4), epsilon]
            for s, ell4 in cells
        ]
        either_positive = [
            positive[s, ell4, 0, epsilon] | positive[s, ell4, 1, epsilon]
            for s, ell4 in cells
        ]
        either_units = [
            units[s, ell4, 0, epsilon] | units[s, ell4, 1, epsilon]
            for s, ell4 in cells
        ]
        summaries[epsilon] = (
            tuple(sorted(set.intersection(*map(set, canonical_positive)))),
            tuple(sorted(set.intersection(*map(set, canonical_units)))),
            tuple(sorted(set.intersection(*map(set, either_positive)))),
            tuple(sorted(set.intersection(*map(set, either_units)))),
            tuple(sorted(Counter(map(len, canonical_positive)).items())),
            tuple(sorted(Counter(map(len, canonical_units)).items())),
            tuple(sorted(Counter(map(len, either_positive)).items())),
            tuple(sorted(Counter(map(len, either_units)).items())),
        )
        require(summaries[epsilon] == EXPECTED_SUMMARIES[epsilon],
                "absolute-tooth section atlas changed")
        print(
            f"epsilon={epsilon} "
            f"canonical_positive_common={summaries[epsilon][0]} "
            f"canonical_unit_common={summaries[epsilon][1]}"
        )
        print(
            f"epsilon={epsilon} either_positive_common={summaries[epsilon][2]} "
            f"either_unit_common={summaries[epsilon][3]}"
        )
        print(
            f"epsilon={epsilon} size_hists="
            f"{summaries[epsilon][4:]}"
        )

    h = 3
    absolute_tooth = (1 - h) % P
    reroute = tuple(
        (epsilon, (absolute_tooth + epsilon) % P,
         7 * absolute_tooth % P)
        for epsilon in range(2)
    )
    require(reroute == ((0, 11, 12), (1, 12, 12)),
            "h=3 absolute-tooth/predecessor labels changed")
    print(f"h3_reroute(epsilon,r,j)={reroute}; adjacent_predecessor=4")
    print(f"coefficient_digest={digest.hexdigest()}")
    require(digest.hexdigest() == EXPECTED_DIGEST,
            "absolute-tooth coefficient digest changed")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
