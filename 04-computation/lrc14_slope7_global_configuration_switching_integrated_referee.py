#!/usr/bin/env python3
"""Exact full-union slope-seven service referee after THM-2672.

This rebuilds all 162 THM-2640 rails and independently classifies every raw
seven-clock row by Euclidean gcd with Phi_7 in F_13[z].  It then forgets the
configuration completely and unions all rails, both delayed sectors, and both
deep half-edges at each physical future half-digit.
"""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from math import gcd
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))
import lrc14_predecessor_carry_private_root_atlas_thm2640 as m

P = 13
PHI7 = (1, 1, 1, 1, 1, 1, 1)
FULL = (1 << P) - 1


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def trim(poly):
    out = [value % P for value in poly]
    while out and out[-1] == 0:
        out.pop()
    return tuple(out)


def polynomial_remainder(dividend, divisor):
    out = list(trim(dividend))
    divisor = trim(divisor)
    require(divisor, "division by zero polynomial")
    inverse = pow(divisor[-1], -1, P)
    while len(out) >= len(divisor):
        coefficient = out[-1] * inverse % P
        shift = len(out) - len(divisor)
        for index, value in enumerate(divisor):
            out[index + shift] = (
                out[index + shift] - coefficient * value
            ) % P
        while out and out[-1] == 0:
            out.pop()
    return tuple(out)


def manual_phi7_unit(values, root, content):
    """Classify the normalized row without the canonical determinant."""
    if root == 0:
        return False
    inverse_root = pow(root, -1, P)
    coefficients = tuple(
        (value // content) * inverse_root % P for value in values
    )
    # Modulo Phi_7, z^6=-(1+z+...+z^5).
    row = trim(tuple(
        (coefficients[index] - coefficients[6]) % P
        for index in range(6)
    ))
    phi = PHI7
    while row:
        phi, row = row, polynomial_remainder(phi, row)
    return len(trim(phi)) == 1


def rebuild_rows():
    with ProcessPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(m.shard, m.core.SHARDS))
    content = 0
    metadata = []
    rows = []
    for bounds, local_content, _, _, _, local_metadata, local_rows in results:
        require(bounds in m.core.SHARDS, "unexpected shard bounds")
        content = gcd(content, local_content)
        metadata.extend(local_metadata)
        rows.extend(local_rows)
    require(content == 26, "global primitive content changed")
    require(len(metadata) == len(rows) == 162, "full rail carrier changed")
    return content, tuple(metadata), tuple(rows)


def expected_missing(digit):
    if digit == 0:
        return (0,)
    if digit <= 11:
        return (0, 6)
    if digit <= 13:
        return (6,)
    if digit <= 24:
        return (6, 12)
    return (12,)


def main():
    content, metadata, rows = rebuild_rows()
    raw_masks = {}
    unit_masks = {}
    configuration_masks = {}
    raw_row_count = 0
    unit_row_count = 0
    manual_gcd_checks = 0

    for h in range(P):
        for kappa in range(2):
            raw_mask = 0
            unit_mask = 0
            masks = []
            for j in range(len(rows)):
                for sector in range(2):
                    for edge in range(2):
                        config_mask = 0
                        for carry in range(P):
                            values = rows[j][sector][edge][carry][kappa][h]
                            root = (
                                2 * carry + (2 * h + kappa) // P
                                + (edge == 0)
                            ) % P
                            raw = any(values)
                            manual = manual_phi7_unit(values, root, content)
                            canonical = m.is_unit(values, root, content)
                            require(manual == canonical,
                                    "manual Phi7 gcd disagrees with determinant")
                            require(not manual or raw,
                                    "unit row appeared without raw support")
                            raw_row_count += int(raw)
                            unit_row_count += int(manual)
                            manual_gcd_checks += 1
                            raw_mask |= int(raw) << carry
                            unit_mask |= int(manual) << carry
                            config_mask |= int(manual) << carry
                        masks.append(config_mask)
            raw_masks[h, kappa] = raw_mask
            unit_masks[h, kappa] = unit_mask
            configuration_masks[h, kappa] = tuple(masks)

    require(manual_gcd_checks == 162 * 2 * 2 * 13 * 2 * 13,
            "seven-clock row universe changed")
    require(raw_masks == unit_masks,
            "raw and unit full-configuration service differ")

    service_table = []
    for digit in range(26):
        h, kappa = divmod(digit, 2)
        mask = unit_masks[h, kappa]
        missing = tuple(carry for carry in range(P)
                        if not (mask >> carry) & 1)
        require(missing == expected_missing(digit),
                "full-configuration missing-carry law changed")
        service_table.append((digit, h, kappa, mask.bit_count(), missing))

    orbit_checks = 0
    for source_carry in range(P):
        orbit = tuple(
            (source_carry + 7 * delta) % P for delta in range(P)
        )
        require(tuple(sorted(orbit)) == tuple(range(P)),
                "slope-seven target orbit stopped being bijective")
        orbit_checks += 1

    # Explicitly audit the proposed THM-2648/2656 two-rainbow bridge.  It
    # cannot start: every pair is contained in the already proper union of
    # all 648 configurations at the fixed half-digit.
    full_pair_count = 0
    best_pair_size_hist = Counter()
    for key, masks in configuration_masks.items():
        best = 0
        for index, first in enumerate(masks):
            for second in masks[index:]:
                union = first | second
                best = max(best, union.bit_count())
                full_pair_count += int(union == FULL)
        best_pair_size_hist[best] += 1
    require(full_pair_count == 0,
            "two physical configurations acquired full carry service")
    require(best_pair_size_hist == Counter({11: 22, 12: 4}),
            "two-configuration best service histogram changed")

    service_size_hist = tuple(sorted(Counter(
        row[3] for row in service_table
    ).items()))
    require(service_size_hist == ((11, 22), (12, 4)),
            "full-union service histogram changed")

    print("THM2672 HETEROGENEOUS FULL-UNION SERVICE EXACT REFEREE")
    print(f"content={content} rails={len(metadata)} sectors=2 edges=2")
    print(f"raw_rows={raw_row_count} unit_rows={unit_row_count} "
          f"manual_phi7_gcd_checks={manual_gcd_checks}")
    print(f"service_table_d_h_kappa_size_missing={tuple(service_table)}")
    print("compressed_missing_law=d0:{0}; d1..11:{0,6}; "
          "d12..13:{6}; d14..24:{6,12}; d25:{12}")
    print(f"service_size_hist={service_size_hist}")
    print(f"raw_missing_equals_unit_missing={raw_masks == unit_masks}")
    print(f"slope7_orbit_checks={orbit_checks}")
    print(f"two_configuration_full_cover_count={full_pair_count} "
          f"best_pair_size_hist={tuple(sorted(best_pair_size_hist.items()))}")
    print("mechanism=the slope-seven lift fixes the physical future half-digit "
          "but permutes all predecessor carries; every fixed half-digit has "
          "an all-configuration seam puncture among carries {0,6,12}")
    print("first_failed_predicate=raw THM2640 row service itself; "
          "rail/present common components, delayed gains, and Helly balance "
          "are never reached")
    print("scope=no arbitrary per-label THM2640 slope-seven 13-fold chart; "
          "with THM2672's positive 12-fold witness the exact maximum is 12")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()


