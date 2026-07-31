#!/usr/bin/env python3
"""Independent all-D audit of the k3/four-drift weighted shape GF.

The primary scout extracts exact-lcm coefficients by grouped
divisor-lattice Mobius inversion.  This audit does not use that recurrence.
It processes each divisor symbol once and propagates the literal current lcm
as a state, truncated at multiset arity four.  For every resolving
denominator in the exact k3 support universe, it compares the complete
weighted feature distribution

    ((multiplicity d=2, d=3, d=4), capacity from d>4) -> shape count.

This is an algebraically independent path for the load-bearing compression,
not just an unweighted total check.
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIMARY_PATH = Path(__file__).with_name(
    "lrc14_k3_four_drift_divisor_status_gf.py"
)
spec = spec_from_file_location("lrc14_k3_primary_gf", PRIMARY_PATH)
primary = module_from_spec(spec)
spec.loader.exec_module(primary)
support = primary.support


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def independent_distribution(D):
    """Literal divisor-symbol DP with current lcm retained."""
    # State: used, current_lcm, count2, count3, count4, large_capacity.
    states = {(0, 1, 0, 0, 0, 0): 1}
    maximum_states = len(states)
    for divisor in (d for d in primary.divisors_of(D) if d > 1):
        capacity = (D // divisor) * ((divisor + 6) // 7)
        old_states = tuple(states.items())
        additions = Counter()
        for state, multiplicity in old_states:
            used, current_lcm, count2, count3, count4, large = state
            for copies in range(1, 5 - used):
                additions[
                    (
                        used + copies,
                        lcm(current_lcm, divisor),
                        count2 + copies * (divisor == 2),
                        count3 + copies * (divisor == 3),
                        count4 + copies * (divisor == 4),
                        large + copies * capacity * (divisor > 4),
                    )
                ] += multiplicity
        for state, multiplicity in additions.items():
            states[state] = states.get(state, 0) + multiplicity
        maximum_states = max(maximum_states, len(states))

    result = Counter()
    for state, multiplicity in states.items():
        used, current_lcm, count2, count3, count4, large = state
        if used == 4 and current_lcm == D:
            result[((count2, count3, count4), large)] += multiplicity
    return result, maximum_states


def main():
    resolving = set()
    body_count = 0
    body_divisor_rows = 0
    support_rows = 0
    for body in combinations(range(1, 15), 6):
        body_count += 1
        L, ranges = support.safe_cell_ranges(body)
        for D in primary.divisors_of(L):
            body_divisor_rows += 1
            support_count = support.support_size_bitset(D, ranges)
            if Q(support_count, D) > primary.SUPPORT_CUTOFF:
                continue
            support_rows += 1
            resolving.add(D)

    require(body_count == 3003, "body universe changed")
    require(body_divisor_rows == 251536, "body/divisor universe changed")
    require(support_rows == 26970, "k3 support-row universe changed")
    require(len(resolving) == 217, "k3 resolving alphabet changed")

    semantic = sha256()
    total_shapes = 0
    total_features = 0
    maximum_states = 0
    maximum_feature_D = None
    maximum_features = 0
    for D in sorted(resolving):
        independent, state_count = independent_distribution(D)
        mobius = primary.shape_distribution(D)
        require(
            independent == mobius,
            ("weighted exact-lcm distributions disagree", D),
        )
        total_shapes += sum(independent.values())
        total_features += len(independent)
        maximum_states = max(maximum_states, state_count)
        if len(independent) > maximum_features:
            maximum_features = len(independent)
            maximum_feature_D = D
        for feature, multiplicity in sorted(independent.items()):
            semantic.update(f"{D}|{feature}|{multiplicity}\n".encode())

    require(total_shapes == 694921995, "weighted shape total changed")
    require(total_features == 121260, "weighted feature total changed")
    print("LRC14 k3/four-drift weighted GF independent exact-lcm audit")
    print(f"resolving_D={len(resolving)}")
    print(f"support_rows={support_rows}")
    print(f"total_shapes={total_shapes}")
    print(f"total_features={total_features}")
    print(f"maximum_lcm_dp_states={maximum_states}")
    print(f"maximum_feature_D={maximum_feature_D}")
    print(f"maximum_features={maximum_features}")
    print(f"semantic_sha256={semantic.hexdigest()}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
