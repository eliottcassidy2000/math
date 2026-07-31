#!/usr/bin/env python3
"""Independent all-D current-lcm audit of the expected-spike GF.

The primary scout uses grouped divisor-lattice Mobius inversion.  This
audit processes each divisor symbol literally, retains the current lcm as
a state, and truncates only at multiset arity four.  It compares the full

    (m2,m3,m4,c,large_capacity)

distribution for every resolving D.  A second pruned current-lcm DP keeps
at most one spike symbol and independently recovers the denominator-resolved
c=3 distribution, including its exact d=7a+r allowance.
"""

from collections import Counter
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path


PRIMARY_PATH = Path(__file__).with_name(
    "lrc14_k3_four_drift_full_activity_q7_gf.py"
)
spec = spec_from_file_location("lrc14_k3_expected_spike", PRIMARY_PATH)
primary = module_from_spec(spec)
spec.loader.exec_module(primary)
base = primary.base
q7 = primary.q7


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def independent_distribution(D):
    """Literal divisor-symbol/current-lcm DP for the complete all-c GF."""
    # State: used,lcm,m2,m3,m4,c,large_capacity.
    states = {(0, 1, 0, 0, 0, 0, 0): 1}
    maximum_states = 1
    for divisor in (d for d in base.divisors_of(D) if d > 1):
        pattern, unit_c, unit_capacity = primary.denominator_feature(D, divisor)
        additions = Counter()
        for state, multiplicity in tuple(states.items()):
            used, current_lcm, m2, m3, m4, c, capacity = state
            for copies in range(1, 5 - used):
                additions[
                    (
                        used + copies,
                        lcm(current_lcm, divisor),
                        m2 + copies * pattern[0],
                        m3 + copies * pattern[1],
                        m4 + copies * pattern[2],
                        c + copies * unit_c,
                        capacity + copies * unit_capacity,
                    )
                ] += multiplicity
        for state, multiplicity in additions.items():
            states[state] = states.get(state, 0) + multiplicity
        maximum_states = max(maximum_states, len(states))

    result = Counter()
    for state, multiplicity in states.items():
        used, current_lcm, m2, m3, m4, c, capacity = state
        if used == 4 and current_lcm == D:
            result[((m2, m3, m4), c, capacity)] += multiplicity
    return result, maximum_states


def independent_c3_distribution(D):
    """Current-lcm DP pruned to exactly three uniform and one spike mask."""
    q = D // 7
    # State: used,lcm,m2,m3,m4,c,large_capacity,unique_spike_d.
    states = {(0, 1, 0, 0, 0, 0, 0, 0): 1}
    maximum_states = 1
    for divisor in (d for d in base.divisors_of(D) if d > 1):
        pattern, unit_c, unit_capacity = primary.denominator_feature(D, divisor)
        additions = Counter()
        for state, multiplicity in tuple(states.items()):
            used, current_lcm, m2, m3, m4, c, capacity, spike_d = state
            spike_count = used - c
            for copies in range(1, 5 - used):
                new_c = c + copies * unit_c
                new_spikes = used + copies - new_c
                if new_c > 3 or new_spikes > 1:
                    continue
                require(
                    not (unit_c == 0 and spike_d),
                    ("second spike escaped pruning", D, divisor, state),
                )
                additions[
                    (
                        used + copies,
                        lcm(current_lcm, divisor),
                        m2 + copies * pattern[0],
                        m3 + copies * pattern[1],
                        m4 + copies * pattern[2],
                        new_c,
                        capacity + copies * unit_capacity,
                        divisor if unit_c == 0 else spike_d,
                    )
                ] += multiplicity
        for state, multiplicity in additions.items():
            states[state] = states.get(state, 0) + multiplicity
        maximum_states = max(maximum_states, len(states))

    result = Counter()
    for state, multiplicity in states.items():
        used, current_lcm, m2, m3, m4, c, capacity, spike_d = state
        if used != 4 or current_lcm != D or c != 3:
            continue
        require(spike_d and q % spike_d == 0, ("bad unique spike", D, state))
        result[
            (
                (m2, m3, m4),
                spike_d,
                capacity,
                primary.one_spike_allowance(q, spike_d),
            )
        ] += multiplicity
    return result, maximum_states


def main():
    by_divisor, body_count, body_divisor_rows = q7.build_rows()
    require(body_count == 3003, "body universe changed")
    require(body_divisor_rows == 251536, "body/divisor universe changed")
    require(sum(map(len, by_divisor.values())) == 26970, "row universe changed")
    require(len(by_divisor) == 217, "resolving-D universe changed")

    all_c_semantic = sha256()
    c3_semantic = sha256()
    total_features = 0
    total_c3_features = 0
    total_shapes = 0
    total_c3_shapes = 0
    maximum_states = 0
    maximum_c3_states = 0
    for D in sorted(by_divisor):
        direct, state_count = independent_distribution(D)
        direct_c3, c3_state_count = independent_c3_distribution(D)
        require(
            direct == primary.expected_distribution(D),
            ("all-c Mobius/current-lcm disagreement", D),
        )
        require(
            direct_c3 == primary.c3_denominator_distribution(D),
            ("c3 Mobius/current-lcm disagreement", D),
        )
        total_features += len(direct)
        total_c3_features += len(direct_c3)
        total_shapes += sum(direct.values())
        total_c3_shapes += sum(direct_c3.values())
        maximum_states = max(maximum_states, state_count)
        maximum_c3_states = max(maximum_c3_states, c3_state_count)
        for feature, multiplicity in sorted(direct.items()):
            all_c_semantic.update(f"{D}|{feature}|{multiplicity}\n".encode())
        for feature, multiplicity in sorted(direct_c3.items()):
            c3_semantic.update(f"{D}|{feature}|{multiplicity}\n".encode())

    require(total_shapes == 694921995, "all-c shape total changed")
    print("LRC14 k3 expected-spike GF independent current-lcm audit")
    print(f"resolving_D={len(by_divisor)}")
    print(f"all_c_features={total_features}")
    print(f"c3_denominator_features={total_c3_features}")
    print(f"all_c_shapes={total_shapes}")
    print(f"c3_shapes={total_c3_shapes}")
    print(f"maximum_all_c_lcm_states={maximum_states}")
    print(f"maximum_c3_lcm_states={maximum_c3_states}")
    print(f"all_c_semantic_sha256={all_c_semantic.hexdigest()}")
    print(f"c3_semantic_sha256={c3_semantic.hexdigest()}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
