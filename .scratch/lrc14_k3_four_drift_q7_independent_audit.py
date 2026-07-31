#!/usr/bin/env python3
"""Independent coefficientwise audit of the all-D q=D/7 compression.

The primary q7 scout obtains exact-lcm coefficients by grouped Mobius
inversion.  This audit instead processes every allowed divisor symbol once,
retains the literal current lcm, and permits zero through four copies of that
symbol.  It compares every (uniform baseline, seven-spike cardinality)
coefficient for all 217 resolving denominators and independently rebuilds
the three default-stage counts and digests.
"""

from collections import Counter
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
Q7_PATH = ROOT / ".scratch" / "lrc14_k3_four_drift_q7_all_D_gf.py"
spec = spec_from_file_location("lrc14_k3_q7_under_audit", Q7_PATH)
q7 = module_from_spec(spec)
spec.loader.exec_module(q7)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def direct_distribution(D):
    """Literal divisor-symbol/current-lcm DP, with no Mobius inversion."""
    q = D // 7
    # State: used, current_lcm, uniform baseline b, total spike count R.
    states = {(0, 1, 0, 0): 1}
    maximum_states = 1
    for divisor in (d for d in q7.base.divisors_of(D) if d > 4):
        ell = (divisor + 6) // 7
        if q % divisor == 0:
            unit_b = 0
            unit_R = (q // divisor) * ell
        else:
            unit_b = 1
            unit_R = 0
        additions = Counter()
        for state, multiplicity in tuple(states.items()):
            used, current_lcm, baseline, spike_count = state
            for copies in range(1, 5 - used):
                additions[
                    (
                        used + copies,
                        lcm(current_lcm, divisor),
                        baseline + copies * unit_b,
                        spike_count + copies * unit_R,
                    )
                ] += multiplicity
        for state, multiplicity in additions.items():
            states[state] = states.get(state, 0) + multiplicity
        maximum_states = max(maximum_states, len(states))

    result = Counter()
    for state, multiplicity in states.items():
        used, current_lcm, baseline, spike_count = state
        if used == 4 and current_lcm == D:
            result[(baseline, spike_count)] += multiplicity
    return result, maximum_states


def digest_row(digest, D, record, count):
    support_count, body, L, _arcs = record
    digest.update(f"{body}|{L}|{D}|{support_count}|{count}\n".encode())


def main():
    by_divisor, body_count, body_divisor_rows = q7.build_rows()
    require(body_count == 3003, "body universe changed")
    require(body_divisor_rows == 251536, "body/divisor universe changed")
    require(sum(map(len, by_divisor.values())) == 26970, "support rows changed")
    require(len(by_divisor) == 217, "resolving-D universe changed")

    feature_digest = sha256()
    stage_digest = {
        name: sha256()
        for name in ("raw_no_small", "scalar_no_small", "q7_and_scalar")
    }
    occurrences = Counter()
    shape_features = {
        name: set()
        for name in ("raw_no_small", "scalar_no_small", "q7_and_scalar")
    }
    feature_multiplicity = {}
    rows = {
        name: set()
        for name in ("raw_no_small", "scalar_no_small", "q7_and_scalar")
    }
    divisors = {
        name: set()
        for name in ("raw_no_small", "scalar_no_small", "q7_and_scalar")
    }
    total_features = 0
    maximum_states = 0

    for D in sorted(by_divisor):
        q = D // 7
        direct, state_count = direct_distribution(D)
        require(
            direct == q7.canonical_q7_distribution(D),
            ("direct/Mobius q7 coefficients disagree", D),
        )
        maximum_states = max(maximum_states, state_count)
        total_features += len(direct)
        for feature, multiplicity in sorted(direct.items()):
            feature_digest.update(f"{D}|7|{feature}|{multiplicity}\n".encode())
            shape_features["raw_no_small"].add((D, feature))
            feature_multiplicity[(D, feature)] = multiplicity
        divisors["raw_no_small"].add(D)

        for record in by_divisor[D]:
            support_count, body, _L, arcs = record
            histogram = q7.combined.residue_load_histogram(arcs, q)
            stage_counts = Counter()
            for (baseline, spike_count), multiplicity in direct.items():
                ambient = baseline * q + 7 * spike_count
                stage_counts["raw_no_small"] += multiplicity
                if ambient >= support_count:
                    stage_counts["scalar_no_small"] += multiplicity
                    shape_features["scalar_no_small"].add(
                        (D, (baseline, spike_count))
                    )
                    if q7.canonical_q7_passes(
                        q, baseline, spike_count, histogram
                    ):
                        stage_counts["q7_and_scalar"] += multiplicity
                        shape_features["q7_and_scalar"].add(
                            (D, (baseline, spike_count))
                        )
            for stage, count in stage_counts.items():
                if not count:
                    continue
                occurrences[stage] += count
                rows[stage].add((body, D))
                divisors[stage].add(D)
                digest_row(stage_digest[stage], D, record, count)

    expected = {
        "raw_no_small": (
            20129326373,
            672632936,
            26970,
            217,
            "dbb12162c4c29bfdb9d3e0a4447b7704965cb070e4db0ea0730d85722992be8c",
        ),
        "scalar_no_small": (
            12852450428,
            672351643,
            18599,
            186,
            "3327619ea1133c4bf17bd11c675c604eef18891fb63f6837e695669cf5592b2b",
        ),
        "q7_and_scalar": (
            1437089787,
            171159117,
            2209,
            112,
            "01a33f9973ddf12cf7175cc8acc1a9bd17f83d1c700fa91decd7c3ac14c44e04",
        ),
    }
    for stage, (occ, shapes, row_count, divisor_count, digest) in expected.items():
        stage_shape_count = sum(
            feature_multiplicity[key] for key in shape_features[stage]
        )
        require(occurrences[stage] == occ, ("occurrence count changed", stage))
        require(
            stage_shape_count == shapes,
            ("shape count changed", stage),
        )
        require(len(rows[stage]) == row_count, ("row count changed", stage))
        require(
            len(divisors[stage]) == divisor_count,
            ("divisor count changed", stage),
        )
        require(
            stage_digest[stage].hexdigest() == digest,
            ("stage digest changed", stage),
        )

    require(
        feature_digest.hexdigest()
        == "d88f0bcc0f29f1a3eccb101c456f775038938c0c13aa02bd031729b6bb38379b",
        "feature digest changed",
    )
    print("LRC14 k3/four-drift q7 independent all-D audit")
    print(f"resolving_D={len(by_divisor)}")
    print(f"total_features={total_features}")
    print(f"maximum_lcm_dp_states={maximum_states}")
    print(f"feature_semantic_sha256={feature_digest.hexdigest()}")
    for stage in expected:
        stage_shape_count = sum(
            feature_multiplicity[key] for key in shape_features[stage]
        )
        print(
            f"stage={stage},occurrences={occurrences[stage]},"
            f"shapes={stage_shape_count},rows={len(rows[stage])},"
            f"divisors={len(divisors[stage])},"
            f"semantic_sha256={stage_digest[stage].hexdigest()}"
        )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
