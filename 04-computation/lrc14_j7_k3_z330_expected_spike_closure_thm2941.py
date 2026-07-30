#!/usr/bin/env python3
"""Close the unique projected k=3 row at z1=330 by expected spike.

The repaired descending scalar atlas and the z1=336 ray/status closure leave
one possible row at the next occupied height::

    E=(1,2,6,8,12,14),  L=2352,  z1=330,  d1=392.

For an actual four-drift packet its resolving denominator D is a multiple of
d1 dividing L.  The THM-2928 support-transfer cutoff 125/143 leaves the four
possibilities D=392,784,1176,2352.  Put q=D/7 and let c be the number of the
four drift denominators that do not divide q.  The fixed d1 does not divide q,
so c>=1.

For the projected support S_D, define

    N_c = #{b mod q : |S_D intersect (b mod q)| > c}.

This file implements the THM-2928 (37lS8) expected-spike predicate directly:

    N_c=0  or  55*N_c < 13*(4-c)*q.

For N_c>0, coverage puts the compact aligned-safe carrier R inside the open
event E={T>=N_c}.  If E is the whole circle then mu(E)=1>55/91; otherwise
compactness gives the strict inclusion-mass inequality mu(E)>55/91.  Markov
and the exact spike mean give the displayed strict necessary condition.  This
wording does not assume that E is always proper.

Literal residue loads and literal later-denominator three-multisets verify
all four rows independently of a generating function.  Every c=1,2,3,4
fails at every D, so the z1=330 row is empty.  The atlas says the next occupied
height is 328, proving the projected k=3 cap z1<=328.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PHYSICAL_PATH = (
    ROOT / "04-computation" / "lrc14_j7_aligned_projected_arc_suffix_thm2941.py"
)
PHYSICAL_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_aligned_projected_arc_suffix_thm2941.out"
)
SUPPORT_PATH = (
    ROOT / "04-computation" / "lrc14_two_drift_body_projection_support_thm2928.py"
)
SUPPORT_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_two_drift_body_projection_support_thm2928.out"
)
COMBINED_PATH = (
    ROOT / "04-computation" /
    "lrc14_three_drift_body_projection_fiber_thm2928.py"
)
COMBINED_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_three_drift_body_projection_fiber_thm2928.out"
)
ATLAS_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_projected_scalar_atlas_thm2941.py"
)
ATLAS_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_projected_scalar_atlas_thm2941.out"
)
Z336_CLOSURE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z336_ray_status_closure_thm2941.py"
)
Z336_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z336_ray_status_closure_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z330_expected_spike_closure_thm2941.out"
)

EXPECTED_INPUTS = (
    ("physical_source", PHYSICAL_PATH,
     "a003d287f618eb301edf6974d0b67dc128c4f380a169e7809ed5b5754e8b8303"),
    ("physical_output", PHYSICAL_OUTPUT_PATH,
     "61e16aab8a368881c574047e576645e6b41837dc9f804f7a78d37230d843612b"),
    ("support_source", SUPPORT_PATH,
     "778842c0e8e7172835ca6ae673fb6156f212d4296e672bce4e7cc2815195bf1a"),
    ("support_output", SUPPORT_OUTPUT_PATH,
     "648327d3b9b5b9a50c7760f0afd89a7a33161f57fa98c1b9e181d6b5b791a25f"),
    ("projection_source", COMBINED_PATH,
     "42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4"),
    ("projection_output", COMBINED_OUTPUT_PATH,
     "2e211620ad7064ea06f7544b5fbac709d6d52d9a0e261b464ae26b595f09b669"),
    ("atlas_source", ATLAS_PATH,
     "ddb2d19c02c4d70cfa74141265ceac585685932e85768baf4cb98aeb3e37935b"),
    ("atlas_output", ATLAS_OUTPUT_PATH,
     "ce6807de6d6b7022c97839d0bf9fc8ba3b90e7b97bc5b0d4069e88563e232be6"),
    ("z336_closure_source", Z336_CLOSURE_PATH,
     "6e991432a1ce5ec9c1cbe97199cb5fd647c5edc89bdcc60ce560382a70835fcf"),
    ("z336_closure_output", Z336_OUTPUT_PATH,
     "e257da724128208a9c80fc5f3e8f0cd4151b2073d3a8afa4814ca1e274f168ac"),
)

BODY = (1, 2, 6, 8, 12, 14)
FIRST = 330
EXPECTED_L = 2352
EXPECTED_D1 = 392
SUPPORT_CUTOFF = Q(125, 143)
EXPECTED_D = (392, 784, 1176, 2352)
EXPECTED_N = {
    392: (56, 56, 56, 56, 56),
    784: (112, 112, 112, 112, 98),
    1176: (144, 144, 144, 116, 18),
    2352: (288, 288, 260, 134, 22),
}
EXPECTED_C_COUNTS = {
    392: ((1, 84), (2, 112), (3, 70), (4, 20)),
    784: ((1, 81), (2, 113), (3, 65), (4, 15)),
    1176: ((1, 596), (2, 848), (3, 470), (4, 100)),
    2352: ((1, 569), (2, 827), (3, 440), (4, 85)),
}

POSITIVE_BODY = (1, 2, 8, 10, 12, 14)
POSITIVE_FIRST = 328
POSITIVE_D = 11760
POSITIVE_PROFILE = (1470, 2, 49, 784)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


INPUT_SHA256 = {}
for _label, _path, _expected in EXPECTED_INPUTS:
    _observed = file_sha256(_path)
    require(_observed == _expected, ("frozen input changed", _label, _observed))
    INPUT_SHA256[_label] = _observed


def require_inputs_unchanged(stage: str) -> None:
    for label, path, _expected in EXPECTED_INPUTS:
        observed = file_sha256(path)
        require(observed == INPUT_SHA256[label], ("input drift", stage, label, observed))


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot load", path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


combined = load_module("z330_projection", COMBINED_PATH)
atlas = load_module("z330_atlas", ATLAS_PATH)
support = combined.support_module
require_inputs_unchanged("after_module_load")


def divisors(number: int) -> tuple[int, ...]:
    return tuple(d for d in range(1, number + 1) if number % d == 0)


def mean_passes(N: int, c: int, q: int) -> bool:
    """Independent direct implementation of THM-2928 (37lS8)."""
    return N == 0 or 55 * N < 13 * (4 - c) * q


def atlas_counts() -> dict[int, int]:
    line = next(
        line for line in ATLAS_OUTPUT_PATH.read_text(encoding="utf-8").splitlines()
        if line.startswith("counts_by_first=")
    )
    return dict(ast.literal_eval(line.split("=", 1)[1]))


def literal_support_record(body: tuple[int, ...], D: int) -> tuple:
    L, ranges = support.safe_cell_ranges(body)
    require(L % D == 0, ("D does not divide body ruler", body, D, L))
    arcs = combined.projected_support_arcs(D, ranges)
    points = tuple(point for left, right in arcs for point in range(left, right))
    require(len(points) == len(set(points)), ("projected arcs overlap", body, D))
    require(
        len(points) == support.support_size_bitset(D, ranges),
        ("literal support size mismatch", body, D),
    )
    q = D // 7
    loads = [0] * q
    for point in points:
        loads[point % q] += 1
    literal_histogram = tuple(sorted(Counter(loads).items()))
    sweep_histogram = combined.residue_load_histogram(arcs, q)
    require(literal_histogram == sweep_histogram, ("load sweep mismatch", body, D))
    N = tuple(sum(load > c for load in loads) for c in range(5))
    return L, len(points), N, literal_histogram


def literal_c_counts(D: int, d1: int) -> Counter:
    q = D // 7
    alphabet = tuple(d for d in divisors(D) if d > 1)
    result = Counter()
    for later in combinations_with_replacement(alphabet, 3):
        if lcm(d1, *later) != D:
            continue
        c = sum(q % d != 0 for d in (d1, *later))
        result[c] += 1
    return result


def main(output: Path) -> None:
    counts = atlas_counts()
    require(counts.get(FIRST) == 1, ("z330 atlas count changed", counts.get(FIRST)))
    require(counts.get(328) == 9, ("z328 atlas count changed", counts.get(328)))
    require(all(not counts.get(first, 0) for first in range(331, 336)),
            "interstitial 331..335 row appeared")

    target_rows = tuple(row for row in atlas.body_rows(BODY) if row[5] == FIRST)
    require(len(target_rows) == 1, ("target scalar row changed", target_rows))
    target_row_text = atlas.row_text(target_rows[0])
    L = target_rows[0][3]
    d1 = L // gcd(L, FIRST)
    require((L, d1) == (EXPECTED_L, EXPECTED_D1), ("target ruler changed", L, d1))

    z336_text = Z336_OUTPUT_PATH.read_text(encoding="utf-8")
    require("totals=states:109;crude_kills:71;status_kills:38;survivors:0" in z336_text,
            "z336 closure total changed")
    require("all_exact_controls=PASS" in z336_text, "z336 closure lost PASS")

    candidate_D = tuple(D for D in divisors(L) if D % d1 == 0)
    require(candidate_D == EXPECTED_D, ("resolving divisor list changed", candidate_D))

    records = []
    total_profiles = 0
    minimum_failure_slack = None
    for D in candidate_D:
        q = D // 7
        row_L, support_count, N, histogram = literal_support_record(BODY, D)
        require(row_L == L, ("body ruler mismatch", D, row_L))
        require(Q(support_count, D) <= SUPPORT_CUTOFF,
                ("resolving D left support-transfer ledger", D, support_count))
        require(N == EXPECTED_N[D], ("target N table changed", D, N))
        require(q % d1 != 0, ("fixed denominator stopped forcing c>=1", D))

        c_counts = literal_c_counts(D, d1)
        require(tuple(sorted(c_counts.items())) == EXPECTED_C_COUNTS[D],
                ("literal c-count changed", D, c_counts))
        total_profiles += sum(c_counts.values())
        tests = []
        for c in range(1, 5):
            lhs = 55 * N[c]
            rhs = 13 * (4 - c) * q
            slack = lhs - rhs
            require(not mean_passes(N[c], c, q),
                    ("z330 expected-spike survivor", D, c, N[c], q))
            require(slack >= 0, ("strict failure arithmetic changed", D, c, slack))
            minimum_failure_slack = (
                slack if minimum_failure_slack is None
                else min(minimum_failure_slack, slack)
            )
            tests.append((c, c_counts[c], N[c], lhs, rhs, slack))
        records.append((D, q, support_count, N, tuple(tests), histogram))
    require(total_profiles == 4495, ("total literal profile count changed", total_profiles))
    require(minimum_failure_slack == 896, ("minimum strict failure changed", minimum_failure_slack))

    # Physical positive control: an explicit z1=328 exact-lcm profile passes.
    positive_rows = tuple(
        row for row in atlas.body_rows(POSITIVE_BODY) if row[5] == POSITIVE_FIRST
    )
    require(len(positive_rows) == 1, "positive scalar row disappeared")
    positive_L, _support_count, positive_N, _histogram = literal_support_record(
        POSITIVE_BODY, POSITIVE_D
    )
    positive_d1 = positive_L // gcd(positive_L, POSITIVE_FIRST)
    require(positive_d1 == POSITIVE_PROFILE[0], "positive first denominator changed")
    require(all(POSITIVE_D % d == 0 for d in POSITIVE_PROFILE),
            "positive profile has nondivisor")
    require(lcm(*POSITIVE_PROFILE) == POSITIVE_D, "positive profile lcm changed")
    positive_q = POSITIVE_D // 7
    positive_c = sum(positive_q % d != 0 for d in POSITIVE_PROFILE)
    require((positive_c, positive_N[positive_c]) == (3, 278),
            ("positive profile data changed", positive_c, positive_N))
    require(mean_passes(positive_N[positive_c], positive_c, positive_q),
            "physical positive control stopped passing")

    # Hostile strict-endpoint and zero-demand controls.
    require(55 * 39 == 13 * 3 * 55, "equality control arithmetic changed")
    require(not mean_passes(39, 1, 55), "strict equality was admitted")
    require(mean_passes(38, 1, 55), "adjacent strict positive control failed")
    require(mean_passes(0, 4, 1), "zero-demand clause failed")

    semantic = hashlib.sha256()
    semantic.update(f"target={target_row_text}\n".encode())
    for record in records:
        semantic.update(f"record={record}\n".encode())
    semantic.update(
        f"positive={POSITIVE_BODY}|{POSITIVE_FIRST}|{POSITIVE_PROFILE}|"
        f"{positive_N}|{positive_c}\n".encode()
    )
    semantic_sha256 = semantic.hexdigest()

    require_inputs_unchanged("end")
    lines = [
        "LRC14 projected k=3 z1=330 expected-spike closure",
        *[f"{label}_sha256={INPUT_SHA256[label]}" for label in sorted(INPUT_SHA256)],
        "dependency=THM-2928-(37lS8)-expected-spike;independent direct predicate",
        "strictness=if open event E is whole circle its mass is 1;otherwise compact R subset E gives strict mass",
        f"atlas_counts=z330:{counts[330]};z328:{counts[328]};z331..335:0",
        f"target_row={target_row_text}",
        f"target_fixed_denominator={d1}",
        f"resolving_D={candidate_D}",
        f"support_cutoff={SUPPORT_CUTOFF.numerator}/{SUPPORT_CUTOFF.denominator}",
        f"literal_profile_total={total_profiles}",
    ]
    for D, q, support_count, N, tests, histogram in records:
        lines.append(
            f"D={D};q={q};support={support_count};N={N};"
            f"tests=(c,multiplicity,N,55N,13(4-c)q,slack)={tests};"
            f"load_histogram={histogram}"
        )
    lines.extend(
        [
            f"minimum_failure_slack={minimum_failure_slack}",
            (
                f"physical_positive_control=E={POSITIVE_BODY};z1={POSITIVE_FIRST};"
                f"D={POSITIVE_D};profile={POSITIVE_PROFILE};c={positive_c};"
                f"N_c={positive_N[positive_c]};lhs={55*positive_N[positive_c]};"
                f"rhs={13*(4-positive_c)*positive_q}:PASS"
            ),
            "hostile_equality_control=(N,c,q)=(39,1,55):REJECTED",
            "adjacent_strict_control=(N,c,q)=(38,1,55):ACCEPTED",
            "zero_demand_control=(N,c,q)=(0,4,1):ACCEPTED",
            "literal_residue_vs_event_sweep=PASS",
            "literal_three_multiset_lcm_counts=PASS",
            "input_snapshot_guard=start/after_module_load/end:PASS",
            "conclusion=unique projected k=3 z1=330 row is empty;next occupied atlas height is 328",
            "projected_k3_first_drift_cap<=328",
            f"semantic_sha256={semantic_sha256}",
            "all_exact_controls=PASS",
        ]
    )
    text = "\n".join(lines) + "\n"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    main(arguments.output)
