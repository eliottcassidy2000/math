#!/usr/bin/env python3
"""Exact first-edge geometry-cache and Walsh-projector audit for THM-862.

This is deliberately not a depth-two component recursion.  It reuses the
frozen c=3 sheet bank, enumerates its exact cap-admissible first edges, and
checks precisely which data descend from a logical lane to the literal
geometric child (R,x_1).  It also verifies the Fourier projector attached to
every canonical affine matching code.
"""

from collections import Counter, defaultdict
from contextlib import redirect_stdout
from fractions import Fraction as F
from hashlib import sha256
from io import StringIO
from itertools import combinations, product
from pathlib import Path
from runpy import run_path


HERE = Path(__file__).resolve().parent
BASE = HERE / "lrc13_scale_three_hamming_six_sheet_classification_codex_S16.py"

# The base verifier is an executable proof artifact and prints its own report.
# Suppress that report here; every object below is recomputed by its assertions.
with redirect_stdout(StringIO()):
    bank = run_path(str(BASE))

CONTEXTS = bank["CONTEXTS"]
ROOT_LENGTH = bank["ROOT_LENGTH"]
MATCHING_CODES = bank["MATCHING_CODES"]
crt_base = bank["crt_base"]
order_type = bank["order_type"]


def ray_base(label: int, order: int, unit: int) -> int:
    """Least proper member of the exact step-39 replacement ray."""
    return 3 * label + 39 if order == 1 else crt_base(label, 3, unit)


def first_edges(row):
    labels, orders, units = row
    real_cap = F(132, 13) / ROOT_LENGTH[labels]
    cap = real_cap.numerator // real_cap.denominator
    for index, (label, order, unit) in enumerate(zip(labels, orders, units)):
        base = ray_base(label, order, unit)
        for speed in range(base, cap + 1, 39):
            yield speed, index


# A logical lane remembers its remaining actual ray languages.  Since every
# ray is {base+39k:k>=0}, the sorted base tuple is a complete, label-free
# description; residues modulo 13 recover the labels if desired.
GEOMETRY_FIBRES = defaultdict(Counter)
LANE_KEYS = set()
LOGICAL_EDGES = 0
for row in CONTEXTS:
    labels, orders, units = row
    kind = order_type(row)
    bases = tuple(ray_base(*triple) for triple in zip(labels, orders, units))
    for speed, chosen in first_edges(row):
        GEOMETRY_FIBRES[labels, speed][kind, orders[chosen]] += 1
        future_language = tuple(sorted(base for i, base in enumerate(bases) if i != chosen))
        lane_key = labels, speed, future_language
        if lane_key in LANE_KEYS:
            raise RuntimeError("two logical first edges have the same exact future language")
        LANE_KEYS.add(lane_key)
        LOGICAL_EDGES += 1


FIBRE_MULTIPLICITIES = Counter(
    sum(contributions.values()) for contributions in GEOMETRY_FIBRES.values()
)
FIBRE_TYPE_VECTORS = Counter(
    tuple(sorted(contributions.items())) for contributions in GEOMETRY_FIBRES.values()
)

EXPECTED_VECTORS = {
    ((((0, 6), 3), 4),): 5_800,
    ((((0, 6), 3), 8), (((1, 5), 3), 4)): 1_100,
    ((((0, 6), 3), 8), (((1, 5), 3), 8), (((2, 4), 3), 2)): 2_201,
    ((((1, 5), 1), 8), (((2, 4), 1), 4)): 1_275,
    ((((1, 5), 3), 4), (((2, 4), 3), 2)): 6_479,
    ((((1, 5), 3), 4),): 1_618,
    ((((2, 4), 1), 4),): 1_335,
    ((((2, 4), 3), 2),): 2_454,
}

assert LOGICAL_EDGES == 146_912
assert len(LANE_KEYS) == LOGICAL_EDGES
assert len(GEOMETRY_FIBRES) == 22_262
assert FIBRE_MULTIPLICITIES == {2: 2_454, 4: 8_753, 6: 6_479, 12: 2_375, 18: 2_201}
assert FIBRE_TYPE_VECTORS == EXPECTED_VECTORS


# The least-live-ray relation is the only canonical tournament supplied by a
# numerical recursion state.  At the root its vertices are the six rays and
# its Hamiltonian path is their increasing least-base order.  Distinct bases
# (already distinct modulo 13) make every tournament transitive.
def numeration_clock(row):
    labels, orders, units = row
    return tuple(
        labels[index]
        for index in sorted(
            range(6),
            key=lambda i: ray_base(labels[i], orders[i], units[i]),
        )
    )


CLOCKS = Counter(numeration_clock(row) for row in CONTEXTS)
assert len(CLOCKS) == 1_151
assert all(
    len(
        {
            ray_base(label, order, unit)
            for label, order, unit in zip(row[0], row[1], row[2])
        }
    )
    == 6
    for row in CONTEXTS
)


# If M is a signed matching and e_i in {+1,-1}, then
#   1_C(e)=2^(-|M|) product_(ij,s in M)(1+s e_i e_j).
# Verify directly that its Walsh support is exactly the unions of matching
# edges, with the predicted signs and no other characters.
WALSH_ROWS = []
for labels, (equations, free_labels) in MATCHING_CODES.items():
    positions = {label: i for i, label in enumerate(labels)}
    expected = {}
    for chosen_count in range(len(equations) + 1):
        for chosen in combinations(range(len(equations)), chosen_count):
            support = frozenset(
                positions[label]
                for edge_index in chosen
                for label in equations[edge_index][:2]
            )
            sign = 1
            for edge_index in chosen:
                sign *= equations[edge_index][2]
            expected[support] = F(sign, 2 ** len(equations))

    actual = {}
    valid_words = 0
    for support_bits in product((0, 1), repeat=len(labels)):
        support = frozenset(i for i, bit in enumerate(support_bits) if bit)
        coefficient = F(0)
        for word in product((-1, 1), repeat=len(labels)):
            valid = all(
                word[positions[a]] * word[positions[b]] == sign
                for a, b, sign in equations
            )
            if valid:
                character = 1
                for i in support:
                    character *= word[i]
                coefficient += character
        coefficient /= 2 ** len(labels)
        if coefficient:
            actual[support] = coefficient
    valid_words = 2 ** (len(labels) - len(equations))
    assert actual == expected
    degree_histogram = Counter(len(support) for support in actual)
    WALSH_ROWS.append(
        (
            labels,
            len(equations),
            len(free_labels),
            valid_words,
            len(actual),
            tuple(sorted(degree_histogram.items())),
        )
    )

assert Counter((row[1], row[3], row[4], row[5]) for row in WALSH_ROWS) == {
    (2, 4, 4, ((0, 1), (2, 2), (4, 1))): 1,
    (2, 8, 4, ((0, 1), (2, 2), (4, 1))): 1,
    (2, 16, 4, ((0, 1), (2, 2), (4, 1))): 2,
    (3, 8, 8, ((0, 1), (2, 3), (4, 3), (6, 1))): 3,
}


GEOMETRY_PAYLOAD = (
    "\n".join(
        f"{','.join(map(str, labels))}:{speed}:"
        + ";".join(
            f"{kind[0]}{kind[1]}D{order}={count}"
            for (kind, order), count in sorted(contributions.items())
        )
        for (labels, speed), contributions in sorted(GEOMETRY_FIBRES.items())
    )
    + "\n"
).encode()

print("SCALE_THREE_HAMMING_SIX_GEOMETRY_CACHE_WALSH")
print("scope=exact_first_edges+matching_projectors depth_two_components=not_run")
print(f"logical_first_edges={LOGICAL_EDGES} exact_lane_keys={len(LANE_KEYS)}")
print(
    f"geometric_children={len(GEOMETRY_FIBRES)} "
    f"cache_factor={F(LOGICAL_EDGES, len(GEOMETRY_FIBRES))}"
)
print(f"geometry_fibre_multiplicities={dict(sorted(FIBRE_MULTIPLICITIES.items()))}")
for vector, count in sorted(FIBRE_TYPE_VECTORS.items(), key=str):
    print(f"geometry_fibre_vector count={count} vector={vector}")
print(f"geometry_payload_sha256={sha256(GEOMETRY_PAYLOAD).hexdigest()}")
print(
    f"numeration_clocks={len(CLOCKS)} tournament=transitive "
    "score_histogram=0:1,1:1,2:1,3:1,4:1,5:1 triangles=0 scc=1^6 hamiltonian_paths=1"
)
for row in WALSH_ROWS:
    print(
        f"walsh C={row[0]} matching_edges={row[1]} free={row[2]} "
        f"valid_words={row[3]} support={row[4]} degrees={row[5]}"
    )
print("verdict=geometry_cache_exact_future_lane_quotient_false")
