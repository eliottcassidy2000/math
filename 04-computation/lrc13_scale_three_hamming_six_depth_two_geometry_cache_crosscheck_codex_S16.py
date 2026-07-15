#!/usr/bin/env python3
"""Independent exact Python depth-two crosscheck for THM-862's c=3 bank.

The primary depth-two engine is C++.  This deliberately different evaluator
materializes the 110 root safe sets as complements of merged closed danger
combs, subtracts one further closed comb once per geometric key (R,x_1), and
then batches every context-specific future ray over that cached child.  It
counts logical depth-two nodes without estimating or sampling heights.
"""

from collections import Counter, defaultdict
from contextlib import redirect_stdout
from fractions import Fraction as F
from hashlib import sha256
from io import StringIO
from pathlib import Path
from runpy import run_path


P = 13
STEP = 39
HERE = Path(__file__).resolve().parent
BASE = HERE / "lrc13_scale_three_hamming_six_sheet_classification_codex_S16.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


with redirect_stdout(StringIO()):
    bank = run_path(str(BASE))

CONTEXTS = tuple(sorted(bank["CONTEXTS"]))
ROOT_LENGTH = bank["ROOT_LENGTH"]
crt_base = bank["crt_base"]
order_type = bank["order_type"]


def ray_base(label: int, order: int, unit: int) -> int:
    return 3 * label + STEP if order == 1 else crt_base(label, 3, unit)


def closed_danger_intervals(speed: int) -> tuple[tuple[F, F], ...]:
    radius = F(1, P * speed)
    answer = []
    for tooth in range(speed):
        lo = F(tooth, speed) - radius
        hi = F(tooth, speed) + radius
        if lo < 0:
            answer.extend(((F(0), hi), (lo + 1, F(1))))
        elif hi > 1:
            answer.extend(((lo, F(1)), (F(0), hi - 1)))
        else:
            answer.append((lo, hi))
    return tuple(sorted(answer))


DANGER_CACHE: dict[int, tuple[tuple[F, F], ...]] = {}


def danger(speed: int) -> tuple[tuple[F, F], ...]:
    if speed not in DANGER_CACHE:
        DANGER_CACHE[speed] = closed_danger_intervals(speed)
    return DANGER_CACHE[speed]


def strict_safe_components(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    intervals = sorted(interval for speed in speeds for interval in danger(speed))
    merged: list[list[F]] = []
    for lo, hi in intervals:
        if not merged or lo > merged[-1][1]:
            merged.append([lo, hi])
        elif hi > merged[-1][1]:
            merged[-1][1] = hi
    answer = []
    last = F(0)
    for lo, hi in merged:
        if lo > last:
            answer.append((last, lo))
        if hi > last:
            last = hi
    if last < 1:
        answer.append((last, F(1)))
    return tuple(answer)


ROOT_COMPONENTS = {
    labels: strict_safe_components(
        tuple(3 * label for label in range(1, P) if label not in labels)
    )
    for labels in ROOT_LENGTH
}

for labels, components in ROOT_COMPONENTS.items():
    longest = max(hi - lo for lo, hi in components)
    require(longest == ROOT_LENGTH[labels], "root component reconstruction mismatch")
    require(len(components) % 3 == 0, "root failed its Z/3 deck census")


def subtract_one_comb(
    components: tuple[tuple[F, F], ...], speed: int
) -> tuple[F, int]:
    """Return the exact longest child component and number of components."""
    teeth = danger(speed)
    tooth_index = 0
    longest = F(0)
    child_count = 0
    for component_lo, component_hi in components:
        while tooth_index < len(teeth) and teeth[tooth_index][1] <= component_lo:
            tooth_index += 1
        position = component_lo
        scan = tooth_index
        while scan < len(teeth) and teeth[scan][0] < component_hi:
            tooth_lo, tooth_hi = teeth[scan]
            if tooth_lo > position:
                gap_hi = min(tooth_lo, component_hi)
                if gap_hi > position:
                    longest = max(longest, gap_hi - position)
                    child_count += 1
            if tooth_hi > position:
                position = tooth_hi
            if position >= component_hi:
                break
            scan += 1
        if position < component_hi:
            longest = max(longest, component_hi - position)
            child_count += 1
    require(longest > 0 and child_count > 0, "lower-runner nonemptiness failed")
    return longest, child_count


def first_edges(row):
    labels, orders, units = row
    real_cap = F(132, 13) / ROOT_LENGTH[labels]
    cap = real_cap.numerator // real_cap.denominator
    for index, (label, order, unit) in enumerate(zip(labels, orders, units)):
        base = ray_base(label, order, unit)
        for speed in range(base, cap + 1, STEP):
            yield speed, index


GEOMETRY_FIBRES = defaultdict(list)
for row in CONTEXTS:
    for speed, chosen in first_edges(row):
        GEOMETRY_FIBRES[row[0], speed].append((row, chosen))

require(sum(map(len, GEOMETRY_FIBRES.values())) == 146_912, "first-edge mismatch")
require(len(GEOMETRY_FIBRES) == 22_262, "geometry-cache mismatch")


CHILD_DATA = {}
for key in sorted(GEOMETRY_FIBRES):
    labels, speed = key
    longest, component_count = subtract_one_comb(ROOT_COMPONENTS[labels], speed)
    real_cap = F(110, 39) / longest
    cap = real_cap.numerator // real_cap.denominator
    CHILD_DATA[key] = longest, component_count, real_cap, cap

# A deterministic endpoint-for-endpoint reconstruction sample exercises the
# cache subtraction against a fresh union of all seven prefix combs.
keys = sorted(CHILD_DATA)
SAMPLE_INDICES = tuple(sorted({0, len(keys) - 1, *range(173, len(keys), 347)}))
for index in SAMPLE_INDICES:
    labels, speed = keys[index]
    prefix = tuple(
        sorted(
            (3 * label for label in range(1, P) if label not in labels),
        )
    ) + (speed,)
    rebuilt = strict_safe_components(tuple(sorted(prefix)))
    longest = max(hi - lo for lo, hi in rebuilt)
    require(longest == CHILD_DATA[labels, speed][0], "sample child rebuild mismatch")
    require(len(rebuilt) == CHILD_DATA[labels, speed][1], "sample component-count mismatch")


def least_strictly_above(base: int, threshold: int) -> int:
    if base > threshold:
        return base
    return base + STEP * ((threshold - base) // STEP + 1)


LOGICAL_BY_TYPE = Counter()
TRANSITION_COUNTS = Counter()
DEAD_DEPTH_ONE_BY_TYPE = Counter()
DEAD_DEPTH_ONE_BY_TRANSITION = Counter()
GEOMETRIC_DEPTH_TWO = 0
GEOMETRY_DEPTH_TWO_MULTIPLICITIES = Counter()
GEOMETRY_ROWS = []

for (labels, first_speed), lanes in sorted(GEOMETRY_FIBRES.items()):
    longest, component_count, real_cap, cap = CHILD_DATA[labels, first_speed]
    local_geometries = Counter()
    logical_at_geometry = 0
    for row, chosen in lanes:
        orders = row[1]
        units = row[2]
        kind = order_type(row)
        first_order = orders[chosen]
        lane_nodes = 0
        for future, (label, order, unit) in enumerate(zip(labels, orders, units)):
            if future == chosen:
                continue
            base = ray_base(label, order, unit)
            first_future = least_strictly_above(base, first_speed)
            if first_future > cap:
                continue
            count = (cap - first_future) // STEP + 1
            lane_nodes += count
            LOGICAL_BY_TYPE[kind] += count
            TRANSITION_COUNTS[kind, first_order, order] += count
            for second_speed in range(first_future, cap + 1, STEP):
                local_geometries[second_speed] += 1
        if lane_nodes == 0:
            DEAD_DEPTH_ONE_BY_TYPE[kind] += 1
            DEAD_DEPTH_ONE_BY_TRANSITION[kind, first_order] += 1
        logical_at_geometry += lane_nodes
    GEOMETRIC_DEPTH_TWO += len(local_geometries)
    GEOMETRY_DEPTH_TWO_MULTIPLICITIES.update(local_geometries.values())
    local_payload = ",".join(
        f"{speed}:{multiplicity}" for speed, multiplicity in sorted(local_geometries.items())
    )
    GEOMETRY_ROWS.append(
        f"{','.join(map(str, labels))}:{first_speed}:"
        f"L={longest.numerator}/{longest.denominator}:C={component_count}:"
        f"B={real_cap.numerator}/{real_cap.denominator}:floor={cap}:"
        f"lanes={len(lanes)}:nodes={logical_at_geometry}:"
        f"children={sha256(local_payload.encode()).hexdigest()}"
    )


LOGICAL_DEPTH_TWO = sum(LOGICAL_BY_TYPE.values())
DEAD_DEPTH_ONE = sum(DEAD_DEPTH_ONE_BY_TYPE.values())
MIN_LONGEST = min(row[0] for row in CHILD_DATA.values())
MAX_LONGEST = max(row[0] for row in CHILD_DATA.values())
MAX_REAL_CAP = max(row[2] for row in CHILD_DATA.values())
MIN_REAL_CAP = min(row[2] for row in CHILD_DATA.values())
MAX_FLOOR_CAP = max(row[3] for row in CHILD_DATA.values())
MIN_FLOOR_CAP = min(row[3] for row in CHILD_DATA.values())
MIN_LONGEST_KEYS = sum(row[0] == MIN_LONGEST for row in CHILD_DATA.values())
MIN_LONGEST_LANES = sum(
    len(GEOMETRY_FIBRES[key])
    for key, row in CHILD_DATA.items()
    if row[0] == MIN_LONGEST
)

ROOT_COMPONENT_HISTOGRAM = Counter(map(len, ROOT_COMPONENTS.values()))
CHILD_COMPONENT_HISTOGRAM = Counter(row[1] for row in CHILD_DATA.values())
GEOMETRY_PAYLOAD = ("\n".join(GEOMETRY_ROWS) + "\n").encode()
GEOMETRY_PAYLOAD_HASH = sha256(GEOMETRY_PAYLOAD).hexdigest()

require(
    LOGICAL_BY_TYPE == {(2, 4): 3_408_353, (1, 5): 6_469_464, (0, 6): 5_114_446},
    "depth-two stratum census mismatch",
)
require(
    TRANSITION_COUNTS
    == {
        ((2, 4), 1, 1): 225_672,
        ((2, 4), 1, 3): 901_620,
        ((2, 4), 3, 1): 911_740,
        ((2, 4), 3, 3): 1_369_321,
        ((1, 5), 1, 3): 1_069_716,
        ((1, 5), 3, 1): 1_078_764,
        ((1, 5), 3, 3): 4_320_984,
        ((0, 6), 3, 3): 5_114_446,
    },
    "depth-two transition census mismatch",
)
require(LOGICAL_DEPTH_TWO == 14_992_263, "logical depth-two total mismatch")
require(DEAD_DEPTH_ONE == 0, "unexpected depth-one dead lane")
require(GEOMETRIC_DEPTH_TWO == 4_307_561, "geometric depth-two total mismatch")
require(
    GEOMETRY_DEPTH_TWO_MULTIPLICITIES
    == {1: 212_990, 2: 2_123_879, 3: 535_281, 4: 337_242,
        6: 892_407, 9: 164_606, 18: 41_156},
    "depth-two geometry multiplicity mismatch",
)
require(
    sum(multiplicity * count for multiplicity, count in GEOMETRY_DEPTH_TWO_MULTIPLICITIES.items())
    == LOGICAL_DEPTH_TWO,
    "depth-two geometry fibres do not recover the logical census",
)
require((MIN_LONGEST, MAX_LONGEST) == (F(11, 15_431), F(1, 36)), "L1 extrema mismatch")
require((MIN_REAL_CAP, MAX_REAL_CAP) == (F(1_320, 13), F(11_870, 3)), "cap extrema mismatch")
require((MIN_FLOOR_CAP, MAX_FLOOR_CAP) == (101, 3_956), "floor-cap extrema mismatch")
require((MIN_LONGEST_KEYS, MIN_LONGEST_LANES) == (1, 4), "minimum-L multiplicity mismatch")
require(
    GEOMETRY_PAYLOAD_HASH == "acebdeb61d6971b83282f6cc5113e13a5fd7d084d2ce79b9e73abf8953689097",
    "canonical geometry payload mismatch",
)

# Independently rebuild the primary scout's sole depth-two pruning witness.
# This verifies its geometry and cap without claiming that this lightweight
# script searched all 4,307,561 second-child geometries for further deaths.
DEAD_CONTEXT_INDEX = 1_448
DEAD_CONTEXT = CONTEXTS[DEAD_CONTEXT_INDEX]
require(
    DEAD_CONTEXT
    == (
        (4, 6, 7, 8, 9, 12),
        (3, 3, 3, 3, 3, 1),
        (2, 1, 1, 1, 2, 0),
    ),
    "primary dead context index mismatch",
)
DEAD_PREFIX_INSERTIONS = (14, 38)
DEAD_PREFIX = tuple(
    sorted(
        (3 * label for label in range(1, P) if label not in DEAD_CONTEXT[0]),
    )
) + DEAD_PREFIX_INSERTIONS
DEAD_COMPONENTS = strict_safe_components(tuple(sorted(DEAD_PREFIX)))
DEAD_LONGEST_INTERVAL = max(DEAD_COMPONENTS, key=lambda interval: interval[1] - interval[0])
DEAD_LONGEST = DEAD_LONGEST_INTERVAL[1] - DEAD_LONGEST_INTERVAL[0]
DEAD_REAL_CAP = F(88, 65) / DEAD_LONGEST
DEAD_CAP = DEAD_REAL_CAP.numerator // DEAD_REAL_CAP.denominator
DEAD_USED_LABELS = {9, 4}
DEAD_REMAINING = []
for label, order, unit in zip(*DEAD_CONTEXT):
    if label in DEAD_USED_LABELS:
        continue
    next_speed = least_strictly_above(ray_base(label, order, unit), 38)
    DEAD_REMAINING.append((label, order, unit, next_speed))
DEAD_LEAST_FUTURE = min(row[3] for row in DEAD_REMAINING)
require(DEAD_LONGEST_INTERVAL == (F(183, 494), F(56, 143)), "dead interval mismatch")
require(DEAD_LONGEST == F(115, 5_434), "dead component length mismatch")
require((DEAD_CAP, DEAD_LEAST_FUTURE) == (63, 70), "dead prefix cap mismatch")
require(
    DEAD_REMAINING
    == [(6, 3, 1, 70), (7, 3, 1, 73), (8, 3, 1, 76), (12, 1, 0, 75)],
    "dead prefix future-language mismatch",
)


print("SCALE_THREE_HAMMING_SIX_DEPTH_TWO_GEOMETRY_CACHE_CROSSCHECK")
print("arithmetic=integer+rational workers=1 height_sampling=none height_cutoff=none")
print(
    f"roots={len(ROOT_COMPONENTS)} first_logical=146912 "
    f"first_geometries={len(GEOMETRY_FIBRES)} rebuild_samples={len(SAMPLE_INDICES)}"
)
print(f"root_component_histogram={dict(sorted(ROOT_COMPONENT_HISTOGRAM.items()))}")
print(
    f"depth1_longest_range={MIN_LONGEST}..{MAX_LONGEST} "
    f"real_cap_range={MIN_REAL_CAP}..{MAX_REAL_CAP} "
    f"floor_cap_range={MIN_FLOOR_CAP}..{MAX_FLOOR_CAP}"
)
print(
    f"depth1_min_longest={MIN_LONGEST} geometric_multiplicity={MIN_LONGEST_KEYS} "
    f"logical_lane_multiplicity={MIN_LONGEST_LANES}"
)
print(f"depth1_component_histogram={dict(sorted(CHILD_COMPONENT_HISTOGRAM.items()))}")
print(
    f"depth2_logical_nodes={LOGICAL_DEPTH_TWO} "
    f"by_type={dict(sorted(LOGICAL_BY_TYPE.items()))}"
)
for key, count in sorted(TRANSITION_COUNTS.items()):
    print(f"depth2_transition type={key[0]} D{key[1]}->D{key[2]} nodes={count}")
print(
    f"depth1_dead_no_depth2={DEAD_DEPTH_ONE} "
    f"by_type={dict(sorted(DEAD_DEPTH_ONE_BY_TYPE.items()))}"
)
print(f"depth1_dead_by_first_order={dict(sorted(DEAD_DEPTH_ONE_BY_TRANSITION.items()))}")
print(
    f"depth2_geometric_children={GEOMETRIC_DEPTH_TWO} "
    f"multiplicities={dict(sorted(GEOMETRY_DEPTH_TWO_MULTIPLICITIES.items()))}"
)
print(f"geometry_payload_sha256={GEOMETRY_PAYLOAD_HASH}")
print(
    "primary_dead_rebuild context=1448 R=4,6,7,8,9,12 "
    "orders=3,3,3,3,3,1 units=2,1,1,1,2,0 insertions=14,38"
)
print(
    "primary_dead_longest_interval=183/494,56/143 length=115/5434 "
    "four_comb_cap=63 least_future=70 "
    "remaining=6:3:1:70,7:3:1:73,8:3:1:76,12:1:0:75"
)
print(
    "primary_dead_verdict=cap_below_least_future "
    "geometry_independently_reconstructed=true uniqueness_owned_by_primary=true"
)
print("scope=depth_two_count_only depth_two_components_not_materialized")
print(f"source_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}")
