#!/usr/bin/env python3
"""Close the exceptional projected k=2 HIGH-TAIL row at z1=1736.

The exact scalar atlas leaves the body E=(1,8,10,12,13,14) with a
positive forced-high scalar gap.  This is not a geometric survivor.  If

    T = h/91 - delta(1736)

is the suffix excess required from the four later drift labels and M3 is
the sum of the unrestricted top three singleton excesses, then every
label in a scalar-eligible packet has excess at least T-M3.  Indeed, the
other three labels contribute at most M3.  Here T-M3 is positive, so the
periodic ray law delta(r+mL)=A(r)/(r+mL) makes the literal universe finite
without a label horizon.

There are 501 labels above this exact cutoff.  Only z=13260 reaches the
forced projected wall floor(13L/150)+1=13250, and exact descending-order
triple enumeration leaves one scalar packet:

    (1736, 1836, 2004, 2340, 13260).

Its lossless projected residual has mass 14/17 after fourteen body cells,
strictly above the two-aligned union cap 25/91.  A direct full-carrier
calculation independently gives projected mass one.  Thus the row is
empty uniformly over all distinct later nonaligned labels.
"""

from __future__ import annotations

import argparse
import os
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
UNIFORM = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k3_uniform_ray_status_closure_thm2941.py"
)
PROJECTED = (
    ROOT
    / "04-computation"
    / "lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py"
)
BAND_SOURCE = (
    ROOT / "04-computation" / "lrc14_j7_k2_scalar_band_1680_1742_thm2941.py"
)
BAND_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1680_1742_thm2941.out"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1736_exceptional_high_packet_closure_thm2941.out"
)

EXPECTED_UNIFORM_SHA256 = (
    "34ab29162ed33d90093e6d2bf781def36c420a1cd6596158b5d6579a3a8f3f46"
)
EXPECTED_PROJECTED_SHA256 = (
    "76f891edfcc029a08202481304a809e03e8bd81f247afaeabab685825c4d3662"
)
EXPECTED_BAND_SOURCE_SHA256 = (
    "89016f939c961fa979ec5b30812981456df5bfb2af3066f1f1b38e5a83f1a412"
)
EXPECTED_BAND_OUTPUT_SHA256 = (
    "4a36611b26585964e185bbaa3d583be3f1c67a7b608cca785920266bc217a779"
)
EXPECTED_CANDIDATE_SHA256 = (
    "0e5c14b04f2e7e345db548d7dcb0b8582c0f543b475f21b01f4baad851a4c8fd"
)
EXPECTED_PACKET_SHA256 = (
    "05d7c674460f95c286b8bf6c9145e1ab1b2dd4057b7c3fdb6615ab207de32634"
)
EXPECTED_PROFILE_SHA256 = (
    "105c65782ef21ba4880a92257f8cd6e1e84d6453c9b807e6252c4af1492a6a26"
)
EXPECTED_SEMANTIC_SHA256 = (
    "2debf1b852cb393a0ed999eada1e2e1aac645848a59b8724459a6cb9132df441"
)

QUANTIFIER = "distinct later nonaligned labels"
BODY = (1, 8, 10, 12, 13, 14)
FIRST = 1736
SUFFIX_SLOTS = 4
HIGH_WALL_RATIO = F(13, 150)
SCALAR_ETA = F(1, 91)
ALIGNED_TWO_UNION_CAP = F(25, 91)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(UNIFORM) == EXPECTED_UNIFORM_SHA256, "uniform engine changed")
require(file_sha256(PROJECTED) == EXPECTED_PROJECTED_SHA256, "projected engine changed")
require(file_sha256(BAND_SOURCE) == EXPECTED_BAND_SOURCE_SHA256, "scalar atlas changed")
require(file_sha256(BAND_OUTPUT) == EXPECTED_BAND_OUTPUT_SHA256, "scalar output changed")
U = load_module("k2_z1736_exception_uniform", UNIFORM)
P = load_module("k2_z1736_exception_projected", PROJECTED)
require(P.A.RULER == U.suffix.A.RULER, "master rulers disagree")


def atlas_control():
    rows = tuple(
        line
        for line in BAND_OUTPUT.read_text(encoding="utf-8").splitlines()
        if line.startswith("SURVIVOR;") and f"z1={FIRST};" in line
    )
    high_rows = tuple(line for line in rows if "HIGH-TAIL" in line)
    target = tuple(
        line
        for line in high_rows
        if f"E={','.join(map(str, BODY))};" in line
    )
    require(len(rows) == 15 and len(high_rows) == 6, "z1736 atlas partition changed")
    require(len(target) == 1, "exceptional atlas row changed")
    require("largest_floor=13250;" in target[0], "atlas high wall changed")
    return sha256(target[0].encode()).hexdigest()


def literal_candidates(amplitudes, L, threshold):
    candidates = []
    for residue in range(1, L):
        amplitude = amplitudes[residue]
        if amplitude <= 0:
            continue
        label = residue
        if label <= FIRST:
            label += ((FIRST + 1 - label + L - 1) // L) * L
        while amplitude / label >= threshold:
            candidates.append((amplitude / label, label, residue, amplitude))
            label += L
    candidates.sort(key=lambda row: (-row[0], row[1:]))
    return tuple(candidates)


def scalar_packets(candidates, high_floor, target):
    high = tuple(row for row in candidates if row[1] >= high_floor)
    require(len(high) == 1, ("high candidate universe changed", high))
    forced = high[0]
    others = tuple(row for row in candidates if row != forced)
    values = tuple(row[0] for row in others)
    packets = []
    nonempty_pairs = 0
    for left in range(len(others) - 2):
        for middle in range(left + 1, len(others) - 1):
            needed = target - forced[0] - values[left] - values[middle]
            low = middle + 1
            high_index = len(others)
            # Values are descending.  The half-open interval
            # [middle+1, low) is exactly the set with value >= needed.
            while low < high_index:
                pivot = (low + high_index) // 2
                if values[pivot] >= needed:
                    low = pivot + 1
                else:
                    high_index = pivot
            cutoff = low
            if cutoff > middle + 1:
                nonempty_pairs += 1
            for right in range(middle + 1, cutoff):
                chosen = (forced, others[left], others[middle], others[right])
                labels = tuple(sorted(row[1] for row in chosen))
                require(len(set(labels)) == SUFFIX_SLOTS, "duplicate literal label")
                value = sum((row[0] for row in chosen), F())
                require(value >= target, "binary-search packet is scalar-ineligible")
                packets.append((value - target, labels, tuple(chosen)))
    packets.sort()
    return high, nonempty_pairs, tuple(packets)


def projected_lower(cells, L, labels, cell_limit=None):
    common_danger = ((F(0), F(1)),)
    cells_used = 0
    for cell in cells:
        local_union = P.merge_fraction(
            [
                interval
                for label in labels
                for interval in P.phase_danger(cell, label, L)
            ]
        )
        common_danger = P.intersect_fraction(common_danger, local_union)
        cells_used += 1
        if cell_limit is not None and cells_used == cell_limit:
            break
        if cell_limit is None and P.interval_mass(common_danger) <= 1 - ALIGNED_TWO_UNION_CAP:
            break
    return 1 - P.interval_mass(common_danger), cells_used, common_danger


def direct_projected_mass(L, labels):
    carrier = tuple(
        (F(left, P.A.RULER), F(right, P.A.RULER))
        for left, right in P.A.carrier_for(BODY)
    )
    removed = P.merge_fraction(
        [interval for label in labels for interval in P.danger_fraction(label)]
    )
    residual = P.subtract_fraction(carrier, removed)
    projected = []
    for left, right in residual:
        scaled_left = L * left
        scaled_right = L * right
        for integer in range(P.floor_fraction(scaled_left), P.ceil_fraction(scaled_right)):
            piece_left = max(scaled_left, F(integer)) - integer
            piece_right = min(scaled_right, F(integer + 1)) - integer
            if piece_left < piece_right:
                projected.append((piece_left, piece_right))
    return P.interval_mass(P.merge_fraction(projected))


def profile():
    atlas_row_digest = atlas_control()
    carrier = U.suffix.A.carrier_for(BODY)
    require(P.A.carrier_for(BODY) == carrier, "carrier engines disagree")
    h = F(sum(right - left for left, right in carrier), U.suffix.A.RULER)
    L = 14 * lcm(*BODY)
    require(U.suffix.A.RULER % L == 0, "body ruler left master ruler")
    high_wall = HIGH_WALL_RATIO * L
    high_floor = max(15, high_wall.numerator // high_wall.denominator + 1)
    require((h, len(carrier), L, high_floor) == (F(811, 2548), 38, 152880, 13250), "body data changed")

    amplitudes = [F(0)]
    signs = {-1: 0, 0: 0, 1: 0}
    for residue in range(1, L):
        amplitude = residue * U.delta(carrier, h, residue)
        require(
            (residue + L) * U.delta(carrier, h, residue + L) == amplitude,
            ("ray recurrence changed", residue),
        )
        amplitudes.append(amplitude)
        signs[(amplitude > 0) - (amplitude < 0)] += 1
    require(
        all(amplitudes[L - residue] == -amplitudes[residue] for residue in range(1, L)),
        "ray antipode changed",
    )
    require(signs[1] == signs[-1] and sum(signs.values()) == L - 1, "ray signs changed")

    first_delta = U.delta(carrier, h, FIRST)
    lower = h * SCALAR_ETA
    target = lower - first_delta
    ray_heads = []
    for residue in range(1, L):
        amplitude = amplitudes[residue]
        if amplitude <= 0:
            continue
        label = residue
        if label <= FIRST:
            label += ((FIRST + 1 - label + L - 1) // L) * L
        ray_heads.append((amplitude / label, label, residue, amplitude))
    ray_heads.sort(key=lambda row: (-row[0], row[1:]))
    top3 = tuple(ray_heads[:3])
    threshold = target - sum((row[0] for row in top3), F())
    require(
        (first_delta, lower, target, threshold)
        == (F(215, 276458), F(811, 231868), F(57, 20956), F(4876247, 30609706218)),
        "scalar cutoff changed",
    )
    require(threshold > 0, "literal cutoff is not finite")

    candidates = literal_candidates(tuple(amplitudes), L, threshold)
    candidate_digest = sha256(repr(candidates).encode()).hexdigest()
    if EXPECTED_CANDIDATE_SHA256 is not None:
        require(candidate_digest == EXPECTED_CANDIDATE_SHA256, "candidate digest changed")
    require(len(candidates) == 501, "candidate count changed")
    require(len({row[2] for row in candidates}) == 501, "a ray contributes twice")
    excluded = max(
        (
            (amplitudes[residue] / (label + L), label + L, residue)
            for _value, label, residue, _amplitude in candidates
            if amplitudes[residue] > 0
        ),
        key=lambda row: (row[0], -row[1]),
    )
    require(excluded[0] < threshold, "cutoff omitted an eligible ray point")

    high, nonempty_pairs, packets = scalar_packets(candidates, high_floor, target)
    packet_digest = sha256(repr(packets).encode()).hexdigest()
    if EXPECTED_PACKET_SHA256 is not None:
        require(packet_digest == EXPECTED_PACKET_SHA256, "packet digest changed")
    require(
        high == ((F(149, 909636), 13260, 13260, F(745, 343)),),
        "forced-high candidate changed",
    )
    require(nonempty_pairs == 1 and len(packets) == 1, "packet universe changed")
    scalar_gap, suffix_labels, _chosen = packets[0]
    require(
        (scalar_gap, suffix_labels)
        == (F(91785, 20406470812), (1836, 2004, 2340, 13260)),
        "unique scalar packet changed",
    )

    labels = (FIRST, *suffix_labels)
    cells = P.body_cells(carrier, L)
    hostile_lower, hostile_cells, hostile_common = projected_lower(cells, L, labels, 13)
    projected, cells_used, common = projected_lower(cells, L, labels)
    margin = projected - ALIGNED_TWO_UNION_CAP
    direct = direct_projected_mass(L, labels)
    require(
        (len(cells), hostile_lower, hostile_cells, hostile_common)
        == (48660, F(0), 13, ((F(0), F(1)),)),
        "thirteen-cell hostile control changed",
    )
    require(
        (projected, cells_used, common, margin, direct)
        == (F(14, 17), 14, ((F(0), F(3, 17)),), F(849, 1547), F(1)),
        "projected closure changed",
    )

    result = (
        atlas_row_digest,
        h,
        len(carrier),
        L,
        high_floor,
        first_delta,
        lower,
        target,
        tuple(signs[index] for index in (1, -1, 0)),
        sha256(repr(tuple(amplitudes)).encode()).hexdigest(),
        top3,
        threshold,
        len(candidates),
        candidate_digest,
        excluded,
        high,
        nonempty_pairs,
        len(packets),
        packet_digest,
        scalar_gap,
        labels,
        len(cells),
        hostile_lower,
        projected,
        cells_used,
        common,
        margin,
        direct,
    )
    profile_digest = sha256(repr(result).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_digest == EXPECTED_PROFILE_SHA256, "profile digest changed")
    return result, profile_digest


def render(result, profile_digest):
    (
        atlas_row_digest,
        h,
        components,
        L,
        high_floor,
        first_delta,
        lower,
        target,
        signs,
        ray_digest,
        top3,
        threshold,
        candidate_count,
        candidate_digest,
        excluded,
        high,
        nonempty_pairs,
        packet_count,
        packet_digest,
        scalar_gap,
        labels,
        body_cells,
        hostile_lower,
        projected,
        cells_used,
        common,
        margin,
        direct,
    ) = result
    semantic_payload = (
        EXPECTED_UNIFORM_SHA256,
        EXPECTED_PROJECTED_SHA256,
        EXPECTED_BAND_SOURCE_SHA256,
        EXPECTED_BAND_OUTPUT_SHA256,
        QUANTIFIER,
        HIGH_WALL_RATIO,
        SUFFIX_SLOTS,
        BODY,
        FIRST,
        high_floor,
        threshold,
        candidate_count,
        candidate_digest,
        packet_count,
        packet_digest,
        projected,
        margin,
        direct,
        profile_digest,
    )
    semantic_digest = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 exceptional HIGH packet closure at z1=1736",
        f"uniform_engine_sha256={file_sha256(UNIFORM)}",
        f"projected_engine_sha256={file_sha256(PROJECTED)}",
        f"scalar_band_source_sha256={file_sha256(BAND_SOURCE)}",
        f"scalar_band_output_sha256={file_sha256(BAND_OUTPUT)}",
        f"scope=body:{','.join(map(str, BODY))};first:{FIRST};all {QUANTIFIER};no finite label horizon",
        (
            f"high_tail=wall_ratio:{ftext(HIGH_WALL_RATIO)};high_floor:{high_floor};"
            f"suffix_slots:{SUFFIX_SLOTS};obligation:at least one later label >= high_floor"
        ),
        (
            f"body=h:{ftext(h)};components:{components};L:{L};"
            f"delta1:{ftext(first_delta)};lower:{ftext(lower)};suffix_target:{ftext(target)}"
        ),
        (
            "cutoff_lemma=each chosen suffix delta >= target minus unrestricted_top3;"
            f"top3:{top3};threshold:{ftext(threshold)};positive:yes"
        ),
        (
            f"rays=signs:+{signs[0]}/-{signs[1]}/0:{signs[2]};"
            f"ray_sha256:{ray_digest};candidate_count:{candidate_count};"
            f"candidate_sha256:{candidate_digest};best_excluded:{excluded}"
        ),
        (
            f"high_candidates={high};nonempty_triple_pairs:{nonempty_pairs};"
            f"scalar_packets:{packet_count};packet_sha256:{packet_digest}"
        ),
        (
            f"unique_packet=labels:{labels};scalar_gap:{ftext(scalar_gap)};"
            "forced_high_scalar_test:SURVIVES"
        ),
        (
            f"hostile_control=first_13_cells_lower:{ftext(hostile_lower)};"
            f"body_cells:{body_cells};projected_cells_used:{cells_used}"
        ),
        (
            f"projected_lower:{ftext(projected)};common_danger:{common};"
            f"two_aligned_cap:{ftext(ALIGNED_TWO_UNION_CAP)};margin:{ftext(margin)};"
            f"direct_projected_mass:{ftext(direct)};conclusion:EMPTY"
        ),
        f"atlas_row_sha256={atlas_row_digest}",
        f"profile_sha256={profile_digest}",
        f"semantic_sha256={semantic_digest}",
        "all_exact_controls=PASS",
    ]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    result, profile_digest = profile()
    output = render(result, profile_digest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
