#!/usr/bin/env python3
"""Exact projected k=2 closure of the integer band 1743 <= z1 <= 1749.

The global stage reconstructs the rigorous scalar upper envelope on all
C(14,6)=3003 literal six-body carriers and all seven first-label heights.  Its
finite horizon is used only with the proved four-slot omitted-tail bound from
the pinned THM-2941 scalar atlas.  Exactly two rows survive: one at z1=1746
and one at z1=1743.

The all-label denominator-ray quotient kills the z1=1746 row by crude fibre
capacity.  At z1=1743 it leaves seven common-status states.  Positive scalar
slack gives a finite literal gate on every residue ray; the resulting seven
packets all have lossless projected residual mass at least 25/91.  Direct
global interval subtraction independently checks the projected cell
calculation on the minimum-margin packet.

This is a scoped exact closure of a contiguous seven-height band.  It does not
by itself close z1=1750 or any height above 1749 and therefore does not change
the current global k=2 cap.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement
from math import comb, gcd, lcm
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_z1750_literal_packet_projected_closure_thm2941.py"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_band_1743_1749_literal_packet_closure_thm2941.out"
)
EXPECTED_SUPPORT_SHA256 = (
    "70c9faa37a0524673e8178ed82cc6abc040438fff043b944a7ee0227d48c8997"
)

START = 1743
END = 1749
TARGET_ROWS = (
    ((1, 4, 8, 10, 12, 14), 1743),
    ((1, 6, 9, 10, 12, 14), 1746),
)
EXPECTED_SURVIVORS = (
    (
        (1, 4, 8, 10, 12, 14),
        1743,
        F(239, 683256),
        (
            (F(3583, 3148740), 1836, "EXACT"),
            (F(1867, 2119740), 2060, "EXACT"),
            (F(263, 310415), 2172, "EXACT"),
            (F(403, 514500), 1750, "EXACT"),
        ),
        F(974213104567, 243613133253000),
        F(1049, 267540),
        F(247366239721, 3166970732289000),
        F(1049, 2940),
        34,
        11760,
        1020,
    ),
    (
        (1, 6, 9, 10, 12, 14),
        1746,
        F(851, 998130),
        (
            (F(1369, 1574370), 1836, "EXACT"),
            (F(2267, 2793735), 1810, "EXACT"),
            (F(8333, 11174940), 2172, "EXACT"),
            (F(3743, 5155290), 2004, "EXACT"),
        ),
        F(7395599147, 1846433101212),
        F(1601, 401310),
        F(1909377817, 120018151578780),
        F(1601, 4410),
        32,
        17640,
        1529,
    ),
)

# scalar, scalar digest, crude kills, crude digest, status kills, status
# digest, survivors, survivor digest.
EXPECTED_STAGES = {
    TARGET_ROWS[0]: (
        11,
        "2af40e3595b2bd431e441ef48f06d9c7f80b3e1b9d18d32369154b7c0eb53185",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        4,
        "19b2fc11b7b411b910aa25c9efb7a27f23787fe9ceafdef08322ee2ceff67568",
        7,
        "f8510fddde7fe3b6cb43a6bab1e996b37782b5a4f3cc3799ee1467001805d960",
    ),
    TARGET_ROWS[1]: (
        1,
        "021628300f498327505c8e12fadaa662a557190358ae39283bc7d0d86a85ff94",
        1,
        "e191526d2e1779a56eb276945530af9a0c5069f750f0673b9b5b2dce4be917da",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
    ),
}
EXPECTED_PACKET_SUMMARY = (
    7,
    F(1026, 16471),
    2,
    F(89330944589, 126678829291560),
    3,
    1,
    (
        F(1026, 16471),
        (1743, 1836, 2060, 2172, 2534),
        F(61, 181),
        1,
    ),
    F(1),
    "2995b91c4bd0e5119db76a417023f1524f1cef7623822a8c96bc78339aaefefb",
    "bf6df9ab5355a8826166f6647691979480a4c53f75f92589432b3bf31c921d2c",
)
EXPECTED_GLOBAL_PROFILE_SHA256 = (
    "adc43465d860ef3cd9744bf382ebb370dbe3c504697b4ebfa062cfd5ad316c59"
)
EXPECTED_GLOBAL_SURVIVOR_SHA256 = (
    "33a833973a66b4a2a98a6c371b7e96f30e3db4264c3ec1c77b6a9b13ffac4ded"
)
EXPECTED_PROFILE_SHA256 = (
    "8c7caacc8a1e9c894eb9b0f2ac6c402c27c328c4924ab7384fede69fa4a527fd"
)
EXPECTED_SEMANTIC_SHA256 = (
    "e07845eaf675a561995164f61e75d9234abbb45ab53cf0ce54b5fb64679ba174"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(SUPPORT) == EXPECTED_SUPPORT_SHA256, "support changed")
spec = spec_from_file_location("k2_z1743_z1749_support", SUPPORT)
require(spec is not None and spec.loader is not None, "cannot load support")
S = module_from_spec(spec)
spec.loader.exec_module(S)
U = S.U
K = S.K
P = S.P
B = S.B
B.START = START
B.END = END
B.LABELS = np.arange(START, B.HORIZON + 1, dtype=np.int64)


def global_scalar_profile(body: tuple[int, ...]):
    return B.profile(body)


def exact_frontier(body: tuple[int, ...], first: int):
    S.FIRST = first
    carrier = U.suffix.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), U.suffix.A.RULER)
    lower = h * U.suffix.ETAS[2]
    period = 14 * lcm(*body)
    first_delta = U.delta(carrier, h, first)
    first_d = period // gcd(period, first)
    amplitudes, signs, ray_digest = S.ray_data(body, carrier, h, period)
    divisors = tuple(d for d in U.support.divisors(period) if d > 1)

    @lru_cache(maxsize=None)
    def top_values(d: int, multiplicity: int):
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (period // d) * direction
            amplitude = amplitudes[residue]
            if amplitude < 0:
                continue
            label = S.first_strictly_after(residue, first, period)
            candidates.extend(
                (
                    amplitude / (label + offset * period),
                    label + offset * period,
                    residue,
                )
                for offset in range(S.SUFFIX_SLOTS)
            )
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        require(len(candidates) >= multiplicity, (body, first, "missing ray", d))
        return tuple(candidates[:multiplicity])

    scalar = []
    for tail in combinations_with_replacement(divisors, S.SUFFIX_SLOTS):
        upper = first_delta
        labels = []
        for d, multiplicity in Counter(tail).items():
            chosen = top_values(d, multiplicity)
            upper += sum((value for value, _label, _residue in chosen), F())
            labels.extend(label for _value, label, _residue in chosen)
        if upper >= lower:
            scalar.append(
                (tuple(sorted((first_d, *tail))), upper, (first, *sorted(labels)))
            )
    require(len({row[0] for row in scalar}) == len(scalar), "duplicate state")

    actual_period, ranges = U.support.safe_cell_ranges(body)
    require(actual_period == period, (body, first, "safe-cell ruler changed"))
    arcs_cache = {}
    crude_kills = []
    crude_survivors = []
    for ds, upper, labels in scalar:
        D = lcm(*ds)
        arcs = arcs_cache.setdefault(D, U.fibre.projected_support_arcs(D, ranges))
        witness = None
        for q in U.support.divisors(D):
            histogram = U.fibre.residue_load_histogram(arcs, q)
            target = max(load for load, count in histogram if count)
            capacity = sum(U.fibre_cap(D, d, q) for d in ds)
            if target > capacity:
                witness = (q, D // q, target, capacity)
                break
        row = (ds, upper, labels, witness)
        (crude_survivors if witness is None else crude_kills).append(row)

    status_kills = []
    states = []
    for ds, upper, labels, _witness in crude_survivors:
        D = lcm(*ds)
        arcs = arcs_cache[D]
        witness = None
        for M in U.support.divisors(D):
            q = D // M
            marginals, capacities = K.hunter_status_data5(D, ds, q)
            histogram = U.fibre.residue_load_histogram(arcs, q)
            feasible, certificate = K.common_status_feasible5(
                q, marginals, capacities, histogram
            )
            if not feasible:
                require(certificate is not None, "missing Farkas witness")
                witness = (
                    q,
                    M,
                    marginals,
                    tuple(sorted(set(capacities))),
                    histogram,
                    certificate,
                )
                break
        row = (ds, upper, labels, witness)
        (states if witness is None else status_kills).append(row)

    stage_rows = tuple(
        tuple(rows) for rows in (scalar, crude_kills, status_kills, states)
    )
    stage_digests = tuple(
        sha256(repr(rows).encode()).hexdigest() for rows in stage_rows
    )
    actual = (
        len(scalar),
        stage_digests[0],
        len(crude_kills),
        stage_digests[1],
        len(status_kills),
        stage_digests[2],
        len(states),
        stage_digests[3],
    )
    require(actual == EXPECTED_STAGES[(body, first)], (body, first, actual))

    packet_summary = None
    cell_count = 0
    if states:
        require((body, first) == TARGET_ROWS[0], "unexpected projected residual")
        cells = P.body_cells(P.A.carrier_for(body), period)
        cell_count = len(cells)
        packet_summary = S.close_literal_packets(
            body, period, amplitudes, cells, lower, first_delta, states
        )
        require(packet_summary[:-1] == EXPECTED_PACKET_SUMMARY, "packet summary changed")
    else:
        require((body, first) == TARGET_ROWS[1], "expected residual vanished")

    return (
        body,
        first,
        h,
        len(carrier),
        period,
        lower,
        first_delta,
        first_d,
        signs,
        ray_digest,
        len(divisors),
        comb(len(divisors) + S.SUFFIX_SLOTS - 1, S.SUFFIX_SLOTS),
        tuple(len(rows) for rows in stage_rows),
        stage_digests,
        cell_count,
        packet_summary,
    )


def exact_profile(row):
    body, first = row
    return exact_frontier(body, first)


def pool_map(function, rows, workers: int, chunksize: int):
    if workers == 1:
        return [function(row) for row in rows]
    with mp.get_context("spawn").Pool(workers) as pool:
        return list(pool.imap(function, rows, chunksize=chunksize))


def render(global_profiles, exact_profiles) -> str:
    global_profiles = tuple(sorted(global_profiles, key=lambda row: row[0]))
    require(len(global_profiles) == comb(14, 6) == 3003, "body universe changed")
    candidate_rows = sum(row[7] for row in global_profiles)
    require(
        candidate_rows == (END - START + 1) * comb(14, 6),
        "candidate-row count changed",
    )
    survivors = tuple(
        sorted(
            (survivor for row in global_profiles for survivor in row[10]),
            key=lambda row: (row[1], row[0]),
        )
    )
    require(survivors == EXPECTED_SURVIVORS, ("global scalar survivors changed", survivors))
    require(
        all(label is not None and kind == "EXACT" for row in survivors for _v, label, kind in row[3]),
        "a surviving scalar suffix is not exact",
    )
    global_profile_hash = sha256(repr(global_profiles).encode()).hexdigest()
    global_survivor_hash = sha256(repr(survivors).encode()).hexdigest()
    if EXPECTED_GLOBAL_PROFILE_SHA256 is not None:
        require(global_profile_hash == EXPECTED_GLOBAL_PROFILE_SHA256, "global profile changed")
    if EXPECTED_GLOBAL_SURVIVOR_SHA256 is not None:
        require(global_survivor_hash == EXPECTED_GLOBAL_SURVIVOR_SHA256, "survivor digest changed")

    exact_profiles = tuple(sorted(exact_profiles, key=lambda row: (row[1], row[0])))
    require(tuple((row[0], row[1]) for row in exact_profiles) == TARGET_ROWS, "exact universe changed")
    packet_profiles = tuple(row for row in exact_profiles if row[15] is not None)
    require(len(packet_profiles) == 1, "projected residual universe changed")
    scalar_states = sum(row[12][0] for row in exact_profiles)
    crude_kills = sum(row[12][1] for row in exact_profiles)
    status_kills = sum(row[12][2] for row in exact_profiles)
    residual_states = sum(row[12][3] for row in exact_profiles)
    literal_packets = sum(row[15][0] for row in packet_profiles)
    minimum_margin = min(row[15][1] for row in packet_profiles)
    maximum_prefix = max(row[15][2] for row in packet_profiles)

    profile_payload = (
        EXPECTED_SUPPORT_SHA256,
        START,
        END,
        global_profile_hash,
        global_survivor_hash,
        exact_profiles,
    )
    profile_hash = sha256(repr(profile_payload).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "complete profile changed")
    semantic_payload = (
        tuple(range(1, 15)),
        6,
        START,
        END,
        5,
        2,
        B.HORIZON,
        S.ALIGNED_TWO_UNION_CAP,
        candidate_rows,
        TARGET_ROWS,
        scalar_states,
        crude_kills,
        status_kills,
        residual_states,
        literal_packets,
        minimum_margin,
        maximum_prefix,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 z1=1743..1749 literal-packet closure",
        f"support_source_sha256={file_sha256(SUPPORT)}",
        (
            "universe=all C(14,6)=3003 bodies;first nonaligned integer labels "
            "1743..1749;four distinct later drifts;two aligned labels"
        ),
        (
            f"global_scalar_horizon={B.HORIZON};omitted_slot_bound=6r/[49(H+1)];"
            f"candidate_rows={candidate_rows};global_scalar_survivors={len(survivors)}"
        ),
        (
            "ray_law=delta(r+mL)=A(r)/(r+mL);"
            "scalar_necessary=sum(delta)>=h/91"
        ),
        (
            "projected_identity=P_(E,Z)=phi_L(C_E minus union_z D_z);"
            "two_aligned_open_union_cap=25/91"
        ),
        (
            f"exact_rows={len(exact_profiles)};scalar_states={scalar_states};"
            f"crude_kills={crude_kills};status_kills={status_kills};"
            f"residual_states={residual_states}"
        ),
        (
            f"all_scalar_eligible_literal_packets={literal_packets};"
            f"projected_kills={literal_packets};survivors=0;"
            f"global_minimum_margin={ftext(minimum_margin)};"
            f"max_prefix_cells={maximum_prefix}"
        ),
    ]
    for row in exact_profiles:
        packet = row[15]
        line = (
            f"EXACT;z1={row[1]};E={','.join(map(str, row[0]))};"
            f"h={ftext(row[2])};r={row[3]};L={row[4]};lower={ftext(row[5])};"
            f"delta1={ftext(row[6])};first_d={row[7]};ray_signs={row[8]};"
            f"ray_sha256={row[9]};denominators={row[10]};trials={row[11]};"
            f"scalar={row[12][0]};crude={row[12][1]};status={row[12][2]};"
            f"residual_states={row[12][3]};stage_sha256={row[13]}"
        )
        if packet is None:
            line += ";conclusion=CRUDE-EMPTY"
        else:
            line += (
                f";body_cells={row[14]};literal_packets={packet[0]};"
                f"projected_kills={packet[0]};min_margin={ftext(packet[1])};"
                f"max_prefix_cells={packet[2]};min_ray_cutoff={ftext(packet[3])};"
                f"max_ray_candidates={packet[4]};max_group_choices={packet[5]};"
                f"minimum_row={packet[6]};direct_control_mass={ftext(packet[7])};"
                f"packet_sha256={packet[8]};state_audit_sha256={packet[9]};"
                "conclusion=PROJECTED-EMPTY"
            )
        lines.append(line)
    lines.extend(
        (
            "independent_control=direct interval subtraction/projection agrees with full cell De Morgan identity on the minimum-margin projected packet",
            "conclusion=the complete projected k=2 first-drift integer band 1743..1749 is empty;no global-cap change claimed",
            f"global_profile_sha256={global_profile_hash}",
            f"global_survivor_sha256={global_survivor_hash}",
            f"profile_sha256={profile_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(4, mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    bodies = tuple(combinations(range(1, 15), 6))
    global_profiles = pool_map(global_scalar_profile, bodies, args.workers, 4)
    exact_profiles = pool_map(exact_profile, TARGET_ROWS, args.workers, 1)
    output = render(global_profiles, exact_profiles)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
