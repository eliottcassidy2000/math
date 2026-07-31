#!/usr/bin/env python3
"""Close the three exact-suffix projected k=2 rows at z1=1836.

The global scalar slice leaves five bodies.  Three have L=11760 and exact
suffix envelopes.  For each of these three bodies this verifier reconstructs
all denominator states over every later label by the exact residue-ray law,
then applies the crude all-divisor and common K5 Hunter-status screens.

Fifty-eight states survive those screens.  The scalar slack on each state
makes its remaining literal label universe finite: a label below the worst
top-ray value minus the total slack cannot occur in a viable packet.  The
result is only 84 scalar-eligible packets.  For every packet we compute the
lossless projected residual

    P_(E,Z) = phi_L(C_E minus union_(z in Z) D_z).

Two aligned combs have open union of measure at most 25/91.  A prefix of at
most two body cells already proves mu(P_(E,Z)) >= 25/91 for all 84 packets,
so compact containment is impossible.  No finite label horizon is used.
"""

import argparse
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROJECTED = (
    ROOT
    / "04-computation"
    / "lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py"
)
PRIMARY = ROOT / "04-computation" / "lrc14_j7_k2_z1932_ray_status_closure_thm2941.py"
BAND_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1836_1836_thm2941.out"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1836_exact_body_projected_closure_thm2941.out"
)
EXPECTED_PROJECTED_SHA256 = (
    "76f891edfcc029a08202481304a809e03e8bd81f247afaeabab685825c4d3662"
)
EXPECTED_PRIMARY_SHA256 = (
    "0db812e766eb25ffeed24e7f242edea9bff330514f47b1a0fb459f89cc9c8ead"
)
EXPECTED_BAND_OUTPUT_SHA256 = (
    "fdbeee8d14057edc7b6d1a5a24c660899846aef5eb2a7c6a3004367caf3aa73b"
)
EXPECTED_PROFILE_SHA256 = (
    "956197a829ac189232695e7ffc65dfc36ee6ab3f67b02a05b790cfd1dd5583e9"
)
EXPECTED_SEMANTIC_SHA256 = (
    "f776c4dc7806047dbc146461ab419e2005c5a600a165c31816b832507a4aca08"
)

BODIES = (
    (1, 2, 8, 10, 12, 14),
    (1, 4, 8, 10, 12, 14),
    (2, 4, 8, 10, 12, 14),
)
FIRST = 1836
ALIGNED_TWO_UNION_CAP = F(25, 91)
EXPECTED = {
    BODIES[0]: (
        8,
        "ebc7f30029d142737ac786dcac22a57ef087a243b932e8301e2311ce69ee0c3c",
        1,
        "cc7d09e5dde9a53a4c977c5ca636bb571ba683060cad078ab8b6cc772d378233",
        2,
        "f5adf490aa2d12dbc506736ca27f37c94a6dd78d2cdff877671bd87a5f22a95f",
        5,
        "b04f724b0d97f0f970bd274725d958b53aca490f43f8c5da187308c0449bf048",
        8,
        F(566453, 1302847),
    ),
    BODIES[1]: (
        807,
        "9c008c0d287aea7c3b353a28c6a3dfb5b828d3f1f5fff46193d34ecf896da469",
        178,
        "75c6775879cf8fa6ae7b7a6f1d0eb85564b7d2725557ec1b99b6b74e80d100fc",
        617,
        "d96a4c55611f74be4f82890ff2834fd6c09fd59e43801a39fc6af5b1f8bf9695",
        12,
        "89b65b61c0a0f4463fe32eec67290b158f6b307c0c3c7ac3eb0edc73d4ec28d5",
        12,
        F(1026, 16471),
    ),
    BODIES[2]: (
        72,
        "a76cb3c461ba43b3b59b24d01599c2c468b73c8471ab987c9cb114503e074291",
        1,
        "20435ec408888fda377c2614e00e1e926cc4f16ebcca5489423e6146c2335a9d",
        30,
        "c3a8392a1de5ee23f7054e24d9b2be93706ea686d6c148cb5ee63952cfc66bba",
        41,
        "2cbbc9cf2fb5933d204643c330129de283e6ba0ba7991cda2ea655551e87f4ad",
        64,
        F(4085, 54691),
    ),
}
EXPECTED_SIGNS = {
    BODIES[0]: Counter({-1: 5877, 0: 5, 1: 5877}),
    BODIES[1]: Counter({-1: 5879, 0: 1, 1: 5879}),
    BODIES[2]: Counter({-1: 5879, 0: 1, 1: 5879}),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(PROJECTED) == EXPECTED_PROJECTED_SHA256, "projected dependency changed")
require(file_sha256(PRIMARY) == EXPECTED_PRIMARY_SHA256, "ray/status dependency changed")
require(file_sha256(BAND_OUTPUT) == EXPECTED_BAND_OUTPUT_SHA256, "z1836 scalar slice changed")

PSPEC = spec_from_file_location("z1836_projected_uniform_support", PROJECTED)
require(PSPEC is not None and PSPEC.loader is not None, "cannot load projected support")
P = module_from_spec(PSPEC)
PSPEC.loader.exec_module(P)
KSPEC = spec_from_file_location("z1836_projected_ray_support", PRIMARY)
require(KSPEC is not None and KSPEC.loader is not None, "cannot load ray/status support")
K = module_from_spec(KSPEC)
KSPEC.loader.exec_module(K)
U = K.U


def projected_safe_lower_bound(cells, L, labels, stop_at_cap=True):
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
        if (
            stop_at_cap
            and P.interval_mass(common_danger) <= 1 - ALIGNED_TWO_UNION_CAP
        ):
            break
    return 1 - P.interval_mass(common_danger), cells_used, common_danger


def direct_projected_mass(body, L, labels):
    carrier = tuple(
        (F(left, P.A.RULER), F(right, P.A.RULER))
        for left, right in P.A.carrier_for(body)
    )
    removed = P.merge_fraction(
        [
            interval
            for label in labels
            for interval in P.danger_fraction(label)
        ]
    )
    residual = P.subtract_fraction(carrier, removed)
    projected = []
    for left, right in residual:
        scaled_left = L * left
        scaled_right = L * right
        for integer in range(
            P.floor_fraction(scaled_left),
            P.ceil_fraction(scaled_right),
        ):
            piece_left = max(scaled_left, F(integer)) - integer
            piece_right = min(scaled_right, F(integer + 1)) - integer
            if piece_left < piece_right:
                projected.append((piece_left, piece_right))
    return P.interval_mass(P.merge_fraction(projected))


def scalar_status_frontier(body, carrier, h, lower, L, first_delta, first_d):
    amplitudes = [F(0)]
    signs = Counter()
    for residue in range(1, L):
        amplitude = residue * U.delta(carrier, h, residue)
        require(
            (residue + L) * U.delta(carrier, h, residue + L) == amplitude,
            (body, "ray recurrence failed", residue),
        )
        amplitudes.append(amplitude)
        signs[(amplitude > 0) - (amplitude < 0)] += 1
    require(
        all(
            amplitudes[L - residue] == -amplitudes[residue]
            for residue in range(1, L)
        ),
        (body, "ray antipode failed"),
    )
    require(
        signs == EXPECTED_SIGNS[body],
        (body, "ray sign census changed", signs),
    )
    ray_digest = sha256(repr(tuple(amplitudes)).encode()).hexdigest()
    divisors = tuple(d for d in U.support.divisors(L) if d > 1)

    @lru_cache(maxsize=None)
    def top_values(d, multiplicity):
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (L // d) * direction
            amplitude = amplitudes[residue]
            if amplitude < 0:
                continue
            label = residue
            if label <= FIRST:
                label += ((FIRST + 1 - label + L - 1) // L) * L
            candidates.extend(
                (
                    amplitude / (label + offset * L),
                    label + offset * L,
                    residue,
                )
                for offset in range(multiplicity)
            )
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        require(len(candidates) >= multiplicity, (body, "missing ray", d))
        return tuple(candidates[:multiplicity])

    trials = comb(len(divisors) + 3, 4)
    scalar = []
    for tail in combinations_with_replacement(divisors, 4):
        upper = first_delta
        labels = []
        for d, multiplicity in Counter(tail).items():
            chosen = top_values(d, multiplicity)
            upper += sum((value for value, _label, _residue in chosen), F())
            labels.extend(label for _value, label, _residue in chosen)
        if upper >= lower:
            scalar.append(
                (
                    tuple(sorted((first_d, *tail))),
                    upper,
                    (FIRST, *sorted(labels)),
                )
            )
    require(
        len({ds for ds, _upper, _labels in scalar}) == len(scalar),
        (body, "duplicate denominator state"),
    )

    actual_L, ranges = U.support.safe_cell_ranges(body)
    require(actual_L == L, (body, "safe-cell ruler changed"))
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
    status_survivors = []
    for ds, upper, labels, _crude_witness in crude_survivors:
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
                require(certificate is not None, (body, ds, "missing Farkas witness"))
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
        (status_survivors if witness is None else status_kills).append(row)

    # A HiGHS dual basis is a solver-selected proof witness, not a canonical
    # status instance (MISTAKE-331).  Verify it through the imported engine,
    # but remove only that final coordinate from the replay digest.
    canonical_status_kills = tuple(
        (ds, upper, labels, witness[:-1])
        for ds, upper, labels, witness in status_kills
    )
    require(
        not status_kills
        or sha256(repr(tuple(status_kills)).encode()).hexdigest()
        != sha256(repr(canonical_status_kills).encode()).hexdigest(),
        (body, "raw LP certificate entered status digest"),
    )
    digests = tuple(
        sha256(repr(tuple(rows)).encode()).hexdigest()
        for rows in (scalar, crude_kills, canonical_status_kills, status_survivors)
    )
    expected = EXPECTED[body]
    actual = (
        len(scalar), digests[0], len(crude_kills), digests[1],
        len(status_kills), digests[2], len(status_survivors), digests[3],
    )
    require(
        all(want is None or got == want for got, want in zip(actual, expected[:8])),
        (body, "scalar/status frontier changed", actual, expected[:8]),
    )
    return (
        tuple(amplitudes),
        ray_digest,
        len(divisors),
        trials,
        tuple(scalar),
        tuple(crude_kills),
        tuple(status_kills),
        tuple(status_survivors),
        digests,
    )


def profile(body):
    carrier = U.suffix.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), U.suffix.A.RULER)
    lower = h * U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    require(L == 11760, (body, "unexpected body ruler"))
    first_delta = U.delta(carrier, h, FIRST)
    first_d = L // gcd(L, FIRST)
    require(first_d == 980, (body, "first denominator changed"))
    cells = P.body_cells(P.A.carrier_for(body), L)
    (
        amplitudes,
        ray_digest,
        divisor_count,
        trials,
        scalar,
        crude_kills,
        status_kills,
        states,
        stage_digests,
    ) = scalar_status_frontier(body, carrier, h, lower, L, first_delta, first_d)

    @lru_cache(maxsize=None)
    def ray_heads(d, multiplicity):
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (L // d) * direction
            amplitude = amplitudes[residue]
            if amplitude < 0:
                continue
            label = residue
            if label <= FIRST:
                label += ((FIRST + 1 - label + L - 1) // L) * L
            candidates.extend(
                (amplitude / (label + offset * L), label + offset * L, residue)
                for offset in range(multiplicity)
            )
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        return tuple(candidates[:multiplicity])

    @lru_cache(maxsize=None)
    def eligible_group(d, multiplicity, slack):
        heads = ray_heads(d, multiplicity)
        top_sum = sum((row[0] for row in heads), F())
        threshold = heads[-1][0] - slack
        require(
            threshold > 0,
            (body, "scalar slack did not make ray finite", d, multiplicity, slack),
        )
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (L // d) * direction
            amplitude = amplitudes[residue]
            if amplitude <= 0:
                continue
            label = residue
            if label <= FIRST:
                label += ((FIRST + 1 - label + L - 1) // L) * L
            while amplitude / label >= threshold:
                candidates.append((amplitude / label, label, residue))
                label += L
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        choices = []
        for chosen in combinations(candidates, multiplicity):
            labels = tuple(sorted(row[1] for row in chosen))
            if len(set(labels)) != multiplicity:
                continue
            value = sum((row[0] for row in chosen), F())
            deficit = top_sum - value
            if deficit <= slack:
                choices.append((deficit, labels, value))
        choices.sort()
        require(
            choices and choices[0][0] == 0,
            (body, "eligible group lost its maximizing choice", d, multiplicity),
        )
        return top_sum, threshold, tuple(candidates), tuple(choices)

    packet_count = 0
    killed_count = 0
    survivor_packets = []
    minimum_margin = None
    maximum_cells_used = 0
    minimum_row = None
    state_audit = []
    for state_index, (ds, upper, maximizing_labels, witness) in enumerate(states, start=1):
        require(witness is None, (body, "status survivor carries a witness", ds))
        slack = upper - lower
        tail = list(ds)
        tail.remove(first_d)
        groups = []
        for d, multiplicity in sorted(Counter(tail).items()):
            top_sum, threshold, candidates, choices = eligible_group(
                d, multiplicity, slack
            )
            groups.append(
                (d, multiplicity, top_sum, threshold, len(candidates), choices)
            )
        expected_upper = first_delta + sum((group[2] for group in groups), F())
        require(
            expected_upper == upper,
            (body, "state upper envelope changed", state_index, expected_upper, upper),
        )
        state_packets = state_kills = 0
        state_minimum = None
        for selection in product(*(group[5] for group in groups)):
            total_deficit = sum((choice[0] for choice in selection), F())
            if total_deficit > slack:
                continue
            labels = (FIRST, *sorted(label for choice in selection for label in choice[1]))
            require(
                len(labels) == 5 and len(set(labels)) == 5,
                (body, "literal label packet is not distinct", labels),
            )
            actual_upper = first_delta + sum((choice[2] for choice in selection), F())
            require(
                actual_upper == upper - total_deficit and actual_upper >= lower,
                (body, "scalar-eligible packet arithmetic failed", labels),
            )
            projected_lower, cells_used, _common = projected_safe_lower_bound(
                cells, L, labels
            )
            margin = projected_lower - ALIGNED_TWO_UNION_CAP
            packet_count += 1
            state_packets += 1
            maximum_cells_used = max(maximum_cells_used, cells_used)
            minimum_margin = margin if minimum_margin is None else min(minimum_margin, margin)
            state_minimum = margin if state_minimum is None else min(state_minimum, margin)
            if minimum_row is None or (margin, labels) < (minimum_row[0], minimum_row[1]):
                minimum_row = (margin, labels, projected_lower, cells_used)
            if margin >= 0:
                killed_count += 1
                state_kills += 1
            else:
                survivor_packets.append(
                    (state_index, ds, labels, actual_upper, projected_lower, cells_used)
                )
        require(state_packets > 0, (body, "status state has no literal packet", ds))
        state_audit.append(
            (
                state_index,
                ds,
                upper,
                maximizing_labels,
                slack,
                tuple(
                    (d, multiplicity, top_sum, threshold, candidate_count, len(choices))
                    for d, multiplicity, top_sum, threshold, candidate_count, choices in groups
                ),
                state_packets,
                state_kills,
                state_minimum,
            )
        )
    expected = EXPECTED[body]
    require(
        packet_count == killed_count == expected[8]
        and not survivor_packets
        and minimum_margin == expected[9]
        and maximum_cells_used <= 2,
        (
            body,
            "projected closure changed",
            packet_count,
            killed_count,
            survivor_packets,
            minimum_margin,
            maximum_cells_used,
        ),
    )
    require(minimum_row is not None, (body, "minimum projected row missing"))
    _margin, control_labels, prefix_mass, prefix_cells = minimum_row
    full_cell_mass, full_cells, _full_common = projected_safe_lower_bound(
        cells, L, control_labels, stop_at_cap=False
    )
    direct_mass = direct_projected_mass(body, L, control_labels)
    require(
        full_cells == len(cells)
        and direct_mass == full_cell_mass
        and direct_mass >= prefix_mass,
        (body, "independent projected-residual control failed"),
    )
    state_digest = sha256(repr(tuple(state_audit)).encode()).hexdigest()
    return (
        body,
        h,
        len(carrier),
        L,
        lower,
        first_delta,
        first_d,
        len(cells),
        ray_digest,
        divisor_count,
        trials,
        len(scalar),
        len(crude_kills),
        len(status_kills),
        len(states),
        stage_digests,
        packet_count,
        killed_count,
        minimum_margin,
        maximum_cells_used,
        minimum_row,
        direct_mass,
        state_digest,
        tuple(state_audit),
    )


def render(profiles):
    require(
        tuple(row[0] for row in profiles) == BODIES,
        "exact-body profile universe changed",
    )
    scalar_total = sum(row[11] for row in profiles)
    crude_total = sum(row[12] for row in profiles)
    status_total = sum(row[13] for row in profiles)
    residual_states = sum(row[14] for row in profiles)
    packet_total = sum(row[16] for row in profiles)
    projected_total = sum(row[17] for row in profiles)
    global_margin = min(row[18] for row in profiles)
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(
            profile_hash == EXPECTED_PROFILE_SHA256,
            (
                "profile digest changed",
                profile_hash,
                tuple((row[0], row[15]) for row in profiles),
            ),
        )
    semantic_payload = (
        EXPECTED_PROJECTED_SHA256,
        EXPECTED_PRIMARY_SHA256,
        EXPECTED_BAND_OUTPUT_SHA256,
        BODIES,
        FIRST,
        ALIGNED_TWO_UNION_CAP,
        scalar_total,
        crude_total,
        status_total,
        residual_states,
        packet_total,
        projected_total,
        global_margin,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_hash == EXPECTED_SEMANTIC_SHA256,
            ("projected-closure semantic digest changed", semantic_hash),
        )
    require(
        (
            scalar_total,
            crude_total,
            status_total,
            residual_states,
            packet_total,
            projected_total,
            global_margin,
        )
        == (887, 180, 649, 58, 84, 84, F(1026, 16471)),
        "global exact-body closure totals changed",
    )

    lines = [
        "LRC14 projected k=2 z1=1836 exact-body projected closure",
        f"projected_residual_source_sha256={file_sha256(PROJECTED)}",
        f"ray_status_source_sha256={file_sha256(PRIMARY)}",
        f"scalar_slice_output_sha256={file_sha256(BAND_OUTPUT)}",
        (
            "scope=three L=11760 exact-suffix bodies;all later labels;"
            "no finite label horizon"
        ),
        (
            "ray_law=(z+L)delta(z+L)=zdelta(z);"
            "scalar slack gives a positive per-state ray cutoff"
        ),
        (
            "projected_identity=P_(E,Z)=phi_L(C_E minus union_z D_z);"
            "two-aligned open-union cap=25/91"
        ),
        (
            f"scalar_states={scalar_total};crude_kills={crude_total};"
            f"common_K5_status_kills={status_total};status_survivors={residual_states}"
        ),
        (
            f"all_scalar_eligible_literal_packets={packet_total};"
            f"projected_residual_kills={projected_total};survivors=0;"
            f"global_minimum_margin={ftext(global_margin)}"
        ),
    ]
    for row in profiles:
        (
            body,
            h,
            components,
            L,
            lower,
            first_delta,
            first_d,
            cell_count,
            ray_digest,
            divisor_count,
            trials,
            scalar_count,
            crude_count,
            status_count,
            survivor_count,
            stage_digests,
            packet_count,
            killed_count,
            minimum_margin,
            maximum_cells,
            minimum_row,
            direct_mass,
            state_digest,
            _state_audit,
        ) = row
        lines.append(
            f"BODY;E={','.join(map(str, body))};h={ftext(h)};r={components};"
            f"L={L};lower={ftext(lower)};first={FIRST};"
            f"first_delta={ftext(first_delta)};first_d={first_d};"
            f"body_cells={cell_count};ray_sha256={ray_digest};"
            f"denominators={divisor_count};trials={trials};"
            f"scalar={scalar_count};crude={crude_count};status={status_count};"
            f"residual_states={survivor_count};literal_packets={packet_count};"
            f"projected_kills={killed_count};min_margin={ftext(minimum_margin)};"
            f"max_prefix_cells={maximum_cells};direct_control_mass={ftext(direct_mass)};"
            f"stage_sha256={stage_digests};state_audit_sha256={state_digest};"
            f"minimum_row={minimum_row}"
        )
    lines.extend(
        (
            "independent_control=direct interval subtraction/projection agrees with full cell De Morgan identity",
            "conclusion=the three exact-suffix projected k=2 z1=1836 rows are empty uniformly",
            f"profile_sha256={profile_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(3, mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    if args.workers == 1:
        profiles = [profile(body) for body in BODIES]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            profiles = list(pool.imap(profile, BODIES))
    profiles.sort(key=lambda row: row[0])
    output = render(profiles)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
