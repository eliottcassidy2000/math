#!/usr/bin/env python3
"""Close the projected k=2 z1=1992 row for all later labels.

Body E=(1,4,8,10,12,14), first drift 1992, and four later drifts.  The
periodic/antipodal ray law gives exact all-label maxima in every denominator
multiset.  The exact quotient leaves ten denominator states.  The crude
all-q screen rejects three; the all-arity common-table transport of THM-2928,
with five-needle Hunter caps, rejects the remaining seven by exactly
checked rational Farkas certificates.  No finite label horizon occurs.
The inherited projected suffix is one of the seven status kills.
"""

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement
from math import comb, gcd, lcm
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
UNIFORM = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k3_uniform_ray_status_closure_thm2941.py"
)
PRIOR_ENGINE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_z2004_ray_status_closure_thm2941.py"
)
BAND_OUTPUT = (
    ROOT
    / "05-knowledge" / "results"
    / "lrc14_j7_k2_scalar_band_1992_2003_thm2941.out"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1992_ray_status_closure_thm2941.out"
)
EXPECTED_UNIFORM_SHA256 = (
    "dfa4788297b8c31fc9b5dce1afadf29d20b267cb4159fa95dadb9346b1980b36"
)
EXPECTED_PRIOR_ENGINE_SHA256 = (
    "a5a7246391e7edb8e60ba3d58a72221edd545b70d14104359f20ae06362d9f24"
)
EXPECTED_BAND_OUTPUT_SHA256 = (
    "477e262c5c1246d78787dde0f9e3b9ef71641aae0e9e7f97cf0694f074531888"
)
EXPECTED_SEMANTIC_SHA256 = (
    "1f217116ab2b461e9f18b8e9b710651b8429908d9e7f742322e60b9c405922ee"
)
EXPECTED_FRONTIER_ROW = (
    "z1=1992;delta1=121/170814;"
    "suffix=2060:1867/2119740,2172:263/310415,"
    "2142:821/1049580,2534:2911/3724980;"
    "lower=1049/267540;upper=6496500403/1624087555020;"
    "gap=417952777/5278284553815"
)
BODY = (1, 4, 8, 10, 12, 14)
FIRST = 1992
K = 2
SUFFIX_SLOTS = 4
EXPECTED_FIXED_SCALAR = (
    (280, 490, 588, 840, 980),
    F(6496500403, 1624087555020),
    (1992, 2060, 2142, 2172, 2534),
)
EXPECTED_SCALAR_COUNT = 10
EXPECTED_SCALAR_SHA256 = (
    "e9fc76697bf185118aa3fc26f30c57f189d3e6fbb0c3f60aa819afa929cf1ee3"
)
EXPECTED_CRUDE_KILLS = (
    (
        (196, 490, 588, 980, 980),
        F(223526089697, 56371760586204),
        (1992, 2004, 2060, 2172, 3180),
        (420, 7, 6, 5),
    ),
    (
        (196, 490, 588, 980, 2940),
        F(802278653347, 202195716114588),
        (1992, 2060, 2172, 2396, 3180),
        (420, 7, 6, 5),
    ),
    (
        (490, 588, 980, 980, 2940),
        F(4206940285153, 1061845427394220),
        (1992, 2004, 2060, 2172, 2396),
        (420, 7, 6, 5),
    ),
)
EXPECTED_STATUS_WITNESS = (
    840,
    7,
    (120, 0, 0, 120, 0),
    (3, 4, 5, 6, 7),
    ((0, 120), (3, 420), (4, 112), (5, 168), (6, 20)),
    (
        (4, 5, 6),
        (F(1, 60), F(0), F(0)),
        (F(0), F(1, 60), F(1, 60), F(1, 60), F(1, 60), F(1, 60)),
    ),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(
    file_sha256(UNIFORM) == EXPECTED_UNIFORM_SHA256,
    "uniform ray/status dependency changed",
)
require(
    file_sha256(BAND_OUTPUT) == EXPECTED_BAND_OUTPUT_SHA256,
    "projected k=2 scalar-band output dependency changed",
)
require(
    file_sha256(PRIOR_ENGINE) == EXPECTED_PRIOR_ENGINE_SHA256,
    "prior K5/full32 referee dependency changed",
)
spec = spec_from_file_location("k2_uniform_support", UNIFORM)
require(spec is not None and spec.loader is not None, "cannot load uniform")
U = module_from_spec(spec)
spec.loader.exec_module(U)
zspec = spec_from_file_location("k2_prior_referee_control", PRIOR_ENGINE)
require(zspec is not None and zspec.loader is not None, "cannot load prior referee")
Z = module_from_spec(zspec)
zspec.loader.exec_module(Z)


def spanning_trees(vertex_count):
    edges = tuple(combinations(range(vertex_count), 2))
    result = []
    for chosen in combinations(edges, vertex_count - 1):
        seen = {0}
        changed = True
        while changed:
            changed = False
            for left, right in chosen:
                if left in seen and right not in seen:
                    seen.add(right)
                    changed = True
                elif right in seen and left not in seen:
                    seen.add(left)
                    changed = True
        if len(seen) == vertex_count:
            result.append(chosen)
    require(
        len(result) == vertex_count ** (vertex_count - 2),
        ("Cayley tree count changed", vertex_count, len(result)),
    )
    return tuple(result)


TREES5 = spanning_trees(5)


@lru_cache(maxsize=None)
def hunter_cap5(M, ds, es):
    sizes = tuple((M // d) * e for d, e in zip(ds, es))
    overlaps = {
        (i, j): U.pair_lower(M, ds[i], es[i], ds[j], es[j])
        for i, j in combinations(range(5), 2)
    }
    invoice = max(
        sum(overlaps[edge] for edge in tree) for tree in TREES5
    )
    return min(M, sum(sizes) - invoice)


def hunter_status_data5(D, ds, q):
    M = D // q
    inner_ds = []
    lows = []
    marginals = []
    for d in ds:
        common = gcd(d, q)
        low, remainder = divmod((d + 6) // 7, common)
        inner_ds.append(d // common)
        lows.append(low)
        marginals.append((q // common) * remainder)
    capacities = tuple(
        hunter_cap5(
            M,
            tuple(inner_ds),
            tuple(
                lows[index] + ((pattern >> index) & 1)
                for index in range(5)
            ),
        )
        for pattern in range(32)
    )
    return tuple(marginals), capacities


@lru_cache(maxsize=None)
def common_status_feasible5(q, marginals, capacities, histogram):
    """The single-table all-threshold condition (37tg+), with exact kills."""
    tail_rows = []
    tail_rhs = []
    thresholds = []
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = sum(
            count for load, count in histogram if load >= threshold
        )
        good = tuple(
            int(capacity >= threshold) for capacity in capacities
        )
        if all(good):
            continue
        thresholds.append(threshold)
        if not any(good):
            return False, ((threshold,), (F(1),), (F(0),) * 6)
        tail_rows.append(good)
        tail_rhs.append(demand)

    equality_rows = [
        (1,) * 32,
        *[
            tuple((pattern >> index) & 1 for pattern in range(32))
            for index in range(5)
        ],
    ]
    equality_rhs = (q, *marginals)
    if not tail_rows:
        return True, None
    primal = linprog(
        np.zeros(32),
        A_ub=-np.asarray(tail_rows, dtype=float),
        b_ub=-np.asarray(tail_rhs, dtype=float),
        A_eq=np.asarray(equality_rows, dtype=float),
        b_eq=np.asarray(equality_rhs, dtype=float),
        bounds=(0, None),
        method="highs",
    )
    if primal.success:
        return True, None
    require(primal.status == 2, ("unexpected primal status", primal.status))

    tail_count = len(tail_rows)
    dual_rows = []
    dual_rhs = []
    for pattern in range(32):
        dual_rows.append(
            tuple(tail_rows[row][pattern] for row in range(tail_count))
            + tuple(-equality_rows[row][pattern] for row in range(6))
        )
        dual_rhs.append(0)
    dual_rows.append(
        tuple(-value for value in tail_rhs) + tuple(equality_rhs)
    )
    dual_rhs.append(-1)
    dual = linprog(
        np.zeros(tail_count + 6),
        A_ub=np.asarray(dual_rows, dtype=float),
        b_ub=np.asarray(dual_rhs, dtype=float),
        bounds=[(0, None)] * tail_count + [(None, None)] * 6,
        method="highs",
    )
    if not dual.success:
        return True, None
    alpha = tuple(
        F(float(value)).limit_denominator(1_000_000)
        for value in dual.x[:tail_count]
    )
    z = tuple(
        F(float(value)).limit_denominator(1_000_000)
        for value in dual.x[tail_count:]
    )
    slacks = tuple(
        sum(z[row] * equality_rows[row][pattern] for row in range(6))
        - sum(
            alpha[row] * tail_rows[row][pattern]
            for row in range(tail_count)
        )
        for pattern in range(32)
    )
    contradiction = (
        sum(z[row] * equality_rhs[row] for row in range(6))
        - sum(alpha[row] * tail_rhs[row] for row in range(tail_count))
    )
    if (
        all(value >= 0 for value in alpha)
        and all(value >= 0 for value in slacks)
        and contradiction < 0
    ):
        return False, (tuple(thresholds), alpha, z)
    return True, None


def main():
    require(
        EXPECTED_FRONTIER_ROW in BAND_OUTPUT.read_text(),
        "inherited projected k=2 z1=1992 row changed",
    )
    pair_rows, control_marginals, control_caps, control_certificate = (
        U.controls()
    )
    require(len(TREES5) == 125, "K5 spanning-tree count changed")
    carrier = U.suffix.A.carrier_for(BODY)
    h = F(
        sum(right - left for left, right in carrier),
        U.suffix.A.RULER,
    )
    L = 14 * lcm(*BODY)
    lower = h * U.suffix.ETAS[K]
    first_delta = U.delta(carrier, h, FIRST)
    first_d = L // gcd(L, FIRST)
    wall = U.suffix.PROJECTED_RATIOS[K] * L
    high_floor = max(15, wall.numerator // wall.denominator + 1)
    require(
        (
            h,
            len(carrier),
            L,
            lower,
            first_delta,
            first_d,
            high_floor,
        )
        == (
            F(1049, 2940),
            34,
            11760,
            F(1049, 267540),
            F(121, 170814),
            490,
            1020,
        ),
        "frontier geometry changed",
    )
    require(FIRST >= high_floor, "unexpected high-label obligation")

    amplitudes = [F(0)]
    signs = Counter()
    for residue in range(1, L):
        amplitude = residue * U.delta(carrier, h, residue)
        require(
            (residue + L)
            * U.delta(carrier, h, residue + L)
            == amplitude,
            ("ray recurrence failed", residue),
        )
        amplitudes.append(amplitude)
        signs[(amplitude > 0) - (amplitude < 0)] += 1
    require(
        all(
            amplitudes[L - residue] == -amplitudes[residue]
            for residue in range(1, L)
        ),
        "ray antipode failed",
    )
    require(
        signs == Counter({-1: 5879, 0: 1, 1: 5879}),
        ("ray sign counts changed", signs),
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
            first = residue
            if first <= FIRST:
                first += ((FIRST + 1 - first + L - 1) // L) * L
            candidates.extend(
                (
                    amplitude / (first + offset * L),
                    first + offset * L,
                    residue,
                )
                for offset in range(multiplicity)
            )
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        require(
            len(candidates) >= multiplicity,
            ("missing nonnegative ray", d),
        )
        return tuple(candidates[:multiplicity])

    trials = comb(len(divisors) + SUFFIX_SLOTS - 1, SUFFIX_SLOTS)
    scalar = []
    for tail in combinations_with_replacement(divisors, SUFFIX_SLOTS):
        upper = first_delta
        labels = []
        for d, multiplicity in Counter(tail).items():
            chosen = top_values(d, multiplicity)
            upper += sum((value for value, _label, _residue in chosen), F(0))
            labels.extend(label for _value, label, _residue in chosen)
        if upper >= lower:
            scalar.append(
                (
                    tuple(sorted((first_d, *tail))),
                    upper,
                    (FIRST, *sorted(labels)),
                )
            )
    require(len({ds for ds, _upper, _labels in scalar}) == len(scalar), "duplicate state")
    scalar_digest = sha256(repr(tuple(scalar)).encode()).hexdigest()
    require(
        len(scalar) == EXPECTED_SCALAR_COUNT
        and scalar_digest == EXPECTED_SCALAR_SHA256
        and EXPECTED_FIXED_SCALAR in scalar,
        ("all-label scalar quotient changed", scalar_digest, scalar),
    )

    actual_L, ranges = U.support.safe_cell_ranges(BODY)
    require(actual_L == L, "safe-cell ruler changed")
    arcs_cache = {}
    crude_kills = []
    crude_survivors = []
    for ds, upper, labels in scalar:
        D = lcm(*ds)
        arcs = arcs_cache.setdefault(
            D, U.fibre.projected_support_arcs(D, ranges)
        )
        witness = None
        for q in U.support.divisors(D):
            histogram = U.fibre.residue_load_histogram(arcs, q)
            target = max(load for load, count in histogram if count)
            capacity = sum(U.fibre_cap(D, d, q) for d in ds)
            if target > capacity:
                witness = (q, D // q, target, capacity)
                break
        row = (ds, upper, labels, witness)
        if witness is None:
            crude_survivors.append(row)
        else:
            crude_kills.append(row)
    require(
        tuple(crude_kills) == EXPECTED_CRUDE_KILLS
        and len(crude_survivors) == EXPECTED_SCALAR_COUNT - len(EXPECTED_CRUDE_KILLS),
        ("crude all-q stage changed", crude_kills, crude_survivors),
    )

    status_kills = []
    status_survivors = []
    for ds, upper, labels, _crude_witness in crude_survivors:
        D = lcm(*ds)
        arcs = arcs_cache[D]
        witness = None
        for M in U.support.divisors(D):
            q = D // M
            marginals, capacities = hunter_status_data5(D, ds, q)
            histogram = U.fibre.residue_load_histogram(arcs, q)
            feasible, certificate = common_status_feasible5(
                q, marginals, capacities, histogram
            )
            if not feasible:
                require(certificate is not None, "status kill lacks certificate")
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
        if witness is None:
            status_survivors.append(row)
        else:
            status_kills.append(row)
    fixed_row = next(
        (row for row in status_kills if row[:3] == EXPECTED_FIXED_SCALAR),
        None,
    )
    require(
        len(status_kills) == EXPECTED_SCALAR_COUNT - len(EXPECTED_CRUDE_KILLS)
        and not status_survivors
        and fixed_row is not None
        and fixed_row[3] == EXPECTED_STATUS_WITNESS,
        ("K5 common-status closure changed", status_kills, status_survivors),
    )

    # A transparent positive/hostile control at the decisive quotient.
    require(fixed_row is not None, "fixed projected state disappeared")
    ds, upper, labels, witness = fixed_row
    D = lcm(*ds)
    q, M, marginals, capacity_values, histogram, certificate = witness
    control_marginals5, control_capacities5 = hunter_status_data5(D, ds, q)
    z_marginals5, z_capacities5 = Z.hunter_status_data5(D, ds, q)
    z_feasible5, z_certificate5 = Z.common_status_feasible5(
        q,
        z_marginals5,
        z_capacities5,
        histogram,
    )
    require(
        (D, q, M, control_marginals5)
        == (5880, 840, 7, marginals),
        "decisive quotient changed",
    )
    require(
        TREES5 == Z.TREES5
        and (z_marginals5, z_capacities5)
        == (control_marginals5, control_capacities5)
        and not z_feasible5
        and z_certificate5 == certificate,
        "prior referee engine disagrees on the fixed z1992 state",
    )
    require(
        (
            control_capacities5[0],
            control_capacities5[2],
            control_capacities5[6],
            control_capacities5[22],
            control_capacities5[1],
        )
        == (3, 4, 5, 6, 7),
        "decisive status capacities changed",
    )
    feasible_control, _ = common_status_feasible5(
        q,
        marginals,
        control_capacities5,
        ((0, 600), (4, 240)),
    )
    require(feasible_control, "240-heavy-fibre positive control rejected")
    heavy_four = sum(count for load, count in histogram if load >= 4)
    require(
        heavy_four == 300
        and sum(marginals) == 240
        and heavy_four > sum(marginals),
        "transparent threshold-four obstruction changed",
    )

    semantic_payload = (
        BODY,
        FIRST,
        h,
        len(carrier),
        L,
        lower,
        first_delta,
        first_d,
        high_floor,
        tuple(sorted(signs.items())),
        ray_digest,
        len(divisors),
        trials,
        tuple(scalar),
        tuple(crude_kills),
        tuple(crude_survivors),
        tuple(status_kills),
        tuple(status_survivors),
        len(TREES5),
        pair_rows,
        control_marginals,
        control_caps,
        control_certificate,
        heavy_four,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_hash == EXPECTED_SEMANTIC_SHA256,
            "k=2 frontier semantic digest changed",
        )

    lines = [
        "LRC14 projected k=2 z1=1992 all-label ray/status closure",
        f"uniform_ray_status_source_sha256={file_sha256(UNIFORM)}",
        f"prior_referee_source_sha256={file_sha256(PRIOR_ENGINE)}",
        f"scalar_band_output_sha256={file_sha256(BAND_OUTPUT)}",
        (
            "scope=inherited projected k=2 z1=1992 body row;"
            "four distinct later nonaligned labels;no finite horizon"
        ),
        (
            f"body={BODY};h={ftext(h)};r={len(carrier)};L={L};"
            f"first={FIRST};first_d={first_d};high_floor={high_floor};"
            f"lower={ftext(lower)};first_delta={ftext(first_delta)}"
        ),
        (
            "ray_law=(z+L)delta(z+L)=zdelta(z);"
            "A(L-b)=-A(b);all denominator maxima attained"
        ),
        (
            f"ray_checks={L-1};ray_signs={dict(sorted(signs.items()))};"
            f"ray_amplitude_sha256={ray_digest}"
        ),
        (
            f"nontrivial_denominators={len(divisors)};"
            f"denominator_trials={trials};all_n_scalar_states={len(scalar)};"
            f"scalar_state_sha256={scalar_digest}"
        ),
        f"fixed_projected_scalar_state={EXPECTED_FIXED_SCALAR}",
        (
            f"crude_all_q_kills={len(crude_kills)};"
            f"crude_all_q_survivors={len(crude_survivors)}"
        ),
        (
            f"K5_spanning_trees={len(TREES5)};"
            f"pair_overlap_exhaustive_controls={pair_rows};"
            "status_cells=32;equality_rows=6"
        ),
        (
            "fixed_state_transparent_obstruction=300 load>=4 fibres;"
            "only 240 exceptional-status events;EXACT_FARKAS_INFEASIBLE"
        ),
        (
            f"common_K5_status_kills={len(status_kills)};"
            f"common_K5_status_survivors={len(status_survivors)}"
        ),
        (
            f"decisive_witness=D={D};q={q};M={M};marginals={marginals};"
            f"capacity_values={capacity_values};histogram={histogram};"
            f"farkas={certificate}"
        ),
        (
            f"transparent_y4_obstruction=heavy_fibres:{heavy_four};"
            f"available_exception_events:{sum(marginals)}"
        ),
        "conclusion=projected k=2 z1=1992 row is empty uniformly",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ]
    lines[10:10] = [
        *[
            (
                f"state[{index}]={crude_ds};upper={ftext(crude_upper)};"
                f"labels={crude_labels};kill=CRUDE;"
                f"witness_q_M_target_capacity={crude_witness}"
            )
            for index, (crude_ds, crude_upper, crude_labels, crude_witness)
            in enumerate(crude_kills, start=1)
        ],
        *[
            (
                f"state[{index}]={state_ds};upper={ftext(state_upper)};"
                f"labels={state_labels};kill=COMMON_STATUS;"
                f"kill_q={state_witness[0]};kill_M={state_witness[1]};"
                f"marginals={state_witness[2]};farkas={state_witness[5]}"
            )
            for index, (state_ds, state_upper, state_labels, state_witness)
            in enumerate(status_kills, start=1 + len(crude_kills))
        ],
    ]
    output = "\n".join(lines) + "\n"
    OUTPUT_PATH.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
