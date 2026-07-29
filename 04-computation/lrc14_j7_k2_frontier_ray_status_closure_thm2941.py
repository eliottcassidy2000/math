#!/usr/bin/env python3
"""Close the unique projected k=2 frontier row for all later labels.

Body E=(1,4,8,10,12,14), first drift 2142, and four later drifts.  The
periodic/antipodal ray law gives exact all-label maxima in each denominator
multiset.  The exact quotient leaves one denominator state.  It survives
the crude all-q capacity screen, but the all-arity common-table transport
of THM-2928, with five-needle Hunter caps, rejects it by an exactly checked
rational Farkas certificate.  No finite label horizon occurs.
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
SUFFIX_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_aligned_projected_arc_suffix_thm2941.out"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_frontier_ray_status_closure_thm2941.out"
)
EXPECTED_UNIFORM_SHA256 = (
    "dfa4788297b8c31fc9b5dce1afadf29d20b267cb4159fa95dadb9346b1980b36"
)
EXPECTED_SUFFIX_OUTPUT_SHA256 = (
    "61e16aab8a368881c574047e576645e6b41837dc9f804f7a78d37230d843612b"
)
EXPECTED_SEMANTIC_SHA256 = (
    "fddb52f665e6c974e3e74839bd377907e21936bdae083a7f79ad9cf1905a63b7"
)
EXPECTED_FRONTIER_ROW = (
    "    FRONTIER;E=1,4,8,10,12,14;h=1049/2940;r=34;L=11760;"
    "analytic_cap=5309;largest_floor=1020;z1=2142;"
    "delta1=821/1049580;"
    "suffix=2172:263/310415,2534:2911/3724980,"
    "3180:279/363580,2396:9419/12327420;"
    "lower=1049/267540;upper=5944239779/1507775985765;"
    "gap=336828043/15680870251956"
)
BODY = (1, 4, 8, 10, 12, 14)
FIRST = 2142
K = 2
SUFFIX_SLOTS = 4
EXPECTED_SCALAR = (
    (
        (196, 280, 840, 980, 2940),
        F(5944239779, 1507775985765),
        (2142, 2172, 2396, 2534, 3180),
    ),
)
EXPECTED_STATUS_WITNESS = (
    840,
    7,
    (0, 120, 120, 0, 0),
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
    file_sha256(SUFFIX_OUTPUT) == EXPECTED_SUFFIX_OUTPUT_SHA256,
    "projected suffix output dependency changed",
)
spec = spec_from_file_location("k2_uniform_support", UNIFORM)
require(spec is not None and spec.loader is not None, "cannot load uniform")
U = module_from_spec(spec)
spec.loader.exec_module(U)


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
        EXPECTED_FRONTIER_ROW in SUFFIX_OUTPUT.read_text(),
        "inherited unique projected k=2 frontier row changed",
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
            F(821, 1049580),
            280,
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
    require(
        tuple(scalar) == EXPECTED_SCALAR,
        ("all-label scalar quotient changed", scalar),
    )
    scalar_digest = sha256(repr(tuple(scalar)).encode()).hexdigest()

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
        not crude_kills
        and tuple(
            (ds, upper, labels)
            for ds, upper, labels, _witness in crude_survivors
        )
        == EXPECTED_SCALAR,
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
    require(
        len(status_kills) == 1
        and not status_survivors
        and status_kills[0][3] == EXPECTED_STATUS_WITNESS,
        ("K5 common-status closure changed", status_kills, status_survivors),
    )

    # A transparent positive/hostile control at the decisive quotient.
    ds, upper, labels, witness = status_kills[0]
    D = lcm(*ds)
    q, M, marginals, capacity_values, histogram, certificate = witness
    control_marginals5, control_capacities5 = hunter_status_data5(D, ds, q)
    require(
        (D, q, M, control_marginals5)
        == (5880, 840, 7, marginals),
        "decisive quotient changed",
    )
    require(
        (
            control_capacities5[0],
            control_capacities5[2],
            control_capacities5[4],
            control_capacities5[6],
        )
        == (3, 7, 7, 7),
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
        "LRC14 projected k=2 z1=2142 all-label ray/status closure",
        f"uniform_ray_status_source_sha256={file_sha256(UNIFORM)}",
        f"suffix_output_sha256={file_sha256(SUFFIX_OUTPUT)}",
        (
            "scope=unique inherited projected k=2 frontier body row;"
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
        f"scalar_state={scalar[0]}",
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
            "common_table_control=240 threshold-four fibres:FEASIBLE;"
            "actual 300 threshold-four fibres:EXACT_FARKAS_INFEASIBLE"
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
        "conclusion=unique projected k=2 z1=2142 frontier row is empty",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ]
    output = "\n".join(lines) + "\n"
    OUTPUT_PATH.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
