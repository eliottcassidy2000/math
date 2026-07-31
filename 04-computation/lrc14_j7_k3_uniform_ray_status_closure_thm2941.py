#!/usr/bin/env python3
"""Uniform z1=250 closure from periodic drift rays and CRT status transport.

Fix the literal body E=(1,4,8,10,12,14), its endpoint ruler L=11760,
carrier C, carrier mass h, and singleton danger coverage c(z).  Every
carrier endpoint is a multiple of R/L, where R is the integration ruler.
The danger primitive satisfies P(x+mR)=P(x)+mR/7.  Consequently, for
0<b<L and a>=0,

    (La+b)(c(La+b)-h/7) = b(c(b)-h/7).                 (RAY)

For b=0, c(La)=h/7.  If d=L/gcd(L,b), then
b=(L/d)u with u a unit modulo d.  Thus denominator d has phi(d) rays.
Reversal gives A_{L-b}=-A_b while preserving the denominator class.
Consequently every class has either a positive ray or an attained zero
ray, and the exact maximum over m distinct labels above a cutoff is the
finite merge of the first m points on each nonnegative ray.

For z1=250, d1=1176.  We enumerate every multiset consisting of d1 and
three arbitrary nontrivial divisors of L, apply the ray-supremum excess
condition, then the crude all-divisor fibre capacity.  The residual is
closed by a four-needle common-status transport.  Within each outer fibre,
sharp CRT-forced pair intersections are inserted into Hunter's maximum
spanning-tree union cap.  One real 16-cell status table must dominate all
target-load tails simultaneously.  Every infeasibility returned here is
independently checked by an exact rational Farkas certificate.

The order of later labels and the projected high-wall condition are
relaxed, so their omission cannot create a false obstruction.  There is
no finite label horizon anywhere in this computation.
"""

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache, reduce
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement
from math import gcd, lcm
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parents[1]
SUFFIX_PATH = (
    ROOT
    / "04-computation"
    / "lrc14_j7_aligned_projected_arc_suffix_thm2941.py"
)
FIBRE_PATH = (
    ROOT
    / "04-computation"
    / "lrc14_three_drift_body_projection_fiber_thm2928.py"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k3_uniform_ray_status_closure_thm2941.out"
)
EXPECTED_SUFFIX_SHA256 = (
    "a003d287f618eb301edf6974d0b67dc128c4f380a169e7809ed5b5754e8b8303"
)
EXPECTED_FIBRE_SHA256 = (
    "42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4"
)
EXPECTED_SEMANTIC_SHA256 = (
    "c3be93a08ba7315136834ccb8dc15ce68421d3b4f102884ac6383bb0768b3948"
)

BODY = (1, 4, 8, 10, 12, 14)
FIRST = 250
K = 3


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


require(
    file_sha256(SUFFIX_PATH) == EXPECTED_SUFFIX_SHA256,
    "THM-2941 suffix dependency changed",
)
require(
    file_sha256(FIBRE_PATH) == EXPECTED_FIBRE_SHA256,
    "THM-2928 fibre dependency changed",
)
suffix = load_module("uniform_ray_suffix", SUFFIX_PATH)
fibre = load_module("uniform_ray_fibre", FIBRE_PATH)
support = fibre.support_module


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def delta(carrier, h, label):
    return suffix.A.singleton_coverage(carrier, label) - h / 7


def spanning_trees():
    edges = tuple(combinations(range(4), 2))
    result = []
    for chosen in combinations(edges, 3):
        vertices = {vertex for edge in chosen for vertex in edge}
        if len(vertices) == 4:
            result.append(chosen)
    require(len(result) == 16, "K4 spanning-tree count changed")
    return tuple(result)


TREES = spanning_trees()


def pair_lower(M, d, e, f, z):
    """Sharp phase-free minimum intersection of two lifted unit needles."""
    common = gcd(d, f)
    ae, re = divmod(e, common)
    az, rz = divmod(z, common)
    base = common * ae * az + ae * rz + az * re
    return (M // lcm(d, f)) * (
        base + max(0, re + rz - common)
    )


@lru_cache(maxsize=None)
def hunter_cap(M, ds, es):
    """Four-needle union cap with the largest forced-overlap tree invoice."""
    sizes = tuple((M // d) * e for d, e in zip(ds, es))
    overlaps = {
        (i, j): pair_lower(M, ds[i], es[i], ds[j], es[j])
        for i, j in combinations(range(4), 2)
    }
    invoice = max(
        sum(overlaps[edge] for edge in tree) for tree in TREES
    )
    return min(M, sum(sizes) - invoice)


def fibre_cap(D, d, q):
    ell = (d + 6) // 7
    common = gcd(d, q)
    height = D // lcm(d, q)
    return height * ((ell + common - 1) // common)


def hunter_status_data(D, ds, q):
    """Marginals and 16 Hunter trace caps at outer modulus q."""
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
        hunter_cap(
            M,
            tuple(inner_ds),
            tuple(
                lows[index] + ((pattern >> index) & 1)
                for index in range(4)
            ),
        )
        for pattern in range(16)
    )
    return tuple(marginals), capacities


@lru_cache(maxsize=None)
def common_status_feasible(q, marginals, capacities, histogram):
    """Common-tail LP; false is backed by an exact rational Farkas witness."""
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
            # B=0, alpha=1, z=0 gives contradiction -demand<0.
            return False, ((threshold,), (F(1),), (F(0),) * 5)
        tail_rows.append(good)
        tail_rhs.append(demand)

    equality_rows = [
        (1,) * 16,
        *[
            tuple((pattern >> index) & 1 for pattern in range(16))
            for index in range(4)
        ],
    ]
    equality_rhs = (q, *marginals)
    if not tail_rows:
        return True, None

    primal = linprog(
        np.zeros(16),
        A_ub=-np.asarray(tail_rows, dtype=float),
        b_ub=-np.asarray(tail_rhs, dtype=float),
        A_eq=np.asarray(equality_rows, dtype=float),
        b_eq=np.asarray(equality_rhs, dtype=float),
        bounds=(0, None),
        method="highs",
    )
    if primal.success:
        # A failed rational reconstruction only weakens the screen.
        rational = tuple(
            F(float(value)).limit_denominator(1_000_000)
            for value in primal.x
        )
        verified = (
            all(value >= 0 for value in rational)
            and all(
                sum(
                    coefficient * value
                    for coefficient, value in zip(row, rational)
                )
                == rhs
                for row, rhs in zip(equality_rows, equality_rhs)
            )
            and all(
                sum(
                    coefficient * value
                    for coefficient, value in zip(row, rational)
                )
                >= rhs
                for row, rhs in zip(tail_rows, tail_rhs)
            )
        )
        return True, rational if verified else None

    require(primal.status == 2, ("unexpected primal status", primal.status))

    # Farkas: A^T z-B^T alpha>=0 and b.z-H.alpha<0.
    tail_count = len(tail_rows)
    dual_rows = []
    dual_rhs = []
    for pattern in range(16):
        dual_rows.append(
            tuple(tail_rows[row][pattern] for row in range(tail_count))
            + tuple(-equality_rows[row][pattern] for row in range(5))
        )
        dual_rhs.append(0)
    dual_rows.append(
        tuple(-value for value in tail_rhs) + tuple(equality_rhs)
    )
    dual_rhs.append(-1)
    dual = linprog(
        np.zeros(tail_count + 5),
        A_ub=np.asarray(dual_rows, dtype=float),
        b_ub=np.asarray(dual_rhs, dtype=float),
        bounds=[(0, None)] * tail_count + [(None, None)] * 5,
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
    column_slacks = tuple(
        sum(z[row] * equality_rows[row][pattern] for row in range(5))
        - sum(
            alpha[row] * tail_rows[row][pattern]
            for row in range(tail_count)
        )
        for pattern in range(16)
    )
    contradiction = (
        sum(z[row] * equality_rhs[row] for row in range(5))
        - sum(alpha[row] * tail_rhs[row] for row in range(tail_count))
    )
    if (
        all(value >= 0 for value in alpha)
        and all(value >= 0 for value in column_slacks)
        and contradiction < 0
    ):
        return False, (tuple(thresholds), alpha, z)
    return True, None


@lru_cache(maxsize=None)
def needle_family(M, d, e):
    """All lifted unit needles, used only by exhaustive positive controls."""
    require(M % d == 0 and 0 <= e <= d, ("bad needle family", M, d, e))
    masks = set()
    for step in range(d):
        if gcd(step, d) != 1:
            continue
        for phase in range(d):
            residues = {
                (phase + step * index) % d for index in range(e)
            }
            mask = 0
            for residue in residues:
                for lift in range(M // d):
                    mask |= 1 << (residue + d * lift)
            masks.add(mask)
    require(masks, ("empty needle family", M, d, e))
    return tuple(masks)


def controls():
    """Independent small exact controls for overlap and common-table logic."""
    pair_rows = 0
    for M in range(1, 9):
        divisors = support.divisors(M)
        for d in divisors:
            for f in divisors:
                for e in range(d + 1):
                    for z in range(f + 1):
                        actual = min(
                            (left & right).bit_count()
                            for left in needle_family(M, d, e)
                            for right in needle_family(M, f, z)
                        )
                        require(
                            actual == pair_lower(M, d, e, f, z),
                            ("pair lower control failed", M, d, e, f, z),
                        )
                        pair_rows += 1

    # Four-needle strict common-table control at D=28, q=2.
    marginals, capacities = hunter_status_data(
        28,
        (4, 4, 28, 28),
        2,
    )
    feasible_histogram = ((9, 2),)
    feasible, primal = common_status_feasible(
        2,
        marginals,
        capacities,
        feasible_histogram,
    )
    require(feasible and primal is not None, "coverable status control rejected")
    hostile_histogram = ((8, 1), (14, 1))
    hostile, certificate = common_status_feasible(
        2,
        marginals,
        capacities,
        hostile_histogram,
    )
    require(
        not hostile and certificate is not None,
        "incompatible-threshold control survived",
    )
    # The alpha/z dual representative is noncanonical.  Return only the
    # deterministic control instances and their active threshold set after
    # exact certificate verification.
    control_instances = (
        feasible_histogram,
        hostile_histogram,
        certificate[0],
    )
    return pair_rows, marginals, capacities, control_instances


def main():
    pair_rows, control_marginals, control_caps, control_instances = controls()
    carrier = suffix.A.carrier_for(BODY)
    h = F(
        sum(right - left for left, right in carrier),
        suffix.A.RULER,
    )
    L = 14 * lcm(*BODY)
    cell = suffix.A.RULER // L
    require(
        L == 11760
        and suffix.A.RULER % L == 0
        and all(
            endpoint % cell == 0
            for interval in carrier
            for endpoint in interval
        ),
        "carrier is not on the endpoint grid",
    )
    lower = h * suffix.ETAS[K]
    first_delta = delta(carrier, h, FIRST)
    require(
        (h, lower, first_delta)
        == (F(1049, 2940), F(1049, 89180), F(484, 128625)),
        "z1=250 constants changed",
    )

    # Exhaustive recurrence and reversal controls on one complete period.
    ray_amplitudes = [F(0)]
    for residue in range(1, L):
        amplitude = residue * delta(carrier, h, residue)
        require(
            (residue + L) * delta(carrier, h, residue + L)
            == amplitude,
            ("periodic ray law failed", residue),
        )
        ray_amplitudes.append(amplitude)
    require(
        all(
            ray_amplitudes[L - residue] == -ray_amplitudes[residue]
            for residue in range(1, L)
        ),
        "ray reversal antisymmetry failed",
    )
    require(delta(carrier, h, L) == delta(carrier, h, 2 * L) == 0, "zero ray")
    require(
        ray_amplitudes[260] == F(1483, 1029)
        and delta(carrier, h, 12020) == F(1483, 12368580),
        "positive ray control changed",
    )
    require(
        ray_amplitudes[980] == F(-2, 3)
        and delta(carrier, h, 12740) == F(-1, 19110),
        "negative ray control changed",
    )
    require(
        ray_amplitudes[5880] == 0
        and delta(carrier, h, 17640) == 0,
        "zero denominator-two ray changed",
    )
    ray_digest = sha256(repr(tuple(ray_amplitudes)).encode()).hexdigest()

    divisors = tuple(d for d in support.divisors(L) if d > 1)

    @lru_cache(maxsize=None)
    def top_ray_values(d, multiplicity):
        candidates = []
        direction_count = 0
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            direction_count += 1
            residue = (L // d) * direction
            amplitude = ray_amplitudes[residue]
            if amplitude < 0:
                continue
            first = residue if residue > FIRST else residue + L
            candidates.extend(
                (
                    amplitude / (first + step * L),
                    first + step * L,
                    residue,
                )
                for step in range(multiplicity)
            )
        require(
            direction_count
            == sum(1 for direction in range(1, d) if gcd(direction, d) == 1),
            "direction count changed",
        )
        # Reversal preserves d and negates amplitude.  Thus each class has
        # a positive direction or an attained zero direction, providing
        # at least ``multiplicity`` literal candidates.
        require(
            len(candidates) >= multiplicity,
            ("missing attained nonnegative ray", d, multiplicity),
        )
        candidates.sort(
            key=lambda row: (
                -row[0],
                L + 1 if row[1] is None else row[1],
            )
        )
        return tuple(candidates[:multiplicity])

    shapes = tuple(
        sorted(
            {
                tuple(sorted((1176, *tail)))
                for tail in combinations_with_replacement(divisors, 3)
            }
        )
    )
    require(len(divisors) == 59 and len(shapes) == 35990, "shape universe changed")

    excess_survivors = []
    excess_digest = sha256()
    for ds in shapes:
        remaining = list(ds)
        remaining.remove(1176)
        counts = Counter(remaining)
        top = tuple(
            (d, top_ray_values(d, multiplicity))
            for d, multiplicity in sorted(counts.items())
        )
        upper = first_delta + sum(
            (
                sum((value for value, _label, _residue in rows), F(0))
                for _d, rows in top
            ),
            F(0),
        )
        if upper >= lower:
            excess_survivors.append(ds)
            excess_digest.update(f"{ds}|{upper}\n".encode())
    require(
        len(excess_survivors) == 1965,
        ("ray-excess survivor count changed", len(excess_survivors)),
    )

    actual_L, ranges = support.safe_cell_ranges(BODY)
    require(actual_L == L, "safe-cell ruler changed")
    body_by_D = {}

    def arcs_for(D):
        if D not in body_by_D:
            body_by_D[D] = fibre.projected_support_arcs(D, ranges)
        return body_by_D[D]

    crude_survivors = []
    crude_kills = 0
    for ds in excess_survivors:
        D = lcm(*ds)
        arcs = arcs_for(D)
        killed = False
        for q in support.divisors(D):
            histogram = fibre.residue_load_histogram(arcs, q)
            maximum = max(load for load, count in histogram if count)
            capacity = sum(fibre_cap(D, d, q) for d in ds)
            if maximum > capacity:
                killed = True
                break
        if killed:
            crude_kills += 1
        else:
            crude_survivors.append(ds)
    require(
        (crude_kills, len(crude_survivors)) == (699, 1266),
        ("crude stage changed", crude_kills, len(crude_survivors)),
    )

    status_kills = []
    status_survivors = []
    kills_by_M = Counter()
    verified_instance_digest = sha256()
    for ds in crude_survivors:
        D = lcm(*ds)
        arcs = arcs_for(D)
        witness = None
        for M in support.divisors(D):
            q = D // M
            marginals, capacities = hunter_status_data(D, ds, q)
            histogram = fibre.residue_load_histogram(arcs, q)
            feasible, certificate = common_status_feasible(
                q,
                marginals,
                capacities,
                histogram,
            )
            if not feasible:
                require(certificate is not None, "kill lacks exact certificate")
                witness = (
                    q,
                    M,
                    marginals,
                    tuple(sorted(set(capacities))),
                    histogram,
                    certificate,
                )
                break
        if witness is None:
            status_survivors.append(ds)
        else:
            status_kills.append((ds, D, witness))
            kills_by_M[witness[1]] += 1
            # The exact certificate has already been verified above.  HiGHS
            # may return a different valid dual basis on another replay, so
            # freeze the deterministic infeasible instance, not the chosen
            # certificate representative.
            verified_instance_digest.update(
                f"{ds}|{D}|{witness[:-1]}\n".encode()
            )
    require(
        len(status_kills) == 1266 and not status_survivors,
        ("status closure changed", len(status_kills), status_survivors),
    )
    require(
        kills_by_M == Counter({3: 44, 4: 128, 5: 743, 6: 59, 7: 292}),
        ("decisive modulus histogram changed", kills_by_M),
    )

    semantic_payload = (
        BODY,
        L,
        h,
        FIRST,
        first_delta,
        lower,
        len(shapes),
        len(excess_survivors),
        crude_kills,
        len(crude_survivors),
        tuple(sorted(kills_by_M.items())),
        ray_digest,
        excess_digest.hexdigest(),
        verified_instance_digest.hexdigest(),
        (
            status_kills[0][0],
            status_kills[0][1],
            status_kills[0][2][:-1],
        ),
        pair_rows,
        control_marginals,
        control_caps,
        control_instances,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_hash == EXPECTED_SEMANTIC_SHA256,
            "uniform closure semantic digest changed",
        )

    first_shape, first_D, first_witness = status_kills[0]
    q, M, marginals, capacities, histogram, _certificate = first_witness
    lines = [
        "LRC14 k=3 z1=250 uniform ray/status closure",
        f"suffix_source_sha256={file_sha256(SUFFIX_PATH)}",
        f"fibre_source_sha256={file_sha256(FIBRE_PATH)}",
        "scope=body (1,4,8,10,12,14), z1=250, three arbitrary later "
        "nonaligned drifts; no finite horizon",
        f"carrier_h={ftext(h)};L={L};cell={cell};"
        f"lower={ftext(lower)};delta250={ftext(first_delta)}",
        "ray_law=(La+b)(c(La+b)-h/7)=b(c(b)-h/7);"
        "b=0 amplitude 0;A(L-b)=-A(b)",
        f"positive_control=b260;amplitude={ftext(ray_amplitudes[260])};"
        f"delta12020={ftext(delta(carrier,h,12020))}",
        f"negative_control=b980;amplitude={ftext(ray_amplitudes[980])};"
        f"delta12740={ftext(delta(carrier,h,12740))}",
        f"zero_control=b5880;amplitude=0;delta17640=0",
        f"ray_amplitude_sha256={ray_digest}",
        f"pair_overlap_exhaustive_controls={pair_rows}",
        "common_table_controls=loads(9,9):FEASIBLE;"
        "loads(14,8):EXACT_FARKAS_INFEASIBLE",
        f"nontrivial_denominators={len(divisors)};"
        f"all_denominator_multisets={len(shapes)}",
        f"ray_excess_kills={len(shapes)-len(excess_survivors)};"
        f"ray_excess_survivors={len(excess_survivors)};"
        f"ray_excess_survivor_sha256={excess_digest.hexdigest()}",
        f"crude_all_q_kills={crude_kills};"
        f"crude_all_q_survivors={len(crude_survivors)}",
        f"common_hunter_status_kills={len(status_kills)};"
        f"common_hunter_status_survivors={len(status_survivors)}",
        f"decisive_inner_M_histogram={dict(sorted(kills_by_M.items()))}",
        f"first_kill=ds{first_shape};D={first_D};q={q};M={M};"
        f"marginals={marginals};capacities={capacities};"
        f"target_histogram={histogram};exact_farkas=VERIFIED",
        (
            "verified_farkas_instance_sha256="
            f"{verified_instance_digest.hexdigest()}"
        ),
        "conclusion=uniform projected z1=250 row is empty",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ]
    output = "\n".join(lines) + "\n"
    OUTPUT_PATH.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
