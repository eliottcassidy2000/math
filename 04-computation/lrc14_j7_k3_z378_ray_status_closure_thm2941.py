#!/usr/bin/env python3
"""All-label closure of the nine projected k=3 scalar rows at z1=378.

The inherited finite-envelope audit leaves nine possible six-body rows at
first drift 378.  For a body ruler L and singleton excess delta, the exact
periodic primitive identity gives

    (z+L) delta(z+L) = z delta(z).

Moreover A(L-b)=-A(b) for A(b)=b delta(b).  Reversal preserves the exact
denominator d=L/gcd(L,b), so every denominator class has a positive ray or
an attained zero ray.  The exact three-label maximum in a fixed denominator
multiset is therefore a finite merge of the first three admissible points
on its nonnegative rays.  A two-state max-plus DP retains the projected
high-label obligation.

The exact ray quotient eliminates four of the nine bodies.  On every
denominator state of the other five, this referee applies the crude all-q
fibre bound and then the common 16-status Hunter transport imported from
the independently controlled z1=250 referee.  Every status rejection has
an exactly verified rational Farkas certificate.  No label horizon or
omitted-tail estimate occurs in this computation.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAY_STATUS_PATH = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k3_uniform_ray_status_closure_thm2941.py"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k3_z378_ray_status_closure_thm2941.out"
)
EXPECTED_RAY_STATUS_SHA256 = (
    "dfa4788297b8c31fc9b5dce1afadf29d20b267cb4159fa95dadb9346b1980b36"
)
EXPECTED_SEMANTIC_SHA256 = (
    "4cc1d3528844818857a3ed4c1c08ca03405cb3857e102c95acae5a95462c4669"
)

FIRST = 378
ORIGINAL_BODIES = (
    (1, 2, 6, 10, 12, 14),
    (1, 6, 8, 10, 12, 14),
    (1, 6, 8, 10, 13, 14),
    (1, 6, 9, 10, 12, 14),
    (1, 6, 10, 12, 13, 14),
    (2, 6, 8, 10, 12, 14),
    (2, 6, 8, 10, 13, 14),
    (2, 6, 10, 12, 13, 14),
    (2, 8, 10, 12, 13, 14),
)
EXPECTED_SURVIVING_BODIES = (
    (1, 2, 6, 10, 12, 14),
    (1, 6, 8, 10, 12, 14),
    (1, 6, 9, 10, 12, 14),
    (2, 6, 8, 10, 12, 14),
    (2, 8, 10, 12, 13, 14),
)
EXPECTED_COUNTS = {
    (1, 2, 6, 10, 12, 14): (11, 3, 8),
    (1, 6, 8, 10, 12, 14): (2, 1, 1),
    (1, 6, 8, 10, 13, 14): (0, 0, 0),
    (1, 6, 9, 10, 12, 14): (15, 10, 5),
    (1, 6, 10, 12, 13, 14): (0, 0, 0),
    (2, 6, 8, 10, 12, 14): (4, 0, 4),
    (2, 6, 8, 10, 13, 14): (0, 0, 0),
    (2, 6, 10, 12, 13, 14): (0, 0, 0),
    (2, 8, 10, 12, 13, 14): (6, 4, 2),
}


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
    file_sha256(RAY_STATUS_PATH) == EXPECTED_RAY_STATUS_SHA256,
    "uniform ray/status dependency changed",
)
local = load_module("k3_uniform_ray_status", RAY_STATUS_PATH)
suffix = local.suffix
fibre = local.fibre
support = local.support


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


class Stream:
    """Exact geometry and inherited fibre screen for one fixed body."""

    def __init__(self, body):
        self.first = FIRST
        self.body = tuple(body)
        self.carrier = suffix.A.carrier_for(self.body)
        self.h = F(
            sum(right - left for left, right in self.carrier),
            suffix.A.RULER,
        )
        self.L = 14 * lcm(*self.body)
        self.lower = self.h * suffix.ETAS[3]
        wall = suffix.PROJECTED_RATIOS[3] * self.L
        self.high_floor = max(
            15, wall.numerator // wall.denominator + 1
        )
        self.first_delta = (
            suffix.A.singleton_coverage(self.carrier, FIRST)
            - self.h / 7
        )
        self.first_d = self.L // gcd(self.L, FIRST)
        actual_L, ranges = support.safe_cell_ranges(self.body)
        require(actual_L == self.L, "body ruler changed")
        self.ranges = tuple(ranges)

    @lru_cache(maxsize=None)
    def target_data(self, D):
        arcs = fibre.projected_support_arcs(D, self.ranges)
        return tuple(
            (
                q,
                max(
                    load
                    for load, multiplicity
                    in fibre.residue_load_histogram(arcs, q)
                    if multiplicity
                ),
            )
            for q in support.divisors(D)
        )

    @lru_cache(maxsize=None)
    def pattern_screen(self, ds):
        """Inherited crude individual all-q fibre-capacity screen."""
        D = lcm(*ds)
        best = None
        for q, target in self.target_data(D):
            capacity = sum(local.fibre_cap(D, d, q) for d in ds)
            candidate = (
                target - capacity,
                q,
                D // q,
                target,
                capacity,
            )
            if best is None or candidate > best:
                best = candidate
        require(best is not None, "missing pattern witness")
        return best[0] > 0, best


def first_on_ray(residue, modulus, threshold):
    """Least positive integer congruent to residue mod modulus, >= threshold."""
    if residue >= threshold:
        return residue
    return residue + ((threshold - residue + modulus - 1) // modulus) * modulus


def ray_class_tables(stream):
    """Top-three all-n candidates in every denominator/high class.

    The exact grid law

        delta(z+L) = z/(z+L) delta(z)

    says that one nonzero residue ray has the form K/z.  Positive rays
    decrease and zero rays are constant.
    Three entries per ray therefore suffice for a three-slot optimization.
    Moreover the centered primitive is odd, so antipodal residues have
    amplitudes ``A_(L-b)=-A_b`` inside the same denominator class.  Every
    class consequently has a positive ray or an attained zero ray, making
    the top-three envelope an actual maximum rather than a mere supremum.
    """
    arbitrary_lists = {}
    high_lists = {}
    amplitudes = {}
    recurrence_checks = 0
    signs = Counter()
    for residue in range(1, stream.L):
        z = first_on_ray(residue, stream.L, stream.first + 1)
        delta = (
            suffix.A.singleton_coverage(stream.carrier, z)
            - stream.h / 7
        )
        constant = z * delta
        amplitudes[residue] = constant
        d = stream.L // gcd(stream.L, residue)
        signs[(d, (constant > 0) - (constant < 0))] += 1

        # Audit the exact functional equation on every residue ray used.
        next_delta = (
            suffix.A.singleton_coverage(
                stream.carrier, z + stream.L
            )
            - stream.h / 7
        )
        require(
            next_delta == F(z, z + stream.L) * delta,
            ("residue-ray law failed", stream.body, residue, z),
        )
        recurrence_checks += 1

        if constant >= 0:
            entries = arbitrary_lists.setdefault(d, [])
            for offset in range(3):
                label = z + offset * stream.L
                entries.append((constant / label, label))
        high_z = first_on_ray(
            residue,
            stream.L,
            max(stream.first + 1, stream.high_floor),
        )
        if constant >= 0:
            entries = high_lists.setdefault(d, [])
            for offset in range(3):
                label = high_z + offset * stream.L
                entries.append((constant / label, label))

    for residue, amplitude in amplitudes.items():
        require(
            amplitude + amplitudes[stream.L - residue] == 0,
            ("antipodal ray law failed", stream.body, residue),
        )

    tables = {}
    for d in support.divisors(stream.L):
        if d == 1:
            continue
        any_candidates = sorted(
            arbitrary_lists.get(d, ()),
            key=lambda item: (item[0], -item[1]),
            reverse=True,
        )[:3]
        high_candidates = sorted(
            high_lists.get(d, ()),
            key=lambda item: (item[0], -item[1]),
            reverse=True,
        )[:3]
        # Antipodal amplitudes satisfy A_(L-b)=-A_b within the same
        # denominator class.  Hence every class has a positive ray or an
        # attained zero ray (including the self-antipodal d=2 class), so its
        # top three values are genuine labels rather than remote suprema.
        require(
            len(any_candidates) == len(high_candidates) == 3,
            ("nonnegative ray missing", stream.body, d),
        )
        tables[d] = (
            tuple(any_candidates),
            tuple(high_candidates),
        )
    return tables, recurrence_checks, signs


def ray_class_choice(stream, tables, d, multiplicity, require_high):
    """Exact all-n scalar maximum for one denominator class."""
    arbitrary, high = tables[d]
    chosen = arbitrary[:multiplicity]
    total = sum((value for value, _label in chosen), F(0))

    def is_high(label):
        return label >= stream.high_floor

    if not require_high or any(is_high(label) for _value, label in chosen):
        return total, tuple(label for _value, label in chosen)
    high_value, high_label = high[0]
    rest = [
        item for item in arbitrary if item[1] != high_label
    ][: multiplicity - 1]
    require(len(rest) == multiplicity - 1, ("ray class too small", d))
    return (
        high_value + sum((value for value, _label in rest), F(0)),
        (high_label, *(label for _value, label in rest)),
    )


def ray_type_upper(stream, tables, suffix_ds, need_high):
    """Two-state all-n max-plus DP for one denominator multiset."""
    dp = {False: (F(0), ())}

    def is_high(label):
        return label >= stream.high_floor

    for d, multiplicity in sorted(Counter(suffix_ds).items()):
        arbitrary = ray_class_choice(
            stream, tables, d, multiplicity, False
        )
        forced = ray_class_choice(stream, tables, d, multiplicity, True)
        options = (arbitrary,) if forced == arbitrary else (arbitrary, forced)
        next_dp = {}
        for seen_high, (total, labels) in dp.items():
            for value, selected in options:
                next_high = seen_high or any(map(is_high, selected))
                candidate = (total + value, (*labels, *selected))
                incumbent = next_dp.get(next_high)
                if incumbent is None or candidate[0] > incumbent[0]:
                    next_dp[next_high] = candidate
        dp = next_dp
    allowed = [
        value for seen_high, value in dp.items() if seen_high or not need_high
    ]
    return max(allowed, key=lambda item: item[0]) if allowed else None


def ray_quotient_states(stream):
    """Exact attained all-n denominator quotient, with no horizon."""
    tables, checks, signs = ray_class_tables(stream)
    alphabet = tuple(sorted(tables))
    trials = comb(len(alphabet) + 2, 3)
    states = {}
    for suffix_ds in combinations_with_replacement(alphabet, 3):
        choice = ray_type_upper(
            stream,
            tables,
            suffix_ds,
            stream.first < stream.high_floor,
        )
        if choice is None:
            continue
        suffix_upper, labels = choice
        total = stream.first_delta + suffix_upper
        if total < stream.lower:
            continue
        full_ds = tuple(sorted((stream.first_d, *suffix_ds)))
        require(full_ds not in states, ("duplicate all-n state", full_ds))
        states[full_ds] = {
            "labels": (stream.first, *labels),
            "total": total,
            "excess": total - stream.lower,
        }
    return trials, states, checks, signs


def common_status_screen(stream, states):
    """Apply crude all-q, then the rigorous common Hunter-status screen."""
    crude_kills = {}
    status_kills = {}
    survivors = []
    arcs_cache = {}
    for ds in sorted(states):
        killed, witness = stream.pattern_screen(ds)
        if killed:
            crude_kills[ds] = witness
            continue
        D = lcm(*ds)
        arcs = arcs_cache.setdefault(
            D,
            fibre.projected_support_arcs(D, stream.ranges),
        )
        status_witness = None
        for M in support.divisors(D):
            q = D // M
            marginals, capacities = local.hunter_status_data(D, ds, q)
            histogram = fibre.residue_load_histogram(arcs, q)
            feasible, certificate = local.common_status_feasible(
                q,
                marginals,
                capacities,
                histogram,
            )
            if not feasible:
                status_witness = (
                    q,
                    M,
                    marginals,
                    tuple(sorted(set(capacities))),
                    histogram,
                    certificate,
                )
                break
        if status_witness is None:
            survivors.append(ds)
        else:
            status_kills[ds] = status_witness
    return crude_kills, status_kills, tuple(survivors)


def main():
    require(
        suffix.ETAS[3] == F(3, 91)
        and suffix.PROJECTED_RATIOS[3] == F(13, 132),
        "projected k=3 constants changed",
    )
    pair_rows, control_marginals, control_caps, control_certificate = (
        local.controls()
    )
    total_states = total_crude = total_status = total_survivors = 0
    surviving_bodies = []
    records = []
    certificate_digest = sha256()
    for body in ORIGINAL_BODIES:
        stream = Stream(body)
        trials, states, checks, signs = ray_quotient_states(stream)
        crude, status, survivors = common_status_screen(stream, states)
        require(
            (len(states), len(crude), len(status))
            == EXPECTED_COUNTS[body],
            (
                "body counts changed",
                body,
                len(states),
                len(crude),
                len(status),
            ),
        )
        if states:
            surviving_bodies.append(body)
        total_states += len(states)
        total_crude += len(crude)
        total_status += len(status)
        total_survivors += len(survivors)
        maximum = max(
            states.items(),
            key=lambda item: (item[1]["excess"], item[0]),
            default=None,
        )
        sign_totals = {
            sign: sum(
                count
                for (_d, candidate_sign), count in signs.items()
                if candidate_sign == sign
            )
            for sign in (-1, 0, 1)
        }
        require(
            sign_totals[-1] == sign_totals[1],
            ("ray sign imbalance", body, sign_totals),
        )
        status_histogram = Counter(
            witness[1] for witness in status.values()
        )
        for ds, witness in sorted(status.items()):
            certificate_digest.update(
                f"{body}|{ds}|{witness}\n".encode()
            )
        records.append(
            (
                body,
                stream.h,
                len(stream.carrier),
                stream.L,
                stream.high_floor,
                stream.first_d,
                checks,
                tuple(sorted(sign_totals.items())),
                len(support.divisors(stream.L)) - 1,
                trials,
                tuple(sorted(states.items())),
                tuple(sorted(crude.items())),
                tuple(sorted(status.items())),
                tuple(survivors),
                tuple(sorted(status_histogram.items())),
                maximum,
            )
        )
    require(
        tuple(surviving_bodies) == EXPECTED_SURVIVING_BODIES,
        ("nine-to-five scalar reduction changed", surviving_bodies),
    )
    require(
        (total_states, total_crude, total_status, total_survivors)
        == (38, 18, 20, 0),
        (
            "global closure counts changed",
            total_states,
            total_crude,
            total_status,
            total_survivors,
        ),
    )

    semantic_payload = (
        FIRST,
        ORIGINAL_BODIES,
        tuple(surviving_bodies),
        tuple(records),
        pair_rows,
        control_marginals,
        control_caps,
        control_certificate,
        certificate_digest.hexdigest(),
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_hash == EXPECTED_SEMANTIC_SHA256,
            "z1=378 closure semantic digest changed",
        )

    eliminated = tuple(
        body for body in ORIGINAL_BODIES if body not in surviving_bodies
    )
    lines = [
        "LRC14 projected k=3 z1=378 all-label ray/status closure",
        f"ray_status_source_sha256={file_sha256(RAY_STATUS_PATH)}",
        (
            "scope=nine inherited projected scalar body rows at z1=378;"
            "three distinct later nonaligned labels;no finite horizon"
        ),
        (
            "ray_law=(z+L)delta(z+L)=zdelta(z);"
            "A(L-b)=-A(b);all denominator maxima attained"
        ),
        f"pair_overlap_exhaustive_controls={pair_rows}",
        (
            "common_table_controls=loads(9,9):FEASIBLE;"
            "loads(14,8):EXACT_FARKAS_INFEASIBLE"
        ),
        (
            f"candidate_bodies={len(ORIGINAL_BODIES)};"
            f"ray_scalar_surviving_bodies={len(surviving_bodies)};"
            f"eliminated_bodies={eliminated}"
        ),
    ]
    for record in records:
        (
            body,
            h,
            components,
            L,
            high_floor,
            first_d,
            checks,
            sign_totals,
            denominator_classes,
            trials,
            states,
            crude,
            status,
            survivors,
            status_histogram,
            maximum,
        ) = record
        state_ds = tuple(ds for ds, _state in states)
        crude_ds = tuple(ds for ds, _witness in crude)
        status_ds = tuple(ds for ds, _witness in status)
        lines.extend(
            (
                (
                    f"E={body};h={ftext(h)};r={components};L={L};"
                    f"high={high_floor};d1={first_d};"
                    f"ray_checks={checks};ray_signs={dict(sign_totals)};"
                    f"denominator_classes={denominator_classes};"
                    f"denominator_trials={trials};"
                    f"scalar_states={len(states)};crude_kills={len(crude)};"
                    f"status_kills={len(status)};survivors={len(survivors)};"
                    f"status_M_histogram={dict(status_histogram)}"
                ),
                f"  scalar_state_denominators={state_ds}",
                f"  crude_kill_denominators={crude_ds}",
                f"  status_kill_denominators={status_ds}",
                f"  max_state={maximum}",
            )
        )
    lines.extend(
        (
            (
                f"totals=states:{total_states};crude_kills:{total_crude};"
                f"status_kills:{total_status};survivors:{total_survivors}"
            ),
            f"farkas_certificate_sha256={certificate_digest.hexdigest()}",
            "conclusion=all nine z1=378 candidate body rows are empty",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    output = "\n".join(lines) + "\n"
    OUTPUT_PATH.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
