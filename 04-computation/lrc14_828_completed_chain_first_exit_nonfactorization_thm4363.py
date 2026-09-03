#!/usr/bin/env python3
"""Finite-exact verifier for the THM-4363 quotient nonfactorization row.

The script instantiates the h=420, u=3, 13u=39 construction with four actual
odd-tail controls.  It checks the inherited necessary gates, deletes exactly
the twelve forced collars, reconstructs every remaining component by two
independent enumerators, and tests whether the common completed-chain quotient
determines the first residual-exit consumer.  This is not an LRC(14) decrement.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import reduce
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import chain
from math import ceil, floor, gcd, isqrt
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]


def require(condition: bool, message: str) -> None:
    """Keep every finite-exact check live under both normal Python and -O."""
    if not condition:
        raise AssertionError(message)


def load(name: str, rel: str):
    spec = spec_from_file_location(name, ROOT / rel)
    require(spec is not None and spec.loader is not None, f"cannot import {rel}")
    mod = module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


renewal = load(
    "renewal4348_lrc828",
    "04-computation/lrc14_third_tooth_competition_probe_third_tooth_20260902.py",
)
current = load(
    "current4345_lrc828",
    "04-computation/lrc14_anchor_current_strip_probe_anchor_current_20260902.py",
)


H = 420
ANCHOR = 2 * H
U = 3
MULTIPLE = 13 * U
FIXED_NINE = (11, 1691, 3371, 5051, 6731, 8411, 10091, 525, 945)
PARAMETERS = (241, 255, 761, 1015)
GLOBAL_PAIR = (241, 255)
ADDRESS_PAIR = (761, 1015)
ROLE_BY_FIXED_SPEED = {
    U: "U",
    MULTIPLE: "M",
    11: "D0",
    1691: "D1",
    3371: "D2",
    5051: "D3",
    6731: "D4",
    8411: "D5",
    10091: "D6",
    525: "C0",
    945: "C1",
}


def tails(parameter: int) -> tuple[int, ...]:
    return (U, MULTIPLE) + FIXED_NINE + (parameter,)


def role(parameter: int, speed: int) -> str:
    return "P" if speed == parameter else ROLE_BY_FIXED_SPEED[speed]


def distance(x: F) -> F:
    r = x % 1
    return min(r, 1 - r)


def divisors(n: int) -> set[int]:
    answer = {n}
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            answer.add(d)
            answer.add(n // d)
    return answer


def sign_class(a: int, modulus: int) -> int:
    a %= modulus
    return min(a, (-a) % modulus)


def inherited_gate_record(parameter: int) -> dict[str, object]:
    odd = tails(parameter)
    speeds = (ANCHOR,) + odd
    require(len(odd) == len(set(odd)) == 12, "tail bank is not twelve distinct speeds")
    require(all(v > 0 and v % 2 for v in odd), "tail bank is not positive odd")
    require(reduce(gcd, speeds) == 1, "row is not primitive")

    denominator_missing = tuple(
        d for d in range(2, 15) if not any(v % d == 0 for v in speeds)
    )
    require(not denominator_missing, "inherited denominator gate failed")

    halfturn = tuple(F(1, 2) + sign * F(1, 28 * H) for sign in (-1, 1))
    halfturn_records = []
    for t in halfturn:
        values = tuple((v, distance(v * t)) for v in speeds)
        minimum = min(x[1] for x in values)
        blockers = tuple(v for v, d in values if d == minimum)
        require(minimum < F(1, 14), "one of the two fixed half-turn clocks is safe")
        halfturn_records.append((t, minimum, blockers))
    require(
        all(item[1:] == (F(829, 11760), (5051,)) for item in halfturn_records),
        "fixed half-turn blocker record changed",
    )

    # This complete translated grid is a hostile scope control, not a gate:
    # unlike the two named clocks above, it contains many safe points.
    modulus = 28 * H
    complete_grid_safe = []
    for j in range(modulus):
        t = F(modulus // 2 + j, modulus)
        values = tuple((v, distance(v * t)) for v in speeds)
        if min(d for _, d in values) >= F(1, 14):
            complete_grid_safe.append((j, t))
    expected_grid_counts = {241: 214, 255: 222, 761: 218, 1015: 170}
    require(
        len(complete_grid_safe) == expected_grid_counts[parameter],
        "complete translated half-turn-grid safe count changed",
    )

    unit_banks = []
    unit_sign_counts = []
    for q in range(8, 15):
        modulus = 2 * q
        units = tuple(a for a in range(1, modulus) if gcd(a, modulus) == 1)
        survivors = tuple(
            a
            for a in units
            if min(distance(v * F(a, modulus)) for v in speeds) >= F(1, 14)
        )
        represented = tuple(
            sorted(
                {
                    sign_class(v, modulus)
                    for v in odd
                    if gcd(v, modulus) == 1
                }
            )
        )
        if H % q:
            universe = tuple(sorted({sign_class(a, modulus) for a in units}))
            require(represented == universe, f"unit sign bank q={q} is incomplete")
        require(not survivors, f"unit bank q={q} contains a safe point")
        unit_banks.append((q, survivors))
        unit_sign_counts.append((q, len(represented)))

    represented_scales = set(chain.from_iterable(divisors(v) for v in speeds))
    capacities = []
    for a in sorted(represented_scales):
        value = sum(
            (
                F(
                    ((a // gcd(a, v)) + 6) // 7,
                    a // gcd(a, v),
                )
                for v in speeds
                if v % a != 0
            ),
            F(),
        )
        require(value >= 1, f"adaptive divisor capacity closes at a={a}")
        capacities.append((value, a))

    r_counts = tuple((p, sum(v % p != 0 for v in odd)) for p in (3, 5, 7))
    require(dict(r_counts)[3] >= 3, "R_3 count gate failed")
    require(dict(r_counts)[5] >= 5, "R_5 count gate failed")
    require(dict(r_counts)[7] >= 7, "R_7 count gate failed")

    A = sum(v % 3 != 0 for v in odd)
    B = sum(v % 3 == 0 and v % 9 != 0 for v in odd)
    d9 = sum(v % 9 == 0 for v in odd)
    d11 = sum(v % 11 == 0 for v in odd)
    d13 = sum(v % 13 == 0 for v in odd)
    require(d9 >= 1 and 2 * A + 3 * B >= 6, "mod-9 mixed-count gate failed")
    require(1 <= d11 <= 7 and 1 <= d13 <= 6, "mod-11/13 divisor gate failed")

    return {
        "gcd": reduce(gcd, speeds),
        "denominator_missing": denominator_missing,
        "halfturn": tuple(halfturn_records),
        "complete_grid_safe_count": len(complete_grid_safe),
        "complete_grid_first": complete_grid_safe[0],
        "unit_banks": tuple(unit_banks),
        "unit_sign_counts": tuple(unit_sign_counts),
        "capacity_count": len(capacities),
        "capacity_min": min(capacities),
        "capacity_minimizers": tuple(a for value, a in capacities if value == min(capacities)[0]),
        "r_counts": r_counts,
        "mc9": (A, B, d9, d11, d13),
    }


@dataclass(frozen=True)
class Tooth:
    speed: int
    n: int
    left: F
    right: F


def fresh_anchor(k: int) -> tuple[F, F]:
    return F(14 * k + 1, 28 * H), F(14 * k + 13, 28 * H)


def fresh_teeth(speed: int, k: int) -> tuple[Tooth, ...]:
    left, right = fresh_anchor(k)
    lo = floor(speed * left - F(1, 14)) - 1
    hi = ceil(speed * right + F(1, 14)) + 1
    answer = []
    for n in range(lo, hi + 1):
        a = F(14 * n - 1, 14 * speed)
        b = F(14 * n + 1, 14 * speed)
        if a < right and left < b:
            answer.append(Tooth(speed, n, a, b))
    return tuple(answer)


def tooth_key(t: Tooth) -> tuple[F, F, int, int]:
    return t.right, -t.left, -t.speed, -t.n


@dataclass(frozen=True)
class Step:
    frontier: F
    winner: Tooth
    equal_right_competitors: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class Trace:
    status: str
    chain: tuple[Tooth, ...]
    exit_or_end: F
    steps: tuple[Step, ...]


def fresh_trace(parameter: int, k: int) -> Trace:
    left, right = fresh_anchor(k)
    bank = tuple(chain.from_iterable(fresh_teeth(w, k) for w in tails(parameter)))
    frontier = left
    chosen: list[Tooth] = []
    steps: list[Step] = []
    while True:
        active = tuple(t for t in bank if t.left < frontier < t.right)
        if not active:
            return Trace("missing", tuple(chosen), frontier, tuple(steps))
        winner = max(active, key=tooth_key)
        ties = tuple(sorted((t.speed, t.n) for t in active if t.right == winner.right))
        steps.append(Step(frontier, winner, ties))
        chosen.append(winner)
        frontier = winner.right
        if frontier > right:
            return Trace("span" if len(chosen) == 1 else "renew", tuple(chosen), frontier, tuple(steps))


def canonical_trace(parameter: int, k: int) -> Trace:
    left, right = renewal.anchor_component(H, k)
    bank = tuple(
        chain.from_iterable(renewal.meeting_teeth(w, left, right) for w in tails(parameter))
    )
    frontier = left
    chosen: list[Tooth] = []
    steps: list[Step] = []
    while True:
        active = tuple(t for t in bank if t.left < frontier < t.right)
        if not active:
            return Trace("missing", tuple(chosen), frontier, tuple(steps))
        raw = max(active, key=renewal.selection_key)
        winner = Tooth(raw.w, raw.n, raw.left, raw.right)
        ties = tuple(sorted((t.w, t.n) for t in active if t.right == raw.right))
        steps.append(Step(frontier, winner, ties))
        chosen.append(winner)
        frontier = winner.right
        if frontier > right:
            return Trace("span" if len(chosen) == 1 else "renew", tuple(chosen), frontier, tuple(steps))


WALL_RESIDUES = (1, 13, 15, 27, 29, 41)
COLLARS = tuple(
    sorted(
        {
            k
            for r in WALL_RESIDUES
            for k in ((20 * r - 1) % (2 * H), (20 * r) % (2 * H))
        }
    )
)
RESIDUAL = tuple(k for k in range(2 * H) if k not in set(COLLARS))


def audit_collars() -> tuple[tuple[int, tuple[tuple[int, int], ...]], ...]:
    require(
        COLLARS == (19, 20, 259, 260, 299, 300, 539, 540, 579, 580, 819, 820),
        "collar address set changed",
    )
    require(len(RESIDUAL) == 828, "residual body does not have 828 components")
    records = []
    for k in COLLARS:
        left, right = fresh_anchor(k)
        covers = tuple(
            (w, t.n)
            for w in (U, MULTIPLE)
            for t in fresh_teeth(w, k)
            if t.left < left and right < t.right
        )
        require(bool(covers), f"collar k={k} has no declared strict cover")
        records.append((k, covers))
    return tuple(records)


@dataclass(frozen=True)
class Row:
    parameter: int
    traces: tuple[tuple[int, Trace], ...]
    status_map: tuple[tuple[int, str], ...]
    completed: tuple[tuple[int, str, tuple[tuple[int, int], ...]], ...]
    coarse: tuple[tuple[tuple[int, str, tuple[str, ...]], int], ...]
    first_exit: tuple[object, ...]
    attachments: tuple[tuple[object, ...], ...]


def safety_record(parameter: int, k: int, trace: Trace) -> tuple[object, ...]:
    require(trace.status == "missing", "safety record requested for a completed trace")
    x = trace.exit_or_end
    values = tuple((v, distance(v * x)) for v in (ANCHOR,) + tails(parameter))
    clearance = min(d for _, d in values)
    binding = tuple(v for v, d in values if d == clearance)
    require(clearance >= F(1, 14), f"reported exit is not safe for P={parameter}, k={k}")
    return k, k % H, int(k >= H), x, clearance, binding


def row(parameter: int) -> Row:
    trace_rows = []
    statuses = []
    completed = []
    coarse_counter: Counter[tuple[int, str, tuple[str, ...]]] = Counter()
    attachments = []
    occurrence: Counter[tuple[int, str, str]] = Counter()
    first_exit = None

    for k in RESIDUAL:
        primary = canonical_trace(parameter, k)
        referee = fresh_trace(parameter, k)
        require(primary == referee, f"two-enumerator disagreement P={parameter}, k={k}")
        trace_rows.append((k, primary))
        statuses.append((k, primary.status))
        addressed = tuple((t.speed, t.n) for t in primary.chain)
        roles = tuple(role(parameter, t.speed) for t in primary.chain)
        if primary.status != "missing":
            completed.append((k, primary.status, addressed))
        elif first_exit is None:
            first_exit = safety_record(parameter, k, primary)
        coarse_counter[(int(k >= H), primary.status, roles if primary.status != "missing" else ())] += 1

        for index, step in enumerate(primary.steps):
            source = step.winner
            source_role = role(parameter, source.speed)
            key = (int(k >= H), source_role, "R")
            ordinal = occurrence[key]
            occurrence[key] += 1
            nxt = primary.chain[index + 1] if index + 1 < len(primary.chain) else None
            owner_nearest_integer = 2 * source.n - int(k >= H) * source.speed
            owner_bit = owner_nearest_integer % 2
            require(owner_bit == int(k >= H), "owner bit/reflection label mismatch")
            attachments.append(
                (
                    k,
                    k % H,
                    int(k >= H),
                    source_role,
                    source.speed,
                    source.n,
                    "R",
                    ordinal,
                    owner_bit,
                    step.frontier,
                    source.right,
                    step.equal_right_competitors,
                    None if nxt is None else role(parameter, nxt.speed),
                    None if nxt is None else nxt.speed,
                    None if nxt is None else nxt.n,
                    primary.status,
                )
            )
    require(first_exit is not None, f"row P={parameter} has no residual exit")
    return Row(
        parameter,
        tuple(trace_rows),
        tuple(statuses),
        tuple(completed),
        tuple(sorted(coarse_counter.items(), key=repr)),
        first_exit,
        tuple(attachments),
    )


def local_current(t: F, odd_speeds: tuple[int, ...]) -> tuple[int, int, int]:
    lower = sum(distance(w * t) < F(1, 14) for w in odd_speeds)
    upper = sum(distance(w * (t + F(1, 2))) < F(1, 14) for w in odd_speeds)
    return lower, upper, lower - upper


def one_sided_current(x: F, banks: tuple[tuple[int, ...], ...]) -> tuple[object, ...]:
    walls = set()
    for bank in banks:
        walls.update(current.tail_events(bank))
    previous = max(w for w in walls if w < x)
    following = min(w for w in walls if w > x)
    left_test = (previous + x) / 2
    right_test = (x + following) / 2
    return (
        previous,
        following,
        tuple(local_current(left_test, bank) for bank in banks),
        tuple(local_current(x, bank) for bank in banks),
        tuple(local_current(right_test, bank) for bank in banks),
    )


def trace_repr(parameter: int, k: int, trace: Trace) -> tuple[tuple[object, ...], ...]:
    return tuple(
        (
            (k, k % H, int(k >= H), index),
            (2 * step.winner.n - int(k >= H) * step.winner.speed) % 2,
            (role(parameter, step.winner.speed), step.winner.speed, step.winner.n, "R"),
            step.frontier,
            step.winner.right,
            step.equal_right_competitors,
            (
                "EXIT"
                if index + 1 == len(trace.chain)
                else (
                    role(parameter, trace.chain[index + 1].speed),
                    trace.chain[index + 1].speed,
                    trace.chain[index + 1].n,
                )
            ),
            trace.status,
        )
        for index, step in enumerate(trace.steps)
    )


def label_signature(parameter: int, k: int, trace: Trace) -> tuple[tuple[object, ...], ...]:
    """Role-level origin/owner/tie/frontier/continuation data, without addresses."""
    epsilon = int(k >= H)
    return tuple(
        (
            (k, k % H, epsilon, index),
            (2 * step.winner.n - epsilon * step.winner.speed) % 2,
            role(parameter, step.winner.speed),
            "R",
            step.frontier,
            tuple(role(parameter, speed) for speed, _ in step.equal_right_competitors),
            (
                "EXIT"
                if index + 1 == len(trace.chain)
                else role(parameter, trace.chain[index + 1].speed)
            ),
            trace.status,
        )
        for index, step in enumerate(trace.steps)
    )


def determinant_records(trace: Trace) -> tuple[tuple[int, int], ...]:
    """Return exact (determinant, overlap numerator) for consecutive teeth."""
    answer = []
    for one, two in zip(trace.chain, trace.chain[1:]):
        determinant = two.n * one.speed - one.n * two.speed
        overlap_numerator = one.speed + two.speed - 14 * determinant
        require(determinant > 0, "transition determinant is not positive")
        require(overlap_numerator > 0, "successive open teeth do not overlap")
        require(
            one.right - two.left
            == F(overlap_numerator, 14 * one.speed * two.speed),
            "transition overlap identity failed",
        )
        answer.append((determinant, overlap_numerator))
    return tuple(answer)


def stable_status_bytes(record: Row) -> bytes:
    """Declared ASCII representation of every labelled residual status."""
    code = {"missing": "M", "span": "S", "renew": "R"}
    traces = dict(record.traces)
    return "".join(f"{k}|{code[traces[k].status]}\n" for k in RESIDUAL).encode("ascii")


def stable_completed_bytes(record: Row) -> bytes:
    """Declared ASCII representation of all statuses and completed chains."""
    traces = dict(record.traces)
    lines = []
    for k in RESIDUAL:
        trace = traces[k]
        completed_chain = "" if trace.status == "missing" else ",".join(
            f"{tooth.speed}@{tooth.n}" for tooth in trace.chain
        )
        lines.append(f"{k}|{trace.status}|{completed_chain}\n")
    return "".join(lines).encode("ascii")


def audit_address_formulas(rows: dict[int, Row]) -> dict[str, object]:
    """Explain the sharp P=761/1015 address change by centered residues."""
    traces = {parameter: dict(record.traces) for parameter, record in rows.items()}
    x0 = F(14 * 26 + 1, 14 * 945)
    require(
        x0 == F(73, 2646) == F(1, 36) - F(1, 5292),
        "common-frontier identity failed",
    )
    formula_rows = {}
    for parameter, address, remainder in ((761, 21, 5), (1015, 28, 7)):
        require(
            parameter == 36 * address + remainder and 0 < remainder < 18,
            "centered quotient/remainder data failed",
        )
        activity = abs(73 * parameter - 2646 * address)
        require(
            activity == abs(18 * address - 73 * remainder),
            "Euclidean-remainder activity identity failed",
        )
        require(activity < 189, "declared parameter tooth is not active at x0")
        require(
            abs(parameter - 147 * remainder) == 2 * activity < 378,
            "centered activity equivalence failed",
        )
        outgoing = F(14 * address + 1, 14 * parameter)
        wall_formula = F(1, 36) + F(18 - 7 * remainder, 252 * parameter)
        require(outgoing == wall_formula, "outgoing-wall formula failed")
        trace = traces[parameter][23]
        require(
            trace.chain[-1] == Tooth(parameter, address, F(14 * address - 1, 14 * parameter), outgoing),
            "declared parameter tooth is not the terminal selected state",
        )
        require(trace.exit_or_end == outgoing, "terminal wall is not the first exit")
        formula_rows[parameter] = (
            address,
            remainder,
            activity,
            abs(parameter - 147 * remainder),
            outgoing,
        )
    determinant = 21 * 1015 - 28 * 761
    wall_cross = (14 * 21 + 1) * 1015 - (14 * 28 + 1) * 761
    require(determinant == 7, "sharp address determinant changed")
    require(
        wall_cross == 14 * determinant + (1015 - 761) == 352,
        "sharp wall-cross numerator changed",
    )
    require(
        formula_rows[761][-1] - formula_rows[1015][-1]
        == F(wall_cross, 14 * 761 * 1015),
        "sharp wall-difference identity failed",
    )
    return {
        "x0": x0,
        "rows": formula_rows,
        "determinant": determinant,
        "wall_cross": wall_cross,
    }


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")
    print("THM4363_LRC14_H420_U3_COMPLETED_CHAIN_FIRST_EXIT_NONFACTORIZATION=FINITE_EXACT")
    print("CLAIM=DECLARED_FINITE_QUOTIENT_NONFACTORIZATION;NOT_AN_LRC14_DECREMENT")
    print(f"anchor={ANCHOR};u={U};13u={MULTIPLE};fixed_nine={FIXED_NINE};P={PARAMETERS}")
    collar_records = audit_collars()
    print(f"collars={COLLARS};count={len(COLLARS)};residual_count={len(RESIDUAL)}")
    print(f"collar_cover_records={collar_records}")

    gate_records = tuple((p, inherited_gate_record(p)) for p in PARAMETERS)
    for p, record in gate_records:
        print(f"gate_P{p}={record}")
        print(
            f"halfturn_scope_P{p}=fixed_clocks_blocked:2;"
            f"complete_translated_grid_safe_count:{record['complete_grid_safe_count']};"
            f"complete_translated_grid_first:{record['complete_grid_first']}"
        )

    rows = tuple(row(p) for p in PARAMETERS)
    row_by_parameter = {record.parameter: record for record in rows}
    prototype = rows[0]
    require(
        all(record.status_map == prototype.status_map for record in rows),
        "componentwise residual status maps differ",
    )
    require(
        all(record.completed == prototype.completed for record in rows),
        "completed physical chains differ",
    )
    require(
        all(record.coarse == prototype.coarse for record in rows),
        "phase/status/completed-role multisets differ",
    )
    require(
        len({record.first_exit for record in rows}) == len(rows),
        "first residual-exit records are not pairwise distinct",
    )

    status_counts = Counter(status for _, status in prototype.status_map)
    require(
        status_counts == Counter({"missing": 546, "span": 276, "renew": 6}),
        "residual status census changed",
    )
    require(len(prototype.completed) == 282, "completed-chain count changed")
    status_serialized = stable_status_bytes(prototype)
    completed_serialized = stable_completed_bytes(prototype)
    require(len(status_serialized) == 4860, "status ASCII byte length changed")
    require(len(completed_serialized) == 10940, "completed ASCII byte length changed")
    status_digest = sha256(status_serialized).hexdigest()
    completed_digest = sha256(completed_serialized).hexdigest()
    require(
        status_digest
        == "4d34447b9eca8c8a9302a0f799a56300b4b96135cd0c3a245a97960975f9a347",
        "status ASCII digest changed",
    )
    require(
        completed_digest
        == "64704c712a0e3e70ca0ecb3264834b3f61c4d473ddc20ea8c44ccc5c2c616d11",
        "completed ASCII digest changed",
    )
    require(
        all(stable_status_bytes(record) == status_serialized for record in rows),
        "declared status representation differs across rows",
    )
    require(
        all(stable_completed_bytes(record) == completed_serialized for record in rows),
        "declared completed representation differs across rows",
    )

    print("declared_common_quotient")
    print(f"status_counts={tuple(sorted(status_counts.items()))}")
    print(
        f"completed_count={len(prototype.completed)};"
        f"completed_exact_equal_all={all(record.completed == prototype.completed for record in rows)}"
    )
    completed_contains_parameter = tuple(
        any(row.parameter in tuple(speed for speed, _ in record[2]) for record in row.completed)
        for row in rows
    )
    require(
        completed_contains_parameter == tuple(False for _ in rows),
        "a varying parameter appears in a completed chain",
    )
    print(f"completed_contains_parameter={completed_contains_parameter}")
    print(f"coarse={prototype.coarse}")
    print(
        f"status_map_equal_all={all(record.status_map == prototype.status_map for record in rows)};"
        f"coarse_equal_all={all(record.coarse == prototype.coarse for record in rows)}"
    )
    print(f"status_ascii_bytes={len(status_serialized)};status_ascii_sha256={status_digest}")
    print(
        f"completed_ascii_bytes={len(completed_serialized)};"
        f"completed_ascii_sha256={completed_digest}"
    )

    expected_pair_statistics = {
        "global": (134, 132, 0, 1661),
        "address": (160, 154, 0, 1529),
    }
    for pair_name, pair in (("global", GLOBAL_PAIR), ("address", ADDRESS_PAIR)):
        a = row_by_parameter[pair[0]]
        b = row_by_parameter[pair[1]]
        trace_a = dict(a.traces)
        trace_b = dict(b.traces)
        changed = tuple(k for k in RESIDUAL if trace_a[k] != trace_b[k])
        changed_nonmissing = tuple(
            k
            for k in changed
            if trace_a[k].status != "missing" or trace_b[k].status != "missing"
        )
        require(not changed_nonmissing, f"{pair_name} pair changes a completed trace")
        role_changed = tuple(
            k
            for k in changed
            if tuple(role(a.parameter, t.speed) for t in trace_a[k].chain)
            != tuple(role(b.parameter, t.speed) for t in trace_b[k].chain)
        )
        attachment_difference = len(set(a.attachments) ^ set(b.attachments))
        require(
            (len(changed), len(role_changed), len(changed_nonmissing), attachment_difference)
            == expected_pair_statistics[pair_name],
            f"{pair_name} pair-difference census changed",
        )
        print(
            f"{pair_name}_pair={pair};changed_partial_traces={len(changed)};"
            f"role_changed={len(role_changed)};changed_nonmissing={len(changed_nonmissing)};"
            f"attachment_symmetric_difference={attachment_difference};"
            f"first_changed_components={changed[:20]}"
        )

    interval = fresh_anchor(23)
    print("first_residual_exit_consumer")
    print(f"component=23;I23={interval}")
    expected_first = {
        241: (((945, 26), (3371, 93)), ((239, 970),), F(1303, 47194), (3371,)),
        255: (((255, 7), (5051, 140)), ((343, 504),), F(1961, 70714), (5051,)),
        761: (((945, 26), (761, 21)), ((59, 880),), F(295, 10654), (761,)),
        1015: (((945, 26), (1015, 28)), ((70, 980),), F(393, 14210), (1015,)),
    }
    for record in rows:
        k, j, epsilon, exit_point, clearance, binding = record.first_exit
        trace = dict(record.traces)[k]
        physical_chain = tuple((tooth.speed, tooth.n) for tooth in trace.chain)
        expected_chain, expected_det_q, expected_exit, expected_binding = expected_first[
            record.parameter
        ]
        require((k, j, epsilon) == (23, 23, 0), "first residual component changed")
        require(physical_chain == expected_chain, "first residual physical chain changed")
        require(determinant_records(trace) == expected_det_q, "first residual det/q changed")
        require(exit_point == expected_exit, "first residual exit point changed")
        require(clearance == F(1, 14), "first residual exit clearance changed")
        require(binding == expected_binding, "first residual exit binding set changed")

    for pair_name, pair in (("global", GLOBAL_PAIR), ("address", ADDRESS_PAIR)):
        a = row_by_parameter[pair[0]]
        b = row_by_parameter[pair[1]]
        trace_a = dict(a.traces)
        trace_b = dict(b.traces)
        k_a, _, _, x_a, _, _ = a.first_exit
        k_b, _, _, x_b, _, _ = b.first_exit
        require(k_a == k_b == 23 and x_a != x_b, "paired first-exit consumer did not separate")
        roles_a = tuple(role(a.parameter, t.speed) for t in trace_a[23].chain)
        roles_b = tuple(role(b.parameter, t.speed) for t in trace_b[23].chain)
        if pair_name == "address":
            require(roles_a == roles_b == ("C1", "P"), "sharp role chains differ")
            require(
                label_signature(a.parameter, 23, trace_a[23])
                == label_signature(b.parameter, 23, trace_b[23]),
                "sharp origin/owner/tie/frontier/continuation labels differ",
            )
        print(f"{pair_name}_pair={pair};local_roles=({roles_a},{roles_b})")
        print(
            f"P{a.parameter}_det_q={determinant_records(trace_a[23])};"
            f"P{b.parameter}_det_q={determinant_records(trace_b[23])}"
        )
        print(f"P{a.parameter}_trace={trace_repr(a.parameter, 23, trace_a[23])}")
        print(f"P{b.parameter}_trace={trace_repr(b.parameter, 23, trace_b[23])}")
        print(f"P{a.parameter}_attachments_k23={tuple(x for x in a.attachments if x[0] == 23)}")
        print(f"P{b.parameter}_attachments_k23={tuple(x for x in b.attachments if x[0] == 23)}")
        print(f"P{a.parameter}_exit={a.first_exit}")
        print(f"P{b.parameter}_exit={b.first_exit}")
        banks = (tails(a.parameter), tails(b.parameter))
        print(f"P{a.parameter}_exit_current={one_sided_current(x_a, banks)}")
        print(f"P{b.parameter}_exit_current={one_sided_current(x_b, banks)}")

    sharp = audit_address_formulas(row_by_parameter)
    print("sharp_address_formulas")
    print(f"x0={sharp['x0']};P_equals_36n_plus_r_rows={sharp['rows']}")
    print(
        "activity_iff=abs(73P-2646n)<189 iff abs(18n-73r)<189 "
        "iff abs(P-147r)<378"
    )
    print("wall_formula=b(P,n)=1/36+(18-7r)/(252P)")
    print(
        f"address_determinant={sharp['determinant']};"
        f"wall_cross_numerator={sharp['wall_cross']}"
    )

    banks = tuple(tails(p) for p in PARAMETERS)
    profiles = tuple(current.exact_profile(H, bank) for bank in banks)
    print("global_current_controls")
    for p, profile in zip(PARAMETERS, profiles):
        print(
            f"P{p}=core_energy:{profile['core_energy']};"
            f"core_q_cubic:{profile['core_q_cubic']};"
            f"max_abs_current:{profile['core_max_abs_current']}"
        )
    print(f"all_core_energies_distinct={len({profile['core_energy'] for profile in profiles}) == len(profiles)}")
    print(f"all_core_q_cubics_distinct={len({profile['core_q_cubic'] for profile in profiles}) == len(profiles)}")
    print(
        "VERDICT=PASS_DECLARED_COMPLETED_CHAIN_QUOTIENT_DOES_NOT_FACTOR_FIRST_RESIDUAL_EXIT;"
        "LOCAL_ORIGIN_OWNER_TIE_FRONTIER_CONTINUATION_LABELS_DO_NOT_FACTOR_EXIT_ADDRESS"
    )
    print(
        "SCOPE=FINITE_FOUR_SAFE_CONTROLS;ONLY_TWO_FIXED_HALF_TURN_CLOCKS_BLOCKED;"
        "COMPLETE_TRANSLATED_GRID_HAS_SAFE_POINTS_214_222_218_170;LRC14_OPEN;NO_DECREMENT"
    )


if __name__ == "__main__":
    main()
