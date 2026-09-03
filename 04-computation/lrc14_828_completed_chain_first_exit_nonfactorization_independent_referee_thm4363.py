#!/usr/bin/env python3
"""Clean-room exact audit for proposed THM-4363.

This verifier intentionally does not import the scout implementation.  It
reconstructs the open-tooth cover directly from Fraction intervals, audits
the inherited finite gates, deletes the twelve u=3/39 collars, and compares
the four declared P rows on the remaining physical components.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt
import sys


sys.stdout.reconfigure(newline="\n")

H = 420
ANCHOR = 2 * H
U = 3
CU = 39
FIXED = (11, 1691, 3371, 5051, 6731, 8411, 10091, 525, 945)
PARAMETERS = (241, 255, 761, 1015)
ROLES = {
    U: "u",
    CU: "13u",
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


@dataclass(frozen=True, order=True)
class Tooth:
    speed: int
    address: int
    left: Fraction
    right: Fraction


@dataclass(frozen=True)
class Step:
    cursor: Fraction
    chosen: Tooth
    right_ties: tuple[Tooth, ...]


@dataclass(frozen=True)
class Trace:
    k: int
    status: str
    chain: tuple[Tooth, ...]
    steps: tuple[Step, ...]
    exit: Fraction | None


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def floor_fraction(x: Fraction) -> int:
    return x.numerator // x.denominator


def ceil_fraction(x: Fraction) -> int:
    return -((-x.numerator) // x.denominator)


def component(k: int) -> tuple[Fraction, Fraction]:
    require(0 <= k < 2 * H, "component address out of range")
    return Fraction(14 * k + 1, 28 * H), Fraction(14 * k + 13, 28 * H)


def tooth(speed: int, address: int) -> Tooth:
    return Tooth(
        speed,
        address,
        Fraction(14 * address - 1, 14 * speed),
        Fraction(14 * address + 1, 14 * speed),
    )


def meeting_teeth(speed: int, left: Fraction, right: Fraction) -> tuple[Tooth, ...]:
    # A deliberately loose exact integer range followed by the literal open
    # interval intersection predicate avoids relying on a derived bound.
    lo = floor_fraction(speed * left) - 2
    hi = ceil_fraction(speed * right) + 2
    return tuple(
        z
        for n in range(lo, hi + 1)
        for z in (tooth(speed, n),)
        if z.left < right and left < z.right
    )


def select(active: tuple[Tooth, ...]) -> tuple[Tooth, tuple[Tooth, ...]]:
    farthest = max(z.right for z in active)
    ties = tuple(sorted((z for z in active if z.right == farthest), key=lambda z: z.speed))
    # Equal right endpoints: wider tooth, equivalently smaller speed, wins.
    chosen = min(ties, key=lambda z: (z.speed, z.address))
    return chosen, ties


def trace_component(speeds: tuple[int, ...], k: int) -> Trace:
    left, right = component(k)
    bank = tuple(z for speed in speeds for z in meeting_teeth(speed, left, right))
    cursor = left
    chain: list[Tooth] = []
    steps: list[Step] = []
    while cursor <= right:
        active = tuple(z for z in bank if z.left < cursor < z.right)
        if not active:
            return Trace(k, "missing", tuple(chain), tuple(steps), cursor)
        chosen, ties = select(active)
        require(not chain or chosen.right > chain[-1].right, "greedy frontier did not advance")
        steps.append(Step(cursor, chosen, ties))
        chain.append(chosen)
        cursor = chosen.right
    return Trace(k, "span" if len(chain) == 1 else "renew", tuple(chain), tuple(steps), None)


def directly_covered(speeds: tuple[int, ...], k: int) -> bool:
    """Literal wall-cell audit, independent of prefix/greedy iteration."""
    left, right = component(k)
    bank = tuple(z for speed in speeds for z in meeting_teeth(speed, left, right))
    walls = {left, right}
    for z in bank:
        if left < z.left < right:
            walls.add(z.left)
        if left < z.right < right:
            walls.add(z.right)
    ordered = tuple(sorted(walls))
    tests = ordered + tuple((a + b) / 2 for a, b in zip(ordered, ordered[1:]))
    return all(any(z.left < x < z.right for z in bank) for x in tests)


def circle_distance(speed: int, t: Fraction) -> Fraction:
    x = speed * t
    residue = x - floor_fraction(x)
    return min(residue, 1 - residue)


def clearance(speeds: tuple[int, ...], t: Fraction) -> tuple[Fraction, tuple[int, ...]]:
    values = tuple((v, circle_distance(v, t)) for v in (ANCHOR,) + speeds)
    m = min(x for _, x in values)
    return m, tuple(v for v, x in values if x == m)


def divisors(n: int) -> set[int]:
    ans = set()
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            ans.add(d)
            ans.add(n // d)
    return ans


def capacity(row: tuple[int, ...], a: int) -> Fraction:
    ans = Fraction(0)
    for v in row:
        if v % a:
            q = a // gcd(a, v)
            ans += Fraction((q + 6) // 7, q)
    return ans


def audit_gates(p: int, speeds: tuple[int, ...]) -> dict[str, object]:
    row = (ANCHOR,) + speeds
    require(len(speeds) == len(set(speeds)) == 12, "tail bank is not twelve distinct speeds")
    require(all(v > 0 and v % 2 for v in speeds), "tail bank is not positive odd")
    require(gcd(*row) == 1, "row is not primitive")
    missing_denominators = tuple(
        m for m in range(2, 15) if not any(v % m == 0 for v in row)
    )
    require(not missing_denominators, "THM-366 denominator gate missing")

    fixed_halfturn = (
        Fraction(1, 2) - Fraction(1, 28 * H),
        Fraction(1, 2) + Fraction(1, 28 * H),
    )
    fixed_halfturn_data = tuple((t,) + clearance(speeds, t) for t in fixed_halfturn)
    require(
        all(item[1] < Fraction(1, 14) for item in fixed_halfturn_data),
        "one of the two fixed half-turn clocks is a witness",
    )

    # Complete integer half-turn grid used by the inherited hostile control.
    modulus = 28 * H
    halfturn_survivors = []
    for j in range(modulus):
        t = Fraction(modulus // 2 + j, modulus)
        if clearance(speeds, t)[0] >= Fraction(1, 14):
            halfturn_survivors.append((j, t))
    # The full translated half-turn grid is a separate hostile-control bank,
    # not one of the analytic gates asserted below.  Retain its result so the
    # theorem's scope can be audited without silently strengthening it.
    expected_grid_counts = {241: 214, 255: 222, 761: 218, 1015: 170}
    require(
        len(halfturn_survivors) == expected_grid_counts[p],
        "complete translated half-turn-grid survivor count changed",
    )

    unit_survivors: dict[int, tuple[int, ...]] = {}
    sign_counts: dict[int, int] = {}
    for q in range(8, 15):
        mod = 2 * q
        units = tuple(a for a in range(1, mod) if gcd(a, mod) == 1)
        survivors = tuple(
            a
            for a in units
            if clearance(speeds, Fraction(a, mod))[0] >= Fraction(1, 14)
        )
        require(not survivors, f"unit bank q={q} contains a witness")
        unit_survivors[q] = survivors
        signs = {
            min(v % mod, mod - v % mod)
            for v in speeds
            if gcd(v, mod) == 1
        }
        sign_counts[q] = len(signs)
        if H % q:
            target = {min(a, mod - a) for a in units}
            require(signs == target, f"unit sign bank q={q} is incomplete")

    represented = sorted({d for v in row for d in divisors(v) if d >= 2})
    cap = {a: capacity(row, a) for a in represented}
    cap_min = min(cap.values())
    cap_argmin = tuple(a for a in represented if cap[a] == cap_min)
    require(cap_min >= 1, "adaptive divisor capacity closes the row")
    return {
        "p": p,
        "missing_denominators": missing_denominators,
        "fixed_halfturn_data": fixed_halfturn_data,
        "halfturn_survivors": tuple(halfturn_survivors),
        "unit_survivors": unit_survivors,
        "unit_sign_counts": sign_counts,
        "capacity_count": len(represented),
        "capacity_min": cap_min,
        "capacity_argmin": cap_argmin,
    }


def component_from_endpoint(x: Fraction) -> tuple[int, str]:
    for k in range(2 * H):
        left, right = component(k)
        if x == left:
            return k, "L"
        if x == right:
            return k, "R"
    raise AssertionError(("not an anchor component endpoint", x))


def containing_tooth(speed: int, k: int) -> Tooth:
    left, right = component(k)
    hits = tuple(
        z
        for z in meeting_teeth(speed, left, right)
        if z.left < left and right < z.right
    )
    require(len(hits) == 1, "collar does not have one strict containing tooth")
    return hits[0]


def audit_collars() -> tuple[tuple[int, ...], tuple[tuple[object, ...], ...]]:
    delta = Fraction(1, 28 * H)
    records = []
    collar_components = set()
    for n in range(U):
        for sigma in (-1, 1):
            raw = Fraction(14 * n + sigma, 14 * U)
            t0 = raw - floor_fraction(raw)
            inward = t0 - sigma * delta
            outward = t0 + sigma * delta
            inward -= floor_fraction(inward)
            outward -= floor_fraction(outward)
            kin, side_in = component_from_endpoint(inward)
            kout, side_out = component_from_endpoint(outward)
            zin = containing_tooth(U, kin)
            zout = containing_tooth(CU, kout)
            lin, rin = component(kin)
            lout, rout = component(kout)
            in_margins = (lin - zin.left, zin.right - rin)
            out_margins = (lout - zout.left, zout.right - rout)
            require(min(in_margins) > 0 and min(out_margins) > 0, "collar containment is not strict")
            # The endpoint adjacent to t0 has the universal delta margin.
            # The other margin is the full-collar slack quoted in the source
            # calculation; reflection reverses left and right.
            require(
                sorted(in_margins) == [Fraction(1, 11760), Fraction(547, 11760)],
                "u=3 collar margins changed",
            )
            require(
                sorted(out_margins) == [Fraction(1, 11760), Fraction(391, 152880)],
                "13u=39 collar margins changed",
            )
            collar_components.update((kin, kout))
            records.append(
                (
                    n,
                    sigma,
                    t0,
                    kin,
                    side_in,
                    zin.address,
                    in_margins,
                    kout,
                    side_out,
                    zout.address,
                    out_margins,
                )
            )
    collars = tuple(sorted(collar_components))
    require(len(collars) == 4 * U == 12, "collar count changed")
    require(
        collars == (19, 20, 259, 260, 299, 300, 539, 540, 579, 580, 819, 820),
        "collar address set changed",
    )
    return collars, tuple(records)


def physical_chain(trace: Trace) -> tuple[tuple[int, int], ...]:
    return tuple((z.speed, z.address) for z in trace.chain)


def role(speed: int, p: int) -> str:
    return "P" if speed == p else ROLES[speed]


def role_chain(trace: Trace, p: int) -> tuple[str, ...]:
    return tuple(role(z.speed, p) for z in trace.chain)


def selection_signature(trace: Trace, p: int) -> tuple[tuple[object, ...], ...]:
    epsilon = int(trace.k >= H)
    ans = []
    for step in trace.steps:
        chosen = step.chosen
        ans.append(
            (
                step.cursor,
                role(chosen.speed, p),
                len(step.right_ties),
                tuple(role(z.speed, p) for z in step.right_ties),
                (2 * chosen.address - epsilon * chosen.speed) % 2,
            )
        )
    return tuple(ans)


def determinant_records(trace: Trace) -> tuple[tuple[int, int], ...]:
    ans = []
    for one, two in zip(trace.chain, trace.chain[1:]):
        det = two.address * one.speed - one.address * two.speed
        q = one.speed + two.speed - 14 * det
        overlap = one.right - two.left
        require(overlap == Fraction(q, 14 * one.speed * two.speed), "overlap identity failed")
        require(det > 0 and q > 0, "transition is not a proper positive crossing")
        ans.append((det, q))
    return tuple(ans)


def stable_completed_bytes(
    residual: tuple[int, ...], traces: dict[int, Trace]
) -> bytes:
    lines = []
    for k in residual:
        tr = traces[k]
        chain = "" if tr.status == "missing" else ",".join(
            f"{z.speed}@{z.address}" for z in tr.chain
        )
        lines.append(f"{k}|{tr.status}|{chain}\n")
    return "".join(lines).encode("ascii")


def stable_status_bytes(residual: tuple[int, ...], traces: dict[int, Trace]) -> bytes:
    code = {"missing": "M", "span": "S", "renew": "R"}
    return "".join(f"{k}|{code[traces[k].status]}\n" for k in residual).encode("ascii")


def audit_sharp_pair(all_traces: dict[int, dict[int, Trace]]) -> dict[str, object]:
    x0 = tooth(945, 26).right
    require(
        x0 == Fraction(73, 2646) == Fraction(1, 36) - Fraction(1, 5292),
        "x0 identity failed",
    )
    rows = {}
    for p, n, r in ((761, 21, 5), (1015, 28, 7)):
        require(p == 36 * n + r and 0 < r < 18 and r % 2 == 1, "P=36n+r data failed")
        active_lhs = abs(73 * p - 2646 * n)
        centered_lhs = abs(p - 147 * r)
        require(active_lhs == abs(18 * n - 73 * r), "Euclidean-remainder activity form failed")
        require(active_lhs < 189, "P tooth is not active at x0")
        require(centered_lhs == 2 * active_lhs < 378, "centered remainder equivalence failed")
        z = tooth(p, n)
        require(z.left < x0 < z.right, "literal open-tooth activity failed")
        wall_formula = Fraction(1, 36) + Fraction(18 - 7 * r, 252 * p)
        require(z.right == wall_formula, "outgoing wall formula failed")
        tr = all_traces[p][23]
        require(tr.chain[-1] == z and tr.exit == z.right, "P tooth is not the terminal selected state")
        rows[p] = (n, r, active_lhs, centered_lhs, wall_formula)
    det = 21 * 1015 - 28 * 761
    require(det == 7, "address determinant changed")
    wall_cross = (14 * 21 + 1) * 1015 - (14 * 28 + 1) * 761
    require(wall_cross == 14 * det + (1015 - 761) == 352, "wall cross numerator changed")
    require(
        tooth(761, 21).right - tooth(1015, 28).right
        == Fraction(wall_cross, 14 * 761 * 1015),
        "wall difference determinant identity failed",
    )
    require(
        role_chain(all_traces[761][23], 761) == role_chain(all_traces[1015][23], 1015),
        "sharp pair role chains differ",
    )
    require(
        selection_signature(all_traces[761][23], 761)
        == selection_signature(all_traces[1015][23], 1015),
        "sharp pair local role/tie/owner-bit signatures differ",
    )
    require(
        all_traces[761][23].exit != all_traces[1015][23].exit,
        "sharp pair metric exits unexpectedly agree",
    )
    return {"x0": x0, "rows": rows, "det": det, "wall_cross": wall_cross}


def main() -> None:
    collars, collar_records = audit_collars()
    residual = tuple(k for k in range(2 * H) if k not in set(collars))
    require(len(residual) == 828, "residual body does not have 828 components")

    all_traces: dict[int, dict[int, Trace]] = {}
    all_gates: dict[int, dict[str, object]] = {}
    for p in PARAMETERS:
        speeds = (U, CU) + FIXED + (p,)
        all_gates[p] = audit_gates(p, speeds)
        traces = {k: trace_component(speeds, k) for k in residual}
        for k, tr in traces.items():
            require(
                (tr.status != "missing") == directly_covered(speeds, k),
                f"greedy/direct cover disagreement P={p}, k={k}",
            )
            if tr.status == "missing":
                require(tr.exit is not None, "missing trace has no exit")
                left, right = component(k)
                bank = tuple(z for speed in speeds for z in meeting_teeth(speed, left, right))
                require(left <= tr.exit <= right, "missing exit lies outside component")
                require(
                    not any(z.left < tr.exit < z.right for z in bank),
                    "reported missing exit is still strictly covered",
                )
        all_traces[p] = traces

    base = all_traces[PARAMETERS[0]]
    base_status = tuple(base[k].status for k in residual)
    for p in PARAMETERS[1:]:
        require(
            tuple(all_traces[p][k].status for k in residual) == base_status,
            f"componentwise status map differs for P={p}",
        )
        for k in residual:
            if base[k].status != "missing":
                require(
                    physical_chain(all_traces[p][k]) == physical_chain(base[k]),
                    f"completed physical chain differs for P={p}, k={k}",
                )
    require(
        all(
            PARAMETERS[0] not in (z.speed for z in base[k].chain)
            for k in residual
            if base[k].status != "missing"
        ),
        "the base varying speed appears in a completed chain",
    )

    census = Counter(base_status)
    require(sum(census.values()) == 828, "residual census total changed")
    completed_serialized = stable_completed_bytes(residual, base)
    completed_digest = sha256(completed_serialized).hexdigest()
    status_serialized = stable_status_bytes(residual, base)
    status_digest = sha256(status_serialized).hexdigest()

    first_missing = {}
    for p in PARAMETERS:
        tr = next(all_traces[p][k] for k in residual if all_traces[p][k].status == "missing")
        require(tr.k == 23 and tr.exit is not None, f"first missing component changed for P={p}")
        speeds = (U, CU) + FIXED + (p,)
        clear, binding = clearance(speeds, tr.exit)
        require(clear == Fraction(1, 14), f"first-exit clearance changed for P={p}")
        require(binding == (tr.chain[-1].speed,), f"first-exit binding owner changed for P={p}")
        first_missing[p] = (tr, clear, binding)
    require(len({tr.exit for tr, _, _ in first_missing.values()}) == 4, "first exits are not pairwise distinct")
    sharp = audit_sharp_pair(all_traces)

    print("THM4363 CLEANROOM EXACT REFEREE")
    print("fixed", FIXED)
    print("parameters", PARAMETERS)
    print("collars", collars)
    print("residual_components", len(residual))
    for rec in collar_records:
        print("collar", rec)
    for p in PARAMETERS:
        g = all_gates[p]
        print(
            "gates",
            p,
            "missing_denominators",
            g["missing_denominators"],
            "fixed_halfturn_data",
            g["fixed_halfturn_data"],
            "complete_grid_safe_count",
            len(g["halfturn_survivors"]),
            "complete_grid_first",
            g["halfturn_survivors"][0],
            "unit_sign_counts",
            g["unit_sign_counts"],
            "capacity_count",
            g["capacity_count"],
            "capacity_min",
            g["capacity_min"],
            "capacity_argmin",
            g["capacity_argmin"],
        )
    print("residual_census", dict(sorted(census.items())))
    print("completed_count", census["span"] + census["renew"])
    print("status_ascii_sha256", status_digest)
    print("status_ascii_bytes", len(status_serialized))
    print("completed_ascii_sha256", completed_digest)
    print("completed_ascii_bytes", len(completed_serialized))
    for p in PARAMETERS:
        tr, clear, binding = first_missing[p]
        print(
            "first_missing",
            p,
            "k",
            tr.k,
            "component",
            component(tr.k),
            "physical_chain",
            physical_chain(tr),
            "role_chain",
            role_chain(tr, p),
            "det_q",
            determinant_records(tr),
            "exit",
            tr.exit,
            "clearance",
            clear,
            "binding",
            binding,
        )
        print("selection_signature", p, selection_signature(tr, p))
    a = all_traces[761][23]
    b = all_traces[1015][23]
    print("sharp_role_equal_761_1015", role_chain(a, 761) == role_chain(b, 1015))
    print("sharp_selection_signature_761", selection_signature(a, 761))
    print("sharp_selection_signature_1015", selection_signature(b, 1015))
    print("sharp_x0", sharp["x0"])
    print("sharp_rows", sharp["rows"])
    print("sharp_determinant", sharp["det"], "wall_cross", sharp["wall_cross"])
    print("PASS")


if __name__ == "__main__":
    main()
