#!/usr/bin/env python3
"""HYP-2679 / T918: true-wide two-far boundary-curvature scout.

After THM-547, the LRC(14) sector crux is the true-wide branch
`second-largest(E)>14`.  This script tests the boundary-function reframe:
when two far speeds are added to a bounded or dilated core, the missing datum is
not only the endpoint `p0(E)` but the exact two-far curvature of the approach.

All LRC quantities are exact Fractions.  Decimal values are printed only as
readable guides.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations
from math import gcd

from lrc14_wide_branch_ridge_codex_s47 import (
    CAP,
    additive_energy,
    has_ap,
    longest_run,
    p0,
    primitive,
    squarefree_profile,
    state_entropy,
    sumset_excess,
)


F = Fraction
Row = tuple[int, ...]


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def row_gcd(row: Row) -> int:
    nonzero = [abs(x) for x in row if x]
    return reduce(gcd, nonzero) if nonzero else 0


def scale_down(row: Row) -> Row:
    g = row_gcd(row)
    if g <= 1:
        return row
    return tuple(x // g for x in row)


def sumset_size(row: Row) -> int:
    return len({a + b for a in row for b in row})


def k2(row: Row) -> F:
    return F(sumset_size(row), len(row))


def longest_ap(row: Row) -> int:
    present = set(row)
    best = 1
    span = row[-1] - row[0]
    for a in row:
        for d in range(1, span + 1):
            length = 0
            x = a
            while x in present:
                length += 1
                x += d
            best = max(best, length)
    return best


def bucket_excess(exc: int) -> str:
    if exc == 0:
        return "AP"
    if exc <= 4:
        return "near_AP"
    if exc <= 10:
        return "small"
    if exc <= 18:
        return "medium"
    return "large"


@lru_cache(maxsize=None)
def p0_cached(row: Row) -> F:
    return p0(tuple(sorted(row)))


@dataclass(frozen=True)
class CurvReport:
    core: Row
    u: int
    v: int
    p_core: F
    p_u: F
    p_v: F
    p_uv: F

    @property
    def row(self) -> Row:
        return tuple(sorted(self.core + (self.u, self.v)))

    @property
    def k(self) -> int:
        return len(self.row)

    @property
    def cap(self) -> F:
        return CAP[self.k]

    @property
    def margin(self) -> F:
        return self.cap - self.p_uv

    @property
    def risk_ratio(self) -> F:
        return self.p_uv / self.cap

    @property
    def du(self) -> F:
        return self.p_u - self.p_core

    @property
    def dv(self) -> F:
        return self.p_v - self.p_core

    @property
    def dv_after_u(self) -> F:
        return self.p_uv - self.p_u

    @property
    def du_after_v(self) -> F:
        return self.p_uv - self.p_v

    @property
    def curvature(self) -> F:
        return self.p_uv - self.p_u - self.p_v + self.p_core


def make_report(core: Row, u: int, v: int) -> CurvReport:
    core = tuple(sorted(core))
    if u > v:
        u, v = v, u
    return CurvReport(
        core=core,
        u=u,
        v=v,
        p_core=p0_cached(core),
        p_u=p0_cached(tuple(sorted(core + (u,)))),
        p_v=p0_cached(tuple(sorted(core + (v,)))),
        p_uv=p0_cached(tuple(sorted(core + (u, v)))),
    )


def push_top(store: list[CurvReport], rep: CurvReport, key, keep: int) -> None:
    store.append(rep)
    store.sort(key=key, reverse=True)
    del store[keep:]


def structure_line(rep: CurvReport) -> str:
    row = rep.row
    core_scaled = scale_down(rep.core)
    support, entropy = state_entropy(row)
    return (
        f"core_gcd={row_gcd(rep.core)} core_scaled={core_scaled} "
        f"row_exc={sumset_excess(row)} core_exc={sumset_excess(rep.core)} "
        f"K2={k2(row)} energy={additive_energy(row)} run={longest_run(row)} "
        f"ap7={has_ap(row, 7)} lap={longest_ap(row)} "
        f"states={support} H={entropy:.4f} sqfree={squarefree_profile(row)}"
    )


def print_report(label: str, rep: CurvReport) -> None:
    print(f"{label}: E={rep.row}, core={rep.core}, far=({rep.u},{rep.v})")
    print(f"  p0(core)={fmt(rep.p_core)}")
    print(f"  p0(core+{rep.u})={fmt(rep.p_u)}  increment={fmt(rep.du)}")
    print(f"  p0(core+{rep.v})={fmt(rep.p_v)}  increment={fmt(rep.dv)}")
    print(f"  p0(full)={fmt(rep.p_uv)} cap={fmt(rep.cap)} margin={fmt(rep.margin)} risk={fmt(rep.risk_ratio)}")
    print(f"  curvature I_B({rep.u},{rep.v})={fmt(rep.curvature)}")
    print(f"  ordered increments: add {rep.v} after {rep.u} -> {fmt(rep.dv_after_u)}; add {rep.u} after {rep.v} -> {fmt(rep.du_after_v)}")
    print(f"  {structure_line(rep)}")


def scan_two_far_k9(far_max: int = 24, keep: int = 16) -> tuple[dict[str, list[CurvReport]], dict[str, object]]:
    tops: dict[str, list[CurvReport]] = {
        "risk": [],
        "positive_curvature": [],
        "negative_curvature_abs": [],
        "abs_curvature": [],
    }
    stats = {
        "raw": 0,
        "primitive": 0,
        "curv_sign": Counter(),
        "by_bucket": defaultdict(lambda: {"count": 0, "max_p0": F(0), "argmax": None, "max_curv": F(-10**9), "argcurv": None}),
        "by_core_gcd": defaultdict(lambda: {"count": 0, "max_p0": F(0), "argmax": None}),
        "far_pair_pressure": defaultdict(F),
        "far_pair_count": Counter(),
    }

    for core_rest in combinations(range(1, 15), 6):
        core = (0,) + core_rest
        for u, v in combinations(range(15, far_max + 1), 2):
            row = tuple(sorted(core + (u, v)))
            stats["raw"] += 1
            if not primitive(row):
                continue
            stats["primitive"] += 1
            rep = make_report(core, u, v)
            curv = rep.curvature
            stats["curv_sign"]["positive" if curv > 0 else "negative" if curv < 0 else "zero"] += 1
            bucket = bucket_excess(sumset_excess(row))
            b = stats["by_bucket"][bucket]
            b["count"] += 1
            if rep.p_uv > b["max_p0"]:
                b["max_p0"] = rep.p_uv
                b["argmax"] = rep
            if rep.curvature > b["max_curv"]:
                b["max_curv"] = rep.curvature
                b["argcurv"] = rep
            g = row_gcd(core)
            bg = stats["by_core_gcd"][g]
            bg["count"] += 1
            if rep.p_uv > bg["max_p0"]:
                bg["max_p0"] = rep.p_uv
                bg["argmax"] = rep
            stats["far_pair_pressure"][(u, v)] += rep.du - rep.dv
            stats["far_pair_count"][(u, v)] += 1

            push_top(tops["risk"], rep, key=lambda r: (r.p_uv, r.curvature, -sumset_excess(r.row), r.row), keep=keep)
            push_top(tops["positive_curvature"], rep, key=lambda r: (r.curvature, r.p_uv, r.row), keep=keep)
            push_top(tops["negative_curvature_abs"], rep, key=lambda r: (-r.curvature, r.p_uv, r.row), keep=keep)
            push_top(tops["abs_curvature"], rep, key=lambda r: (abs(r.curvature), r.p_uv, r.row), keep=keep)
    return tops, stats


def row_tournament(reports: list[CurvReport]) -> None:
    verts = reports[:12]
    score = [0] * len(verts)
    wins: set[tuple[int, int]] = set()
    for i, a in enumerate(verts):
        for j, b in enumerate(verts):
            if i >= j:
                continue
            ka = (a.risk_ratio, a.curvature, -sumset_excess(a.row), tuple(-x for x in a.row))
            kb = (b.risk_ratio, b.curvature, -sumset_excess(b.row), tuple(-x for x in b.row))
            if ka >= kb:
                score[i] += 1
                wins.add((i, j))
            else:
                score[j] += 1
                wins.add((j, i))
    cycles = 0
    for i, j, k in combinations(range(len(verts)), 3):
        if ((i, j) in wins and (j, k) in wins and (k, i) in wins) or (
            (i, k) in wins and (k, j) in wins and (j, i) in wins
        ):
            cycles += 1
    path = sorted(range(len(verts)), key=lambda i: (-score[i], -verts[i].risk_ratio, -verts[i].curvature, verts[i].row))
    print("Row-risk Tournament Analysis")
    print("  vertices: true-wide two-far rows / proof obligations, not runners")
    print("  pairwise observable: larger exact p0/cap wins; ties use curvature then additive excess")
    print("  switch/gauge: direct sector-cover risk; tie Hamiltonian path is lexicographic")
    print(f"  score_hist={dict(sorted(Counter(score).items()))}, directed_3cycles={cycles}")
    print("  Hamiltonian path:")
    for rank, idx in enumerate(path, 1):
        rep = verts[idx]
        print(
            f"    {rank:2d}. E={rep.row} risk={rep.risk_ratio} "
            f"p0={rep.p_uv} curv={rep.curvature} exc={sumset_excess(rep.row)}"
        )


def far_speed_tournament(stats: dict[str, object], far_min: int = 15, far_max: int = 24) -> None:
    pressure: dict[tuple[int, int], F] = stats["far_pair_pressure"]  # type: ignore[assignment]
    speeds = list(range(far_min, far_max + 1))
    score = {v: 0 for v in speeds}
    wins: set[tuple[int, int]] = set()
    ties = 0
    for u, v in combinations(speeds, 2):
        p = pressure[(u, v)]
        if p > 0 or (p == 0 and u < v):
            score[u] += 1
            wins.add((u, v))
            if p == 0:
                ties += 1
        else:
            score[v] += 1
            wins.add((v, u))
    cycles = 0
    for a, b, c in combinations(speeds, 3):
        if ((a, b) in wins and (b, c) in wins and (c, a) in wins) or (
            (a, c) in wins and (c, b) in wins and (b, a) in wins
        ):
            cycles += 1
    path = sorted(speeds, key=lambda x: (-score[x], x))
    print("Far-speed first-increment Tournament Analysis")
    print("  vertices: far speeds 15..24")
    print("  pairwise observable: aggregate over cores of [p0(B+u)-p0(B)] - [p0(B+v)-p0(B)]")
    print("  switch/gauge: larger first increment wins; lower speed wins exact zero ties")
    print(f"  score_hist={dict(sorted(Counter(score.values()).items()))}, directed_3cycles={cycles}, ties={ties}")
    print(f"  Hamiltonian path={path}")
    for u, v in list(combinations(speeds, 2))[:8]:
        print(f"    edge {u}->{v} pressure={pressure[(u, v)]}")


def path_report(label: str, row: Row, core_max: int = 14) -> None:
    core = tuple(x for x in row if x <= core_max)
    fars = tuple(x for x in row if x > core_max)
    print(f"{label}: row={row}")
    print(f"  core<={core_max}: {core}, fars={fars}")
    current = core
    prev = p0_cached(current)
    print(f"  p0(core)={fmt(prev)}")
    for far in fars:
        nxt = tuple(sorted(current + (far,)))
        value = p0_cached(nxt)
        print(f"  +{far}: p0={fmt(value)} increment={fmt(value - prev)}")
        current = nxt
        prev = value
    if len(fars) >= 2:
        for u, v in combinations(fars, 2):
            rep = make_report(core, u, v)
            print(f"  pair curvature over base core ({u},{v}) = {fmt(rep.curvature)}")
    print()


def main() -> None:
    print("HYP-2679 / T918 LRC14 true-wide boundary-curvature scout")
    print("Arithmetic: exact Fractions; decimal displays are guides only.")
    print("Assumption challenge: vertices below are rows, far speeds, or approach orders, not runners.")
    print("Preserved predicate: direct p0(E)<=cap_k. Destroyed data: internal wall endpoints unless state words are reattached.")
    print()

    print("Named calibration rows")
    named = [
        ("HYP2675 true-wide leader", make_report((0, 4, 6, 8, 10, 12, 14), 15, 16)),
        ("k8 true-wide leader", make_report((0, 3, 6, 9, 12, 14), 15, 18)),
        ("even AP core plus two extras", make_report((0, 2, 4, 6, 8, 10, 12), 14, 15)),
    ]
    for label, rep in named:
        print_report(label, rep)
    print()

    print("Multi-far approach paths")
    path_report("HYP2675 true-wide leader sorted path", (0, 4, 6, 8, 10, 12, 14, 15, 16))
    path_report("KPS third-pocket sorted path", (0, 3, 5, 16, 28, 30, 33, 35))

    far_max = 24
    print(f"Exact k=9 two-far scan: core=(0)+6-subsets of [1,14], far pairs 15..{far_max}")
    tops, stats = scan_two_far_k9(far_max=far_max, keep=16)
    print(f"  raw_rows={stats['raw']}, primitive_rows={stats['primitive']}")
    print(f"  curvature signs={dict(stats['curv_sign'])}")
    cache_info = p0_cached.cache_info()
    print(f"  p0 cache: hits={cache_info.hits}, misses={cache_info.misses}, size={cache_info.currsize}")
    print()

    for title, key in [
        ("Top direct-risk rows", "risk"),
        ("Top positive-curvature rows", "positive_curvature"),
        ("Top negative-curvature-overlap rows", "negative_curvature_abs"),
        ("Top absolute-curvature rows", "abs_curvature"),
    ]:
        print(title)
        for idx, rep in enumerate(tops[key][:10], 1):
            print(
                f"  {idx:2d}. E={rep.row} p0={rep.p_uv} margin={rep.margin} "
                f"curv={rep.curvature} du={rep.du} dv={rep.dv} "
                f"exc={sumset_excess(rep.row)} core_gcd={row_gcd(rep.core)}"
            )
        print()

    print("Bucket maxima")
    for bucket in ["AP", "near_AP", "small", "medium", "large"]:
        b = stats["by_bucket"][bucket]
        if not b["count"]:
            continue
        arg = b["argmax"]
        carg = b["argcurv"]
        print(
            f"  {bucket}: count={b['count']}, max_p0={b['max_p0']} at {arg.row if arg else None}, "
            f"max_curv={b['max_curv']} at {carg.row if carg else None}"
        )
    print()

    print("Core-gcd maxima")
    for g in sorted(stats["by_core_gcd"]):
        b = stats["by_core_gcd"][g]
        arg = b["argmax"]
        print(f"  gcd={g}: count={b['count']}, max_p0={b['max_p0']} at {arg.row if arg else None}")
    print()

    row_tournament(tops["risk"])
    print()
    far_speed_tournament(stats, far_min=15, far_max=far_max)
    print()

    leader = tops["risk"][0]
    print("Proof-route reading")
    print("  1. The exact two-far curvature is the finite analogue of curvilinear boundary approach data.")
    print("  2. Positive curvature means the far speeds are complementary; negative curvature means one-far overlap.")
    print("  3. If the risk leaders are confined to low-excess/dilated-core rows, HYP-2678's d=1 finite check is the right next theorem.")
    print("  4. If high-excess rows only appear with small or negative curvature, the d>=2 branch should be attacked by a signed dimension-penalty bound.")
    print(
        f"  5. Scan leader: E={leader.row}, p0={leader.p_uv}, margin={leader.margin}, "
        f"curvature={leader.curvature}, core_scaled={scale_down(leader.core)}."
    )
    print("PASS: exact boundary-curvature scout complete; no LRC(14) proof claimed.")


if __name__ == "__main__":
    main()
