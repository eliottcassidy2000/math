#!/usr/bin/env python3
"""HYP-2708 / T941: live-depth kernel for two-far survival deviation.

HYP-2701 reduced the true-wide two-far problem to:

    actual survival currency >= decorrelated death-chain boundary - margin.

This script decomposes the actual-minus-boundary deviation at the exact
wall-atom level.  The useful surprise is algebraic: for the currency

    C = p1+p2+p3+p4-4p6,

and two far hits, before-depths 3 and 4 contribute no deviation at all.
Only depths 1, 2, 5, and 6 can spend the death-chain boundary margin.

The script verifies that identity exactly and audits the tight two-far rows.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, gcd

from lrc14_twofar_survival_currency_boundary_codex_s64 import (
    Row,
    boundary_currency,
    currency,
    floor_required,
    primitive,
    profile,
    relation_distance,
)
from lrc14_wide_branch_ridge_codex_s47 import (
    additive_energy,
    fmt,
    squarefree_profile,
    sumset_excess,
    wall_breakpoints,
)


INNER_BITS = sum(1 << j for j in range(1, 7))
LIVE_DEPTHS = (1, 2, 5, 6)


def currency_coeff(t: int) -> F:
    if 1 <= t <= 4:
        return F(1)
    if t == 6:
        return F(-4)
    return F(0)


def transition_prob(t: int, s: int, r: int = 2) -> F:
    if not 0 <= s <= t <= 6:
        return F(0)
    need_hit = t - s
    total = F(0)
    for j in range(need_hit + 1):
        total += ((-1) ** j) * comb(need_hit, j) * F(7 - s - j, 7) ** r
    return F(comb(t, s)) * total


def boundary_after_two(t: int) -> F:
    return sum(transition_prob(t, s, 2) * currency_coeff(s) for s in range(t + 1))


def hit_boundary_prob(t: int, h: int) -> F:
    """Boundary probability that two iid colors hit h distinct missed colors."""

    if h == 0:
        return F((7 - t) ** 2, 49)
    if h == 1:
        return F(t * (2 * (7 - t) + 1), 49)
    if h == 2:
        return F(t * (t - 1), 49)
    return F(0)


def sector_at(speed: int, midnum: int, den2: int) -> int:
    return (speed * midnum // den2) % 7


@dataclass(frozen=True)
class KernelAudit:
    row: Row
    core: Row
    far: Row
    mass_by_th: dict[tuple[int, int], F]
    dev_by_t: dict[int, F]
    actual: F
    boundary: F
    required: F

    @property
    def k(self) -> int:
        return len(self.row)

    @property
    def deviation(self) -> F:
        return self.actual - self.boundary

    @property
    def slack(self) -> F:
        return self.actual - self.required

    @property
    def boundary_margin(self) -> F:
        return self.boundary - self.required

    @property
    def relation(self) -> tuple[int, tuple[int, ...], int]:
        return relation_distance(self.far)

    @property
    def live_mass(self) -> F:
        return sum(mass for (t, _h), mass in self.mass_by_th.items() if t in LIVE_DEPTHS)


def audit_row(row: Row) -> KernelAudit:
    row = tuple(sorted(row))
    core = tuple(x for x in row if x <= 14)
    far = tuple(x for x in row if x > 14)
    assert len(far) == 2

    d, bps = wall_breakpoints(row)
    l = d // 7
    den2 = 2 * l
    mass_by_th: dict[tuple[int, int], F] = defaultdict(F)
    actual_sum = F(0)
    boundary_sum = F(0)
    dev_by_t: dict[int, F] = defaultdict(F)

    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mass = F(hi - lo, d)
        midnum = lo + hi
        core_mask = 0
        for speed in core:
            if speed:
                core_mask |= 1 << sector_at(speed, midnum, den2)
        missed = INNER_BITS & ~core_mask
        t = missed.bit_count()
        hit = 0
        for speed in far:
            sec = sector_at(speed, midnum, den2)
            if 1 <= sec <= 6:
                hit |= 1 << sec
        h = (hit & missed).bit_count()
        after_t = t - h
        actual_c = currency_coeff(after_t)
        boundary_c = boundary_after_two(t)
        mass_by_th[(t, h)] += mass
        actual_sum += mass * actual_c
        boundary_sum += mass * boundary_c
        dev_by_t[t] += mass * (actual_c - boundary_c)

    assert actual_sum == currency(row)
    assert boundary_sum == boundary_currency(core, 2)
    assert sum(dev_by_t.values(), F(0)) == actual_sum - boundary_sum
    return KernelAudit(
        row=row,
        core=core,
        far=far,
        mass_by_th=dict(mass_by_th),
        dev_by_t=dict(dev_by_t),
        actual=actual_sum,
        boundary=boundary_sum,
        required=floor_required(len(row)),
    )


def live_depth_formula_checks() -> None:
    print("LIVE-DEPTH FORMULA CHECK")
    print("  h = number of distinct before-missed sectors hit by the two far colors")
    print("  dev_t = E_actual[C(t-h)] - E_iid[C(t-h)]")
    for t in range(7):
        iid = sum(hit_boundary_prob(t, h) * currency_coeff(t - h) for h in range(3))
        assert iid == boundary_after_two(t)
        coeffs = [currency_coeff(t - h) - boundary_after_two(t) for h in range(3)]
        live = any(c != 0 for c in coeffs)
        print(
            f"  t={t}: K2={str(boundary_after_two(t)):>8} "
            f"coeffs(h0,h1,h2)={[str(c) for c in coeffs]} "
            f"{'LIVE' if live else 'silent'}"
        )
    print("  conclusion: only t=1,2,5,6 can contribute; t=3,4 are exactly silent.")
    print()


def row_line(label: str, audit: KernelAudit) -> None:
    rel_val, rel_coeff, rel_l1 = audit.relation
    prof = profile(audit.row)
    print(f"{label}: E={audit.row}, core={audit.core}, far={audit.far}")
    print(
        f"  actual={fmt(audit.actual)} boundary={fmt(audit.boundary)} "
        f"required={fmt(audit.required)} slack={fmt(audit.slack)}"
    )
    print(
        f"  deviation={fmt(audit.deviation)} boundary_margin={fmt(audit.boundary_margin)} "
        f"live_mass={fmt(audit.live_mass)} relation={rel_val} coeff={rel_coeff} l1={rel_l1}"
    )
    print(
        f"  dist={[str(x) for x in prof]} exc={sumset_excess(audit.row)} "
        f"energy={additive_energy(audit.row)} sqfree={squarefree_profile(audit.row)}"
    )
    print("  deviation by before-depth:")
    for t in range(7):
        if audit.dev_by_t.get(t, F(0)):
            print(f"    t={t}: {fmt(audit.dev_by_t[t])}")
    print("  actual hit kernel mass by (t,h), only nonzero:")
    for key in sorted(audit.mass_by_th):
        t, h = key
        if t not in LIVE_DEPTHS:
            continue
        mass = audit.mass_by_th[key]
        if mass:
            iid = profile(audit.core)[t] * hit_boundary_prob(t, h)
            print(
                f"    t={t}, h={h}: actual={fmt(mass)} "
                f"iid_boundary={fmt(iid)} diff={fmt(mass - iid)}"
            )
    print()


def scan_twofar(k: int, bound: int, keep: int = 8) -> list[KernelAudit]:
    rows: list[KernelAudit] = []
    for combo in combinations(range(1, bound + 1), k - 1):
        row = (0,) + combo
        if not primitive(row):
            continue
        if sum(1 for x in row if x > 14) != 2:
            continue
        audit = audit_row(row)
        rows.append(audit)
    rows.sort(key=lambda a: (a.slack, a.deviation, a.row))
    print(f"SCAN k={k}, bound={bound}: two-far rows={len(rows)}")
    print("  tightest rows:")
    for i, audit in enumerate(rows[:keep], 1):
        print(
            f"    {i}: slack={audit.slack} dev={audit.deviation} "
            f"live_mass={audit.live_mass} row={audit.row}"
        )
    depth_totals: dict[int, F] = defaultdict(F)
    dev_totals: dict[int, F] = defaultdict(F)
    for audit in rows[: min(250, len(rows))]:
        for (t, _h), mass in audit.mass_by_th.items():
            depth_totals[t] += mass
        for t, dev in audit.dev_by_t.items():
            dev_totals[t] += dev
    print("  aggregate over top min(250, rows) by actual slack:")
    for t in range(7):
        if depth_totals.get(t, F(0)) or dev_totals.get(t, F(0)):
            print(
                f"    t={t}: mass={fmt(depth_totals.get(t, F(0)))} "
                f"dev={fmt(dev_totals.get(t, F(0)))}"
            )
    print()
    return rows


def tournament(rows: list[KernelAudit]) -> None:
    verts = rows[:10]
    edges: set[tuple[int, int]] = set()
    scores = [0] * len(verts)
    for i, a in enumerate(verts):
        for j, b in enumerate(verts):
            if i >= j:
                continue
            ka = (a.slack, a.deviation, -a.live_mass, a.row)
            kb = (b.slack, b.deviation, -b.live_mass, b.row)
            if ka <= kb:
                edges.add((i, j))
                scores[i] += 1
            else:
                edges.add((j, i))
                scores[j] += 1
    cycles = 0
    for a, b, c in combinations(range(len(verts)), 3):
        if (a, b) in edges and (b, c) in edges and (c, a) in edges:
            cycles += 1
        if (a, c) in edges and (c, b) in edges and (b, a) in edges:
            cycles += 1
    print("TOURNAMENT ANALYSIS")
    print("  vertices: two-far row obligations")
    print("  pairwise observable: smaller cap slack, then more negative live-depth deviation")
    print("  switch/gauge: live-depth kernel before scalar row invariants")
    print(f"  score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  Hamiltonian risk path:")
    for i, audit in enumerate(sorted(verts, key=lambda a: (a.slack, a.deviation, -a.live_mass, a.row)), 1):
        print(
            f"    {i}. k={audit.k} row={audit.row} slack={audit.slack} "
            f"dev={audit.deviation} live_mass={audit.live_mass}"
        )
    print()


def main() -> None:
    print("HYP-2708 / T941 -- two-far live-depth survival kernel")
    print("Exact arithmetic: Fraction wall atoms, imported S64/S47 engines.\n")
    live_depth_formula_checks()

    detailed = [
        ("k8 floor-fail", (0, 3, 6, 9, 12, 14, 15, 18)),
        ("k9 leader", (0, 4, 6, 8, 10, 12, 14, 15, 16)),
        ("k10 leader", (0, 2, 4, 6, 8, 10, 12, 14, 15, 16)),
    ]
    print("NAMED LIVE-DEPTH DECOMPOSITIONS")
    for label, row in detailed:
        row_line(label, audit_row(row))

    bank_rows = [
        (0, 3, 6, 9, 12, 14, 15, 18),
        (0, 3, 6, 7, 9, 12, 15, 18),
        (0, 2, 3, 6, 9, 12, 15, 18),
        (0, 4, 6, 8, 10, 12, 14, 15, 16),
        (0, 3, 5, 7, 9, 11, 13, 15, 17),
        (0, 3, 6, 9, 11, 12, 14, 15, 18),
        (0, 2, 4, 6, 8, 10, 12, 14, 15, 16),
        (0, 3, 5, 6, 9, 12, 14, 15, 18),
        (0, 9, 10, 11, 12, 13, 14, 15, 16),
    ]
    risk_rows = [audit_row(row) for row in bank_rows]
    print("COMPACT RISK BANK FROM S64 LEADERS")
    depth_mass: dict[int, F] = defaultdict(F)
    depth_dev: dict[int, F] = defaultdict(F)
    for audit in risk_rows:
        print(
            f"  k={audit.k} row={audit.row} slack={fmt(audit.slack)} "
            f"dev={fmt(audit.deviation)} live_mass={fmt(audit.live_mass)}"
        )
        for (t, _h), mass in audit.mass_by_th.items():
            if t in LIVE_DEPTHS:
                depth_mass[t] += mass
        for t, dev in audit.dev_by_t.items():
            if t in LIVE_DEPTHS:
                depth_dev[t] += dev
    print("  aggregate live-depth mass/deviation over risk bank:")
    for t in LIVE_DEPTHS:
        print(f"    t={t}: mass={fmt(depth_mass[t])} dev={fmt(depth_dev[t])}")
    print()
    tournament(sorted(risk_rows, key=lambda a: (a.slack, a.deviation, a.row)))

    print("SYNTHESIS")
    print("  The two-far actual-boundary deviation has only four live before-depths: 1,2,5,6.")
    print("  Depths 3 and 4 are exactly silent because C(t), C(t-1), and C(t-2) all equal 1.")
    print("  Therefore the two-far proof can ignore the bulk middle mass after the boundary step.")
    print("  The remaining analytic lemma is a signed four-depth kernel bound, not a seven-depth bound.")


if __name__ == "__main__":
    main()
