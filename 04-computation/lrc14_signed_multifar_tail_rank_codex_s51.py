#!/usr/bin/env python3
"""HYP-2680: signed multi-far tail and bounded relation-rank scout.

The previous S51 script isolates the individual three-far residual

    Delta_S(B) - Phi_|S|(B).

This scout uses the simultaneous-peel architecture of THM-548.  For a bounded
core B and a whole far block F it decomposes

    p0(B union F) - P_|F|(B)

by Newton order.  The point is to see whether positive two-far mass is cancelled
by signed three-/four-far mass, and whether the obstruction is measured by
bounded low-height relations among the far runners.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, gcd

from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import (
    boundary_value_direct,
    delta_s,
    phi_s,
    p0_cached,
)
from lrc14_wide_branch_ridge_codex_s47 import (
    CAP,
    additive_energy,
    primitive,
    squarefree_profile,
    sumset_excess,
)


Row = tuple[int, ...]


def fmt(q: F) -> str:
    return f"{q} ({float(q):+.9f})"


def cap_for(k: int) -> F | None:
    if k == 13:
        return F(1)
    return CAP.get(k)


def sign(q: F) -> str:
    if q > 0:
        return "+"
    if q < 0:
        return "-"
    return "0"


def primitive_coeff(coeff: tuple[int, ...]) -> bool:
    g = 0
    for c in coeff:
        g = gcd(g, abs(c))
    return g == 1


def rank_vectors(vectors: list[tuple[int, ...]], width: int) -> int:
    rows = [[F(x) for x in v] for v in vectors if any(v)]
    rank = 0
    col = 0
    while rank < len(rows) and col < width:
        pivot = None
        for r in range(rank, len(rows)):
            if rows[r][col]:
                pivot = r
                break
        if pivot is None:
            col += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pv = rows[rank][col]
        rows[rank] = [x / pv for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][col]:
                factor = rows[r][col]
                rows[r] = [a - factor * b for a, b in zip(rows[r], rows[rank])]
        rank += 1
        col += 1
    return rank


def relation_profile(vals: Row, height: int = 3) -> dict[str, object]:
    exact: list[tuple[int, ...]] = []
    near_counts: Counter[int] = Counter()
    best: tuple[int, int, tuple[int, ...]] | None = None

    for coeff in product(range(-height, height + 1), repeat=len(vals)):
        if all(c == 0 for c in coeff) or not primitive_coeff(coeff):
            continue
        dot = sum(c * v for c, v in zip(coeff, vals))
        ad = abs(dot)
        l1 = sum(abs(c) for c in coeff)
        if ad <= 3:
            near_counts[ad] += 1
        if ad == 0:
            exact.append(coeff)
        item = (ad, l1, coeff)
        if best is None or item < best:
            best = item

    basis: list[tuple[int, ...]] = []
    for coeff in sorted(exact, key=lambda c: (sum(abs(x) for x in c), c)):
        if rank_vectors(basis + [coeff], len(vals)) > rank_vectors(basis, len(vals)):
            basis.append(coeff)

    assert best is not None
    return {
        "height": height,
        "best": best,
        "exact_count": len(exact),
        "exact_rank": rank_vectors(exact, len(vals)),
        "basis": basis,
        "near_counts": dict(sorted(near_counts.items())),
    }


@dataclass(frozen=True)
class OrderStats:
    order: int
    count: int
    total: F
    abs_total: F
    max_abs: F
    max_subset: Row
    positive: int
    negative: int
    zero: int

    @property
    def average(self) -> F:
        return self.total / self.count if self.count else F(0)

    @property
    def scaled(self) -> F:
        return self.total * (7 ** (self.order + 1))


def order_stats(core: Row, far: Row, order: int) -> OrderStats:
    total = F(0)
    abs_total = F(0)
    max_abs = F(-1)
    max_subset: Row = ()
    signs = Counter()
    count = 0
    for sub in combinations(far, order):
        dev = delta_s(core, sub) - phi_s(core, order)
        total += dev
        abs_total += abs(dev)
        signs[sign(dev)] += 1
        count += 1
        if abs(dev) > max_abs:
            max_abs = abs(dev)
            max_subset = sub
    return OrderStats(
        order=order,
        count=count,
        total=total,
        abs_total=abs_total,
        max_abs=max_abs,
        max_subset=max_subset,
        positive=signs["+"],
        negative=signs["-"],
        zero=signs["0"],
    )


@dataclass(frozen=True)
class BlockReport:
    label: str
    core: Row
    far: Row
    direct: F
    boundary: F
    orders: tuple[OrderStats, ...]

    @property
    def row(self) -> Row:
        return tuple(sorted(self.core + self.far))

    @property
    def residual(self) -> F:
        return self.direct - self.boundary

    @property
    def cap(self) -> F | None:
        return cap_for(len(self.row))

    @property
    def margin(self) -> F | None:
        cap = self.cap
        return None if cap is None else cap - self.direct

    @property
    def sign_word(self) -> str:
        return "".join(sign(o.total) for o in self.orders)


def make_block(label: str, core: Row, far: Row) -> BlockReport:
    core = tuple(sorted(core))
    far = tuple(sorted(far))
    direct = p0_cached(tuple(sorted(core + far)))
    boundary = boundary_value_direct(core, len(far))
    orders = tuple(order_stats(core, far, s) for s in range(1, len(far) + 1))
    assert direct == boundary + sum((o.total for o in orders), F(0))
    return BlockReport(label, core, far, direct, boundary, orders)


def print_block(rep: BlockReport) -> None:
    rel = relation_profile(rep.far, 3)
    cap_text = "n/a" if rep.cap is None else fmt(rep.cap)
    margin_text = "n/a" if rep.margin is None else fmt(rep.margin)
    print(f"CASE {rep.label}")
    print(f"  core={rep.core}")
    print(f"  far={rep.far} row={rep.row}")
    print(
        f"  relation_height={rel['height']} best={rel['best']} "
        f"exact_rank={rel['exact_rank']} exact_count={rel['exact_count']} "
        f"basis={rel['basis']} near_counts={rel['near_counts']}"
    )
    print(
        f"  row_exc={sumset_excess(rep.row)} energy={additive_energy(rep.row)} "
        f"primitive={primitive(rep.row)} sqfree={squarefree_profile(rep.row)}"
    )
    print(f"  P_r={fmt(rep.boundary)}")
    print(f"  p0={fmt(rep.direct)} cap={cap_text} margin={margin_text}")
    print(f"  total_residual=p0-P_r={fmt(rep.residual)} order_sign_word={rep.sign_word}")
    print("  order residuals:")
    for o in rep.orders:
        print(
            f"    s={o.order} count={o.count} sum={fmt(o.total)} "
            f"avg={fmt(o.average)} abs_sum={fmt(o.abs_total)} "
            f"max_abs={fmt(o.max_abs)} at={o.max_subset} "
            f"signs=+{o.positive}/-{o.negative}/0{o.zero} "
            f"apex_scaled=7^{o.order + 1}*sum={fmt(o.scaled)}"
        )
    tournament_analysis(rep)
    print()


def tournament_analysis(rep: BlockReport) -> None:
    """Tournament on Newton orders as proof-obligation vertices."""

    orders = list(rep.orders)
    wins = [0] * len(orders)

    def key(o: OrderStats) -> tuple[F, F, int]:
        positive_stress = o.total if o.total > 0 else F(0)
        return positive_stress, abs(o.total), -o.order

    def beats(a: OrderStats, b: OrderStats) -> bool:
        return key(a) > key(b)

    cycles = 0
    for i, j in combinations(range(len(orders)), 2):
        if beats(orders[i], orders[j]):
            wins[i] += 1
        else:
            wins[j] += 1
    for i, j, k in combinations(range(len(orders)), 3):
        eij = beats(orders[i], orders[j])
        ejk = beats(orders[j], orders[k])
        eki = beats(orders[k], orders[i])
        if (eij and ejk and eki) or ((not eij) and (not ejk) and (not eki)):
            cycles += 1

    order = sorted(range(len(orders)), key=lambda i: wins[i], reverse=True)
    print("  tournament(vertices=Newton orders): observable=larger positive cap-stress, tie by |sum|")
    print(f"    score_hist={dict(Counter(wins))} directed_3cycles={cycles}")
    print(
        "    Hamiltonian_path="
        + " > ".join(f"R{orders[i].order}({sign(orders[i].total)}{abs(orders[i].total)})" for i in order)
    )


def named_blocks() -> list[BlockReport]:
    cases = [
        ("dilated true-wide, three consecutive far", (0, 4, 6, 8, 10, 12, 14), (15, 16, 17)),
        ("dilated true-wide, four consecutive far", (0, 4, 6, 8, 10, 12, 14), (15, 16, 17, 18)),
        ("dilated true-wide, six consecutive far", (0, 4, 6, 8, 10, 12, 14), (15, 16, 17, 18, 19, 20)),
        ("dilated true-wide, three separated far", (0, 4, 6, 8, 10, 12, 14), (17, 23, 31)),
        ("dilated true-wide, five separated far", (0, 4, 6, 8, 10, 12, 14), (17, 23, 31, 43, 59)),
        ("consec8, four consecutive far", (0, 1, 2, 3, 4, 5, 6, 7), (15, 16, 17, 18)),
        ("consec8, six consecutive far", (0, 1, 2, 3, 4, 5, 6, 7), (15, 16, 17, 18, 19, 20)),
        ("third-pocket active, three far", (0, 3, 5, 16, 28), (30, 33, 35)),
        ("Newton high-order witness", (0, 1, 2, 3), (15, 16, 17, 18, 19, 20, 21)),
    ]
    reports: list[BlockReport] = []
    for label, core, far in cases:
        rep = make_block(label, core, far)
        reports.append(rep)
        print_block(rep)
    return reports


def small_four_far_bank() -> None:
    print("SMALL FOUR-FAR BANK: all primitive 7-cores in [0,14], far=(15,16,17,18)")
    far = (15, 16, 17, 18)
    stats: Counter[str] = Counter()
    top_direct: list[BlockReport] = []
    top_total: list[BlockReport] = []
    top_r4: list[BlockReport] = []
    primitive_count = 0

    def push(store: list[BlockReport], rep: BlockReport, key, keep: int = 8) -> None:
        store.append(rep)
        store.sort(key=key, reverse=True)
        del store[keep:]

    for rest in combinations(range(1, 15), 6):
        core = (0,) + rest
        row = tuple(sorted(core + far))
        if not primitive(row):
            continue
        primitive_count += 1
        rep = make_block("bank", core, far)
        r2, r3, r4 = rep.orders[1], rep.orders[2], rep.orders[3]
        stats["sign_" + sign(r2.total) + sign(r3.total) + sign(r4.total)] += 1
        stats["R2_R3_opposite"] += int(sign(r2.total) != "0" and sign(r3.total) != "0" and sign(r2.total) != sign(r3.total))
        stats["R3_R4_opposite"] += int(sign(r3.total) != "0" and sign(r4.total) != "0" and sign(r3.total) != sign(r4.total))
        stats["total_positive"] += int(rep.residual > 0)
        stats["total_negative"] += int(rep.residual < 0)
        push(top_direct, rep, lambda r: (r.direct, r.residual, r.row))
        push(top_total, rep, lambda r: (abs(r.residual), r.direct, r.row))
        push(top_r4, rep, lambda r: (abs(r.orders[3].total), r.direct, r.row))

    print(f"  primitive_count={primitive_count}")
    print(f"  stats={dict(stats)}")
    print("  top direct p0:")
    for rep in top_direct[:6]:
        print(
            f"    row={rep.row} p0={rep.direct} margin={rep.margin} "
            f"P4={rep.boundary} residual={rep.residual} sign_word={rep.sign_word}"
        )
    print("  top |total residual|:")
    for rep in top_total[:6]:
        print(
            f"    row={rep.row} residual={rep.residual} p0={rep.direct} "
            f"margin={rep.margin} sign_word={rep.sign_word}"
        )
    print("  top |R4|:")
    for rep in top_r4[:6]:
        r4 = rep.orders[3]
        print(
            f"    row={rep.row} R4={r4.total} abs={abs(r4.total)} "
            f"p0={rep.direct} margin={rep.margin} sign_word={rep.sign_word}"
        )
    print()


def synthesis(reports: list[BlockReport]) -> None:
    print("SYNTHESIS")
    print("  Simultaneous peel should be read as a signed order ledger, not as an absolute tail.")
    print("  Consecutive/AP far blocks have low-height relation rank and large individual order terms,")
    print("  but their order signs often oppose; separated far blocks are already close to P_r.")
    print("  The sharp proof obligation is a rank-sensitive packet estimate for each order sum")
    print("  Sigma_{|S|=s}(Delta_S-Phi_s), with exact low-height relations routed to finite scale atlases.")
    print("  High orders do not vanish; the useful claim is geometric apex-prime suppression plus signs.")
    print()
    print("  Residual leaderboard:")
    for rep in sorted(reports, key=lambda r: (abs(r.residual), r.direct), reverse=True)[:8]:
        rel = relation_profile(rep.far, 3)
        print(
            f"    {rep.label}: residual={rep.residual} p0={rep.direct} "
            f"orders={rep.sign_word} exact_rank={rel['exact_rank']} far={rep.far}"
        )


def main() -> None:
    print("HYP-2680 / S51 -- signed multi-far tail rank scout\n")
    reports = named_blocks()
    small_four_far_bank()
    synthesis(reports)


if __name__ == "__main__":
    main()
