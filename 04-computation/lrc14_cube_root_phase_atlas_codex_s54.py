#!/usr/bin/env python3
"""HYP-2682 / T921: cube-root phase/support atlas for LRC(14).

S52 showed that the seven three-far packets

    A,B,C ; D,E,F ; G

should not be collapsed to a single scalar.  This scout keeps the exact
Eisenstein imbalance from S52, but asks a sharper finite-atlas question:

    if the far triple has the same rank-one relation u - 2v + w = 0,
    what phase/support data still changes the signs and the direct risk?

The answer is tested on consecutive triples (m,m+1,m+2), where relation rank is
fixed but the mod-7 phase and core support address move.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations

from lrc14_cube_root_order_filter_codex_s52 import TriplePackets, make_packets, sgn
from lrc14_signed_multifar_tail_rank_codex_s51 import relation_profile
from lrc14_wide_branch_ridge_codex_s47 import CAP, additive_energy, primitive, sumset_excess


Row = tuple[int, ...]


def fmt(q: F | None) -> str:
    if q is None:
        return "n/a"
    return f"{q} ({float(q):+.9f})"


def sign_tuple(pkt: TriplePackets) -> tuple[str, str, str, str]:
    """actual residual, pair-tax shadow, pair layer, triple packet."""

    return (sgn(pkt.total), sgn(pkt.recursion), sgn(pkt.pairs), sgn(pkt.g))


def residue_word(row: Row, mod: int = 7) -> tuple[int, ...]:
    return tuple(sorted(x % mod for x in row))


def residue_hist(row: Row, mod: int = 7) -> tuple[int, ...]:
    counts = [0] * mod
    for x in row:
        counts[x % mod] += 1
    return tuple(counts)


def chamber(pkt: TriplePackets) -> str:
    """A2/Weyl chamber of the Eisenstein imbalance a+b*omega.

    The chamber walls are a=0, b=0, and a-b=0, i.e. the six natural sectors in
    the triangular lattice.  Wall labels are kept because they are useful
    finite-atlas addresses rather than numerical annoyances.
    """

    z = pkt.imbalance
    return f"a{sgn(z.a)} b{sgn(z.b)} d{sgn(z.a - z.b)}"


@dataclass(frozen=True)
class AtlasRow:
    pkt: TriplePackets
    chamber: str
    signs: tuple[str, str, str, str]
    far_residue: tuple[int, ...]
    core_hist: tuple[int, ...]
    row_exc: int
    energy: int

    @property
    def row(self) -> Row:
        return self.pkt.row

    @property
    def p0(self) -> F:
        return self.pkt.p0

    @property
    def margin(self) -> F | None:
        cap = CAP.get(len(self.row))
        return None if cap is None else cap - self.p0


def make_row(core: Row, far: Row) -> AtlasRow:
    pkt = make_packets(core, far)
    return AtlasRow(
        pkt=pkt,
        chamber=chamber(pkt),
        signs=sign_tuple(pkt),
        far_residue=residue_word(far),
        core_hist=residue_hist(core),
        row_exc=sumset_excess(pkt.row),
        energy=additive_energy(pkt.row),
    )


def push_top(store: list[AtlasRow], row: AtlasRow, key, keep: int = 10) -> None:
    store.append(row)
    store.sort(key=key, reverse=True)
    del store[keep:]


def scan_consecutive_triples(m_start: int = 15, m_stop: int = 35) -> list[AtlasRow]:
    rows: list[AtlasRow] = []
    per_far: dict[Row, dict[str, object]] = {}

    for m in range(m_start, m_stop + 1):
        far = (m, m + 1, m + 2)
        rel = relation_profile(far, height=3)
        stats: Counter[str] = Counter()
        sign_stats: Counter[tuple[str, str, str, str]] = Counter()
        chamber_stats: Counter[str] = Counter()
        top_p0: list[AtlasRow] = []
        top_imbalance: list[AtlasRow] = []
        top_pair_tax: list[AtlasRow] = []

        for rest in combinations(range(1, 15), 6):
            core = (0,) + rest
            full = tuple(sorted(core + far))
            if not primitive(full):
                continue
            row = make_row(core, far)
            rows.append(row)
            stats["primitive"] += 1
            stats[f"actual_{row.signs[0]}"] += 1
            stats[f"pairtax_{row.signs[1]}"] += 1
            stats[f"same_actual_pairtax_{row.signs[0] == row.signs[1]}"] += 1
            sign_stats[row.signs] += 1
            chamber_stats[row.chamber] += 1
            push_top(top_p0, row, lambda r: (r.p0, abs(r.pkt.total), r.row))
            push_top(top_imbalance, row, lambda r: (r.pkt.imbalance.norm(), r.p0, r.row))
            push_top(top_pair_tax, row, lambda r: (abs(r.pkt.recursion), r.p0, r.row))

        per_far[far] = {
            "relation": rel,
            "stats": stats,
            "sign_stats": sign_stats,
            "chamber_stats": chamber_stats,
            "top_p0": top_p0[0],
            "top_imbalance": top_imbalance[0],
            "top_pair_tax": top_pair_tax[0],
        }

    print("CONSECUTIVE RANK-ONE FAR TRIPLES")
    print(f"  far triples: ({m_start}, {m_start + 1}, {m_start + 2}) through ({m_stop}, {m_stop + 1}, {m_stop + 2})")
    print("  every far triple has exact relation rank 1 at height 3 (basis includes -u+2v-w=0)")
    for far, info in per_far.items():
        rel = info["relation"]
        stats = info["stats"]
        top_p0 = info["top_p0"]
        top_imb = info["top_imbalance"]
        top_tax = info["top_pair_tax"]
        print(
            f"  far={far} res7={residue_word(far)} best_relation={rel['best']} "
            f"rank={rel['exact_rank']} count={stats['primitive']} "
            f"actual(+/-)={stats['actual_+']}/{stats['actual_-']} "
            f"pairtax(+/-)={stats['pairtax_+']}/{stats['pairtax_-']} "
            f"same_sign={stats['same_actual_pairtax_True']}"
        )
        print(
            f"    top_p0 row={top_p0.row} p0={top_p0.p0} margin={top_p0.margin} "
            f"signs={''.join(top_p0.signs)} chamber={top_p0.chamber} "
            f"Nimb={top_p0.pkt.imbalance.norm()}"
        )
        print(
            f"    top_imb row={top_imb.row} p0={top_imb.p0} "
            f"signs={''.join(top_imb.signs)} chamber={top_imb.chamber} "
            f"Nimb={top_imb.pkt.imbalance.norm()}"
        )
        print(
            f"    top_pairtax row={top_tax.row} p0={top_tax.p0} "
            f"signs={''.join(top_tax.signs)} chamber={top_tax.chamber} "
            f"abs_pairtax={abs(top_tax.pkt.recursion)}"
        )

    return rows


def aggregate_phase_report(rows: list[AtlasRow]) -> None:
    print()
    print("AGGREGATE PHASE/SUPPORT ATLAS")
    print(f"  rows={len(rows)}")

    sign_stats: Counter[tuple[str, str, str, str]] = Counter(r.signs for r in rows)
    chamber_stats: Counter[str] = Counter(r.chamber for r in rows)
    phase_stats: Counter[tuple[int, ...]] = Counter(r.far_residue for r in rows)

    print("  most common sign ledgers (actual,pairtax,pair,triple):")
    for key, count in sign_stats.most_common(10):
        print(f"    {''.join(key)}: {count}")

    print("  A2/Eisenstein chamber counts:")
    for key, count in chamber_stats.most_common():
        print(f"    {key}: {count}")

    print("  mod-7 far phase counts:")
    for key, count in sorted(phase_stats.items()):
        print(f"    {key}: {count}")

    top_direct = sorted(rows, key=lambda r: (r.p0, abs(r.pkt.total), r.row), reverse=True)[:12]
    top_total = sorted(rows, key=lambda r: (abs(r.pkt.total), r.p0, r.row), reverse=True)[:12]
    top_pairtax = sorted(rows, key=lambda r: (abs(r.pkt.recursion), r.p0, r.row), reverse=True)[:12]
    top_imb = sorted(rows, key=lambda r: (r.pkt.imbalance.norm(), r.p0, r.row), reverse=True)[:12]

    def print_top(title: str, selected: list[AtlasRow], value) -> None:
        print(f"  top {title}:")
        for r in selected[:6]:
            print(
                f"    row={r.row} far={r.pkt.far} res7={r.far_residue} "
                f"{title}={value(r)} p0={r.p0} margin={r.margin} "
                f"signs={''.join(r.signs)} chamber={r.chamber} "
                f"core_hist7={r.core_hist} exc={r.row_exc} energy={r.energy}"
            )

    print_top("direct_p0", top_direct, lambda r: r.p0)
    print_top("|actual_residual|", top_total, lambda r: abs(r.pkt.total))
    print_top("|pair_tax_shadow|", top_pairtax, lambda r: abs(r.pkt.recursion))
    print_top("Eisenstein_norm", top_imb, lambda r: r.pkt.imbalance.norm())

    print("  top-12 overlap matrix:")
    groups = {
        "direct": {r.row for r in top_direct},
        "actual": {r.row for r in top_total},
        "pairtax": {r.row for r in top_pairtax},
        "eisenstein": {r.row for r in top_imb},
    }
    names = list(groups)
    for a in names:
        line = []
        for b in names:
            line.append(str(len(groups[a] & groups[b])))
        print(f"    {a}: {' '.join(line)}")

    signature_leaders(rows)
    tournament_analysis(top_direct, top_total, top_pairtax, top_imb)


def signature_leaders(rows: list[AtlasRow]) -> None:
    print("  signature leaders: (far_residue, signs, chamber) -> max direct p0")
    groups: dict[tuple[tuple[int, ...], tuple[str, str, str, str], str], list[AtlasRow]] = defaultdict(list)
    for row in rows:
        groups[(row.far_residue, row.signs, row.chamber)].append(row)

    leaders = []
    for key, bucket in groups.items():
        leader = max(bucket, key=lambda r: (r.p0, abs(r.pkt.total), r.row))
        leaders.append((leader, len(bucket), key))
    leaders.sort(key=lambda item: (item[0].p0, item[1], item[0].row), reverse=True)

    for leader, count, key in leaders[:12]:
        far_residue, signs, ch = key
        print(
            f"    res7={far_residue} signs={''.join(signs)} chamber={ch} "
            f"count={count} leader={leader.row} p0={leader.p0} margin={leader.margin} "
            f"total={leader.pkt.total} recur={leader.pkt.recursion} Nimb={leader.pkt.imbalance.norm()}"
        )

    print(
        "  reading: direct risk is not a function of relation rank alone; it is a finite "
        "address in (mod-7 phase, order signs, A2 chamber, bounded-core support)."
    )


def tournament_analysis(
    top_direct: list[AtlasRow],
    top_total: list[AtlasRow],
    top_pairtax: list[AtlasRow],
    top_imb: list[AtlasRow],
) -> None:
    print()
    print("TOURNAMENT ANALYSIS ON PHASE-ATLAS LENSES")
    lenses = {
        "direct_p0": top_direct[0],
        "actual_residual": top_total[0],
        "pair_tax_shadow": top_pairtax[0],
        "eisenstein_norm": top_imb[0],
    }
    wins = {name: 0 for name in lenses}

    def lens_key(row: AtlasRow) -> tuple[F, F, F, Row]:
        # Direct cap risk is the first proof obligation; ties keep the residual ledger.
        return (row.p0, abs(row.pkt.total), row.pkt.imbalance.norm(), row.row)

    lens_names = list(lenses)
    cycles = 0
    for a, b in combinations(lens_names, 2):
        wins[a if lens_key(lenses[a]) > lens_key(lenses[b]) else b] += 1
    for a, b, c in combinations(lens_names, 3):
        ab = lens_key(lenses[a]) > lens_key(lenses[b])
        bc = lens_key(lenses[b]) > lens_key(lenses[c])
        ca = lens_key(lenses[c]) > lens_key(lenses[a])
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            cycles += 1

    print(f"  leaders={{ {', '.join(f'{name}: {row.row}' for name, row in lenses.items())} }}")
    print(f"  score_hist={dict(Counter(wins.values()))} directed_3cycles={cycles}")
    print("  Hamiltonian_path=" + " > ".join(sorted(lens_names, key=lambda name: wins[name], reverse=True)))
    print(
        "  challenged assumption: the natural vertices are not far runners; they are "
        "phase/support proof lenses and packet signatures."
    )


def synthesis() -> None:
    print()
    print("SYNTHESIS")
    print("  Holding the rank-one relation u-2v+w=0 fixed does not freeze the packet signs.")
    print("  The mod-7 far phase and A2/Eisenstein chamber move large mass between actual residual,")
    print("  pair-tax shadow, and direct p0.  This explains why S51 saw AP-relation sign flips")
    print("  under shifts and why S52's cube-root coordinate is useful but not a proof scalar.")
    print("  Proof target: a finite resonant phase/support atlas for low-height relation tuples,")
    print("  followed by signed Abel/Koksma bounds only after the atlas key has been retained.")
    print("  No LRC(14) proof is claimed.")


def main() -> None:
    print("HYP-2682 / T921 -- LRC14 cube-root phase/support atlas\n")
    rows = scan_consecutive_triples()
    aggregate_phase_report(rows)
    synthesis()


if __name__ == "__main__":
    main()
