#!/usr/bin/env python3
"""HYP-2682 / T921: AP-triple phase atlas for the LRC(14) true-wide route.

The S51/S52 scouts showed that far triples F_m=(m,m+1,m+2) all carry the same
rank-one relation u-2v+w=0, but their signed Newton/cube-root phase changes
with m.  This script keeps the exact seven-packet data and asks whether the
rank-one branch looks like a finite small-offset atlas plus a direct-p0-safe
tail.

Incoming KPS S19 refutes the scalar C(k)=sup w|Delta_w| route.  Therefore this
scout reports direct p0 and cap margin first; residual packets are treated as
proof-obligation coordinates, not as a standalone danger scalar.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations

from lrc14_cube_root_order_filter_codex_s52 import TriplePackets, make_packets, sgn
from lrc14_wide_branch_ridge_codex_s47 import CAP, additive_energy, primitive, sumset_excess


Row = tuple[int, ...]


def fmt(q: F | None) -> str:
    if q is None:
        return "n/a"
    return f"{q} ({float(q):+.9f})"


def order_word(pkt: TriplePackets) -> str:
    """Signs of R1,R2,R3,total,pair-tax shadow."""

    return "".join(sgn(x) for x in (pkt.singles, pkt.pairs, pkt.g, pkt.total, pkt.recursion))


def packet_word(pkt: TriplePackets) -> str:
    return "".join(sgn(x) for x in (pkt.a, pkt.b, pkt.c, pkt.d, pkt.e, pkt.f, pkt.g))


def band(m: int) -> str:
    if m <= 30:
        return "15..30"
    if m <= 60:
        return "31..60"
    if m <= 90:
        return "61..90"
    return "91+"


@dataclass(frozen=True)
class PhaseRecord:
    label: str
    core: Row
    m: int
    pkt: TriplePackets

    @property
    def far(self) -> Row:
        return (self.m, self.m + 1, self.m + 2)

    @property
    def cap(self) -> F | None:
        return CAP.get(len(self.pkt.row))

    @property
    def margin(self) -> F | None:
        cap = self.cap
        return None if cap is None else cap - self.pkt.p0


def top_records(rows: list[PhaseRecord], key, keep: int = 5) -> list[PhaseRecord]:
    return sorted(rows, key=key, reverse=True)[:keep]


def print_record(prefix: str, rec: PhaseRecord) -> None:
    pkt = rec.pkt
    print(
        f"{prefix} label={rec.label} m={rec.m} mod7={rec.m % 7} "
        f"far={rec.far} row={pkt.row}"
    )
    print(
        f"    p0={fmt(pkt.p0)} margin={fmt(rec.margin)} "
        f"total={fmt(pkt.total)} recur={fmt(pkt.recursion)} "
        f"R1={fmt(pkt.singles)} R2={fmt(pkt.pairs)} R3={fmt(pkt.g)} "
        f"Nimb={pkt.imbalance.norm()}"
    )
    print(
        f"    order_word={order_word(pkt)} packet_word={packet_word(pkt)} "
        f"sumset_excess={sumset_excess(pkt.row)} energy={additive_energy(pkt.row)}"
    )


def tournament_from_values(names: list[str], values: dict[str, tuple[F, F, str]]) -> None:
    """Tournament where larger primary value wins; tie by secondary value/name."""

    n = len(names)
    adj = [[False] * n for _ in range(n)]
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i == j:
                continue
            adj[i][j] = values[a] > values[b]

    scores = {names[i]: sum(adj[i]) for i in range(n)}
    hist = Counter(scores.values())
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        edges = (adj[i][j], adj[j][k], adj[k][i])
        redges = (adj[j][i], adj[k][j], adj[i][k])
        if all(edges) or all(redges):
            c3 += 1

    # Hamiltonian-path count by DP; n is intentionally tiny here.
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            cur = dp.get((mask, last), 0)
            if not cur:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + cur
    hp_count = sum(dp.get(((1 << n) - 1, i), 0) for i in range(n))

    path = sorted(names, key=lambda name: values[name], reverse=True)
    print("  TOURNAMENT ANALYSIS")
    print("    vertices are proof lenses/core families, not runners.")
    print("    pairwise observable=(primary exact value, secondary exact value, stable label)")
    print("    gauge: larger observable wins; ties are lexicographic through the stored label")
    print(f"    scores={scores}")
    print(f"    score_hist={dict(hist)} directed_3cycles={c3} Hamiltonian_path_count={hp_count}")
    print(f"    tie Hamiltonian path={' > '.join(path)}")


def named_core_records(max_m: int = 120) -> list[PhaseRecord]:
    print("NAMED CORE AP-TRIPLE PHASE SCAN")
    cores: list[tuple[str, Row, int]] = [
        ("dilated", (0, 4, 6, 8, 10, 12, 14), 15),
        ("consec8", (0, 1, 2, 3, 4, 5, 6, 7), 15),
        ("direct_p0_leader_core", (0, 9, 10, 11, 12, 13, 14), 15),
        ("s51_top_dev_core", (0, 5, 6, 7, 11, 13, 14), 15),
        ("third_pocket_mixed_core", (0, 3, 5, 16, 28), 30),
    ]

    all_records: list[PhaseRecord] = []
    family_values: dict[str, tuple[F, F, str]] = {}

    for label, core, start in cores:
        rows: list[PhaseRecord] = []
        for m in range(start, max_m + 1):
            far = (m, m + 1, m + 2)
            if set(core) & set(far):
                continue
            pkt = make_packets(core, far)
            rec = PhaseRecord(label, core, m, pkt)
            rows.append(rec)
            all_records.append(rec)

        print(f"\nCORE {label}: core={core}, m={start}..{max_m}, rows={len(rows)}")
        print(f"  order_word_counts={dict(Counter(order_word(r.pkt) for r in rows))}")
        print(f"  packet_word_top={Counter(packet_word(r.pkt) for r in rows).most_common(8)}")

        by_mod7: dict[int, Counter[str]] = defaultdict(Counter)
        by_band: dict[str, list[PhaseRecord]] = defaultdict(list)
        for rec in rows:
            by_mod7[rec.m % 7][order_word(rec.pkt)] += 1
            by_band[band(rec.m)].append(rec)
        print("  order words by m mod 7:")
        for r in range(7):
            print(f"    {r}: {dict(by_mod7[r])}")
        print("  band maxima:")
        for name in ("15..30", "31..60", "61..90", "91+"):
            chunk = by_band.get(name, [])
            if not chunk:
                continue
            best_p0 = max(chunk, key=lambda r: r.pkt.p0)
            best_total = max(chunk, key=lambda r: abs(r.pkt.total))
            print(
                f"    {name}: max_p0 m={best_p0.m} p0={best_p0.pkt.p0} "
                f"margin={best_p0.margin}; max_abs_total m={best_total.m} "
                f"|total|={abs(best_total.pkt.total)}"
            )

        print("  top direct p0:")
        for rec in top_records(rows, lambda r: (r.pkt.p0, abs(r.pkt.total), -r.m)):
            print_record("   ", rec)
        print("  top |actual residual|:")
        for rec in top_records(rows, lambda r: (abs(r.pkt.total), r.pkt.p0, -r.m)):
            print_record("   ", rec)
        print("  top Eisenstein imbalance norm:")
        for rec in top_records(rows, lambda r: (r.pkt.imbalance.norm(), r.pkt.p0, -r.m)):
            print_record("   ", rec)

        best = max(rows, key=lambda r: (r.pkt.p0, abs(r.pkt.total), -r.m))
        family_values[label] = (best.pkt.p0, abs(best.pkt.total), label)

    print("\nNAMED CORE FAMILY TOURNAMENT")
    tournament_from_values([name for name, _core, _start in cores], family_values)
    return all_records


def bank_for_m(m: int) -> list[PhaseRecord]:
    rows: list[PhaseRecord] = []
    far = (m, m + 1, m + 2)
    for rest in combinations(range(1, 15), 6):
        core = (0,) + rest
        row = tuple(sorted(core + far))
        if not primitive(row):
            continue
        rows.append(PhaseRecord(f"bank_m{m}", core, m, make_packets(core, far)))
    return rows


def selected_all_core_banks(ms: tuple[int, ...] = (15, 16, 22, 28, 42)) -> None:
    print("\nSELECTED ALL-CORE BANKS")
    print("  cores=(0)+6-subsets of [1,14]; far=(m,m+1,m+2)")
    for m in ms:
        rows = bank_for_m(m)
        print(f"\nBANK m={m}, far={(m, m + 1, m + 2)}, primitive rows={len(rows)}")
        print(f"  total_signs={dict(Counter(sgn(r.pkt.total) for r in rows))}")
        print(f"  recursion_signs={dict(Counter(sgn(r.pkt.recursion) for r in rows))}")
        print(f"  order_word_top={Counter(order_word(r.pkt) for r in rows).most_common(10)}")
        print(f"  packet_word_top={Counter(packet_word(r.pkt) for r in rows).most_common(10)}")
        print(f"  same_total_recursion_sign={sum(sgn(r.pkt.total) == sgn(r.pkt.recursion) for r in rows)}")

        leaders = {
            "direct_p0": max(rows, key=lambda r: (r.pkt.p0, abs(r.pkt.total), r.core)),
            "actual_residual": max(rows, key=lambda r: (abs(r.pkt.total), r.pkt.p0, r.core)),
            "pair_tax_shadow": max(rows, key=lambda r: (abs(r.pkt.recursion), r.pkt.p0, r.core)),
            "triple_packet": max(rows, key=lambda r: (abs(r.pkt.g), r.pkt.p0, r.core)),
            "eisenstein_imbalance": max(rows, key=lambda r: (r.pkt.imbalance.norm(), r.pkt.p0, r.core)),
        }
        for name, rec in leaders.items():
            print_record(f"  leader {name}", rec)

        tournament_from_values(
            list(leaders),
            {
                name: (rec.pkt.p0, abs(rec.pkt.total), name)
                for name, rec in leaders.items()
            },
        )


def synthesis() -> None:
    print("\nSYNTHESIS")
    print("  Incoming KPS S19 refutes the scalar C(k)=sup w|Delta_w| route.")
    print("  Therefore AP-triple packets are not being used to bound C(k).")
    print("  They are phase coordinates for the direct wide=>small-p0 coverage lemma.")
    print("  Rank-one relation u-2v+w=0 is preserved by every AP triple, but packet")
    print("  signs, pair-tax shadow, and Eisenstein imbalance vary with offset and core.")
    print("  A finite resonant atlas must retain this phase/support address before")
    print("  applying any signed relation-lattice or Freiman/GAP reduction.")
    print("  No LRC(14) proof is claimed.")


def main() -> None:
    print("HYP-2682 / T921 -- AP-triple phase atlas for LRC14")
    print()
    named_core_records()
    selected_all_core_banks()
    synthesis()


if __name__ == "__main__":
    main()
