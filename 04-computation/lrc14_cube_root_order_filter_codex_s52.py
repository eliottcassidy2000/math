#!/usr/bin/env python3
"""HYP-2681 / T920: cube-root order filter for LRC(14) multi-far packets.

For a far triple (u,v,w), write the seven nonempty simultaneous-peel residuals

    A,B,C = one-far packets
    D,E,F = two-far packets
    G     = three-far packet.

The actual correction is A+B+C+D+E+F+G.  The user's recursion
A+B+C-D-E-F+G is the pair-tax shadow.  This scout keeps the exact C3 cyclic
modes

    A + omega B + omega^2 C,   D + omega E + omega^2 F

as Eisenstein pairs a+b*omega with omega^2+omega+1=0.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import sqrt

from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import (
    delta_s,
    phi_s,
    p0_cached,
)
from lrc14_wide_branch_ridge_codex_s47 import CAP, primitive, sumset_excess, additive_energy


Row = tuple[int, ...]


def fmt(q: F) -> str:
    return f"{q} ({float(q):+.9f})"


def sgn(q: F) -> str:
    if q > 0:
        return "+"
    if q < 0:
        return "-"
    return "0"


@dataclass(frozen=True)
class Eisenstein:
    """Exact a + b*omega, omega^2 + omega + 1 = 0."""

    a: F
    b: F

    def __add__(self, other: "Eisenstein") -> "Eisenstein":
        return Eisenstein(self.a + other.a, self.b + other.b)

    def __sub__(self, other: "Eisenstein") -> "Eisenstein":
        return Eisenstein(self.a - other.a, self.b - other.b)

    def norm(self) -> F:
        return self.a * self.a - self.a * self.b + self.b * self.b

    def approx_abs(self) -> float:
        return sqrt(float(self.norm()))

    def __str__(self) -> str:
        return f"({self.a}) + ({self.b})*omega ; N={self.norm()} ({self.approx_abs():.9f})"


def omega_mode(x0: F, x1: F, x2: F) -> Eisenstein:
    """Return x0 + omega*x1 + omega^2*x2 in basis 1, omega."""

    return Eisenstein(x0 - x2, x1 - x2)


def residual(core: Row, far: Row) -> F:
    return delta_s(core, far) - phi_s(core, len(far))


@dataclass(frozen=True)
class TriplePackets:
    core: Row
    far: Row
    a: F
    b: F
    c: F
    d: F
    e: F
    f: F
    g: F

    @property
    def row(self) -> Row:
        return tuple(sorted(self.core + self.far))

    @property
    def singles(self) -> F:
        return self.a + self.b + self.c

    @property
    def pairs(self) -> F:
        return self.d + self.e + self.f

    @property
    def total(self) -> F:
        return self.singles + self.pairs + self.g

    @property
    def recursion(self) -> F:
        return self.singles - self.pairs + self.g

    @property
    def pair_tax_checksum(self) -> F:
        return self.total - 2 * self.pairs

    @property
    def single_mode(self) -> Eisenstein:
        return omega_mode(self.a, self.b, self.c)

    @property
    def pair_mode(self) -> Eisenstein:
        # Cyclic edge order: uv -> vw -> wu.
        return omega_mode(self.d, self.e, self.f)

    @property
    def imbalance(self) -> Eisenstein:
        return self.single_mode - self.pair_mode

    @property
    def p0(self) -> F:
        return p0_cached(self.row)

    @property
    def margin(self) -> F | None:
        cap = CAP.get(len(self.row))
        return None if cap is None else cap - self.p0


def make_packets(core: Row, far: Row) -> TriplePackets:
    core = tuple(sorted(core))
    u, v, w = tuple(sorted(far))
    return TriplePackets(
        core=core,
        far=(u, v, w),
        a=residual(core, (u,)),
        b=residual(core, (v,)),
        c=residual(core, (w,)),
        d=residual(core, (u, v)),
        e=residual(core, (v, w)),
        f=residual(core, (w, u)),
        g=residual(core, (u, v, w)),
    )


def print_packets(label: str, pkt: TriplePackets) -> None:
    print(f"CASE {label}")
    print(f"  core={pkt.core} far={pkt.far} row={pkt.row}")
    print(f"  p0={fmt(pkt.p0)} margin={fmt(pkt.margin) if pkt.margin is not None else 'n/a'}")
    print(f"  row_exc={sumset_excess(pkt.row)} energy={additive_energy(pkt.row)} primitive={primitive(pkt.row)}")
    print("  packets:")
    print(f"    A={fmt(pkt.a)}  B={fmt(pkt.b)}  C={fmt(pkt.c)}")
    print(f"    D={fmt(pkt.d)}  E={fmt(pkt.e)}  F={fmt(pkt.f)}")
    print(f"    G={fmt(pkt.g)}")
    print(f"  actual H(1)=A+B+C+D+E+F+G = {fmt(pkt.total)} sign={sgn(pkt.total)}")
    print(f"  pair_tax_shadow A+B+C-D-E-F+G = {fmt(pkt.recursion)} sign={sgn(pkt.recursion)}")
    print(f"  checksum H(1)-2(D+E+F) equal recursion: {pkt.pair_tax_checksum == pkt.recursion}")
    print(f"  singles_sum={fmt(pkt.singles)} pairs_sum={fmt(pkt.pairs)} triple={fmt(pkt.g)}")
    print(f"  S_omega={pkt.single_mode}")
    print(f"  P_omega={pkt.pair_mode}")
    print(f"  S_omega-P_omega={pkt.imbalance}")
    print()


def named_reports() -> list[TriplePackets]:
    named = [
        ("dilated consecutive triple", (0, 4, 6, 8, 10, 12, 14), (15, 16, 17)),
        ("dilated shifted AP triple", (0, 4, 6, 8, 10, 12, 14), (22, 23, 24)),
        ("dilated positive phase triple", (0, 4, 6, 8, 10, 12, 14), (43, 44, 45)),
        ("dilated separated triple", (0, 4, 6, 8, 10, 12, 14), (17, 23, 31)),
        ("consec8 consecutive triple", (0, 1, 2, 3, 4, 5, 6, 7), (15, 16, 17)),
        ("third-pocket active triple", (0, 3, 5, 16, 28), (30, 33, 35)),
    ]
    reports: list[TriplePackets] = []
    for label, core, far in named:
        pkt = make_packets(core, far)
        reports.append(pkt)
        print_packets(label, pkt)
    return reports


def push_top(store: list[TriplePackets], pkt: TriplePackets, key, keep: int = 8) -> None:
    store.append(pkt)
    store.sort(key=key, reverse=True)
    del store[keep:]


def small_core_bank() -> None:
    print("SMALL CORE BANK: cube-root packets for far=(15,16,17)")
    far = (15, 16, 17)
    stats: Counter[str] = Counter()
    top_total: list[TriplePackets] = []
    top_recursion: list[TriplePackets] = []
    top_imbalance: list[TriplePackets] = []
    top_p0: list[TriplePackets] = []

    for rest in combinations(range(1, 15), 6):
        core = (0,) + rest
        row = tuple(sorted(core + far))
        if not primitive(row):
            continue
        pkt = make_packets(core, far)
        stats[f"total_{sgn(pkt.total)}"] += 1
        stats[f"recursion_{sgn(pkt.recursion)}"] += 1
        stats[f"G_{sgn(pkt.g)}"] += 1
        stats[f"total_recursion_same_{pkt.total == 0 or pkt.recursion == 0 or sgn(pkt.total) == sgn(pkt.recursion)}"] += 1
        stats[f"pair_tax_{sgn(pkt.pairs)}"] += 1
        push_top(top_total, pkt, lambda p: (abs(p.total), p.p0, p.row))
        push_top(top_recursion, pkt, lambda p: (abs(p.recursion), p.p0, p.row))
        push_top(top_imbalance, pkt, lambda p: (p.imbalance.norm(), p.p0, p.row))
        push_top(top_p0, pkt, lambda p: (p.p0, abs(p.total), p.row))

    print(f"  stats={dict(stats)}")

    def print_top(title: str, rows: list[TriplePackets], value) -> None:
        print(f"  top {title}:")
        for pkt in rows[:6]:
            print(
                f"    row={pkt.row} {title}={value(pkt)} p0={pkt.p0} margin={pkt.margin} "
                f"total={pkt.total} recur={pkt.recursion} pair={pkt.pairs} "
                f"G={pkt.g} Nimb={pkt.imbalance.norm()}"
            )

    print_top("|actual residual|", top_total, lambda p: abs(p.total))
    print_top("|A+B+C-D-E-F+G|", top_recursion, lambda p: abs(p.recursion))
    print_top("Eisenstein imbalance norm", top_imbalance, lambda p: p.imbalance.norm())
    print_top("direct p0", top_p0, lambda p: p.p0)
    tournament_analysis(top_total, top_recursion, top_imbalance, top_p0)
    print()


def tournament_analysis(
    top_total: list[TriplePackets],
    top_recursion: list[TriplePackets],
    top_imbalance: list[TriplePackets],
    top_p0: list[TriplePackets],
) -> None:
    print("  TOURNAMENT ANALYSIS ON CUBE-ROOT PROOF LENSES")
    lenses = ["actual_residual", "pair_tax_shadow", "eisenstein_imbalance", "direct_p0"]
    leaders = {
        "actual_residual": top_total[0].row,
        "pair_tax_shadow": top_recursion[0].row,
        "eisenstein_imbalance": top_imbalance[0].row,
        "direct_p0": top_p0[0].row,
    }
    wins = {lens: 0 for lens in lenses}
    # Lens i beats lens j if its leader has larger direct p0; tie by larger actual residual.
    leader_pkt = {
        "actual_residual": top_total[0],
        "pair_tax_shadow": top_recursion[0],
        "eisenstein_imbalance": top_imbalance[0],
        "direct_p0": top_p0[0],
    }
    for a, b in combinations(lenses, 2):
        ka = (leader_pkt[a].p0, abs(leader_pkt[a].total), leaders[a])
        kb = (leader_pkt[b].p0, abs(leader_pkt[b].total), leaders[b])
        wins[a if ka > kb else b] += 1
    print(f"    leaders={leaders}")
    print(f"    score_hist={dict(Counter(wins.values()))}")
    print("    Hamiltonian_path=" + " > ".join(sorted(lenses, key=lambda x: wins[x], reverse=True)))
    print("    challenged assumption: cube-root imbalance is not automatically direct risk.")


def synthesis() -> None:
    print("SYNTHESIS")
    print("  A+B+C-D-E-F+G is exact, but it is the pair-tax shadow H(1)-2*R2, not the cap residual.")
    print("  The cube-root modes split cyclic imbalance among singles and pairs without numeric approximation.")
    print("  The useful proof lens is a three-coordinate packet: actual residual, pair-tax shadow,")
    print("  and Eisenstein imbalance.  These can rank different rows, so none should be collapsed early.")


def main() -> None:
    print("HYP-2681 / T920 -- cube-root order filter for LRC14 multi-far packets\n")
    named_reports()
    small_core_bank()
    synthesis()


if __name__ == "__main__":
    main()
