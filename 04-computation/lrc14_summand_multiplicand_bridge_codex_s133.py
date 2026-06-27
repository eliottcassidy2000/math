#!/usr/bin/env python3
"""S133: summand and multiplicand graph bridge for the LRC14 Farey split.

This is the follow-up to S131.  S131 read the Farey product payload

    M(S)=p/q  ->  K_{p,q},     p*q = |E(K_{p,q})|.

Here we reconnect that payload to the older summand/product graph work.
There are three nearby, easily conflated graph layers:

1. Farey pair ledger: (p,q) is an additive pair at node p+q.
2. LRC shell ledger: for n=14, q=27 is the odd summand node C=2n-1.
3. Multiplicand ledger: (p,q) is a factor pair at product node p*q.

The useful new distinction is that the LRC14 row 2/27 belongs to the
C=27 summand-unit quotient (HYP-2083/HYP-2161), while 3/41 is the first
unit-excess row whose K_{p,q} expansion has a K_{3,3} incidence wall.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from math import gcd, isqrt
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
N = 14
C_SECOND = 2 * N - 1


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s130 = load_module(
    "s130_mutated_farey_for_s133",
    REPO / "04-computation" / "lrc14_mutated_farey_tournament_codex_s130.py",
)
s124 = s130.s124
s127 = s130.s127


@dataclass(frozen=True)
class BipartiteLedger:
    p: int
    q: int

    @property
    def edges(self) -> int:
        return self.p * self.q

    @property
    def rank(self) -> str:
        if min(self.p, self.q) >= 3:
            return "K33-wall"
        if min(self.p, self.q) == 2:
            return "two-block"
        return "star"


def summand_pairs(m: int) -> list[tuple[int, int]]:
    return [(a, m - a) for a in range(1, (m + 1) // 2) if a < m - a]


def factor_pairs(m: int) -> list[tuple[int, int]]:
    out = []
    for a in range(2, isqrt(m) + 1):
        if m % a == 0:
            b = m // a
            if a < b:
                out.append((a, b))
    return out


def factorization(m: int) -> str:
    n = m
    parts = []
    p = 2
    while p * p <= n:
        if n % p:
            p += 1
            continue
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        parts.append(f"{p}^{e}" if e > 1 else str(p))
        p += 1
    if n > 1:
        parts.append(str(n))
    return " * ".join(parts) if parts else "1"


def fmt_pairs(pairs: list[tuple[int, int]], limit: int = 4) -> str:
    if not pairs:
        return "-"
    shown = ", ".join(f"{a}x{b}" for a, b in pairs[:limit])
    if len(pairs) > limit:
        shown += f", ...(+{len(pairs) - limit})"
    return shown


def antipodal_shell_strata(C: int) -> dict[int, list[tuple[int, int]]]:
    strata: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for a in range(1, (C + 1) // 2):
        strata[gcd(a, C)].append((a, C - a))
    return dict(sorted(strata.items()))


def line(title: str) -> None:
    print()
    print(title)
    print("-" * len(title))


def operation_shadow_contrast() -> None:
    line("[1] Operation-graph contrast from the older repo thread")
    print("  summand graph:")
    print("    labelled fiber at node m keeps pairs a+b=m.")
    print("    after forgetting labels, addition collapses to x<m, the natural order.")
    print("  multiplicand/product graph:")
    print("    labelled fiber at node m keeps pairs a*b=m.")
    print("    after forgetting labels, multiplication still leaves the divisor DAG.")
    print("  S133 rule:")
    print("    for M=p/q, keep BOTH the additive pair (p,q) at p+q and the")
    print("    multiplicand pair (p,q) at p*q; then expand the latter to K_{p,q}.")
    print("  guardrail:")
    print("    the LRC C=2n-1 summand graph is the denominator shell q=27, not p+q.")


def unit_excess_chain(limit: int = 8) -> None:
    line("[2] LRC14 unit-excess chain with summand and multiplicand shadows")
    print("  e=14p-q=1, q=14p-1.  The q column remains the theorem scale.")
    header = (
        f"  {'p/q':>7s} {'p+q':>5s} {'sum pairs':>9s} {'p*q':>5s} "
        f"{'factor fibers':>24s} {'defect':>7s} {'K-rank':>10s}  reading"
    )
    print(header)
    for p in range(1, limit + 1):
        q = N * p - 1
        M = F(p, q)
        led = BipartiteLedger(p, q)
        add_node = p + q
        prod_node = p * q
        defect = prod_node - add_node
        factors = factor_pairs(prod_node)
        if p == 1:
            reading = "q-threshold parent; additive star"
        elif p == 2:
            reading = "q=C=27 second-gap shell; planar two-block"
        elif p == 3:
            reading = "first K33 wall; near-miss 12->36"
        else:
            reading = "higher K33 wall"
        print(
            f"  {str(M):>7s} {add_node:5d} {len(summand_pairs(add_node)):9d} "
            f"{prod_node:5d} {fmt_pairs(factors):>24s} {defect:7d} "
            f"{led.rank:>10s}  {reading}"
        )
    print()
    print("  product-sum defect = p*q - (p+q).  For p>=2 this is the number")
    print("  of 1s needed to turn the two-element core {p,q} into a")
    print("  product-sum equality.  The p=1 parent is the identity-degenerate")
    print("  case, so its defect is negative rather than a core extension count.")


def denominator_shell_focus() -> None:
    line("[3] The denominator shell: why 2/27 is different from 3/41")
    print(f"  For n=14, the second-gap/shell modulus is C=2n-1={C_SECOND}.")
    strata = antipodal_shell_strata(C_SECOND)
    print("  C=27 antipodal summand pairs split by gcd:")
    for g, pairs in strata.items():
        tag = "unit-visible" if g == 1 else f"nonunit gcd={g}"
        print(f"    {tag:16s}: count={len(pairs):2d}  {fmt_pairs(pairs, limit=5).replace('x', '+')}")
    print()
    print("  Chain interpretation:")
    print("    1/13: right Farey parent; not the C=27 shell.")
    print("    2/27: exactly the S553 value 2/(2n-1).  This is the summand")
    print("          graph at C=27 acted on by units; composite 27 creates")
    print("          nonunit holes (gcd 3 and 9 strata).")
    print("    3/41: leaves the C=27 shell and enters the first K_{3,3}")
    print("          product-incidence wall.  It should be handled by a")
    print("          three-owner/K33 packet, not by the second-gap unit clock.")


def selected_lrc_rows() -> None:
    line("[4] Named LRC14 rows under the bridge")
    ap = tuple(range(1, 14))
    rows = [
        s130.Row("AP", ap, "tight"),
        s130.Row("GW 12->24", tuple(list(range(1, 12)) + [13, 24]), "tight"),
        s130.Row("residue-liar 12->26", tuple(list(range(1, 12)) + [13, 26]), "loose"),
        s130.Row("petal 10->20", tuple(sorted((set(ap) - {10}) | {20})), "loose"),
        s130.Row("petal 13->26", tuple(sorted((set(ap) - {13}) | {26})), "loose"),
        s130.Row("near-miss 12->36", tuple(list(range(1, 12)) + [13, 36]), "loose"),
    ]
    print(
        f"  {'row':24s} {'M':>7s} {'q(S)':>4s} {'p+q':>5s} {'p*q':>5s} "
        f"{'K-rank':>10s}  bridge tag"
    )
    for row in rows:
        M = s130.exact_M(row.speeds)
        p, q = M.numerator, M.denominator
        led = BipartiteLedger(p, q)
        if M == F(2, C_SECOND):
            tag = "C27 second-gap / petal two-block"
        elif M == F(3, 41):
            tag = "Farey-41 K33 wall"
        elif M == F(1, N):
            tag = "apex floor"
        elif p == 1:
            tag = "star/q-threshold"
        else:
            tag = "off-shell"
        print(
            f"  {row.label:24s} {str(M):>7s} {s124.q_threshold(row.speeds):4d} "
            f"{p + q:5d} {p * q:5d} {led.rank:>10s}  {tag}"
        )


def product_fiber_closeup() -> None:
    line("[5] Multiplicand graph closeup")
    for p, q in [(2, 27), (3, 41), (4, 55)]:
        prod = p * q
        led = BipartiteLedger(p, q)
        print(f"  {p}/{q}: product node {prod} = {factorization(prod)}")
        print(f"    product-graph factor fibers: {fmt_pairs(factor_pairs(prod), limit=8)}")
        print(f"    K_{{{p},{q}}}: edges={led.edges}, rank={led.rank}")
    print()
    print("  Reading: the multiplicand graph and K_{p,q} are not the same object.")
    print("  The product graph records factor fibers of the integer p*q; K_{p,q}")
    print("  is the incidence blow-up of the particular fiber (p,q).  The p=2")
    print("  row has rich 3-adic product fibers but no K33 blow-up; p=3 has a")
    print("  simple product node but enough incidence rank to cross K33.")


@dataclass(frozen=True)
class Carrier:
    name: str
    theorem_safety: int
    additive_shell: int
    multiplicand_fiber: int
    k33_signal: int
    lrc_specificity: int
    false_positive_resistance: int

    def score(self) -> tuple[int, int, int, int, int, int]:
        return (
            self.theorem_safety,
            self.additive_shell,
            self.multiplicand_fiber,
            self.k33_signal,
            self.lrc_specificity,
            self.false_positive_resistance,
        )


def carrier_tournament() -> None:
    line("[6] Tournament Analysis on bridge carriers")
    carriers = [
        Carrier("q-binding", 5, 3, 1, 1, 5, 5),
        Carrier("C27-shell", 4, 5, 1, 1, 5, 4),
        Carrier("Farey p+q", 3, 4, 1, 0, 3, 3),
        Carrier("product fiber", 2, 1, 5, 2, 2, 2),
        Carrier("Kpq blowup", 2, 1, 4, 3, 3, 2),
        Carrier("K33 packet", 1, 0, 3, 5, 4, 3),
        Carrier("raw apex iso", 0, 1, 0, 0, 1, 0),
    ]
    wins: dict[int, set[int]] = {i: set() for i in range(len(carriers))}
    bitmask = 0
    bit = 0
    flips_vs_product_only = 0
    for i, j in combinations(range(len(carriers)), 2):
        ci, cj = carriers[i], carriers[j]
        if ci.score() >= cj.score():
            winner, loser = i, j
        else:
            winner, loser = j, i
        wins[winner].add(loser)
        if winner == i:
            bitmask |= 1 << bit
        product_pref = i if ci.multiplicand_fiber >= cj.multiplicand_fiber else j
        if product_pref != winner:
            flips_vs_product_only += 1
        bit += 1

    fp = s127.tournament_fingerprint(bitmask, len(carriers))
    ham = 0
    first_path: tuple[str, ...] | None = None
    for perm in permutations(range(len(carriers))):
        if all(perm[t + 1] in wins[perm[t]] for t in range(len(perm) - 1)):
            ham += 1
            if first_path is None:
                first_path = tuple(carriers[i].name for i in perm)
    score_hist = Counter(len(v) for v in wins.values())
    print("  vertices considered/challenged:")
    print("    runners, gaps, denominator shells, Farey pairs, factor fibers,")
    print("    Kpq incidence blowups, K33 packets, and proof obligations.")
    print("  chosen vertices: proof carriers, not runners.")
    print("  pair observable: role score")
    print("    (theorem safety, additive shell retention, multiplicand fiber retention,")
    print("     K33 signal, LRC specificity, false-positive resistance)")
    print("  switch/gauge: lexicographically larger role score wins.")
    print("  preserves: which carrier should handle which branch of the LRC14 split.")
    print("  destroys: exact safe arcs, runner labels, and continuous phase geometry.")
    print(
        f"  fingerprint: score_hist={dict(sorted(score_hist.items()))} "
        f"c3={fp['c3']} scc={fp['scc']} hp={ham}"
    )
    print(f"  first Hamiltonian path: {' > '.join(first_path or ())}")
    print(f"  product-only gauge would flip {flips_vs_product_only} carrier pairs.")


def proof_readout() -> None:
    line("[7] Proof readout")
    print("  New split refined from S131:")
    print("    p=2 is not just planar.  It is the exact C=27 second-gap")
    print("    summand-unit quotient, so the right tools are HYP-2083/HYP-2161,")
    print("    petal rigidity, nonunit gcd strata, and lift/CRT conservativity.")
    print("    p>=3 is a different branch.  It starts at 3/41, where the")
    print("    denominator leaves C=27 and the K_{p,q} incidence blow-up has")
    print("    a K_{3,3} minor.  That branch should feed a finite three-owner")
    print("    obstruction packet, possibly the HYP-2908 tournament-state lift.")
    print("  Salient theorem target:")
    print("    prove every remaining q=14 non-AP/GW atom either reduces to the")
    print("    C=27 p=2 shell/petal branch or crosses the p>=3 K33 packet wall.")
    print("    Do not scalarize product fibers or apex iso classes before the")
    print("    q-shell, owner, and incidence labels have done their work.")


def main() -> None:
    print("S133 LRC14 SUMMAND x MULTIPLICAND FAREY BRIDGE")
    print("=" * 78)
    operation_shadow_contrast()
    unit_excess_chain()
    denominator_shell_focus()
    selected_lrc_rows()
    product_fiber_closeup()
    carrier_tournament()
    proof_readout()


if __name__ == "__main__":
    main()
