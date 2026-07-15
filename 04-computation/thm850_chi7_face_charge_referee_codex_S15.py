#!/usr/bin/env python3
"""Exact referee for THM-850, the chi-seven B3 face-charge law.

The arithmetic is integral except for the quadratic cyclotomic field.  That
field is represented exactly as pairs a+b*eta with eta^2+eta+2=0; no floating
point or external computer algebra is used.

Tournament Analysis uses the seven nonempty B3 strata as vertices.  The
pairwise observable is charge difference in F_7, switched to an orientation
by the Legendre character.  The increasing residue order is the tie
Hamiltonian path (there are no actual ties).

Assumption challenge: vertices could have been runners, chords, gap arcs,
wall events, residues, or proof obligations.  Face strata are used here
because this theorem preserves the B3 Mobius sign and operation congruence.
It does not preserve a tiling, merged node, blue/black colour, or LRC data.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from math import comb


QR = frozenset({1, 2, 4})
NQR = frozenset({3, 5, 6})
CHARGE = {"A": 1, "B": 2, "C": 4}


def chi7(q: int) -> int:
    q %= 7
    if q == 0:
        return 0
    return 1 if q in QR else -1


def augmented_chi7(q: int) -> int:
    q %= 7
    return 1 if q == 0 else chi7(q)


@dataclass(frozen=True)
class Quad:
    """a+b*eta in Z[eta], where eta^2+eta+2=0."""

    a: int
    b: int = 0

    def __add__(self, other: "Quad") -> "Quad":
        return Quad(self.a + other.a, self.b + other.b)

    def __sub__(self, other: "Quad") -> "Quad":
        return Quad(self.a - other.a, self.b - other.b)

    def __neg__(self) -> "Quad":
        return Quad(-self.a, -self.b)

    def __mul__(self, other: "Quad") -> "Quad":
        # eta^2=-eta-2
        return Quad(
            self.a * other.a - 2 * self.b * other.b,
            self.a * other.b + self.b * other.a - self.b * other.b,
        )

    def __str__(self) -> str:
        names = {
            (1, 0): "1",
            (0, 1): "eta",
            (-1, 0): "-1",
            (-1, -1): "eta_bar",
            (0, 0): "0",
        }
        return names.get((self.a, self.b), f"{self.a:+d}{self.b:+d}*eta")


ZERO = Quad(0)
ONE = Quad(1)
ETA = Quad(0, 1)
ETA_BAR = Quad(-1, -1)


def charge_census(r: int) -> Counter[int]:
    """Counts x+4y+2z mod 7 on x+y+z=r (A,C,B address order)."""
    out: Counter[int] = Counter()
    for x in range(r + 1):
        for y in range(r - x + 1):
            z = r - x - y
            out[(x + 4 * y + 2 * z) % 7] += 1
    return out


def predicted_census(r: int) -> dict[int, int]:
    m = comb(r + 2, 2)
    s = r % 7
    if s in (0, 4):
        return {q: (m + 6) // 7 if q == 0 else (m - 1) // 7 for q in range(7)}
    if s == 1:
        return {
            q: (m + 4) // 7 if q in QR else (m - 3) // 7
            for q in range(7)
        }
    if s == 2:
        return {q: (m - 6) // 7 if q == 0 else (m + 1) // 7 for q in range(7)}
    if s == 3:
        return {
            q: (m + 4) // 7 if q in NQR else (m - 3) // 7
            for q in range(7)
        }
    return {q: m // 7 for q in range(7)}


def weighted_census(census: Counter[int]) -> Quad:
    # sum_q c(q) zeta^q.  The formula first checks that the values are
    # constant on QR and NQR, as the theorem predicts.
    qr_values = {census[q] for q in QR}
    nqr_values = {census[q] for q in NQR}
    assert len(qr_values) == len(nqr_values) == 1
    c_qr = next(iter(qr_values))
    c_nqr = next(iter(nqr_values))
    return Quad(census[0]) + Quad(c_qr) * ETA + Quad(c_nqr) * ETA_BAR


def cyclotomic_sequence(limit: int) -> list[Quad]:
    out = [ONE, ETA, -ONE]
    while len(out) <= limit:
        out.append(ETA * out[-1] - ETA_BAR * out[-2] + out[-3])
    return out[: limit + 1]


def directed_3cycles(adj: list[list[int]]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        total += int(
            (adj[a][b] and adj[b][c] and adj[c][a])
            or (adj[a][c] and adj[c][b] and adj[b][a])
        )
    return total


def hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            if not dp[mask][v]:
                continue
            for w in range(n):
                if not (mask & (1 << w)) and adj[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[-1])


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach(start: int, reverse: bool) -> set[int]:
        seen = {start}
        todo = [start]
        while todo:
            v = todo.pop()
            for w in range(n):
                edge = adj[w][v] if reverse else adj[v][w]
                if edge and w not in seen:
                    seen.add(w)
                    todo.append(w)
        return seen

    unseen = set(range(n))
    sizes: list[int] = []
    while unseen:
        v = min(unseen)
        component = reach(v, False) & reach(v, True)
        sizes.append(len(component))
        unseen -= component
    return sorted(sizes, reverse=True)


def main() -> None:
    print("THM-850 exact chi-seven face-charge referee")
    print("=" * 72)

    print("\nA. Exact B3 sign convention")
    sign_word = []
    for mask in range(1, 8):
        subset = tuple(role for bit, role in enumerate(("A", "B", "C")) if mask & (1 << bit))
        q = sum(CHARGE[role] for role in subset) % 7
        mobius = 1 if len(subset) % 2 else -1
        assert augmented_chi7(q) == mobius
        sign_word.append("+" if mobius > 0 else "-")
        print(
            f"mask={mask} subset={''.join(subset):3s} charge={q} "
            f"chi7={chi7(q):+d} augmented={augmented_chi7(q):+d} mobius={mobius:+d}"
        )
    assert "".join(sign_word) == "++-+--+"
    print(f"mask-order sign word={''.join(sign_word)}")

    print("\nB. Cyclotomic recurrence and seven-periodicity")
    seq = cyclotomic_sequence(140)
    period = [ONE, ETA, -ONE, ETA_BAR, ONE, ZERO, ZERO]
    for r, value in enumerate(seq):
        assert value == period[r % 7]
    assert ETA * ETA + ETA + Quad(2) == ZERO
    assert ETA + ETA_BAR == -ONE
    assert ETA * ETA_BAR == Quad(2)
    print("F_r for r mod 7 = [1, eta, -1, eta_bar, 1, 0, 0]")
    print("exact recurrence and period checked through r=140")

    print("\nC. Address charge census")
    print(" n   r    M     c(0)   per-QR   per-NQR   vector q=0..6")
    for n in range(3, 22):
        r = n - 3
        census = charge_census(r)
        predicted = predicted_census(r)
        assert dict(census) == {q: predicted[q] for q in range(7) if predicted[q]}
        assert sum(census.values()) == comb(r + 2, 2)
        assert max(census[q] for q in range(7)) - min(census[q] for q in range(7)) <= 1
        assert weighted_census(census) == seq[r]
        vector = [census[q] for q in range(7)]
        print(
            f"{n:2d} {r:3d} {sum(vector):4d} {census[0]:8d} "
            f"{census[1]:8d} {census[3]:10d}   {vector}"
        )

    for r in range(0, 400):
        census = charge_census(r)
        predicted = predicted_census(r)
        assert all(census[q] == predicted[q] for q in range(7))
    print("formula checked exactly through r=399 (80,200 addresses on the last row)")

    print("\nD. First repeated-role neutral layer")
    for r in (1, 2):
        assert not [
            (a, b, r - a - b)
            for a in range(r + 1)
            for b in range(r - a + 1)
            if (a + 2 * b + 4 * (r - a - b)) % 7 == 0
        ]
    neutral3 = [(1, 1, 1)]
    neutral4 = [(0, 1, 3), (1, 3, 0), (3, 0, 1)]
    for r, expected in ((3, neutral3), (4, neutral4)):
        actual = sorted(
            (a, b, r - a - b)
            for a in range(r + 1)
            for b in range(r - a + 1)
            if (a + 2 * b + 4 * (r - a - b)) % 7 == 0
        )
        assert actual == expected
        print(f"depth={r} neutral (A,B,C) exponent triples={actual}")

    print("\nE. Fixed-depth operation-plane fibre")
    kernel = (2, 4, 1)
    assert sum(kernel) % 7 == 0
    assert sum(CHARGE[r] * d for r, d in zip(("A", "B", "C"), kernel)) % 7 == 0
    print("kernel direction in alpha+beta+gamma=constant is (2,4,1) mod 7")
    print("the charge records one parallel-line pencil, not all Kakeya directions")

    print("\nF. Tournament Analysis on the seven B3 strata")
    adj = [[0] * 7 for _ in range(7)]
    for u in range(7):
        for v in range(7):
            if u != v and (v - u) % 7 in QR:
                adj[u][v] = 1
    scores = [sum(row) for row in adj]
    assert scores == [3] * 7
    assert all(adj[q][(q + 1) % 7] for q in range(6))
    c3 = directed_3cycles(adj)
    scc = scc_sizes(adj)
    hp = hamiltonian_paths(adj)
    assert (c3, scc, hp) == (14, [7], 189)
    print(f"score_hist={{3: 7}}, directed_3cycles={c3}, SCC_sizes={scc}, Hamiltonian_paths={hp}")
    print("tie path=0->1->2->3->4->5->6 (no ties occur)")

    print("\nVERDICT: ALL EXACT")


if __name__ == "__main__":
    main()
