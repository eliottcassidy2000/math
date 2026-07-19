#!/usr/bin/env python3
"""Exact referee for the four-torsion centre-hit criterion (codex-S79).

For an integer direction d, write g = gcd(d) and e = d/g.  The closed
geodesic gamma_d(u) = ({-d_i u}) hits the four-torsion centre r/4, where r
is a permutation of (1,2,3), iff e is congruent to r or -r modulo 4.

This refutes the proposed assertion that only directions proportional to
(1,2,3) can hit a centre.  The proof in THM-1211 is algebraic; this script
only replays the equivalence exhaustively and freezes explicit witnesses.

Tournament analysis is intentionally not imposed here.  The preserved
predicate is membership of one four-torsion point in one cyclic subtorus,
and its essential datum is the *single common gauge* a modulo 4 across all
three coordinates.  Pairwise orientations (whether on frequencies, centre
slots, residues, or proof obligations) forget that common gauge, so score
histograms, cycles, SCCs, and Hamiltonian-path counts would be artificial
telemetry rather than evidence for the claim.
"""

from fractions import Fraction
from itertools import combinations, permutations
from math import gcd


CENTRES = tuple(permutations((1, 2, 3)))


def gcd3(d: tuple[int, int, int]) -> int:
    return gcd(gcd(d[0], d[1]), d[2])


def primitive(d: tuple[int, int, int]) -> tuple[int, int, int]:
    g = gcd3(d)
    return tuple(x // g for x in d)


def frac(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def direct_hit(
    d: tuple[int, int, int], r: tuple[int, int, int], u: Fraction
) -> bool:
    return all(frac(-di * u) == Fraction(ri, 4) for di, ri in zip(d, r))


def common_gauge(
    d: tuple[int, int, int], r: tuple[int, int, int]
) -> int | None:
    """Return a mod 4 with -a(d/g)=r mod 4, if one exists."""
    e = primitive(d)
    for a in range(4):
        if all((-a * ei - ri) % 4 == 0 for ei, ri in zip(e, r)):
            return a
    return None


def signed_residue_criterion(
    d: tuple[int, int, int], r: tuple[int, int, int]
) -> bool:
    e = primitive(d)
    return all((ei - ri) % 4 == 0 for ei, ri in zip(e, r)) or all(
        (ei + ri) % 4 == 0 for ei, ri in zip(e, r)
    )


def main() -> None:
    explicit_d = (1, 2, 7)
    explicit_r = (1, 2, 3)
    explicit_u = Fraction(3, 4)
    assert direct_hit(explicit_d, explicit_r, explicit_u)
    assert primitive(explicit_d) != (1, 2, 3)

    rows = 0
    hits = 0
    nonproportional_hits = 0
    for d in combinations(range(1, 65), 3):
        g = gcd3(d)
        e = primitive(d)
        for r in CENTRES:
            rows += 1
            a = common_gauge(d, r)
            signed = signed_residue_criterion(d, r)
            assert (a is not None) == signed
            if a is not None:
                # The proof shows that any centre hit has a preimage of this
                # form; conversely this u is an exact witness.
                u = Fraction(a, 4 * g)
                assert direct_hit(d, r, u)
                assert a in (1, 3)  # r has odd coordinates.
                hits += 1
                if e != r and e != tuple(ri * e[0] for ri in (1, 2, 3)):
                    nonproportional_hits += 1

    family_rows = 64
    for m in range(1, family_rows + 1):
        d = (1, 2, 4 * m + 3)
        assert direct_hit(d, (1, 2, 3), Fraction(3, 4))
        assert d != (1, 2, 3)

    print("THM-1211 four-torsion centre-hit congruence referee")
    print(f"complete sorted triples 1<=d1<d2<d3<=64, labelled centres: {rows}")
    print(f"exact centre hits: {hits}")
    print(f"hits outside the literal proportional (1,2,3) ray: {nonproportional_hits}")
    print("criterion mismatches: 0")
    print("explicit refutation: d=(1,2,7), r=(1,2,3), u=3/4")
    print("gamma_d(u)=(1/4,1/2,3/4) exactly")
    print(f"infinite family prefix checked: (1,2,4m+3), 1<=m<={family_rows}")
    print("tournament analysis: inapplicable; the common Z/4 gauge is ternary/global")


if __name__ == "__main__":
    main()
