#!/usr/bin/env python3
"""Finite-exact stopping probe for within-fibre owner/carry data at (1,2,5).

This uses the sharp local phase from THM-2246.  It is not asserted to extend
to a global covering row.  It only tests whether the 169-sheet, horizon-two
owner/carry/minimal-nonface package can itself force an obstruction.
"""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import comb


P = 13
Z = Fraction(325007, 700000)
TAU = Fraction(1, 14)
H, C1, C2, C3 = 1, 2 * P, 2 * P**2, 2 * P**5
QS = tuple(1 + P**2 * 700000 * i for i in range(1, 6))


def require(test, message):
    if not test:
        raise RuntimeError(message)


def norm(value):
    residue = value % 1
    return min(residue, 1 - residue)


def v13(value):
    answer = 0
    while value % P == 0:
        value //= P
        answer += 1
    return answer


def point(n):
    return (Z + n) / P**2


def obligation(n):
    x = point(n)
    return (
        norm(H * x) > Fraction(1, 7)
        and norm(C1 * x) >= TAU
        and all(norm(q * x) >= TAU for q in QS)
    )


def owners(n):
    x = point(n)
    return tuple(c for c in (C2, C3) if norm(c * x) < TAU)


def floor_fraction(value):
    return value.numerator // value.denominator


def potential(c, n):
    return floor_fraction(c * point(n))


def parent_potential(c, h):
    return floor_fraction((c // P) * (Z + h) / P)


def terminal_potential(c):
    return floor_fraction((c // P**2) * Z)


def marked_increments(c, marked):
    increments = []
    gaps = []
    for i, n in enumerate(marked):
        nxt = marked[i + 1] if i + 1 < len(marked) else marked[0] + P**2
        increments.append(potential(c, nxt) - potential(c, n))
        gaps.append(nxt - n)
    return tuple(increments), tuple(gaps)


def main():
    require(tuple(map(v13, (C1, C2, C3))) == (1, 2, 5), "profile drift")
    require(H % 2 and H % P and all(q > 0 and q % P for q in QS), "scope")
    marked = tuple(n for n in range(P**2) if obligation(n))
    require(len(marked) == 112, "THM-2246 hostile occupancy drift")

    for n in range(P**2):
        x = point(n)
        require(norm(H * x) != Fraction(1, 7), "guard boundary")
        require(norm(C1 * x) != TAU, "peeled boundary")
        require(all(norm(q * x) != TAU for q in QS), "unit boundary")
        require(all(norm(c * x) != TAU for c in (C2, C3)), "owner boundary")

    owner_words = tuple(owners(n) for n in marked)
    require(set(owner_words) == {(C2,)}, "persistent owner lost")
    require(norm(2 * Z) < TAU < norm(2 * P**3 * Z), "terminal gate drift")

    # The witness stalk uses a single common root word n=h+13k.  Carries for
    # the two owners are not chosen independently at individual obligations.
    stalk = []
    for n in marked:
        h, k = n % P, n // P
        row = [n, h, k, C2]
        for c in (C2, C3):
            leaf = potential(c, n)
            parent = parent_potential(c, h)
            terminal = terminal_potential(c)
            slope = c // P**2
            require(parent == terminal + slope * h, "parent carry drift")
            require(leaf == parent + P * slope * k, "leaf carry drift")
            row.extend((leaf, leaf % P))
        stalk.append(tuple(row))

    # This one legal local stalk services every marked obligation by c2.
    # Therefore any complex formed as a union over legal stalk facets already
    # contains the full marked simplex; enumeration of 2^112 faces is needless.
    require(all(word == (C2,) for word in owner_words), "facet drift")
    minimal_nonfaces = ()
    obstruction_fraction = Fraction(0, len(marked))

    holonomy = {}
    for c in (C2, C3):
        slope = c // P**2
        increments, gaps = marked_increments(c, marked)
        defects = tuple(a - slope * b for a, b in zip(increments, gaps))
        require(sum(increments) == c, "integer winding drift")
        require(set(defects) == {0}, "nonzero holonomy defect")
        holonomy[c] = (slope, slope % P, sum(increments), sum(defects))

    switches = sum(
        owner_words[i] != owner_words[(i + 1) % len(owner_words)]
        for i in range(len(owner_words))
    )
    require(switches == 0, "owner switch drift")
    first_digits = Counter(n % P for n in marked)
    carry_pairs = Counter((row[5], row[7]) for row in stalk)
    serialization = "\n".join(",".join(map(str, row)) for row in stalk)

    print("scope=FINITE-EXACT local hostile control; not a global cover")
    print("probe=profile_(1,2,5)_within_fibre_horizon_two")
    print(f"z={Z} H={H} c1={C1} c2={C2} c3={C3}")
    print(f"unit_speeds={QS}")
    print(
        "universe=n in [0,168], x_n=(z+n)/169; "
        "filters=||x||>1/7, ||26x||>=1/14, five ||q_i x||>=1/14"
    )
    print(f"obligations={len(marked)} indices={marked}")
    print(f"by_first_digit={dict(sorted(first_digits.items()))}")
    print(
        f"terminal_norms=owner2:{norm(2 * Z)},"
        f"owner3:{norm(2 * P**3 * Z)}"
    )
    print(
        f"owner_words={dict(Counter(owner_words))} "
        f"common_stalks=1 switches={switches}"
    )
    print(f"carry_pair_distribution={dict(sorted(carry_pairs.items()))}")
    print(f"stalk_sha256={sha256(serialization.encode()).hexdigest()}")
    print(
        f"complex=full_{len(marked)-1}_simplex faces={2**len(marked)} "
        f"minimal_nonfaces={minimal_nonfaces}"
    )
    print(
        f"all_pairs={comb(len(marked), 2)} "
        f"all_triples={comb(len(marked), 3)} "
        f"obstruction_fraction={obstruction_fraction}"
    )
    for c, (slope, residue, winding, defect) in holonomy.items():
        print(
            f"carry_c={c} sheet_slope={slope} "
            f"residue_increment={residue} "
            f"raw_integer_winding={winding} holonomy_defect={defect}"
        )
    print("verdict=full-simplex collapse at this local horizon")


if __name__ == "__main__":
    main()
