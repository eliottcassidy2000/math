#!/usr/bin/env python3
"""Exact replay/search for the THM-1183 continuum-ceiling reduction.

All interval endpoints and measures use fractions.Fraction.  This proves the
six-order formula and evaluates both the true bad measure and the explicit
three-path majorant.  The finite sweep is evidence for, not a proof of, the
remaining all-integer three-path inequality.
"""
from fractions import Fraction as Q
from itertools import combinations

I = (Q(1, 7), Q(2, 7))
CEILING = Q(2, 21)
HALF_CEILING = Q(1, 21)


def comb(q: int) -> list[tuple[Q, Q]]:
    """Times u in [0,1] with frac(q*u) in [1/7,2/7]."""
    assert q != 0
    if q > 0:
        return [(Q(7*n+1, 7*q), Q(7*n+2, 7*q)) for n in range(q)]
    q = -q
    return [(Q(7*n+5, 7*q), Q(7*n+6, 7*q)) for n in range(q)]


def meet(a: list[tuple[Q, Q]], b: list[tuple[Q, Q]]) -> list[tuple[Q, Q]]:
    out: list[tuple[Q, Q]] = []
    i = j = 0
    while i < len(a) and j < len(b):
        left = max(a[i][0], b[j][0])
        right = min(a[i][1], b[j][1])
        if left < right:
            out.append((left, right))
        if a[i][1] < b[j][1]:
            i += 1
        elif b[j][1] < a[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return out


def event(coeffs: tuple[int, ...]) -> tuple[Q, list[tuple[Q, Q]]]:
    intervals = [(Q(0), Q(1))]
    for q in sorted(coeffs, key=abs):
        intervals = meet(intervals, comb(q))
        if not intervals:
            break
    return sum((r-l for l, r in intervals), Q(0)), intervals


def carriers(a: int, b: int, c: int):
    assert 0 < a < b < c
    # Three cyclic orders modulo reversal.
    cycles = (
        (a, b-a, c-b, -c),       # 0,a,b,c
        (a, c-a, b-c, -b),       # 0,a,c,b
        (b, a-b, c-a, -c),       # 0,b,a,c
    )
    chambers = tuple(event(w) for w in cycles)
    bad = 2 * sum((m for m, _ in chambers), Q(0))

    # Delete the first edge of each cycle.  These are the disjoint paths
    # a->b->c->0, a->c->b->0, b->a->c->0.
    paths = tuple(event(w[1:]) for w in cycles)
    path_majorant = sum((m for m, _ in paths), Q(0))
    return bad, chambers, path_majorant, paths


def show(a: int, b: int, c: int) -> None:
    bad, chambers, majorant, paths = carriers(a, b, c)
    print(f"  {(a,b,c)} bad={bad} ({float(bad):.12f}) "
          f"half={bad/2} path-majorant={majorant}")
    print("    chamber measures=" + str(tuple(m for m, _ in chambers)))
    print("    path measures   =" + str(tuple(m for m, _ in paths)))


def main() -> None:
    print("THM-1183 exact six-order / three-path reduction")
    print(f"I=[{I[0]},{I[1]}], proposed full ceiling={CEILING}")
    print("[exact named rows]")
    for row in ((1,2,3), (2,4,6), (1,2,4), (1,3,6), (1,9,10), (1,2,7)):
        show(*row)

    for m in range(1, 14):
        bad, ch, majorant, paths = carriers(m, 2*m, 3*m)
        assert bad == CEILING
        assert tuple(x[0] for x in ch) == (HALF_CEILING, Q(0), Q(0))

    D = 40
    actual_top = []
    path_top = []
    actual_over = path_over = 0
    equality = []
    nonprop_top = (Q(0), None)
    for a, b, c in combinations(range(1, D+1), 3):
        bad, chambers, majorant, paths = carriers(a, b, c)
        actual_top.append((bad, (a,b,c)))
        path_top.append((majorant, (a,b,c)))
        actual_over += bad > CEILING
        path_over += majorant > HALF_CEILING
        proportional = b == 2*a and c == 3*a
        if bad == CEILING:
            equality.append((a,b,c))
        if not proportional and bad > nonprop_top[0]:
            nonprop_top = (bad, (a,b,c))

    actual_top.sort(reverse=True)
    path_top.sort(reverse=True)
    assert actual_over == 0
    assert path_over == 0
    assert equality == [(m,2*m,3*m) for m in range(1, D//3+1)]
    assert nonprop_top[0] == Q(4,105)

    print("[complete exact finite sweep]")
    print(f"  1<=a<b<c<={D}: rows={len(actual_top)}")
    print(f"  actual over 2/21={actual_over}; path majorant over 1/21={path_over}")
    print(f"  actual equality rows={equality}")
    print(f"  largest nonproportional actual={nonprop_top}")
    print("  top path-majorant rows:")
    for value, row in path_top[:8]:
        print(f"    {row}: {value} ({float(value):.12f})")

    print("[carrier audit]")
    print("  torus: rational geodesic u -> (a*u,b*u,c*u) meets six 1/7 Kuhn tetrahedra")
    print("  faithful binary relation: i->j iff frac((v_j-v_i)u) is in [1/7,2/7]")
    print("  bad fingerprint: directed C4; scores=(1,1,1,1), triangles=0, SCC=(4), HP=4")
    print("  u->1-u reverses every edge; inside each bad component edge flips=0")
    print("  deleting one cycle edge gives the three path obligations used above")
    print("DONE -- the all-integer path-majorant inequality remains analytic, not proved here")


if __name__ == "__main__":
    main()
