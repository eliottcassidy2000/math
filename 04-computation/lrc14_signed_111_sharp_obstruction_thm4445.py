#!/usr/bin/env python3
"""Exact audit and threshold anatomy for signed (1,1,1) LRC14 tails.

The only sorted positive signed (1,1,1) relation is a+b=c.  This script
implements the resulting deleted-third carrier ray without importing repo
producers, checks it against a literal raw carrier scan, reproduces the sharp
finite heads, and audits the uniform 6/77 obstruction derived from two-sided
deleted-third quadrature.
"""

from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from math import gcd

R = Q(3, 14)
OLD = Q(6, 77)
SHARP = Q(6, 55)
CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def eligible(w):
    a, b, c = w
    return (0 < a < b < c and a + b == c and gcd(gcd(a, b), c) == 1
            and all(x % 3 for x in w))


def ray_row(w):
    """Exact normalized ray profile, projections, physical mass, carriers."""
    a, b, c = w
    need(eligible(w), ("typed row", w))
    t = Q(a, c)
    addresses = tuple(k for k in range(1, (3*c - 1)//14 + 1) if k % 3)
    samples = []
    for k in addresses:
        s = Q(k, c)
        vals = (
            min(2*R, (R*(2-t)-s)/(1-t)),
            min(2*R, (R*(1+t)-s)/t),
            min(2*R, (R-s)/(t*(1-t))),
        )
        need(min(vals) > 0, ("positive strict ray sample", w, k, vals))
        samples.append(vals)
    E = tuple(Q(2, c)*sum((v[i] for v in samples), Q(0))
              for i in range(3))
    physical = Q(2, c)*sum((min(v) for v in samples), Q(0))
    carriers = {(sgn*k, sgn*k, -sgn*k)
                for k in addresses for sgn in (-1, 1)}
    return E, physical, carriers, addresses


def raw_row(w):
    """Direct roof/congruence scan, independent of the one-ray assertion."""
    a, b, c = w
    bounds = tuple((3*(sum(w)-w[i])-1)//14 for i in range(3))
    live = set()
    for x in range(-bounds[0], bounds[0]+1):
        if x % 3 == 0:
            continue
        for y in range(-bounds[1], bounds[1]+1):
            if y % 3 == 0:
                continue
            q = -(a*x+b*y)
            if q % c:
                continue
            z = q//c
            if z % 3 and abs(z) <= bounds[2]:
                live.add((x, y, z))
    E = [Q(0)]*3
    physical = Q(0)
    for C in live:
        vals = []
        for i in range(3):
            j, k = [h for h in range(3) if h != i]
            v = min(Q(3, 7*c),
                    Q(3*(w[j]+w[k])-14*abs(C[i]), 14*w[j]*w[k]))
            need(v > 0, ("positive raw carrier", w, C, i, v))
            vals.append(v)
            E[i] += v
        physical += min(vals)
    return tuple(E), physical, live


def parity_mask(w):
    return sum((x % 2) << i for i, x in enumerate(w))


def update(table, key, value, payload):
    old = table.get(key)
    if old is None or value > old[0]:
        table[key] = (value, [payload])
    elif value == old[0]:
        old[1].append(payload)


def update_min(table, key, value, payload):
    old = table.get(key)
    if old is None or value < old[0]:
        table[key] = (value, [payload])
    elif value == old[0]:
        old[1].append(payload)


def qstr(x):
    return f"{x.numerator}/{x.denominator}"


def tail_entry(bulk, candidate):
    """First integer c for which bulk+4/(7c) is strictly below candidate."""
    c = 1
    while bulk + Q(4, 7*c) >= candidate:
        c += 1
    return c


def lower_tail_entry(candidate):
    """First integer c where 9/98-4/(7c) is strictly above candidate."""
    c = 1
    while Q(9, 98) - Q(4, 7*c) <= candidate:
        c += 1
    return c


def main():
    digest = sha256()
    leaders = {}
    old_threshold = defaultdict(lambda: {"lt": 0, "eq": [], "gt": 0})
    last_not_gt = {"N": None, "P": None}
    head_rows = 0
    old_boundary_rows = 0
    lower_boundary_rows = 0
    raw_rows = 0
    family_rows = []
    finite_old_boundary = {"N": [], "P": []}
    low_nonexception = {}

    # c<=223 is enough not only for the old-bound boundary, but also for
    # every parity-stratum maximum: the universal continuum/error tail is
    # already below each exact head leader after the cutoffs checked below.
    for c in range(1, 224):
        if c % 3 == 0:
            continue
        for a in range(1, (c+1)//2):
            b = c-a
            w = (a, b, c)
            if not eligible(w):
                continue
            E, physical, carriers, addresses = ray_row(w)
            Er, pr, cr = raw_row(w)
            need((E, physical, carriers) == (Er, pr, cr),
                 ("ray/raw mismatch", w, E, physical, carriers, Er, pr, cr))
            raw_rows += 1
            need(E[1] >= E[0] and physical <= min(E),
                 ("projection/physical order", w, E, physical))
            t = Q(a, c)
            A = Q(9, 98)+Q(3, 98)*t
            C = Q(12, 98)*(1-t+t*t)
            L = min(A, C)
            P = Q(9, 98)
            err = Q(4, 7*c)
            need(L-err <= min(E) < L+err,
                 ("two-sided network quadrature", w, min(E), L, err))
            need(P-err <= physical < P+err,
                 ("two-sided physical quadrature", w, physical, P, err))
            if c >= 43:
                need(P-err > OLD and min(E) > OLD and physical > OLD,
                     ("uniform old-bound obstruction", w, E, physical, P-err))
            for kind, value in (("N", min(E)), ("P", physical)):
                rel = "lt" if value < OLD else "eq" if value == OLD else "gt"
                if rel == "eq":
                    old_threshold[(kind, parity_mask(w))][rel].append(w)
                else:
                    old_threshold[(kind, parity_mask(w))][rel] += 1
                if value <= OLD:
                    last_not_gt[kind] = w
                if w != (1, 4, 5):
                    update_min(low_nonexception, kind, value, w)
            if c <= 41:
                old_boundary_rows += 1
                for kind, value in (("N", min(E)), ("P", physical)):
                    if value <= OLD:
                        finite_old_boundary[kind].append((w, value))
            if c <= 44:
                lower_boundary_rows += 1
            # Candidate maxima over a head long enough for every tail cutoff.
            if c <= 223:
                head_rows += 1
                update(leaders, ("N", "all"), min(E), w)
                update(leaders, ("N", parity_mask(w)), min(E), w)
                update(leaders, ("P", "all"), physical, w)
                update(leaders, ("P", parity_mask(w)), physical, w)
            if a % 3 == 1 and b == 3*a+1 and c == 4*a+1:
                family_rows.append((w, min(E), physical, L, err))
            digest.update(repr((w, E, physical, tuple(sorted(carriers)))).encode())

    need(finite_old_boundary == {
        "N": [((1, 4, 5), Q(1, 28))],
        "P": [((1, 4, 5), Q(1, 28))],
    }, ("complete low-height old-bound boundary", finite_old_boundary))
    need(leaders[("N", "all")] == (SHARP, [(1, 10, 11)]),
         ("sharp network equality", leaders[("N", "all")]))
    need(leaders[("P", "all")] == (SHARP, [(1, 10, 11)]),
         ("sharp physical equality", leaders[("P", "all")]))
    expected_leaders = {
        ("N", "all"): (SHARP, [(1, 10, 11)]),
        ("N", 3): (Q(2946, 28861), [(7, 31, 38)]),
        ("N", 5): (SHARP, [(1, 10, 11)]),
        ("N", 6): (Q(223, 2156), [(4, 7, 11)]),
        ("P", "all"): (SHARP, [(1, 10, 11)]),
        ("P", 3): (Q(222, 2275), [(1, 25, 26)]),
        ("P", 5): (SHARP, [(1, 10, 11)]),
        ("P", 6): (Q(102, 1001), [(2, 11, 13)]),
    }
    need(leaders == expected_leaders, ("all parity-stratum maxima", leaders))
    need(low_nonexception == {
        "N": (Q(31, 392), [(1, 7, 8)]),
        "P": (Q(31, 392), [(1, 7, 8)]),
    }, ("sharp nonexception lower boundary", low_nonexception))
    need(Q(9, 98)-Q(4, 7*43) > OLD,
         ("uniform c=43 margin", Q(9, 98)-Q(4, 7*43)-OLD))
    for key, (candidate, witnesses) in leaders.items():
        bulk = Q(39, 392) if key[0] == "N" else Q(9, 98)
        cutoff = tail_entry(bulk, candidate)
        need(cutoff <= 224, ("parity-tail cutoff outside head", key, cutoff))
    for kind, (candidate, witnesses) in low_nonexception.items():
        need(lower_tail_entry(candidate) <= 224,
             ("lower-tail cutoff outside head", kind, candidate,
              lower_tail_entry(candidate)))

    print("SIGNED_111_EQUIVALENCE sorted positive signs force a+b=c")
    print("RAW_RAY_CROSSCHECK c<=223 rows", raw_rows)
    print("SHARP_ALL_HEIGHT network=physical_upper", SHARP,
          "unique_equality", (1, 10, 11))
    print("CONTINUUM network_range", Q(9, 98), Q(39, 392),
          "physical", Q(9, 98), "two_sided_error", "4/(7c)")
    print("UNIFORM_OLD_BOUND_FAILURE c>=43 margin_at_43",
          Q(9, 98)-Q(4, 7*43)-OLD)
    print("COMPLETE_C_LE_41 rows", old_boundary_rows,
          "AT_OR_BELOW_6_77", finite_old_boundary)
    print("LAST_ROWS_AT_OR_BELOW_6_77", last_not_gt)
    print("SHARP_NONEXCEPTION_LOWER_BOUNDS")
    for kind in sorted(low_nonexception):
        candidate, witnesses = low_nonexception[kind]
        print(kind, candidate, witnesses, "finite_rows_c_le_44", lower_boundary_rows,
              "tail_entry", lower_tail_entry(candidate))
    print("SHARP_PARITY_STRATUM_LEADERS_WITH_TAIL_CUTOFFS")
    for key in sorted(leaders, key=str):
        bulk = Q(39, 392) if key[0] == "N" else Q(9, 98)
        print(key, leaders[key], "tail_entry", tail_entry(bulk, leaders[key][0]))
    print("THRESHOLD_CENSUS_C_LE_223")
    for key in sorted(old_threshold, key=str):
        print(key, old_threshold[key])
    print("COFINAL_RATIO_ONE_QUARTER_FAMILY")
    for row in family_rows[:8]:
        print(row)
    if family_rows:
        print("last", family_rows[-1])
    print("SEMANTIC_SHA256", digest.hexdigest())
    print("EXPLICIT_CHECKS", CHECKS)
    print("PASS")


if __name__ == "__main__":
    main()
