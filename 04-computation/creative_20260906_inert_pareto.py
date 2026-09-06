#!/usr/bin/env python3
"""Exact all-height inert-pair envelope: finite heads plus proved tails.

No producer imports. Fractions and always-active checks remain valid under -O.
The proof of the height cutoffs is in the matching result note.
"""
from fractions import Fraction as F
from math import gcd, isqrt


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def factors(n):
    answer = []
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            answer.append((p, e))
        p += 1
    if n > 1:
        answer.append((n, 1))
    return answer


def eligible(p, q):
    return (0 < p < q and gcd(p, q) == 1 and gcd(p * q, 6) == 1
            and all(r % 3 == 2 and e <= 2 for r, e in factors(p + q)))


def bernoulli_mass(p, q):
    u = F((7 + q - p) % 14, 14)
    v = F((7 + q + p) % 14, 14)
    return F(2, 49) + 2 * (u*u - u - v*v + v) / (p*q)


def intervals(p, q):
    """Intersect literal opposite-parity danger intervals in y=2x.

    Endpoints are excluded. Distinct teeth of either family are separated,
    so every positive intersection is an entire connected component.
    """
    def teeth(s):
        return [(max(F(0), F(7*k-1, 7*s)),
                 min(F(1), F(7*k+1, 7*s)), k % 2)
                for k in range(s+1)]
    a, b = teeth(p), teeth(q)
    i = j = 0
    out = []
    while i < len(a) and j < len(b):
        left, right = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if left < right and a[i][2] != b[j][2]:
            out.append((left, right))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def profile(p, q):
    arcs = intervals(p, q)
    widths = [b-a for a, b in arcs]
    return sum(widths, F(0)), max(widths, default=F(0))


FRONTIER = (
    (7, 13, F(2, 49), F(1, 49)),
    (19, 25, F(138, 3325), F(37, 3325)),
    (5, 41, F(12, 287), F(2, 287)),
    (5, 53, F(78, 1855), F(2, 371)),
    (1, 67, F(20, 469), F(2, 469)),
)


def main():
    cap = F(20, 469)
    require(cap - F(2, 49) == F(6, 3283), "mass gap")
    require(F(2, 49) + F(1, 2*274) < cap, "all-height mass cutoff")
    mass_head = [(p, q) for p in range(1, isqrt(273)+1)
                 for q in range(p+1, 273//p+1) if eligible(p, q)]
    require(len(mass_head) == 21, "complete product head")
    maxima = [(p, q) for p, q in mass_head if bernoulli_mass(p, q) == cap]
    require(maxima == [(1, 67)], "unique mass maximizer")
    require(all(bernoulli_mass(p, q) <= cap for p, q in mass_head), "head cap")
    for p, q in mass_head:
        require(profile(p, q)[0] == bernoulli_mass(p, q), "literal mass head")

    require(F(2, 7*15) < F(1, 49), "component cutoff")
    width_head = [(p, q) for q in range(2, 15) for p in range(1, q)
                  if eligible(p, q)]
    require(width_head == [(7, 13)], "complete width head")
    require(intervals(7, 13) == [(F(1, 7), F(8, 49)),
                               (F(41, 49), F(6, 7))], "width equality intervals")

    pareto_head = [(p, q, *profile(p, q)) for q in range(2, 68)
                   for p in range(1, q) if eligible(p, q)]
    maximal = [r for r in pareto_head
               if not any(s[2] >= r[2] and s[3] >= r[3]
                          and s[2:] != r[2:] for s in pareto_head)]
    require(set(maximal) == set(FRONTIER), "complete Pareto frontier")
    require(all(any(m <= a and b <= c for _, _, a, c in FRONTIER)
                for _, _, m, b in pareto_head), "all head profiles dominated")
    require(all(profile(p, q) == (m, b) for p, q, m, b in FRONTIER), "attainment")

    # Formula versus literal geometry, including excluded arithmetic classes.
    controls = 0
    for q in range(3, 76, 2):
        for p in range(1, q, 2):
            if gcd(p, q) != 1:
                continue
            m, b = profile(p, q)
            require(m == bernoulli_mass(p, q), "formula/geometry control")
            require(b <= F(2, 7*q), "component upper bound")
            controls += 1
    require(profile(1, 9)[0] == F(4, 63) > cap, "3-unit hypothesis hostile")
    require(profile(5, 11)[0] == F(18, 385) > cap, "exponent hypothesis hostile")
    require(profile(1, 11)[0] == F(4, 77) > cap, "inertness hostile")
    require(profile(1, 3) == (F(0), F(0)), "empty-carrier control")
    require(2*cap == F(40, 469) < F(8, 91), "q4 gate")
    require(cap+F(8, 63) == F(716, 4221) < F(20, 117), "q2 gate")

    # The actual decoder atlas is only a subset of this all-height theorem.
    atlas = [(p, q) for q in range(2, 356) for p in range(1, q)
             if p+q <= 356 and eligible(p, q)]
    require(len(atlas) == 548, "odd-3-unit inert atlas")
    require(all((p, q) in atlas for p, q, _, _ in FRONTIER), "atlas attainment")

    print("INERT_PAIR_ALL_HEIGHT_PARETO")
    print(f"mass_head={len(mass_head)} product<=273")
    for p, q in mass_head:
        print(f"mass_head_pair={p},{q} mass={bernoulli_mass(p,q)}")
    print(f"pareto_head={len(pareto_head)} q<=67")
    for p, q, m, b in FRONTIER:
        print(f"frontier={p},{q} mass={m} component={b}")
    print(f"literal_formula_controls={controls} q<=75")
    print(f"actual_atlas_subclass={len(atlas)} p+q<=356")
    print("sharp_mass=20/469 unique_pair=1,67")
    print("sharp_component=1/49 unique_pair=7,13")
    print("inclusive_q4_body_gate=40/469 inclusive_q2_one_even_gate=716/4221")
    print("hypothesis_hostiles=1,9;5,11;1,11 empty_control=1,3")
    print("PASS")


if __name__ == "__main__":
    main()
