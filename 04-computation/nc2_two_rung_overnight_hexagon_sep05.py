#!/usr/bin/env python3
"""Exact referee controls for the all-carry two-first-channel certificate.

No repository mathematics or floating point is imported. The proof universe
is declared in main; original charge equations independently reconstruct both
moment rows before the parameter/carry formulas and Bezout identity are used.
Every gate remains active under Python optimization.
"""

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import comb, factorial, gcd, prod


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def multinomial(v):
    return factorial(sum(v)) // prod(factorial(a) for a in v)


def channels(charges, mass, weights=True):
    """Original equations: scan negative count, solve the two positive counts."""
    a, b, c = charges
    result = {}
    for x in range(mass + 1):
        numerator = a*x - b*(mass-x)
        if numerator % (c-b):
            continue
        z = numerator // (c-b)
        y = mass-x-z
        if min(y, z) >= 0:
            result[x, y, z] = multinomial((x, y, z)) if weights else 1
    return result


def plus(t, L):
    return F(comb(2*t+L, t), comb(2*t, t))


def minus(t, L):
    return F(comb(2*t+L, t), comb(2*t+2*L, t+L))


def parameter_audit(A, B, r, z, x):
    delta = B-A
    a = A*(B+r)+B*z
    g = x+B+r+z
    b, c = g*A-a, g*B-a
    require(b > 0 and c > b and gcd(g, a) == 1, "eligible primitive parameter row")
    require(gcd(gcd(a, b), c) == 1 and gcd(a+b, a+c) == g, "charge content and return modulus")
    require(all(a*m % g for m in range(1, g)), "every earlier mass ruled out")
    v, w = (x, B+r, z), (x+delta, r, z+A)
    first = channels((a, b, c), g)
    require(set(first) == {v, w}, ("complete first row", (a, b, c), first))
    p, q = first[v], first[w]
    second = channels((a, b, c), 2*g)
    eps_y, eps_z = 2*r//B, 2*z//A
    indices = range(-eps_z, 3+eps_y)
    expected = {tuple(2*v[i]+j*(w[i]-v[i]) for i in range(3)): j for j in indices}
    require(set(second) == set(expected), ("complete second row including carries", (a, b, c)))
    by_index = {j: second[u] for u, j in expected.items()}

    Nv = prod(comb(2*t, t) for t in v)
    Nw = prod(comb(2*t, t) for t in w)
    mixed = prod(comb(v[i]+w[i], v[i]) for i in range(3))
    require(F(mixed, Nv) == plus(x, delta)*minus(r, B)*plus(z, A), "first ratio decomposition")
    require(F(mixed, Nw) == minus(x, delta)*plus(r, B)*minus(z, A), "second ratio decomposition")
    require(2*mixed < Nv and 2*mixed <= Nw, "strict and weak half bounds")
    require((2*mixed == Nw) == (A == 1 and B == 2 and r == 1 and z == 0),
            "exact equality boundary for the weak ratio")

    C, D, E = (by_index[j] for j in (0, 1, 2))
    Delta = C*q*q-D*p*q+E*p*p
    require(F(Delta, p*p*q*q) == comb(2*g, g)*(F(1, Nv)-F(1, mixed)+F(1, Nw)),
            "normalized factorial determinant")
    require(Delta < 0, "generated three-channel determinant is strictly negative")
    left, right = by_index.get(-1, 0), by_index.get(3, 0)
    require((left > 0) == bool(eps_z) and (right > 0) == bool(eps_y), "both carry labels")
    t = -F(p, q)
    normalized_second = sum(weight*t**j for j, weight in by_index.items())
    require(normalized_second < 0, "all channels at the tuned first zero")
    require(all(weight*t**j < 0 for j, weight in by_index.items() if j in (-1, 3)),
            "both extra carries have favorable sign")
    K = -p*q*Delta+left*q**4+right*p**4
    require(K > 0 and q**4*t*normalized_second == K, "positive exact remainder")

    # Coefficients in increasing powers of t of the cubic quotient Q(1,t).
    Q = [q**3*C-p*q*q*D+p*p*q*E-p**3*right,
         q**3*D-p*q*q*E+p*p*q*right,
         q**3*E-p*q*q*right,
         q**3*right]
    G = [left, C, D, E, right]
    difference = [q**4*c for c in G]
    for j, coefficient in enumerate(Q):
        difference[j] -= p*coefficient
        difference[j+1] -= q*coefficient
    require(difference == [K, 0, 0, 0, 0], "full integer polynomial Bezout identity")
    require(2*g <= a+c, "sharp-width consumer")
    require((2*g == a+c) == (A == 1 and B == 2), "width equality boundary")
    return ((a, b, c), g, (eps_y, eps_z), v, w, p, q, by_index, Delta,
            normalized_second, K)


def recover(charges):
    a, b, c = charges
    g = gcd(a+b, a+c)
    A, B = (a+b)//g, (a+c)//g
    z = next(j for j in range(A) if (a-B*j) % A == 0)
    y = (a-B*z)//A
    h, r = divmod(y, B)
    require(h == 1, ("two-channel control hypothesis", charges, h))
    return A, B, r, z, g-B-r-z


def main():
    # Independent exact arithmetic for the analytic scalar inequalities.
    ratio_rows = 0
    for t in range(41):
        for L in range(1, 62):
            require(plus(t, L) < 2**L, "unrestricted plus bound")
            require(minus(t, L) <= F(1, 2**L), "minus bound")
            if L > t:
                require(plus(t, L) <= 2**(L-1), "canonical-residue plus bound")
                require((plus(t, L) == 2**(L-1)) == ((t, L) in ((0, 1), (1, 2))),
                        "scalar equality boundary")
            ratio_rows += 1
    require(plus(1, 1) > 1, "hostile to dropping strict L>t")

    kinds, nonfree, rows = Counter(), 0, 0
    digest = sha256()
    for B in range(2, 10):
        for A in range(1, B):
            if gcd(A, B) != 1:
                continue
            for r in range(B):
                for z in range(A):
                    a = A*(B+r)+B*z
                    for x in range(1, 11):
                        g = x+B+r+z
                        if A*x-(B-A)*z <= 0 or gcd(a, g) != 1:
                            continue
                        data = parameter_audit(A, B, r, z, x)
                        rows += 1
                        nonfree += bool(r or z)
                        kinds[data[2]] += 1
                        digest.update(repr(data).encode())
    require(rows == 3665, "complete parameter universe")
    require(kinds == {(0, 0): 1187, (1, 0): 999, (0, 1): 794, (1, 1): 685},
            "all four doubling-carry classes")
    print("PARAMETER_UNIVERSE 1<=A<B<=9; gcd(A,B)=1; 0<=r<B; 0<=z<A; 1<=x<=10; b>0; gcd(a,g)=1")
    print("ROWS", rows, "NONFREE", nonfree, "CARRY_CLASSES", sorted(kinds.items()))
    print("SCALAR_RATIO_ROWS", ratio_rows)
    for charges in ((3, 1, 9), (4, 1, 11), (3, 1, 5), (13, 1, 8)):
        data = parameter_audit(*recover(charges))
        print("CONTROL", charges, "g", data[1], "carries", data[2],
              "v,w", data[3:5], "weights", (data[5], data[6], sorted(data[7].items())),
              "Delta", data[8], "tuned_second", str(data[9]), "K", data[10])

    # Infinite families are proved in the note; these are explicit wide controls.
    family_rows = 0
    for g in range(5, 102, 2):
        data = parameter_audit(*recover((4, g-4, 3*g-4)))
        require(data[2] == (0, 0), "nonfree family has no doubling carry")
        require(len(channels(data[0], 3*g)) == 5, "first later carry at the third rung")
        family_rows += 1
    unbounded_endpoint_rows = 0
    for n in range(1, 21):
        g = (n+1)**2
        support = (n*(n+2), n*(n*n+n-1), n*(n+1)**2+1)
        data = parameter_audit(*recover(support))
        require(data[1] == g, "unbounded-endpoint nonfree family first mass")
        require(len(channels(support, n*g, weights=False)) == n+1,
                "no extra channel before the first later carry")
        require(len(channels(support, (n+1)*g, weights=False)) == n+3,
                "first later carry in the unbounded-endpoint family")
        unbounded_endpoint_rows += 1
    sharp_rows = 0
    for g in range(4, 81):
        if gcd(3, g) != 1:
            continue
        data = parameter_audit(*recover((3, g-3, 2*g-3)))
        require(data[2] == (1, 0) and 2*g == sum((data[0][0], data[0][2])),
                "primitive width-three sharp family")
        sharp_rows += 1
    hostile = channels((4, 1, 6), 5)
    require(len(hostile) == 3, "higher first-row multiplicity is not a two-channel row")
    print("NONFREE_LATER_CARRY_FAMILY_ROWS", family_rows,
          "UNBOUNDED_ENDPOINT_FAMILY_ROWS", unbounded_endpoint_rows,
          "SHARP_WIDTH_THREE_FAMILY_ROWS", sharp_rows)
    print("HIGHER_MULTIPLICITY_HOSTILE", sorted(hostile.items()))
    print("SEMANTIC_SHA256", digest.hexdigest())
    print("EXPLICIT_GATES", CHECKS, "RESULT PASS; all-first-row-multiplicity>=3 remains OPEN")


if __name__ == "__main__":
    main()
