"""Exact toric/conic controls, not an RH or tournament-connectivity test."""
from fractions import Fraction
from itertools import product
from math import comb, gcd


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def phi(t, r):
    x, y, z = t
    a, b = (z+x)**r, (z-x)**r
    need((a-b) % 2 == 0 and (a+b) % 2 == 0, "integral half")
    return ((a-b)//2, y**r, (a+b)//2)


def euclid(m, n):
    return (m*m-n*n, 2*m*n, m*m+n*n)


def canonical(values, p):
    for v in values:
        if v % p:
            inv = pow(v, -1, p)
            return tuple(x*inv % p for x in values)
    raise RuntimeError("zero projective coordinates")


def main():
    field_points = 0
    for p in (2, 3, 5, 7, 11, 13):
        source = {(1, n) for n in range(p)} | {(0, 1)}
        image = {canonical((m*m, n*n, m*n), p) for m, n in source}
        target = {canonical((u, v, w), p)
                  for u in range(p) for v in range(p) for w in range(p)
                  if (u or v or w) and (u*v-w*w) % p == 0}
        need(image == target and len(target) == p+1, "Veronese bijection")
        field_points += len(target)
        for t in target:
            for r in (2, 3, 5):
                got = canonical(tuple(x**r for x in t), p)
                need((got[0]*got[1]-got[2]**2) % p == 0, "power stability")
        # Original conic has p+1 geometric F_p points even at the bad p=2;
        # this count alone does NOT certify reducedness/model equality.
    print("Projective toric bijection: primes 2,3,5,7,11,13; points", field_points)

    samples = 0
    for m in range(2, 21):
        for n in range(1, m):
            if gcd(m, n) != 1 or (m+n) % 2 == 0:
                continue
            t = euclid(m, n)
            for r in range(1, 8):
                u = phi(t, r)
                expected = tuple(2**(r-1)*x for x in euclid(m**r, n**r))
                need(u == expected, "exact scale law")
                need(u[0]**2+u[1]**2 == u[2]**2, "conic equation")
                need(gcd(gcd(u[0], u[1]), u[2]) == 2**(r-1), "primitive scale")
                for s in (2, 3):
                    need(phi(u, s) == phi(t, r*s), "commutation")
                samples += 1
    print("Primitive parameter/power controls:", samples)

    for p in (3, 5, 7, 11, 13, 17, 19):
        for j in range(p+1):
            # x coordinate uses odd x powers; z coordinate even x powers.
            if j % 2:
                need(comb(p, j) % p == (1 if j == p else 0), "x Frobenius")
            else:
                need(comb(p, j) % p == (1 if j == 0 else 0), "z Frobenius")
    need(phi((3, 4, 5), 3) == (252, 64, 260), "first power witness")
    need(euclid(8, 1) == (63, 16, 65), "normalized witness")
    need((63, (63**2-1)//2, (63**2+1)//2) == (63, 1984, 1985), "blue witness")
    for n in range(1, 21):
        for r in range(2, 10):
            need((n+1)**r-n**r > 1, "consecutive parameter obstruction")
    mod2_point = (1, 1, 0)
    need((mod2_point[0]**2+mod2_point[1]**2-mod2_point[2]**2) % 2 == 0,
         "original-model F2 point")
    need(tuple(x % 2 for x in phi(mod2_point, 2)) != mod2_point,
         "Phi2 differs from F2 coordinate Frobenius")
    need(tuple(x % 2 for x in phi((1, 0, 1), 2)) == (0, 0, 0),
         "original-model Phi2 base point")
    print("Odd-prime coefficient lifts PASS; characteristic-two original-model hostile retained")
    print("psi_3: (3,4,5) -> (63,16,65); blue recipe at 63 -> (63,1984,1985)")

    set_cases = 0
    for m in range(2, 5):
        for n in range(1, m):
            for r in range(1, 4):
                aa = list(product(range(m), repeat=r))
                bb = set(product(range(n), repeat=r))
                xs = [(0, t) for t in aa] + [(1, t) for t in sorted(bb)]
                def kappa(x):
                    side, t = x
                    return (1-side, t) if t in bb else x
                fixed = [x for x in xs if kappa(x) == x]
                all_orbits, nontrivial = set(), set()
                for x in xs:
                    need(kappa(kappa(x)) == x, "involution")
                    for y in xs:
                        pair, other = (x, y), (kappa(x), kappa(y))
                        rep = min(pair, other)
                        all_orbits.add(rep)
                        if pair != other:
                            nontrivial.add(rep)
                triple = (len(xs)*len(fixed), len(nontrivial), len(all_orbits))
                need(triple == euclid(m**r, n**r), "actual orbit cardinalities")
                need(len(fixed) == (m-n)*sum(m**(r-1-j)*n**j for j in range(r)),
                     "exceptional-state growth")
                set_cases += 1
    print("Nested-set power carriers: explicit involution orbits in", set_cases, "cases PASS")

    for q in (2, 3, 4, 5, 8, 9, 11):
        coeff = [Fraction(1)]
        for n in range(1, 11):
            coeff.append(sum(Fraction(q**m+1)*coeff[n-m]
                             for m in range(1, n+1))/n)
            need(coeff[n] == sum(q**j for j in range(n+1)), "zeta coefficient")
    print("Local zeta: fixed-count exponential recursion equals two-factor product through degree 10")
    print("PASS: exact carrier and failed subfamily transport; no RH implication established")


if __name__ == "__main__":
    main()
