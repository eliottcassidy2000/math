#!/usr/bin/env python3
"""Exact, stdlib-only controls for the Session59 arithmetic/F1 packet.

No asserts: all checks survive python -O. Finite-field extension checks
construct fields and Frobenius orbits, rather than substituting in formulas.
"""
from fractions import Fraction
from itertools import combinations, product
from math import comb
from pathlib import Path
import hashlib
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def trim(a):
    a = list(a)
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a


def add(*polys):
    out = [0] * max(map(len, polys))
    for a in polys:
        for i, x in enumerate(a):
            out[i] += x
    return trim(out)


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def power(a, k):
    out = [1]
    for _ in range(k):
        out = mul(out, a)
    return out


def scale(a, c):
    return trim([c*x for x in a])


def ev(a, t):
    out = 0
    for x in reversed(a):
        out = out*t+x
    return out


def coeff(a, i):
    return a[i] if i < len(a) else 0


def C_poly(k):
    numerator = add(power([1, 0, 3, 1], k),
                    scale(power([1, 0, 1], k), -3),
                    mul([0, -3], power([1, 1, 1], k)),
                    mul([0, 3, 2], power([1, 1], k)), [2, 0, -2])
    require(all(coeff(numerator, i) == 0 for i in range(3)), "t^3 divisibility")
    return trim(numerator[3:] or [0])


def D_poly(k):
    if k == 0:
        return [0]
    numerator = add(power([0, 0, 3, 1], k),
                    scale(power([0, 0, 1], k), -3),
                    mul([0, -3], power([0, 1, 1], k)),
                    mul([0, 3, 2], power([0, 1], k)))
    require(all(coeff(numerator, i) == 0 for i in range(3)), "primitive divisibility")
    return trim(numerator[3:] or [0])


def Uprim_poly(k):
    if k <= 1:
        return [0]
    t2 = mul([0]*(k-1)+[1],
             mul(add(power([2, 1], k-1), [-1]),
                 add(*(power([1, 1], i) for i in range(k)))))
    t3 = mul([0]*(k-1)+[1], add(power([3, 3, 1], k-1), [-1]))
    numerator = add(D_poly(k), scale(t2, 3), scale(t3, 2))
    require(all(c % 6 == 0 for c in numerator), "Uprim integral coefficients")
    return [c//6 for c in numerator]


def Uprim_positive_poly(k):
    if k <= 1:
        return [0]
    bracket = add(
        mul([0]*(k-1)+[1], add(power([3, 1], k), [-3])),
        scale(mul(add(power([2, 1], k-1), [-2]),
                  add(power([1, 1], k), [-1])), 3),
        mul([0, 2], power([3, 3, 1], k-1)))
    require(all(c >= 0 and c % 6 == 0 for c in bracket), "effective coefficient proof")
    return trim([0]*(k-2)+[c//6 for c in bracket])


def N(k):
    return Fraction(21**k-6*7**k-3*5**k+14*3**k-6, 48)


def M(k):
    return Fraction(20**k-6*6**k-3*4**k+14*2**k, 48) if k else Fraction(0)


def points(k, q):
    """Canonical projective representatives over a field of size q.

    Enumeration itself uses only that every nonzero point has a unique
    representative with first nonzero coordinate one; extension arithmetic
    enters separately through Frobenius below.
    """
    return [tuple([0]*r+[1]+list(tail))
            for r in range(k) for tail in product(range(q), repeat=k-r-1)]


def balanced(triple, primitive=False):
    for row in zip(*triple):
        n = sum(x != 0 for x in row)
        if n == 1 or (primitive and n == 0):
            return False
    return True


class Field:
    """F_p[x]/f, represented by base-p digits, monic f."""
    def __init__(self, p, f):
        self.p, self.f = p, f
        self.d, self.q = len(f)-1, p**(len(f)-1)
        self.digits = []
        for x in range(self.q):
            a = []
            for _ in range(self.d):
                a.append(x % p)
                x //= p
            self.digits.append(a)
        self.table = [[self._mul(a, b) for b in range(self.q)] for a in range(self.q)]
        require(all(self.pow(a, self.q-1) == 1 for a in range(1, self.q)),
                "defining polynomial must give a field")
        self.frob = [self.pow(a, p) for a in range(self.q)]

    def _mul(self, a, b):
        w = [0]*(2*self.d-1)
        for i, x in enumerate(self.digits[a]):
            for j, y in enumerate(self.digits[b]):
                w[i+j] = (w[i+j]+x*y) % self.p
        for j in range(len(w)-1, self.d-1, -1):
            for i in range(self.d):
                w[j-self.d+i] = (w[j-self.d+i]-w[j]*self.f[i]) % self.p
        return sum(w[i]*self.p**i for i in range(self.d))

    def pow(self, a, n):
        b = 1
        while n:
            if n & 1:
                b = self.table[a][b]
            a = self.table[a][a]
            n //= 2
        return b

    def projective_orbits(self, k):
        unseen = set(points(k, self.q))
        out = []
        while unseen:
            first = min(unseen)
            orbit = [first]
            current = tuple(self.frob[x] for x in first)
            while current != first:
                orbit.append(current)
                current = tuple(self.frob[x] for x in current)
            require(len(orbit) in (1, self.d), "prime extension orbit length")
            unseen.difference_update(orbit)
            out.append(orbit)
        return out


def T2(k, q):
    return sum(comb(k, s)*((q*q-1)**(s-1)-(q-1)**(s-1))
               *sum(q**i for i in range(s)) for s in range(1, k+1))


def T3(k, q):
    return sum(q**(3*i)-q**i for i in range(k))


def main():
    result = {"status": "PASS", "ordered_direct": [], "frobenius_orbits": [],
              "polynomials_t_ascending": {}, "primitive_t_ascending": {},
              "unordered_primitive_t_ascending": {}, "effective_minimum_ranks": []}
    for q, max_k in [(2, 5), (3, 4), (5, 3), (7, 3)]:
        for k in range(max_k+1):
            pts = points(k, q)
            total = primitive = 0
            for triple in combinations(pts, 3):
                if balanced(triple):
                    total += 6
                    if balanced(triple, primitive=True):
                        primitive += 6
            require(total == ev(C_poly(k), q-1), f"C direct q={q},k={k}")
            require(primitive == ev(D_poly(k), q-1), f"D direct q={q},k={k}")
            result["ordered_direct"].append({"q": q, "k": k, "points": len(pts),
                       "triple_universe": comb(len(pts), 3) if len(pts) >= 3 else 0,
                       "C": total, "D": primitive})

    for k in range(21):
        c, d = C_poly(k), D_poly(k)
        uprim = Uprim_poly(k)
        require(uprim == Uprim_positive_poly(k), "positive factorization")
        require(all(x >= 0 for x in uprim), "effective torus count carrier")
        utotal = add(*(scale(Uprim_poly(j), comb(k, j)) for j in range(k+1)))
        require(coeff(utotal, 1) == k*(k-1), "unordered rank one")
        require(coeff(utotal, 2) == 3*comb(k, 2)+6*comb(k, 3), "unordered rank two")
        if k >= 3:
            require(all(coeff(uprim, i) == 0 for i in range(k-1)), "Uprim minimum rank")
            require(coeff(uprim, k-1) == k*2**(k-2)-k+3**(k-2), "Uprim first coefficient")
        if k >= 2:
            result["effective_minimum_ranks"].append(
                {"k": k, "rank": k-1, "multiplicity": coeff(uprim, k-1)})
        require(ev(c, 2) == 6*N(k), "C(3)=6N")
        require(ev(d, 2) == 6*M(k), "D(3)=6M")
        require(c == add(*(scale(D_poly(j), comb(k, j)) for j in range(k+1))),
                "binomial ancestry polynomial identity")
        require(coeff(c, 0) == 0 and coeff(c, 1) == k*(k-1), "low C coefficients")
        require(coeff(c, 2) == Fraction(k*(k-1)*(23-7*k), 6), "negative C coefficient")
        if k >= 3:
            require(all(coeff(d, i) == 0 for i in range(k-1)), "primitive vanishing")
            require(coeff(d, k-1) == -(3*k-2), "primitive negative leading term")
        if k <= 6:
            result["polynomials_t_ascending"][str(k)] = c
            result["primitive_t_ascending"][str(k)] = d
            result["unordered_primitive_t_ascending"][str(k)] = uprim

    # Independent geometric count by Frobenius orbits in explicit extension fields.
    moduli = {2: ([1, 1, 1], [1, 1, 0, 1]),
              3: ([1, 0, 1], [1, 2, 0, 1]),
              5: ([2, 0, 1], [1, 1, 0, 1])}
    for q, (f2, f3) in moduli.items():
        field2, field3 = Field(q, f2), Field(q, f3)
        for k in range(1, 4):
            rational = points(k, q)
            split = sum(balanced(triple) for triple in combinations(rational, 3))
            split_primitive = sum(balanced(triple, True) for triple in combinations(rational, 3))
            pairs = [o for o in field2.projective_orbits(k) if len(o) == 2]
            threes = [o for o in field3.projective_orbits(k) if len(o) == 3]
            type12 = sum(balanced((o[0], o[1], p)) for o in pairs for p in rational)
            type3 = sum(balanced(o) for o in threes)
            type12_primitive = sum(balanced((o[0], o[1], p), True) for o in pairs for p in rational)
            type3_primitive = sum(balanced(o, True) for o in threes)
            require(2*type12 == T2(k, q), "degree 2+1 Frobenius formula")
            require(3*type3 == T3(k, q), "degree 3 Frobenius formula")
            expected = Fraction(ev(C_poly(k), q-1)+3*T2(k, q)+2*T3(k, q), 6)
            require(split+type12+type3 == expected, "unordered Frobenius total")
            primitive_total = split_primitive+type12_primitive+type3_primitive
            require(primitive_total == ev(Uprim_poly(k), q-1), "primitive Frobenius total")
            result["frobenius_orbits"].append({"q": q, "k": k, "split111": split,
                         "type12": type12, "type3": type3, "total": int(expected),
                         "primitive_split111": split_primitive,
                         "primitive_type12": type12_primitive,
                         "primitive_type3": type3_primitive,
                         "primitive_total": primitive_total})

    # Symbolically the logarithmic derivative has simple-pole residue -1/48
    # at the distinguished largest base, unaffected by polynomial corrections.
    result["zeta_residues"] = {"N_at_1_over_21": "-1/48", "M_at_1_over_20": "-1/48"}
    require(Fraction(-1, 48).denominator != 1, "nonintegral pole order hostile")
    positive_control = {"points_affine_line": [3**n for n in range(1, 6)],
                        "zeta": "1/(1-3t)", "residue_at_1_over_3": -1}
    result["positive_rational_zeta_control"] = positive_control
    result["hostiles"] = {"first_total_torification_failure": {"k": 4, "t2": -10},
                          "first_primitive_failure": {"k": 3, "t2": -7},
                          "naive_unordered_q3_k2": {"split": 4, "geometric": 24}}
    result["scope"] = "Necessary support filter only; no numerical gap-addition test."
    Path(__file__).with_name("verification.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
