"""Exact, bounded fruit-curve audit; Python standard library only.

Run from any directory. Output is the sibling JSON. No assert statements,
floating-point tests, rank computation, or global-minimality claim.
"""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import permutations
import json
from math import gcd, lcm
from pathlib import Path


CHECKS = []


def require(condition, label):
    if not condition:
        raise RuntimeError(label)
    CHECKS.append(label)


class Poly:
    """Sparse Z[a,b,c,N], used to compare every coefficient exactly."""
    def __init__(self, terms=None):
        self.terms = {k: v for k, v in (terms or {}).items() if v}

    @staticmethod
    def coerce(value):
        if isinstance(value, Poly):
            return value
        return Poly({(0, 0, 0, 0): value})

    @staticmethod
    def variable(i):
        powers = [0] * 4
        powers[i] = 1
        return Poly({tuple(powers): 1})

    def __add__(self, other):
        terms = dict(self.terms)
        for k, v in Poly.coerce(other).terms.items():
            terms[k] = terms.get(k, 0) + v
        return Poly(terms)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + -Poly.coerce(other)

    def __rsub__(self, other):
        return Poly.coerce(other) + -self

    def __mul__(self, other):
        terms = {}
        for p, v in self.terms.items():
            for q, w in Poly.coerce(other).terms.items():
                key = tuple(i + j for i, j in zip(p, q))
                terms[key] = terms.get(key, 0) + v * w
        return Poly(terms)

    __rmul__ = __mul__

    def __pow__(self, n):
        result = Poly.coerce(1)
        for _ in range(n):
            result = result * self
        return result

    def __eq__(self, other):
        return self.terms == Poly.coerce(other).terms


def cubic(v, n=4):
    a, b, c = v
    s = a + b + c
    return s**3 - (n + 2) * s * (a*b + a*c + b*c) + (n + 3)*a*b*c


def pair_product(v):
    a, b, c = v
    return (a+b)*(a+c)*(b+c)


def involution(v):
    a, b, c = v
    e = a*b + a*c + b*c
    return (-a*a+b*b+c*c+e, a*a-b*b+c*c+e, a*a+b*b-c*c+e)


def forward(v, n=4):
    a, b, c = v
    return (-4*(n+3)*(a+b+2*c), 4*(n+3)*(2*n+5)*(a-b),
            (n+2)*(a+b)-c)


def inverse(v, n=4):
    x, y, z = v
    return (8*(n+3)*z-x+y, 8*(n+3)*z-x-y,
            -8*(n+3)*z-2*(n+2)*x)


def symbolic_checks():
    a, b, c, n = (Poly.variable(i) for i in range(4))
    v = (a, b, c)
    s = a+b+c
    numerator = a*(a+b)*(a+c)+b*(a+b)*(b+c)+c*(a+c)*(b+c)
    require(numerator-n*pair_product(v) == cubic(v, n),
            "all-N cleared fruit equation coefficient identity")
    x, y, z = forward(v, n)
    residual = y*y*z-x**3-(4*n*n+12*n-3)*x*x*z-32*(n+3)*x*z*z
    require(residual == 64*(n+3)**2*(2*n+5)**2*cubic(v, n),
            "all-N projective curve identity")
    for i, item in enumerate(inverse((x, y, z), n)):
        require(item == 8*(n+3)*(2*n+5)*v[i], f"all-N inverse matrix coordinate {i}")
    j = involution(v)
    require(cubic(j, n) == 8*pair_product(v)*cubic(v, n),
            "all-N reciprocal involution preserves cubic")
    for i, item in enumerate(involution(j)):
        require(item == 8*pair_product(v)*v[i], f"reciprocal involution square coordinate {i}")
    x, y, z = forward(v)
    u = forward(j)
    w = (224*x*z, -224*y*z, x*x)
    require(u[0]*w[1]-u[1]*w[0] == 0, "translation-U polynomial minor XY")
    require(u[0]*w[2]-u[2]*w[0] == -81536*(a+b+2*c)*cubic(v),
            "translation-U polynomial minor XZ modulo cubic")
    require(u[1]*w[2]-u[2]*w[1] == -1059968*(a-b)*cubic(v),
            "translation-U polynomial minor YZ modulo cubic")


def normalize(v):
    v = tuple(Q(z) for z in v)
    d = lcm(*(z.denominator for z in v))
    integers = tuple(int(z*d) for z in v)
    g = gcd(gcd(integers[0], integers[1]), integers[2])
    if not g:
        raise RuntimeError("zero projective vector")
    integers = tuple(z//g for z in integers)
    if next(z for z in integers if z) < 0:
        integers = tuple(-z for z in integers)
    return integers


def positive(v):
    return all(z > 0 for z in v) or all(z < 0 for z in v)


def fruit(v):
    a, b, c = v
    return Q(a, b+c)+Q(b, a+c)+Q(c, a+b)


def add(p, q):
    if p is None:
        return q
    if q is None:
        return p
    x, y = p
    u, v = q
    if x == u and y == -v:
        return None
    m = (v-y)/(u-x) if x != u else (3*x*x+218*x+224)/(2*y)
    w = m*m-109-x-u
    return w, m*(x-w)-y


def neg(p):
    return None if p is None else (p[0], -p[1])


def from_point(p):
    return normalize(inverse((0, 1, 0) if p is None else (*p, 1)))


def to_point(v):
    x, y, z = forward(v)
    return None if z == 0 else (Q(x, z), Q(y, z))


def on_curve(p):
    return p is None or p[1]**2 == p[0]**3+109*p[0]**2+224*p[0]


def point_json(p):
    return "O" if p is None else [str(z) for z in p]


def counts(p, a, b, c):
    brute = 1+sum(y*y % p == (x**3+a*x*x+b*x+c) % p
                  for x in range(p) for y in range(p))
    legendre = 1
    for x in range(p):
        r = (x**3+a*x*x+b*x+c) % p
        legendre += 1 if r == 0 else 2 if pow(r, (p-1)//2, p) == 1 else 0
    require(brute == legendre, f"independent F{p} counts for ({a},{b},{c})")
    return brute


def main():
    symbolic_checks()
    a = 154476802108746166441951315019919837485664325669565431700026634898253202035277999
    b = 368751317941299998271978115652254748254929799689719709962831374716372246340555790
    c = 43736126779286972578612526023713901528165375581616136186214379933784234677720360
    comma_c = "43,736,126,779,286,972,578,612,526,023,713,901,528,165,375,581,616,136,186,214,379,933,784,234,677,720,360"
    require(c == int(comma_c.replace(",", "")), "literal comma transcription")
    require(fruit((a,b,c)) != 4 and cubic((a,b,c)) != 0, "literal pasted triple is hostile")
    repaired = (a, b//10, c//10)
    require(b % 10 == c % 10 == 0 and fruit(repaired) == 4,
            "explicit decimal repair satisfies fruit equation")
    require(gcd(gcd(*repaired[:2]), repaired[2]) == 1, "repaired triple primitive")
    g = (Q(-4), Q(28))
    t = (Q(56), Q(728))
    u = (Q(0), Q(0))
    torsion = [None]
    for _ in range(5):
        torsion.append(add(torsion[-1], t))
    require(len(set(torsion)) == 6 and add(torsion[-1], t) is None,
            "explicit T has exact order six")
    require(torsion[2] == (Q(4), Q(52)) and torsion[3] == u,
            "torsion coordinate identification")
    for j, p in enumerate(torsion):
        require(on_curve(p) and to_point(from_point(p)) == p, f"torsion map roundtrip {j}")
        require(pair_product(from_point(p)) == 0, f"torsion fruit denominators excluded {j}")
    records = []
    p = None
    ninth = None
    positive_indices = []
    for n in range(1, 13):
        p = add(p, g)
        v = from_point(p)
        require(on_curve(p) and to_point(v) == p, f"multiple {n} roundtrip")
        require(fruit(v) == 4, f"multiple {n} exact fruit identity")
        x = p[0]
        gate = x < Q(-14, 3) and x*x+112*x+784 > 0
        require(gate == positive(v), f"multiple {n} exact positivity gate")
        if positive(v):
            positive_indices.append(n)
        if n == 9:
            ninth = p
            require(v == repaired, "ninth multiple equals repaired triple exactly")
        records.append({"n": n, "positive": positive(v),
                        "coordinate_digits": [len(str(abs(z))) for z in v],
                        "torsion_translates_positive": []})
        for j, tj in enumerate(torsion):
            q = add(p, tj)
            w = from_point(q)
            require(on_curve(q) and to_point(w) == q, f"n={n},t={j} roundtrip")
            require(fruit(w) == 4, f"n={n},t={j} fruit")
            require(to_point((w[1], w[2], w[0])) == add(q, torsion[2]),
                    f"n={n},t={j} cyclic translation")
            require(to_point((w[1], w[0], w[2])) == neg(q),
                    f"n={n},t={j} transposition negation")
            require(to_point(involution(w)) == add(q, u),
                    f"n={n},t={j} central involution translation")
            require(normalize(involution(involution(w))) == normalize(w),
                    f"n={n},t={j} involution roundtrip")
            if positive(w):
                records[-1]["torsion_translates_positive"].append(j)
                w = tuple(abs(z) for z in w)
                iw = involution(w)
                require(sum(z < 0 for z in iw) == 1 and sum(z > 0 for z in iw) == 2,
                        f"n={n},t={j} positive input becomes one-negative")
    require(positive_indices[0] == 9, "first positive in specified first twelve multiples is ninth")
    v = from_point(g)
    orbit = {normalize(w) for w in permutations(v)}
    orbit |= {normalize(involution(w)) for w in permutations(v)}
    require(len(orbit) == 12, "C2 times S3 action is faithful on explicit orbit")
    fibonacci = [0, 1]
    for _ in range(2, 33):
        fibonacci.append(fibonacci[-1]+fibonacci[-2])
    sharpness = []
    for n in range(3, 32, 2):
        q, p = fibonacci[n], fibonacci[n+1]
        v = (2*p*q, 1, 2*q*q-1)
        jv = involution(v)
        value = fruit(v)
        require(q*q+p*q-p*p == 1, f"Fibonacci n={n} Cassini identity")
        require(v[0] > v[2] > v[1] > 0 and gcd(gcd(*v[:2]),v[2]) == 1,
                f"Fibonacci n={n} primitive positive ordered triple")
        require(jv[0] == 2*q*q+1 and all(z > 0 for z in jv),
                f"Fibonacci n={n} positive reciprocal image")
        require(value > 0 and value*value < 5,
                f"Fibonacci n={n} exact strict sqrt-five bound")
        sharpness.append({"n":n,"p":p,"q":q,"triple":list(v),
                          "fruit_sum":str(value),"J_a":jv[0]})
    finite_counts = {str(p): {"fruit": counts(p,109,224,0), "other": counts(p,4,0,4)}
                     for p in (3,11,17,19)}
    require(finite_counts["11"] == {"fruit":12,"other":16}, "isogeny hostile at prime eleven")
    require(finite_counts["17"]["fruit"] == 18, "fruit torsion bound second prime")
    source = Path(__file__)
    result = {
        "status": "PASS", "source_sha256_lf": sha256(source.read_bytes().replace(b"\r\n",b"\n")).hexdigest(),
        "scope": "Polynomial coefficient identities over Z[a,b,c,N]; exact rational finite controls only. No rank or global minimum computed.",
        "universe": {"multiples_of_G": [1,12], "torsion_translates_per_multiple":6,
                     "faithful_group_orbit_size":12, "good_primes":[3,11,17,19],
                     "fibonacci_sharpness_odd_indices":list(range(3,32,2)),
                     "finite_field_paths":"all (x,y) pairs and Legendre-symbol count"},
        "literal": {"triple":list(map(str,(a,b,c))),"fruit_value":str(fruit((a,b,c))),
                    "cubic_residual":str(cubic((a,b,c))),"equals_four":False},
        "decimal_repair": {"operation":"leave aa unchanged; divide bb and cc by ten", "triple":list(map(str,repaired)),
                           "equals_four":True,"point":point_json(ninth)},
        "ninth_plus_U":list(map(str,from_point(add(ninth,u)))),
        "torsion_points":list(map(point_json,torsion)),"multiples":records,"finite_field_counts":finite_counts,
        "fibonacci_sharpness":sharpness,
        "checks_passed":len(CHECKS),"checks":CHECKS
    }
    target = source.with_suffix(".json")
    target.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":"PASS","checks_passed":len(CHECKS),"output":target.name,
                      "first_positive_multiple_in_universe":positive_indices[0]}))


if __name__ == "__main__":
    main()
