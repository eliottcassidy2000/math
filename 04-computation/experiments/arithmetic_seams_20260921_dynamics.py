"""Exact controls for arithmetic seams: three-cycles, sign lifts and trace quotients.

No third-party dependencies. Universe and non-consequences are frozen in JSON.
Run normally and with python -O; require() keeps every check active.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import gcd
from pathlib import Path
import hashlib
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


I = (F(1), F(0), F(0), F(1))


def mm(a, b):
    return (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3])


def mv(a, v):
    return (a[0]*v[0]+a[1]*v[1], a[2]*v[0]+a[3]*v[1])


def det(a):
    return a[0]*a[3]-a[1]*a[2]


def inverse(a):
    d = det(a)
    require(d != 0, "singular matrix")
    return tuple(t/d for t in (a[3], -a[1], -a[2], a[0]))


def mpow(a, n):
    out = I
    for _ in range(n):
        out = mm(out, a)
    return out


def slope(v):
    require(v[1] != 0, "unexpected infinite cycle point")
    return v[0]/v[1]


def closure(b, r):
    group = {mm(mpow(b, k), mpow(r, e)) for k in range(6) for e in (0, 1)}
    require(len(group) == 12, "dihedral order")
    require(all(mm(x, y) in group for x in group for y in group), "closure")
    require(mm(mm(r, b), r) == inverse(b), "reflection reverses rotation")
    central = {I, mpow(b, 3)}
    s3 = {mm(mpow(b, 2*k), mpow(r, e)) for k in range(3) for e in (0, 1)}
    require(len(s3) == 6 and central & s3 == {I}, "direct-product factors")
    require({mm(x, y) for x in central for y in s3} == group, "C2 times S3")
    return group


def marked_cycle(t):
    d = 2*t*(t+1)
    require(d != 0, "excluded chart point")
    p = ((t**3+2*t*t+t+1)/d, (t**3-t-1)/d,
         -(t**3+2*t*t+3*t+1)/d)
    c = p[1]-p[0]**2
    require(len(set(p)) == 3, "exact three-cycle")
    require(all(p[j]**2+c == p[(j+1) % 3] for j in range(3)), "chart dynamics")
    return p, c


# Dense integer polynomial arithmetic, constant coefficient first.
def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)


def padd(a, b):
    out = [0]*max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    return trim(out)


def pscale(a, s):
    return trim([s*x for x in a])


def pmul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def pdiv(a, b):
    r = list(trim(a))
    b = trim(b)
    require(b != (0,), "polynomial zero divisor")
    q = [0]*max(1, len(r)-len(b)+1)
    while len(r) >= len(b) and r != [0]:
        shift = len(r)-len(b)
        require(r[-1] % b[-1] == 0, "integral polynomial quotient")
        value = r[-1] // b[-1]
        q[shift] = value
        for j in range(len(b)):
            r[shift+j] -= value*b[j]
        r = list(trim(r))
    require(r == [0], "exact polynomial division")
    return trim(q)


def pcompose(a, b):
    out = (0,)
    for x in reversed(a):
        out = padd(pmul(out, b), (x,))
    return out


@lru_cache(None)
def cyclotomic(n):
    p = [-1]+[0]*(n-1)+[1]
    for d in range(1, n):
        if n % d == 0:
            p = pdiv(p, cyclotomic(d))
    return trim(p)


def trace_cyclotomic(n):
    phi = cyclotomic(n)
    d = (len(phi)-1)//2
    require(n > 2 and len(phi) == 2*d+1 and phi == phi[::-1], "reciprocal Phi")
    traces = [(2,), (0, 1)]
    for k in range(1, d):
        traces.append(padd(pmul((0, 1), traces[-1]), pscale(traces[-2], -1)))
    out = (phi[d],)
    for k in range(1, d+1):
        out = padd(out, pscale(traces[k], phi[d+k]))
    return out


def multiplicative_order_two(n):
    require(n > 1 and n % 2 == 1, "odd conductor")
    value, k = 2 % n, 1
    while value != 1:
        value = 2*value % n
        k += 1
        require(k <= n, "order termination")
    return k


def trace_period(n):
    value, k = 2 % n, 1
    while value not in (1, n-1):
        value = 2*value % n
        k += 1
        require(k <= n, "trace period termination")
    return k, 1 if value == 1 else -1


def main():
    b = (F(1, 4), F(-13, 4), F(1, 4), F(3, 4))
    r = (F(-1), F(-2), F(0), F(1))
    require(det(b) == 1 and mpow(b, 3) == tuple(-x for x in I), "central sign")
    group = closure(b, r)
    v = (F(-7), F(1))
    hexagon = [mv(mpow(b, k), v) for k in range(6)]
    norm = lambda w: w[0]**2+2*w[0]*w[1]+13*w[1]**2
    integer_level = {(F(x), F(y)) for x in range(-9, 10) for y in range(-2, 3)
                     if x*x+2*x*y+13*y*y == 48}
    require(set(hexagon) == integer_level, "complete integral level Q=48")
    require(all(norm(w) == 48 for w in hexagon), "invariant norm")
    require(all({mv(g, w) for w in hexagon} == set(hexagon) for g in group), "hexagon action")
    require([slope(w) for w in hexagon] == [F(-7), F(5), F(-1)]*2, "projective quotient")
    # Full linear stabilizer: an invertible map is determined by images of a basis.
    frame = (hexagon[0][0], hexagon[1][0], hexagon[0][1], hexagon[1][1])
    full_stabilizer = set()
    for u in hexagon:
        for w in hexagon:
            candidate_frame = (u[0], w[0], u[1], w[1])
            if det(candidate_frame):
                g = mm(candidate_frame, inverse(frame))
                if {mv(g, h) for h in hexagon} == set(hexagon):
                    full_stabilizer.add(g)
    require(full_stabilizer == group, "full rational linear hexagon stabilizer")
    c = (F(1), F(1), F(0), F(2))
    require(mm(mm(c, b), inverse(c)) == (F(1, 2), F(-3, 2), F(1, 2), F(1, 2)),
            "Eisenstein rotation coordinates")
    require(mm(mm(c, r), inverse(c)) == (F(-1), F(0), F(0), F(1)),
            "negative conjugation coordinates")

    generic_count = 0
    parameter_set = {F(a, d) for d in range(1, 9) for a in range(-12, 13)}-{F(0), F(-1)}
    for t in sorted(parameter_set):
        points, parameter = marked_cycle(t)
        sigma = sum(points)
        lift = (sigma+1, -(sigma*sigma+sigma+1), F(1), -sigma)
        start = (points[0], F(1))
        nxt = mv(lift, start)
        frame = (start[0], nxt[0], start[1], nxt[1])
        model_b = (F(0), F(-1), F(1), F(1))
        model_j = (F(1), F(1), F(0), F(-1))
        reflection = mm(mm(frame, model_j), inverse(frame))
        require(mm(mm(inverse(frame), lift), frame) == model_b, "universal lift basis")
        lifts = [mv(mpow(lift, k), start) for k in range(6)]
        require([slope(w) for w in lifts] == list(points)*2, "generic cycle projection")
        require(len(closure(lift, reflection)) == 12, "generic D6 action")
        affine_reverse = (F(-1), F(-1, 2), F(0), F(1))
        invariant = {slope(mv(affine_reverse, (x, F(1)))) for x in points} == set(points)
        require(invariant == (parameter == F(-29, 16)), "AP affine-reflection rigidity")
        generic_count += 1

    # The all-Q graph proof reduces to exactly these eight quarter-integers.
    graph = {a: (a*a-29)//4 for a in (-7, -5, -3, -1, 1, 3, 5, 7)}
    require(graph == {-7: 5, -5: -1, -3: -5, -1: -7, 1: -7, 3: -5, 5: -1, 7: 5},
            "complete bounded rational graph")
    require(F(8)*F(-7, 4)*F(5, 4)*F(-1, 4) == F(35, 8), "quadratic cycle multiplier")
    critical, x = [], F(0)
    for _ in range(4):
        x = x*x-F(7, 4)
        critical.append(x)
    require(critical[:3] == [F(-7, 4), F(21, 16), F(-7, 256)], "critical third exception")

    # Independent finite primitive-support check by stripping all old prime powers.
    old_product, exceptional = 1, []
    for n in range(1, 41):
        value = 2**n-1
        fresh = value
        while gcd(fresh, old_product) > 1:
            fresh //= gcd(fresh, old_product)
        if fresh == 1:
            exceptional.append(n)
        old_product *= value
    require(exceptional == [1, 6], "Mersenne finite exceptional indices")
    require(multiplicative_order_two(3) == 2 and multiplicative_order_two(7) == 3,
            "old prime orders")
    require(multiplicative_order_two(9) == 6 and multiplicative_order_two(63) == 6,
            "new order without new prime support")

    trace_rows = []
    for n in range(3, 130, 2):
        period, sign = trace_period(n)
        order = multiplicative_order_two(n)
        require(order == period*(1 if sign == 1 else 2), "trace quotient period law")
        if period in (3, 6):
            degree = sum(gcd(a, n) == 1 for a in range(1, n))//2
            trace_rows.append({"conductor": n, "trace_period": period,
                               "squaring_period": order, "multiplier": sign*2**period,
                               "trace_points": degree, "trace_cycles": degree//period})
    require([v["conductor"] for v in trace_rows if v["trace_period"] == 3] == [7, 9],
            "period-three conductors")
    require([v["conductor"] for v in trace_rows if v["trace_period"] == 6] == [13, 21, 63, 65],
            "period-six conductors")
    xpoly, cheb = (0, 1), (-2, 0, 1)
    iterates = [xpoly]
    for _ in range(6):
        iterates.append(pcompose(cheb, iterates[-1]))
    difference = lambda k: padd(iterates[k], (0, -1))
    dyn3 = pdiv(difference(3), difference(1))
    dyn6 = pdiv(pmul(difference(6), difference(1)), pmul(difference(3), difference(2)))
    require(dyn3 == pmul(trace_cyclotomic(7), trace_cyclotomic(9)), "third trace factorization")
    factor6 = (1,)
    for n in (13, 21, 63, 65):
        factor6 = pmul(factor6, trace_cyclotomic(n))
    require(dyn6 == factor6 and len(dyn6)-1 == 54, "complete sixth trace factorization")

    output = {
        "status": "PASS", "new_convergence_or_sphere_theorem": False,
        "universe": {"generic_rational_chart": "distinct t=a/d, |a|<=12, 1<=d<=8; exclude0,-1",
                     "generic_chart_count": generic_count, "mersenne_indices": [1, 40],
                     "odd_conductors": [3, 129], "cyclotomic_factorization": "exact polynomial identities",
                     "Q48_integer_bounds": {"X": [-9, 9], "Y": [-2, 2]}},
        "hexagon": [[str(x) for x in w] for w in hexagon],
        "linear_stabilizer_order": len(full_stabilizer), "projective_group_order": 6,
        "generic_lift_group": "C2 x S3; reflections reverse arrows",
        "rational_graph_y_equals_4x": graph,
        "critical_orbit_minus_seven_quarters": [str(x) for x in critical],
        "finite_mersenne_exceptions": exceptional, "trace_conductors": trace_rows,
        "trace_cyclotomic_polynomials_constant_first": {str(n): trace_cyclotomic(n) for n in (7, 9, 13, 21, 63, 65)},
        "dynatomic_six_degree": 54, "chebyshev_exact_six_cycles": 9,
        "anthropic_exact_S2_times_S3_reference": "UNRESOLVED; no substitution",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    Path(__file__).with_suffix(".json").write_text(json.dumps(output, indent=2, sort_keys=True)+"\n",
                                                encoding="utf-8", newline="\n")
    print(f"PASS: {generic_count} rational marked cycles, full12-element lift, exact trace factors and controls")


if __name__ == "__main__":
    main()
