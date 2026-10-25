#!/usr/bin/env python3
"""Exact supporting controls for the shape/spectral bridge; standard library only.

The unrestricted claims are proved in shapes.md. This finite certificate checks
named positive and hostile controls; it is not an RH or Jacobian computation.
"""
from fractions import Fraction as Q
import json


checks = 0


def require(predicate, label):
    global checks
    checks += 1
    if not predicate:
        raise RuntimeError(label)


def transpose(a):
    return tuple(zip(*a))


def mul(a, b):
    return tuple(tuple(sum(x*y for x, y in zip(row, col))
                       for col in zip(*b)) for row in a)


def scale(c, a):
    return tuple(tuple(c*x for x in row) for row in a)


def power(a, n):
    r = tuple(tuple(int(i == j) for j in range(len(a)))
              for i in range(len(a)))
    for _ in range(n):
        r = mul(r, a)
    return r


def det2(a):
    return a[0][0]*a[1][1]-a[0][1]*a[1][0]


def act(a, v):
    return tuple(sum(x*y for x, y in zip(row, v)) for row in a)


def clean(p):
    return {e: c for e, c in p.items() if c}


def add(p, q, factor=1):
    out = p.copy()
    for e, c in q.items():
        out[e] = out.get(e, 0)+factor*c
    return clean(out)


def product(p, q):
    out = {}
    for (a,b), c in p.items():
        for (d,e), f in q.items():
            exponent = (a+d,b+e)
            out[exponent] = out.get(exponent, 0)+c*f
    return clean(out)


def derivative(p, axis):
    out = {}
    for exponent, coefficient in p.items():
        if exponent[axis]:
            target = list(exponent)
            target[axis] -= 1
            out[tuple(target)] = coefficient*exponent[axis]
    return clean(out)


def jac(p, q):
    return add(product(derivative(p,0), derivative(q,1)),
               product(derivative(p,1), derivative(q,0)), -1)


I = ((1,0),(0,1))
A3 = ((0,-1),(1,-1))
A4 = ((0,-1),(1,0))
triangle = {(1,0),(0,1),(-1,-1)}
square = {(x,y) for x in (-1,1) for y in (-1,1)}
require({act(A3,v) for v in triangle} == triangle, "triangle permutation")
require({act(A4,v) for v in square} == square, "square permutation")
for a, order in ((A3,3),(A4,4)):
    require(power(a,order) == I, "finite-order positive control")
    for k in range(1,order):
        require(power(a,k) != I, "minimal order")

triangle_cov = tuple(tuple(sum(Q(v[i]*v[j],12) for v in triangle)
                           for j in range(2)) for i in range(2))
require(triangle_cov == ((Q(1,6),Q(1,12)),(Q(1,12),Q(1,6))),
        "uniform triangle covariance formula")
require(mul(mul(A3,triangle_cov),transpose(A3)) == triangle_cov,
        "covariance equivariance")
require(det2(triangle_cov) > 0 and triangle_cov[0][0] > 0,
        "positive covariance")
require(mul(mul(transpose(A3),((2,-1),(-1,2))),A3) == ((2,-1),(-1,2)),
        "triangle positive form")

# Every rational trace below is nonintegral. If any matrix had finite order,
# its eigenvalues and trace would be algebraic integers, a contradiction.
rotations = []
for t in (Q(1,2),Q(2,3),Q(2),Q(3)):
    c, s = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    r = ((c,-s),(s,c))
    require(mul(transpose(r),r) == I, "rational rotation orthogonality")
    require((2*c).denominator > 1, "infinite-order trace certificate")
    rotations.append({"parameter":str(t), "trace":str(2*c)})

# Rescaled purity does not require normalized eigenvalues to be roots of unity.
T, G = ((0,-2),(1,1)), ((2,1),(1,4))
require(mul(mul(transpose(T),G),T) == scale(2,G), "weight-one similitude")
require(G[0][0] > 0 and det2(G) > 0, "positive similitude form")
require(det2(T) == 2 and T[0][0]+T[1][1] == 1,
        "characteristic polynomial z^2-z+2")

H, indefinite = ((2,1),(1,1)), ((2,-1),(-1,-2))
require(mul(mul(transpose(H),indefinite),H) == indefinite,
        "indefinite invariant form")
require(det2(indefinite) == -5, "indefinite signature")
require(det2(H) == 1 and H[0][0]+H[1][1] == 3,
        "hyperbolic reciprocal polynomial z^2-3z+1")
J = ((1,1),(0,1))
for k in (0,1,2,5,20):
    require(power(J,k) == ((1,k),(0,1)), "Jordan growth hostile")

# Same full supports, not just the same Newton polygons, can hide Keller failure.
P = {(1,0):1,(0,2):1}
target_support = {(0,1),(2,0),(1,2),(0,4)}
jacobian_rows = []
for t in (-1,1,2,3):
    Qt = {(0,1):1,(2,0):1,(1,2):t,(0,4):1}
    require(set(Qt) == target_support, "unchanged full support")
    expected = clean({(0,0):1,(1,1):2*(t-2),(0,3):-2*(t-2)})
    actual = jac(P,Qt)
    require(actual == expected, "exact Keller coefficient obstruction")
    require((actual == {(0,0):1}) == (t == 2), "unit Jacobian iff t=2")
    jacobian_rows.append({"t":t, "coefficients":{
        f"x^{i}y^{j}":c for (i,j),c in sorted(actual.items())}})

# An independent determinant expansion at nine rational points per parameter.
for t in (-1,1,2,3):
    for x in (Q(-1),Q(0),Q(2,3)):
        for y in (Q(-1,2),Q(0),Q(2)):
            by_derivatives = (1+2*t*x*y+4*y**3)-2*y*(2*x+t*y*y)
            require(by_derivatives == 1+2*(t-2)*y*(x-y*y),
                    "independent pointwise determinant path")

print(json.dumps({
    "status":"FINITE-EXACT",
    "checks":checks,
    "universe":"named 2x2 matrix controls; t=-1,1,2,3 polynomial pairs; 36 rational point evaluations",
    "all_order_claims":"proved separately in shapes.md; no extrapolation from this finite run",
    "rational_infinite_order_rotations":rotations,
    "jacobian_support_hostiles":jacobian_rows,
    "passed":True
}, indent=2, sort_keys=True))
