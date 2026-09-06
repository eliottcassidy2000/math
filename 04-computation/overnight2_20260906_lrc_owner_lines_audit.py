#!/usr/bin/env python3
"""Independent owner-line proof audit: solve first coordinate; new line chart."""
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import gcd
import sys

CHECKS = 0


def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def cross(p, q):
    return (p[1]*q[2]-p[2]*q[1], p[2]*q[0]-p[0]*q[2], p[0]*q[1]-p[1]*q[0])


def dictionary(w):
    a, b, c = w
    xroof, yroof, zroof = (3*(sum(w)-v) for v in w)
    found = set()
    # Different raw enumeration from the producer: enumerate C2,C3,
    # solve C1, and keep literal strict integer roof comparisons.
    for y in range(-(yroof//14), yroof//14+1):
        if y % 3 == 0 or 14*abs(y) >= yroof:
            continue
        for z in range(-(zroof//14), zroof//14+1):
            if z % 3 == 0 or 14*abs(z) >= zroof:
                continue
            numerator = -b*y-c*z
            if numerator % a:
                continue
            x = numerator//a
            if x % 3 and 14*abs(x) < xroof:
                found.add((x,y,z))
    return found


def classify(w, points):
    colors = {1:[], 2:[]}
    for p in points:
        values = {(wi*pi) % 3 for wi, pi in zip(w, p)}
        need(len(values) == 1 and 0 not in values, "literal owner color")
        colors[values.pop()].append(p)
    A = sorted(colors[1])
    need({tuple(-x for x in p) for p in A} == set(colors[2]), "opposite owner classes")
    if len(A) < 2:
        return None
    p, q = A[0], A[-1]
    if cross(p, q) == (0,0,0):
        # An owner-collinear set with these endpoints would have rank one.
        return None
    delta = tuple(qi-pi for pi, qi in zip(p, q))
    if any(cross(delta, tuple(xi-pi for xi, pi in zip(x,p))) != (0,0,0) for x in A):
        return None
    divisor = gcd(gcd(abs(delta[0]), abs(delta[1])), abs(delta[2]))
    r = tuple(x//divisor for x in delta)
    pivot = next(j for j in range(3) if r[j])
    if r[pivot] < 0:
        r = tuple(-x for x in r)
    need(cross(p,r) in (w,tuple(-x for x in w)), "saturated transverse determinant")
    need(any(x % 3 == 0 for x in r), "invisible primitive line direction")
    coordinates = [Q(x[pivot]-p[pivot], r[pivot]) for x in A]
    need(sorted(coordinates) == [3*j for j in range(len(A))], "all same-color points have step 3r")
    N, c = len(points), w[2]
    need(Q(N) < Q(c,7)+2, "uniform strict diameter/chord bound")
    if max(map(abs,r)) <= 3:
        need(sorted(map(abs,r)) == [1,2,3], "complete small invisible coefficient shape")
        need(Q(N) < Q(2*c,21)+2, "small-direction affine chord bound")
    return r


def ratio_interval(r):
    # Intersect the relation with the closed ratio triangle by its edges.
    # This avoids the producer's one-variable inequality interval compiler.
    triangle = ((Q(0),Q(0)), (Q(0),Q(1)), (Q(1),Q(1)))
    evaluate = lambda p:r[0]*p[0]+r[1]*p[1]+r[2]
    points = set()
    for j in range(3):
        p, q = triangle[j], triangle[(j+1) % 3]
        fp, fq = evaluate(p), evaluate(q)
        if fp == 0:
            points.add(p)
        if fp*fq < 0:
            t = fp/(fp-fq)
            points.add(tuple(p[i]+t*(q[i]-p[i]) for i in range(2)))
    if len(points) != 2:
        return None
    p, q = sorted(points)
    alpha, beta = ((p[j]+q[j])/2 for j in range(2))
    if not 0 < alpha < beta < 1:
        return None
    return min(p[1],q[1]), max(p[1],q[1])


def interval_roofs(r, beta):
    alpha = Q(-r[1]*beta-r[2], r[0])
    w = (alpha,beta,Q(1))
    # A different affine section: first coordinate zero, instead of third.
    u = (Q(0),Q(-1,r[0]),beta/r[0])
    need(cross(u,r) == w, "independent normalized affine section")
    lower, upper = [], []
    for j in range(3):
        roof = Q(3,14)*(sum(w)-w[j])
        ends = sorted(((-roof-u[j])/r[j], (roof-u[j])/r[j]))
        lower.append(ends[0])
        upper.append(ends[1])
    return lower, upper


def chord_certificates():
    rows = []
    for magnitudes in permutations((1,2,3)):
        for signs in product((-1,1),repeat=2):
            r = (magnitudes[0],signs[0]*magnitudes[1],signs[1]*magnitudes[2])
            bounds = ratio_interval(r)
            if bounds is None:
                continue
            lo, hi = bounds
            l0, u0 = interval_roofs(r,lo)
            l1, u1 = interval_roofs(r,hi)
            options = []
            for i,j in product(range(3),repeat=2):
                v0, v1 = u0[i]-l0[j], u1[i]-l1[j]
                slope = (v1-v0)/(hi-lo)
                intercept = v0-slope*lo
                options.append((max(v0,v1),i+1,j+1,slope,intercept))
            best = min(options)
            need(best[0] <= Q(1,7), "entire ratio interval has a 1/7 chord certificate")
            rows.append((r,tuple(map(str,bounds)),tuple(map(str,best))))
    need(len(rows) == 7, "all feasible sign and order cases")
    return sorted(rows)


def main():
    sys.stdout.reconfigure(newline="\n")
    certs = chord_certificates()
    speeds = [n for n in range(1,52) if n % 2 and n % 3]
    universe, hits, sizes, maximizers, maximum = 0, 0, {}, [], Q(0)
    for w in combinations(speeds,3):
        if gcd(gcd(w[0],w[1]),w[2]) != 1:
            continue
        universe += 1
        points = dictionary(w)
        r = classify(w,points)
        if r is None:
            continue
        hits += 1
        sizes[len(points)] = sizes.get(len(points),0)+1
        ratio = Q(len(points),w[2])
        need(ratio < Q(2,11), ("complete finite head strictly count-safe",w))
        if ratio > maximum:
            maximum, maximizers = ratio, [w]
        elif ratio == maximum:
            maximizers.append(w)
    need((universe,hits,sizes) == (678,139,{4:137,6:2}), "complete independent finite-head classification")
    named = []
    for w in ((17,23,25),(19,23,29),(23,29,37),(41,47,49),(1,137,277),(1,499,1001)):
        points = dictionary(w)
        r = classify(w,points)
        named.append((w,len(points),r))
        need((r is None) == (w in ((19,23,29),(23,29,37))), "included and excluded dictionary controls")
    print("status=INDEPENDENT_AUDIT_PASS; finite base and symbolic affine chord certificates")
    print("finite_universe="+str(universe)+"; owner_line_rows="+str(hits)+"; N_counts="+str(sorted(sizes.items())))
    print("finite_max_N_over_c="+str(maximum)+"; all_maximizers="+str(maximizers))
    print("seven_chord_certificates=(r,beta_interval,(bound,upper_i,lower_j,slope,constant))")
    for row in certs:
        print(row)
    print("named_controls="+str(named))
    print("strict_survivor=N<2c/11; all_three_projections<6/77; full_rank_two_owner_line_scope_only")
    print("checks="+str(CHECKS))
    print("LRC14=OPEN; RESULT=PASS")


if __name__ == "__main__":
    main()
