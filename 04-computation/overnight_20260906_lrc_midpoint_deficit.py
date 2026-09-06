#!/usr/bin/env python3
"""Exact controls for midpoint curvature payment and its affine kernel.

Named controls and a declared family sample, not a new height census.
The all-height assertions are proved in the companion report.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from hashlib import sha256
import json
import sys

CHECKS = 0


def need(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def pos(x):
    return max(Q(0), x)


def carriers(w):
    a, b, c = w
    bx, by, bz = ((3*(sum(w)-x)-1)//14 for x in w)
    out = set()
    for x in range(-bx, bx+1):
        if not x % 3:
            continue
        for y in range(-by, by+1):
            if not y % 3:
                continue
            znum = -a*x-b*y
            if znum % c:
                continue
            z = znum//c
            if abs(z) <= bz and z % 3:
                out.add((x, y, z))
    return out


def params(w, i):
    j, k = [j for j in range(3) if j != i]
    alpha = Q(7*w[2], 3*w[j]*w[k])
    threshold = Q(3*(w[j]+w[k]), 14)-Q(3*w[j]*w[k], 7*w[2])
    need(threshold > 0, ("positive kink threshold", w, i))
    return alpha, threshold


def deficit(w, p, i):
    j, k = [j for j in range(3) if j != i]
    return pos(1-Q(w[2]*(3*(w[j]+w[k])-14*abs(p[i])), 6*w[j]*w[k]))


def bucket(w, p):
    colors = tuple(a*b % 3 for a, b in zip(w, p))
    need(len(set(colors)) == 1 and colors[0] != 0, ("owner color", w, p))
    return colors[0], tuple(x % 2 for x in p)


def payment(w, points):
    groups = {}
    for p in points:
        groups.setdefault(bucket(w, p), []).append(p)
    Ds, Bs, Rs = [], [], []
    all_curvatures = []
    for i in range(3):
        alpha, t = params(w, i)
        total = sum((deficit(w, p, i) for p in points), Q(0))
        bound = residual = Q(0)
        curvatures = []
        for group in groups.values():
            for p, q in combinations(group, 2):
                mid = tuple((a+b)//2 for a, b in zip(p, q))
                step = Q(q[i]-p[i], 2)
                need(mid in points, ("live AP midpoint", w, p, q, mid))
                curve = deficit(w, p, i)+deficit(w, q, i)-2*deficit(w, mid, i)
                hinge = alpha*(pos(abs(step)-abs(mid[i]-t))+pos(abs(step)-abs(mid[i]+t)))
                need(curve == hinge and curve >= 0, ("exact hinge curvature", w, i, p, q))
                curvatures.append(curve)
            ordered = sorted(group, key=lambda p:(p[i], p))
            for j in range(len(ordered)//2):
                p, q = ordered[j], ordered[-j-1]
                mid = tuple((a+b)//2 for a, b in zip(p, q))
                bound += deficit(w, p, i)+deficit(w, q, i)-2*deficit(w, mid, i)
                residual += 2*deficit(w, mid, i)
            if len(ordered) % 2:
                residual += deficit(w, ordered[len(ordered)//2], i)
        need(bound+residual == total and 0 <= bound <= total,
             ("payment plus affine remainder", w, i, total, bound, residual))
        Ds.append(total)
        Bs.append(bound)
        Rs.append(residual)
        all_curvatures.append(curvatures)
    return tuple(Ds), tuple(Bs), tuple(Rs), all_curvatures


def family(m):
    lo, hi = Q(5*m-9, 42), Q(3*m-1, 14)
    ks = tuple(k for k in range(int(lo)+1, int(hi)+1) if k < hi and k % 3 == 2)
    positive = {(m-3*k, -2*k-1, k) for k in ks}
    return ks, positive | {tuple(-x for x in p) for p in positive}


def matching_test(xs, t):
    f = lambda x:pos(abs(x)-t)
    weights = {(i, j):f(xs[i])+f(xs[j])-2*f(Q(xs[i]+xs[j], 2))
               for i in range(len(xs)) for j in range(i+1, len(xs))}
    @lru_cache(None)
    def solve(mask):
        if not mask:
            return Q(0)
        i = (mask & -mask).bit_length()-1
        rest = mask ^ (1 << i)
        best = solve(rest)
        for j in range(i+1, len(xs)):
            if rest & (1 << j):
                best = max(best, weights[i, j]+solve(rest ^ (1 << j)))
        return best
    greedy = sum((weights[j, len(xs)-1-j] for j in range(len(xs)//2)), Q(0))
    need(greedy == solve((1 << len(xs))-1), ("independent optimal matching", xs, t))


def main():
    sys.stdout.reconfigure(newline="\n")
    scalar_cases = 0
    for length in range(2, 8):
        for xs in combinations(range(-3, 4), length):
            for t in (Q(0), Q(1, 2), Q(3, 2), Q(3)):
                matching_test(xs, t)
                scalar_cases += 1
    for xs in ((-2,-2,0,1,1), (-1,0,0,0,1), (0,0,0,0)):
        matching_test(xs, Q(1, 2))
        scalar_cases += 1
    records = []
    for w in ((19,23,29), (23,29,37), (29,35,41), (49,59,61), (1,137,277), (1,499,1001)):
        points = carriers(w)
        ds, bs, residual, curvatures = payment(w, points)
        invoice = Q(len(points))-Q(2*w[2], 11)
        records.append((w, len(points), tuple(map(str,ds)), tuple(map(str,bs)),
                        tuple(map(str,residual)), str(invoice), max(bs) >= invoice))
    sample = (5,7,11,13,17,31,43,61,79,101,137,199,499,1001)
    for m in sample:
        w = (1,m,2*m+3)
        ks, predicted = family(m)
        actual = carriers(w)
        need(actual == predicted, ("independent complete two-line dictionary", m))
        ds, bs, residual, curvatures = payment(w, actual)
        expected = Q((11*m-9)*len(ks)-42*sum(ks), 3*m)
        need(ds == (expected,Q(0),Q(0)) and bs == (Q(0),)*3,
             ("flat family deficit", m, ds, bs))
        need(all(kappa == 0 for row in curvatures for kappa in row),
             ("every eligible AP is curvature invisible", m))
    ks, points = family(137)
    ds, bs, _, _ = payment((1,137,277), points)
    need(ks == (17,20,23,26,29) and len(points) == 10 and ds[0] == Q(2660,411),
         "named ten-carrier affine-baseline hostile")
    print("status=PROVED_ANALYTICALLY local payment and complete flat family; FINITE_EXACT controls")
    print("scalar_optimal_matching_cases="+str(scalar_cases))
    print("named_rows=(w,N,D,B,remainder,invoice,criterion_passes)")
    for row in records:
        print(row)
    print("family_controls="+str(sample))
    print("family_w=(1,m,2m+3), odd m>=5, 3 does not divide m")
    print("family_carriers=+- (m-3k,-2k-1,k), k=2 mod3, (5m-9)/42<k<(3m-1)/14")
    print("all_same_owner_midpoint_curvatures=0; positive_D1_is_affine_baseline")
    print("unbounded_cardinality=N=2(2m+3)/63+O(1)")
    print("N10_hostile_m137_k="+str(ks)+"; D1=2660/411; B=(0,0,0)")
    print("checks="+str(CHECKS))
    print("semantic_sha256="+sha256(json.dumps(records, sort_keys=True).encode()).hexdigest())
    print("all_remaining_count_dense_full_rank_dictionaries_satisfy_payment_condition=OPEN; LRC14=OPEN")


if __name__ == "__main__":
    main()
