"""Independent exact rational checks of the phase-two weighted stability bound.

No producer imports. Fraction Gaussian determinants and Hamilton subset DP
check the claimed lower bounds, including the augmented-skew cofactor route.
"""
from fractions import Fraction as F
from itertools import combinations, product
from hashlib import sha256
import json, random
from audit_feedback2 import determinant

def require(test, message):
    if not test:
        raise RuntimeError(message)

def uv(a):
    n = len(a)
    return (determinant([[a[i][j] + (i == j) for j in range(n)] for i in range(n)]),
            (-1)**n * determinant(a))

def hamilton(a):
    n = len(a)
    if not n:
        return F(1)
    dp = [[F(0)] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = F(1)
    for mask in range(1, 1 << n):
        for i, value in enumerate(dp[mask]):
            if value:
                for j in range(n):
                    if not mask >> j & 1:
                        dp[mask | 1 << j][j] += value * a[i][j]
    return sum(dp[-1])

def matrix(n, values):
    a = [[F(0)] * n for _ in range(n)]
    for (i, j), value in zip(combinations(range(n), 2), values):
        a[i][j], a[j][i] = value, 1-value
    return a

def run():
    cases = remainders = root_cases = 0
    digest = sha256()
    by_order = {}
    def audit(a):
        nonlocal cases, remainders, root_cases
        n = len(a)
        u, v = uv(a)
        gap = hamilton(a) - u - abs(v)
        k = [[a[i][j]-a[j][i] for j in range(n)] for i in range(n)]
        det_m = determinant([[k[i][j] + (i == j) for j in range(n)] for i in range(n)])
        augmented = [row+[F(1)] for row in k] + [[F(-1)]*n+[F(0)]]
        det_aug = determinant([[augmented[i][j]+(i == j) for j in range(n+1)] for i in range(n+1)])
        rooted = det_aug-det_m
        require(det_m >= 1 and rooted >= n, 'principal skew-minor lower bounds')
        require(u+v == det_m/F(2**(n-1)), 'D identity')
        require(u-v == rooted/F(2**(n-1)), 'sD cofactor identity')
        require(u-abs(v) >= F(1, 2**(n-1)), 'uniform determinant margin')
        if n >= 3:
            deletion_margins = {}
            for i, j in combinations(range(n), 2):
                vs = [t for t in range(n) if t not in (i,j)]
                b = [[a[s][t] for t in vs] for s in vs]
                ur, vr = uv(b)
                margin = ur-abs(vr)
                require(margin >= F(1, 2**(n-3)), 'pair-deletion uniform margin')
                deletion_margins[i,j] = margin
                remainders += 1
            root_variances = [sum(a[x][i]*a[i][x] for i in range(n) if i != x) for x in range(n)]
            root_lower = [sum(a[x][i]*a[i][x]*deletion_margins[tuple(sorted((x,i)))]
                              for i in range(n) if i != x) for x in range(n)]
            require(gap >= max(root_lower), 'exact pair-deletion stability bound')
            require(max(root_lower) >= F(1, 2**(n-3))*max(root_variances), 'uniform root stability bound')
            total = sum(a[i][j]*a[j][i] for i,j in combinations(range(n),2))
            require(F(1, 2**(n-3))*max(root_variances) >= F(2, n*2**(n-3))*total,
                    'averaged stability bound')
            root_cases += n
            digest.update(str((n,gap,root_lower,root_variances,total)).encode())
        cases += 1
        by_order[n] = by_order.get(n,0)+1

    for n in range(1,5):
        for values in product((F(0),F(1,2),F(1)), repeat=n*(n-1)//2):
            audit(matrix(n,values))
    rng = random.Random(1372026)
    for n in range(5,9):
        for _ in range(24):
            values = []
            for _ in combinations(range(n),2):
                d = rng.choice((2,3,4,5,7))
                values.append(F(rng.randrange(d+1),d))
            audit(matrix(n,values))
    # One fractional edge, both endpoints beating the sink: sharp at n=3.
    sharp = matrix(3,(F(1,3),F(1),F(1)))
    audit(sharp)
    u,v = uv(sharp)
    require(hamilton(sharp)-u-abs(v) == F(2,9), 'n3 sharp row-variance control')
    # n=2 is necessarily excluded: positive variance with zero full gap.
    hostile = matrix(2,(F(1,3),))
    u,v = uv(hostile)
    require(hamilton(hostile)-u-abs(v) == 0 and hostile[0][1]*hostile[1][0] == F(2,9),
            'n2 hostile exception')
    print(json.dumps({'scope':'all ternary weighted matrices n1..4;24 fixed-seed rational matrices each n5..8;one n3 sharp control',
                      'matrices':cases,'by_order':by_order,'pair_deletion_margins':remainders,
                      'root_stability_checks':root_cases,'sharp_n3_gap':'2/9','n2_hostile_variance':'2/9',
                      'semantic_sha256':digest.hexdigest(),
                      'verdict':'PASS exact pair-deletion, uniform root, averaged stability, augmented-skew cofactor lower bounds'},sort_keys=True,indent=2))

if __name__ == '__main__':
    run()
