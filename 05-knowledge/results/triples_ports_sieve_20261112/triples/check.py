"""Exact matrix checks for theta Cartan shadows and simplex gasket elimination."""
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
from math import lcm
from pathlib import Path
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def determinant(matrix):
    a = [[F(x) for x in row] for row in matrix]
    d = F(1)
    for i in range(len(a)):
        j = next((j for j in range(i, len(a)) if a[j][i]), None)
        if j is None:
            return F(0)
        if i != j:
            a[i], a[j] = a[j], a[i]
            d = -d
        pivot = a[i][i]
        d *= pivot
        for j in range(i+1, len(a)):
            c = a[j][i]/pivot
            for k in range(i+1, len(a)):
                a[j][k] -= c*a[i][k]
    return d


def theta(arms):
    n = sum(arms)-2
    c = [[int(i == j)*2 for j in range(n)] for i in range(n)]
    arm_vertices = []
    nxt = 1
    for s in arms:
        vs = list(range(nxt, nxt+s-1))
        for x, y in zip([0]+vs, vs):
            c[x][y] = c[y][x] = -1
        nxt += s-1
        arm_vertices.append(vs)
    return c, arm_vertices


def star_ldl(c, arm_vertices):
    # Leaf-to-root elimination is an independent exact inertia calculation.
    a = [[F(x) for x in row] for row in c]
    order = [v for vs in arm_vertices for v in reversed(vs)]+[0]
    pivots = []
    for t, v in enumerate(order):
        pivot = a[v][v]
        pivots.append(pivot)
        if pivot == 0:
            need(t == len(order)-1, 'unexpected zero pivot')
            continue
        active = order[t+1:]
        for i in active:
            for j in active:
                a[i][j] -= a[i][v]*a[v][j]/pivot
    return tuple(sum((x > 0, x < 0, x == 0)[k] for x in pivots) for k in range(3)), pivots[-1]


def gasket(q):
    inners = list(combinations(range(q), 2))
    labels = list(range(q))+inners
    idx = {x: i for i, x in enumerate(labels)}
    n = len(labels)
    lap = [[F(0) for _ in labels] for _ in labels]
    for i in range(q):
        cell = [idx[i]]+[idx[tuple(sorted((i,j)))] for j in range(q) if j != i]
        for u, v in combinations(cell, 2):
            lap[u][u] += 1
            lap[v][v] += 1
            lap[u][v] -= 1
            lap[v][u] -= 1
    # Solve all boundary basis data using the claimed harmonic extension,
    # then independently Schur-eliminate the literal Laplacian.
    for j in range(q):
        potentials = [F(int(i == j)) for i in range(q)]
        potentials += [F(1+int(i == j)+int(k == j), q+2) for i,k in inners]
        current = [sum(lap[i][k]*potentials[k] for k in range(n)) for i in range(n)]
        need(all(x == 0 for x in current[q:]), (q,'harmonic'))
        need(current[:q] == [F(q,q+2)*(q*int(i == j)-1) for i in range(q)], (q,'current'))
    for k in range(n-1,q-1,-1):
        pivot = lap[k][k]
        need(pivot > 0, 'interior positive')
        for i in range(k):
            for j in range(k):
                lap[i][j] -= lap[i][k]*lap[k][j]/pivot
    need(all(lap[i][j] == F(q,q+2)*(q*int(i == j)-1)
             for i in range(q) for j in range(q)), (q,'Schur'))


def main():
    count = 0
    zero = []
    for arms in combinations_with_replacement(range(1,13),3):
        p,q,r = arms
        c, vs = theta(arms)
        chi = F(1,p)+F(1,q)+F(1,r)-1
        need(determinant(c) == p*q+p*r+q*r-p*q*r, (arms,'determinant'))
        inertia, schur = star_ldl(c,vs)
        n = len(c)
        expected = (n,0,0) if chi>0 else ((n-1,1,0) if chi<0 else (n-1,0,1))
        need((inertia,schur) == (expected,chi), (arms,'inertia'))
        if not chi:
            root = lcm(*arms)
            vector = [root]+[root*(s-j)//s for s in arms for j in range(1,s)]
            need(all(sum(x*y for x,y in zip(row,vector)) == 0 for row in c), (arms,'null'))
            zero.append({'arms':arms,'null_vector':vector})
        count += 1
    for q in range(2,10):
        gasket(q)
    report = {'status':'FINITE-EXACT','sorted_triples':count,'max_arm':12,
              'zero_cases':zero,'gasket_q':list(range(2,10)),
              'theta_parity_hostile':[[2,3,r] for r in (4,6,8)]}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,sort_keys=True))


if __name__ == '__main__':
    main()
