#!/usr/bin/env python3
"""Small exact controls for frame_and_design.md; standard library only."""
from fractions import Fraction
from itertools import product


def require(ok, message):
    if not ok:
        raise ValueError(message)


def mm(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def cyclic(n, steps):
    require(n >= 3 and n % 2 == 1, 'odd order required')
    require(len(steps) == (n-1)//2 and set(steps) | {n-x for x in steps} == set(range(1,n)), 'bad steps')
    return [[0 if i == j else (1 if (j-i) % n in steps else -1) for j in range(n)] for i in range(n)]


def check(k):
    n = len(k)
    require(all(len(row) == n for row in k), 'not square')
    require(all(k[i][j] == -k[j][i] and (k[i][j] in (-1,1) if i != j else k[i][j] == 0)
                for i in range(n) for j in range(n)), 'not skew tournament')
    require(all(sum(row) == 0 for row in k), 'regularity required')
    k2 = mm(k, k)
    d = [[k2[i][j]+n*(i == j)-1 for j in range(n)] for i in range(n)]
    require(all(d[i][i] == 0 and sum(d[i]) == 0 for i in range(n)), 'defect identities')
    a = [[int(x == 1) for x in row] for row in k]
    gram = mm(a, list(zip(*a)))
    for i in range(n):
        for j in range(n):
            require(4*gram[i][j] == (n+1)*(i == j)+(n-3)-d[i][j], 'design identity')
    energy = sum(x*x for row in d for x in row)
    if energy:
        i,j = next((i,j) for i in range(n) for j in range(n) if d[i][j] > 0)
        x = [int(t == i)-int(t == j) for t in range(n)]
        y = [sum(v*w for v,w in zip(row,x)) for row in k]
        require(sum(v*v for v in y) == 2*n+2*d[i][j], 'negative-vector identity')
        require(Fraction(n*sum(v*v for v in x)-sum(v*v for v in y),2) == -d[i][j], 'wrong witness value')
    return d, energy, gram


def main():
    count = 0
    for n in range(3,16,2):
        for signs in product((0,1), repeat=(n-1)//2):
            steps = {r if bit else n-r for r,bit in zip(range(1,(n+1)//2),signs)}
            check(cyclic(n,steps)); count += 1
    paley = cyclic(11,{1,3,4,5,9})
    d,e,g = check(paley)
    require(e == 0 and g[0] == [5]+[2]*10, 'Paley11 design/frame')
    complement = [[1-int(x == 1) for x in row] for row in paley]
    cg = mm(complement,list(zip(*complement)))
    require(all(cg[i][j] == 3*(i == j)+3 for i in range(11) for j in range(11)), 'complement design')
    d,e,g = check(cyclic(11,{1,2,3,4,5}))
    require(e == 3520 and d[0][4] == 4, 'cyclic11 exact obstruction')
    require(Fraction(e,16*121) == Fraction(20,11), 'projection defect normalization')
    bad = [row[:] for row in paley]
    bad[0][1] *= -1; bad[1][0] *= -1
    try:
        check(bad)
    except ValueError:
        pass
    else:
        raise ValueError('broken regularity not rejected')
    print('PASS:',count,'cyclic regular inputs; Paley11, complement design, cyclic11 negative witness; malformed input rejected')


if __name__ == '__main__':
    main()
