"""Exact tests of quadratic branch labels, extremal incidences, and their boundaries."""
from fractions import Fraction
from itertools import combinations
from math import comb
import subprocess


def decode(line, n):
    rows = [0] * n
    for bit, (i, j) in zip(line.strip(), combinations(range(n), 2)):
        a, b = (i, j) if bit == "1" else (j, i)
        rows[a] |= 1 << b
    return rows


def valid(rows):
    n = len(rows)
    incoming = [sum(1 << i for i in range(n) if rows[i] >> j & 1) for j in range(n)]
    cap = max((row & col).bit_count() for row in rows for col in incoming)
    for i, row in enumerate(rows):
        out = [k for k in range(n) if row >> k & 1]
        neighbors = {k: [l for l in out if k != l and rows[k] & rows[l]] for k in out}
        color = {}
        for start in out:
            if start in color:
                continue
            color[start] = 0
            queue = [start]
            while queue:
                u = queue.pop()
                for v in neighbors[u]:
                    if v not in color:
                        color[v] = color[u] ^ 1
                        queue.append(v)
                    elif color[v] == color[u]:
                        return cap, False
    return cap, True


def tournament_test():
    for n in range(3, 8):
        words = subprocess.check_output(["gentourng", "-q", str(n)],text=True).splitlines()
        data = [valid(decode(line,n)) for line in words]
        print("branches", n, "types", len(data), "two-route test", sum(k<=2 for k,_ in data),
              "bipartite test", sum(b for _,b in data))
        if n == 7:
            assert len(data) == 456 and sum(k<=2 for k,_ in data) == 1 and not any(b for _,b in data)
            paley = [sum(1<<((v+d)%7) for d in (1,2,4)) for v in range(7)]
            assert valid(paley) == (2,False)


def sharp_examples():
    # Rational parametrization of the unit circle. At j=3 it gives (4/5,3/5).
    X = set()
    for j in range(3,31):
        X.update((Fraction(j*j-1,j*j+1), Fraction(2*j,j*j+1)))
    first = {x:x**4 for x in X}
    second = {x:(x*x-1)**2 for x in X}
    assert len(X) == len(set(first.values())) == 56
    assert set(first.values()) == set(second.values()) and all(first[x]!=second[x] for x in X)

    # A skew translate of the projective trace-zero difference set from Session 35.
    D = {(x+3)%57 for x in (8,15,19,38,48,50,51,56)}
    assert D.isdisjoint({(-x)%57 for x in D})
    rows = [{(i+d)%57 for d in D} for i in range(57)]
    assert sum(map(len,rows)) == 456
    assert all(len(a&b)==1 for a,b in combinations(rows,2))
    # At the quadratic bound, the two convex minorants are simultaneously active.
    minimum = lambda m, f: (57-m%57)*f(m//57) + (m%57)*f(m//57+1)
    f = lambda v: comb(v,2) - (v*v)//4
    assert minimum(386,f) + minimum(386,lambda v: comb(v,2)) == 1593
    assert minimum(387,f) + minimum(387,lambda v: comb(v,2)) == 1602
    assert minimum(632,lambda v: comb(v,2)) <= 3192 < minimum(633,lambda v: comb(v,2))
    for cap, bound in ((1,456),(7,1157),(8,1235)):
        assert minimum(bound,lambda v: comb(v,2)) <= cap*1596 < minimum(bound+1,lambda v: comb(v,2))
    print("two disjoint rational quartic transports on 56 sheets; 456 oriented Singer incidences")


if __name__ == "__main__":
    tournament_test()
    sharp_examples()
