#!/usr/bin/env python3
"""Finite independent controls for the signed fold and central difference."""
from itertools import combinations, permutations, product
from math import factorial


def require(condition, text):
    if not condition:
        raise RuntimeError(text)


def acyclic_paths_profile(P, n, q):
    h = [0] * (1 << n)
    h[0] = 1
    for walk in permutations(range(n)):
        mask = 1 << walk[0]
        h[mask] += 1
        for a, b in zip(walk, walk[1:]):
            if not P[a][b]:
                break
            mask |= 1 << b
            h[mask] += 1
    for mask in range(1, 1 << n):
        h[mask] //= factorial(n - mask.bit_count())
    def profile(slots):
        d = [0] * (1 << n)
        d[0] = 1
        for _ in range(slots):
            out = [0] * (1 << n)
            for mask in range(1 << n):
                sub = mask
                while True:
                    out[mask] += h[sub] * d[mask ^ sub]
                    if not sub:
                        break
                    sub = (sub - 1) & mask
            d = out
        return d[-1]
    return [profile(a) for a in range(q + 1)]


def adjacency(m, owners, label, b, odd):
    n = 2*m + odd
    P = [[0]*n for _ in range(n)]
    def arc(a, z):
        P[a][z] = 1
    C = n - 1
    for i in range(m):
        win = b[i] if not odd else 0
        arc(2*i+win, 2*i+(1^win))
        if odd:
            arc(2*i+b[i], C)
            arc(C, 2*i+(1^b[i]))
    for idx, (i,j) in enumerate(owners):
        p = 2*i+label[idx]
        z = p ^ 1
        for v in (2*j, 2*j+1):
            arc(p, v)
            arc(v, z)
    return P


def direct_fixed(P,m,b):
    return sum(
        signs[-1] == b[order[-1]]
        and all(P[2*i+s][2*j+t] for i,j,s,t in
                zip(order,order[1:],signs,signs[1:]))
        for order in permutations(range(m))
        for signs in product((0,1), repeat=m))


def coefficient(owners, m, mask):
    e = len(owners)
    rows = [[] for _ in range(m)]
    graph = [set() for _ in range(m)]
    terms = [i for i in range(m) if (mask >> (e+i)) & 1]
    if len(terms)>1:
        return 0
    for idx,(i,j) in enumerate(owners):
        if (mask >> idx) & 1:
            rows[i].append(j)
            graph[i].add(j)
            graph[j].add(i)
    s = 0
    for i, row in enumerate(rows):
        chosen = len(row)+(i in terms)
        if chosen not in (0,2) or len(graph[i])>2 or (i in terms and len(graph[i])!=1):
            return 0
        s += chosen==2
    visited = set()
    c = p = 0
    for i in range(m):
        if i in visited:
            continue
        unused = {i}
        component = set()
        while unused:
            j = unused.pop()
            component.add(j)
            unused.update(graph[j] - component)
        if sum(len(graph[j]) for j in component)//2 != len(component)-1:
            return 0
        visited.update(component)
        c += 1
        p += len(component)>1
    return (-1)**s * 2**(p-len(terms)) * factorial(c-len(terms))


def check_walsh():
    count = 0
    for m in range(1,5):
        pairs = list(combinations(range(m),2))
        ds = range(1 << len(pairs)) if m < 4 else [0,1,7,18,39,63]
        for orient_bits in ds:
            owners = [e[::(-1 if (orient_bits>>i)&1 else 1)]
                      for i,e in enumerate(pairs)]
            e = len(owners)
            values = []
            for mask in range(1<<(e+m)):
                label = [(mask>>i)&1 for i in range(e)]
                b = [(mask>>(e+i))&1 for i in range(m)]
                P = adjacency(m,owners,label,b,True)
                f = direct_fixed(P,m,b)
                values.append(f)
                require(f % 2 == 1, "even fixed count")
            norm = len(values)
            for i in range(m+e):
                for start in range(0,norm,2<<i):
                    for j in range(start,start+(1<<i)):
                        a,z = values[j], values[j+(1<<i)]
                        values[j],values[j+(1<<i)] = a+z,a-z
            for mask,v in enumerate(values):
                require(v == norm * coefficient(owners,m,mask), "Walsh coefficient")
            count += norm
    return count


def check_binary():
    count = 0
    for m in range(1,5):
        D = [(i,j) for i in range(m) for j in range(i)]
        for bits in range(1<<(m-1)):
            b = [0]+[(bits>>i)&1 for i in range(m-1)]
            for odd in (False,True):
                P = adjacency(m,D,[0]*len(D),b,odd)
                require(direct_fixed(P,m,b) == 2*bits+1, "binary fixed count")
                if odd or m==1:
                    continue
                old = [r[:-2] for r in P[:-2]]
                base = acyclic_paths_profile(old,2*m-2,4)
                lifted = acyclic_paths_profile(P,2*m,3)
                for q in (1,2,3):
                    expect = q*q*base[q]+b[-1]*q*(base[q+1]-base[q-1])
                    require(lifted[q]==expect, "central-difference profile")
                count += 1
    return count


if __name__ == "__main__":
    print(f"signed assignments checked: {check_walsh()}")
    print(f"active/inactive seam profiles checked: {check_binary()}")
