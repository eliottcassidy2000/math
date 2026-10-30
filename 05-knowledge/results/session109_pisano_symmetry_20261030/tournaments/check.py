#!/usr/bin/env python3
"""Stdlib exact independent signed-constraint vs order/Burnside audit."""
import math, json
def require(condition, detail):
    if not condition:
        raise RuntimeError(f'Exact verification failed: {detail}')

def mm(A, B, n):
    return tuple((sum((A[2 * i + k] * B[2 * k + j] for k in range(2))) % n for i in range(2) for j in range(2)))

def mp(A, t, n):
    B = (1, 0, 0, 1)
    while t:
        if t & 1:
            B = mm(B, A, n)
        A = mm(A, A, n)
        t //= 2
    return B

def act(A, x, n):
    return ((A[0] * x[0] + A[1] * x[1]) % n, (A[2] * x[0] + A[3] * x[1]) % n)

def period(n):
    A = (1, 0, 0, 1)
    q = (1, 1, 1, 0)
    t = 0
    while True:
        A = mm(A, q, n)
        t += 1
        if A == (1, 0, 0, 1):
            return t

def constraint_exponent(n, A, reverse):
    X = {(a, b) for a in range(n) for b in range(n)} - {(0, 0)}
    color = {}
    components = 0
    for x in X:
        if x in color:
            continue
        components += 1
        color[x] = 0
        stack = [x]
        while stack:
            y = stack.pop()
            for z, e in [(act(A, y, n), reverse), ((-y[0] % n, -y[1] % n), 1)]:
                c = color[y] ^ e
                if z in color:
                    if color[z] != c:
                        return None
                else:
                    color[z] = c
                    stack.append(z)
    return components

def fib(t):
    a, b = (0, 1)
    for _ in range(t):
        a, b = (b, a + b)
    return (a, b)

def fixed_snf(n, t):
    if t == 0:
        return n * n
    f, fp = fib(t)
    fm = fp - f
    s = math.gcd(math.gcd(fp - 1, f), fm - 1)
    det = abs((fp - 1) * (fm - 1) - f * f)
    return math.gcd(n, s) * math.gcd(n, det // s)

def burnside(n, r, d):
    total = sum((fixed_snf(n, r * j) - 1 for j in range(d)))
    require(total % d == 0, 'total % d == 0')
    return total // d

def predicted(n, r, reverse):
    P = period(n)
    d = P // math.gcd(P, r)
    if reverse:
        if d % 4 != 2:
            return None
        if mp((1, 1, 1, 0), r * (d // 2), n) != (n - 1, 0, 0, n - 1):
            return None
        N = burnside(n, 2 * r, d // 2)
    else:
        if d % 2 == 0:
            return None
        N = burnside(n, r, d)
    require(N % 2 == 0, 'N % 2 == 0')
    return N // 2

def orbit_partition(n, A):
    unseen = {(a, b) for a in range(n) for b in range(n)} - {(0, 0)}
    out = []
    while unseen:
        x = next(iter(unseen))
        orb = set()
        y = x
        while y not in orb:
            orb.add(y)
            y = act(A, y, n)
        unseen -= orb
        out.append(orb)
    return out

def main():
    cases = 0
    examples = []
    for n in range(3, 26, 2):
        P = period(n)
        for r in range(P):
            A = mp((1, 1, 1, 0), r, n)
            for rev in [0, 1]:
                a = constraint_exponent(n, A, rev)
                b = predicted(n, r, rev)
                require(a == b, (n, r, rev, a, b))
                cases += 1
        examples.append({'n': n, 'period': P, 'Q_invariant': predicted(n, 1, 0), 'Q_reversing': predicted(n, 1, 1), 'Q2_invariant': predicted(n, 2, 0), 'Q2_reversing': predicted(n, 2, 1)})
    ram = []
    for k in range(1, 4):
        n = 5 ** k
        a = constraint_exponent(n, mp((1, 1, 1, 0), 4, n), 0)
        require(a == 5 ** k - 1, 'a == 5 ** k - 1')
        ram.append({'k': k, 'exponent': a})
    mixed = (1, 0, 0, 2)
    require(mp(mixed, 1, 3) != (2, 0, 0, 2), 'mp(mixed, 1, 3) != (2, 0, 0, 2)')
    require(constraint_exponent(3, mixed, 0) is None, 'constraint_exponent(3, mixed, 0) is None')
    require(constraint_exponent(3, mixed, 1) is None, 'constraint_exponent(3, mixed, 1) is None')
    n = 7
    H = (2, 0, 0, 4)
    E = constraint_exponent(n, H, 0)
    require(E == 8, 'E == 8')
    O = orbit_partition(7, H)
    lookup = {x: i for i, o in enumerate(O) for x in o}
    pairs = []
    used = set()
    for i, o in enumerate(O):
        if i in used:
            continue
        x = next(iter(o))
        j = lookup[-x[0] % 7, -x[1] % 7]
        pairs.append((i, j))
        used |= {i, j}
    for mask in range(1 << len(pairs)):
        S = set().union(*(O[p[mask >> j & 1]] for j, p in enumerate(pairs)))
        CS = {(b, a) for a, b in S}
        NS = {(-a % 7, -b % 7) for a, b in S}
        require(S != CS and S != NS and (CS != NS), 'S != CS and S != NS and (CS != NS)')
    attachments = []
    for p, a in [(5, 1), (13, 3)]:
        for k in [2, 3] if p == 5 else [2]:
            n = p ** k
            low = n // p
            Q = (a, 1, 1, 0)
            lo = orbit_partition(low, mp(Q, 4, low))
            lookup = {x: i for i, o in enumerate(lo) for x in o}
            hi = orbit_partition(n, mp(Q, 4, n))
            children = [0] * len(lo)
            kernel = 0
            for o in hi:
                reduced = {(x % low, y % low) for x, y in o}
                if reduced == {(0, 0)}:
                    kernel += 1
                    continue
                i = lookup[next(iter(reduced))]
                require(reduced == lo[i], 'reduced == lo[i]')
                require(len(o) == p * len(lo[i]), 'len(o) == p * len(lo[i])')
                children[i] += 1
            require(set(children) == {p} and kernel == 2 * (p - 1), 'set(children) == {p} and kernel == 2 * (p - 1)')
            require(len(hi) == 2 * (p ** k - 1), 'len(hi) == 2 * (p ** k - 1)')
            attachments.append({'p': p, 'a': a, 'k': k, 'parents': len(lo), 'children_each': p, 'kernel_pairs': kernel // 2})
    out = {'status': 'PASS', 'observer_choices_checked': 256, 'attachment_checks': attachments, 'signed_constraint_cases': cases, 'examples': examples, 'ramified_checks': ram, 'hostile_mixed_involution': 'PASS', 'split_observer_exponent': E}
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    main()
