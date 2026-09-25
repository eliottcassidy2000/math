#!/usr/bin/env python3
"""Independent audit of ternary_berggren section 5: consecutive nondegenerate
plus edge-triangles are never Berggren-comparable (ancestor/descendant at any
distance).  Own implementation of the parent map (B3) with accelerated B1-runs.
Also counts comparable consecutive pairs on the MINUS sheet (sign control) and
records where they occur.
"""
import sys


def U(n, s):
    m = 3 * n + s
    while m % 2 == 0:
        m //= 2
    return m


def ancestor_runs(s, t):
    """Yield runs (t, lo, hi): the set {(u, t): lo <= u <= hi, u == hi mod 2t}
    covering the whole ancestor chain of (s,t) INCLUDING (s,t) itself, up to the
    root (3,1).  Parent rules: (s-2t,t) if s>3t; (t,s-2t) if 2t<s<3t;
    (t,2t-s) if t<s<2t."""
    while True:
        if s > 3 * t:
            # B1^{-1} run: s, s-2t, ..., down to s_f in (t, 3t]
            m = (s - t - 1) // (2 * t)  # number of steps keeping value > t
            sf = s - 2 * t * m
            # ensure sf > t and sf <= 3t
            while sf <= t:
                m -= 1
                sf = s - 2 * t * m
            while sf > 3 * t:
                m += 1
                sf = s - 2 * t * m
            yield (t, sf, s)
            s = sf
        else:
            yield (t, s, s)
        if s == 3 * t:  # root (3,1) (coprime forces t=1)
            return
        if 2 * t < s < 3 * t:
            s, t = t, s - 2 * t
        elif t < s < 2 * t:
            s, t = t, 2 * t - s
        else:
            raise ValueError("bad pair", s, t)
        # after a non-B1 parent, loop continues


def is_ancestor(A, D):
    """True iff A is a (non-strict) ancestor of D (A may equal D)."""
    sa, ta = A
    for (t, lo, hi) in ancestor_runs(*D):
        if t == ta and lo <= sa <= hi and (hi - sa) % (2 * t) == 0:
            return True
        if hi < sa:
            # upper roots only decrease going up the chain
            return False
    return False


def edge_pair(x, y):
    return (max(x, y), min(x, y))


def run(sign, xmax, verbose_limit=40):
    counts = {}
    comparable = []
    total = 0
    for x in range(3, xmax + 1, 2):
        if sign == -1 and x == 1:
            continue
        y = U(x, sign)
        if y == 1 or y == x:
            continue
        z = U(y, sign)
        if sign == -1 and (y == 1):
            continue
        P1 = edge_pair(x, y)
        P2 = edge_pair(y, z)
        if P1[0] == P1[1] or P2[0] == P2[1]:
            continue
        total += 1
        kind = ("rise" if y > x else "fall") + "-" + ("rise" if z > y else "fall")
        counts[kind] = counts.get(kind, 0) + 1
        comp = (P1 == P2) or is_ancestor(P1, P2) or is_ancestor(P2, P1)
        if comp:
            comparable.append((x, y, z, P1, P2))
    return total, counts, comparable


def selftest():
    # (27,5) = B1^2 (7,5)
    assert is_ancestor((7, 5), (27, 5))
    assert not is_ancestor((27, 5), (7, 5))
    # generate a Berggren tree to depth 7 and check ancestor relation exactly
    def children(s, t):
        return [(s + 2 * t, t), (2 * s + t, s), (2 * s - t, s)]
    nodes = {(3, 1): None}
    frontier = [(3, 1)]
    for _ in range(7):
        nf = []
        for v in frontier:
            for c in children(*v):
                nodes[c] = v
                nf.append(c)
        frontier = nf
    import random
    keys = list(nodes)
    rnd = random.Random(1)
    for _ in range(20000):
        A = rnd.choice(keys)
        D = rnd.choice(keys)
        # true ancestor via parent pointers
        v = D
        truth = False
        while v is not None:
            if v == A:
                truth = True
                break
            v = nodes[v]
        assert is_ancestor(A, D) == truth, (A, D)
    print("selftest: ancestor test agrees with explicit tree (depth 7,", len(keys), "nodes) PASS")


if __name__ == "__main__":
    selftest()
    xmax = int(sys.argv[1]) if len(sys.argv) > 1 else 100001
    for sign in (1, -1):
        total, counts, comp = run(sign, xmax)
        print("sign %+d odd x in [3,%d]: consecutive nondegenerate pairs %d; kinds %s; comparable %d"
              % (sign, xmax, total, dict(sorted(counts.items())), len(comp)))
        for c in comp[:40]:
            print("   comparable:", c)
