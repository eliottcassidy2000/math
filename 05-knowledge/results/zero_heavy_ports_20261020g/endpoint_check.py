"""Exact marked-path tests for the (n-3,2) two-port endpoint layer.

Universe: every fixed-path blue middle of orders 4..8, every allowed port;
no isomorphism quotient and no implicit assertion-dependent checks.
"""
from collections import Counter
from itertools import combinations, permutations


def need(test, message):
    if not test:
        raise RuntimeError(message)


def arc(t, a, b):
    return (t[a] >> b) & 1


def path(t, p):
    return len(p) == len(t) and len(set(p)) == len(t) and all(
        arc(t, a, b) for a, b in zip(p, p[1:]))


def blue(t, p):
    n = len(p)
    return all(arc(t, p[i], p[j]) == arc(t, p[n-1-j], p[n-1-i])
               for i, j in combinations(range(n), 2))


def flip(t, edges):
    result = list(t)
    for a, b in edges:
        result[a] ^= 1 << b
        result[b] ^= 1 << a
    return tuple(result)


def opposite(t):
    n = len(t)
    return tuple(((1 << n)-1) ^ row ^ (1 << a) for a, row in enumerate(t))


def K(t, p):
    need(path(t, p), ("nonpath", t, p))
    return flip(opposite(t), zip(p, p[1:]))


def blue_presentations(n):
    free = {(a, b) for a in range(n) for b in range(a+2, n)}
    orbits = []
    while free:
        edge = min(free)
        a, b = edge
        orbit = {edge, (n-1-b, n-1-a)}
        orbits.append(orbit)
        free -= orbit
    for mask in range(1 << len(orbits)):
        rows = [0] * n
        for a in range(n-1):
            rows[a] |= 1 << (a+1)
        for bit, orbit in enumerate(orbits):
            for a, b in orbit:
                if (mask >> bit) & 1:
                    rows[a] |= 1 << b
                else:
                    rows[b] |= 1 << a
        yield tuple(rows)


def wrap(middle, x):
    r = len(middle)
    t = (sum(1 << v for v in range(1, r+1)),) + tuple(
        (row << 1) | (1 << (r+1)) for row in middle) + (1,)
    return flip(t, ((0, x), (r+1, r+1-x)))


def ledger(t, words, name):
    original = t
    for i, q in enumerate(words):
        need(path(t, q), (name, "direction", i, original, q))
        need(not blue(t, q), (name, "blue", i, original, q))
        t = K(t, q)
    need(t == K(original, tuple(range(len(original)))), (name, "wrong labelled target"))


def paths(t, allowed_starts=None, end=None):
    n = len(t)
    def visit(p, used):
        if len(p) == n:
            if end is None or p[-1] == end:
                yield p
            return
        for v in range(n):
            if used >> v & 1 or (v == end and len(p) != n-1):
                continue
            if (not p and (allowed_starts is None or v in allowed_starts)) or (
                    p and arc(t, p[-1], v)):
                yield from visit(p+(v,), used | (1 << v))
    yield from visit((), 0)


def strong(t):
    n = len(t)
    for start in range(n):
        seen = 1 << start
        while True:
            expanded = seen
            for v in range(n):
                if seen >> v & 1:
                    expanded |= t[v]
            if expanded == seen:
                break
            seen = expanded
        if seen != (1 << n)-1:
            return False
    return True


def twins(t, a, b):
    return all(arc(t, a, w) == arc(t, b, w)
               for w in range(len(t)) if w not in (a, b))


def full_tiling(n):
    return tuple(sum(1 << b for b in range(n) if b != a and
                     (b == a+1 or b < a-1)) for a in range(n))


def cycle_lengths(g):
    seen = set()
    lengths = []
    for v in range(len(g)):
        if v in seen:
            continue
        count = 0
        while v not in seen:
            seen.add(v)
            count += 1
            v = g[v]
        lengths.append(count)
    return sorted(lengths)


def compose_reflections(q):
    n = len(q)
    reflection = {q[i]: q[-1-i] for i in range(n)}
    return [reflection[n-1-i] for i in range(n)]


def main():
    counts = Counter()
    sidecars = Counter()
    first_missing = None
    for n in range(6, 11):
        a, b, c, z = 0, 1, n-2, n-1
        W = tuple(range(2, c))
        M = tuple(range(1, z))
        u, v = W[0], W[-1]
        for middle in blue_presentations(n-2):
            for x in range(2, z):
                y = n-1-x
                t = wrap(middle, x)
                p = tuple(range(n))
                need(path(t, p) and blue(t, p) and arc(t, z, a), "normal form")
                need((t[a].bit_count(), t[z].bit_count()) == (n-3, 2), "scores")
                need({i for i in range(n) if arc(t, i, a)} == {z, x}, "first port")
                need({i for i in range(n) if arc(t, z, i)} == {a, y}, "last port")
                counts[(n, "sources")] += 1
                if n >= 7:
                    need(not blue(t, M+(z,a)), "one-cut blackness")
                opened = arc(t, b, c)
                epsilon = arc(t, b, v)
                if opened and n >= 7 and x != v:
                    Vtail = tuple(range(2, z))
                    words = (M+(z,a), (z,)+Vtail+(b,a), (a,z)+Vtail+(b,))
                    ledger(t, words, "open")
                    counts[(n, "open")] += 1
                if not opened and not epsilon and n >= 8 and x != c:
                    words = ((c,)+W+(b,z,a), (z,c)+W+(b,a), (a,z,b)+W+(c,))
                    ledger(t, words, "closed0")
                    counts[(n, "closed0")] += 1
                if not opened and epsilon and n >= 8 and x not in (c,u):
                    words = ((c,b,z,a)+W, (z,c,b,a)+W, (a,z,b)+W+(c,))
                    ledger(t, words, "closed1")
                    counts[(n, "closed1")] += 1
                if not opened and epsilon and n >= 7 and x == c:
                    tail = W[1:]
                    words = ((a,)+W+(c,z,b), (z,)+tail+(c,u,b,a),
                             (b,z)+tail+(c,u,a))
                    ledger(t, words, "closed1-port-c")
                    counts[(n, "closed1-port-c")] += 1
                if not opened and x == c:
                    bd = all(arc(t,b,w) for w in W)
                    cz = all(arc(t,w,c) for w in W)
                    need(bd == cz == twins(t,a,b) == twins(t,c,z), "clone equivalence")
                    counts[(n, "clone-controls")] += 1
                if opened and x == v and 7 <= n <= 8:
                    C = K(middle, tuple(range(n-2)))
                    allowed = set(range(n-2)) - {c-1,y-1}
                    H = next(paths(C, allowed_starts=allowed, end=b-1), None)
                    if H is None:
                        sidecars[(n,"unavailable")] += 1
                        if first_missing is None:
                            first_missing = (t, x)
                    else:
                        H = tuple(v+1 for v in H)
                        ledger(t, (M+(z,a), (z,)+H+(a,), (a,z)+H), "open-sidecar")
                        sidecars[(n,"available")] += 1
        print("order", n, " ".join(name+"="+str(counts[(n,name)]) for name in
              ("sources","open","closed0","closed1","closed1-port-c","clone-controls")))

    print("open residual endpoint-sidecar", sorted(sidecars.items()))
    print("first unavailable open residual", first_missing)

    # Independent permutation controls for the symbolic all-order cycle tables.
    for n in range(6, 101):
        q = (n-2,1,n-1,0)+tuple(range(2,n-2))
        g = compose_reflections(q)
        if n % 2:
            expected = [n]
        elif n % 4 == 0:
            expected = [n//2,n//2]
        else:
            h = (n-2)//4
            expected = sorted([h,h,h+1,h+1])
        need(cycle_lengths(g) == expected, ("epsilon1 cycle table",n))
        if n >= 8:
            formula = [n-5,n-3,n-2,1,n-1,0]+list(range(2,n-6))+[n-4,n-6]
            need(g == formula, ("epsilon1 piecewise permutation",n,g,formula))
        q = (0,n-1)+tuple(range(1,n-1))
        g = compose_reflections(q)
        need(cycle_lengths(g) == ([n] if n%2 else [n//2,n//2]), "final cycle table")
        if n % 2 == 0:
            orbit = {0}
            a = g[0]
            while a != 0:
                orbit.add(a)
                a = g[a]
            need(orbit == {0} | set(range(1,n-2,2)), "final first orbit")
    print("reflection-composition controls: all orders 6..100 PASS")

    # Infinite obstruction to the particular endpoint-preserving two-block repair.
    for r in range(4, 25):
        n = r+2
        t = wrap(full_tiling(r), r)
        R = tuple(range(2,n))
        induced = tuple(sum(1 << j for j,w in enumerate(R) if arc(t,v,w)) for v in R)
        reversed_R = flip(induced, zip(range(len(R)-1), range(1,len(R))))
        need(reversed_R[0] == 1 << (len(R)-1), "full-tiling sole exit")
        need(next(paths(reversed_R, allowed_starts={0}, end=len(R)-1),None) is None,
             "full-tiling endpoint obstruction")
    print("full-tiling endpoint-repair obstruction: all r=4..24 PASS")

    # Exact small boundary: keep W=(2,3) contiguous and forward, permute five blocks.
    blocks = ((0,),(1,),(2,3),(4,),(5,))
    words = tuple(tuple(v for block in order for v in block) for order in permutations(blocks))
    controls = (("closed0",(14,4,40,50,39,3)),
                ("closed1-clones",(14,12,56,48,35,3)))
    for name,t in controls:
        p = tuple(range(6))
        target = K(t,p)
        need(blue(t,p) and strong(t) and strong(target), "small boundary type")
        target_edges = {K(target,q) for q in words if path(target,q) and not blue(target,q)}
        count = 0
        for q1 in words:
            if not path(t,q1) or blue(t,q1):
                continue
            U1 = K(t,q1)
            for q2 in words:
                if path(U1,q2) and not blue(U1,q2) and K(U1,q2) in target_edges:
                    count += 1
        need(count == 0, "restricted-template hostile")
        witness = next(q for q in paths(t) if not blue(t,q))
        witness_target = next(q for q in paths(target) if not blue(target,q))
        print("small boundary",name,"rows",t,"restricted triples",count,
              "black paths",witness,witness_target)
    clone = controls[1][1]
    need(blue(clone,(1,2,3,4,5,0)), "six-vertex blue rotation hostile")
    print("six-vertex shift-by-two/clone hostile PASS")
    print("ALL CHECKS PASS; finite controls support the stated proofs, not global connectivity.")


if __name__ == "__main__":
    main()
