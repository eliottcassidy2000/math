"""Exact labelled complement game census; no random or external graph package."""
from itertools import combinations, permutations
from math import comb, gcd


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def adjacency(n, edges):
    a = [0] * n
    for u, v in edges:
        a[u] |= 1 << v
        a[v] |= 1 << u
    return tuple(a)


def connected(a):
    seen, todo = 0, 1
    while todo:
        bit = todo & -todo
        todo ^= bit
        if seen & bit:
            continue
        seen |= bit
        todo |= a[bit.bit_length() - 1] & ~seen
    return seen == (1 << len(a)) - 1


def traceable(a):
    n = len(a)
    ends = [0] * (1 << n)
    for i in range(n):
        ends[1 << i] = 1 << i
    for mask in range(1, 1 << n):
        vertices = mask
        while vertices:
            bit = vertices & -vertices
            vertices ^= bit
            if ends[mask ^ bit] & a[bit.bit_length() - 1]:
                ends[mask] |= bit
    return bool(ends[-1])


def literal_traceable(a):
    return any(all(a[u] & (1 << v) for u, v in zip(p, p[1:]))
               for p in permutations(range(len(a))))


def minimum_sw_hp(n):
    m = n // 2
    return min(2 * n - 5, comb(m + 1, 2) if n % 2 == 0
               else comb(m + 2, 2) - 1)


def census(n):
    pairs = list(combinations(range(n), 2))
    size = 1 << len(pairs)
    hp, sw = bytearray(size), bytearray(size)
    for mask in range(size):
        a = adjacency(n, (e for i, e in enumerate(pairs) if mask >> i & 1))
        hp[mask], sw[mask] = traceable(a), connected(a)
        if n <= 5:
            require(hp[mask] == literal_traceable(a), (n, mask, "DP/permutation"))
    counts, minima = [0] * 4, [None] * 4
    for mask in range(size):
        co = size - 1 - mask
        predicates = (not hp[mask] and not hp[co],
                      not hp[mask] and not sw[co],
                      not sw[mask] and not hp[co],
                      not sw[mask] and not sw[co])
        for i, win in enumerate(predicates):
            if win:
                counts[i] += 1
                m = mask.bit_count()
                minima[i] = m if minima[i] is None else min(minima[i], m)
    expected = [n - 1, n - 1, minimum_sw_hp(n), None] if n >= 4 else [None] * 4
    require(minima == expected, (n, minima, expected))
    require(counts[1] == counts[2], "complement reciprocity")
    print("n=", n, "graphs=", size, "wins HH/HW/WH/WW=", counts,
          "minimum_red=", minima)


def extremal_family(n, t):
    # Blue: K_t join (I_(t+1) disjoint-union K_(n-2t-1)).
    c = set(range(t))
    r = set(range(2 * t + 1, n))
    return adjacency(n, ((u, v) for u, v in combinations(range(n), 2)
                         if u in c or v in c or (u in r and v in r)))


def constructions():
    for n in range(4, 101):
        candidates = []
        for t in range(1, n // 2):
            b = extremal_family(n, t)
            red = tuple(((1 << n) - 1) ^ (1 << u) ^ b[u] for u in range(n))
            require(connected(b) and not connected(red), (n, t, "connectivity"))
            cost = sum(row.bit_count() for row in red) // 2
            require(cost == comb(t + 1, 2) + (t + 1) * (n - 2 * t - 1),
                    "construction edge count")
            # Deleting universal blue vertices leaves t+2 nonempty components.
            require(n - 2 * t - 1 >= 1, "separator boundary")
            if n <= 12:
                require(not traceable(b), (n, t, "separator/DP"))
            candidates.append(cost)
        require(min(candidates) == minimum_sw_hp(n), (n, "endpoint minimum"))
    print("all extremal families n=4..100: counts/connectivity; n<=12: HP DP PASS")
    print("WH minimum n=4..20:", [minimum_sw_hp(n) for n in range(4, 21)])


def directed_hp_count(a):
    n = len(a)
    return sum(all(a[u] & (1 << v) for u, v in zip(p, p[1:]))
               for p in permutations(range(n)))


def ports_and_carries():
    # Undirected line graph loses which endpoint must be used at each transition.
    star_edges = [(0, 1), (0, 2), (0, 3)]
    line = adjacency(3, ((i, j) for i, j in combinations(range(3), 2)
                         if set(star_edges[i]) & set(star_edges[j])))
    require(traceable(line), "star line graph has a Hamiltonian path")
    require(all(sum(v in e for e in star_edges) % 2 for v in range(4)),
            "star has four odd vertices")
    count = 0
    for p in range(2, 31):
        for q in range(1, 101):
            nxt = [(i + q) % p for i in range(p)]
            carry = [(i + q) // p for i in range(p)]
            require(sum(carry) == q, "total winding")
            unseen = set(range(p))
            cycles = []
            while unseen:
                start = min(unseen)
                cycle, v = [], start
                while v in unseen:
                    unseen.remove(v)
                    cycle.append(v)
                    v = nxt[v]
                require(v == start, "permutation cycle")
                cycles.append(cycle)
            d = gcd(p, q)
            require(len(cycles) == d, "gcd component count")
            require(all(len(c) == p // d and sum(carry[v] for v in c) == q // d
                        for c in cycles), "cycle bidegree")
            if p <= 7 and q <= 14:
                hp = directed_hp_count([1 << v for v in nxt])
                require(hp == (p if d == 1 else 0), "permutation HP count")
            count += 1
    print("port hostile K1,3 -> K3 PASS; modular/carry pairs=", count, "PASS")


if __name__ == "__main__":
    for n in range(1, 7):
        census(n)
    constructions()
    ports_and_carries()
    print("PASS")
