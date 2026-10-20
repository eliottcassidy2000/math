"""Small exact controls for the ordinary Hamiltonian-path spectrum audit.

No large-order census, random sampling, optimizer, or third-party dependency.
Normal and -O must have identical output.
"""
from collections import Counter
from itertools import combinations, permutations
from math import prod


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def rows(n, code):
    a = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        if code >> k & 1:
            a[i] |= 1 << j
        else:
            a[j] |= 1 << i
    return tuple(a)


def hamiltonian(a):
    n = len(a)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for end, count in enumerate(dp[mask]):
            todo = a[end] & ~mask if count else 0
            while todo:
                bit = todo & -todo
                todo -= bit
                dp[mask | bit][bit.bit_length() - 1] += count
    return sum(dp[-1])


def literal(a):
    return sum(all(a[u] >> v & 1 for u, v in zip(p, p[1:]))
               for p in permutations(range(len(a))))


def strong(a):
    n = len(a)
    full = (1 << n) - 1
    for start in range(n):
        seen = 1 << start
        frontier = seen
        while frontier:
            bit = frontier & -frontier
            frontier -= bit
            added = a[bit.bit_length() - 1] & ~seen
            seen |= added
            frontier |= added
        if seen != full:
            return False
    return True


def minimum(n):
    if n % 3 == 0:
        return 3 * 5 ** ((n - 3) // 3)
    if n % 3 == 1:
        return 5 ** ((n - 1) // 3)
    return 9 * 5 ** ((n - 5) // 3)


def cutoff(h):
    if h == 1:
        return 1
    n = 0
    while 3 ** (n + 1) <= h ** 3:
        n += 1
    return n


def strong_cutoff(h):
    n = 3
    while minimum(n) <= h:
        n += 1
    return n - 1


def factors(n):
    out = []
    p = 3
    while p * p <= n:
        while n % p == 0:
            out.append(p)
            n //= p
        p += 2
    if n > 1:
        out.append(n)
    return out


def prime(n):
    return n >= 3 and n % 2 and factors(n) == [n]


def lane_recipe(h):
    need(h > 0 and h % 2 and h not in (7, 21), "allowed target")
    a, m = 0, h
    while m % 7 == 0:
        a += 1
        m //= 7
    if not a:
        return factors(m)
    if a % 2 == 0:
        return [49] * (a // 2) + factors(m)
    if a >= 3:
        return [343] + [49] * ((a - 3) // 2) + factors(m)
    fs = factors(m)
    p = next((p for p in fs if p != 3), None)
    if p is not None:
        fs.remove(p)
        return [7 * p] + fs
    need(len(fs) >= 2, "excluded 7/21 reached")
    return [63] + fs[2:]


def circulant(n, steps):
    return tuple(sum(1 << ((i + d) % n) for d in steps) for i in range(n))


def reverse_arcs(a, arcs):
    b = list(a)
    for u, v in arcs:
        need(a[u] >> v & 1, "arc must be present before reversal")
        b[u] ^= 1 << v
        b[v] ^= 1 << u
    return tuple(b)


def automorphisms(a):
    return sum(all(bool(a[i] >> j & 1) == bool(a[p[i]] >> p[j] & 1)
                   for i, j in combinations(range(len(a)), 2))
               for p in permutations(range(len(a))))


def score_transport(a, target):
    need(tuple(x.bit_count() for x in a) == tuple(x.bit_count() for x in target),
         "transport needs the same labelled score vector")
    steps = 0
    while a != target:
        difference = tuple(a[i] & ~target[i] for i in range(len(a)))
        left = sum(row.bit_count() for row in difference)
        start = next(i for i, row in enumerate(difference) if row)
        path, position, v = [], {}, start
        while v not in position:
            position[v] = len(path)
            path.append(v)
            need(difference[v], "balanced disagreement cannot dead-end")
            v = (difference[v] & -difference[v]).bit_length() - 1
        cycle = path[position[v]:]
        arcs = tuple(zip(cycle, cycle[1:] + cycle[:1]))
        a = reverse_arcs(a, arcs)
        need(tuple(x.bit_count() for x in a) == tuple(x.bit_count() for x in target),
             "cycle transport preserves labelled scores")
        need(sum((a[i] & ~target[i]).bit_count() for i in range(len(a)))
             == left - len(cycle), "cycle removes precisely its disagreement edges")
        steps += 1
    return steps


def main():
    for n in range(1, 7):
        all_values, strong_values = set(), set()
        strong_number = 0
        for code in range(1 << (n * (n - 1) // 2)):
            a = rows(n, code)
            h = hamiltonian(a)
            need(h % 2 == 1 and h not in (7, 21), "small parity/exclusion")
            all_values.add(h)
            if n <= 5:
                need(h == literal(a), "independent permutation count")
            if strong(a):
                strong_number += 1
                strong_values.add(h)
        if n >= 3:
            need(min(strong_values) == minimum(n), "small Busch minimum")
        print("n=", n, "labelled=", 1 << (n * (n - 1) // 2),
              "strong_labelled=", strong_number,
              "ordinary_spectrum=", sorted(all_values),
              "strong_spectrum=", sorted(strong_values), sep="")
    need(strong_cutoff(7) == 4 and strong_cutoff(21) == 6,
         "Busch converts global exclusions to the checked small bank")
    need(all(minimum(n + 3) == 5 * minimum(n) for n in range(3, 100)),
         "correct recurrence")
    need(minimum(10) != 3 * minimum(8), "stale recurrence hostile")
    need(all(minimum(n) ** 3 >= 3 ** n for n in range(3, 100))
         and all(minimum(n) ** 3 > 3 ** n for n in range(4, 100)),
         "integer reduced-witness floor control")
    print("Busch minima n=3..14:", [minimum(n) for n in range(3, 15)])
    print("False 3*f(n-2) at n10:", minimum(10), "!=", 3 * minimum(8))
    print("Exclusion proof cutoff: strong H7 requires n<=4; strong H21 requires n<=6.")
    for h in (3, 7, 21, 49, 63, 343, 80557):
        print("h=", h, "SCC-reduced cutoff=", cutoff(h),
              "strong-lane cutoff=", strong_cutoff(h))
    for k in range(1, 20):
        need(cutoff(3 ** k) == 3 * k, "C3-chain sharp deletion-only cutoff")

    for h in range(1, 10002, 2):
        if h in (7, 21):
            continue
        fs = lane_recipe(h)
        need(prod(fs) == h, "symbolic lane recipe product")
        need(all(v in (49, 63, 343) or (prime(v) and v != 7)
                 or (v % 7 == 0 and prime(v // 7) and v // 7 != 3)
                 for v in fs), "only permitted conditional atoms")
    print("Conditional two-lane arithmetic recipe checked on all allowed odd h<=10001.")
    print("This checks arithmetic only; it does not construct the unproved lane witnesses.")

    c5 = circulant(5, (1, 2))
    need(hamiltonian(c5) == 15 and automorphisms(c5) == 5, "regular5 count types")
    regular5 = [a for a in (rows(5, code) for code in range(1 << 10))
                if all(row.bit_count() == 2 for row in a)]
    need(len(regular5) == 24, "labelled regular5 bank")
    distances = Counter(score_transport(a, b) for a in regular5 for b in regular5)
    print("Directed-cycle score-fibre transport: all 576 ordered regular5 pairs; step histogram=",
          sorted(distances.items()))
    for u, v in combinations(range(5), 2):
        arc = (u, v) if c5[u] >> v & 1 else (v, u)
        flipped = reverse_arcs(c5, (arc,))
        need(sorted(x.bit_count() for x in flipped) == [1, 2, 2, 2, 3],
             "single-edge reversal exits regularity")
    c7 = circulant(7, (1, 2, 3))
    base_h = hamiltonian(c7)
    jumps, hostile = Counter(), None
    for triple in combinations(range(7), 3):
        u = triple[0]
        for v in triple[1:]:
            w = next(x for x in triple if x not in (u, v))
            if all(c7[x] >> y & 1 for x, y in ((u, v), (v, w), (w, u))):
                arcs = ((u, v), (v, w), (w, u))
                b = reverse_arcs(c7, arcs)
                need(all(row.bit_count() == 3 for row in b), "triangle preserves regularity")
                delta = hamiltonian(b) - base_h
                jumps[delta] += 1
                if abs(delta) > 2 and hostile is None:
                    hostile = (arcs, b, hamiltonian(b))
    need(hostile is not None, "cycle reversal has nonunit path-count jump")
    print("Regular carousel7 H=", base_h, "; directed-triangle-reversal deltas=", sorted(jumps.items()))
    print("Explicit regular nonunit jump: source rows=", c7, "triangle=", hostile[0],
          "target rows=", hostile[1], "H changes", base_h, "to", hostile[2])
    paley7 = circulant(7, (1, 2, 4))
    need(hamiltonian(paley7) == 189 and automorphisms(paley7) == 21,
         "ordinary count versus automorphism quotient")
    print("Type controls: regular5 H=15, Aut=5, H/Aut=3; Paley7 H=189, Aut=21, H/Aut=9.")
    print("PASS small all-labelled exclusions, independent counts, exact cutoff arithmetic, and regular-move hostiles.")


if __name__ == "__main__":
    main()
