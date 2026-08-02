#!/usr/bin/env python3
"""Exact controls for THM-3134 endpoint jets and C3 Newton substitution.

Repository files are not imported.  The two engines are:
  (1) coefficient powering of the proved C3 quotient-walk content kernel;
  (2) direct vertex-subset Hamiltonian-path/set-partition DP for n <= 9.
"""

from __future__ import annotations

from collections import defaultdict
from math import comb, factorial


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def walk_terms(bound: int) -> dict[tuple[int, int, int], int]:
    out: dict[tuple[int, int, int], int] = {}
    for a in range(bound + 1):
        for b in range(max(0, a - 1), min(bound, a + 1) + 1):
            for c in range(max(0, a - 1), min(bound, a + 1) + 1):
                if (a, b, c) != (0, 0, 0) and max(a, b, c) - min(a, b, c) <= 1:
                    out[(a, b, c)] = 3 if a == b == c else 1
    return out


def compose_equal_f_profile(old_f: tuple[int, ...]) -> tuple[int, ...]:
    """Return F_new(d)=d!*pc_new(d) for C3[S,S,S]."""
    bound = len(old_f) - 1
    terms = walk_terms(bound)
    power: dict[tuple[int, int, int], int] = {(0, 0, 0): 1}
    result = [0] * (3 * bound + 1)
    for d in range(1, 3 * bound + 1):
        nxt: dict[tuple[int, int, int], int] = defaultdict(int)
        for (a, b, c), value in power.items():
            for (u, v, w), coeff in terms.items():
                aa, bb, cc = a + u, b + v, c + w
                if aa <= bound and bb <= bound and cc <= bound:
                    nxt[(aa, bb, cc)] += value * coeff
        power = nxt
        result[d] = sum(
            value * old_f[a] * old_f[b] * old_f[c]
            for (a, b, c), value in power.items()
        )
    return tuple(result)


def compose_equal_newton(old_f: tuple[int, ...]) -> tuple[int, ...]:
    """Independent 1D Gregory--Newton recurrence derived from rational W_C3."""
    n = len(old_f) - 1
    max_d = 3 * n

    # B_s(k) = L(x^s(1+x)^k), where L(x^a)=F(a).
    b = [[0] * (max_d + 1) for _ in range(n + 1)]
    for s in range(n + 1):
        for k in range(max_d + 1):
            b[s][k] = sum(
                comb(k, r) * old_f[s + r]
                for r in range(min(k, n - s) + 1)
            )

    # delta[s][m] = Delta^m (B_s(k)^3)|_{k=0}.
    delta = [[0] * (max_d + 1) for _ in range(n + 1)]
    for s in range(n + 1):
        row = [value**3 for value in b[s]]
        delta[s][0] = row[0]
        for m in range(1, max_d + 1):
            row = [row[k + 1] - row[k] for k in range(len(row) - 1)]
            delta[s][m] = row[0]

    result = [0] * (max_d + 1)
    for d in range(1, max_d + 1):
        value = 0
        for j in range(min(d, n) + 1):
            m = d - j
            outer = comb(d, j) * 2**j
            value += outer * sum(
                comb(d + t - 1, t) * delta[t + j][m]
                for t in range(n - j + 1)
            )
        result[d] = value
    return tuple(result)


def diagonal_h(f: tuple[int, ...]) -> int:
    total = 0
    n = len(f) - 1
    for r in range(1, n + 1):
        total += 3 * f[r] ** 3
        if r + 1 <= n:
            total += 3 * (f[r + 1] * f[r] ** 2 + f[r + 1] ** 2 * f[r])
    return total


def c3_power(level: int) -> tuple[tuple[int, ...], ...]:
    # vertices are ternary words; first differing digit controls the arc.
    words = []
    for value in range(3**level):
        digits = []
        x = value
        for _ in range(level):
            digits.append(x % 3)
            x //= 3
        words.append(tuple(reversed(digits)))

    n = len(words)
    rows = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            k = next(k for k in range(level) if words[i][k] != words[j][k])
            arc = (words[j][k] - words[i][k]) % 3 == 1
            rows[i][j] = int(arc)
            rows[j][i] = int(not arc)
    return tuple(tuple(row) for row in rows)


def direct_profile(adj: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    """Independent subset DP: induced Hamiltonian paths then set partitions."""
    n = len(adj)
    size = 1 << n
    end = [[0] * n for _ in range(size)]
    h = [0] * size
    for v in range(n):
        end[1 << v][v] = 1
    for mask in range(1, size):
        for last in range(n):
            value = end[mask][last]
            if not value:
                continue
            h[mask] += value
            for nxt in range(n):
                if not (mask >> nxt) & 1 and adj[last][nxt]:
                    end[mask | (1 << nxt)][nxt] += value

    covers = [[0] * (n + 1) for _ in range(size)]
    covers[0][0] = 1
    for mask in range(1, size):
        root_bit = mask & -mask
        sub = mask
        while sub:
            if sub & root_bit:
                rest = mask ^ sub
                hs = h[sub]
                if hs:
                    for d in range(n):
                        if covers[rest][d]:
                            covers[mask][d + 1] += hs * covers[rest][d]
            sub = (sub - 1) & mask
    return tuple(covers[-1])


def main() -> None:
    profiles_f = [(0, 1)]
    h_values = [1]
    print("level=0 vertices=1 F=(0,1) H=1")
    for level in range(1, 5):
        old = profiles_f[-1]
        new = compose_equal_f_profile(old)
        if level <= 4:
            closed = compose_equal_newton(old)
            require(closed == new, f"Newton/kernel mismatch at level {level}")
            print(f"newton_closed_recurrence_level_{level}=PASS")
        profiles_f.append(new)
        h_values.append(new[1])
        group_order = 3 ** ((3**level - 1) // 2)
        require(
            all(value % group_order == 0 for value in new),
            f"wreath divisibility failure at level {level}",
        )
        print(f"wreath_divisibility_full_profile_level_{level}=PASS")
        require(new[1] == diagonal_h(old), f"diagonal H mismatch at level {level}")
        require(
            all(new[d] % factorial(d) == 0 for d in range(len(new))),
            f"factorial divisibility failure at level {level}",
        )
        print(f"level={level} vertices={3**level} H={new[1]}")
        print(f"pc={tuple(new[d] // factorial(d) for d in range(len(new)))}")
        if level <= 2:
            direct = direct_profile(c3_power(level))
            predicted = tuple(new[d] // factorial(d) for d in range(len(new)))
            require(direct == predicted, f"direct profile mismatch at level {level}")
            print(f"direct_profile_control_level_{level}=PASS")
    f5 = compose_equal_newton(profiles_f[-1])
    require(f5[1] == diagonal_h(profiles_f[-1]), "level-5 diagonal mismatch")
    profiles_f.append(f5)
    h_values.append(f5[1])
    require(
        all(value % (3 ** ((3**5 - 1) // 2)) == 0 for value in f5),
        "level-5 wreath divisibility failure",
    )
    print("wreath_divisibility_full_profile_level_5=PASS")
    print(f"level=5 vertices=243 H={f5[1]}")
    h6 = diagonal_h(f5)
    h_values.append(h6)
    print(f"level=6 vertices=729 H={h6}")
    for level, h in enumerate(h_values):
        group_exponent = (3**level - 1) // 2
        x, valuation = h, 0
        while x % 3 == 0:
            x //= 3
            valuation += 1
        print(
            f"valuation level={level} v3(H)={valuation} "
            f"wreath_lower_bound={group_exponent} "
            f"wreath_quotient_mod3={(h // 3**group_exponent) % 3}"
        )


if __name__ == "__main__":
    main()

# --- independent endpoint-jet and hostile controls ---

"""Independent hostile audit of the path-cover/Taylor-jet bridge.

Uses direct successor maps for path covers and direct permutation enumeration
for forward-edge counts.  It does not import the THM-3121 companion or the
iterated-profile probe.
"""

from itertools import permutations, product
from math import comb, factorial


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def tournaments(n: int):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for mask in range(1 << len(pairs)):
        a = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            a[i][j] = (mask >> k) & 1
            a[j][i] = 1 - a[i][j]
        yield tuple(tuple(row) for row in a), mask


def path_cover_profile(a: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    n = len(a)
    none = n
    choices = [tuple(v for v in range(n) if a[u][v]) + (none,) for u in range(n)]
    counts = [0] * (n + 1)
    for succ in product(*choices):
        indeg = [0] * n
        ok = True
        for v in succ:
            if v < n:
                indeg[v] += 1
                if indeg[v] > 1:
                    ok = False
                    break
        if not ok:
            continue
        for start in range(n):
            seen = set()
            v = start
            while succ[v] < n:
                if v in seen:
                    ok = False
                    break
                seen.add(v)
                v = succ[v]
            if not ok:
                break
        if ok:
            counts[sum(v == none for v in succ)] += 1
    return tuple(counts)


def forward_counts(a: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    n = len(a)
    counts = [0] * n
    for p in permutations(range(n)):
        good = sum(a[p[i]][p[i + 1]] for i in range(n - 1))
        counts[good] += 1
    return tuple(counts)


def jet_from_forward(forward: tuple[int, ...]) -> tuple[int, ...]:
    n = len(forward)
    # Palindromy permits indexing by the number b of bad adjacencies.
    return (0,) + tuple(
        sum(forward[b] * comb(n - 1 - b, d - 1 - b) for b in range(d))
        for d in range(1, n + 1)
    )


def forward_from_jet(jet: tuple[int, ...]) -> tuple[int, ...]:
    n = len(jet) - 1
    low = []
    for b in range(n):
        low.append(sum(
            (-1) ** (b - j) * comb(n - 1 - j, b - j) * jet[j + 1]
            for j in range(b + 1)
        ))
    return tuple(low)


def substitute_c3_two_singletons(a: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    n = len(a)
    sizes = (n, 1, 1)
    offsets = (0, n, n + 1, n + 2)
    out = [[0] * (n + 2) for _ in range(n + 2)]
    for i in range(n):
        for j in range(n):
            out[i][j] = a[i][j]
    c3 = ((0, 1, 0), (0, 0, 1), (1, 0, 0))
    for bi in range(3):
        for bj in range(3):
            if c3[bi][bj]:
                for u in range(offsets[bi], offsets[bi + 1]):
                    for v in range(offsets[bj], offsets[bj + 1]):
                        out[u][v] = 1
    return tuple(tuple(row) for row in out)


def held_karp_h(a: tuple[tuple[int, ...], ...]) -> int:
    n = len(a)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            value = dp[mask][last]
            if value:
                for nxt in range(n):
                    if not (mask >> nxt) & 1 and a[last][nxt]:
                        dp[mask | 1 << nxt][nxt] += value
    return sum(dp[-1])


def stirling_second(n: int, d: int) -> int:
    table = [[0] * (d + 1) for _ in range(n + 1)]
    table[0][0] = 1
    for i in range(1, n + 1):
        for j in range(1, min(i, d) + 1):
            table[i][j] = table[i - 1][j - 1] + j * table[i - 1][j]
    return table[n][d]


def main() -> None:
    checked = 0
    witness = {}
    for n in range(1, 6):
        for a, mask in tournaments(n):
            pc = path_cover_profile(a)
            f = tuple(factorial(d) * pc[d] for d in range(n + 1))
            forward = forward_counts(a)
            require(forward == tuple(reversed(forward)), f"palindromy n={n}, mask={mask}")
            require(jet_from_forward(forward) == f, f"jet mismatch n={n}, mask={mask}")
            require(forward_from_jet(f) == forward, f"inverse mismatch n={n}, mask={mask}")
            if n == 5 and mask in (40, 76):
                witness[mask] = (a, pc, held_karp_h(substitute_c3_two_singletons(a)))
            checked += 1

    require(witness[40][1] == (0, 15, 39, 33, 10, 1), "mask-40 profile")
    require(witness[76][1] == (0, 15, 45, 35, 10, 1), "mask-76 profile")
    require(witness[40][2] == 123, "mask-40 composite H")
    require(witness[76][2] == 135, "mask-76 composite H")

    for n in range(1, 9):
        transitive = tuple(tuple(int(i < j) for j in range(n)) for i in range(n))
        pc = path_cover_profile(transitive)
        require(
            pc == (0,) + tuple(stirling_second(n, d) for d in range(1, n + 1)),
            f"transitive Stirling boundary n={n}",
        )

    print(f"labelled_tournaments_jet_and_inverse={checked}:PASS")
    print("same_order_same_H_hostile=mask40:pc(15,39,33,10,1)->123;")
    print("                          mask76:pc(15,45,35,10,1)->135:PASS")
    print("transitive_boundary_pc_equals_Stirling2=n1..8:PASS")


if __name__ == "__main__":
    main()
