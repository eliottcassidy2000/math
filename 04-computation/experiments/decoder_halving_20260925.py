"""Canonical tournament halving and guarded ordinary-Collatz certificates.

Standard library; all finite universes are stated in main().
"""
from itertools import combinations


def odd_address(n):
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k, n


def odd_cores(max_q):
    cores = {1: [[False]]}
    for q in range(1, max_q, 2):
        old = cores[q]
        new = [[False] * (q + 2) for _ in range(q + 2)]
        for u in range(q):
            new[u][:q] = old[u]
        a, b = q, q + 1
        new[a][b] = True
        for v in range(q):
            if v < (q - 1) // 2:
                new[a][v] = new[v][b] = True
            else:
                new[b][v] = new[v][a] = True
        assert all(sum(row) == (q + 1) // 2 for row in new)
        cores[q + 2] = new
    return cores


def encode(n, cores):
    k, q = odd_address(n)
    size = 2**k
    H = cores[q]
    matrix = [[False] * n for _ in range(n)]
    for u, v in combinations(range(n), 2):
        bu, iu = divmod(u, size)
        bv, iv = divmod(v, size)
        matrix[u][v] = iu < iv if bu == bv else H[bu][bv]
        matrix[v][u] = not matrix[u][v]
    return matrix


def is_pair_module(H, u, v):
    return all(H[w][u] == H[w][v] for w in range(len(H)) if w not in (u, v))


def four_core(mask, H):
    q = len(H)
    base = [[False] * 4 for _ in range(4)]
    for bit, (u, v) in enumerate(combinations(range(4), 2)):
        base[u][v] = bool(mask & (1 << bit))
        base[v][u] = not base[u][v]
    N = 3 * q + 1
    F = [[False] * N for _ in range(N)]
    for u, v in combinations(range(N), 2):
        bu, iu = divmod(u, q)
        bv, iv = divmod(v, q)
        F[u][v] = H[iu][iv] if bu == bv else base[bu][bv]
        F[v][u] = not F[u][v]
    return F


def verify_certificate(target, word):
    state = 4
    for letter in word:
        if letter == "D":
            state *= 2
        elif letter == "O" and state % 6 == 4:
            state = (state - 1) // 3
        else:
            return False
    return state == target


def discover_odd_certificate(q):
    """Finite discovery only, not an assumed terminating universal routine."""
    state, seen, reverse_letters = q, set(), []
    for _ in range(10000):
        if state == 4:
            return "".join(reversed(reverse_letters))
        if state in seen:
            raise RuntimeError(f"Unexpected cycle from {q}")
        seen.add(state)
        if state % 2:
            reverse_letters.append("O")
            state = 3 * state + 1
        else:
            reverse_letters.append("D")
            state //= 2
    raise RuntimeError(f"Undecided start {q}; never silently omitted")


def main():
    cores = odd_cores(127)
    for q, H in cores.items():
        assert len(H) == q and all(sum(row) == (q - 1) // 2 for row in H)
        assert all(H[u][v] != H[v][u] for u, v in combinations(range(q), 2))
        if q > 1:
            assert [row[:q - 2] for row in H[:q - 2]] == cores[q - 2]
    print("Odd +2 cores: all 64 odd orders 1..127 are regular; each extends the previous core induced.")

    pair_count = 0
    for n in range(1, 129):
        E = encode(n, cores)
        k, q = odd_address(n)
        size = 2**k
        actual = {(u, v) for u, v in combinations(range(n), 2) if is_pair_module(E, u, v)}
        expected = {(u, u + 1) for u in range(n - 1) if u // size == (u + 1) // size}
        assert actual == expected
        assert len(actual) == n - q
        # Recover the unique matching intrinsically, by removing forced leaves.
        remaining = set(range(n))
        matched = []
        if k:
            while remaining:
                neighbors = {u: {v for v in remaining if tuple(sorted((u, v))) in actual}
                             for u in remaining}
                leaf = next(u for u in sorted(remaining) if len(neighbors[u]) == 1)
                mate = next(iter(neighbors[leaf]))
                matched.append(tuple(sorted((leaf, mate))))
                remaining.difference_update((leaf, mate))
            assert sorted(matched) == [(u, u + 1) for u in range(0, n, 2)]
        else:
            assert not actual
        pair_count += len(actual)
    print(f"Intrinsic pairs: all n=1..128; pair-module graph is q disjoint paths of order 2^k ({pair_count} pairs).")
    print("  Odd cores have no pairs; first doubles have disjoint edges; higher rows have unique perfect matchings.")

    for n in range(2, 129, 2):
        E, half = encode(n, cores), encode(n // 2, cores)
        for u in range(0, n, 2):
            assert is_pair_module(E, u, u + 1)
        quotient = [[E[2 * u][2 * v] for v in range(n // 2)] for u in range(n // 2)]
        assert quotient == half
    print("Canonical halving: every even n=2..128 has intrinsic pair modules and quotient exactly E(n/2).")

    checked = 0
    for q in range(3, 16, 2):
        for mask in range(64):
            F = four_core(mask, cores[q])
            assert not any(is_pair_module(F, u, v) for u, v in combinations(range(len(F)), 2))
            checked += 1
    print(f"Tripling hostile: {checked} substitutions, all 64 labeled four-cores and regular odd q=3..15;")
    print("  no two-vertex module in any; cannot equal canonical E(3q+1), which is pairable.")

    for n in range(1, 4097):
        k, q = odd_address(n)
        new_k, _ = odd_address(n + 2)
        assert (new_k == 0) if k == 0 else ((new_k >= 2) if k == 1 else (new_k == 1))
        assert 2 * (n + 2) == 2 * n + 4
    print("Literal +2 seam: all n=1..4096; odd stays odd, v2=1 enters v2>=2, v2>=2 enters v2=1.")

    certificates = {q: discover_odd_certificate(q) for q in range(1, 2001, 2)}
    lengths = []
    for n in range(1, 2001):
        k, q = odd_address(n)
        word = certificates[q] + "D" * k
        assert verify_certificate(n, word)
        # Independent forward check of every decoded edge, using the reversed word.
        x = n
        for symbol in reversed(word):
            if symbol == "D":
                assert x % 2 == 0
                x //= 2
            else:
                assert x % 2 == 1
                x = 3 * x + 1
        assert x == 4
        lengths.append(len(word))
    assert not verify_certificate(3, "DO")  # O at8 is forbidden
    assert not verify_certificate(4, "D")  # wrong endpoint
    assert not verify_certificate(4, "X")  # unknown transition
    print("Guarded decoder: all 2000 targets 1..2000 certified; 1000 odd witnesses then append doublings.")
    print(f"  Maximum certificate length={max(lengths)}; discovery limit 10000 steps, no filtered starts.")
    print(f"  Example 3: {certificates[3]}; invalid O-at-8, wrong endpoint, and unknown symbol rejected.")
    print("FINITE-EXACT discovery only. Tournament halving is universal; odd root-certificate existence stays OPEN.")


if __name__ == "__main__":
    main()
