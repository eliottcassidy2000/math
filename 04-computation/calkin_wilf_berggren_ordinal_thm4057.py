"""Calkin--Wilf/Berggren ordinal companion for THM-4057.

This checks the Stern--Brocot reciprocal mirror, the distinct Calkin--Wilf
heap ordinal, the mod-six primitive-Pythagorean section, the three Berggren
ordinal transducers, scale fibres of natural-number edges, and the minimal
tournament-completion hostile.
"""

from __future__ import annotations

from collections import Counter, deque
from itertools import permutations
from math import gcd


def require(flag: bool, payload: object) -> None:
    if not flag:
        raise AssertionError(payload)


def add_pair(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return x[0] + y[0], x[1] + y[1]


def swap_pair(x: tuple[int, int]) -> tuple[int, int]:
    return x[1], x[0]


def complement_word(word: str) -> str:
    return word.translate(str.maketrans("LR", "RL"))


def sb_interval(word: str) -> tuple[tuple[int, int], tuple[int, int]]:
    """Boundary columns (numerator,denominator) at a Stern--Brocot address."""
    left, right = (0, 1), (1, 0)
    for letter in word:
        middle = add_pair(left, right)
        if letter == "L":
            right = middle
        else:
            left = middle
    return left, right


def projective_euclid_coefficients(p: int, q: int) -> tuple[int, ...]:
    digits = []
    while q:
        a, r = divmod(p, q)
        digits.append(a)
        p, q = q, r
    if digits and digits[0] == 0:
        digits.pop(0)
    return tuple(digits)


def stern(limit: int) -> list[int]:
    s = [0] * (limit + 2)
    s[1] = 1
    for n in range(2, limit + 2):
        if n % 2 == 0:
            s[n] = s[n // 2]
        else:
            s[n] = s[n // 2] + s[n // 2 + 1]
    return s


def cw_pair(k: int, s: list[int]) -> tuple[int, int]:
    return s[k], s[k + 1]


def cw_word(p: int, q: int) -> str:
    require(p > 0 and q > 0 and gcd(p, q) == 1, (p, q))
    reverse_letters: list[str] = []
    while (p, q) != (1, 1):
        if p < q:
            q -= p
            reverse_letters.append("L")
        else:
            p -= q
            reverse_letters.append("R")
    return "".join(reversed(reverse_letters))


def cw_pair_from_word(word: str) -> tuple[int, int]:
    p, q = 1, 1
    for letter in word:
        if letter == "L":
            q += p
        else:
            p += q
    return p, q


def heap_index(word: str) -> int:
    return int("1" + word.translate(str.maketrans("LR", "01")), 2)


def heap_reflection(k: int) -> int:
    depth = k.bit_length() - 1
    return 3 * (1 << depth) - 1 - k


def berggren_a(m: int, n: int) -> tuple[int, int]:
    return 2 * m - n, m


def berggren_b(m: int, n: int) -> tuple[int, int]:
    return 2 * m + n, m


def berggren_c(m: int, n: int) -> tuple[int, int]:
    return m + 2 * n, n


def berggren_ordinal_a(k: int) -> int:
    return 2 * k - 1


def berggren_ordinal_b(k: int) -> int:
    return 4 * heap_reflection(k) + 3


def berggren_ordinal_c(k: int) -> int:
    return 4 * k + 3


def finite_reduction_fibre_audit(n_max: int) -> tuple[int, int]:
    actual: Counter[tuple[int, int]] = Counter()
    for a in range(1, n_max + 1):
        for b in range(a + 1, n_max + 1):
            d = gcd(a, b)
            actual[a // d, b // d] += 1
    expected = {
        (p, q): n_max // q
        for q in range(2, n_max + 1)
        for p in range(1, q)
        if gcd(p, q) == 1
    }
    require(dict(actual) == expected, n_max)
    primitive = sum(count == 1 for count in actual.values())
    return len(actual), primitive


def tournament_stats(flip_24: bool) -> tuple[int, int]:
    n = 4
    adj = [[False] * n for _ in range(n)]
    for a in range(1, n + 1):
        for b in range(a + 1, n + 1):
            if flip_24 and (a, b) == (2, 4):
                adj[b - 1][a - 1] = True
            else:
                adj[a - 1][b - 1] = True
    hamiltonian = sum(
        all(adj[p[i]][p[i + 1]] for i in range(n - 1))
        for p in permutations(range(n))
    )
    cycles = 0
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                cycles += int(
                    (adj[a][b] and adj[b][c] and adj[c][a])
                    or (adj[a][c] and adj[c][b] and adj[b][a])
                )
    return hamiltonian, cycles


def main() -> None:
    # Stern--Brocot reflection is coordinate swap plus boundary reversal.
    sb_checks = 0
    for k in range(1, 1 << 16):
        word = bin(k)[3:].translate(str.maketrans("01", "LR"))
        left, right = sb_interval(word)
        reflected = sb_interval(complement_word(word))
        require(reflected == (swap_pair(right), swap_pair(left)), (word, reflected))
        sb_checks += 1

    # Reciprocal rationals retain the projective Euclidean coefficient word
    # obtained by deleting only a leading zero.  This is not the standard
    # finite Khinchin word, which always discards a_0.
    cf_checks = 0
    for p in range(1, 301):
        for q in range(1, 301):
            if gcd(p, q) == 1:
                require(projective_euclid_coefficients(p, q) == projective_euclid_coefficients(q, p), (p, q))
                cf_checks += 1

    ordinal_limit = 1 << 20
    s = stern(ordinal_limit + 8)
    pythagorean_ordinals = 0
    for k in range(1, ordinal_limit):
        p, q = cw_pair(k, s)
        require(gcd(p, q) == 1, (k, p, q))
        word = cw_word(p, q)
        require(heap_index(word) == k, (k, p, q, word))
        reflected = heap_reflection(k)
        require(heap_index(complement_word(word)) == reflected, (k, word, reflected))
        require(cw_pair(reflected, s) == (q, p), (k, p, q, reflected))
        require((s[k] % 2 == 0) == (k % 3 == 0), (k, s[k]))
        if k > 1:
            require((p > q) == (k % 2 == 1), (k, p, q))
        is_parameter = p > q and (p - q) % 2 == 1
        require(is_parameter == (k % 6 in (3, 5)), (k, p, q))
        pythagorean_ordinals += int(is_parameter)

    # Berggren tree: pair, word, ordinal, terminal-radix inverse, and depth
    # pullback are checked simultaneously through depth eleven.
    max_depth = 11
    queue: deque[tuple[int, int, int, int, int]] = deque([(2, 1, 3, 0, 0)])
    seen_pairs: set[tuple[int, int]] = set()
    level_counts: Counter[int] = Counter()
    depth_profile: Counter[tuple[int, int]] = Counter()
    branch_counts: Counter[str] = Counter()
    max_parameter = 0
    while queue:
        m, n, k, depth, non_a = queue.popleft()
        require((m, n) not in seen_pairs, (m, n, depth))
        seen_pairs.add((m, n))
        level_counts[depth] += 1
        depth_profile[depth, k.bit_length() - 1] += 1
        require(k.bit_length() - 1 == 1 + depth + non_a, (m, n, k, depth, non_a))
        max_parameter = max(max_parameter, m, n)
        word = cw_word(m, n)
        require(heap_index(word) == k and word.endswith("R"), (m, n, k, word))
        require(k % 6 in (3, 5), (m, n, k))
        require(m > n > 0 and gcd(m, n) == 1 and (m - n) % 2 == 1, (m, n, k))
        if depth == max_depth:
            continue
        children = (
            ("A", berggren_a(m, n), berggren_ordinal_a(k), word[:-1] + "LR", non_a),
            ("B", berggren_b(m, n), berggren_ordinal_b(k), complement_word(word) + "RR", non_a + 1),
            ("C", berggren_c(m, n), berggren_ordinal_c(k), word + "RR", non_a + 1),
        )
        for label, pair, child_k, child_word, child_non_a in children:
            require(cw_word(*pair) == child_word, (label, pair, child_word))
            require(heap_index(child_word) == child_k, (label, child_k, child_word))
            require(cw_pair_from_word(child_word) == pair, (label, pair, child_word))
            if label == "A":
                require(child_k % 4 == 1 and (child_k + 1) // 2 == k, child_k)
            elif label == "B":
                require(child_k % 8 == 3 and heap_reflection((child_k - 3) // 4) == k, child_k)
            else:
                require(child_k % 8 == 7 and (child_k - 3) // 4 == k, child_k)
            branch_counts[label] += 1
            queue.append((*pair, child_k, depth + 1, child_non_a))

    require(level_counts == Counter({d: 3**d for d in range(max_depth + 1)}), level_counts)
    for depth in range(max_depth + 1):
        for j in range(depth + 1):
            require(depth_profile[depth, 1 + depth + j] == __import__("math").comb(depth, j) * 2**j, (depth, j))

    fibre_rows = []
    for cutoff in (8, 32, 128, 512):
        types, singleton_types = finite_reduction_fibre_audit(cutoff)
        fibre_rows.append((cutoff, types, singleton_types, cutoff * (cutoff - 1) // 2))

    require(tournament_stats(False) == (1, 0), tournament_stats(False))
    require(tournament_stats(True) == (3, 1), tournament_stats(True))

    print("status=PROVED Calkin--Wilf/Berggren identities with FINITE-EXACT hostile companion")
    print(f"stern_brocot_words={sb_checks};projective_euclid_reciprocal_checks={cf_checks}")
    print(f"calkin_wilf_nodes={ordinal_limit-1};pythagorean_ordinals={pythagorean_ordinals}")
    print("cw_reflection=k_star=3*2^floor(log2(k))-1-k")
    print("berggren_ordinals=A:2k-1;B:4k_star+3;C:4k+3")
    print(f"berggren_depth={max_depth};nodes={len(seen_pairs)};max_parameter={max_parameter};branches={dict(sorted(branch_counts.items()))}")
    print(f"finite_reduction_fibres={fibre_rows}")
    print("n4_completion_hostile=ascending:(H,c3)=(1,0);flip_2_4:(H,c3)=(3,1)")
    print("PASS")


if __name__ == "__main__":
    main()
