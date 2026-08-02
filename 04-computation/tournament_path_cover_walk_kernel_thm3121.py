#!/usr/bin/env python3
"""Exact referee for THM-3121 (path-cover substitution kernels).

The proof is a direct bijection.  This companion independently compares the
closed transfer formula with vertex-level enumeration on small hostile and
positive controls.  It uses integer arithmetic only and has no dependencies.
"""

from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from itertools import product
from math import factorial


Adj = tuple[tuple[int, ...], ...]
Profile = tuple[int, ...]  # entry d is the number of d-path covers


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def transitive(n: int) -> Adj:
    return tuple(tuple(int(i < j) for j in range(n)) for i in range(n))


C3: Adj = (
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 0),
)


def labelled_tournaments(n: int):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for mask in range(1 << len(pairs)):
        rows = [[0] * n for _ in range(n)]
        for bit, (i, j) in enumerate(pairs):
            rows[i][j] = (mask >> bit) & 1
            rows[j][i] = 1 - rows[i][j]
        yield tuple(tuple(row) for row in rows)


def substitute(quotient: Adj, blocks: tuple[Adj, ...]) -> Adj:
    sizes = [len(block) for block in blocks]
    offsets = [0]
    for size in sizes:
        offsets.append(offsets[-1] + size)
    total = offsets[-1]
    rows = [[0] * total for _ in range(total)]
    for i, block in enumerate(blocks):
        for u in range(sizes[i]):
            for v in range(sizes[i]):
                rows[offsets[i] + u][offsets[i] + v] = block[u][v]
        for j in range(len(blocks)):
            if i != j and quotient[i][j]:
                for u in range(sizes[i]):
                    for v in range(sizes[j]):
                        rows[offsets[i] + u][offsets[j] + v] = 1
    return tuple(tuple(row) for row in rows)


def hamiltonian_paths(adj: Adj) -> int:
    """Held--Karp count, independent of the path-cover enumerator below."""
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            value = dp[mask][last]
            if not value:
                continue
            for nxt in range(n):
                if not (mask >> nxt) & 1 and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += value
    return sum(dp[-1])


@lru_cache(maxsize=None)
def path_cover_profile(adj: Adj) -> Profile:
    """Count spanning covers by unordered directed paths via successor maps."""
    n = len(adj)
    none = n
    choices = [tuple(v for v in range(n) if adj[u][v]) + (none,) for u in range(n)]
    counts = [0] * (n + 1)
    for successor in product(*choices):
        indegree = [0] * n
        valid = True
        for nxt in successor:
            if nxt < n:
                indegree[nxt] += 1
                if indegree[nxt] > 1:
                    valid = False
                    break
        if not valid:
            continue
        for start in range(n):
            seen: set[int] = set()
            vertex = start
            while successor[vertex] < n:
                if vertex in seen:
                    valid = False
                    break
                seen.add(vertex)
                vertex = successor[vertex]
            if not valid:
                break
        if valid:
            counts[sum(int(nxt == none) for nxt in successor)] += 1
    return tuple(counts)


def add_count(counts: tuple[int, ...], vertex: int) -> tuple[int, ...]:
    out = list(counts)
    out[vertex] += 1
    return tuple(out)


@lru_cache(maxsize=None)
def walk_content(quotient: Adj, bounds: tuple[int, ...]) -> dict[tuple[int, ...], int]:
    """Coefficient table of W_Q: directed quotient walks by vertex content."""
    q = len(quotient)
    zero = (0,) * q
    ending: dict[tuple[tuple[int, ...], int], int] = {}
    for start in range(q):
        if bounds[start]:
            content = add_count(zero, start)
            ending[(content, start)] = 1

    # Every extension raises total content, so a total-degree sweep is acyclic.
    for degree in range(1, sum(bounds)):
        layer = [item for item in ending.items() if sum(item[0][0]) == degree]
        for (content, last), value in layer:
            for nxt in range(q):
                if quotient[last][nxt] and content[nxt] < bounds[nxt]:
                    new_content = add_count(content, nxt)
                    key = (new_content, nxt)
                    ending[key] = ending.get(key, 0) + value

    totals: dict[tuple[int, ...], int] = defaultdict(int)
    for (content, _last), value in ending.items():
        totals[content] += value
    return dict(totals)


def predicted_h(quotient: Adj, profiles: tuple[Profile, ...]) -> int:
    bounds = tuple(len(profile) - 1 for profile in profiles)
    total = 0
    for content, words in walk_content(quotient, bounds).items():
        if 0 in content:
            continue
        term = words
        for profile, count in zip(profiles, content):
            term *= factorial(count) * profile[count]
        total += term
    return total


def convolve(
    left: dict[tuple[int, ...], int],
    right: dict[tuple[int, ...], int],
    bounds: tuple[int, ...],
) -> dict[tuple[int, ...], int]:
    out: dict[tuple[int, ...], int] = defaultdict(int)
    for alpha, a_value in left.items():
        for beta, b_value in right.items():
            gamma = tuple(a + b for a, b in zip(alpha, beta))
            if all(g <= bound for g, bound in zip(gamma, bounds)):
                out[gamma] += a_value * b_value
    return dict(out)


def predicted_profile(quotient: Adj, profiles: tuple[Profile, ...]) -> Profile:
    """Full path-cover composition from exp(y W_Q)."""
    bounds = tuple(len(profile) - 1 for profile in profiles)
    q = len(bounds)
    zero = (0,) * q
    walk = walk_content(quotient, bounds)
    power: dict[tuple[int, ...], int] = {zero: 1}
    result = [0] * (sum(bounds) + 1)
    for paths in range(1, sum(bounds) + 1):
        power = convolve(power, walk, bounds)
        for content, ordered_words in power.items():
            if 0 in content:
                continue
            numerator = ordered_words
            block_choices = 1
            for profile, count in zip(profiles, content):
                numerator *= factorial(count)
                block_choices *= profile[count]
            require(
                numerator % factorial(paths) == 0,
                f"nonintegral set-of-walks kernel at d={paths}, content={content}",
            )
            result[paths] += numerator // factorial(paths) * block_choices
    return tuple(result)


def c3_start_count(content: tuple[int, int, int]) -> int:
    """Closed coefficient [x^content] W_C3."""
    low = min(content)
    high = max(content)
    if low <= 0 or high - low > 1:
        return 0
    return 3 if low == high else 1


def c3_diagonal_h(profiles: tuple[Profile, Profile, Profile]) -> int:
    transformed = [
        tuple(factorial(c) * profile[c] for c in range(len(profile)))
        for profile in profiles
    ]
    bound = min(len(profile) - 1 for profile in profiles)
    total = 0
    for r in range(1, bound + 1):
        total += 3 * transformed[0][r] * transformed[1][r] * transformed[2][r]
        for high in range(3):
            counts = [r, r, r]
            counts[high] += 1
            if all(counts[i] < len(transformed[i]) for i in range(3)):
                total += transformed[0][counts[0]] * transformed[1][counts[1]] * transformed[2][counts[2]]
        for low in range(3):
            counts = [r + 1, r + 1, r + 1]
            counts[low] -= 1
            if all(counts[i] < len(transformed[i]) for i in range(3)):
                total += transformed[0][counts[0]] * transformed[1][counts[1]] * transformed[2][counts[2]]
    return total


def verify_c3_kernel() -> int:
    bounds = (6, 6, 6)
    dynamic = walk_content(C3, bounds)
    checked = 0
    for content in product(range(1, 7), repeat=3):
        require(
            dynamic.get(content, 0) == c3_start_count(content),
            f"C3 kernel mismatch at content={content}",
        )
        checked += 1
    return checked


def verify_c3_blocks() -> int:
    blocks = (transitive(1), transitive(2), transitive(3), C3)
    checked = 0
    for triple in product(blocks, repeat=3):
        profiles = tuple(path_cover_profile(block) for block in triple)
        composite = substitute(C3, triple)
        direct = hamiltonian_paths(composite)
        require(direct == predicted_h(C3, profiles), "C3 general-kernel mismatch")
        require(direct == c3_diagonal_h(profiles), "C3 diagonal-law mismatch")
        checked += 1
    return checked


def verify_general_h() -> int:
    singleton = transitive(1)
    edge = transitive(2)
    checked = 0
    # All labelled four-vertex quotients and all 1/2-vertex block-size words.
    for quotient in labelled_tournaments(4):
        for blocks in product((singleton, edge), repeat=4):
            profiles = tuple(path_cover_profile(block) for block in blocks)
            require(
                hamiltonian_paths(substitute(quotient, blocks)) == predicted_h(quotient, profiles),
                "general quotient Hamiltonian-path mismatch",
            )
            checked += 1
    return checked


def verify_full_profiles() -> int:
    singleton = transitive(1)
    edge = transitive(2)
    controls = []
    for quotient in (transitive(2), transitive(3), C3):
        controls.extend((quotient, blocks) for blocks in product((singleton, edge), repeat=len(quotient)))
    controls.extend(
        [
            (C3, (C3, edge, singleton)),
            (C3, (C3, C3, singleton)),
            (transitive(3), (C3, edge, singleton)),
        ]
    )
    for quotient, blocks in controls:
        profiles = tuple(path_cover_profile(block) for block in blocks)
        direct = path_cover_profile(substitute(quotient, blocks))
        predicted = predicted_profile(quotient, profiles)
        require(direct == predicted, "full path-cover profile mismatch")
        require(
            direct[1] == hamiltonian_paths(substitute(quotient, blocks)),
            "path-cover/Held--Karp mismatch",
        )
    return len(controls)


def main() -> None:
    kernel_cells = verify_c3_kernel()
    c3_cases = verify_c3_blocks()
    general_cases = verify_general_h()
    profile_cases = verify_full_profiles()

    c3_profile = path_cover_profile(C3)
    cube = substitute(C3, (C3, C3, C3))
    direct_cube = hamiltonian_paths(cube)
    formula_cube = c3_diagonal_h((c3_profile, c3_profile, c3_profile))
    require(c3_profile == (0, 3, 3, 1), "unexpected C3 path-cover profile")
    require(direct_cube == formula_cube == 3159, "C3 substitution-cube mismatch")

    print("THM-3121 exact referee")
    print(f"C3 walk-kernel cells checked: {kernel_cells}")
    print("C3 kernel: K(a,b,c)=a!b!c! times (3 if a=b=c; 1 if max-min=1; 0 otherwise)")
    print(f"C3 ordered block triples checked against Held-Karp: {c3_cases}")
    print(f"general quotient/block controls checked against Held-Karp: {general_cases}")
    print(f"full path-cover profiles checked against exp(y W_Q): {profile_cases}")
    print(f"pc(C3)={c3_profile[1:]}; factorial transform=(3, 6, 6)")
    print(f"H(C3[C3,C3,C3]) direct=formula={direct_cube}")
    print("all checks passed")


if __name__ == "__main__":
    main()
