#!/usr/bin/env python3
"""THM-3181 independent subset audit with no permutation enumeration.

The backward distribution is computed by endpoint deletion.  The path-cover
profile is computed independently by Held--Karp path counts plus a canonical
set-partition recurrence.
"""

from functools import lru_cache
from math import comb, factorial


def require(condition, payload):
    if not condition:
        raise AssertionError(payload)


def tournament_from_code(n, code):
    beats = [[False] * n for _ in range(n)]
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> bit) & 1:
                beats[i][j] = True
            else:
                beats[j][i] = True
            bit += 1
    return tuple(tuple(row) for row in beats)


def transitive(n):
    return tuple(tuple(i < j for j in range(n)) for i in range(n))


def directed_triangle():
    return (
        (False, True, False),
        (False, False, True),
        (True, False, False),
    )


def order_join(left, right):
    a, b = len(left), len(right)
    out = [[False] * (a + b) for _ in range(a + b)]
    for i in range(a):
        for j in range(a):
            out[i][j] = left[i][j]
    for i in range(b):
        for j in range(b):
            out[a + i][a + j] = right[i][j]
    for i in range(a):
        for j in range(b):
            out[i][a + j] = True
    return tuple(tuple(row) for row in out)


def substitute(quotient, blocks):
    offsets = []
    total = 0
    for block in blocks:
        offsets.append(total)
        total += len(block)
    out = [[False] * total for _ in range(total)]
    for quotient_vertex, block in enumerate(blocks):
        offset = offsets[quotient_vertex]
        for i in range(len(block)):
            for j in range(len(block)):
                out[offset + i][offset + j] = block[i][j]
    for i in range(len(blocks)):
        for j in range(len(blocks)):
            if i == j or not quotient[i][j]:
                continue
            for u in range(len(blocks[i])):
                for v in range(len(blocks[j])):
                    out[offsets[i] + u][offsets[j] + v] = True
    return tuple(tuple(row) for row in out)


def backward_subset_dp(tournament, cap=None):
    """Return A_T(b), optionally only 0<=b<=cap, by last-vertex deletion."""
    n = len(tournament)
    if cap is None:
        cap = n - 1
    cap = min(cap, n - 1)
    size = 1 << n
    dp = [[None] * n for _ in range(size)]
    for v in range(n):
        polynomial = [0] * (cap + 1)
        polynomial[0] = 1
        dp[1 << v][v] = polynomial
    for mask in range(1, size):
        if mask & (mask - 1) == 0:
            continue
        for v in range(n):
            if not ((mask >> v) & 1):
                continue
            previous = mask ^ (1 << v)
            polynomial = [0] * (cap + 1)
            for u in range(n):
                source = dp[previous][u]
                if source is None:
                    continue
                shift = 1 if tournament[v][u] else 0
                for backward in range(cap + 1 - shift):
                    polynomial[backward + shift] += source[backward]
            dp[mask][v] = polynomial
    answer = [0] * (cap + 1)
    full = size - 1
    for v in range(n):
        for backward, value in enumerate(dp[full][v]):
            answer[backward] += value
    return tuple(answer)


def hamilton_counts(tournament):
    """Hamiltonian path count of every nonempty induced vertex subset."""
    n = len(tournament)
    size = 1 << n
    ending = [[0] * n for _ in range(size)]
    hamilton = [0] * size
    for v in range(n):
        ending[1 << v][v] = 1
        hamilton[1 << v] = 1
    for mask in range(1, size):
        if mask & (mask - 1) == 0:
            continue
        total = 0
        for v in range(n):
            if not ((mask >> v) & 1):
                continue
            previous = mask ^ (1 << v)
            value = 0
            for u in range(n):
                if ((previous >> u) & 1) and tournament[u][v]:
                    value += ending[previous][u]
            ending[mask][v] = value
            total += value
        hamilton[mask] = total
    return tuple(hamilton)


def path_cover_profile(tournament):
    """Unordered path-cover profile via a canonical least-vertex component."""
    n = len(tournament)
    hamilton = hamilton_counts(tournament)

    @lru_cache(maxsize=None)
    def covers(mask):
        if mask == 0:
            return (1,)
        least_bit = mask & -mask
        remainder = mask ^ least_bit
        out = [0] * (mask.bit_count() + 1)
        sub = remainder
        while True:
            component = sub | least_bit
            rest_profile = covers(mask ^ component)
            weight = hamilton[component]
            for d, count in enumerate(rest_profile):
                out[d + 1] += weight * count
            if sub == 0:
                break
            sub = (sub - 1) & remainder
        return tuple(out)

    return covers((1 << n) - 1)


def falling(t, d):
    value = 1
    for j in range(d):
        value *= t - j
    return value


def q_from_profile(profile, t):
    return sum(profile[d] * falling(t, d) for d in range(1, len(profile)))


def q_from_truncated_backward(tournament, k):
    """Compute Q_T(k) from only A_T(0),...,A_T(k-1)."""
    n = len(tournament)
    low = backward_subset_dp(tournament, k - 1)
    return sum(
        low[backward] * comb(n + k - backward - 1, n)
        for backward in range(k)
    )


def half_grid_recover(n, q_values):
    half = (n + 1) // 2
    answer = [None] * n
    for k in range(1, half + 1):
        value = sum(
            (-1) ** (k - m) * comb(n + 1, k - m) * q_values[m]
            for m in range(k + 1)
        )
        answer[k - 1] = value
        answer[n - k] = value
    return tuple(answer)


def selected_codes(n, count):
    universe = 1 << comb(n, 2)
    if count >= universe:
        return tuple(range(universe))
    values = {0, universe - 1}
    state = 0x6A09E667
    while len(values) < count:
        state = (1664525 * state + 1013904223) & 0xFFFFFFFF
        values.add(state % universe)
    return tuple(sorted(values))


def product(values):
    value = 1
    for item in values:
        value *= item
    return value


def main():
    cases = []
    samples = (
        (1, 1),
        (2, 2),
        (3, 8),
        (4, 64),
        (5, 128),
        (6, 48),
        (7, 24),
        (8, 8),
    )
    for n, count in samples:
        cases.extend((n, code) for code in selected_codes(n, count))

    independent_cases = 0
    truncated_checks = 0
    for n, code in cases:
        tournament = tournament_from_code(n, code)
        backward = backward_subset_dp(tournament)
        profile = path_cover_profile(tournament)
        require(
            sum(backward) == factorial(n),
            ("subset distribution", n, code, backward),
        )
        require(
            backward == tuple(reversed(backward)),
            ("palindrome", n, code, backward),
        )

        q_values = [
            q_from_profile(profile, m) for m in range((n + 1) // 2 + 1)
        ]
        require(q_values[0] == 0, (n, code, q_values))
        for m in range(1, len(q_values)):
            require(
                q_from_profile(profile, -m) == (-1) ** n * q_values[m],
                ("parity", n, code, m),
            )
            require(
                q_values[m] == q_from_truncated_backward(tournament, m),
                ("truncated deletion DP", n, code, m, q_values[m]),
            )
            truncated_checks += 1
        require(
            half_grid_recover(n, q_values) == backward,
            ("half recovery", n, code, q_values, backward),
        )
        independent_cases += 1

    # SCC-chain control through two independent computations.
    c3 = directed_triangle()
    k1 = transitive(1)
    chain = order_join(order_join(c3, k1), c3)
    chain_profile = path_cover_profile(chain)
    n = len(chain)
    half = (n + 1) // 2
    factors = (c3, k1, c3)
    product_values = [0] + [
        product(
            q_from_profile(path_cover_profile(factor), m) for factor in factors
        )
        for m in range(1, half + 1)
    ]
    direct_values = [q_from_profile(chain_profile, m) for m in range(half + 1)]
    require(
        product_values == direct_values,
        ("SCC Q product", product_values, direct_values),
    )
    require(
        half_grid_recover(n, product_values) == backward_subset_dp(chain),
        ("SCC half recovery", product_values),
    )

    # Cyclic-substitution hostile: scalar plethysm is false.
    lifted = substitute(c3, (c3, c3, c3))
    lifted_profile = path_cover_profile(lifted)
    lifted_h = lifted_profile[1]
    base_profile = path_cover_profile(c3)
    q_c3_at_1 = q_from_profile(base_profile, 1)
    naive_plethysm = q_from_profile(base_profile, q_c3_at_1)
    require(lifted_h == 3159, ("C3 cube direct H", lifted_h))
    require(
        naive_plethysm == 33 and lifted_h != naive_plethysm,
        ("cyclic scalar hostile", lifted_h, naive_plethysm),
    )
    lifted_backward = backward_subset_dp(lifted)
    lifted_half = (len(lifted) + 1) // 2
    lifted_q = [
        q_from_profile(lifted_profile, m) for m in range(lifted_half + 1)
    ]
    require(
        half_grid_recover(len(lifted), lifted_q) == lifted_backward,
        "C3 cube half-grid recovery",
    )

    print("TOURNAMENT HALF-GRID RECIPROCITY -- THM-3181 INDEPENDENT AUDIT")
    print("status=VERIFIED-EXACT;no-permutation-enumeration")
    print(
        f"selected_tournament_cases={independent_cases};"
        "orders=1..8;fixed_code_sampler=LCG"
    )
    print(
        "engines=last_vertex_backward_polynomial_DP + "
        "Held_Karp_paths + canonical_set_partition_profile"
    )
    print(
        f"truncated_deletion_DP_checks={truncated_checks};"
        "Q(k)_from_A(0..k-1);PASS"
    )
    print(
        "half_grid_parity_and_recovery=PASS;"
        "positive_and_negative_profile_evaluations=independent"
    )
    print(
        "SCC_chain=C3_join_K1_join_C3;"
        "factor_Q_jets=direct_profile_Q_jet;PASS"
    )
    print(
        f"cyclic_substitution_hostile=C3[C3,C3,C3];"
        f"direct_H={lifted_h};naive_Q_plethysm={naive_plethysm};REFUTED"
    )
    print(
        f"cyclic_positive_control=order_9_half_grid_values={lifted_q};"
        "direct_subset_A_recovered;PASS"
    )


if __name__ == "__main__":
    main()
