#!/usr/bin/env python3
"""P7_burnside.py — Compute P(7) via Burnside formula by cycle type.

Session: kind-pasteur-2026-03-20-S5

P(n) = (1/n!) * sum_{sigma in S_n} c_T(sigma) * fix(sigma)

where c_T(sigma) = # tournaments fixed by sigma, fix(sigma) = # fixed vertices.

By grouping by cycle type, this becomes:
P(n) = (1/n!) * sum_{partitions lambda of n} N(lambda) * c_T(lambda) * a_1(lambda)

where N(lambda) = n! / (prod k^{a_k} * prod a_k!) is the number of permutations
with cycle type lambda, and a_1 = number of 1-cycles.

c_T(lambda) = 2^{d(lambda)} where d(lambda) = number of "free" edge orbits.

Also verify the formula P(n) = n*T(n) - 2*C(2(n-3), n-3).
"""

from math import factorial, gcd, comb
from itertools import combinations
from functools import reduce

def partitions(n, max_part=None):
    """Generate all partitions of n."""
    if max_part is None:
        max_part = n
    if n == 0:
        yield []
        return
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            yield [first] + rest

def cycle_type_to_a(partition, n):
    """Convert partition to a_k counts. a_k = number of k-cycles."""
    a = [0] * (n + 1)
    for part in partition:
        a[part] += 1
    return a

def num_perms_with_cycle_type(a, n):
    """Number of permutations in S_n with cycle type a."""
    result = factorial(n)
    for k in range(1, n + 1):
        result //= (k ** a[k]) * factorial(a[k])
    return result

def count_fixed_tournaments(a, n):
    """Count tournaments on n vertices fixed by a permutation with cycle type a.

    A tournament is fixed by sigma iff the arc orientation is consistent
    across each orbit of sigma on ordered pairs.

    For each orbit of sigma on UNORDERED pairs {u,v}:
    - If the orbit of (u,v) and (v,u) are the SAME orbit on ordered pairs:
      the orientation is FORCED (no free choice). This happens when sigma
      maps some element to its "reverse" within the orbit.
    - If the orbit of (u,v) and (v,u) are DIFFERENT orbits:
      we have a free choice (one determines the other). +1 free bit.

    For two vertices in the same k-cycle: {u, sigma^j(u)} for j=1,...,k-1.
    The orbit of (u, sigma^j(u)) under sigma has size k/gcd(j,k).
    The reverse (sigma^j(u), u) = sigma^j applied to (u, sigma^{k-j}(u)).
    This is in the orbit of (u, sigma^{k-j}(u)). So (u,sigma^j(u)) and
    (u,sigma^{k-j}(u)) are in "reverse" orbits.
    They're in the SAME orbit iff j = k-j mod k iff 2j = 0 mod k iff k|2j.
    For j=k/2 (when k is even): self-reverse orbit (forced).
    For j != k/2 (or k odd): paired with a different orbit (1 free choice per pair).

    Within a k-cycle: floor(k/2) unordered pair orbits.
    Of these: 1 self-reverse if k even (j=k/2), 0 if k odd.
    Free choices: (floor(k/2) - [k even]) / 2 + ... hmm, this needs care.

    Actually, let me think of it differently.
    Within a k-cycle on vertices v_0, v_1, ..., v_{k-1}:
    - Unordered pairs: {v_i, v_j} for 0 <= i < j < k
    - Orbit of {v_i, v_j} = orbit of {v_0, v_{j-i}} (by shifting)
    - So orbits are determined by the "distance" d = j-i, for d=1,...,floor(k/2)
    - For d < k/2: orbit of {v_0, v_d} has size k, contains k unordered pairs
      The reverse pair {v_d, v_0} has distance k-d, which is a DIFFERENT orbit
      → these two orbits are paired → 1 free choice per pair
    - For d = k/2 (k even): orbit has size k/2, and is its own reverse
      → orientation forced → 0 free choices

    Total free choices within a k-cycle:
    - k odd: floor(k/2) distances, all paired → floor(k/2)/2 = (k-1)/4? NO
      Wait: distances 1,...,(k-1)/2. Pairs: {d, k-d} for d=1,...,(k-1)/2.
      Number of such pairs: (k-1)/2 distances, BUT each pair {d,k-d} counts
      distances d and k-d. Since d < k/2 < k-d, each pair is unique.
      So floor((k-1)/2 / 2) = ... hmm, (k-1)/2 distances, paired into
      ((k-1)/2) / 2 free choices. But (k-1)/2 might be odd!
      Actually: the distances are 1, 2, ..., (k-1)/2. There are (k-1)/2 of them.
      They pair as {d, k-d}. Since d != k-d (because k is odd), each distance
      pairs with a different one. So we have (k-1)/2 / 2 pairs... but (k-1)/2
      might be odd. Wait, for k odd: (k-1)/2 distances. These pair as
      {1,k-1}, {2,k-2}, ..., {(k-1)/2, (k+1)/2}. But (k+1)/2 > (k-1)/2,
      so the last "pair" is {(k-1)/2, (k+1)/2}. ALL distances pair up,
      giving (k-1)/2 / 2 complete pairs. But (k-1)/2 is always even for
      k odd (since k-1 is even)? No: k=3: (k-1)/2 = 1 distance (d=1).
      Pairs: {1, 2}. That's 1 pair, 1 free choice. So (k-1)/2 = 1 pair? No,
      there's 1 distance and it pairs with k-1=2, but d=2 > (k-1)/2=1.
      So the orbit of distance 1 pairs with the orbit of distance 2.
      Free choices from 3-cycle: 1. But floor(3/2) = 1 unordered pair orbit.
      OK let me just compute directly.

    For a single k-cycle, the number of FREE binary choices for tournament
    orientations is floor(k/2) if k is odd, floor(k/2)-1 + 0 = (k-2)/2 if k even.
    Wait, I'll just count by example:
    - k=1: 0 pairs, 0 free (single vertex)
    - k=2: 1 pair {v_0,v_1}, distance 1=k/2, self-reverse → 0 free
    - k=3: 1 pair orbit (dist 1), paired with dist 2. 1 free choice.
    - k=4: 2 pair orbits (dist 1,2). Dist 2=k/2 is self-reverse (0 free).
            Dist 1 pairs with dist 3. 1 free choice. Total: 1.
    - k=5: 2 pair orbits (dist 1,2). {1,4} paired, {2,3} paired. 2 free? No,
            1 free per pair. 2 pairs → but we need to count free choices.
            Each pair {d, k-d} gives 1 free choice. There are floor((k-1)/2) / 2
            such pairs... actually for k=5: distances 1,2. Pairs: {1,4}, {2,3}.
            That's 2 unordered-pair-orbits, forming 1 reverse-pair. Wait, no:
            {1,4} is the pairing of orbit-of-dist-1 with orbit-of-dist-4=orbit-of-dist-1
            ... I'm confusing myself.

    Let me just compute d(lambda) directly by constructing the cycle representation.
    """
    # Build a concrete permutation with cycle type a
    perm = list(range(n))
    pos = 0
    for k in range(1, n + 1):
        for _ in range(a[k]):
            # Create a k-cycle starting at pos
            for i in range(k - 1):
                perm[pos + i] = pos + i + 1
            perm[pos + k - 1] = pos
            pos += k

    # Find orbits on ordered pairs
    # An ordered pair (i,j) maps to (perm[i], perm[j])
    visited = [[False]*n for _ in range(n)]
    free_choices = 0

    for i in range(n):
        for j in range(n):
            if i == j or visited[i][j]:
                continue
            # Trace the orbit of (i,j)
            orbit = []
            ci, cj = i, j
            while not visited[ci][cj]:
                visited[ci][cj] = True
                orbit.append((ci, cj))
                ci, cj = perm[ci], perm[cj]

            # Check if the reverse orbit (j,i) is the same orbit
            # The reverse of (i,j) is (j,i). Is (j,i) in this orbit?
            reverse_in_orbit = any(a == j and b == i for a, b in orbit)

            if not reverse_in_orbit:
                # Mark the reverse orbit too
                ci, cj = j, i
                while not visited[ci][cj]:
                    visited[ci][cj] = True
                    ci, cj = perm[ci], perm[cj]
                free_choices += 1
            else:
                # Self-reverse orbit: orientation is forced
                # But for tournaments, this means the orbit constrains
                # which orientation is consistent. If the orbit has
                # no consistent orientation, c_T = 0 for this sigma.
                # For tournaments, a self-reverse orbit of ordered pairs
                # means sigma maps (u,v) to (v,u) at some step, which
                # reverses the arc. For a tournament, arc(u,v) = 1 implies
                # arc(v,u) = 0. If sigma maps (u,v) to (v,u), then
                # consistency requires arc(u,v) = arc(v,u), contradiction!
                # So sigma fixes 0 tournaments if any self-reverse orbit exists?
                # NO — only if the ORBIT itself is self-reverse (contains both
                # (u,v) and (v,u)). If the orbit size is odd, the orbit can't
                # be self-reverse (since reversing changes the orbit structure).
                # Self-reverse means: orbit size is even, and the orbit contains
                # paired arcs.

                # For a self-reverse orbit: the orbit has (u,v) and (v,u).
                # In a fixed tournament: arc(u,v) and arc(v,u) must be
                # consistent along the orbit. Since arc(u,v) + arc(v,u) = 1,
                # this means the orbit forces one of two choices (arc(u,v)=1
                # or arc(u,v)=0). BUT: since the orbit maps (u,v) to (v,u),
                # we need arc(u,v) = arc(v,u) = contradiction! Unless...
                # Actually, sigma maps (u,v) to (sigma(u), sigma(v)). If
                # at some point sigma^k maps (u,v) to (v,u), then we need
                # arc(u,v) = arc(sigma(u), sigma(v)) = ... = arc(v,u) = 1-arc(u,v).
                # This means 0 = 1, contradiction. So c_T(sigma) = 0 when
                # there's a self-reverse orbit!

                # WAIT — this means c_T(sigma) = 0 for any sigma that has
                # a self-reverse ordered-pair orbit. Is this always the case
                # for permutations with even cycles?

                # For a 2-cycle (i,j): sigma swaps i and j.
                # (i,j) maps to (j,i) in one step. Self-reverse!
                # So c_T(sigma) = 0 for any sigma with a 2-cycle.

                # This means: only sigma with ALL ODD cycle lengths can fix
                # any tournament. This is a known result!

                free_choices = -1  # signal: c_T = 0
                break
        if free_choices == -1:
            break

    if free_choices == -1:
        return 0
    return 2 ** free_choices


def main():
    print("=" * 70)
    print("COMPUTING P(n) VIA BURNSIDE BY CYCLE TYPE")
    print("=" * 70)

    T_known = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}

    for n in range(1, 9):
        print(f"\n  --- n={n} ---")
        total = 0

        for partition in partitions(n):
            a = cycle_type_to_a(partition, n)
            N_lambda = num_perms_with_cycle_type(a, n)
            fix = a[1]  # number of fixed points = number of 1-cycles
            c_T = count_fixed_tournaments(a, n)

            contribution = N_lambda * c_T * fix

            if contribution > 0:
                print(f"    partition {partition}: N={N_lambda}, fix={fix}, "
                      f"c_T={c_T}, contribution={contribution}")

            total += contribution

        P_n = total // factorial(n)
        T_n = T_known.get(n, '?')
        nTn = n * T_known.get(n, 0)
        deficit = nTn - P_n if isinstance(T_n, int) else '?'
        k = n - 3
        central_binom = comb(2*k, k) if k >= 0 else 0
        formula_pred = nTn - 2 * central_binom if n >= 3 else P_n

        print(f"\n    P({n}) = {P_n}")
        print(f"    T({n}) = {T_n}")
        print(f"    n*T(n) = {nTn}")
        print(f"    Deficit = n*T(n) - P(n) = {deficit}")
        if n >= 3:
            print(f"    2*C(2*{k},{k}) = {2*central_binom}")
            print(f"    Formula: n*T(n) - 2*C(2(n-3),n-3) = {formula_pred}")
            print(f"    Match: {formula_pred == P_n}")

    print(f"\n  {'='*60}")
    print(f"  SUMMARY")
    print(f"  {'='*60}")
    print(f"\n  {'n':>4} {'T(n)':>8} {'P(n)':>8} {'n*T(n)':>8} {'D(n)':>8} {'2C(2k,k)':>8} {'Formula':>8}")
    for n in range(1, 9):
        # Recompute P(n) from Burnside
        total = 0
        for partition in partitions(n):
            a = cycle_type_to_a(partition, n)
            N_lambda = num_perms_with_cycle_type(a, n)
            fix = a[1]
            c_T = count_fixed_tournaments(a, n)
            total += N_lambda * c_T * fix
        P_n = total // factorial(n)

        T_n = T_known.get(n, 0)
        nTn = n * T_n
        D_n = nTn - P_n
        k = n - 3
        cb = 2 * comb(2*k, k) if k >= 0 else 0
        formula = nTn - cb if n >= 3 else P_n
        match = '  Y' if formula == P_n else '  N'

        print(f"  {n:4d} {T_n:8d} {P_n:8d} {nTn:8d} {D_n:8d} {cb:8d} {formula:8d}{match}")


if __name__ == "__main__":
    main()
