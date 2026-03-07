#!/usr/bin/env python3
"""
Grinberg-Stanley U_T computation v2 -- corrected understanding.

Key insight from the paper (arXiv:2307.05569):
- For a digraph D on V, U_D = sum_{sigma in S_V} (-1)^{phi_D(sigma)} p_{type(sigma)}
  where phi_D(sigma) = #{i : (sigma_i, sigma_{i+1}) NOT in D}
- For a tournament T, every pair has exactly one arc direction, so
  phi_T(sigma) = #{i : sigma_{i+1} -> sigma_i in T} = #{non-descents w.r.t. T}

Actually, let me reconsider the definition. The issue is what "descent" means
in the Grinberg-Stanley context.

Their Theorem 1.38 states: For a tournament T, U_T can be written as a polynomial
in p_1, 2p_3, 2p_5, ... with NONNEG INTEGER coefficients.

But in the raw power-sum expansion, even-part partitions CAN appear -- they cancel
when you group terms properly.

NEW APPROACH: Instead of trying to match U_T directly, let me:
1. Compute the EXACT q-expansion (eliminate even parts via the involution)
2. Compare q-expansion coefficients with independence polynomial coefficients
3. Verify the principal specialization at specific points

Actually, the simplest check: for the odd-part partitions in U_T,
the coefficient pattern should encode the independence polynomial of Omega(T).

From n=3 data:
  T_3: U_T has p_1^3 + ... -2*p_3  =>  q_1^3 - q_3  (negative! not nonneg)
  C_3: U_T has p_1^3 + ... +2*p_3  =>  q_1^3 + q_3

Wait, T_3 has q_3 coefficient = -1, which is NEGATIVE. This contradicts Theorem 1.38.

RESOLUTION: I think my computation of phi is wrong. Let me re-read the definition.

In Grinberg-Stanley:
- D is a digraph on vertex set V = {v_1, ..., v_n}
- Dbar is the complement (all arcs NOT in D)
- sigma is a permutation of V
- phi(sigma) = #{i in {1,...,n-1} : (sigma(i), sigma(i+1)) is an arc of Dbar}
- BUT WAIT: they might use (sigma_i, sigma_{i+1}) or (sigma_{i+1}, sigma_i)

Let me try BOTH conventions and see which gives the right answer.

Also: The descent statistic for permutations is des(sigma) = #{i : sigma_i > sigma_{i+1}}.
For a tournament, a "T-descent" at position i means sigma_i -> sigma_{i+1} is an arc of T
and sigma_i > sigma_{i+1} in the natural order... no, that's not right either.

Let me just try: phi(sigma) = #{i : (sigma_{i+1}, sigma_i) is an arc of D}
i.e., the number of "ascents of D" = positions where the REVERSE arc is in D.
"""

from itertools import permutations, combinations
from collections import defaultdict
from fractions import Fraction


def transitive_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

def cyclic_tournament(n):
    assert n % 2 == 1
    A = [[0]*n for _ in range(n)]
    half = (n-1) // 2
    for i in range(n):
        for d in range(1, half+1):
            j = (i + d) % n
            A[i][j] = 1
    return A

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(edges):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def cycle_type(sigma):
    n = len(sigma)
    visited = [False]*n
    lengths = []
    for i in range(n):
        if not visited[i]:
            length = 0
            j = i
            while not visited[j]:
                visited[j] = True
                j = sigma[j]
                length += 1
            lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def compute_UT(A, n, phi_convention="forward"):
    """
    Compute U_T with different phi conventions.

    forward: phi(sigma) = #{i : A[sigma[i]][sigma[i+1]] == 0} (original)
    reverse: phi(sigma) = #{i : A[sigma[i+1]][sigma[i]] == 0}
    complement: phi(sigma) = #{i : A[sigma[i]][sigma[i+1]] == 1} (count arcs of D, not Dbar)
    """
    UT = defaultdict(int)
    for sigma in permutations(range(n)):
        if phi_convention == "forward":
            phi = sum(1 for i in range(n-1) if A[sigma[i]][sigma[i+1]] == 0)
        elif phi_convention == "reverse":
            phi = sum(1 for i in range(n-1) if A[sigma[i+1]][sigma[i]] == 0)
        elif phi_convention == "complement":
            phi = sum(1 for i in range(n-1) if A[sigma[i]][sigma[i+1]] == 1)
        elif phi_convention == "reverse_complement":
            phi = sum(1 for i in range(n-1) if A[sigma[i+1]][sigma[i]] == 1)
        ct = cycle_type(sigma)
        UT[ct] += (-1)**phi
    return dict(UT)


def hamiltonian_paths(A, n):
    paths = []
    for perm in permutations(range(n)):
        valid = all(A[perm[i]][perm[i+1]] == 1 for i in range(n-1))
        if valid:
            paths.append(perm)
    return paths

def descent_count(perm):
    return sum(1 for i in range(len(perm)-1) if perm[i] > perm[i+1])

def find_all_odd_cycles(A, n):
    cycles = []
    seen = set()
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                if all(A[perm[i]][perm[(i+1)%length]] == 1 for i in range(length)):
                    min_idx = list(perm).index(min(perm))
                    canonical = perm[min_idx:] + perm[:min_idx]
                    if canonical not in seen:
                        seen.add(canonical)
                        cycles.append(canonical)
    return cycles

def independence_polynomial(cycles):
    m = len(cycles)
    cycle_sets = [set(c) for c in cycles]
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycle_sets[i] & cycle_sets[j]:
                adj[i][j] = adj[j][i] = True

    alpha = defaultdict(int)
    for mask in range(2**m):
        verts = [i for i in range(m) if (mask >> i) & 1]
        indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            alpha[len(verts)] += 1

    max_k = max(alpha.keys()) if alpha else 0
    return [alpha.get(k, 0) for k in range(max_k + 1)]


def format_partition(partition):
    """Format partition as string like p_3*p_1^2."""
    from collections import Counter
    c = Counter(partition)
    parts = []
    for k in sorted(c.keys(), reverse=True):
        if c[k] == 1:
            parts.append(f"p_{k}")
        else:
            parts.append(f"p_{k}^{c[k]}")
    return "*".join(parts)


def analyze_all_conventions(A, n, name):
    print(f"\n{'='*60}")
    print(f"{name} (n={n})")
    print(f"{'='*60}")

    paths = hamiltonian_paths(A, n)
    cycles = find_all_odd_cycles(A, n)
    ip = independence_polynomial(cycles)

    print(f"H(T) = {len(paths)}, I(Omega,2) = {sum(c*2**k for k,c in enumerate(ip))}")
    print(f"I(Omega,x) = {ip}")
    print(f"#odd cycles = {len(cycles)}")

    for conv in ["forward", "complement", "reverse", "reverse_complement"]:
        UT = compute_UT(A, n, conv)
        print(f"\n  phi={conv}:")

        # Show only odd-part entries
        odd_entries = {k: v for k, v in UT.items() if all(p % 2 == 1 for p in k)}
        even_entries = {k: v for k, v in UT.items() if any(p % 2 == 0 for p in k) and v != 0}

        for partition in sorted(UT.keys()):
            coeff = UT[partition]
            if coeff != 0:
                print(f"    {coeff:+d} * {format_partition(partition)}")

        # Check q-basis (q_k = 2*p_k for odd k>=3)
        has_even = bool(even_entries)
        print(f"    Has even-part partitions: {has_even}")

        if odd_entries:
            print(f"    Odd-part q-expansion:")
            for partition in sorted(odd_entries.keys()):
                coeff = odd_entries[partition]
                if coeff == 0:
                    continue
                num_big = sum(1 for k in partition if k >= 3)
                q_coeff = Fraction(coeff, 2**num_big)
                q_parts = "*".join(f"q_{k}" for k in partition)
                print(f"      {q_coeff} * {q_parts}")


def main():
    print("Testing all phi conventions to find the one matching Theorem 1.38")
    print("(nonneg coefficients in q_1, q_3, q_5, ... basis)")

    # n=3
    T3 = transitive_tournament(3)
    C3 = cyclic_tournament(3)
    analyze_all_conventions(T3, 3, "T_3 transitive")
    analyze_all_conventions(C3, 3, "C_3 cyclic")

    # n=5
    T5 = transitive_tournament(5)
    C5 = cyclic_tournament(5)
    analyze_all_conventions(T5, 5, "T_5 transitive")
    analyze_all_conventions(C5, 5, "C_5 cyclic")

    # Now: find the RIGHT convention by checking which one gives ALL nonneg
    # q-coefficients for ALL tournaments of a given n

    print("\n\n" + "="*60)
    print("EXHAUSTIVE CHECK: which convention gives nonneg q-coefficients for ALL n=3 tournaments?")
    print("="*60)

    for conv in ["forward", "complement", "reverse", "reverse_complement"]:
        all_nonneg = True
        for A in all_tournaments(3):
            UT = compute_UT(A, 3, conv)
            for partition, coeff in UT.items():
                if all(p % 2 == 1 for p in partition):
                    num_big = sum(1 for k in partition if k >= 3)
                    q_coeff = Fraction(coeff, 2**num_big)
                    if q_coeff < 0:
                        all_nonneg = False
                        break
            if not all_nonneg:
                break
        print(f"  {conv}: all nonneg = {all_nonneg}")

    print("\n\n" + "="*60)
    print("EXHAUSTIVE CHECK: which convention gives nonneg q-coefficients for ALL n=5 tournaments?")
    print("="*60)

    for conv in ["forward", "complement", "reverse", "reverse_complement"]:
        all_nonneg = True
        violation_count = 0
        for A in all_tournaments(5):
            UT = compute_UT(A, 5, conv)
            for partition, coeff in UT.items():
                if all(p % 2 == 1 for p in partition):
                    num_big = sum(1 for k in partition if k >= 3)
                    q_coeff = Fraction(coeff, 2**num_big)
                    if q_coeff < 0:
                        all_nonneg = False
                        violation_count += 1
        print(f"  {conv}: all nonneg = {all_nonneg}, violations = {violation_count}")

    # Key question: Does ANY convention make the q-coefficients match the
    # independence polynomial of Omega(T)?

    print("\n\n" + "="*60)
    print("KEY TEST: q-coefficient of q_3 vs alpha_1(Omega) for all n=3 tournaments")
    print("="*60)

    for conv in ["forward", "complement", "reverse", "reverse_complement"]:
        print(f"\n  Convention: {conv}")
        match_count = 0
        total = 0
        for A in all_tournaments(3):
            n = 3
            UT = compute_UT(A, n, conv)
            cycles = find_all_odd_cycles(A, n)
            ip = independence_polynomial(cycles)
            alpha_1 = ip[1] if len(ip) > 1 else 0

            p3_coeff = UT.get((3,), 0)
            q3_coeff = Fraction(p3_coeff, 2)

            total += 1
            if q3_coeff == alpha_1:
                match_count += 1
        print(f"    q_3 == alpha_1: {match_count}/{total}")

    # Maybe the relationship is different. Let me check the TOTAL sum
    # sum over all odd-part partitions lambda of q-coeff(lambda)
    # vs I(Omega, 1) or I(Omega, 2)

    print("\n\n" + "="*60)
    print("DEEPER: Sum of ALL q-coefficients vs I(Omega, x) for n=3")
    print("="*60)

    conv = "forward"  # Use the original convention
    for A in all_tournaments(3):
        n = 3
        UT = compute_UT(A, n, conv)
        cycles = find_all_odd_cycles(A, n)
        ip = independence_polynomial(cycles)

        # Sum of q-coefficients for odd-part partitions
        q_sum = Fraction(0)
        q_details = {}
        for partition, coeff in UT.items():
            if all(p % 2 == 1 for p in partition):
                num_big = sum(1 for k in partition if k >= 3)
                q_coeff = Fraction(coeff, 2**num_big)
                q_sum += q_coeff
                q_details[partition] = q_coeff

        ip_at_1 = sum(ip)
        ip_at_2 = sum(c * 2**k for k, c in enumerate(ip))

        A_flat = tuple(tuple(row) for row in A)
        print(f"  T={A_flat}")
        print(f"    q-coeffs: {q_details}")
        print(f"    sum(q-coeffs) = {q_sum}")
        print(f"    I(Omega,1) = {ip_at_1}, I(Omega,2) = {ip_at_2} = H(T)")
        print(f"    q_(1,1,1) + q_(3) = {q_details.get((1,1,1),0)} + {q_details.get((3,),0)} = {q_details.get((1,1,1),0) + q_details.get((3,),0)}")
        print()


    # Actually, maybe I should look at this from the OPPOSITE direction.
    # The Grinberg-Stanley formula for the CHROMATIC symmetric function X_G
    # involves the graph complement. For a tournament T, what's the right object?

    # Their Theorem 1.30 says U_D for a digraph D.
    # For a tournament T:
    #   U_T = sum_{sigma in S_n} (-1)^{phi_T(sigma)} p_{type(sigma)}
    # where phi_T(sigma) = #{i : (sigma_i, sigma_{i+1}) in Tbar}
    #                     = #{i : T[sigma_i][sigma_{i+1]] = 0 and T[sigma_{i+1}][sigma_i] = 0}
    #
    # Wait! For a TOURNAMENT, Tbar has NO arcs between vertices that have arcs in T.
    # In fact, Tbar is the EMPTY digraph (since T is a tournament = complete).
    # So phi_T(sigma) = n-1 for ALL sigma?? That can't be right.
    #
    # I think the issue is: the complement is the complement as a SIMPLE digraph.
    # T has arc (i,j) for each pair. The complement Tbar has arc (i,j) iff
    # T does NOT have arc (i,j). But since T has exactly one of (i,j) or (j,i),
    # Tbar also has exactly one of (i,j) or (j,i) -- namely the OTHER one.
    # So Tbar = T^op (the opposite tournament).
    #
    # Therefore: phi_T(sigma) = #{i : (sigma_i, sigma_{i+1}) is an arc of T^op}
    #                         = #{i : T[sigma_{i+1}][sigma_i] = 1}
    #                         = #{i : the consecutive pair goes BACKWARD in T}
    #
    # For a Hamiltonian path sigma: phi_T(sigma) = 0 (all arcs forward in T)
    # So (-1)^{phi} = 1 for Hamiltonian paths, and the sign alternates otherwise.

    print("\n\n" + "="*60)
    print("CORRECTED: phi_T(sigma) = #{backward arcs} = #{i : T[sigma_{i+1}][sigma_i] = 1}")
    print("This means phi = 'reverse_complement' convention")
    print("="*60)

    conv = "reverse_complement"
    for A in all_tournaments(3):
        n = 3
        UT = compute_UT(A, n, conv)
        cycles = find_all_odd_cycles(A, n)
        ip = independence_polynomial(cycles)
        paths = hamiltonian_paths(A, n)
        ET = defaultdict(int)
        for p in paths:
            ET[descent_count(p)] += 1

        print(f"\n  T = {[row for row in A]}")
        print(f"  H(T)={len(paths)}, I(Omega,2)={sum(c*2**k for k,c in enumerate(ip))}, I(Omega,x)={ip}")
        print(f"  E_T(t) = {dict(ET)}")
        print(f"  U_T power-sum expansion:")
        for partition in sorted(UT.keys()):
            coeff = UT[partition]
            if coeff != 0:
                print(f"    {coeff:+d} * {format_partition(partition)}")

        # q-expansion
        print(f"  q-expansion (odd parts only):")
        for partition in sorted(UT.keys()):
            coeff = UT[partition]
            if coeff == 0 or any(p % 2 == 0 for p in partition):
                continue
            num_big = sum(1 for k in partition if k >= 3)
            q_coeff = Fraction(coeff, 2**num_big)
            q_parts = "*".join(f"q_{k}" for k in partition)
            print(f"    {q_coeff} * {q_parts}")


if __name__ == "__main__":
    main()
