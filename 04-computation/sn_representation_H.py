#!/usr/bin/env python3
"""
sn_representation_H.py -- kind-pasteur-2026-03-14-S68

How does S_n act on tournament space, and how does H decompose into irreps?

FRAMEWORK:
- S_n acts on tournaments by vertex relabeling:
    (sigma.T)[i][j] = T[sigma^{-1}(i)][sigma^{-1}(j)]
- H(T) is constant on orbits: H(sigma.T) = H(T) for all sigma.

- The VECTOR SPACE of functions f: {tournaments on [n]} -> R is a rep of S_n.
  Dimension = 2^{n choose 2} (number of labeled tournaments).
  S_n acts by (sigma.f)(T) = f(sigma^{-1}.T).

- The permutation representation on tournaments decomposes into irreps.
  Character: chi(sigma) = #{T : sigma.T = T} (fixed tournaments).

- H is a vector in this space. Its isotypic decomposition tells us
  which irreps "see" H.

For n=3,4,5: compute orbits, H per orbit, permutation character,
irrep multiplicities, and H's projection onto each isotypic component.

Author: kind-pasteur-2026-03-14-S68
"""

import math
from itertools import permutations, combinations
from collections import defaultdict, Counter
from fractions import Fraction


# ─────────────────── Tournament utilities ───────────────────

def binary_to_tournament(bits, n):
    """Convert integer bits to adjacency matrix. Bit i encodes edge (u,v) for u<v."""
    m = n * (n - 1) // 2
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def tournament_to_bits(A, n):
    """Convert adjacency matrix back to bits."""
    bits = 0
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j]:
                bits |= (1 << pos)
            pos += 1
    return bits


def count_ham_paths(A, n):
    """Count Hamiltonian paths using bitmask DP."""
    if n <= 1:
        return 1
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def apply_perm_to_tournament(A, sigma, n):
    """Apply permutation sigma to tournament A.
    (sigma.T)[i][j] = T[sigma_inv[i]][sigma_inv[j]]
    """
    sigma_inv = [0]*n
    for i in range(n):
        sigma_inv[sigma[i]] = i
    B = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            B[i][j] = A[sigma_inv[i]][sigma_inv[j]]
    return B


# ─────────────────── S_n character theory ───────────────────

def cycle_type(perm):
    """Return the cycle type of a permutation as a sorted tuple (descending)."""
    n = len(perm)
    visited = [False]*n
    cycles = []
    for i in range(n):
        if visited[i]:
            continue
        length = 0
        j = i
        while not visited[j]:
            visited[j] = True
            j = perm[j]
            length += 1
        cycles.append(length)
    return tuple(sorted(cycles, reverse=True))


def partitions(n):
    """Generate all partitions of n in descending order."""
    if n == 0:
        yield ()
        return
    def _partitions(n, max_part):
        if n == 0:
            yield ()
            return
        for part in range(min(n, max_part), 0, -1):
            for rest in _partitions(n - part, part):
                yield (part,) + rest
    yield from _partitions(n, n)


def conjugacy_classes(n):
    """Return conjugacy classes of S_n: list of (cycle_type, size, representative).
    cycle_type = partition of n (descending)."""
    classes = {}
    for perm in permutations(range(n)):
        ct = cycle_type(perm)
        if ct not in classes:
            classes[ct] = [perm, 0]
        classes[ct][1] += 1
    result = []
    for ct in sorted(classes.keys(), reverse=True):
        rep, size = classes[ct]
        result.append((ct, size, rep))
    return result


def character_table_sn(n):
    """Compute the character table of S_n using the Murnaghan-Nakayama rule.
    Returns: (list of partitions, list of conj_classes, 2D table chi[i][j])
    where chi[i][j] = character of irrep labeled by partition i, evaluated at conjugacy class j.
    """
    parts = list(partitions(n))
    conj_cls = conjugacy_classes(n)

    # ── Rim hook (border strip) removal ──
    # A rim hook of size r in partition lam is a connected skew shape
    # along the SE boundary with r cells and no 2x2 block.
    # Leg length = (number of rows it spans) - 1.
    #
    # Standard algorithm: convert partition to its "beta numbers" representation.
    # beta_i = lam_i + (num_parts - 1 - i) for each row i.
    # A rim hook of size r exists iff some beta_i - r is NOT already a beta number
    # and is >= 0. The new partition is obtained by replacing beta_i with beta_i - r.
    # The leg length = number of beta_j strictly between beta_i - r and beta_i.

    def remove_rim_hooks(lam, r):
        """Find all rim hooks of size r removable from partition lam.
        Returns list of (new_partition_tuple, leg_length)."""
        if not lam:
            return []
        k = len(lam)  # number of parts (including trailing zeros if any)
        # Compute beta numbers: beta_i = lam_i + (k - 1 - i)
        beta = [lam[i] + (k - 1 - i) for i in range(k)]
        beta_set = set(beta)

        results = []
        for i in range(k):
            new_b = beta[i] - r
            if new_b < 0:
                continue
            if new_b in beta_set:
                continue  # would create a duplicate beta number
            # Compute leg length = number of beta_j strictly between new_b and beta[i]
            leg = sum(1 for j in range(k) if j != i and new_b < beta[j] < beta[i])
            # Build new beta set
            new_beta = sorted([new_b if j == i else beta[j] for j in range(k)], reverse=True)
            # Convert back to partition
            new_lam = [new_beta[j] - (k - 1 - j) for j in range(k)]
            # Remove trailing zeros
            while new_lam and new_lam[-1] == 0:
                new_lam.pop()
            results.append((tuple(new_lam), leg))
        return results

    # Memoized MN character computation
    _mn_cache = {}

    def mn_char(lam, mu):
        """Compute chi^lam(mu) via Murnaghan-Nakayama rule.
        lam = partition (irrep label), mu = cycle type (conjugacy class label)."""
        key = (lam, mu)
        if key in _mn_cache:
            return _mn_cache[key]

        if not mu:
            result = 1 if (not lam or lam == ()) else 0
            _mn_cache[key] = result
            return result

        # Peel off the last (smallest) part of mu
        mu_list = list(mu)
        r = mu_list.pop()
        mu_rest = tuple(mu_list)

        total = 0
        for new_lam, leg in remove_rim_hooks(lam, r):
            sign = (-1) ** leg
            total += sign * mn_char(new_lam, mu_rest)

        _mn_cache[key] = total
        return total

    # Build character table
    chi = []
    ctypes = [c[0] for c in conj_cls]
    for lam in parts:
        row = []
        for ct in ctypes:
            row.append(mn_char(lam, ct))
        chi.append(row)

    return parts, conj_cls, chi


def verify_character_table(parts, conj_classes, chi, n):
    """Verify orthogonality relations of the character table."""
    n_fact = math.factorial(n)
    num_classes = len(conj_classes)

    # Row orthogonality: sum_g chi_i(g) * chi_j(g)^* = |G| delta_{ij}
    # = sum_c |c| * chi_i(c) * chi_j(c) where |c| is conjugacy class size
    print(f"\n  Verifying character table for S_{n}...")
    all_ok = True
    for i in range(len(parts)):
        for j in range(i, len(parts)):
            inner = sum(conj_classes[k][1] * chi[i][k] * chi[j][k]
                       for k in range(num_classes))
            expected = n_fact if i == j else 0
            if inner != expected:
                print(f"  FAIL: <chi_{parts[i]}, chi_{parts[j]}> = {inner}, expected {expected}")
                all_ok = False
    if all_ok:
        print(f"  Character table VERIFIED (row orthogonality)")

    # Print dimensions
    # Identity permutation has cycle type (1,1,...,1), which is LAST in reverse-sorted order
    id_idx = num_classes - 1  # index of identity conjugacy class
    dims = [chi[i][id_idx] for i in range(len(parts))]  # chi(identity) = dimension
    print(f"  Irrep dimensions: {dict(zip(parts, dims))}")
    print(f"  Sum of d^2 = {sum(d*d for d in dims)} (should be {n_fact})")


# ─────────────────── Main computation ───────────────────

def compute_fixed_tournaments(sigma, n):
    """Count tournaments T on [n] such that sigma.T = T."""
    m = n * (n - 1) // 2
    count = 0
    for bits in range(1 << m):
        A = binary_to_tournament(bits, n)
        B = apply_perm_to_tournament(A, sigma, n)
        if all(A[i][j] == B[i][j] for i in range(n) for j in range(n)):
            count += 1
    return count


def compute_orbits(n):
    """Compute orbits of S_n acting on tournaments.
    Returns list of (representative_bits, orbit_size, H_value)."""
    m = n * (n - 1) // 2
    total = 1 << m
    visited = set()
    orbits = []

    perms = list(permutations(range(n)))

    for bits in range(total):
        if bits in visited:
            continue
        A = binary_to_tournament(bits, n)
        H = count_ham_paths(A, n)
        orbit = set()
        for sigma in perms:
            B = apply_perm_to_tournament(A, sigma, n)
            b = tournament_to_bits(B, n)
            orbit.add(b)
        for b in orbit:
            visited.add(b)
        orbits.append((bits, len(orbit), H))

    return orbits


def compute_permutation_character(n, conj_classes):
    """For each conjugacy class, compute chi(sigma) = # fixed tournaments."""
    char = []
    for ct, size, rep in conj_classes:
        fixed = compute_fixed_tournaments(rep, n)
        char.append(fixed)
    return char


def decompose_character(perm_char, conj_classes, chi, n):
    """Decompose a character into irrep multiplicities.
    mult_i = (1/|G|) * sum_g perm_char(g) * chi_i(g)^*
    = (1/n!) * sum_c |c| * perm_char(c) * chi_i(c)
    """
    n_fact = math.factorial(n)
    mults = []
    for i in range(len(chi)):
        inner = sum(conj_classes[k][1] * perm_char[k] * chi[i][k]
                    for k in range(len(conj_classes)))
        mult = Fraction(inner, n_fact)
        mults.append(mult)
    return mults


def project_H_onto_irreps(n, orbits, conj_classes, chi, parts):
    """Project the function H onto each isotypic component.

    H is a function on tournaments. We compute its projection onto
    each irrep using:

    For each irrep lambda with character chi_lambda:
      ||proj_lambda(H)||^2 = (d_lambda / |G|^2) * sum_{g,h} chi_lambda(g^{-1}h) * H(g.T0) * H(h.T0)
      summed over all T0.

    But simpler: the norm-squared of H's projection onto the lambda-isotypic component is:
      (d_lambda / |G|) * sum_g chi_lambda(g) * <H, g.H>

    where <H, g.H> = sum_T H(T) * H(g.T) = sum_T H(T) * H(g.T).

    Since H(g.T) = H(T) (H is S_n-invariant!), we get:
      <H, g.H> = sum_T H(T)^2 for ALL g.

    So ||proj_lambda(H)||^2 = (d_lambda / |G|) * sum_T H(T)^2 * sum_g chi_lambda(g)

    Now sum_g chi_lambda(g) = |G| if lambda = trivial, 0 otherwise.

    WAIT: This means H lives ENTIRELY in the trivial representation!

    That's because H is S_n-INVARIANT: H(sigma.T) = H(T).
    So sigma.H = H for all sigma, meaning H is a fixed vector,
    which is exactly the trivial representation.

    But the INTERESTING question is about the PERMUTATION REPRESENTATION
    on tournaments, not about H as a single vector.

    Let me reframe: The space of functions {tournaments} -> R is a rep of S_n.
    H is a vector in this space that is FIXED by S_n (it's in the trivial isotypic component).
    The trivial rep appears with multiplicity = number of orbits.

    The more interesting analysis is:
    1. How does the permutation rep on tournaments decompose?
    2. What is the structure of the trivial component (= orbit-constant functions)?
    3. What is the "H landscape" on the orbit space?
    """
    pass  # We handle this in the main analysis


def compute_H_weighted_character(n, conj_classes):
    """Compute the H-weighted character: chi_H(sigma) = sum_T H(T) * [sigma.T = T].
    This is the character of the action weighted by H values on fixed points."""
    m = n * (n - 1) // 2
    total = 1 << m
    char = []

    for ct, size, sigma in conj_classes:
        weighted_sum = 0
        for bits in range(total):
            A = binary_to_tournament(bits, n)
            B = apply_perm_to_tournament(A, sigma, n)
            if all(A[i][j] == B[i][j] for i in range(n) for j in range(n)):
                H = count_ham_paths(A, n)
                weighted_sum += H
        char.append(weighted_sum)

    return char


def compute_H2_weighted_character(n, conj_classes):
    """Compute chi_{H^2}(sigma) = sum_T H(T)^2 * [sigma.T = T].
    Decomposing this gives the squared norm of H's projection, and since H is
    invariant, tells us about ||H||^2 restricted to the trivial subspace."""
    m = n * (n - 1) // 2
    total = 1 << m
    char = []

    for ct, size, sigma in conj_classes:
        weighted_sum = 0
        for bits in range(total):
            A = binary_to_tournament(bits, n)
            B = apply_perm_to_tournament(A, sigma, n)
            if all(A[i][j] == B[i][j] for i in range(n) for j in range(n)):
                H = count_ham_paths(A, n)
                weighted_sum += H * H
        char.append(weighted_sum)

    return char


def analyze_orbit_structure(n, orbits):
    """Analyze the structure of orbits: sizes, H values, etc."""
    print(f"\n{'='*70}")
    print(f"  ORBIT STRUCTURE for n = {n}")
    print(f"{'='*70}")
    print(f"  Number of labeled tournaments: {2**(n*(n-1)//2)}")
    print(f"  Number of orbits (non-isomorphic tournaments): {len(orbits)}")
    print()

    # Sort orbits by H value
    sorted_orbits = sorted(orbits, key=lambda x: x[2])

    print(f"  {'Orbit':>5} {'|orbit|':>8} {'H(T)':>8} {'|orbit|*H':>10} {'|orbit|*H^2':>12}")
    print(f"  {'-'*5:>5} {'-'*8:>8} {'-'*8:>8} {'-'*10:>10} {'-'*12:>12}")

    total_H = 0
    total_H2 = 0
    H_values = []
    for idx, (bits, size, H) in enumerate(sorted_orbits):
        print(f"  {idx+1:>5} {size:>8} {H:>8} {size*H:>10} {size*H*H:>12}")
        total_H += size * H
        total_H2 += size * H * H
        H_values.append(H)

    print(f"\n  Total sum H = {total_H}")
    n_fact = math.factorial(n)
    print(f"  Average H = {total_H} / {2**(n*(n-1)//2)} = {Fraction(total_H, 2**(n*(n-1)//2))}")
    print(f"  n! = {n_fact}")
    print(f"  Average H should be n!/2^(n-1) = {Fraction(n_fact, 2**(n-1))} "
          f"= {n_fact / 2**(n-1):.4f}")
    print(f"  Sum H^2 = {total_H2}")

    # H value distribution
    H_counter = Counter(H_values)
    print(f"\n  H value distribution across orbits:")
    for h_val in sorted(H_counter.keys()):
        print(f"    H = {h_val}: {H_counter[h_val]} orbit(s)")

    # H value distribution across labeled tournaments
    print(f"\n  H value distribution across labeled tournaments:")
    labeled_counter = defaultdict(int)
    for bits, size, H in orbits:
        labeled_counter[H] += size
    for h_val in sorted(labeled_counter.keys()):
        print(f"    H = {h_val}: {labeled_counter[h_val]} tournament(s)")

    return sorted_orbits


def main():
    print("="*70)
    print("  S_n REPRESENTATION THEORY OF TOURNAMENT SPACE")
    print("  Decomposition of H into irreducible representations")
    print("  kind-pasteur-2026-03-14-S68")
    print("="*70)

    for n in [3, 4, 5]:
        m = n * (n - 1) // 2
        total = 1 << m

        print(f"\n{'#'*70}")
        print(f"#  n = {n}: {total} labeled tournaments, |S_{n}| = {math.factorial(n)}")
        print(f"{'#'*70}")

        # Step 1: Compute character table of S_n
        print(f"\n  Step 1: Character table of S_{n}")
        parts, conj_classes, chi = character_table_sn(n)
        verify_character_table(parts, conj_classes, chi, n)

        # Print character table
        ctypes = [c[0] for c in conj_classes]
        csizes = [c[1] for c in conj_classes]
        print(f"\n  Character table of S_{n}:")
        header = f"  {'Partition':<15}"
        for ct in ctypes:
            header += f" {str(ct):>12}"
        print(header)
        print(f"  {'class size':<15}", end="")
        for sz in csizes:
            print(f" {sz:>12}", end="")
        print()
        print(f"  {'-'*15}", end="")
        for _ in ctypes:
            print(f" {'-'*12}", end="")
        print()
        for i, lam in enumerate(parts):
            row = f"  {str(lam):<15}"
            for j in range(len(ctypes)):
                row += f" {chi[i][j]:>12}"
            print(row)

        # Step 2: Compute orbits
        print(f"\n  Step 2: Tournament orbits under S_{n}")
        orbits = compute_orbits(n)
        sorted_orbits = analyze_orbit_structure(n, orbits)

        # Step 3: Permutation character
        print(f"\n  Step 3: Permutation character (# fixed tournaments per conj class)")
        perm_char = compute_permutation_character(n, conj_classes)
        print(f"  Conjugacy classes and fixed-tournament counts:")
        for k, (ct, size, _) in enumerate(conj_classes):
            print(f"    {str(ct):<15} (size {size:>3}): chi = {perm_char[k]}")

        # Step 4: Decompose permutation character
        print(f"\n  Step 4: Decomposition of permutation representation into irreps")
        mults = decompose_character(perm_char, conj_classes, chi, n)
        print(f"  Irrep multiplicities:")
        for i, lam in enumerate(parts):
            id_idx = len(conj_classes) - 1  # identity class is last
            dim = chi[i][id_idx]
            print(f"    {str(lam):<15}: multiplicity = {mults[i]}, "
                  f"dim = {dim}, contributes {int(mults[i])*dim} to total dim")
        id_idx = len(conj_classes) - 1
        total_dim = sum(int(mults[i]) * chi[i][id_idx] for i in range(len(parts)))
        print(f"  Total dimension: {total_dim} (should be {total})")

        # Step 5: H as a vector — it's S_n-invariant, so lives in trivial component
        print(f"\n  Step 5: H as an S_n-invariant vector")
        trivial_mult = int(mults[0])  # First partition is (n), the trivial rep
        print(f"  H is invariant under S_{n}, so H lives in the trivial isotypic component.")
        print(f"  Trivial rep multiplicity = {trivial_mult} = number of orbits = {len(orbits)}")
        print(f"  The trivial isotypic component has dimension {trivial_mult}.")
        print(f"  H is ONE specific vector in this {trivial_mult}-dimensional space.")

        # Step 6: H-weighted character
        print(f"\n  Step 6: H-weighted permutation character")
        print(f"  chi_H(sigma) = sum_{{T: sigma.T=T}} H(T)")
        H_char = compute_H_weighted_character(n, conj_classes)
        for k, (ct, size, _) in enumerate(conj_classes):
            print(f"    {str(ct):<15} (size {size:>3}): chi_H = {H_char[k]}")

        # Decompose H-weighted character
        H_mults = decompose_character(H_char, conj_classes, chi, n)
        print(f"\n  Decomposition of H-weighted character:")
        for i, lam in enumerate(parts):
            if H_mults[i] != 0:
                print(f"    {str(lam):<15}: coefficient = {H_mults[i]}")

        # Step 7: H^2-weighted character (gives ||H||^2 in each component)
        print(f"\n  Step 7: H^2-weighted permutation character")
        print(f"  chi_{{H^2}}(sigma) = sum_{{T: sigma.T=T}} H(T)^2")
        H2_char = compute_H2_weighted_character(n, conj_classes)
        for k, (ct, size, _) in enumerate(conj_classes):
            print(f"    {str(ct):<15} (size {size:>3}): chi_{{H^2}} = {H2_char[k]}")

        H2_mults = decompose_character(H2_char, conj_classes, chi, n)
        print(f"\n  Decomposition of H^2-weighted character:")
        for i, lam in enumerate(parts):
            if H2_mults[i] != 0:
                print(f"    {str(lam):<15}: coefficient = {H2_mults[i]}")

        # Step 8: Alternative analysis — the MARK MATRIX approach
        # For each orbit O, define indicator function 1_O (= 1 on orbit, 0 elsewhere)
        # H = sum_O H(O) * 1_O
        # The 1_O form a basis for the trivial isotypic component
        # The "norm" of H in this basis is sum_O H(O)^2 * |O| / n!
        # (weighted by orbit stabilizer size)
        print(f"\n  Step 8: H in the orbit basis")
        print(f"  H = sum_O H(O) * 1_O where 1_O is the indicator of orbit O")
        print(f"  Norm^2 in L^2(tournaments) = sum_O |O| * H(O)^2")
        norm2 = sum(size * H * H for _, size, H in orbits)
        avg = Fraction(sum(size * H for _, size, H in orbits), total)
        var = Fraction(norm2, total) - avg * avg
        print(f"  ||H||^2 = {norm2}")
        print(f"  Mean H = {avg} = {float(avg):.6f}")
        print(f"  Var(H) = {var} = {float(var):.6f}")

        # Step 9: The REFINED question — representation on H-level sets
        print(f"\n  Step 9: Representation theory of H-level sets")
        print(f"  For each value h, the set {{T : H(T) = h}} is S_n-stable.")
        print(f"  Its indicator function is a vector in the permutation rep.")
        H_levels = defaultdict(list)
        for bits, size, H in orbits:
            H_levels[H].append((bits, size))

        for h_val in sorted(H_levels.keys()):
            level_orbits = H_levels[h_val]
            level_size = sum(size for _, size in level_orbits)
            # The indicator of this level set is a character too
            # Its character: chi_h(sigma) = #{T with H(T)=h and sigma.T=T}
            level_char = []
            for ct, size, sigma in conj_classes:
                count = 0
                for bits in range(total):
                    A = binary_to_tournament(bits, n)
                    B = apply_perm_to_tournament(A, sigma, n)
                    if all(A[i][j] == B[i][j] for i in range(n) for j in range(n)):
                        if count_ham_paths(A, n) == h_val:
                            count += 1
                level_char.append(count)

            level_mults = decompose_character(level_char, conj_classes, chi, n)
            nonzero = [(parts[i], int(level_mults[i])) for i in range(len(parts)) if level_mults[i] != 0]
            print(f"    H = {h_val} ({level_size} tournaments, {len(level_orbits)} orbits): "
                  f"irreps = {nonzero}")

        # Step 10: Burnside's lemma verification
        print(f"\n  Step 10: Burnside's lemma verification")
        burnside = Fraction(sum(csz * perm_char[k] for k, (_, csz, _) in enumerate(conj_classes)),
                           math.factorial(n))
        print(f"  Number of orbits by Burnside = {burnside} (should be {len(orbits)})")

        # Step 11: The CYCLE INDEX of S_n acting on edges
        # Each permutation sigma acts on the m = C(n,2) edges.
        # sigma sends edge {i,j} to {sigma(i),sigma(j)}.
        # The cycle index on edges determines the tournament enumeration.
        print(f"\n  Step 11: Edge action analysis")
        print(f"  S_{n} acts on {m} edges. The cycle structure on edges determines fixed tournaments.")
        for k, (ct, size, sigma) in enumerate(conj_classes):
            # Compute cycle structure of sigma on edges
            edge_list = []
            for i in range(n):
                for j in range(i+1, n):
                    edge_list.append((i, j))

            # sigma acts on edges: {i,j} -> {sigma(i), sigma(j)}
            edge_perm = []
            for i, j in edge_list:
                si, sj = sigma[i], sigma[j]
                target = (min(si, sj), max(si, sj))
                target_idx = edge_list.index(target)
                edge_perm.append(target_idx)

            edge_ct = cycle_type(edge_perm)
            # Number of fixed tournaments = 2^(# cycles that DON'T pair with reverse)
            # Actually for tournaments, each edge has an orientation.
            # sigma fixes T iff for each edge cycle, the orientation is consistent.
            # For a cycle of length L on edges: need T[sigma^L(i)][sigma^L(j)] = T[i][j]
            # which is automatic. But also sigma reversal: if the cycle maps (i,j) to
            # (j,i) at some point, it forces T[i][j] and T[j][i], contradiction.
            # For tournaments: fixed count = 2^(# edge orbits under <sigma>)
            # where we treat ordered pairs, and each orbit is either a "directed" orbit
            # (stays consistently oriented) or gets canceled.

            # Actually: sigma acts on ORIENTED edges (i->j) by sending i->j to sigma(i)->sigma(j).
            # The unoriented edge {i,j} splits into two oriented edges.
            # For a tournament, exactly one of each pair is chosen.
            # An edge cycle of length L under sigma on unoriented edges:
            # Check if sigma preserves or reverses orientation along the cycle.
            # If L is odd and the cycle is a "reversing" cycle, no fixed tournament has this.
            # If the cycle never reverses, we have 2 choices (all one way or all the other).

            n_edge_cycles = len(edge_ct)
            print(f"    {str(ct):<15}: edge cycle type = {edge_ct}, "
                  f"{n_edge_cycles} cycles, fixed = {perm_char[k]} = 2^{n_edge_cycles} ? "
                  f"{'YES' if perm_char[k] == 2**n_edge_cycles else 'NO (orientation constraint)'}")

    # ─────────── Additional analysis: the REPRESENTATION on orbit space ───────────
    print(f"\n{'#'*70}")
    print(f"#  SUMMARY: Representation-theoretic structure of H")
    print(f"{'#'*70}")

    print(f"""
  KEY FINDINGS:

  1. H(T) is S_n-INVARIANT: H(sigma.T) = H(T) for all sigma in S_n.
     Therefore H, viewed as a vector in Fun(Tournaments, R), lives
     ENTIRELY in the trivial isotypic component.

  2. The trivial isotypic component has dimension = number of orbits
     (non-isomorphic tournaments). This is the space of "isomorphism
     class functions."

  3. The PERMUTATION REPRESENTATION on tournaments decomposes into
     irreps with specific multiplicities shown above. The multiplicity
     of the trivial rep equals the number of orbits (Burnside).

  4. The MORE INTERESTING question: the LEVEL SETS {{T : H(T) = h}} are
     S_n-stable subsets. Their indicator functions span the trivial
     component, and we computed their decomposition.

  5. H is NOT equidistributed on orbits — different orbits have different
     H values. The variance of H across orbits measures how "spread out"
     the Hamiltonian path count is.

  6. For each H-level set, the representation decomposes into multiple
     irreps, showing how the symmetry group "sees" the tournament
     structure within each level.
""")

    # ─────────── Bonus: Representation on ORBIT-PAIR correlations ───────────
    # The interesting representation-theoretic content is in the
    # permutation character and its decomposition. Let's do a deeper analysis
    # for n=4 specifically.
    print(f"\n{'#'*70}")
    print(f"#  DEEP DIVE: n=4 — Orbit adjacency and representation structure")
    print(f"{'#'*70}")

    n = 4
    m = n * (n - 1) // 2
    total = 1 << m
    orbits_4 = compute_orbits(4)

    # For each pair of orbits, compute the number of single-arc-flips between them
    print(f"\n  Arc-flip graph on orbits:")
    print(f"  Two orbits are adjacent if a single arc flip connects them.")

    orbit_map = {}  # bits -> orbit index
    for idx, (bits, size, H) in enumerate(orbits_4):
        A = binary_to_tournament(bits, n)
        perms = list(permutations(range(n)))
        for sigma in perms:
            B = apply_perm_to_tournament(A, sigma, n)
            b = tournament_to_bits(B, n)
            orbit_map[b] = idx

    # Build adjacency between orbits via arc flips
    adj = defaultdict(lambda: defaultdict(int))
    for bits in range(total):
        oi = orbit_map[bits]
        A = binary_to_tournament(bits, n)
        # Try flipping each edge
        for pos in range(m):
            flipped = bits ^ (1 << pos)
            oj = orbit_map[flipped]
            adj[oi][oj] += 1

    print(f"\n  Orbit adjacency matrix (# arc flips between orbits):")
    n_orbits = len(orbits_4)
    sorted_idx = list(range(n_orbits))
    sorted_idx.sort(key=lambda i: orbits_4[i][2])

    header = f"  {'Orbit':<12}"
    for j in sorted_idx:
        header += f" O{j}(H={orbits_4[j][2]})"
    print(header)
    for i in sorted_idx:
        row = f"  O{i}(H={orbits_4[i][2]},|{orbits_4[i][1]}|)"
        for j in sorted_idx:
            row += f" {adj[i][j]:>10}"
        print(row)

    # ─────────── Score sequence analysis per orbit ───────────
    print(f"\n  Score sequences per orbit (n=4):")
    for idx, (bits, size, H) in enumerate(orbits_4):
        A = binary_to_tournament(bits, n)
        scores = sorted([sum(A[i]) for i in range(n)])
        print(f"    Orbit {idx}: H={H}, |orbit|={size}, score={scores}")

    # ─────────── Automorphism group analysis ───────────
    print(f"\n  Automorphism groups per orbit (n=4):")
    for idx, (bits, size, H) in enumerate(orbits_4):
        A = binary_to_tournament(bits, n)
        aut_size = math.factorial(n) // size
        # Find actual automorphisms
        auts = []
        for sigma in permutations(range(n)):
            B = apply_perm_to_tournament(A, list(sigma), n)
            if all(A[i][j] == B[i][j] for i in range(n) for j in range(n)):
                auts.append(sigma)
        aut_types = Counter(cycle_type(list(a)) for a in auts)
        print(f"    Orbit {idx}: H={H}, |Aut|={len(auts)}, "
              f"Aut cycle types: {dict(aut_types)}")


if __name__ == "__main__":
    main()
