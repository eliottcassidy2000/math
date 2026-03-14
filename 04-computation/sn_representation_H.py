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
    Returns: (list of partitions, list of cycle_types, 2D table chi[i][j])
    where chi[i][j] = character of irrep labeled by partition i, evaluated at conjugacy class j.
    """
    parts = list(partitions(n))
    conj_classes = conjugacy_classes(n)
    ctypes = [c[0] for c in conj_classes]

    # Murnaghan-Nakayama rule
    def mn_character(lam, mu):
        """Character of irrep lambda at conjugacy class of type mu."""
        if sum(mu) == 0:
            return 1
        # Remove first part of mu, apply MN recursion
        return _mn_recurse(list(lam), list(mu))

    def _mn_recurse(lam, mu):
        """Recursive Murnaghan-Nakayama."""
        if not mu:
            if all(p == 0 for p in lam):
                return 1
            return 0
        # Take the last part of mu
        r = mu[-1]
        mu_rest = mu[:-1]
        total = 0
        # Find all border strips of size r in the Young diagram of lam
        for strip in border_strips(lam, r):
            new_lam, height = strip
            sign = (-1)**height
            total += sign * _mn_recurse(new_lam, mu_rest)
        return total

    def border_strips(lam, r):
        """Find all border strips (rim hooks) of size r that can be removed from lambda.
        Returns list of (new_lambda, height) where height = #rows - 1."""
        # Represent lambda as a list of parts
        n_rows = len(lam)
        # Pad with zeros
        lam_padded = list(lam) + [0]*r

        results = []
        # A border strip is a connected skew shape on the border
        # Try all possible removals
        _find_border_strips(lam_padded, r, 0, [], results)
        return results

    def _find_border_strips(lam, r, start_row, removed, results):
        """Find border strips by removing r cells from the border."""
        if r == 0:
            # Check valid partition and connectivity
            new_lam = list(lam)
            # Remove trailing zeros
            while new_lam and new_lam[-1] == 0:
                new_lam.pop()
            if not new_lam:
                new_lam = []
            # Check it's still a valid partition
            for i in range(len(new_lam)-1):
                if new_lam[i] < new_lam[i+1]:
                    return
            # Height = number of distinct rows touched - 1
            rows_touched = set(row for row, _ in removed)
            height = len(rows_touched) - 1
            results.append((tuple(new_lam), height))
            return
        # This approach is too slow for large n. Use a cleaner method.
        # ... Actually let's use the standard approach for border strip removal.
        pass

    # Actually, let me use a cleaner implementation of the MN rule.
    # Use the standard formulation via rim hooks.

    def remove_rim_hooks(partition, size):
        """Find all rim hooks of given size removable from partition.
        Returns list of (new_partition, leg_length)."""
        lam = list(partition)
        n_rows = len(lam)
        if n_rows == 0:
            return []

        results = []
        # For each possible bottom-right cell of the rim hook
        # A rim hook starts at some cell on the boundary and goes up-left
        # The boundary consists of cells (i, lam[i]-1) for each row i, and
        # cells (n_rows, j) for j = 0, ..., lam[-1]-1 (below the diagram)

        # Better approach: try removing a rim hook ending in each row
        # A rim hook of size r ending in row `end_row` and starting in row `start_row`
        # removes cells from the rightmost of rows start_row through end_row.

        # Standard method: for each row i (bottom of hook), try hooks going up
        for end_row in range(n_rows):
            # The hook must include the rightmost cell of row end_row
            # Going up from end_row, we remove cells from right side
            cells_remaining = size
            new_lam = list(lam)
            valid = True
            start_row = end_row

            for row in range(end_row, -1, -1):
                if cells_remaining <= 0:
                    break
                # How many cells can we remove from this row?
                # Must keep at least lam[row+1] cells (if row+1 exists) to maintain partition shape
                min_width = lam[row+1] if row+1 < n_rows else 0
                # But for the hook to be connected, we must remove from the right
                # and the removal must be contiguous on the border

                # For a rim hook: remove all cells in this row from column `min_width` to `lam[row]-1`
                # that are on the rim, up to cells_remaining
                available = lam[row] - min_width
                if available <= 0:
                    valid = False
                    break
                remove_here = min(available, cells_remaining)
                new_lam[row] = lam[row] - remove_here
                cells_remaining -= remove_here
                start_row = row

                # Check rim hook connectivity: after removing from this row,
                # the new row width must equal min_width (the row below's width)
                # unless this is the bottom row of the hook
                if row < end_row and new_lam[row] != (new_lam[row+1] if row+1 < len(new_lam) else 0):
                    # Not a valid rim hook (gap in connectivity)
                    # Actually, for a rim hook, each row's new width should equal
                    # the ORIGINAL width of the row below (for internal rows)
                    pass

            if cells_remaining > 0:
                valid = False

            if valid and cells_remaining == 0:
                # Verify it's a valid partition
                ok = True
                for row in range(len(new_lam)-1):
                    if new_lam[row] < new_lam[row+1]:
                        ok = False
                        break
                if ok:
                    # Remove trailing zeros
                    while new_lam and new_lam[-1] == 0:
                        new_lam.pop()
                    leg = end_row - start_row  # leg length = height
                    result_tuple = (tuple(new_lam), leg)
                    if result_tuple not in results:
                        results.append(result_tuple)

        return results

    def mn_char(partition, cycle_type_mu):
        """Compute character using Murnaghan-Nakayama rule recursively."""
        if not cycle_type_mu:
            if not partition or all(p == 0 for p in partition):
                return 1
            return 0

        mu = list(cycle_type_mu)
        r = mu.pop()  # Remove last (smallest) cycle length
        total = 0
        for new_lam, leg in remove_rim_hooks(partition, r):
            sign = (-1)**leg
            total += sign * mn_char(new_lam, tuple(mu))
        return total

    # Build character table
    chi = []
    for lam in parts:
        row = []
        for ct in ctypes:
            row.append(mn_char(lam, ct))
        chi.append(row)

    return parts, conj_classes, chi


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
    dims = [chi[i][0] for i in range(len(parts))]  # chi(identity) = dimension
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
            dim = chi[i][0]
            print(f"    {str(lam):<15}: multiplicity = {mults[i]}, "
                  f"dim = {dim}, contributes {int(mults[i])*dim} to total dim")
        total_dim = sum(int(mults[i]) * chi[i][0] for i in range(len(parts)))
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
