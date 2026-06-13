"""
paley_ocf_mechanism.py — WHY Paley maximizes odd cycle count (α₁)

Key observation from dihedral_hp_analysis.out:
  Paley T_7:     α₁=80, α₂=7,  H = 1 + 2*80 + 4*7 = 189
  Non-Paley:     α₁=59, α₂=14, H = 1 + 2*59 + 4*14 = 175

The paradox: Paley has MORE cycles but FEWER disjoint pairs.
Yet it wins because 2*80 = 160 > 2*59+4*14-4*7 = 146.

Questions:
1. What cycle LENGTH distribution gives Paley its edge? (length 3, 5, or 7?)
2. Is the cycle excess explained by the automorphism group?
3. Can we prove a "cycle count inequality" from the dihedral structure?
4. What is the EXACT relationship between eigenvalues and cycle counts?

The spectral connection: for a circulant with eigenvalues λ_k,
the number of directed cycles of length ℓ is:
  C_ℓ = (1/ℓ) Σ_k λ_k^ℓ    (since Tr(A^ℓ) counts closed walks, /ℓ for cycles)

For odd ℓ in a tournament, all closed walks of length ℓ ARE directed cycles
(since tournament has no 2-cycles).

Actually: Tr(A^ℓ) counts closed walks, which include non-simple ones.
For length 3: all closed walks are simple → C_3 = Tr(A^3)/3.
For length 5: need to subtract non-simple...

Let's use the eigenvalue formula directly and compare.

Author: opus-2026-03-12-S60
"""
import sys
import time
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def circulant_eigenvalues(n, S):
    """Eigenvalues of circulant tournament on Z_n with connection set S."""
    omega = np.exp(2j * np.pi / n)
    return [sum(omega ** (k * s) for s in S) for k in range(n)]


def closed_walks(n, S, length):
    """Tr(A^ℓ) for circulant tournament = Σ_k λ_k^ℓ."""
    evals = circulant_eigenvalues(n, S)
    return sum(ev ** length for ev in evals)


def count_directed_3cycles(n, S):
    """Exact count of directed 3-cycles in circulant tournament on Z_n."""
    S_set = set(S)
    count = 0
    for a in S_set:
        for b in S_set:
            # edge 0→a, a→a+b, need a+b→0, i.e., (n-(a+b))%n ∈ S
            if (n - (a + b) % n) % n in S_set:
                count += 1
    # Each 3-cycle is counted 3 times (3 starting vertices in the orbit, but we're fixing vertex 0)
    # Actually: for circulant, each 3-cycle through vertex 0 has a specific structure
    # A 3-cycle is (0, a, a+b) with a∈S, b∈S, (n-a-b)%n ∈ S
    # The total number of directed 3-cycles is n * count / 3
    return n * count // 3


def count_all_directed_cycles(n, S, max_len=None):
    """Count all directed cycles in circulant tournament by DFS from vertex 0.

    By circulant symmetry, total = n * cycles_through_0.
    But each cycle of length ℓ passes through vertex 0 exactly ℓ/n * n = ℓ times...
    No: each cycle of length ℓ visits ℓ vertices, so passes through 0 with probability ℓ/n.
    So total cycles of length ℓ = n/ℓ * (cycles through 0 of length ℓ).

    Actually: for circulant, each cycle through 0 of length ℓ has exactly ℓ rotated copies
    (which form the same cycle). And there are n copies from translation. So:
    total distinct cycles of length ℓ = (n * cycles_from_0) / ℓ.
    """
    if max_len is None:
        max_len = n
    S_set = set(S)

    # Count cycles starting at vertex 0
    cycles_by_len = Counter()

    def dfs(current, visited, length):
        if length > max_len:
            return
        for s in S_set:
            next_v = (current + s) % n
            if next_v == 0 and length >= 3:
                # Found a cycle back to 0
                cycles_by_len[length] += 1
            elif next_v not in visited and next_v != 0:
                visited.add(next_v)
                dfs(next_v, visited, length + 1)
                visited.remove(next_v)

    dfs(0, {0}, 1)

    # Convert to total cycle counts
    total_by_len = {}
    for ell, c in cycles_by_len.items():
        total_by_len[ell] = n * c // ell

    return total_by_len


def independence_number(n, S):
    """Compute α₁ (odd cycles), α₂ (disjoint pairs), etc. for OCF.

    Actually, we need the FULL odd cycle independence polynomial.
    An odd cycle C is a subset of vertices forming a directed cycle of odd length.
    Two cycles are "independent" if they share no vertex.

    I(Ω,x) = 1 + α₁*x + α₂*x² + ... where α_k = number of collections of k
    pairwise vertex-disjoint odd cycles.

    H(T) = I(Ω, 2) = Σ_k α_k * 2^k.
    """
    # Find all odd directed cycles
    cycles = find_all_odd_cycles(n, S)

    # Build conflict graph: two cycles conflict if they share a vertex
    vertex_sets = [frozenset(c) for c in cycles]
    m = len(cycles)

    # Compute independence polynomial I(x) of the conflict graph
    # For small m this is tractable
    if m > 30:
        # Use approximation: just compute α₁, α₂, α₃
        alpha = [0] * 4
        alpha[1] = m
        # α₂ = number of non-conflicting pairs
        for i in range(m):
            for j in range(i+1, m):
                if vertex_sets[i].isdisjoint(vertex_sets[j]):
                    alpha[2] += 1
        # α₃ = number of non-conflicting triples
        for i in range(m):
            for j in range(i+1, m):
                if not vertex_sets[i].isdisjoint(vertex_sets[j]):
                    continue
                for k in range(j+1, m):
                    if vertex_sets[i].isdisjoint(vertex_sets[k]) and vertex_sets[j].isdisjoint(vertex_sets[k]):
                        alpha[3] += 1
        return {k: v for k, v in enumerate(alpha)}

    # Exact: enumerate all independent sets
    alpha = Counter()

    def find_indep_sets2(idx, size, used_vertices):
        alpha[size] += 1
        for i in range(idx, m):
            if vertex_sets[i].isdisjoint(used_vertices):
                find_indep_sets2(i + 1, size + 1, used_vertices | vertex_sets[i])

    find_indep_sets2(0, 0, frozenset())

    return dict(sorted(alpha.items()))


def find_all_odd_cycles(n, S):
    """Find all directed odd cycles in circulant tournament on Z_n."""
    S_set = set(S)
    cycles = set()

    def dfs(start, current, visited_tuple, length):
        if length > n:
            return
        for s in S_set:
            next_v = (current + s) % n
            if next_v == start and length >= 3 and length % 2 == 1:
                # Normalize: use lexicographically smallest rotation
                cycle = visited_tuple
                min_rot = min(cycle[i:] + cycle[:i] for i in range(len(cycle)))
                cycles.add(min_rot)
            elif next_v not in set(visited_tuple) and next_v > start:
                # Only allow next_v > start to avoid duplicates from different starts
                dfs(start, next_v, visited_tuple + (next_v,), length + 1)

    for start in range(n):
        dfs(start, start, (start,), 1)

    return [list(c) for c in cycles]


def main():
    print("PALEY OCF MECHANISM: WHY MORE ODD CYCLES?")
    print("=" * 70)

    for n in [3, 5, 7]:
        print(f"\n{'='*70}")
        print(f"n = {n}")
        print(f"{'='*70}")

        k = (n - 1) // 2
        pairs = [(d, n - d) for d in range(1, k + 1)]

        # All circulant tournaments
        tournaments = []
        for bits in range(1 << k):
            S = set()
            for i in range(k):
                if bits & (1 << i):
                    S.add(pairs[i][0])
                else:
                    S.add(pairs[i][1])
            tournaments.append(frozenset(S))

        for S in tournaments:
            evals = circulant_eigenvalues(n, S)

            # Cycle counts by length
            cycle_counts = count_all_directed_cycles(n, S)
            total_odd_cycles = sum(c for l, c in cycle_counts.items() if l % 2 == 1)

            # Trace formula check: Tr(A^3) = Σ λ_k^3, should equal 3 * C_3
            tr3 = sum(ev**3 for ev in evals)

            # Check if Paley
            is_paley = False
            if all(pow(x, 2, n) != 0 for x in S):
                qr = set(pow(x, 2, n) for x in range(1, n))
                is_paley = (set(S) == qr or set(S) == set(range(1, n)) - qr)

            label = "PALEY" if is_paley else "other"

            print(f"\n  S={sorted(S)} [{label}]")
            print(f"    Cycle counts by length: {dict(sorted(cycle_counts.items()))}")
            print(f"    Total odd cycles: {total_odd_cycles}")
            print(f"    Tr(A^3)/3 = {tr3.real/3:.1f} (should = C_3 = {cycle_counts.get(3, 0)})")

            if n <= 7:
                # Eigenvalue powers
                for ell in [3, 5, 7]:
                    if ell > n:
                        break
                    tr_ell = sum(ev**ell for ev in evals)
                    print(f"    Tr(A^{ell}) = {tr_ell.real:.1f}, Σ|λ|^{ell} = {sum(abs(ev)**ell for ev in evals):.4f}")

            if n <= 7:
                # Full independence polynomial
                indep = independence_number(n, S)
                H = sum(v * (2**k) for k, v in indep.items())
                print(f"    Independence polynomial: {indep}")
                print(f"    H = I(Ω,2) = {H}")

    # Key analysis: eigenvalue formula for cycle counts
    print(f"\n{'='*70}")
    print(f"EIGENVALUE FORMULA FOR CYCLE COUNTS")
    print(f"{'='*70}")

    n = 7
    # Paley: S = {1,2,4}
    S_paley = frozenset([1, 2, 4])
    S_other = frozenset([1, 2, 3])  # interval

    evals_p = circulant_eigenvalues(n, S_paley)
    evals_o = circulant_eigenvalues(n, S_other)

    print(f"\n  Paley S={set(S_paley)}: eigenvalues = {[f'{e:.3f}' for e in evals_p]}")
    print(f"  Other S={set(S_other)}: eigenvalues = {[f'{e:.3f}' for e in evals_o]}")

    # For circulant: Tr(A^ℓ) = Σ_k λ_k^ℓ
    # For Paley with flat |λ_k| = r for k≥1:
    #   Tr(A^ℓ) = λ_0^ℓ + Σ_{k=1}^{n-1} λ_k^ℓ
    #   = ((n-1)/2)^ℓ + Σ_{k=1}^{n-1} r^ℓ e^{iℓθ_k}
    #   where θ_k = arg(λ_k)

    print(f"\n  Power sums comparison:")
    for ell in range(1, 8):
        tr_p = sum(ev**ell for ev in evals_p)
        tr_o = sum(ev**ell for ev in evals_o)
        print(f"    ℓ={ell}: Tr_Paley={tr_p.real:>10.3f}, Tr_Other={tr_o.real:>10.3f}, diff={tr_p.real-tr_o.real:>8.3f}")

    # Key identity: for odd ℓ, Tr(A^ℓ) counts closed directed walks of length ℓ
    # For ℓ=n=7, this includes Hamiltonian cycles.
    # But also non-Hamiltonian walks.

    # The permanent connection: H(T) = sum of Hamiltonian paths
    # = sum over permutations σ of ∏ A[i][σ(i)] where σ has no fixed points... no
    # Actually H(T) = perm(M) where M is the path adjacency matrix... not quite.

    # H(T) = I(Ω(T), 2) where Ω(T) is the odd-cycle complex.
    # The connection to eigenvalues is INDIRECT, through the cycle structure.

    # But: for CIRCULANT tournaments, there's a more direct connection.
    # All eigenvalues have Re(λ_k) = -1/2 for k≥1 (universal for tournaments on Z_p
    # with p prime and |S| = (p-1)/2).

    print(f"\n  Real parts of eigenvalues:")
    for k in range(n):
        print(f"    k={k}: Paley Re={evals_p[k].real:.6f}, Other Re={evals_o[k].real:.6f}")

    # Verification: for any circulant tournament on Z_p,
    # Σ_{s∈S} cos(2πks/p) = -1/2 for k≠0 (since S ∪ (-S) = Z_p\{0})
    print(f"\n  Verification: Re(λ_k) = -1/2 for all k≥1 and ALL circulant tournaments on Z_p?")
    for S in tournaments:
        evals = circulant_eigenvalues(n, S)
        for k in range(1, n):
            if abs(evals[k].real + 0.5) > 1e-10:
                print(f"    FAIL: S={sorted(S)}, k={k}, Re(λ_k)={evals[k].real}")
                break
        else:
            continue
        break
    else:
        print(f"    CONFIRMED: All circulant tournaments on Z_7 have Re(λ_k)=-1/2 for k≥1")

    # So the ONLY difference between circulants is the imaginary parts / phases of eigenvalues!
    # Re(λ_k) = -1/2 universally, |λ_k| varies.
    # Since Re = -1/2 and |λ_k|² = 1/4 + Im(λ_k)², we have |λ_k| = √(1/4 + Im²)
    # Paley: all |λ_k| = √((p+1)/4), so Im(λ_k)² = p/4 for all k

    print(f"\n  Imaginary parts:")
    for k in range(1, n):
        print(f"    k={k}: Paley Im={evals_p[k].imag:.6f} (|Im|={abs(evals_p[k].imag):.6f}), "
              f"Other Im={evals_o[k].imag:.6f} (|Im|={abs(evals_o[k].imag):.6f})")

    # The UNIVERSALITY of Re(λ_k) = -1/2 is key!
    # This means: λ_k = -1/2 + i*y_k where y_k are the only free parameters.
    # Constraint: Σ y_k² = p(p-1)/8 (from Parseval/sum-of-squares)
    # This is the "spectral simplex" mentioned in kind-pasteur's analysis.

    print(f"\n  Sum of Im² (should be p(p-1)/8 = {n*(n-1)/8}):")
    print(f"    Paley: Σ Im² = {sum(evals_p[k].imag**2 for k in range(1,n)):.6f}")
    print(f"    Other: Σ Im² = {sum(evals_o[k].imag**2 for k in range(1,n)):.6f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
