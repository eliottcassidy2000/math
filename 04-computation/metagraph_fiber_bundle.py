"""
metagraph_fiber_bundle.py — opus-2026-04-04-S12

BEYOND TRANSFER MATRICES: The metagraph G_n as a FIBER BUNDLE over G_{n-1}.

A tournament on n vertices = (tournament on {1,...,n-1}) + (extension σ ∈ {0,1}^{n-1}).
σ_i = 1 means vertex n beats vertex i.

Two extensions (T, σ₁) and (T, σ₂) give isomorphic tournaments iff
σ₂ = α(σ₁) for some automorphism α of T acting on {1,...,n-1}.

More generally: (T₁, σ₁) ≅ (T₂, σ₂) iff there exists a permutation π
of {1,...,n} with π(n) = n, T₂ = π(T₁), σ₂ = π(σ₁). This requires T₁ ≅ T₂
AND σ₂ related to σ₁ by the isomorphism.

APPROACH:
1. For each iso class C at n-1, compute Aut(C) acting on {0,1}^{n-1}
2. The orbits of Aut(C) on signatures = the "extensions" of C
3. Each orbit gives one iso class at n (with parent C)
4. The metagraph edges decompose into inner (change C) and boundary (change σ)

THIS IS THE RECURSIVE GENERATION without full 2^m enumeration.
"""
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import make_tiles, tiling_to_tournament, count_hamiltonian_paths

def tournament_canonical(T, n):
    """Canonical form = lexicographically smallest upper triangle under all permutations."""
    best = None
    best_perm = None
    for perm in permutations(range(n)):
        key = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(i+1, n))
        if best is None or key < best:
            best = key
            best_perm = perm
    return best, best_perm

def compute_aut_group(T, n):
    """Compute the automorphism group of tournament T.
    Returns list of permutations (as tuples)."""
    auts = []
    for perm in permutations(range(n)):
        if all(T[perm[i]][perm[j]] == T[i][j] for i in range(n) for j in range(n)):
            auts.append(perm)
    return auts

def extend_tournament(T_sub, n_sub, sigma):
    """Extend tournament T_sub on n_sub vertices by adding vertex n_sub
    with extension signature sigma ∈ {0,1}^{n_sub}.
    sigma[i] = 1 means new vertex beats vertex i."""
    n = n_sub + 1
    T = np.zeros((n, n), dtype=int)
    # Copy subtournament
    T[:n_sub, :n_sub] = T_sub
    # Add new vertex n_sub (0-indexed)
    for i in range(n_sub):
        if sigma[i]:
            T[n_sub][i] = 1  # new vertex beats i
        else:
            T[i][n_sub] = 1  # i beats new vertex
    return T

def sigma_orbit(sigma, aut_group, n_sub):
    """Compute the orbit of sigma under the automorphism group.
    Each aut permutes the positions of sigma."""
    orbit = set()
    for perm in aut_group:
        # Apply perm to sigma: new_sigma[perm[i]] = sigma[i] → new_sigma[j] = sigma[perm^{-1}[j]]
        # Actually: perm relabels vertices. If vertex i maps to perm[i],
        # then sigma[i] (does new vertex beat vertex i?) becomes:
        # in the relabeled tournament, does new vertex beat vertex perm[i]?
        # = sigma_new[perm[i]] = sigma[i]
        # So sigma_new[j] = sigma[perm^{-1}[j]]
        inv = [0]*n_sub
        for i in range(n_sub):
            inv[perm[i]] = i
        new_sigma = tuple(sigma[inv[j]] for j in range(n_sub))
        orbit.add(new_sigma)
    return frozenset(orbit)

# ==============================================================
# MAIN: Build the fiber bundle at n = 5 (from n-1 = 4)
# ==============================================================
def fiber_bundle(n_target):
    """Build G_{n_target} from G_{n_target-1} via the fiber bundle."""
    n_sub = n_target - 1
    print(f"\n{'='*60}")
    print(f"FIBER BUNDLE: G_{n_target} from G_{n_sub}")
    print(f"{'='*60}")

    # Step 1: Enumerate all iso classes at n_sub
    print(f"\n  Step 1: Iso classes at n={n_sub}")
    all_tournaments_sub = []
    for bits in range(1 << (n_sub*(n_sub-1)//2)):
        T = np.zeros((n_sub, n_sub), dtype=int)
        idx = 0
        for i in range(n_sub):
            for j in range(i+1, n_sub):
                if bits & (1 << idx):
                    T[i][j] = 1
                else:
                    T[j][i] = 1
                idx += 1
        all_tournaments_sub.append(T)

    # Group by canonical form
    classes_sub = defaultdict(list)
    for i, T in enumerate(all_tournaments_sub):
        canon, _ = tournament_canonical(T, n_sub)
        classes_sub[canon].append(i)

    print(f"    {len(classes_sub)} iso classes at n={n_sub}")

    # For each class, pick a representative and compute Aut
    class_reps = {}  # canon → (representative T, aut_group)
    for canon, indices in classes_sub.items():
        T = all_tournaments_sub[indices[0]]
        auts = compute_aut_group(T, n_sub)
        h = count_hamiltonian_paths(T, n_sub)
        class_reps[canon] = {
            'T': T, 'auts': auts, 'aut_size': len(auts),
            'H': h, 'count': len(indices)
        }
        if len(classes_sub) <= 20:
            print(f"    Class {canon[:6]}...: |Aut|={len(auts)}, H={h}, "
                  f"labeled copies={len(indices)}")

    # Step 2: For each class at n_sub, compute extension orbits
    print(f"\n  Step 2: Extension orbits (fiber structure)")

    total_classes_at_n = 0
    class_to_parent = {}  # class at n → parent class at n_sub
    class_to_sigma_orbit = {}  # class at n → orbit info
    all_extensions = []  # list of (parent_canon, sigma_orbit, T_extended, H_ext)

    for canon, info in class_reps.items():
        T_sub = info['T']
        auts = info['auts']

        # Enumerate all 2^{n_sub} extensions
        orbits = {}  # orbit representative → orbit info
        seen_sigmas = set()

        for sigma_int in range(1 << n_sub):
            sigma = tuple((sigma_int >> i) & 1 for i in range(n_sub))
            if sigma in seen_sigmas:
                continue

            # Compute orbit
            orbit = sigma_orbit(sigma, auts, n_sub)
            for s in orbit:
                seen_sigmas.add(s)

            # Pick canonical sigma (lex smallest in orbit)
            canon_sigma = min(orbit)

            # Extend tournament
            T_ext = extend_tournament(T_sub, n_sub, canon_sigma)
            h_ext = count_hamiltonian_paths(T_ext, n_target)
            canon_ext, _ = tournament_canonical(T_ext, n_target)

            orbits[canon_sigma] = {
                'orbit_size': len(orbit),
                'H': h_ext,
                'canon_n': canon_ext,
                'sigma': canon_sigma,
                'score_n': sum(canon_sigma),  # score of new vertex
            }
            total_classes_at_n += 1

            all_extensions.append({
                'parent': canon,
                'sigma': canon_sigma,
                'orbit_size': len(orbit),
                'H': h_ext,
                'canon_n': canon_ext,
            })

        n_orbits = len(orbits)
        if len(classes_sub) <= 20:
            print(f"\n    Parent class H={info['H']}, |Aut|={info['aut_size']}:")
            print(f"      {n_orbits} extension orbits (= {n_orbits} new classes at n={n_target})")
            for cs, orb_info in sorted(orbits.items(), key=lambda x: x[1]['H']):
                print(f"        σ={cs}, orbit_size={orb_info['orbit_size']}, "
                      f"H_ext={orb_info['H']}, score(n)={orb_info['score_n']}")

    # Step 3: Verify total count
    print(f"\n  Step 3: Verification")
    print(f"    Total classes at n={n_target} (from fiber bundle): {total_classes_at_n}")

    # Compare with direct computation
    all_T_n = []
    for bits in range(1 << (n_target*(n_target-1)//2)):
        T = np.zeros((n_target, n_target), dtype=int)
        idx = 0
        for i in range(n_target):
            for j in range(i+1, n_target):
                if bits & (1 << idx):
                    T[i][j] = 1
                else:
                    T[j][i] = 1
                idx += 1
        all_T_n.append(T)

    classes_n_direct = defaultdict(list)
    for i, T in enumerate(all_T_n):
        canon, _ = tournament_canonical(T, n_target)
        classes_n_direct[canon].append(i)

    print(f"    Total classes at n={n_target} (direct): {len(classes_n_direct)}")
    print(f"    Match: {total_classes_at_n == len(classes_n_direct)}")

    # NOTE: the fiber bundle count might OVERCOUNT because different
    # parent classes can give ISOMORPHIC extensions (the new vertex might
    # bridge two non-isomorphic subtournaments into isomorphic n-tournaments).

    # Check: do any extension classes coincide?
    ext_canons = Counter()
    for ext in all_extensions:
        ext_canons[ext['canon_n']] += 1

    multi_parent = {c: n for c, n in ext_canons.items() if n > 1}
    print(f"    Classes with MULTIPLE parents: {len(multi_parent)}")
    if multi_parent:
        print(f"    (These are classes reachable from different parent deletions)")
        for canon_n, count in list(multi_parent.items())[:5]:
            parents = [e['parent'] for e in all_extensions if e['canon_n'] == canon_n]
            h = [e['H'] for e in all_extensions if e['canon_n'] == canon_n][0]
            print(f"      Class H={h}: {count} parents")

    # Step 4: The bundle projection map
    print(f"\n  Step 4: Bundle projection (class at n → parent classes at n-1)")

    # For each class at n, find ALL possible parents (delete each vertex)
    # A class at n can have multiple parents depending on which vertex we delete
    if n_target <= 6:
        print(f"    Computing full projection (all vertex deletions)...")
        class_parents_full = defaultdict(set)  # class_n → set of class_{n-1}

        for canon_n, indices in classes_n_direct.items():
            # Pick one representative
            T = all_T_n[indices[0]]
            for v in range(n_target):
                # Delete vertex v
                T_del = np.zeros((n_sub, n_sub), dtype=int)
                remaining = [u for u in range(n_target) if u != v]
                for ii, i in enumerate(remaining):
                    for jj, j in enumerate(remaining):
                        T_del[ii][jj] = T[i][j]
                canon_del, _ = tournament_canonical(T_del, n_sub)
                class_parents_full[canon_n].add(canon_del)

        # Statistics
        parent_counts = [len(v) for v in class_parents_full.values()]
        print(f"    Parent count distribution: {Counter(parent_counts)}")
        print(f"    Classes with 1 parent: {sum(1 for c in parent_counts if c == 1)}")
        print(f"    Classes with max parents: {max(parent_counts)}")

    # Step 5: The FIBER SIZES
    print(f"\n  Step 5: Fiber sizes (extensions per parent)")
    fiber_sizes = Counter()
    for canon, info in class_reps.items():
        auts = info['auts']
        n_orbits = 0
        seen = set()
        for sigma_int in range(1 << n_sub):
            sigma = tuple((sigma_int >> i) & 1 for i in range(n_sub))
            if sigma in seen:
                continue
            orbit = sigma_orbit(sigma, auts, n_sub)
            for s in orbit:
                seen.add(s)
            n_orbits += 1
        fiber_sizes[n_orbits] += 1

    print(f"    Fiber size distribution: {dict(sorted(fiber_sizes.items()))}")
    print(f"    Total fibers (= total classes at n): "
          f"{sum(k*v for k,v in fiber_sizes.items())}")

    # The formula: fiber_size(C) = 2^{n-1} / |Aut(C)| (by orbit-stabilizer)
    # But orbits can have different sizes, so it's not exact...
    # Actually by Burnside: #orbits = (1/|G|) sum_{g∈G} |Fix(g)|
    # where Fix(g) = {σ : g(σ) = σ}

    print(f"\n    Burnside verification:")
    for canon, info in class_reps.items():
        auts = info['auts']
        aut_size = info['aut_size']
        # Burnside: #orbits = (1/|Aut|) sum_{α∈Aut} #{σ : α(σ)=σ}
        total_fixed = 0
        for alpha in auts:
            # σ is fixed by α iff σ[i] = σ[α(i)] for all i
            # This means σ is constant on each cycle of α
            # Number of such σ = 2^(#cycles of α on {0,...,n_sub-1})
            # Count cycles
            visited = [False]*n_sub
            n_cycles = 0
            for i in range(n_sub):
                if not visited[i]:
                    n_cycles += 1
                    j = i
                    while not visited[j]:
                        visited[j] = True
                        j = alpha[j]
            total_fixed += 2**n_cycles

        n_orbits_burnside = total_fixed // aut_size
        if len(classes_sub) <= 20:
            print(f"      |Aut|={aut_size}: #orbits = {total_fixed}/{aut_size} "
                  f"= {n_orbits_burnside}")

    return all_extensions, class_reps, classes_n_direct

# ==============================================================
# Run for small n
# ==============================================================
if __name__ == "__main__":
    for n in [4, 5, 6]:
        fiber_bundle(n)
