#!/usr/bin/env python3
"""perspectives_investigation.py -- The Perspective-Class Coincidence.

Session: kind-pasteur-2026-03-20-S4

CONCEPT: A "perspective" of a tournament class C is a vertex orbit under Aut(T).
It counts how many distinct ROLES a vertex can play in that tournament.

  - Transitive T_3: 3 roles (source, middle, sink) → 3 perspectives
  - Cyclic C_3: 1 role (all vertices equivalent) → 1 perspective
  - Total at n=3: 3 + 1 = 4

OBSERVATION: P(n) = total perspectives at n equals T(n+1) = isomorphism classes
at n+1, at least for small n:
  - P(3) = 4 = T(4)
  - P(4) = 12 = T(5)

QUESTION: When does P(n) diverge from T(n+1)? Why do they agree initially?
What does this reveal about the structure of tournament isomorphism?

COMPUTATION:
  Part 1: Enumerate all iso classes for n=3..8, count vertex orbits
  Part 2: Compare P(n) vs T(n+1) — find the divergence point
  Part 3: If they diverge, study HOW — which perspectives map to which classes?
  Part 4: Natural maps between rooted tournaments at n and classes at n+1
"""

from itertools import permutations, combinations
from collections import defaultdict, Counter
from math import factorial

# ================================================================
# TOURNAMENT ISOMORPHISM CLASSES
# ================================================================

def adj_from_bits(bits, n):
    adj = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    return adj

def canon_form(adj, n):
    """Canonical form of a tournament (for isomorphism testing).
    Returns the lexicographically smallest bit encoding over all relabelings.
    """
    best = None
    for perm in permutations(range(n)):
        bits = 0
        k = 0
        for i in range(n):
            for j in range(i+1, n):
                if adj[perm[i]][perm[j]]:
                    bits |= (1 << k)
                k += 1
        if best is None or bits < best:
            best = bits
    return best

def find_automorphisms(adj, n):
    """Find all automorphisms of tournament."""
    auts = []
    for perm in permutations(range(n)):
        is_aut = True
        for i in range(n):
            for j in range(i+1, n):
                if adj[perm[i]][perm[j]] != adj[i][j]:
                    is_aut = False
                    break
            if not is_aut:
                break
        if is_aut:
            auts.append(perm)
    return auts

def vertex_orbits(adj, n):
    """Find vertex orbits under Aut(T).
    Returns list of orbits, each orbit is a frozenset of vertex indices.
    """
    auts = find_automorphisms(adj, n)
    # Union-find via orbit computation
    orbit_of = list(range(n))

    def find(x):
        while orbit_of[x] != x:
            orbit_of[x] = orbit_of[orbit_of[x]]
            x = orbit_of[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            orbit_of[rx] = ry

    for aut in auts:
        for v in range(n):
            union(v, aut[v])

    # Collect orbits
    orbits = defaultdict(set)
    for v in range(n):
        orbits[find(v)].add(v)
    return list(orbits.values())

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))

def is_self_converse(adj, n):
    """Check if T is isomorphic to T^op."""
    opp = [[1 - adj[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
    c1 = canon_form(adj, n)
    c2 = canon_form(opp, n)
    return c1 == c2

def held_karp_H(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


# ================================================================
# PART 1: ENUMERATE CLASSES AND COUNT PERSPECTIVES
# ================================================================

def compute_perspectives(max_n=8):
    print("=" * 70)
    print("PART 1: TOURNAMENT CLASSES AND PERSPECTIVES")
    print("=" * 70)

    # Known T(n) values (OEIS A000568)
    known_T = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}

    perspectives_by_n = {}
    classes_by_n = {}

    for n in range(1, max_n + 1):
        m = n * (n - 1) // 2
        N = 1 << m

        if N > 5000000:
            print(f"\n  n={n}: too large for exhaustive enumeration ({N} tournaments)")
            classes_by_n[n] = known_T.get(n, '?')
            continue

        print(f"\n  n={n} ({N} labeled tournaments)")

        # Find all isomorphism classes
        seen_canons = {}  # canon_bits -> representative_bits
        class_data = []

        for bits in range(N):
            adj = adj_from_bits(bits, n)
            c = canon_form(adj, n)
            if c not in seen_canons:
                seen_canons[c] = bits
                # Analyze this class representative
                auts = find_automorphisms(adj, n)
                orbits = vertex_orbits(adj, n)
                sc = score_seq(adj, n)
                is_sc = is_self_converse(adj, n)
                H = held_karp_H(adj, n) if n <= 7 else None

                labeled_count = factorial(n) // len(auts)
                class_data.append({
                    'canon': c, 'bits': bits,
                    'score': sc, 'aut_size': len(auts),
                    'num_orbits': len(orbits),
                    'orbits': orbits,
                    'labeled_count': labeled_count,
                    'is_sc': is_sc,
                    'H': H,
                })

        num_classes = len(class_data)
        total_perspectives = sum(d['num_orbits'] for d in class_data)

        classes_by_n[n] = num_classes
        perspectives_by_n[n] = total_perspectives

        print(f"    Classes: {num_classes} (expected: {known_T.get(n, '?')})")
        print(f"    Total perspectives: {total_perspectives}")

        # Detailed table
        if num_classes <= 20:
            print(f"\n    {'Score':>20} {'|Aut|':>6} {'Orbits':>7} {'SC':>4} {'H':>6} {'Labeled':>8}")
            for d in sorted(class_data, key=lambda x: x['score']):
                H_str = str(d['H']) if d['H'] is not None else '?'
                print(f"    {str(d['score']):>20} {d['aut_size']:>6} {d['num_orbits']:>7} "
                      f"{'Y' if d['is_sc'] else 'N':>4} {H_str:>6} {d['labeled_count']:>8}")

        # SC vs non-SC perspective breakdown
        sc_classes = [d for d in class_data if d['is_sc']]
        nsc_classes = [d for d in class_data if not d['is_sc']]
        sc_persp = sum(d['num_orbits'] for d in sc_classes)
        nsc_persp = sum(d['num_orbits'] for d in nsc_classes)
        print(f"\n    SC classes: {len(sc_classes)}, perspectives: {sc_persp}")
        print(f"    Non-SC classes: {len(nsc_classes)} (in {len(nsc_classes)//2} converse pairs), "
              f"perspectives: {nsc_persp}")

    return perspectives_by_n, classes_by_n


# ================================================================
# PART 2: COMPARE P(n) vs T(n+1)
# ================================================================

def compare_sequences(persp, classes):
    print("\n" + "=" * 70)
    print("PART 2: P(n) vs T(n+1) COMPARISON")
    print("=" * 70)

    known_T = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}

    print(f"\n  {'n':>4} {'P(n)':>10} {'T(n+1)':>10} {'Match':>8} {'Ratio':>10}")
    for n in sorted(persp.keys()):
        p = persp[n]
        t_next = known_T.get(n + 1, '?')
        if isinstance(t_next, int):
            match = 'YES' if p == t_next else 'NO'
            ratio = p / t_next
            print(f"  {n:4d} {p:10d} {t_next:10d} {match:>8} {ratio:10.4f}")
        else:
            print(f"  {n:4d} {p:10d} {str(t_next):>10} {'?':>8} {'?':>10}")

    # THE KEY QUESTION: Where do they diverge?
    print(f"\n  Divergence analysis:")
    for n in sorted(persp.keys()):
        p = persp[n]
        t_next = known_T.get(n + 1)
        if t_next and p != t_next:
            diff = p - t_next
            print(f"  DIVERGENCE at n={n}: P(n)={p}, T(n+1)={t_next}, diff={diff}")
            print(f"    P(n)/T(n+1) = {p/t_next:.6f}")
            if p > t_next:
                print(f"    P(n) EXCEEDS T(n+1) by {diff}")
                print(f"    This means: more perspectives than classes at n+1")
                print(f"    Some classes at n+1 must 'share' perspectives")
            else:
                print(f"    P(n) FALLS SHORT of T(n+1) by {-diff}")
                print(f"    This means: some classes at n+1 have no perspective parent")


# ================================================================
# PART 3: ROOTED TOURNAMENT → CLASS EXTENSION MAP
# ================================================================

def study_extension_map(max_n=6):
    print("\n" + "=" * 70)
    print("PART 3: EXTENSION MAP — ROOTED TOURNAMENTS TO CLASSES")
    print("=" * 70)

    for n in range(3, max_n + 1):
        n_ext = n + 1
        m = n * (n - 1) // 2
        m_ext = n_ext * (n_ext - 1) // 2

        if (1 << m_ext) > 5000000:
            print(f"\n  n={n} -> n+1={n_ext}: too large")
            continue

        print(f"\n  === Extension from n={n} to n+1={n_ext} ===")

        # 1. Find all rooted tournaments (T, v) at n
        # A rooted tournament = iso class + vertex orbit
        rooted = []
        seen = {}
        for bits in range(1 << m):
            adj = adj_from_bits(bits, n)
            c = canon_form(adj, n)
            if c in seen:
                continue
            seen[c] = True
            orbits = vertex_orbits(adj, n)
            for orb in orbits:
                v_rep = min(orb)  # representative vertex
                rooted.append({
                    'bits': bits, 'vertex': v_rep,
                    'orbit_size': len(orb),
                    'score_v': sum(adj[v_rep][j] for j in range(n) if j != v_rep),
                    'class_score': score_seq(adj, n),
                })

        print(f"  Rooted tournaments at n={n}: {len(rooted)}")

        # 2. For each rooted tournament (T, v), enumerate all extensions
        # adding vertex w with all 2^n possible arc orientations
        # and find the iso class of each extension
        extension_map = defaultdict(set)  # (rooted_idx) -> set of extension class canons

        for idx, rt in enumerate(rooted):
            bits = rt['bits']
            adj = adj_from_bits(bits, n)
            v = rt['vertex']

            for new_arcs in range(1 << n):
                # Build extended tournament
                adj_ext = [row[:] + [0] for row in adj]
                adj_ext.append([0] * n_ext)
                w = n  # new vertex index
                for i in range(n):
                    if (new_arcs >> i) & 1:
                        adj_ext[w][i] = 1  # w -> i
                    else:
                        adj_ext[i][w] = 1  # i -> w

                c_ext = canon_form(adj_ext, n_ext)
                extension_map[idx].add(c_ext)

        # 3. Find all classes at n+1
        ext_classes = set()
        for bits in range(1 << m_ext):
            adj = adj_from_bits(bits, n_ext)
            c = canon_form(adj, n_ext)
            ext_classes.add(c)

        print(f"  Classes at n+1={n_ext}: {len(ext_classes)}")

        # 4. Analyze the extension map
        # Which classes at n+1 are reachable from each rooted tournament?
        # Which rooted tournaments reach each class?
        class_to_parents = defaultdict(set)  # class canon -> set of rooted indices
        for idx in range(len(rooted)):
            for c in extension_map[idx]:
                class_to_parents[c].add(idx)

        # Statistics
        reach_counts = [len(extension_map[idx]) for idx in range(len(rooted))]
        parent_counts = [len(class_to_parents[c]) for c in ext_classes]

        print(f"\n  Extension map statistics:")
        print(f"    Rooted tournaments: {len(rooted)}")
        print(f"    Target classes: {len(ext_classes)}")
        print(f"    Reach per rooted: min={min(reach_counts)}, max={max(reach_counts)}, "
              f"avg={sum(reach_counts)/len(reach_counts):.1f}")
        print(f"    Parents per class: min={min(parent_counts)}, max={max(parent_counts)}, "
              f"avg={sum(parent_counts)/len(parent_counts):.1f}")

        # KEY: Is there a BIJECTION between rooted tournaments and classes?
        # This would explain the P(n) = T(n+1) coincidence.
        if len(rooted) == len(ext_classes):
            print(f"\n    P(n) = T(n+1) = {len(rooted)}")
            # Check if each rooted tournament reaches a UNIQUE class
            # (some class reachable from only that rooted tournament)
            unique_matches = 0
            for idx in range(len(rooted)):
                for c in extension_map[idx]:
                    if len(class_to_parents[c]) == 1:
                        unique_matches += 1
                        break
            print(f"    Rooted tournaments with a unique reachable class: {unique_matches}/{len(rooted)}")

            # Try to find a bijection greedily
            unmatched_rooted = set(range(len(rooted)))
            unmatched_classes = set(ext_classes)
            bijection = {}

            # First pass: match rooted tournaments that reach a class with only 1 parent
            for c in ext_classes:
                if len(class_to_parents[c]) == 1:
                    idx = list(class_to_parents[c])[0]
                    if idx in unmatched_rooted:
                        bijection[idx] = c
                        unmatched_rooted.discard(idx)
                        unmatched_classes.discard(c)

            print(f"    After unique-parent matching: {len(bijection)} matched, "
                  f"{len(unmatched_rooted)} rooted and {len(unmatched_classes)} classes remaining")

        else:
            print(f"\n    P(n) = {len(rooted)} != T(n+1) = {len(ext_classes)}")
            print(f"    DIVERGENCE! The coincidence breaks down.")

            # Analyze the mismatch
            if len(rooted) > len(ext_classes):
                # More perspectives than classes: multiple rooted tournaments must map to same class
                print(f"    {len(rooted) - len(ext_classes)} extra perspectives need to share classes")

                # Find classes with multiple parents
                multi_parent = [(c, parents) for c, parents in class_to_parents.items()
                               if len(parents) > 1]
                print(f"    Classes with multiple parent perspectives: {len(multi_parent)}")

                if len(multi_parent) <= 10:
                    for c, parents in multi_parent:
                        parent_info = []
                        for idx in parents:
                            rt = rooted[idx]
                            parent_info.append(f"({rt['class_score']}, v_score={rt['score_v']})")
                        print(f"      Class {c}: parents = {parent_info}")

            else:
                # Fewer perspectives than classes: some classes have no perspective parent
                orphan_classes = [c for c in ext_classes if len(class_to_parents[c]) == 0]
                print(f"    Orphan classes (no perspective parent): {len(orphan_classes)}")

        # Show the detailed mapping for small n
        if len(rooted) <= 15 and len(ext_classes) <= 60:
            print(f"\n  Detailed rooted tournament descriptions:")
            for idx, rt in enumerate(rooted):
                reaches = extension_map[idx]
                print(f"    RT[{idx}]: class_score={rt['class_score']}, "
                      f"v_score={rt['score_v']}, orbit_size={rt['orbit_size']}, "
                      f"reaches {len(reaches)} classes")


# ================================================================
# MAIN
# ================================================================

if __name__ == "__main__":
    persp, classes = compute_perspectives(max_n=7)
    compare_sequences(persp, classes)
    study_extension_map(max_n=6)
