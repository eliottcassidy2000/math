#!/usr/bin/env python3
"""
creative_explorations_s24.py — opus-2026-04-05-S24

THREE GENUINELY NOVEL DIRECTIONS:

1. TOURNAMENT CELLULAR AUTOMATON
   Conway's Game of Life on tournament space. Three update rules:
   - Majority rule: each arc polls 2-paths through all mediators
   - Triangle stress: flip arcs in too many 3-cycles
   - H-gradient ascent: flip the single arc that maximizes H increase
   Analyze: fixed points, limit cycles, basin sizes, attractor landscape.

2. RANDOM WALK PHASE TRANSITION
   Start from transitive tournament, flip random arcs one at a time.
   Track H(T), c3(T), score_var(T) as function of flip count.
   Question: is there a sharp phase transition? At what fraction of
   C(n,2) flips does the system become "disordered"?

3. SPECTRAL ANALYSIS OF THE ISO-CLASS METAGRAPH
   Build the adjacency matrix of G_n (wiggly edges between iso classes).
   Compute full eigenvalue spectrum. The spectral gap controls:
   - Expansion (Cheeger inequality)
   - Mixing time of random walk on iso classes
   - Community structure (Fiedler vector = algebraic connectivity)
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
import random
import time
import sys

sys.stdout.reconfigure(encoding='utf-8')

# ================================================================
# CORE TOURNAMENT FUNCTIONS
# ================================================================

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=np.int8)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def adj_to_bits(A, n):
    bits = 0
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j] == 1:
                bits |= (1 << idx)
            idx += 1
    return bits

def count_ham_paths(A, n):
    """Count Hamiltonian paths via DP bitmask."""
    m = 1 << n
    dp = [[0]*n for _ in range(m)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(m):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            if not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u] == 1:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = m - 1
    return sum(dp[full][v] for v in range(n))

def count_3cycles(A, n):
    """Count directed 3-cycles."""
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check both orientations
                if A[i][j] and A[j][k] and A[k][i]:
                    c += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    c += 1
    return c

def score_sequence(A, n):
    return tuple(sorted(int(A[i].sum()) for i in range(n)))

def score_variance(A, n):
    scores = [int(A[i].sum()) for i in range(n)]
    mean = sum(scores) / n
    return sum((s - mean)**2 for s in scores) / n

def canonical_form(A, n):
    """Isomorphism canonical form via brute-force over S_n (fine for n<=7)."""
    best = None
    for perm in permutations(range(n)):
        bits = 0
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if A[perm[i]][perm[j]] == 1:
                    bits |= (1 << idx)
                idx += 1
        if best is None or bits < best:
            best = bits
    return best

def transitive_tournament(n):
    """The transitive tournament: i beats j iff i > j."""
    A = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1  # higher index beats lower? Actually: i beats j for i < j
    return A

def flip_arc(A, i, j):
    """Flip arc between i and j. Returns new matrix."""
    B = A.copy()
    B[i][j] = 1 - B[i][j]
    B[j][i] = 1 - B[j][i]
    return B


# ================================================================
# EXPERIMENT 1: TOURNAMENT CELLULAR AUTOMATON
# ================================================================

def majority_rule_step(A, n):
    """Each arc (i,j) polls all 2-paths: how many k have i->k->j vs j->k->i?
    If majority says i->j, keep; if majority says j->i, flip.
    Ties: keep current. All updates simultaneous (synchronous CA)."""
    B = A.copy()
    for i in range(n):
        for j in range(i+1, n):
            votes_ij = 0  # votes for i->j
            votes_ji = 0  # votes for j->i
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[i][k] and A[k][j]:
                    votes_ij += 1
                if A[j][k] and A[k][i]:
                    votes_ji += 1
            if votes_ij > votes_ji:
                B[i][j] = 1; B[j][i] = 0
            elif votes_ji > votes_ij:
                B[i][j] = 0; B[j][i] = 1
            # ties: keep current
    return B

def triangle_stress_step(A, n):
    """Each arc (i,j) counts how many 3-cycles it's in.
    If more than (n-2)/2, it's "stressed" and flips.
    All updates simultaneous."""
    B = A.copy()
    threshold = (n - 2) / 2
    for i in range(n):
        for j in range(i+1, n):
            cycle_count = 0
            for k in range(n):
                if k == i or k == j:
                    continue
                # 3-cycle through i,j,k: either i->j->k->i or i->k->j->i
                if A[i][j] and A[j][k] and A[k][i]:
                    cycle_count += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    cycle_count += 1
            if cycle_count > threshold:
                # Flip this arc
                B[i][j] = 1 - A[i][j]
                B[j][i] = 1 - A[j][i]
    return B

def h_gradient_step(A, n):
    """Greedy H-ascent: find the single arc flip that maximizes H. Flip it.
    If no flip increases H, we're at a local maximum."""
    current_H = count_ham_paths(A, n)
    best_delta = 0
    best_ij = None
    for i in range(n):
        for j in range(i+1, n):
            B = flip_arc(A, i, j)
            new_H = count_ham_paths(B, n)
            delta = new_H - current_H
            if delta > best_delta:
                best_delta = delta
                best_ij = (i, j)
    if best_ij is not None:
        return flip_arc(A, best_ij[0], best_ij[1])
    return A  # local max

def run_ca(A, n, rule_fn, max_steps=200, name="CA"):
    """Run a CA rule until fixed point or cycle detected."""
    history = []
    seen = {}
    current = A.copy()
    for step in range(max_steps):
        bits = adj_to_bits(current, n)
        if bits in seen:
            cycle_start = seen[bits]
            cycle_len = step - cycle_start
            return {
                'name': name,
                'steps': step,
                'cycle_start': cycle_start,
                'cycle_len': cycle_len,
                'H_history': [count_ham_paths(bits_to_adj(h, n), n) for h in history],
                'c3_history': [count_3cycles(bits_to_adj(h, n), n) for h in history],
                'final_bits': bits,
                'final_H': count_ham_paths(current, n),
                'final_c3': count_3cycles(current, n),
            }
        seen[bits] = step
        history.append(bits)
        current = rule_fn(current, n)
    return {
        'name': name,
        'steps': max_steps,
        'cycle_start': -1,
        'cycle_len': -1,
        'H_history': [count_ham_paths(bits_to_adj(h, n), n) for h in history],
        'c3_history': [count_3cycles(bits_to_adj(h, n), n) for h in history],
        'final_bits': adj_to_bits(current, n),
        'final_H': count_ham_paths(current, n),
        'final_c3': count_3cycles(current, n),
    }

def experiment_1_tournament_ca():
    print("=" * 70)
    print("EXPERIMENT 1: TOURNAMENT CELLULAR AUTOMATON")
    print("=" * 70)

    for n in [4, 5]:
        m = n * (n - 1) // 2
        total = 2**m
        print(f"\n--- n = {n}, {m} arcs, {total} tournaments ---")

        rules = [
            ("Majority", majority_rule_step),
            ("TriStress", triangle_stress_step),
            ("H-Gradient", h_gradient_step),
        ]

        for rule_name, rule_fn in rules:
            print(f"\n  Rule: {rule_name}")
            # Run from multiple starting points
            fixed_points = Counter()
            cycle_lengths = Counter()
            attractor_count = 0
            max_trials = min(total, 200)

            if total <= 200:
                starts = range(total)
            else:
                starts = random.sample(range(total), max_trials)

            t0 = time.time()
            for bits in starts:
                A = bits_to_adj(bits, n)
                result = run_ca(A, n, rule_fn, max_steps=100, name=rule_name)
                if result['cycle_len'] == 1:
                    fixed_points[result['final_bits']] += 1
                elif result['cycle_len'] > 0:
                    cycle_lengths[result['cycle_len']] += 1
                else:
                    attractor_count += 1
            elapsed = time.time() - t0

            print(f"    Trials: {max_trials}, Time: {elapsed:.1f}s")
            print(f"    Fixed points: {len(fixed_points)} distinct ({sum(fixed_points.values())} total convergences)")
            if fixed_points:
                for fp_bits, count in fixed_points.most_common(5):
                    fp_A = bits_to_adj(fp_bits, n)
                    H = count_ham_paths(fp_A, n)
                    c3 = count_3cycles(fp_A, n)
                    sc = score_sequence(fp_A, n)
                    print(f"      bits={fp_bits}: H={H}, c3={c3}, score={sc}, basin={count}")
            if cycle_lengths:
                print(f"    Limit cycles: {dict(cycle_lengths)}")
            if attractor_count:
                print(f"    Non-converged: {attractor_count}")

            # Detailed trajectory from transitive tournament
            T = transitive_tournament(n)
            result = run_ca(T, n, rule_fn, max_steps=50, name=rule_name)
            Hh = result['H_history'][:20]
            print(f"    From transitive: H trajectory = {Hh}")
            print(f"      -> settled to H={result['final_H']}, c3={result['final_c3']}, cycle_len={result['cycle_len']}")


# ================================================================
# EXPERIMENT 2: RANDOM WALK PHASE TRANSITION
# ================================================================

def experiment_2_phase_transition():
    print("\n" + "=" * 70)
    print("EXPERIMENT 2: RANDOM WALK PHASE TRANSITION")
    print("=" * 70)

    for n in [5, 6, 7]:
        m = n * (n - 1) // 2
        print(f"\n--- n = {n}, m = {m} arcs ---")

        num_trials = 200 if n <= 6 else 50
        # For each trial, start from transitive and flip random arcs
        # Track H, c3, score_var at each step
        max_flips = 3 * m  # flip up to 3x the number of arcs

        # Accumulate statistics at each step
        H_by_step = defaultdict(list)
        c3_by_step = defaultdict(list)
        sv_by_step = defaultdict(list)

        t0 = time.time()
        for trial in range(num_trials):
            A = transitive_tournament(n)
            flipped = set()  # track which arcs have been flipped (can re-flip)
            for step in range(max_flips):
                # Record current state
                if step % max(1, m // 10) == 0 or step < m:
                    H = count_ham_paths(A, n)
                    c3 = count_3cycles(A, n)
                    sv = score_variance(A, n)
                    H_by_step[step].append(H)
                    c3_by_step[step].append(c3)
                    sv_by_step[step].append(sv)

                # Pick a random arc to flip
                i = random.randint(0, n-2)
                j = random.randint(i+1, n-1)
                A = flip_arc(A, i, j)

        elapsed = time.time() - t0
        print(f"  {num_trials} trials x {max_flips} flips, time: {elapsed:.1f}s")

        # Print trajectory summary
        steps_sorted = sorted(H_by_step.keys())
        print(f"\n  Step | mean(H)    | std(H)   | mean(c3) | mean(sv)")
        print(f"  -----|------------|----------|----------|--------")
        for step in steps_sorted[:30]:
            Hs = H_by_step[step]
            c3s = c3_by_step[step]
            svs = sv_by_step[step]
            print(f"  {step:4d} | {np.mean(Hs):10.1f} | {np.std(Hs):8.1f} | {np.mean(c3s):8.2f} | {np.mean(svs):6.3f}")

        # Find the "transition point" — steepest H increase
        if len(steps_sorted) >= 3:
            mean_H = [np.mean(H_by_step[s]) for s in steps_sorted]
            max_jump = 0
            jump_step = 0
            for i in range(1, len(mean_H)):
                jump = mean_H[i] - mean_H[i-1]
                if jump > max_jump:
                    max_jump = jump
                    jump_step = steps_sorted[i]
            print(f"\n  Steepest H increase at step {jump_step} (delta_H = {max_jump:.1f})")
            print(f"  H ranges: start={mean_H[0]:.1f}, end={mean_H[-1]:.1f}")
            # Theoretical: random tournament mean H
            if n <= 7:
                # mean H for random tournament = n! / 2^(n-1)
                import math
                mean_random = math.factorial(n) / 2**(n-1)
                print(f"  Theoretical mean(H) for random tournament: {mean_random:.1f}")
                # Transition midpoint: step where H crosses half of (max - min)
                mid_H = (mean_H[0] + mean_random) / 2
                for i, s in enumerate(steps_sorted):
                    if mean_H[i] >= mid_H:
                        print(f"  H crosses midpoint ({mid_H:.1f}) at step {s} = {s/m:.2f}*m")
                        break


# ================================================================
# EXPERIMENT 3: SPECTRAL ANALYSIS OF THE METAGRAPH
# ================================================================

def build_iso_class_metagraph(n):
    """Build the iso class metagraph at level n.
    Returns: classes (list of canonical bits), class_H, adjacency matrix."""
    m = n * (n - 1) // 2
    total = 2**m

    # Enumerate all tournaments, group by canonical form
    canon_to_class = {}
    class_members = defaultdict(list)

    print(f"  Enumerating all {total} tournaments on n={n}...")
    t0 = time.time()
    for bits in range(total):
        A = bits_to_adj(bits, n)
        cf = canonical_form(A, n)
        if cf not in canon_to_class:
            canon_to_class[cf] = len(canon_to_class)
        class_members[canon_to_class[cf]].append(bits)
    elapsed = time.time() - t0
    num_classes = len(canon_to_class)
    print(f"  {num_classes} iso classes found in {elapsed:.1f}s")

    # Map each bits -> class_id
    bits_to_class = {}
    for bits in range(total):
        A = bits_to_adj(bits, n)
        cf = canonical_form(A, n)
        bits_to_class[bits] = canon_to_class[cf]

    # Compute H and c3 for each class
    class_H = {}
    class_c3 = {}
    class_scores = {}
    for cf, cid in canon_to_class.items():
        A = bits_to_adj(cf, n)
        class_H[cid] = count_ham_paths(A, n)
        class_c3[cid] = count_3cycles(A, n)
        class_scores[cid] = score_sequence(A, n)

    # Build wiggly adjacency matrix (d=1 layer: single arc flips)
    print(f"  Building wiggly adjacency matrix...")
    adj = np.zeros((num_classes, num_classes), dtype=int)
    edge_set = set()
    for bits in range(total):
        A = bits_to_adj(bits, n)
        src_class = bits_to_class[bits]
        for i in range(n):
            for j in range(i+1, n):
                B = flip_arc(A, i, j)
                new_bits = adj_to_bits(B, n)
                dst_class = bits_to_class[new_bits]
                if src_class != dst_class:
                    edge_set.add((min(src_class, dst_class), max(src_class, dst_class)))
                    adj[src_class][dst_class] += 1
                    adj[dst_class][src_class] += 1

    # Make symmetric adjacency (binary: is there at least one connection?)
    adj_binary = (adj > 0).astype(int)
    np.fill_diagonal(adj_binary, 0)
    num_edges = len(edge_set)
    print(f"  {num_edges} edges in the metagraph")

    # Also build complement-merged version
    # Two classes are complement-paired if T^op of one is isomorphic to the other
    canon_to_complement = {}
    for cf in canon_to_class:
        A = bits_to_adj(cf, n)
        # Complement: flip all arcs
        Ac = np.zeros_like(A)
        for i in range(n):
            for j in range(n):
                if i != j:
                    Ac[i][j] = 1 - A[i][j]
        cf_comp = canonical_form(Ac, n)
        canon_to_complement[cf] = cf_comp

    # Identify SC classes
    sc_classes = set()
    for cf, cid in canon_to_class.items():
        if cf == canon_to_complement[cf]:
            sc_classes.add(cid)

    return {
        'n': n,
        'num_classes': num_classes,
        'adj_weighted': adj,
        'adj_binary': adj_binary,
        'edge_set': edge_set,
        'class_H': class_H,
        'class_c3': class_c3,
        'class_scores': class_scores,
        'sc_classes': sc_classes,
        'canon_to_class': canon_to_class,
        'bits_to_class': bits_to_class,
        'class_members': class_members,
    }


def experiment_3_spectral():
    print("\n" + "=" * 70)
    print("EXPERIMENT 3: SPECTRAL ANALYSIS OF ISO-CLASS METAGRAPH")
    print("=" * 70)

    for n in [4, 5, 6]:
        print(f"\n{'='*60}")
        print(f"n = {n}")
        print(f"{'='*60}")

        if n == 6:
            # n=6 has 56 iso classes, canonical_form is expensive
            print("  (n=6 may take a while due to canonical form computation...)")

        t0 = time.time()
        G = build_iso_class_metagraph(n)
        build_time = time.time() - t0
        print(f"  Build time: {build_time:.1f}s")

        nc = G['num_classes']
        adj = G['adj_binary']
        adj_w = G['adj_weighted']

        print(f"\n  Basic metagraph stats:")
        print(f"    Vertices (iso classes): {nc}")
        print(f"    Edges: {len(G['edge_set'])}")
        print(f"    SC classes: {len(G['sc_classes'])}")
        degrees = adj.sum(axis=1)
        print(f"    Degree range: {degrees.min()} to {degrees.max()}, mean={degrees.mean():.1f}")
        density = 2 * len(G['edge_set']) / (nc * (nc - 1)) if nc > 1 else 0
        print(f"    Density: {density:.3f}")

        # Eigenvalue spectrum
        eigenvalues = np.sort(np.linalg.eigvalsh(adj.astype(float)))[::-1]
        print(f"\n  Eigenvalue spectrum (binary adj):")
        print(f"    lambda_1 = {eigenvalues[0]:.4f}")
        if nc > 1:
            print(f"    lambda_2 = {eigenvalues[1]:.4f}")
            gap = eigenvalues[0] - eigenvalues[1]
            print(f"    Spectral gap = {gap:.4f}")
            ratio = eigenvalues[1] / eigenvalues[0] if eigenvalues[0] > 0 else 0
            print(f"    lambda_2/lambda_1 = {ratio:.4f}")
        if nc > 2:
            print(f"    lambda_3 = {eigenvalues[2]:.4f}")
        print(f"    lambda_min = {eigenvalues[-1]:.4f}")
        # Algebraic connectivity (second smallest eigenvalue of Laplacian)
        L = np.diag(degrees.astype(float)) - adj.astype(float)
        lap_eigs = np.sort(np.linalg.eigvalsh(L))
        print(f"\n  Laplacian spectrum:")
        print(f"    mu_1 (should be 0) = {lap_eigs[0]:.6f}")
        if nc > 1:
            print(f"    mu_2 (algebraic connectivity) = {lap_eigs[1]:.4f}")
            print(f"    mu_max = {lap_eigs[-1]:.4f}")
            # Cheeger inequality: h(G) >= mu_2 / 2
            print(f"    Cheeger bound: h(G) >= {lap_eigs[1]/2:.4f}")
            # Mixing time bound: t_mix <= log(n)/(mu_2)
            if lap_eigs[1] > 0:
                mix_bound = np.log(nc) / lap_eigs[1]
                print(f"    Mixing time bound: t_mix <= {mix_bound:.2f}")

        # Fiedler vector (eigenvector of mu_2)
        if nc > 1:
            evals, evecs = np.linalg.eigh(L)
            fiedler = evecs[:, 1]  # second eigenvector
            # Sort classes by Fiedler value
            order = np.argsort(fiedler)
            print(f"\n  Fiedler ordering (algebraic 'geometry' of metagraph):")
            print(f"    class_id: H, c3, score, SC?, Fiedler_val")
            for idx in order:
                is_sc = "SC" if idx in G['sc_classes'] else "  "
                print(f"    {idx:3d}: H={G['class_H'][idx]:5d}, c3={G['class_c3'][idx]:3d}, "
                      f"score={G['class_scores'][idx]}, {is_sc}, f={fiedler[idx]:+.4f}")

        # Weighted spectrum (number of labeled tournament connections)
        if nc > 1:
            eigs_w = np.sort(np.linalg.eigvalsh(adj_w.astype(float)))[::-1]
            print(f"\n  Weighted adjacency spectrum (multiplicity of connections):")
            print(f"    lambda_1 = {eigs_w[0]:.1f}")
            print(f"    lambda_2 = {eigs_w[1]:.1f}")
            print(f"    Spectral gap = {eigs_w[0] - eigs_w[1]:.1f}")

        # Chromatic number / clique number from adjacency
        max_clique = 0
        for size in range(nc, 0, -1):
            if size <= max_clique:
                break
            if size > 10 and nc > 20:
                break  # skip expensive computation
            found = False
            for subset in combinations(range(nc), size):
                if all(adj[i][j] for i in subset for j in subset if i != j):
                    max_clique = size
                    found = True
                    break
            if not found and size <= max_clique + 1:
                break
        print(f"\n  Clique number (lower bound): omega >= {max_clique}")

        # Check if metagraph is bipartite
        # (it won't be, but let's verify)
        coloring = [-1] * nc
        coloring[0] = 0
        queue = [0]
        is_bipartite = True
        while queue and is_bipartite:
            v = queue.pop(0)
            for u in range(nc):
                if adj[v][u] and u != v:
                    if coloring[u] == -1:
                        coloring[u] = 1 - coloring[v]
                        queue.append(u)
                    elif coloring[u] == coloring[v]:
                        is_bipartite = False
        print(f"  Bipartite: {is_bipartite}")

        # H-distribution by Fiedler position
        if nc > 1:
            # Correlation between Fiedler value and H
            H_vals = np.array([G['class_H'][i] for i in range(nc)])
            c3_vals = np.array([G['class_c3'][i] for i in range(nc)])
            corr_fiedler_H = np.corrcoef(fiedler, H_vals)[0, 1]
            corr_fiedler_c3 = np.corrcoef(fiedler, c3_vals)[0, 1]
            print(f"\n  Correlations:")
            print(f"    corr(Fiedler, H) = {corr_fiedler_H:.4f}")
            print(f"    corr(Fiedler, c3) = {corr_fiedler_c3:.4f}")

        print()


# ================================================================
# MAIN
# ================================================================

if __name__ == '__main__':
    experiment_1_tournament_ca()
    experiment_2_phase_transition()
    experiment_3_spectral()
