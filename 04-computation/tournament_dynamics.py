"""
tournament_dynamics.py -- kind-pasteur-2026-03-14-S71
Tournament dynamics: arc-flip as discrete evolution.

CREATIVE IDEAS:
1. "Hamiltonian flow" on tournament space: flip arc that maximizes dH
2. Random walk on the H-landscape: what is the mixing time?
3. "Tournament automata": each vertex updates based on neighbors
4. Markov chain on tournaments: H as energy, Boltzmann distribution
5. The "H-potential" landscape and its gradient flow

But the REALLY creative idea:
CONWAY'S GAME OF LIFE ON TOURNAMENTS
Each arc evolves based on the triangle structure around it.
Rule: arc (i,j) flips if the majority of 3-cycles through (i,j) are "stressed"
where a 3-cycle is stressed if it creates an odd cycle.

OR:
THRESHOLD DYNAMICS: each arc flips to match the majority of its
"support triangles." If more triangles containing (i,j) have i→j as the
majority direction, keep it; otherwise flip.

OR even simpler:
MAJORITY RULE ON ARCS: for each pair (i,j), count how many paths
of length 2 go i→k→j vs j→k→i. If i beats j through more mediators,
keep i→j; otherwise flip.
"""

import numpy as np
from collections import Counter, defaultdict
import sys

sys.stdout.reconfigure(encoding='utf-8')

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
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
            if A[j][i] == 1:
                bits |= (1 << idx)
            idx += 1
    return bits

def compute_H_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def majority_step(A, n):
    """
    Majority dynamics: for each pair (i,j), flip to the majority
    direction among all 2-paths through third vertices.

    For each (i,j): count mediators k where i→k→j vs j→k→i.
    If more mediators support i→j, set i→j. Otherwise set j→i.
    Ties: keep current direction.
    """
    A_new = A.copy()
    for i in range(n):
        for j in range(i+1, n):
            # Count 2-paths i→k→j vs j→k→i
            support_ij = 0  # i→k and k→j
            support_ji = 0  # j→k and k→i
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[i][k] == 1 and A[k][j] == 1:
                    support_ij += 1
                if A[j][k] == 1 and A[k][i] == 1:
                    support_ji += 1

            if support_ij > support_ji:
                A_new[i][j] = 1
                A_new[j][i] = 0
            elif support_ji > support_ij:
                A_new[j][i] = 1
                A_new[i][j] = 0
            # Ties: keep current

    return A_new

def triangle_stress_step(A, n):
    """
    Triangle stress dynamics: for each arc (i→j), count how many
    3-cycles contain this arc. If the arc participates in MORE 3-cycles
    than its reverse would, flip it.
    """
    A_new = A.copy()
    for i in range(n):
        for j in range(i+1, n):
            # Count 3-cycles through (i,j) in current direction
            cycles_current = 0
            cycles_flipped = 0

            for k in range(n):
                if k == i or k == j:
                    continue
                # With current i→j:
                if A[i][j] == 1:
                    # 3-cycle i→j→k→i: need A[j][k]=1 and A[k][i]=1
                    if A[j][k] == 1 and A[k][i] == 1:
                        cycles_current += 1
                    # 3-cycle k→i→j→k: need A[k][i]=1 and A[j][k]=0
                    # Wait, this is the same cycle
                else:
                    # j→i: cycle j→i→k→j: need A[i][k]=1 and A[k][j]=1
                    if A[i][k] == 1 and A[k][j] == 1:
                        cycles_current += 1

            # Count cycles if we flip
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[i][j] == 1:
                    # After flip: j→i. Cycle j→i→k→j: A[i][k] and A[k][j]
                    if A[i][k] == 1 and A[k][j] == 1:
                        cycles_flipped += 1
                else:
                    # After flip: i→j. Cycle i→j→k→i: A[j][k] and A[k][i]
                    if A[j][k] == 1 and A[k][i] == 1:
                        cycles_flipped += 1

            # Flip to minimize 3-cycles (transitivize)
            if cycles_flipped < cycles_current:
                A_new[i][j], A_new[j][i] = A_new[j][i], A_new[i][j]

    return A_new

def main():
    print("=" * 70)
    print("TOURNAMENT DYNAMICS")
    print("kind-pasteur-2026-03-14-S71")
    print("=" * 70)

    # ========================================
    # PART 1: Majority dynamics
    # ========================================
    print("\n" + "=" * 70)
    print("PART 1: MAJORITY DYNAMICS")
    print("  Each arc aligns with the majority of 2-path mediators")
    print("  Question: Does this converge? To what?")
    print("=" * 70)

    for n in [5, 6, 7]:
        print(f"\n--- n = {n} ---")
        np.random.seed(42)

        num_trials = 20 if n <= 6 else 10
        fixed_points = Counter()
        cycle_lengths = Counter()

        for trial in range(num_trials):
            bits = np.random.randint(0, 2**(n*(n-1)//2))
            A = bits_to_adj(bits, n)
            H_init = compute_H_dp(A, n)

            # Run majority dynamics
            trajectory = [adj_to_bits(A, n)]
            visited = {trajectory[0]: 0}

            for step in range(50):
                A = majority_step(A, n)
                bits_new = adj_to_bits(A, n)

                if bits_new in visited:
                    # Found cycle
                    cycle_start = visited[bits_new]
                    cycle_len = step + 1 - cycle_start
                    H_final = compute_H_dp(A, n)

                    if cycle_len == 1:
                        fixed_points[H_final] += 1
                    cycle_lengths[cycle_len] += 1
                    break

                visited[bits_new] = step + 1
                trajectory.append(bits_new)
            else:
                cycle_lengths['none'] += 1

            if trial < 3:
                H_traj = [compute_H_dp(bits_to_adj(b, n), n) for b in trajectory[:10]]
                print(f"  Trial {trial}: H trajectory = {H_traj}")

        print(f"\n  Fixed point H values: {dict(fixed_points)}")
        print(f"  Cycle lengths: {dict(cycle_lengths)}")

    # ========================================
    # PART 2: Triangle stress dynamics (transitivization)
    # ========================================
    print("\n" + "=" * 70)
    print("PART 2: TRIANGLE STRESS DYNAMICS (TRANSITIVIZATION)")
    print("  Each arc flips to minimize 3-cycles through it")
    print("  Expected: converges to transitive tournament")
    print("=" * 70)

    for n in [5, 6]:
        print(f"\n--- n = {n} ---")
        np.random.seed(123)

        for trial in range(10):
            bits = np.random.randint(0, 2**(n*(n-1)//2))
            A = bits_to_adj(bits, n)
            H_init = compute_H_dp(A, n)
            c3_init = int(np.trace(A @ A @ A)) // 3

            trajectory_H = [H_init]
            trajectory_c3 = [c3_init]

            for step in range(20):
                A_prev = A.copy()
                A = triangle_stress_step(A, n)

                if np.array_equal(A, A_prev):
                    break  # fixed point

                H = compute_H_dp(A, n)
                c3 = int(np.trace(A @ A @ A)) // 3
                trajectory_H.append(H)
                trajectory_c3.append(c3)

            H_final = compute_H_dp(A, n)
            c3_final = int(np.trace(A @ A @ A)) // 3

            if trial < 5:
                print(f"  Trial {trial}: H: {trajectory_H}, c3: {trajectory_c3}")

    # ========================================
    # PART 3: Simulated annealing for H-maximization
    # ========================================
    print("\n" + "=" * 70)
    print("PART 3: SIMULATED ANNEALING FOR H-MAXIMIZATION")
    print("  Can SA escape the H=37 trap at n=6?")
    print("=" * 70)

    import random as rng

    for n in [6, 7]:
        print(f"\n--- n = {n} ---")
        m = n * (n-1) // 2
        rng.seed(42)

        best_H_list = []

        for trial in range(20 if n == 6 else 10):
            bits = rng.randint(0, 2**m - 1)
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            best_H = H
            best_bits = bits

            temp = 10.0

            for step in range(200):
                # Random arc flip
                arc_bit = rng.randint(0, m-1)
                new_bits = bits ^ (1 << arc_bit)
                A_new = bits_to_adj(new_bits, n)
                H_new = compute_H_dp(A_new, n)

                delta = H_new - H

                if delta > 0 or rng.random() < np.exp(delta / temp):
                    bits = new_bits
                    A = A_new
                    H = H_new

                    if H > best_H:
                        best_H = H
                        best_bits = bits

                temp *= 0.99  # cool down

            best_H_list.append(best_H)

        H_dist = Counter(best_H_list)
        print(f"  SA results ({len(best_H_list)} runs): {dict(sorted(H_dist.items()))}")
        if n == 6:
            print(f"  Fraction reaching H=45: {H_dist.get(45, 0)}/{len(best_H_list)}")
        if n == 7:
            print(f"  Fraction reaching H=189: {H_dist.get(189, 0)}/{len(best_H_list)}")

    # ========================================
    # PART 4: Ergodicity — random walk mixing time
    # ========================================
    print("\n" + "=" * 70)
    print("PART 4: RANDOM WALK MIXING ON H-LANDSCAPE")
    print("  Random walk: flip a random arc at each step")
    print("  How quickly does the H distribution converge to equilibrium?")
    print("=" * 70)

    n = 5
    m = n * (n-1) // 2
    rng.seed(42)

    # Run many random walks and record H trajectory
    num_walks = 200
    walk_length = 50
    H_at_step = defaultdict(list)

    for walk in range(num_walks):
        bits = rng.randint(0, 2**m - 1)
        for step in range(walk_length):
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            H_at_step[step].append(H)
            # Random flip
            arc_bit = rng.randint(0, m-1)
            bits ^= (1 << arc_bit)

    print(f"\n  n={n}, {num_walks} random walks of length {walk_length}")
    print(f"  Mean H at each step:")
    for step in range(0, walk_length, 5):
        vals = H_at_step[step]
        print(f"    step {step:3d}: mean={np.mean(vals):.2f}, std={np.std(vals):.2f}")

    print(f"\n  Equilibrium mean H = n!/2^(n-1) = {np.math.factorial(n)/2**(n-1):.2f}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)

if __name__ == '__main__':
    import math
    np.math = math
    main()
