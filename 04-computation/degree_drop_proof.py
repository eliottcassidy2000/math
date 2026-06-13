"""
degree_drop_proof.py -- kind-pasteur-2026-03-14-S72
Prove the Degree Drop Theorem: deg(H) = 2*floor((n-1)/2).

STRATEGY: Show that for even n, the coefficient of the top-degree monomial
x_{e1}*...*x_{e_{n-1}} vanishes for EVERY set of n-1 arcs.

The coefficient c_S of monomial prod_{e in S} x_e in H(T) equals:
c_S = sum over paths P using arcs in S: (-1)^{# reversed arcs in P relative to S}

More precisely: H(T) = sum_P prod_step T[P[i]][P[i+1]]
where T[a][b] = x_{ab} if a<b, or (1-x_{ab}) if a>b.

For a path P = (v_0,...,v_{n-1}), the expansion gives:
prod_{i=0}^{n-2} T[v_i][v_{i+1}]
= prod_{i: v_i < v_{i+1}} x_{v_i,v_{i+1}} * prod_{i: v_i > v_{i+1}} (1 - x_{v_{i+1},v_i})

The degree-(n-1) term of this product requires selecting the "x" part of
every factor, which means:
- For ascending step (v_i < v_{i+1}): contribute x_{v_i,v_{i+1}} with coeff +1
- For descending step (v_i > v_{i+1}): contribute -x_{v_{i+1},v_i} with coeff -1

So the degree-(n-1) coefficient for the arc set S(P) = {canonical arcs used by P} is:
(-1)^{des(P)} where des(P) = #{i : v_i > v_{i+1}} = #descents

THEREFORE: c_S = sum over paths P with arc set S: (-1)^{des(P)}

For the degree drop at even n, we need: for every valid arc set S of size n-1
that forms a Hamiltonian path, sum (-1)^{des(P)} = 0 over all paths P using S.

But wait: a set of n-1 arcs can form at most ONE Hamiltonian path (since n-1
arcs on n vertices with maximum degree constraints determine a unique path).
So c_S = (-1)^{des(P)} for the unique path P using S, or c_S = 0 if no path
uses S. This would mean c_S = +1 or -1, never 0, and the degree CAN'T drop!

CONTRADICTION with our computation. What am I missing?

Ah wait: the arc set S(P) is the set of CANONICAL arcs {min(v_i,v_{i+1}), max(...)}.
Multiple paths CAN use the same canonical arc set! Because the DIRECTION can differ.

No: the canonical arc set determines which pairs are adjacent in the path,
but NOT the direction. A canonical arc {a,b} can be traversed a->b or b->a.

So: given a set S of n-1 canonical arcs (unordered pairs), the set of paths
using these arcs includes both the original path AND any path that traverses
some edges in reverse. But wait, reversing a single edge in a path might
break the Hamiltonian path property...

Actually, the canonical arc set determines the UNDERLYING GRAPH of the path
(an undirected path graph on n vertices). This graph has exactly 2 orientations
as a directed path: the forward direction and the reverse direction.

So for each canonical arc set S forming a Hamiltonian path, there are exactly
2 directed paths: P and P^rev. And:
des(P) + des(P^rev) = n-1 (since each step is ascending in P iff descending in P^rev).

Therefore: c_S = (-1)^{des(P)} + (-1)^{des(P^rev)} = (-1)^{des(P)} + (-1)^{n-1-des(P)}

If n is even (n-1 is odd): (-1)^{des} + (-1)^{n-1-des} = (-1)^{des}(1 + (-1)^{n-1-2des})
= (-1)^{des}(1 + (-1)^{n-1}) since (-1)^{-2des} = 1
= (-1)^{des}(1 - 1) = 0.  BECAUSE n-1 IS ODD, (-1)^{n-1} = -1.

If n is odd (n-1 is even): (-1)^{des} + (-1)^{n-1-des}
= (-1)^{des}(1 + (-1)^{n-1-2des})
The sign of n-1-2des depends on des. But (-1)^{n-1} = +1 (even),
so this is (-1)^{des} + (-1)^{-des} * (-1)^{n-1} = (-1)^des + (-1)^{-des}
Wait let me redo this.

(-1)^{n-1-des} = (-1)^{n-1} * (-1)^{-des} = (-1)^{n-1} * (-1)^{des}
since (-1)^{-des} = (-1)^{des} (because (-1)^2 = 1).

Wait no: (-1)^{-k} = 1/(-1)^k = (-1)^{-k}. And (-1)^{-1} = -1. So
(-1)^{-k} = (-1)^k if k is even, (-1)^{-k} = (-1)^k for all k actually,
because (-1)^{2k} = 1 so (-1)^{-k} * (-1)^k = (-1)^0 = 1, hence
(-1)^{-k} = (-1)^{-k} = 1/(-1)^k = (-1)^{k} when (-1)^{2k}=1 which is
always true. So (-1)^{-k} = (-1)^k for all integer k. RIGHT!

So: (-1)^{n-1-des} = (-1)^{n-1} * (-1)^{-des} = (-1)^{n-1} * (-1)^{des}

c_S = (-1)^{des} + (-1)^{n-1} * (-1)^{des} = (-1)^{des} * (1 + (-1)^{n-1})

If n even (n-1 odd): 1 + (-1)^{n-1} = 1 + (-1) = 0. So c_S = 0. QED!
If n odd (n-1 even): 1 + (-1)^{n-1} = 1 + 1 = 2. So c_S = 2*(-1)^{des}.

THIS PROVES THE DEGREE DROP THEOREM!

For EVEN n: every top-degree coefficient is 0, so degree <= n-2.
For ODD n: top-degree coefficients are +-2, so degree = n-1.

The involution is PATH REVERSAL: P <-> P^rev.
The key: des(P) + des(P^rev) = n-1 (complementary descents).
For even n: (-1)^{des} + (-1)^{n-1-des} = 0 (opposite signs).
For odd n: (-1)^{des} + (-1)^{n-1-des} = 2*(-1)^{des} (same sign).

BEAUTIFUL! Let me verify computationally.
"""

import numpy as np
from itertools import permutations
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

def main():
    print("=" * 70)
    print("DEGREE DROP THEOREM — PROOF BY PATH REVERSAL INVOLUTION")
    print("kind-pasteur-2026-03-14-S72")
    print("=" * 70)

    print("""
THEOREM: For tournament on n vertices, the multilinear polynomial H(T)
in the C(n,2) arc variables has degree:
  deg(H) = n-1    if n is odd
  deg(H) = n-2    if n is even

PROOF:
Step 1: H(T) = sum over permutations P: prod_{i=0}^{n-2} T[P[i],P[i+1]]
where T[a,b] = x_{ab} if a<b, else (1-x_{ba}) if a>b.

Step 2: The degree-(n-1) coefficient of monomial prod_{e in S} x_e is:
c_S = sum over paths P with canonical arc set = S: (-1)^{des(P)}
where des(P) = #{i : P[i] > P[i+1]} = number of descents.

Step 3: Each canonical arc set S that forms a valid path has exactly 2
directed paths using it: P and P^rev (the reversed path). These satisfy:
  des(P) + des(P^rev) = n-1
because each step is ascending in P iff descending in P^rev.

Step 4: c_S = (-1)^{des(P)} + (-1)^{des(P^rev)}
       = (-1)^{des(P)} + (-1)^{n-1-des(P)}
       = (-1)^{des(P)} * [1 + (-1)^{n-1}]

Step 5: If n is EVEN: n-1 is odd, so 1 + (-1)^{n-1} = 1-1 = 0. QED!
        If n is ODD: n-1 is even, so 1 + (-1)^{n-1} = 1+1 = 2.
        So c_S = 2*(-1)^{des(P)}, which is +-2 (nonzero).

The involution is PATH REVERSAL, which maps des(P) to n-1-des(P).
For even n, path and its reversal have opposite signs.
For odd n, they have the same sign.
""")

    # VERIFICATION
    print("=" * 70)
    print("VERIFICATION")
    print("=" * 70)

    for n in [3, 4, 5, 6]:
        print(f"\n--- n = {n}, n-1 = {n-1} ({'odd' if (n-1)%2==1 else 'even'}) ---")

        # For each valid path, compute (-1)^des and check pairing
        arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
        arc_idx = {a: i for i, a in enumerate(arcs)}

        # Group paths by their canonical arc set
        arc_set_to_paths = defaultdict(list)
        for perm in permutations(range(n)):
            # Check: is this a valid sequence? (All consecutive pairs are arcs)
            # YES — every permutation defines a Hamiltonian path in SOME tournament
            # The canonical arc set is the set of unordered consecutive pairs
            canonical = frozenset(
                (min(perm[i], perm[i+1]), max(perm[i], perm[i+1]))
                for i in range(n-1)
            )
            des = sum(1 for i in range(n-1) if perm[i] > perm[i+1])
            arc_set_to_paths[canonical].append((perm, des))

        print(f"  Distinct canonical arc sets: {len(arc_set_to_paths)}")
        print(f"  Total permutations: {len(list(permutations(range(n))))}")

        # Check: each arc set should have exactly 2 paths (P and P^rev)
        sizes = Counter(len(v) for v in arc_set_to_paths.values())
        print(f"  Arc set sizes: {dict(sizes)}")

        # Check coefficient = (-1)^des(P) + (-1)^des(P^rev) for each arc set
        all_zero = True
        nonzero_coeffs = []
        for S, paths in arc_set_to_paths.items():
            coeff = sum((-1)**des for _, des in paths)
            if coeff != 0:
                all_zero = False
                nonzero_coeffs.append(coeff)

        if n % 2 == 0:
            print(f"  n EVEN: all top-degree coefficients = 0? {all_zero}")
            if not all_zero:
                print(f"    COUNTEREXAMPLE! Some nonzero: {nonzero_coeffs[:5]}")
        else:
            print(f"  n ODD: nonzero top-degree coefficients: {Counter(nonzero_coeffs)}")
            all_pm2 = all(abs(c) == 2 for c in nonzero_coeffs)
            print(f"    All coefficients = +-2? {all_pm2}")

        # Verify des(P) + des(P^rev) = n-1
        complement_check = True
        for S, paths in arc_set_to_paths.items():
            if len(paths) == 2:
                d1, d2 = paths[0][1], paths[1][1]
                if d1 + d2 != n-1:
                    complement_check = False
                    print(f"    FAIL: des={d1}, des_rev={d2}, sum={d1+d2} != {n-1}")

        print(f"  des(P) + des(P^rev) = n-1 for all pairs? {complement_check}")

    # CONSEQUENCES
    print("\n" + "=" * 70)
    print("CONSEQUENCES OF THE DEGREE DROP THEOREM")
    print("=" * 70)

    print("""
1. VASSILIEV TYPE:
   The Vassiliev type of H(T) = degree of multilinear polynomial
   = 2*floor((n-1)/2) for all n.
   Type n-1 for odd n, type n-2 for even n.

2. DELTA VANISHING:
   Delta_{k+1} H = 0 for k = degree = 2*floor((n-1)/2).
   This means the (k+1)-th order finite difference vanishes for ALL
   choices of k+1 arcs (not just special configurations).

3. PARITY STRUCTURE:
   For even n, H(T) mod 2 is a polynomial of degree <= n-2.
   This might give additional mod-2 constraints on H.

4. COEFFICIENT BOUND:
   For odd n, all top-degree coefficients = +-2.
   For even n, all top-degree coefficients = 0.
   The next-to-top coefficients have different structure.

5. INVARIANT THEORY:
   H(T) is invariant under the REVERSAL INVOLUTION at even n
   in a very specific sense: the involution acts as -1 on the
   top-degree part, forcing it to vanish.
""")

    # CHECK: Does degree = n-2 (not just n-2 upper bound) at even n?
    # I.e., are there nonzero degree-(n-2) coefficients?
    print("Checking: degree EXACTLY n-2 at even n?")
    for n in [4, 6]:
        arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
        m = len(arcs)

        # Expand H as polynomial
        poly = Counter()
        for perm in permutations(range(n)):
            factors = []
            for step in range(n-1):
                i, j = perm[step], perm[step+1]
                arc = (min(i,j), max(i,j))
                idx = arcs.index(arc)
                positive = (i < j)
                factors.append((idx, positive))

            for mask in range(1 << (n-1)):
                coeff = 1
                arc_set = set()
                for bit_pos in range(n-1):
                    arc_id, positive = factors[bit_pos]
                    if mask & (1 << bit_pos):
                        if not positive:
                            coeff *= -1
                        arc_set.add(arc_id)
                    else:
                        if positive:
                            coeff = 0
                            break
                if coeff != 0:
                    poly[frozenset(arc_set)] += coeff

        deg_n_minus_2 = sum(1 for S, c in poly.items() if len(S) == n-2 and c != 0)
        max_deg = max((len(S) for S, c in poly.items() if c != 0), default=0)
        print(f"  n={n}: max degree = {max_deg}, nonzero deg-(n-2) monomials = {deg_n_minus_2}")

    print("\n" + "=" * 70)
    print("PROOF COMPLETE")
    print("=" * 70)

if __name__ == '__main__':
    main()
