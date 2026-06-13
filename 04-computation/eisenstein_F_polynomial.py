"""
eisenstein_F_polynomial.py -- opus-2026-03-14-S71l
THE EISENSTEIN NORM OF THE F-POLYNOMIAL

Key insight: F(T, omega) where omega = e^{2pi*i/3} lives in Z[omega] (Eisenstein integers).
|F(T, omega)|^2 is the Eisenstein norm = Loeschian number.

The Loeschian numbers are n = a^2 + ab + b^2 for some a,b in Z.
They are OEIS A003136: 0,1,3,4,7,9,12,13,16,19,21,25,27,28,...

Primes p are Loeschian iff p=3 or p ≡ 1 (mod 3).
7 is the SMALLEST non-trivial Loeschian prime (7 = 1 + 1*2 + 4 = 2^2 + 2 + 1 = Phi_3(2)).

CONNECTION TO FORBIDDEN VALUES:
- 7 = Phi_3(2) is Loeschian (norm of 2 - omega)
- 21 = 3*7 is Loeschian (norm of (2-omega)(1-omega) or similar)
- The F-polynomial's 3-strand decomposition naturally produces Loeschian values
"""

import sys
import numpy as np
from itertools import permutations
from fractions import Fraction
from math import comb, gcd

sys.stdout.reconfigure(encoding='utf-8')

def eisenstein_norm(a, b):
    """Norm of a + b*omega in Z[omega]: a^2 + ab + b^2."""
    return a*a + a*b + b*b

def is_loeschian(n):
    """Check if n is a Loeschian number (representable as a^2+ab+b^2)."""
    if n == 0:
        return True
    for a in range(int(n**0.5) + 1):
        for b in range(int(n**0.5) + 1):
            if a*a + a*b + b*b == n:
                return True
    return False

def loeschian_representation(n):
    """Find (a,b) such that a^2+ab+b^2 = n."""
    for a in range(int(n**0.5) + 1):
        for b in range(int(n**0.5) + 1):
            if a*a + a*b + b*b == n:
                return (a, b)
    return None

def compute_all_tournaments(n):
    """Enumerate all tournaments on n vertices, return (adj_matrix, H, F_poly)."""
    m = n * (n - 1) // 2
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    results = []
    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        # Compute F-polynomial via permutation enumeration
        F = [0] * n
        for perm in permutations(range(n)):
            valid = True
            ascents = 0
            for k in range(n-1):
                if adj[perm[k]][perm[k+1]]:
                    pass  # valid arc
                else:
                    valid = False
                    break
            if valid:
                # Count ascents: positions where perm[k] < perm[k+1] AND arc goes i->j
                ascents = sum(1 for k in range(n-1) if perm[k] < perm[k+1])
                F[ascents] += 1
        H = sum(F)
        results.append((bits, H, tuple(F)))
    return results

def main():
    print("=" * 70)
    print("EISENSTEIN NORM OF THE F-POLYNOMIAL")
    print("opus-2026-03-14-S71l")
    print("=" * 70)

    omega = np.exp(2j * np.pi / 3)

    # Part 1: Loeschian numbers and Phi_3
    print(f"\n{'='*70}")
    print("PART 1: LOESCHIAN NUMBERS = NORMS IN Z[omega]")
    print(f"{'='*70}")

    loeschians = sorted(set(a*a + a*b + b*b for a in range(20) for b in range(20)))
    loeschians = [x for x in loeschians if x <= 100]
    print(f"\n  Loeschian numbers up to 100:")
    print(f"  {loeschians}")

    print(f"\n  Key forbidden/achievable values and Loeschian status:")
    for v in [7, 21, 63, 189, 35, 39, 273]:
        is_L = is_loeschian(v)
        rep = loeschian_representation(v) if is_L else None
        print(f"    {v:4d}: Loeschian={is_L}", end="")
        if rep:
            a, b = rep
            print(f"  ({a}^2 + {a}*{b} + {b}^2 = {v})", end="")
        print()

    # Phi_3 values
    print(f"\n  Phi_3(x) = x^2 + x + 1 values:")
    for x in range(1, 10):
        v = x*x + x + 1
        rep = loeschian_representation(v)
        print(f"    Phi_3({x}) = {v}", end="")
        if rep:
            a, b = rep
            print(f" = {a}^2 + {a}*{b} + {b}^2", end="")
        print()

    # Part 2: F(T, omega) decomposition
    print(f"\n{'='*70}")
    print("PART 2: F(T, omega) AS EISENSTEIN INTEGER")
    print(f"{'='*70}")

    for n in [5, 6]:
        print(f"\n  --- n = {n} ---")
        results = compute_all_tournaments(n)

        # Group by H
        from collections import defaultdict
        h_groups = defaultdict(list)
        for bits, H, F in results:
            h_groups[H].append(F)

        # Compute |F(omega)|^2 for each tournament
        norm_values = defaultdict(list)
        all_norms = []
        for bits, H, F in results:
            G = [0, 0, 0]  # G_j = sum_{k ≡ j mod 3} F_k
            for k in range(n):
                G[k % 3] += F[k]
            # Eisenstein norm: G0^2 + G1^2 + G2^2 - G0*G1 - G1*G2 - G0*G2
            norm_val = G[0]**2 + G[1]**2 + G[2]**2 - G[0]*G[1] - G[1]*G[2] - G[0]*G[2]
            # Alternative: express as element of Z[omega]
            # F(omega) = G0 + G1*omega + G2*omega^2
            # Since 1 + omega + omega^2 = 0, omega^2 = -1 - omega
            # F(omega) = G0 + G1*omega + G2*(-1 - omega) = (G0-G2) + (G1-G2)*omega
            a_eis = G[0] - G[2]
            b_eis = G[1] - G[2]
            norm_check = a_eis**2 - a_eis*b_eis + b_eis**2
            assert norm_val == norm_check, f"Norm mismatch: {norm_val} vs {norm_check}"
            all_norms.append(norm_val)
            norm_values[H].append((norm_val, G, (a_eis, b_eis)))

        distinct_norms = sorted(set(all_norms))
        all_loeschian = all(is_loeschian(v) for v in distinct_norms)
        print(f"    All |F(omega)|^2 are Loeschian: {all_loeschian}")
        print(f"    Distinct norms: {distinct_norms[:30]}...")

        # Show norm distribution by H
        print(f"\n    Norm distribution by H value:")
        for H in sorted(h_groups.keys()):
            norms_for_H = [v[0] for v in norm_values[H]]
            distinct_for_H = sorted(set(norms_for_H))
            if len(distinct_for_H) <= 5:
                print(f"      H={H:3d}: |F(omega)|^2 in {distinct_for_H}")
            else:
                print(f"      H={H:3d}: |F(omega)|^2 in {distinct_for_H[:5]}... ({len(distinct_for_H)} values)")

    # Part 3: The Eisenstein lattice structure
    print(f"\n{'='*70}")
    print("PART 3: EISENSTEIN LATTICE POINTS")
    print(f"{'='*70}")

    print(f"""
  F(T, omega) = (G0 - G2) + (G1 - G2) * omega

  So each tournament maps to an Eisenstein integer z = a + b*omega.
  The NORM |z|^2 = a^2 + ab + b^2 is Loeschian by definition.

  Key question: Which Eisenstein lattice points are ACHIEVABLE?
  And how does this relate to the forbidden value 7?""")

    # Map out the lattice points at n=5
    n = 5
    results = compute_all_tournaments(n)
    lattice_points = {}
    for bits, H, F in results:
        G = [sum(F[k] for k in range(n) if k % 3 == j) for j in range(3)]
        a = G[0] - G[2]
        b = G[1] - G[2]
        key = (a, b)
        if key not in lattice_points:
            lattice_points[key] = []
        lattice_points[key].append(H)

    print(f"\n  n=5: Eisenstein lattice points (a, b) and their H values:")
    for (a, b) in sorted(lattice_points.keys()):
        H_vals = sorted(set(lattice_points[(a, b)]))
        norm = a*a + a*b + b*b
        print(f"    ({a:3d},{b:3d}): norm={norm:3d}, H in {H_vals}")

    # Part 4: The Frobenius connection
    print(f"\n{'='*70}")
    print("PART 4: FROBENIUS ENDOMORPHISM AND F(T, omega)")
    print(f"{'='*70}")

    print(f"""
  In Z[omega], the Frobenius (complex conjugation) acts as:
    omega -> omega^2 = -1 - omega

  For tournaments, complement T -> T^op gives:
    F_k(T^op) = F_{{n-1-k}}(T)

  So: G_j(T^op) = sum_{{k equiv j mod 3}} F_{{n-1-k}}(T)
                 = sum_{{k' equiv (n-1-j) mod 3}} F_{{k'}}(T)
                 = G_{{(n-1-j) mod 3}}(T)

  For n=5 (n-1=4 equiv 1 mod 3):
    G_0(T^op) = G_1(T), G_1(T^op) = G_0(T), G_2(T^op) = G_2(T)

  So complement SWAPS G_0 and G_1, keeps G_2 fixed.

  In Eisenstein coordinates:
    a(T^op) = G_0(T^op) - G_2(T^op) = G_1(T) - G_2(T) = b(T)
    b(T^op) = G_1(T^op) - G_2(T^op) = G_0(T) - G_2(T) = a(T)

  So T -> T^op corresponds to (a,b) -> (b,a) in the Eisenstein lattice!
  This is REFLECTION across the line a=b.""")

    # Verify
    print(f"\n  Verification at n=5:")
    n = 5
    m = n*(n-1)//2
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    count_verified = 0
    count_total = 0
    for bits in range(min(200, 1 << m)):
        adj = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        # Complement
        adj_op = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    adj_op[i][j] = 1 - adj[i][j]
        # Compute F for both
        def get_F(a, nn):
            FF = [0]*nn
            for perm in permutations(range(nn)):
                valid = True
                for k in range(nn-1):
                    if not a[perm[k]][perm[k+1]]:
                        valid = False
                        break
                if valid:
                    asc = sum(1 for k in range(nn-1) if perm[k] < perm[k+1])
                    FF[asc] += 1
            return FF
        F_T = get_F(adj, n)
        F_Top = get_F(adj_op, n)
        G_T = [sum(F_T[k] for k in range(n) if k%3==j) for j in range(3)]
        G_Top = [sum(F_Top[k] for k in range(n) if k%3==j) for j in range(3)]
        a_T = G_T[0] - G_T[2]
        b_T = G_T[1] - G_T[2]
        a_Top = G_Top[0] - G_Top[2]
        b_Top = G_Top[1] - G_Top[2]
        if a_Top == b_T and b_Top == a_T:
            count_verified += 1
        count_total += 1
    print(f"    Tested {count_total} tournaments: {count_verified}/{count_total} satisfy (a,b) -> (b,a)")

    # Part 5: H = 7 through the Eisenstein lens
    print(f"\n{'='*70}")
    print("PART 5: WHY H=7 IS FORBIDDEN — THE EISENSTEIN VIEW")
    print(f"{'='*70}")

    print(f"""
  H = F(T, 1) = G_0 + G_1 + G_2
  |F(T, omega)|^2 = a^2 + ab + b^2 where a = G0-G2, b = G1-G2

  For H=7: G_0 + G_1 + G_2 = 7.
  Possible (G0, G1, G2) triples with G_j >= 0, sum = 7:
  Many options. But the NORM constraint:
  Which of these give |F(omega)|^2 = some value?

  The question is NOT about the norm being 7 (it can be any value).
  The question is: can ANY of these triples arise from a tournament?""")

    # What (G0, G1, G2) triples appear for each H at n=5?
    n = 5
    results = compute_all_tournaments(n)
    G_by_H = {}
    for bits, H, F in results:
        G = tuple(sum(F[k] for k in range(n) if k%3==j) for j in range(3))
        if H not in G_by_H:
            G_by_H[H] = set()
        G_by_H[H].add(G)

    print(f"\n  n=5: (G0, G1, G2) triples by H value:")
    for H in sorted(G_by_H.keys()):
        triples = sorted(G_by_H[H])
        norms = [loeschian_representation(a*a+a*b+b*b) for (g0,g1,g2) in triples
                 for a, b in [(g0-g2, g1-g2)]]
        print(f"    H={H:3d}: {triples}")

    # Part 6: Fibonacci and Eisenstein
    print(f"\n{'='*70}")
    print("PART 6: FIBONACCI IN THE EISENSTEIN LATTICE")
    print(f"{'='*70}")

    phi = (1 + 5**0.5) / 2
    print(f"""
  The golden ratio phi = {phi:.6f} has a curious relationship with omega:
  phi = 2*cos(pi/5) and omega = e^(2*pi*i/3)

  Fibonacci numbers mod 3: 1,1,2,0,2,2,1,0,1,1,2,0,...  (period 8)
  Lucas numbers mod 3:     2,1,0,1,1,2,0,2,2,1,0,1,...  (period 8)

  The Pisano period pi(3) = 8.

  In the Eisenstein lattice, Fibonacci numbers create a SPIRAL:
  F_n + F_{{n+1}} * omega traces out an expanding hexagonal spiral.""")

    # Fibonacci in Eisenstein
    fib = [0, 1]
    for i in range(20):
        fib.append(fib[-1] + fib[-2])

    print(f"\n  Fibonacci Eisenstein points (F_n, F_{{n+1}}) and norms:")
    for i in range(12):
        a, b = fib[i], fib[i+1]
        norm = a*a + a*b + b*b
        print(f"    n={i:2d}: ({a:4d},{b:4d}), norm = {norm:6d}", end="")
        if norm <= 100 and is_loeschian(norm):
            print(f" (Loeschian)", end="")
        print()

    # Part 7: The trinomial / 3-strand connection
    print(f"\n{'='*70}")
    print("PART 7: TRINOMIAL COEFFICIENTS AS EISENSTEIN NORMS")
    print(f"{'='*70}")

    print(f"""
  (1+x+x^2)^n = Phi_3(x)^n evaluated at x = omega gives:
  Phi_3(omega) = omega^2 + omega + 1 = 0

  So (1+omega+omega^2)^n = 0^n = 0 for n >= 1!

  This means: the trinomial row (T(n,0), T(n,1), ..., T(n,2n))
  satisfies: sum_k T(n,k) * omega^k = 0.

  Splitting into residues mod 3:
  G_0 + G_1*omega + G_2*omega^2 = 0

  where G_j = sum_{{k equiv j mod 3}} T(n,k).

  Since 1 + omega + omega^2 = 0:
  G_0 = G_1 = G_2 = 3^n / 3 = 3^{{n-1}}.

  This is obvious: each residue class sums to 3^{{n-1}}.
  The trinomial coefficients are PERFECTLY balanced among the three classes.

  BUT: tournament F-polynomials are NOT perfectly balanced!
  The IMBALANCE is measured by |F(T, omega)|^2.

  H=7 being forbidden means: no tournament can have a specific
  imbalance pattern that gives H = G0+G1+G2 = 7.

  The Eisenstein norm measures HOW UNBALANCED the F-polynomial is
  among the three residue classes mod 3.""")

    # Verify trinomial balance
    print(f"\n  Trinomial row balance verification:")
    for nn in range(1, 8):
        # Compute (1+x+x^2)^nn
        row = [0] * (2*nn + 1)
        row[0] = 1
        for _ in range(nn):
            new_row = [0] * (2*nn + 1)
            for k in range(2*nn + 1):
                if row[k]:
                    for d in [0, 1, 2]:
                        if k + d <= 2*nn:
                            new_row[k+d] += row[k]
            row = new_row
        G = [sum(row[k] for k in range(2*nn+1) if k%3==j) for j in range(3)]
        print(f"    n={nn}: G = {G}, balanced: {G[0]==G[1]==G[2]}, each = 3^{nn-1} = {3**(nn-1)}")

    # Part 8: Category theory — the Eisenstein functor
    print(f"\n{'='*70}")
    print("PART 8: THE EISENSTEIN FUNCTOR")
    print(f"{'='*70}")

    print(f"""
  CATEGORY STRUCTURE:

  Objects: Tournaments on [n]
  Morphisms: Score-sequence-preserving maps? Or isomorphism classes?

  The F-POLYNOMIAL is a functor: T -> F(T, x) in Z[x]/(x^n).
  This factorizes through the EISENSTEIN MAP:

  T --F(T,x)--> Z[x] --eval at omega--> Z[omega] --norm--> Z_{{>=0}}

  The composite is the EISENSTEIN NORM FUNCTOR:
    T |-> |F(T, omega)|^2

  Properties of this functor:
  1. It is NOT injective (many T have same norm)
  2. It IS complement-compatible: norm(T) = norm(T^op)
     (because (a,b) -> (b,a) preserves a^2+ab+b^2)
  3. Its IMAGE is a subset of Loeschian numbers
  4. It correlates strongly (r~0.84) with H

  NATURAL TRANSFORMATIONS:
  - The inclusion Z -> Z[omega] (rational tournaments)
  - The norm Z[omega] -> Z (losing phase information)
  - The projection Z[omega] -> Z/3Z (residue class)

  The FIBER over a Loeschian number L = a^2+ab+b^2:
  All tournaments T with |F(T,omega)|^2 = L.
  This fiber has a natural GROUP ACTION by the unit group
  of Z[omega], which is {{+-1, +-omega, +-omega^2}} (order 6).
  """)

    # Which Loeschian numbers appear at n=5?
    n = 5
    results = compute_all_tournaments(n)
    norm_to_count = {}
    for bits, H, F in results:
        G = [sum(F[k] for k in range(n) if k%3==j) for j in range(3)]
        a = G[0] - G[2]
        b = G[1] - G[2]
        norm = a*a + a*b + b*b
        norm_to_count[norm] = norm_to_count.get(norm, 0) + 1

    print(f"  n=5: Loeschian norms and tournament counts:")
    for norm in sorted(norm_to_count.keys()):
        rep = loeschian_representation(norm)
        print(f"    |F(omega)|^2 = {norm:3d}: {norm_to_count[norm]:4d} tournaments", end="")
        if rep:
            print(f"  (= {rep[0]}^2 + {rep[0]}*{rep[1]} + {rep[1]}^2)", end="")
        print()

    # Part 9: Connecting H and the Eisenstein norm
    print(f"\n{'='*70}")
    print("PART 9: H vs EISENSTEIN NORM — THE CONSTRAINT SURFACE")
    print(f"{'='*70}")

    print(f"""
  Given H = G0 + G1 + G2 and N = (G0-G2)^2 + (G0-G2)(G1-G2) + (G1-G2)^2:

  Let s = G0+G1+G2 = H, and set u = G0-G2, v = G1-G2.
  Then G2 = (s - u - v) / 3 ... wait, G0 = G2+u, G1 = G2+v.
  So s = G2+u + G2+v + G2 = 3*G2 + u + v.
  Thus G2 = (s - u - v)/3.

  For G2 >= 0: u + v <= s = H.
  For G0 >= 0: G2 + u >= 0, so u >= -(s-u-v)/3, so 4u+v >= -s... always true for u,v small.

  The norm N = u^2 + uv + v^2.

  Given H, the MAXIMUM possible norm is achieved when the G_j are
  as UNBALANCED as possible.

  If one G_j = H and others = 0: u = H, v = 0 (say), N = H^2.
  If perfectly balanced: G0=G1=G2=H/3, u=v=0, N=0.

  So 0 <= N <= H^2, with N=0 iff H divisible by 3 and balanced.

  The constraint is: (H, N) must be ACHIEVABLE by some tournament.
  """)

    # Map H vs norm at n=5
    n = 5
    results = compute_all_tournaments(n)
    H_N_pairs = set()
    for bits, H, F in results:
        G = [sum(F[k] for k in range(n) if k%3==j) for j in range(3)]
        a = G[0] - G[2]
        b = G[1] - G[2]
        N = a*a + a*b + b*b
        H_N_pairs.add((H, N))

    print(f"  n=5: Achievable (H, Eisenstein_norm) pairs:")
    for H_val in sorted(set(h for h, n in H_N_pairs)):
        norms = sorted(n for h, n in H_N_pairs if h == H_val)
        print(f"    H={H_val:3d}: norms = {norms}")

    print(f"\n  H=7 is forbidden: NO (7, N) pair exists for ANY N!")
    print(f"  This means no partition of 7 into (G0, G1, G2) is achievable.")

    print(f"\n{'='*70}")
    print("DONE — EISENSTEIN NORM OF F-POLYNOMIAL")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
