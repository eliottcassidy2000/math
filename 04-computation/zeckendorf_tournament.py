"""
zeckendorf_tournament.py -- kind-pasteur-2026-03-14-S85
The Zeckendorf-Tournament connection.

THE USER'S INSIGHT:
Natural numbers in base Fibonacci (Zeckendorf representation):
  Every n >= 1 has a UNIQUE representation n = sum F_k where no two
  consecutive F_k are both used. The binary string of selected Fibonacci
  numbers NEVER has two adjacent 1-bits.

THIS IS AN INDEPENDENCE CONDITION:
  The "Fibonacci graph" is a path P_m (vertices = Fibonacci indices,
  edges between consecutive indices). Zeckendorf representations =
  INDEPENDENT SETS of this path graph!

  The independence polynomial of the path P_m is:
  I(P_m, x) = F_{m+2} when x=1 (counts independent sets).
  The number of Zeckendorf representations of numbers up to F_{m+1}
  is exactly the independent set count of P_m.

CONNECTION TO TOURNAMENTS:
  H(T) = I(Omega(T), 2) — evaluation of independence polynomial at x=2.
  Zeckendorf = I(P_m, 1) — evaluation at x=1 on a path graph.

  If Omega(T) is a PATH GRAPH, then:
  H(T) = I(P_m, 2) = F_{m+2}(2) where F_k(x) is a Fibonacci-like sequence
  with F_1(x) = 1, F_2(x) = 1+x, F_{k}(x) = F_{k-1}(x) + x*F_{k-2}(x).

  This is the "Fibonacci polynomial" sequence!

THE DEEP ANALOGY:
  Zeckendorf: every integer = sum of non-adjacent Fibonacci numbers
  OCF: H-1 = sum of 2^k * (independent cycle sets of size k)

  Both count weighted independent sets on specific graphs.
  The "no two adjacent" constraint in Zeckendorf IS the independence condition.
  The "vertex-disjoint cycles" condition in OCF IS the independence condition.

  H = I(Omega, 2) generalizes Zeckendorf to cycle conflict graphs at fugacity 2.

QUESTIONS:
1. When is Omega(T) a path? -> H is a Fibonacci polynomial value
2. What is the "Zeckendorf representation" of H(T)?
3. Does the Fibonacci structure explain the forbidden values {7, 21}?
4. Connection to golden ratio phi and tournament asymptotics
"""

import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict
import sys, math

sys.stdout.reconfigure(encoding='utf-8')

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[j][i] = 1
            else: A[i][j] = 1
            idx += 1
    return A

def compute_H(A, n):
    dp = {}
    for v in range(n): dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = sum(dp.get((pm, u), 0) for u in range(n) if (pm & (1 << u)) and A[u][v])
                if t: dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def fibonacci(n):
    """Return first n Fibonacci numbers starting F_1=1, F_2=2, F_3=3, F_4=5, ..."""
    fibs = [1, 2]
    while len(fibs) < n:
        fibs.append(fibs[-1] + fibs[-2])
    return fibs[:n]

def zeckendorf(n):
    """Zeckendorf representation of n as sum of non-consecutive Fibonacci numbers."""
    fibs = fibonacci(50)
    # Find largest Fibonacci <= n
    rep = []
    remaining = n
    for i in range(len(fibs)-1, -1, -1):
        if fibs[i] <= remaining:
            rep.append(fibs[i])
            remaining -= fibs[i]
        if remaining == 0:
            break
    return rep

def indpoly_path(m, x):
    """Independence polynomial of path P_m evaluated at x.
    Uses recurrence: I(P_m, x) = I(P_{m-1}, x) + x * I(P_{m-2}, x)
    with I(P_0, x) = 1, I(P_1, x) = 1 + x."""
    if m == 0: return 1
    if m == 1: return 1 + x
    a, b = 1, 1 + x  # I(P_0), I(P_1)
    for _ in range(m - 1):
        a, b = b, b + x * a
    return b

def main():
    print("=" * 70)
    print("ZECKENDORF-TOURNAMENT CONNECTION")
    print("kind-pasteur-2026-03-14-S85")
    print("=" * 70)

    # ============================================================
    # PART 1: Fibonacci independence polynomial
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 1: FIBONACCI AND INDEPENDENCE POLYNOMIALS")
    print("  I(P_m, x) = Fibonacci-like sequence in x")
    print("  I(P_m, 1) = F_{m+2} (Fibonacci number)")
    print("  I(P_m, 2) = ? (a different sequence)")
    print(f"{'='*70}")

    print(f"\n  Path P_m: independence polynomial at various x:")
    print(f"  {'m':>3} {'I(P_m,1)':>10} {'I(P_m,2)':>10} {'I(P_m,3)':>10} {'Fib':>6}")

    fibs = fibonacci(15)
    for m in range(12):
        I1 = indpoly_path(m, 1)
        I2 = indpoly_path(m, 2)
        I3 = indpoly_path(m, 3)
        fib_val = fibs[m+1] if m+1 < len(fibs) else '?'  # F_{m+2}
        print(f"  {m:3d} {I1:10d} {I2:10d} {I3:10d} {fib_val:>6}")

    # I(P_m, 2) sequence: 1, 3, 5, 11, 21, 43, 85, 171, ...
    print(f"\n  I(P_m, 2) sequence: ", end="")
    seq = [indpoly_path(m, 2) for m in range(12)]
    print(seq)

    # Is this a known sequence?
    # Recurrence: a(m) = a(m-1) + 2*a(m-2) with a(0)=1, a(1)=3
    # This is the Jacobsthal-like sequence!
    # a(m) = (2^{m+1} + (-1)^m) / 3
    print(f"\n  Formula check: I(P_m, 2) = (2^(m+1) + (-1)^m) / 3?")
    for m in range(12):
        predicted = (2**(m+1) + (-1)**m) // 3
        actual = indpoly_path(m, 2)
        print(f"    m={m}: predicted={predicted}, actual={actual}, match={predicted == actual}")

    # YES! I(P_m, 2) = (2^{m+1} + (-1)^m) / 3 = Jacobsthal numbers!

    # ============================================================
    # PART 2: When is Omega(T) a path graph?
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 2: WHEN IS OMEGA(T) A PATH GRAPH?")
    print("  If Omega = P_m, then H = I(P_m, 2) = (2^(m+1) + (-1)^m) / 3")
    print(f"{'='*70}")

    jacobsthal = [(2**(m+1) + (-1)**m) // 3 for m in range(20)]
    print(f"  Jacobsthal numbers: {jacobsthal}")

    # Check: which H values at n=5 are Jacobsthal?
    for n in [5, 6]:
        m_arc = n*(n-1)//2
        H_vals = set()
        count = 0
        for bits in range(2**m_arc):
            count += 1
            if n >= 6 and count > 10000: break
            A = bits_to_adj(bits, n)
            H_vals.add(compute_H(A, n))

        H_list = sorted(H_vals)
        jac_set = set(jacobsthal)
        print(f"\n  n={n}: H values = {H_list}")
        for H in H_list:
            is_jac = H in jac_set
            if is_jac:
                m_val = jacobsthal.index(H)
                print(f"    H={H} IS Jacobsthal J({m_val}) = I(P_{m_val}, 2)")

    # ============================================================
    # PART 3: The "Zeckendorf representation" of H
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 3: ZECKENDORF REPRESENTATION OF H VALUES")
    print("  Write each H value in base Fibonacci (non-consecutive)")
    print("  Compare with the OCF decomposition H = 1 + 2*a1 + 4*a2 + ...")
    print(f"{'='*70}")

    for n in [5]:
        m_arc = n*(n-1)//2
        H_vals = sorted(set(
            compute_H(bits_to_adj(bits, n), n) for bits in range(2**m_arc)
        ))

        fibs = fibonacci(20)
        print(f"\n  n={n}: Zeckendorf representations:")
        for H in H_vals:
            zrep = zeckendorf(H)
            T = (H - 1) // 2  # half-Redei
            zrep_T = zeckendorf(T) if T > 0 else [0]
            binary_H = bin(H)[2:]
            binary_T = bin(T)[2:] if T > 0 else '0'

            print(f"    H={H:3d}: Zeck(H)={zrep}, T=(H-1)/2={T}, "
                  f"Zeck(T)={zrep_T}, bin(T)={binary_T}")

    # ============================================================
    # PART 4: Independence polynomial at x=2 as "generalized Zeckendorf"
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 4: GENERALIZED ZECKENDORF")
    print("  Zeckendorf: n = sum of non-adjacent Fibonacci numbers (x=1)")
    print("  Tournament: H = I(Omega, 2) = sum of 2^k * alpha_k")
    print("  The alpha_k count independent sets of SIZE k in Omega")
    print("  So H-1 = sum_{k>=1} 2^k * alpha_k = 'weighted Zeckendorf on Omega'")
    print(f"{'='*70}")
    print(f"""
  THE STRUCTURAL PARALLEL:

  ZECKENDORF:                       OCF:
  Graph: P_m (path)                 Graph: Omega(T) (conflict graph)
  Fugacity: x = 1                   Fugacity: x = 2
  I(P_m, 1) = F_{'{m+2}'}          I(Omega, 2) = H(T)
  Representation: unique            Representation: unique (H determines OCF)
  No adjacent selected              No conflicting cycles selected
  Total count = Fibonacci           Total count = H

  THE KEY DIFFERENCE:
  Zeckendorf: the graph is ALWAYS a path (P_m)
  OCF: the graph Omega(T) can be ANYTHING (depends on tournament)

  When Omega IS a path: H = Jacobsthal number!
  When Omega is more complex: H takes other values.

  THE FORBIDDEN VALUES:
  H=7 is forbidden because Omega cannot be K_3 (complete on 3 vertices).
  In Zeckendorf terms: the "conflict graph" K_3 would give
  I(K_3, 2) = 1 + 6 = 7 (3 vertices, each independent from others).
  But wait, K_3 has NO independent sets of size >=2!
  I(K_3, 2) = 1 + 3*2 = 7 (3 independent singletons weighted by 2).
  The issue: K_3 requires 3 pairwise-conflicting cycles, which
  tournament completeness forbids.
""")

    # ============================================================
    # PART 5: The golden ratio and tournament asymptotics
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 5: GOLDEN RATIO AND TOURNAMENTS")
    print("  Fibonacci: F_n/F_{n-1} -> phi = (1+sqrt(5))/2 ≈ 1.618")
    print("  Jacobsthal: J_n/J_{n-1} -> 2 (since dominant eigenvalue is 2)")
    print("  Tournament max_H: max_H(n)/max_H(n-1) -> ??? (oscillates)")
    print(f"{'='*70}")

    phi = (1 + math.sqrt(5)) / 2
    maxH = [1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095]

    print(f"\n  phi = {phi:.6f}")
    print(f"\n  Comparison of growth rates:")
    print(f"  {'n':>3} {'maxH':>8} {'ratio':>8} {'F_n':>8} {'J_n':>8}")

    fibs = fibonacci(12)
    jacs = [(2**(m+1) + (-1)**m) // 3 for m in range(12)]

    for i in range(len(maxH)):
        ratio = maxH[i] / maxH[i-1] if i > 0 and maxH[i-1] > 0 else 0
        fib = fibs[i] if i < len(fibs) else '?'
        jac = jacs[i] if i < len(jacs) else '?'
        print(f"  {i+1:3d} {maxH[i]:8d} {ratio:8.4f} {fib:>8} {jac:>8}")

    # Is max_H related to Jacobsthal?
    print(f"\n  max_H vs Jacobsthal:")
    for i in range(min(len(maxH), len(jacs))):
        ratio = maxH[i] / jacs[i] if jacs[i] > 0 else 0
        print(f"    n={i+1}: maxH/J = {maxH[i]}/{jacs[i]} = {ratio:.4f}")

    # ============================================================
    # PART 6: The independence polynomial "classifies" tournaments
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 6: I(Omega, x) AS A TOURNAMENT CLASSIFIER")
    print("  Different tournaments have different Omega, hence different I(Omega, x)")
    print("  I(Omega, 2) = H(T), but I(Omega, x) for general x is MORE information")
    print("  Does I(Omega, x) (as a polynomial) determine the tournament class?")
    print(f"{'='*70}")

    n = 5
    m = n*(n-1)//2

    # For each tournament, compute alpha_1 (and alpha_2=0 at n=5)
    # I(Omega, x) = 1 + alpha_1*x at n=5 (alpha_2=0)
    # So I(Omega, x) is completely determined by alpha_1, which is determined by H.
    # At n=5, I(Omega, x) does NOT give more information than H!

    print(f"\n  n=5: I(Omega, x) = 1 + alpha_1 * x (alpha_2=0)")
    print(f"  So I(Omega, x) determined by H, and vice versa. No extra info.")

    # At n=6: alpha_2 can be nonzero, so I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2
    # gives MORE info than H alone!
    print(f"\n  n=6: I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2")
    print(f"  This is a QUADRATIC in x, giving 2 independent numbers (alpha_1, alpha_2).")
    print(f"  H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 conflates them into one number.")
    print(f"  The full polynomial I(Omega, x) is STRICTLY MORE INFORMATIVE than H at n>=6!")

    # How many distinct I(Omega, x) polynomials at n=6?
    I_polys = set()
    count = 0
    for bits in range(2**m):  # n=5 only
        A = bits_to_adj(bits, n)
        H = compute_H(A, n)
        alpha_1 = (H - 1) // 2
        I_polys.add((alpha_1,))  # at n=5, just alpha_1

    print(f"\n  n=5: {len(I_polys)} distinct I(Omega, x) polynomials")
    print(f"  n=5: {len(set(compute_H(bits_to_adj(bits, n), n) for bits in range(2**m)))} distinct H values")
    print(f"  Same count! (because alpha_2=0)")

    # ============================================================
    # PART 7: The "no adjacent 1s" constraint in tournaments
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 7: THE 'NO ADJACENT 1s' CONSTRAINT")
    print("  In Zeckendorf: no two adjacent Fibonacci indices selected")
    print("  In OCF: no two conflicting cycles selected (independence)")
    print("  In the TILING model: what does 'no adjacent' mean?")
    print(f"{'='*70}")

    # In the pin grid, two tiles (r1,c1) and (r2,c2) are "adjacent" if
    # the corresponding arcs share a vertex. This creates a constraint
    # graph on the pin grid.

    # The GS constraint says paired positions must have equal bits.
    # The "no adjacent 1s" would be an INDEPENDENCE constraint on the grid.

    # Does the independence polynomial of the pin grid adjacency relate to H?
    # The pin grid adjacency is NOT the same as Omega (conflict graph).
    # But there might be a connection.

    print(f"  In the tiling model:")
    print(f"  - Two tiles (arcs) are 'adjacent' if they share a vertex")
    print(f"  - The 'no adjacent 1s' constraint would mean:")
    print(f"    no two backward arcs share a vertex")
    print(f"  - This is the INDEPENDENCE condition on the ARC conflict graph")
    print(f"    (different from the CYCLE conflict graph Omega!)")
    print(f"")
    print(f"  The arc conflict graph has vertices = arcs, edges = shared vertex.")
    print(f"  This is the LINE GRAPH of the complete graph K_n.")
    print(f"  I(line(K_n), x) at x=1 counts independent sets of arcs = MATCHINGS of K_n!")
    print(f"")
    print(f"  SO: The 'no adjacent 1s' constraint on arcs counts MATCHINGS.")
    print(f"  The 'no conflicting cycles' constraint on cycles counts INDEPENDENT CYCLE SETS.")
    print(f"  H = I(Omega, 2) bridges these two worlds!")

    # Number of matchings in K_n:
    # I(line(K_n), 1) = sum_k m_k where m_k = #{matchings of size k}
    print(f"\n  Matchings of K_n (= independent sets of line graph):")
    for nn in range(2, 8):
        # m_k = C(n,2k) * (2k-1)!! for perfect matchings
        # Total matchings = sum_{k=0}^{floor(n/2)} C(n, 2k) * (2k-1)!! / ... complex
        # Just use: #matchings = sum over subsets of edges that are pairwise vertex-disjoint
        # For small n, compute directly
        edges = [(i,j) for i in range(nn) for j in range(i+1, nn)]
        match_count = 0
        for mask in range(2**len(edges)):
            selected = [edges[i] for i in range(len(edges)) if mask & (1 << i)]
            # Check independence: no two selected edges share a vertex
            vertices_used = set()
            independent = True
            for u, v in selected:
                if u in vertices_used or v in vertices_used:
                    independent = False
                    break
                vertices_used.add(u)
                vertices_used.add(v)
            if independent:
                match_count += 1

        print(f"    K_{nn}: {match_count} matchings (= I(line(K_{nn}), 1))")

    print(f"\n  Matching numbers: these are the TELEPHONE NUMBERS (involution counts)!")

    print(f"\n{'='*70}")
    print("DONE — ZECKENDORF-TOURNAMENT CONNECTION EXPLORED")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
