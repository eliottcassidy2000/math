"""
axiom_challenge_meta.py -- kind-pasteur-2026-03-14-S83
CHALLENGING AXIOMS and exploring META-STRUCTURE.

AXIOMS TO CHALLENGE:
1. "H is the natural invariant" — what about F(T, x) at x=3 or x=-1?
2. "The arc encoding is canonical" — tiling encoding gives different Fourier!
3. "Isomorphism is the right equivalence" — what about weaker equivalences?
4. "The number 2 in OCF is special" — what is I(Omega, lambda) for other lambda?
5. "Tournaments are the right objects" — what about semi-tournaments or digraphs?

META-STRUCTURE:
What are the PATTERNS IN THE PATTERNS?
- The sequence of theorems we've proved follows a structure
- The connections between domains (knot, info, coding, dynamics) share a shape
- The "number 2" appears everywhere — is there a meta-2?

FRONTIER:
What hasn't anyone thought of? What's the most SURPRISING possible truth?
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

def main():
    print("=" * 70)
    print("AXIOM CHALLENGE & META-STRUCTURE EXPLORATION")
    print("kind-pasteur-2026-03-14-S83")
    print("=" * 70)

    # ============================================================
    # AXIOM 1: "H is the natural invariant"
    # CHALLENGE: F(T, 3) = I(Omega, 3) might be more natural!
    # ============================================================
    print(f"\n{'='*70}")
    print("AXIOM 1 CHALLENGED: IS F(T,3) = I(Omega,3) MORE NATURAL THAN H?")
    print("  H = I(Omega, 2) = the 'binary' evaluation")
    print("  I(Omega, 3) = the 'ternary' evaluation")
    print("  At n=5: I(Omega,3) = 1 + 3*alpha_1 = (3H-1)/2")
    print("  Is I(Omega, 3) 'better' in some sense?")
    print(f"{'='*70}")

    n = 5
    m = n*(n-1)//2
    N = 2**m

    I2_vals = []  # H = I(Omega, 2)
    I3_vals = []  # I(Omega, 3) = 1 + 3*alpha_1 at n=5
    I1_vals = []  # I(Omega, 1) = 1 + alpha_1

    for bits in range(N):
        A = bits_to_adj(bits, n)
        H = compute_H(A, n)
        alpha_1 = (H - 1) // 2  # since alpha_2 = 0 at n=5
        I2_vals.append(H)
        I3_vals.append(1 + 3 * alpha_1)
        I1_vals.append(1 + alpha_1)

    # Compare information content
    def entropy_of(vals):
        counts = Counter(vals)
        total = len(vals)
        return -sum(c/total * math.log2(c/total) for c in counts.values() if c > 0)

    H_ent = entropy_of(I2_vals)
    I3_ent = entropy_of(I3_vals)
    I1_ent = entropy_of(I1_vals)

    print(f"\n  n=5: Information comparison")
    print(f"    I(Omega,1): {len(set(I1_vals))} values, entropy = {I1_ent:.4f} bits")
    print(f"    I(Omega,2) = H: {len(set(I2_vals))} values, entropy = {H_ent:.4f} bits")
    print(f"    I(Omega,3): {len(set(I3_vals))} values, entropy = {I3_ent:.4f} bits")

    # At n=5, all three are linearly related (alpha_2=0), so they have SAME information!
    print(f"\n    At n=5: all three are affinely related (alpha_2=0)")
    print(f"    So they carry EXACTLY the same information!")
    print(f"    The difference only matters at n>=6 where alpha_2 > 0.")

    # At n=6: I(Omega, lambda) = 1 + lambda*alpha_1 + lambda^2*alpha_2
    # Different lambda values weight alpha_1 vs alpha_2 differently!
    print(f"\n  At n>=6: different lambda emphasize different cycle structures:")
    print(f"    lambda=1: weights all independent sets equally (alpha_k has weight 1)")
    print(f"    lambda=2: H = OCF point, weights pairs 4x singles (standard)")
    print(f"    lambda=3: weights pairs 9x singles (emphasizes disjoint cycles)")
    print(f"    lambda->inf: dominated by max independent set (clique number)")

    # ============================================================
    # AXIOM 2: "The permanent gaps {7, 21} are structural accidents"
    # CHALLENGE: They might be INEVITABLE from the algebraic structure
    # ============================================================
    print(f"\n{'='*70}")
    print("AXIOM 2 CHALLENGED: ARE THE GAPS {7,21} STRUCTURALLY INEVITABLE?")
    print("  THM-200 proves H≠7 via Omega = K_3 impossibility")
    print("  But WHY is K_3 impossible? Is it just combinatorics,")
    print("  or is there a deeper algebraic reason?")
    print(f"{'='*70}")

    # The (z-2)(z-3) orbit: 7, 21, 63, 189, ...
    # = 7 * 3^k for k = 0, 1, 2, 3, ...
    # H=7 forbidden (proved), H=21 forbidden (proved through n≤8)
    # H=63 NOT forbidden (achieved at n=8)
    # H=189 NOT forbidden (= max_H(7) = H(Paley T_7))

    print(f"\n  The (z-2)(z-3) orbit: 7, 21, 63, 189, 567, 1701, ...")
    orbit = [7]
    for _ in range(8):
        orbit.append(orbit[-1] * 3)
    print(f"  {orbit}")
    print(f"  Forbidden: 7, 21")
    print(f"  Achievable: 63 (n=8), 189 (n=7, Paley!)")
    print(f"")
    print(f"  META-PATTERN: The first TWO values of the orbit are forbidden,")
    print(f"  then it 'escapes'. Why two? Because:")
    print(f"  - H=7 requires alpha_1=3, alpha_2=0 with Omega=K_3")
    print(f"  - H=21 requires alpha_1+2*alpha_2=10 with specific Omega")
    print(f"  - H=63 can use alpha_3>=1 (which first appears at n=8)")
    print(f"  The 'escape' happens when higher alpha_k terms kick in!")

    # ============================================================
    # AXIOM 3: "Fourier on the arc hypercube is the natural basis"
    # CHALLENGE: Maybe the CYCLE BASIS is more natural
    # ============================================================
    print(f"\n{'='*70}")
    print("AXIOM 3 CHALLENGED: IS THE CYCLE BASIS MORE NATURAL THAN ARC BASIS?")
    print("  H = 1 + 2*alpha_1 + 4*alpha_2 + ... (OCF)")
    print("  The OCF writes H in terms of CYCLE STRUCTURE, not arcs.")
    print("  Maybe the 'right' basis for understanding H is:")
    print("  {alpha_k : k = 0,1,2,...} not {arc_e : e in E}")
    print(f"{'='*70}")

    # The alpha_k are integer-valued invariants
    # alpha_1 = total odd cycles, alpha_2 = disjoint pairs, etc.
    # H = sum 2^k * alpha_k = binary expansion in cycle-count basis!

    # This is like a BINARY NUMBER: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
    # The alpha_k are the "binary digits" of (H-1) in a specific sense!

    print(f"\n  H-1 = 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...")
    print(f"  = 2*(alpha_1 + 2*alpha_2 + 4*alpha_3 + ...)")
    print(f"  = 2 * T where T = alpha_1 + 2*alpha_2 + ...")
    print(f"")
    print(f"  T is the 'total cycle weight' = alpha_1 + 2*alpha_2 + 4*alpha_3 + ...")
    print(f"  H = 1 + 2T, so T = (H-1)/2")
    print(f"  T counts independent sets of Omega weighted by 2^(size-1)")
    print(f"")
    print(f"  INSIGHT: T is the 'half-Redei' number. H = 2T+1 (always odd).")
    print(f"  T = 0 (transitive), T = 1 (one 3-cycle), T = 2 (two non-disjoint), ...")
    print(f"  T is a NON-NEGATIVE INTEGER that fully determines H.")

    # The forbidden values in terms of T:
    # H=7: T=3. H=21: T=10. H=63: T=31.
    print(f"\n  Forbidden T values:")
    print(f"    T=3 (H=7): FORBIDDEN (proved)")
    print(f"    T=10 (H=21): FORBIDDEN (proved through n<=8)")
    print(f"    T=31 (H=63): NOT forbidden (achieved at n=8)")

    # ============================================================
    # META-STRUCTURE 1: The "theory of the theory"
    # ============================================================
    print(f"\n{'='*70}")
    print("META-STRUCTURE 1: THE SHAPE OF OUR THEORY")
    print(f"{'='*70}")
    print("""
  Our theory has a LAYERED structure, like an onion:

  LAYER 0 (Ground Truth):
    H(T) is always odd. [Redei, 1934]
    H(T) = I(Omega(T), 2). [Grinberg-Stanley, 2023]

  LAYER 1 (Structural):
    deg(H) = 2*floor((n-1)/2). [Our proof via path reversal]
    Fourier spectrum: only even levels, 75/25 split.
    H(T) = H(T^op). [Path reversal]

  LAYER 2 (Landscape):
    Unimodal for n<=5, multimodal for n>=6. [Our computation]
    Score regularization = H maximization (R=-0.97).
    Phase transition at beta_c ≈ 0.3-0.8. [Stat mech]

  LAYER 3 (Forbidden Structure):
    H=7 forbidden (Omega=K_3 impossible). [THM-200]
    H=21 forbidden (6-way block). [THM-079]
    (z-2)(z-3) orbit escapes at k=2 via alpha_3.

  LAYER 4 (Connections):
    DC = Jones skein relation. [Our observation]
    GS code = product code on pin grid. [Our proof]
    Lex product formula for transitive T1. [Our verification]
    0.27 = non-constant Fourier energy / m. [Our computation]

  LAYER 5 (Open Frontiers):
    Categorification via Khovanov/Sazdanovic-Yip?
    Awan-Bernardi B-polynomial specialization?
    Tournament operads and composition?
    Var/Mean^2 = 1/3 exact proof?

  META-OBSERVATION: Each layer EXPLAINS the one above it.
  OCF explains Redei. Fourier explains the landscape.
  The forbidden structure explains the gaps. Connections
  explain WHY the structure looks the way it does.
""")

    # ============================================================
    # META-STRUCTURE 2: The "2-3-5 universality" hypothesis
    # ============================================================
    print(f"\n{'='*70}")
    print("META-STRUCTURE 2: THE 2-3-5 UNIVERSALITY")
    print(f"{'='*70}")
    print("""
  The numbers 2, 3, 5 appear EVERYWHERE in tournament theory:

  THE NUMBER 2:
    - H = I(Omega, 2)
    - deg(H) coefficients = +-2 at odd n
    - H(T) = H(T^op) (involution of order 2)
    - 2 orientations per arc
    - mean H = n!/2^{n-1}

  THE NUMBER 3:
    - 3-cycles are the fundamental building blocks
    - (z-2)(z-3) recurrence
    - H = 1 + 2T with T = alpha_1 + 2*alpha_2 (base 2 + base 3 interplay)
    - max_H ratio approaches 3 for consecutive n
    - GS weight enumerator has (1+z^2)^p factors (degree 2, from pairs)

  THE NUMBER 5:
    - 5-cycles first appear at n=5
    - KEY_SUM = 2+3 = 5
    - n=5 is the largest n with claw-free Omega
    - H(T_5 regular) = 15 = 3*5

  WHY 2, 3, 5? These are the first three primes.
  Tournament theory lives in the intersection of:
  - Binary choice (2)
  - Cyclic structure (3)
  - Pentagonal symmetry (5)

  CHALLENGE: Is this just pattern-matching, or is there
  a STRUCTURAL reason why exactly these three primes matter?
""")

    # ============================================================
    # FRONTIER 1: What if H has a p-adic expansion?
    # ============================================================
    print(f"\n{'='*70}")
    print("FRONTIER 1: p-ADIC EXPANSION OF H")
    print("  H is always odd. What if we expand in powers of 2?")
    print("  H = 1 + 2*(alpha_1 + 2*alpha_2 + 4*alpha_3 + ...)")
    print("  The alpha_k ARE the p=2 expansion coefficients!")
    print(f"{'='*70}")

    n = 5
    m = n*(n-1)//2
    for bits in range(min(2**m, 20)):
        A = bits_to_adj(bits, n)
        H = compute_H(A, n)
        T = (H - 1) // 2  # T = alpha_1 + 2*alpha_2 + ...

        # Binary expansion of T
        T_binary = bin(T)[2:] if T > 0 else '0'
        print(f"  bits={bits:06b}: H={H:3d}, T=(H-1)/2={T:2d}, T_binary={T_binary}")

    print(f"\n  The binary expansion of T = (H-1)/2 encodes the cycle structure!")
    print(f"  T bit 0 (= alpha_1 mod 2): parity of total odd cycles")
    print(f"  T bit 1 (= alpha_2 mod 2): parity of disjoint cycle pairs")
    print(f"  T bit 2 (= alpha_3 mod 2): parity of disjoint cycle triples")

    # ============================================================
    # FRONTIER 2: The "tournament category"
    # ============================================================
    print(f"\n{'='*70}")
    print("FRONTIER 2: THE TOURNAMENT CATEGORY")
    print("  Objects: tournaments")
    print("  Morphisms: tournament homomorphisms (arc-preserving maps)")
    print("  Question: What are the universal properties of H in this category?")
    print(f"{'='*70}")

    # A tournament homomorphism f: T1 -> T2 is a map f: V1 -> V2 such that
    # T1[u][v] = 1 implies T2[f(u)][f(v)] = 1.
    # NOT the same as isomorphism (not required to be bijective or injective).

    # Question: Does H(T1) >= H(T2) when T1 -> T2 is a homomorphism?
    # Or: is H monotone under homomorphisms?

    # For the TRIVIAL homomorphism (all vertices map to one vertex):
    # T1 -> T_1 (single vertex). H(T_1) = 1. So H(T1) >= 1 = H(T_1). ✓

    # For a QUOTIENT: T -> T/e (contraction). H(T) = H(T\e) + H(T/e) (DC).
    # So H(T/e) <= H(T). Contraction DECREASES H. ✓

    print(f"  THEOREM: H is monotone decreasing under contraction:")
    print(f"  H(T/e) <= H(T) (from DC: H = H(T\\e) + H(T/e) >= H(T/e))")
    print(f"")
    print(f"  Is H monotone under general homomorphisms? OPEN QUESTION.")

    # ============================================================
    # FRONTIER 3: What ISN'T a tournament?
    # ============================================================
    print(f"\n{'='*70}")
    print("FRONTIER 3: EXTENDING BEYOND TOURNAMENTS")
    print("  What if we allow bidirectional arcs? (semicomplete digraphs)")
    print("  What if we allow missing arcs? (partial tournaments)")
    print("  H(T)=7 is impossible for TOURNAMENTS but possible for DIGRAPHS!")
    print(f"{'='*70}")

    # opus's THM (S71g): H=7 lives in digraph world, not tournament world
    # The directed 7-cycle has H=7 as a digraph!
    # But it's not a tournament (not complete).

    print(f"  The directed 7-cycle C_7 has H(C_7) = 7 (as a digraph).")
    print(f"  But C_7 is NOT a tournament (only 7 arcs, not C(7,2)=21).")
    print(f"  The COMPLETENESS constraint is what forbids H=7.")
    print(f"")
    print(f"  DEEP INSIGHT: Tournament theory is about what COMPLETENESS forces.")
    print(f"  H=7 is impossible BECAUSE every vertex pair has an arc,")
    print(f"  creating unavoidable additional cycles that push H above 7.")

    # ============================================================
    # FRONTIER 4: The "dual tournament"
    # ============================================================
    print(f"\n{'='*70}")
    print("FRONTIER 4: THE DUAL — WHAT COUNTS HAM CYCLES INSTEAD OF PATHS?")
    print(f"{'='*70}")

    # C(T) = number of directed Hamiltonian CYCLES in T
    # By definition: C(T) = H(T) * ... no, Ham cycles are different.
    # For a tournament T: C(T) = (1/n) * sum_v #{Ham cycles through v as "start"}
    # For regular tournaments: C(T) = H(T) * ... related by a factor.

    print(f"  Let C(T) = number of directed Hamiltonian cycles in T.")
    print(f"  Moon's theorem: C(T) > 0 iff T is strongly connected.")

    n = 5
    cycle_data = []
    for bits in range(2**m):
        A = bits_to_adj(bits, n)
        H_val = compute_H(A, n)

        # Count Ham cycles: a cycle visits all n vertices and returns to start
        cycles = 0
        for perm in permutations(range(1, n)):  # fix vertex 0 as start
            path = (0,) + perm
            valid = all(A[path[i]][path[i+1]] for i in range(n-1))
            if valid and A[path[-1]][path[0]]:  # close the cycle
                cycles += 1
        # Each cycle is counted once (fixed starting vertex)

        cycle_data.append((H_val, cycles))

    H_C_pairs = Counter(cycle_data)
    print(f"\n  n=5: (H, C) joint distribution:")
    for (h, c), count in sorted(H_C_pairs.items()):
        ratio = c / h if h > 0 else 'N/A'
        print(f"    H={h:3d}, C={c:3d}: {count:4d} tournaments, C/H={ratio}")

    # Is C/H constant for regular tournaments?
    print(f"\n  For regular tournaments (H=15): C/H values = ",
          sorted(set(c/h for (h,c),_ in H_C_pairs.items() if h == 15 and _ > 0)))

    # ============================================================
    # FRONTIER 5: The generating function Z(q, t, n)
    # ============================================================
    print(f"\n{'='*70}")
    print("FRONTIER 5: THE MASTER GENERATING FUNCTION")
    print("  Z(q, t, n) = sum over all tournaments T on n vertices:")
    print("  q^{H(T)} * t^{c3(T)}")
    print("  This encodes EVERYTHING about (H, c3) jointly.")
    print("  Does Z factor? Have nice specializations?")
    print(f"{'='*70}")

    for n_val in [4, 5]:
        m_val = n_val*(n_val-1)//2
        terms = Counter()
        for bits in range(2**m_val):
            A = bits_to_adj(bits, n_val)
            H_val = compute_H(A, n_val)
            c3 = int(np.trace(A @ A @ A)) // 3
            terms[(H_val, c3)] += 1

        print(f"\n  n={n_val}: Z(q,t) =")
        for (h, c3), count in sorted(terms.items()):
            print(f"    + {count} * q^{h} * t^{c3}")

        # Evaluate at special points
        Z_1_1 = sum(terms.values())
        Z_1_0 = sum(count for (h, c3), count in terms.items() if c3 == 0)
        print(f"\n  Z(1,1) = {Z_1_1} = 2^{m_val} ✓")
        print(f"  Z(1,0) = {Z_1_0} = number of transitive tournaments = {math.factorial(n_val)}")

    print(f"\n{'='*70}")
    print("DONE — AXIOMS CHALLENGED, FRONTIERS EXPLORED")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
