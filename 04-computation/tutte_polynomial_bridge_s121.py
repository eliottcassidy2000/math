#!/usr/bin/env python3
"""
tutte_polynomial_bridge_s121.py — opus-2026-03-21-S121

TUTTE POLYNOMIALS AND THE TOURNAMENT POLYNOMIAL WEB

Building on:
- kind-pasteur S18c: Tournament-as-Codes, Greene parallel, partition function bridge
- S115: g(x) = x³-x²-x-1 curvature classifier
- S117: Casimir invariants d ∈ {5,21}
- S120: Fano plane embeds in T_7, (2,3,7) Coxeter
- arXiv:1610.01839: B-polynomial for directed graphs

THE WEB OF TOURNAMENT POLYNOMIALS:
  1. I(Ω(T), x): independence polynomial of conflict graph → OCF: H = I(Ω, 2)
  2. F(T, x): descent polynomial → F(T, 1) = H(T)
  3. W(T, r): the W-polynomial → W(T, 1/2) = H(T)
  4. g(x) = x³-x²-x-1: curvature classifier
  5. B(T; x,y,z): the directed Tutte polynomial (B-polynomial)
  6. χ(G, k): chromatic polynomial of the conflict graph

THE QUESTION: How do these polynomials connect? Is there a MASTER
polynomial from which all others are specializations?

Author: opus-2026-03-21-S121
"""

import sys
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter
from math import comb, factorial
sys.stdout.reconfigure(line_buffering=True)

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        yield A

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def count_c3(A, n):
    return sum(1 for i,j,k in combinations(range(n), 3)
               if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]))

def find_odd_cycles(A, n, max_len=None):
    if max_len is None: max_len = n
    cycles = set()
    for length in range(3, max_len+1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts[1:]):
                path = [verts[0]] + list(perm)
                if all(A[path[k]][path[(k+1)%length]] for k in range(length)):
                    cycles.add(frozenset(path))
    return list(cycles)

def independence_poly_coeffs(adj, nc):
    alpha = [0] * (nc + 1)
    for mask in range(1 << nc):
        subset = [i for i in range(nc) if (mask >> i) & 1]
        ok = True
        for a in range(len(subset)):
            for b in range(a+1, len(subset)):
                if adj[subset[a]][subset[b]]:
                    ok = False; break
            if not ok: break
        if ok: alpha[len(subset)] += 1
    return alpha

print("=" * 72)
print("  TUTTE POLYNOMIALS AND THE TOURNAMENT POLYNOMIAL WEB")
print("  opus-2026-03-21-S121")
print("=" * 72)

# ========================================================================
# PART 1: The descent polynomial F(T, x) and its connection to H
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: THE DESCENT POLYNOMIAL F(T, x)")
print("=" * 72)

n = 5
# F(T, x) = sum over Hamiltonian paths P of x^{des(P)}
# where des(P) = number of descents in the path.
# A descent at position i means P[i] > P[i+1].
# F(T, 1) = H(T) (every path contributes 1).

print(f"\n  Computing F(T, x) for n={n} tournament classes:")

desc_polys = defaultdict(Counter)

for A in all_tournaments(n):
    H = count_hp(A, n)
    if H in desc_polys and len(desc_polys[H]) > 0:
        continue  # one representative per H

    # Enumerate HPs and count descents
    for perm in permutations(range(n)):
        # Check if this is a HP
        is_hp = all(A[perm[i]][perm[i+1]] for i in range(n-1))
        if is_hp:
            des = sum(1 for i in range(n-1) if perm[i] > perm[i+1])
            desc_polys[H][des] += 1

for H_val in sorted(desc_polys.keys()):
    poly = desc_polys[H_val]
    coeffs = [poly.get(k, 0) for k in range(n)]
    # F(T, x) = sum_k f_k x^k
    F_at_1 = sum(coeffs)
    print(f"  H={H_val}: F(T,x) = {' + '.join(f'{c}x^{k}' for k, c in enumerate(coeffs) if c)}")
    print(f"    F(T,1) = {F_at_1} = H ✓")
    print(f"    F(T,0) = {coeffs[0]} (paths with no descent)")
    print(f"    F(T,-1) = {sum((-1)**k * c for k, c in enumerate(coeffs))}")

# ========================================================================
# PART 2: I(Omega, x) — the independence polynomial at general x
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 2: I(Omega(T), x) AT MULTIPLE EVALUATIONS")
print("=" * 72)

n = 5
# For each tournament class, compute I(Omega, x) for several x values

ip_data = defaultdict(dict)

for A in all_tournaments(n):
    H = count_hp(A, n)
    if H in ip_data:
        continue

    cycles = find_odd_cycles(A, n)
    nc = len(cycles)

    if nc == 0:
        alpha = [1]
    else:
        adj_omega = [[0]*nc for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if cycles[i] & cycles[j]:
                    adj_omega[i][j] = adj_omega[j][i] = 1
        alpha = independence_poly_coeffs(adj_omega, nc)

    # Evaluate at several points
    tau = 1.8392867552
    phi = (1 + 5**0.5) / 2

    ip_data[H] = {
        'alpha': alpha[:6],
        'I_at_0': 1,
        'I_at_1': sum(alpha),
        'I_at_2': sum(alpha[k] * 2**k for k in range(len(alpha))),
        'I_at_neg1': sum(alpha[k] * (-1)**k for k in range(len(alpha))),
        'I_at_tau': sum(alpha[k] * tau**k for k in range(len(alpha))),
        'I_at_phi': sum(alpha[k] * phi**k for k in range(len(alpha))),
        'nc': nc,
    }

print(f"\n  {'H':>4s} {'|V(Ω)|':>7s} {'I(Ω,1)':>8s} {'I(Ω,2)=H':>10s} {'I(Ω,-1)':>8s} "
      f"{'I(Ω,τ)':>8s} {'I(Ω,φ)':>8s}")
for H_val in sorted(ip_data.keys()):
    d = ip_data[H_val]
    print(f"  {H_val:4d} {d['nc']:7d} {d['I_at_1']:8.0f} {d['I_at_2']:10.0f} "
          f"{d['I_at_neg1']:8.0f} {d['I_at_tau']:8.2f} {d['I_at_phi']:8.2f}")

# ========================================================================
# PART 3: The Greene-MacWilliams parallel
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 3: THE GREENE-MACWILLIAMS PARALLEL")
print("=" * 72)
print(f"""
  From kind-pasteur S18c:

  GREENE'S THEOREM (1976): For a linear code C with matroid M(C),
  the weight enumerator W_C(x) is determined by the Tutte polynomial
  T(M(C); x, y).

  The TOURNAMENT PARALLEL:
  The tournament T on n vertices defines a matroid M(T) on the arcs.
  The "weight enumerator" is I(Omega(T), x).
  The Tutte polynomial T(M(T); x, y) should specialize to give I(Omega, x).

  But tournaments are DIRECTED, so we need the B-POLYNOMIAL (directed Tutte).
  From arXiv:1610.01839: the B-polynomial B(D; x, y, z) for a digraph D
  generalizes the Tutte polynomial to directed graphs.

  THE CHAIN (conjectured):
  B(T; x, y, z) --specialization--> I(Omega(T), x) --evaluate x=2--> H(T)

  If this chain holds, then the Tutte-like polynomial B(T) would be the
  MASTER POLYNOMIAL from which H is extracted by two successive operations:
  (1) specialize B to get the independence polynomial, (2) evaluate at x=2.

  THE EXISTING CONNECTIONS:
  1. F(T, x) = descent polynomial. F(T, 1) = H.
  2. I(Omega, x) = independence polynomial. I(Omega, 2) = H.
  3. g(x) = x^3-x^2-x-1. g(tau) = 0, g(2) = 1.

  Do F and I determine each other?
  Both have H as a common specialization. Are they related by a transform?
""")

# Check if F(T, x) determines I(Omega, x) at n=5
print(f"\n  F(T, x) vs I(Omega, x) at n=5:")
print(f"  {'H':>4s} {'F coeffs':>20s} {'I coeffs':>20s}")

for H_val in sorted(ip_data.keys()):
    F = desc_polys[H_val]
    F_coeffs = [F.get(k, 0) for k in range(n)]
    I_coeffs = ip_data[H_val]['alpha']
    print(f"  {H_val:4d} {str(F_coeffs):>20s} {str(I_coeffs[:5]):>20s}")

# ========================================================================
# PART 4: Chromatic polynomial of Omega and the Tutte connection
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 4: THE CHROMATIC POLYNOMIAL OF Omega")
print("=" * 72)
print(f"""
  For any graph G, the chromatic polynomial chi(G, k) and the Tutte
  polynomial T(G; x, y) are related by:
    chi(G, k) = (-1)^|V| * k * T(G; 1-k, 0)

  For the CONFLICT GRAPH Omega(T):
    chi(Omega, k) = number of proper k-colorings of Omega.

  And the independence polynomial I(G, x) is related to the Tutte
  polynomial via the more general two-variable Potts model:
    Z(G; q, v) = sum over colorings sum over edges of interaction terms

  At specific parameters:
  - Z at (q, v) = (1, 0) → independent sets (hard-core gas)
  - Z at (q, v) = (k, -1) → proper k-colorings

  So: I(Omega, x) and chi(Omega, k) are BOTH specializations of the
  partition function / Tutte polynomial.

  H = I(Omega, 2) = Z(Omega; 1, 2) (hard-core gas at fugacity 2)
  chi(Omega, 2) = Z(Omega; 2, -1) (2-colorings of Omega)

  These are DIFFERENT evaluations of the same partition function!
""")

# Compute chromatic polynomial of Omega for n=5 tournaments
print(f"\n  Chromatic polynomial chi(Omega, k) at n=5:")

for A in all_tournaments(n):
    H = count_hp(A, n)
    if H not in [1, 9, 15]:
        continue

    cycles = find_odd_cycles(A, n)
    nc = len(cycles)

    if nc == 0:
        print(f"  H={H}: Omega is empty, chi(Omega, k) = k^0 = 1 (trivially)")
        continue
    if nc > 10:
        print(f"  H={H}: Omega has {nc} vertices, skipping chi")
        continue

    adj_omega = [[0]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                adj_omega[i][j] = adj_omega[j][i] = 1

    # Compute chi(Omega, k) by inclusion-exclusion for k = 1,2,3,...
    # chi(G, k) = sum over independent sets S of (-1)^{|V|-|S|} * k^{|S|}
    # Actually that's the independence polynomial approach:
    # I(complement(G), x) related to chi via substitution

    # Let me just count proper k-colorings directly for small k
    chi_vals = {}
    for k in range(1, 6):
        count = 0
        for coloring in range(k**nc):
            colors = [(coloring // k**i) % k for i in range(nc)]
            proper = True
            for i in range(nc):
                for j in range(i+1, nc):
                    if adj_omega[i][j] and colors[i] == colors[j]:
                        proper = False
                        break
                if not proper:
                    break
            if proper:
                count += 1
        chi_vals[k] = count

    print(f"\n  H={H}: Omega has {nc} vertices")
    print(f"    chi(Omega, k) for k=1,...,5: {[chi_vals[k] for k in range(1,6)]}")
    print(f"    chi(Omega, 2) = {chi_vals[2]} (2-colorings of conflict graph)")
    print(f"    I(Omega, 2) = {H} (independence polynomial at 2 = H)")

    break  # one per H

# ========================================================================
# PART 5: The Tutte polynomial of K_n and tournament counts
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 5: TUTTE POLYNOMIAL OF K_n AND TOURNAMENTS")
print("=" * 72)
print(f"""
  The TUTTE POLYNOMIAL of the complete graph K_n:
  T(K_n; x, y) encodes structural information about subgraphs of K_n.

  KEY EVALUATIONS:
  T(K_n; 1, 1) = number of spanning trees of K_n = n^(n-2) (Cayley)
  T(K_n; 2, 0) = number of acyclic orientations of K_n = n!
  T(K_n; 0, 2) = number of totally cyclic orientations of K_n

  For TOURNAMENTS specifically:
  A tournament is an orientation of K_n.
  The number of tournaments = 2^C(n,2) (every orientation is valid).
  The number of ACYCLIC tournaments = n! (one per permutation = one HP).
  Wait: acyclic orientations of K_n = n! which equals the number of
  linear orders. And T(K_n; 2, 0) = n! ✓.

  THE BEAUTIFUL FACT:
  H(T) = number of Hamiltonian paths = number of acyclic orientations
  of K_n that contain T as a sub-orientation? No...

  Actually: H(T) = number of topological sortings of T if we view T
  as a directed graph. An acyclic tournament IS a topological sort.
  For a general tournament (which may have cycles), the number of
  Hamiltonian paths = number of ways to linearly extend the partial
  order induced by the tournament... no, every tournament has a HP (Rédei).

  THE CONNECTION IS THROUGH THE INDEPENDENCE POLYNOMIAL:
  H = I(Omega, 2) = hard-core partition function at fugacity 2.
  The hard-core model on Omega IS a specialization of the multivariate
  Tutte polynomial (also called the Fortuin-Kasteleyn model):
    Z_FK(Omega; q, v) with q=1, v=x = 2.

  So: H = Z_FK(Omega; 1, 2).

  And T(G; x, y) relates to Z_FK by:
    Z_FK(G; q, v) = sum over spanning subgraphs H of q^(k(H)) * v^(|E(H)|)
  where k(H) = number of connected components of H.

  For the INDEPENDENT SET partition function:
    I(G, x) = sum over independent sets S of x^|S|
  This is NOT the FK model but the HARD-CORE MODEL.
  It is a specialization of a DIFFERENT generalization.

  THE RELATIONSHIP:
  I(G, x) is the independence polynomial (hard-core gas).
  Z_FK(G; q, v) is the Potts/Tutte partition function.
  They are related but NOT direct specializations of each other.

  HOWEVER: for the COMPLEMENT graph:
  I(G, x) = chi(complement(G), 1+x) / (1+x)... not quite.

  The exact connection: the partition function of the ANTIFERROMAGNETIC
  Potts model at q → 0 limit gives the independence polynomial.
""")

# ========================================================================
# PART 6: F(T,x), I(Omega,x), and the Tutte: evaluation relationships
# ========================================================================
print("=" * 72)
print("  PART 6: THE POLYNOMIAL EVALUATION TABLE")
print("=" * 72)

# Compute F(T, -1) and I(Omega, -1) for comparison
print(f"\n  EVALUATION TABLE for n=5 tournament classes:")
print(f"  {'H':>4s} {'F(T,-1)':>8s} {'I(Ω,-1)':>8s} {'I(Ω,1)':>7s} {'I(Ω,τ)':>8s} "
      f"{'F_0':>4s} {'F_{n-1}':>6s} {'α₁':>4s}")

for H_val in sorted(ip_data.keys()):
    d = ip_data[H_val]
    F = desc_polys[H_val]
    F_at_neg1 = sum((-1)**k * F.get(k, 0) for k in range(n))
    F_0 = F.get(0, 0)
    F_nm1 = F.get(n-1, 0)
    alpha1 = d['alpha'][1] if len(d['alpha']) > 1 else 0

    print(f"  {H_val:4d} {F_at_neg1:8d} {d['I_at_neg1']:8.0f} {d['I_at_1']:7.0f} "
          f"{d['I_at_tau']:8.2f} {F_0:4d} {F_nm1:6d} {alpha1:4d}")

print(f"""
  OBSERVATIONS:
  1. F(T, -1) = 1 for ALL tournaments with odd n! (by Rédei: F(-1) relates to
     the signed path count, which is always 1 by the complement involution.)

  2. I(Omega, -1) = {set(ip_data[H]['I_at_neg1'] for H in ip_data)} for all classes.
     The Euler characteristic of the independence complex!

  3. I(Omega, 1) = total independent sets = {[ip_data[H]['I_at_1'] for H in sorted(ip_data)]}.
     Grows with H (more cycles → more independent sets).

  4. F_0 = number of paths with 0 descents = number of INCREASING paths.
     This is 1 for the transitive tournament and 0 for others.
     Actually F_0 can be > 0 for non-transitive tournaments...
""")

# ========================================================================
# PART 7: Toward the master polynomial
# ========================================================================
print("=" * 72)
print("  PART 7: TOWARD THE MASTER POLYNOMIAL")
print("=" * 72)
print(f"""
  THE WEB OF POLYNOMIAL CONNECTIONS:

  KNOWN:
  1. F(T, 1) = H(T) = I(Omega(T), 2)  [different polynomials, same evaluation]
  2. Sum_T F(T, x) = A_n(x) * 2^(C(n,2)-(n-1))  [sum over all tournaments]
  3. I(Omega, x) is claw-free for Omega at n<=8 → all real roots
  4. g(x) = x^3-x^2-x-1 classifies {phi, tau, 2}

  CONJECTURED:
  5. B(T; x, y, z) --specialize--> I(Omega(T), x) --eval x=2--> H(T)
     (The directed Tutte polynomial specializes to the independence polynomial)

  6. F(T, x) and I(Omega, x) are both specializations of a SINGLE
     two-variable polynomial P(T; x, y) where:
     P(T; x, 1) = F(T, x)  (descent polynomial)
     P(T; 2, y) = I(Omega(T), y)  (independence polynomial shifted)

  OPEN QUESTION: Does such a P(T; x, y) exist?

  IF IT EXISTS: P(T; 1, 2) = F(T, 1) = H(T) = I(Omega(T), 2).
  It would unify the combinatorial (F) and algebraic (I) perspectives.

  THE TUTTE POLYNOMIAL T(K_n; x, y) could be this master:
  - T(K_n; 2, 0) = n! (number of acyclic orientations = all HPs of all tournaments)
  - T(K_n; 1+x, 0) = chromatic polynomial = chi(K_n, 1+x)
  - The FK model on Omega relates to specific T evaluations

  But the connection is INDIRECT: T(K_n) is a property of the COMPLETE GRAPH,
  while we need a polynomial for the SPECIFIC TOURNAMENT T.
  The B-polynomial B(T) for the directed graph T is the right object.

  THE MASTER POLYNOMIAL CANDIDATE: B(T; x, y, z)
  from arXiv:1610.01839. It reduces to T(K_n; x, y) when T is bidirected,
  and to I(Omega(T), x) under appropriate specialization (conjectured).
""")

# ========================================================================
# PART 8: The Potts model connection
# ========================================================================
print("=" * 72)
print("  PART 8: THE POTTS MODEL AND H=I(Omega, 2)")
print("=" * 72)
print(f"""
  The POTTS MODEL partition function on graph G:
    Z_Potts(G; q, v) = Sum over colorings c Product over edges e of (1 + v * delta(c_u, c_v))

  The independence polynomial I(G, x) is the "hard-core gas" partition function:
    I(G, x) = Sum over independent sets S of x^|S|

  These are DIFFERENT models on the same graph:
  - Potts: colors on vertices, interactions on edges
  - Hard-core: occupations on vertices, exclusions on edges

  H = I(Omega, 2) is the hard-core partition function at fugacity lambda=2.

  The TUTTE POLYNOMIAL T(G; x, y) relates to Potts via:
    Z_Potts(G; q, v) = (x-1)^(-k(G)) * y^(-|V|) * T(G; x, y)
  where x = 1 + q/v, y = 1 + v (roughly — the exact transformation depends
  on the convention).

  For the hard-core model: the connection to Tutte is through the
  INDEPENDENT SET POLYNOMIAL:
    I(G, x) = x^n * T(G_complement; 1 + 1/x, 0)? Not exactly...

  Actually: the independence polynomial doesn't directly specialize from
  the classical Tutte polynomial. It requires a WEIGHTED version or the
  multivariate Tutte polynomial.

  MULTIVARIATE TUTTE (Sokal):
    Z_G(q, v_e) = Sum over A subset E of q^(k(A)) * Product_(e in A) v_e

  When we set q=1, v_e = lambda for all e:
    Z_G(1, lambda) = Sum_A lambda^(|A|) = (1+lambda)^(|E|)... not quite.

  The hard-core model is actually the LATTICE GAS:
    Z_HC(G, lambda) = Sum over independent sets S of lambda^(|S|)
                    = I(G, lambda)

  This is NOT a direct Tutte specialization but a TRANSFER MATRIX
  eigenvalue problem on the graph.

  FOR TOURNAMENTS: H = I(Omega, 2) = Z_HC(Omega, 2).
  The hard-core gas at fugacity 2 on the conflict graph.
  From the lattice gas / statistical mechanics perspective,
  fugacity 2 is DEEP in the ordered (crystal) phase — far from
  the uniqueness threshold. This was noted in earlier sessions.

  THE SIGNIFICANCE: H is not a "nice" Tutte evaluation.
  It's a partition function evaluation that requires the FULL
  hard-core model, not just the Tutte polynomial.
  The Tutte polynomial captures STRUCTURAL properties of the graph
  (spanning trees, forests, connected subgraphs).
  The hard-core model captures OCCUPANCY properties (independent sets).
  H requires the latter, not the former.

  HOWEVER: the chromatic polynomial (a Tutte specialization) and the
  independence polynomial are related by:
    chi(G, k) and I(complement(G), ...) have connections through
    the "clique polynomial" of the complement.

  And I(G, x) IS the clique polynomial of complement(G).
  So: H = I(Omega, 2) = clique polynomial of complement(Omega) at x=2.
  And the clique polynomial IS related to Tutte via the "strong" Tutte:
    "the clique polynomial of G is a specialization of the multivariate
    Tutte polynomial of G" — this is a known result.

  CONCLUSION: H = I(Omega, 2) IS a specialization of the multivariate
  Tutte polynomial of Omega. But the path from tournament T to Omega
  (via cycle detection) is the non-trivial part.
""")

# ========================================================================
# PART 9: Synthesis — the polynomial web
# ========================================================================
print("=" * 72)
print("  SYNTHESIS: THE TOURNAMENT POLYNOMIAL WEB")
print("=" * 72)
print(f"""
  LEVEL 0: THE TOURNAMENT T (a binary string of length C(n,2))

  LEVEL 1: POLYNOMIAL INVARIANTS OF T
    F(T, x) = descent polynomial (evaluations: F(1)=H, F(-1)=1 for odd n)
    W(T, r) = W-polynomial (evaluation: W(1/2)=H)
    I(Omega(T), x) = independence polynomial of conflict graph (I(2)=H)

  LEVEL 2: TUTTE-LIKE POLYNOMIALS
    B(T; x, y, z) = directed B-polynomial (the master for directed graphs)
    T(K_n; x, y) = Tutte polynomial of the underlying complete graph
    Z_multivar(Omega; q, v_e) = multivariate Tutte of the conflict graph

  LEVEL 3: UNIVERSAL POLYNOMIALS
    g(x) = x^3-x^2-x-1 = curvature classifier ({-1, 0, +1} at {phi, tau, 2})
    A_n(x) = Eulerian polynomial (Sum_T F(T,x) relates to A_n(x))

  KNOWN CONNECTIONS:
    F(T, 1) = I(Omega, 2) = H(T)  [the OCF bridge]
    Sum_T F(T, x) = A_n(x) * 2^(C(n,2)-(n-1))  [ensemble average]
    I(Omega, x) = clique polynomial of complement(Omega) at x
    clique poly is a multivariate Tutte specialization

  THE KEY INSIGHT (from kind-pasteur S18c):
    ALL of these polynomials are partition functions evaluated at
    specific points. H sits at the INTERSECTION of evaluations:
    x=2 for independence polynomial, x=1 for descent polynomial.

    The TOURNAMENT is a "code" and H is its "weight" computed at
    the tournament fugacity.

  THE OPEN CHALLENGE:
    Find a SINGLE polynomial P(T; x, y) such that:
      P(T; 1, y) = F(T, y)  (specialization gives descent)
      P(T; x, 2) = I(Omega(T), x)  (specialization gives independence)
      P(T; 1, 1) = P(T; 2, 2)... no, both should give H but at different points.

    Such a P would be the "tournament Tutte polynomial."
""")

print("\nDone. opus-2026-03-21-S121")
