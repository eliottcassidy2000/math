#!/usr/bin/env python3
"""petersen_tricks_s116k.py — Tricks from the Cayley-Paley-Petersen connection.

The Petersen graph sits at the intersection of {2,3,5,7,10}.
Its fractional chromatic number = Q(3/7) = 5/2.
Its vertices = C(5,2) = T_4 = 10 = the reset cost.
Its girth = 5 = the pentagon. Its regularity = 3 = curvature quantum.

What can we DO with this?
"""
from math import log, sqrt, comb
from fractions import Fraction
from itertools import combinations

print()
print("  PETERSEN GRAPH TRICKS")
print()
print("="*70)
print()

# ============================================================
print("  TRICK 1: THE FRACTIONAL CHROMATIC NUMBER VIA CAYLEY")
print("  " + "-"*40)
print()
print("  For the Kneser graph K(n,k): vertices = k-subsets of [n],")
print("  edges between disjoint subsets.")
print("  Fractional chromatic number chi_f(K(n,k)) = n/k.")
print()
print("  In Cayley: chi_f = n/k. And Q((n-2k)/(n)) = n/(2k)... hmm.")
print("  Actually: chi_f(K(n,k)) = n/k.")
print("  For Petersen K(5,2): chi_f = 5/2.")
print()
print("  Now: 5/2 = Q(3/7). The Cayley transform of the curvature/threshold ratio.")
print("  Is there a GENERAL pattern?")
print()
print("  For K(n,k): chi_f = n/k = Q((n-2k)/(n+2k))? Let's check:")
for n in range(3, 10):
    for k in range(1, n//2 + 1):
        chi_f = Fraction(n, k)
        # What x gives Q(x) = n/k?
        # Q(x) = (1+x)/(1-x) = n/k
        # x = (n/k - 1)/(n/k + 1) = (n-k)/(n+k)
        x = Fraction(n-k, n+k)
        q_check = (1+x)/(1-x)
        print(f"  K({n},{k}): chi_f = {chi_f}, Q^{{-1}}({chi_f}) = {x}, Q({x}) = {q_check}, match: {q_check == chi_f}")

print()
print("  YES. The Cayley address of chi_f(K(n,k)) is always (n-k)/(n+k).")
print("  For Petersen K(5,2): address = 3/7. Curvature over threshold.")
print("  For K(7,3): address = 4/10 = 2/5. Q(2/5) = 7/3.")
print("  For K(7,2): address = 5/9. Q(5/9) = 7/2.")
print()
print("  The Cayley address of the fractional chromatic number of K(n,k)")
print("  = (n-k)/(n+k) = the Cayley address of (n-k) relative to (n+k).")
print("  = the 'imbalance ratio' of the Kneser parameters.")
print()

# ============================================================
print()
print("  TRICK 2: INDEPENDENT SETS OF THE PETERSEN GRAPH")
print("  " + "-"*40)
print()
print("  The independence number alpha(Petersen) = 4.")
print("  The independence polynomial I(Petersen, x):")
print()

# Build the Petersen graph
V = list(combinations(range(5), 2))  # 10 vertices = 2-subsets of {0,1,2,3,4}
adj = [[0]*10 for _ in range(10)]
for i in range(10):
    for j in range(i+1, 10):
        if len(set(V[i]) & set(V[j])) == 0:  # disjoint = adjacent in Kneser
            adj[i][j] = adj[j][i] = 1

# Count independent sets by size
alpha = [0] * 11
for size in range(11):
    for subset in combinations(range(10), size):
        independent = True
        for a in range(len(subset)):
            for b in range(a+1, len(subset)):
                if adj[subset[a]][subset[b]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            alpha[size] += 1

print(f"  I(Petersen, x) = ", end="")
terms = []
for k in range(11):
    if alpha[k] > 0:
        terms.append(f"{alpha[k]}*x^{k}")
print(" + ".join(terms))
print()

# Evaluate at key points
print(f"  I(Petersen, 1) = {sum(alpha)} (total independent sets)")
print(f"  I(Petersen, 2) = {sum(alpha[k]*2**k for k in range(11))}")
print(f"  I(Petersen, -1) = {sum(alpha[k]*(-1)**k for k in range(11))}")
print()

# The OCF connection: if this were a tournament conflict graph,
# H = I(Omega, 2). What H would the Petersen graph give?
H_petersen = sum(alpha[k] * 2**k for k in range(11))
print(f"  If Petersen were a conflict graph Omega(T):")
print(f"  H = I(Petersen, 2) = {H_petersen}")
print(f"  Factor: {H_petersen}", end="")
temp = H_petersen
factors = []
d = 2
while d*d <= temp:
    while temp % d == 0:
        factors.append(d)
        temp //= d
    d += 1
if temp > 1: factors.append(temp)
print(f" = {'*'.join(str(f) for f in factors)}")
print()

# ============================================================
print()
print("  TRICK 3: PETERSEN AS A TOURNAMENT CONFLICT GRAPH?")
print("  " + "-"*40)
print()
print("  Can the Petersen graph arise as Omega(T) for some tournament T?")
print("  Omega(T) has vertices = directed odd cycles of T,")
print("  edges = pairs sharing a vertex.")
print()
print("  The Petersen graph is vertex-transitive, 3-regular, girth 5.")
print("  Tournament conflict graphs at n=5: all cycles are 3-cycles or 5-cycles.")
print("  At n=5: max cycles = 4 (three-cycles) + 2 (five-cycles) = 6.")
print("  The Petersen graph has 10 vertices. That's too many for n=5.")
print()
print("  At n=7: the max number of 3-cycles = C(7,3) - ... could be up to 14.")
print("  Petersen has 10 vertices. Could 10 of the 14 three-cycles form")
print("  a Petersen-like conflict graph? Unlikely (Petersen is 3-regular,")
print("  but conflict graphs of 3-cycles at n=7 have varying degree).")
print()
print("  The Petersen graph probably does NOT arise as Omega(T).")
print("  But I(Petersen, 2) = {H_petersen} is the H-value it WOULD give.")
print()

# ============================================================
print()
print("  TRICK 4: THE CAYLEY-PETERSEN IDENTITY")
print("  " + "-"*40)
print()
print("  Q(3/7) = 5/2 = chi_f(Petersen).")
print("  3 + 7 = 10 = |V(Petersen)|.")
print("  3 * 7 = 21 = FORBIDDEN H-value.")
print("  C(3+7, 2) = C(10, 2) = 45 = |E(K_{10})| = edges of the complete graph on Petersen's vertices.")
print()
print("  These are not independent facts. They all flow from:")
print("  Petersen = K(5, 2) where 5 = Q(2/3) and 2 is the Cayley derivative.")
print()
print("  The number 5 enters via Q: it is the image of the simplest")
print("  non-trivial address (2/3). Then K(5,2) builds the Petersen graph")
print("  by selecting 2-subsets of [5]. The parameter 2 is the Cayley")
print("  derivative Q'(0) = 2. So Petersen = K(Q(2/3), Q'(0)/Q'(0)).")
print("  = K(Q(2/3), 1)... no, K(5,2) not K(5,1).")
print()
print("  More directly: the Petersen graph's structure constants ARE")
print("  the Cayley framework's fundamental numbers:")
print("  |V| = 10 = T_4 = the reset cost.")
print("  |E| = 15 = T_5 = the next triangular.")
print("  degree = 3 = the curvature quantum.")
print("  girth = 5 = the pentagon.")
print("  chi_f = 5/2 = Q(3/7).")
print("  chi = 3 (chromatic number) = the curvature quantum.")
print("  alpha = 4 = the Cayley period.")
print("  |Aut| = 120 = 5! = the factorial of the pentagon.")
print()
print("  EVERY invariant of the Petersen graph is one of our numbers.")
print()

# ============================================================
print()
print("  TRICK 5: USING CHORD TO PREDICT KNESER INDEPENDENCE POLYNOMIALS")
print("  " + "-"*40)
print()
print("  For K(n,k), the independence number alpha = C(n-1, k-1).")
print("  (An independent set in K(n,k) = k-subsets that pairwise intersect")
print("  = k-subsets all containing a common element = stars.)")
print()
for n in range(4, 8):
    for k in range(1, n//2 + 1):
        alpha_val = comb(n-1, k-1)
        chi_f = Fraction(n, k)
        cayley_addr = Fraction(n-k, n+k)
        print(f"  K({n},{k}): alpha = C({n-1},{k-1}) = {alpha_val}, "
              f"chi_f = {chi_f} = Q({cayley_addr})")

print()
print("  The independence number alpha(K(n,k)) = C(n-1, k-1).")
print("  The fractional chromatic number chi_f = n/k = Q((n-k)/(n+k)).")
print("  Their product: alpha * chi_f = C(n-1,k-1) * n/k.")
print()
for n in range(4, 8):
    for k in range(1, n//2 + 1):
        alpha_val = comb(n-1, k-1)
        prod = Fraction(alpha_val * n, k)
        total_v = comb(n, k)
        print(f"  K({n},{k}): alpha * chi_f = {prod} = {float(prod):.1f}. "
              f"|V| = {total_v}. Ratio |V|/(alpha*chi_f) = {Fraction(total_v, 1)/prod}")

print()
print("  alpha * chi_f = |V| always! (This is a known identity for vertex-transitive graphs.)")
print("  |V| = alpha * chi_f.")
print("  For Petersen: 10 = 4 * 5/2. CHECK.")
print()
print("  In Cayley terms: |V(K(n,k))| = alpha(K(n,k)) * Q((n-k)/(n+k)).")
print("  The vertex count = independence number * Cayley transform of the imbalance.")
print()

# ============================================================
print()
print("  TRICK 6: THE SPECTRUM OF THE PETERSEN GRAPH")
print("  " + "-"*40)
print()
print("  The adjacency matrix eigenvalues of the Petersen graph:")
print("  lambda = 3 (multiplicity 1), 1 (multiplicity 5), -2 (multiplicity 4).")
print()
print("  In our framework:")
print("  3 = the curvature quantum. The dominant eigenvalue.")
print("  1 = the identity. Multiplicity 5 = the pentagon.")
print("  -2 = negative doubling. Multiplicity 4 = the Cayley period.")
print()
print("  The eigenvalues ARE our fundamental numbers: 3, 1, -2.")
print("  Their Cayley addresses:")
print(f"  Q^{{-1}}(3) = 1/2. The midpoint.")
print(f"  Q^{{-1}}(1) = 0. The center.")
print(f"  Q^{{-1}}(-2): Q(x) = -2 => (1+x)/(1-x) = -2 => 1+x = -2+2x => x = 3.")
print(f"  Address of -2 would be 3, which is > 1. On the upper sheet.")
print(f"  arctanh(3) = arctanh(1/3) + i*pi/2 = ln(2)/2 + i*pi/2.")
print()
print("  The eigenvalue -2 lives on the UPPER SHEET of the Cayley helix!")
print("  Its rapidity: (ln(2)/2, pi/2). The octave, rotated by a quarter turn.")
print("  -2 is the SHADOW of 2 on the imaginary sheet.")
print("  The negative eigenvalue carries the imaginary phase of the helix.")
print()
print("  The Petersen spectrum in rapidity:")
print(f"  lambda=3:  rapidity ln(3)/2 = {log(3)/2:.4f} (real, ground floor)")
print(f"  lambda=1:  rapidity 0 (the center)")
print(f"  lambda=-2: rapidity ln(2)/2 + i*pi/2 (complex, upper floor)")
print()
print("  The spectrum lives on TWO FLOORS of the Cayley helix.")
print("  The positive eigenvalues on the ground floor.")
print("  The negative eigenvalue on the upper floor.")
print("  The Petersen graph SPANS the helix.")
print()

# ============================================================
print()
print("  TRICK 7: PETERSEN AND THE FORBIDDEN")
print("  " + "-"*40)
print()
print("  The Petersen graph is the SMALLEST graph that is:")
print("  - Not 3-edge-colorable (obstruction to 3-coloring)")
print("  - A snark (bridgeless, 3-regular, non-3-edge-colorable)")
print()
print("  3-edge-coloring fails because 3 = curvature quantum.")
print("  The curvature quantum OBSTRUCTS the Petersen graph")
print("  just as it OBSTRUCTS H=7 in tournaments.")
print()
print("  In both cases: 3 is the source of the impossibility.")
print("  3 pairwise-intersecting cycles cascade at H=7.")
print("  3 colors can't cover the Petersen edges.")
print()
print("  The NUMBER 3 is the fundamental obstruction.")
print("  It is too small to tile certain structures (Petersen, H=7)")
print("  but too large to be trivial (unlike 1 or 2).")
print()
print("  3 is the curvature quantum because it is the SMALLEST")
print("  number that creates non-trivial topology:")
print("  - The smallest cycle has 3 vertices.")
print("  - The smallest non-planar complete graph K_{3,3} uses 3+3.")
print("  - The smallest non-trivial group Z_3 has 3 elements.")
print("  - The smallest odd prime is 3.")
print()
print("  Everything that goes wrong in combinatorics goes wrong at 3.")
print("  The Petersen graph is WHERE it goes wrong for edge-coloring.")
print("  H=7 is WHERE it goes wrong for Hamiltonian path counts.")
print("  Both are obstructed by the same 3.")
print()

# ============================================================
print()
print("  SUMMARY OF TRICKS")
print("  " + "-"*40)
print()
print("  1. chi_f(K(n,k)) = Q((n-k)/(n+k)). Fractional chromatic via Cayley.")
print("  2. I(Petersen, 2) = {0}. Independence polynomial at tournament fugacity.".format(H_petersen))
print("  3. Petersen probably NOT a tournament conflict graph.")
print("  4. Every Petersen invariant is one of our fundamental numbers.")
print("  5. |V| = alpha * chi_f = alpha * Q(imbalance). Vertex-transitive identity.")
print("  6. Spectrum {3, 1, -2} spans two floors of the Cayley helix.")
print("  7. Petersen and H=7 share the same obstruction: the number 3.")
