#!/usr/bin/env python3
"""
cumulant_ihara.py - Explore connection between cumulant hierarchy and Ihara zeta.

The Ihara zeta function of a digraph D is:
  zeta_D(u) = prod_{prime cycles C} (1 - u^{|C|})^{-1}

For a tournament T:
  log zeta_T(u) = sum_{k>=1} N_k * u^k / k
where N_k = #{directed closed walks of length k} = tr(A^k).

The cumulant hierarchy involves:
  t_3 = c_3/2 = directed 3-cycle count / 2 (since each undirected 3-cycle gives 2 directed)
  t_5, t_7 etc.

Actually, t_{2k+1} in our convention counts directed (2k+1)-cycles (one per direction?).
Let me clarify our convention from THM-089 and the computation scripts.

Our convention from count_3cycles: t3 counts UNDIRECTED 3-cycles
(each vertex triple with a cyclic orientation counted once).
So t3 = c_3 / 2 where c_3 = #directed 3-cycles.

For t5: count_5cycles in our scripts counts directed 5-cycles / 5 (dividing by start vertex).
A directed 5-cycle has 5 starting points but 2 directions.
An undirected 5-cycle gives 2 directed 5-cycles, each counted 5 times by the permutation loop.
So our t5 = #directed 5-cycles / 5 = 2 * #undirected 5-cycles * 5 / 5 = ... hmm.

Let me just check from the scripts.

Author: opus-2026-03-07-S46d
"""
from itertools import permutations, combinations
from fractions import Fraction

# At n=5: the doubly-regular tournament (Paley T5) has all vertices equal.
# T5 has score (2,2,2,2,2) and is unique up to isomorphism.
# It has t3 = 5 (there are exactly C(5,3)/? 3-cycles).
# Wait, C(5,3) = 10 triples. In a regular tournament on 5 vertices,
# each triple is either transitive (no cycle) or cyclic (one cycle).
# Score (2,2,2,2,2) => each vertex has outdeg 2.
# By a well-known identity: t3 = C(n,3) - sum C(s_i, 2)
# = 10 - 5*C(2,2) = 10 - 5 = 5. So t3=5. ✓

# For t5: the Paley T5 has 2*5!/5 = 48 directed 5-cycles? No...
# Actually in a tournament on 5 vertices, a directed 5-cycle uses all 5 vertices.
# The number of directed Hamiltonian cycles = (n-1)!/2 * (fraction that are cycles).
# For T5 with t3=5, H(T5)=15 by OCF: 1 + 2*5 + 2*t5 = 15 => t5=2.

# So t5=2 for the regular T5. Let me verify with the counting function.
n = 5
# Build Paley T5: 0->1, 1->2, 2->3, 3->4, 4->0 (pentagons are the QR edges)
# Actually Paley T5 has connection set {1,4} in Z/5 (QR mod 5 = {1,4}).
# Wait, 5 ≡ 1 mod 4 so Paley doesn't give a tournament. Let me just use the cyclic tournament.
# C5 tournament: i->j iff j-i in {1,2} mod 5.
adj = [[0]*n for _ in range(n)]
for i in range(n):
    for d in [1, 2]:
        j = (i + d) % n
        adj[i][j] = 1

t3 = sum(1 for i,j,k in combinations(range(n), 3)
         if (adj[i][j] and adj[j][k] and adj[k][i]) or
            (adj[i][k] and adj[k][j] and adj[j][i]))

# t5 from our counting function
t5_count = 0
for combo in combinations(range(n), 5):
    for perm in permutations(combo):
        if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
            t5_count += 1
t5 = t5_count // 5

print(f"C5 tournament at n=5: t3={t3}, t5={t5}")
print(f"OCF: H = 1 + 2*{t3} + 2*{t5} = {1 + 2*t3 + 2*t5}")

# t5_count = total directed 5-cycles * starting points / 5
# Actually t5_count / 5 gives: each directed cycle counted once.
# But a 5-vertex set can have 2 directed 5-cycles (CW and CCW).
# So t5 = number of directed 5-cycles / 1 = 2 * (undirected 5-cycles)?
# No: t5_count // 5 = number of directed 5-cycles counting each direction separately.

# Let me check: at n=5 with the C5 tournament (all 5 vertices):
# The C5 tournament has the Hamiltonian cycle 0->1->2->3->4->0 (if adj allows).
# adj[0][1]=1, adj[1][2]=1, adj[2][3]=1, adj[3][4]=1, adj[4][0]=1. Yes!
# And the reverse: 0->4->3->2->1->0. adj[0][4]? d=4, (0+1)%5=1, (0+2)%5=2. Not 4.
# adj[0][4]: 4-0=4 mod 5 = 4. Is 4 in {1,2}? No! So adj[0][4]=0.
# The reverse cycle doesn't work. So there's only 1 directed 5-cycle direction?
# Let me check all 24 = 4! orderings as directed cycles starting from 0.

print("\nDirected 5-cycles in C5 tournament:")
for perm in permutations(range(5)):
    if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
        print(f"  {'->'.join(str(p) for p in perm)}->{perm[0]}")

# Actually permutations gives all 120 orderings. Each directed cycle is counted 5 times.
# So t5 = t5_count // 5 counts directed cycles (both directions if both exist).

# Now: the connection to Ihara zeta.
# The Ihara zeta function counts PRIME cycles (non-backtracking, not a power of shorter).
# For tournaments, all cycles are non-backtracking (since each edge has one direction).
# A directed k-cycle is prime iff it's not a power of a shorter cycle.
# For odd k (prime cycle length), all directed k-cycles are prime.

# So: log zeta_T(u) = sum_{k>=1} N_k * u^k / k
# where N_k = tr(A^k) = sum of directed closed walks of length k.
# N_3 = 3*c3 = 6*t3 (each 3-cycle gives 3 closed walks * 2 directions? No...)
# Actually N_3 = tr(A^3) = sum_{i,j,k} A[i][j]*A[j][k]*A[k][i]
# = 3 * (#directed 3-cycles). Wait, each directed 3-cycle i->j->k->i
# gives closed walks starting at i, j, k. So 3 walks per directed cycle.
# And there are 2*t3 directed 3-cycles.
# So N_3 = 3 * 2*t3 = 6*t3.

# Actually: N_3 = tr(A^3) = sum_i (A^3)[i][i] = sum of all directed triangles through each vertex.
# Each directed triangle i->j->k->i contributes 1 to (A^3)[i][i], plus j and k start versions.
# So each directed 3-cycle contributes 3 to N_3.
# Number of directed 3-cycles = 2*t3 (each undirected gives CW and CCW).
# N_3 = 3*2*t3 = 6*t3.

# This is the classical identity: t3 = tr(A^3)/6. ✓

# For 5-cycles: N_5 = tr(A^5) = ... complex (includes non-simple walks).
# But directed simple 5-cycles contribute 5*2*t5_undirected to N_5.
# Plus non-simple walks of length 5.

# INSIGHT: The cumulant hierarchy coefficients 2/C(n,2k) have a nice interpretation.
# kappa_{2k} probes the (2k+1)-cycle structure.
# The NEW contribution from t_{2k+1} is: 2/C(n,2k) * t_{2k+1}.
# This is 2*t_{2k+1} / C(n,2k).

# Now t_{2k+1} * 2 = contribution of (2k+1)-cycles to OCF at the FIRST LEVEL.
# In OCF: H(T) = 1 + 2*alpha_1 + 4*alpha_2 + ...
# where alpha_1 = t3 + t5 + t7 + ... (all directed odd cycles).
# The "2*t_{2k+1}" is the alpha_1-level contribution from (2k+1)-cycles.

# So kappa_{2k} has NEW cycle contribution = (alpha_1-level OCF contribution of new cycles) / C(n,2k).

# And C(n,2k) = number of (2k)-element subsets = dimension of the moment space at level 2k.

# This is the cumulant-OCF bridge: cumulants "unpack" the OCF contributions
# by filtering through moment level, with binomial coefficient normalization.

print("\n" + "=" * 60)
print("CUMULANT-OCF BRIDGE")
print("=" * 60)
print()
print("OCF: H(T) = 1 + 2*sum_k t_{2k+1} + 4*alpha_2 + ...")
print()
print("Cumulant hierarchy:")
print("  kappa_{2k} = Bernoulli_constant + 2*t_{2k+1}/C(n,2k) + lower-order terms")
print()
print("Each even cumulant captures the OCF contribution of (2k+1)-cycles,")
print("normalized by C(n,2k) = dim of the 2k-th moment interaction space.")
print()
print("This means: the full set of even cumulants {kappa_2, kappa_4, kappa_6, ...}")
print("ENCODES the OCF data {t3, t5, t7, ...} in a graded fashion.")
print()
print("In principle: knowing all kappa_{2k} determines H(T) via OCF.")
print("But the nonlinear cross terms (t3^2 in kappa_4, t3*t5 in kappa_6)")
print("make the recovery non-trivial.")
