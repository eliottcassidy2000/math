#!/usr/bin/env python3
"""
coefficient_proof.py - Prove the 2/C(n,2k) coefficient pattern.

The claim: in kappa_{2k}(T), the coefficient of t_{2k+1}
(count of directed (2k+1)-cycles) is 2/C(n,2k).

Proof sketch for kappa_2 (k=1):
  kappa_2 = Var[fwd] = sum_i Var(X_i) + 2*sum_{i<j} Cov(X_i, X_j)
  = (n-1)/4 + 2*(n-2)*Cov(X_i, X_{i+1})  [only adjacent pairs correlated]

  Cov(X_i, X_{i+1}) = E[X_i X_{i+1}] - 1/4
  E[X_i X_{i+1}] = #{2-paths a->b->c} / (n(n-1)(n-2))
  = [C(n,3) + 2*t3] / (n(n-1)(n-2))

  So: Cov = 2*t3/(n(n-1)(n-2))
  Var = (n-1)/4 + 2*(n-2)*2*t3/(n(n-1)(n-2))
      = (n-1)/4 + 4*t3/(n(n-1))

  Coefficient of t3 = 4/(n(n-1)) = 2/C(n,2). Done for k=1.

Now for kappa_4 (k=2):
  The t5 dependence enters through 4th-order correlations.
  E[X_i X_{i+1} X_{i+2} X_{i+3}] for consecutive quadruple.
  This involves the 5-vertex substructure {P[i],...,P[i+4]}.

  The number of directed 4-paths a->b->c->d->e is:
  #fwd5path = C(n,5) * #{4-paths in random 5-tournament}
  + correction from cycle structure

  Let me compute: among all 5-vertex subtournaments, how many
  directed 4-paths (a0->a1->a2->a3->a4) does each have?

For a 5-vertex tournament, the number of directed 4-paths
through all 5 vertices (hamiltonian paths) equals H(T[5-subset]).

So sum over 5-vertex subsets of H(T[S]) counts 5-step directed paths.

Actually, #fwd5path = #{ordered 5-tuples (a,b,c,d,e) distinct : a->b->c->d->e}
        = sum_{S in C(V,5)} #{hamiltonian paths of T[S]}
        = sum_{S} H(T[S])

Now H(T[S]) = I(Omega(T[S]), 2) by OCF.
For a 5-vertex tournament: H(T[S]) = 1 + 2*t3(S) + 2*t5(S)
(at n=5: alpha_2 = 0 always, so just t3 and t5 contributions)

So: sum_S H(T[S]) = C(n,5) + 2*sum_S t3(S) + 2*sum_S t5(S)

sum_S t3(S): each 3-cycle is in C(n-3,2) five-vertex subsets.
So sum_S t3(S) = t3 * C(n-3, 2).

sum_S t5(S): each directed 5-cycle is in exactly 1 five-vertex subset (uses all 5).
Wait, at n>5, a 5-cycle uses exactly 5 vertices, so it appears in C(n-5, 0)=1 subset.
So sum_S t5(S) = t5 * 1 = t5.

Therefore: #fwd5path = C(n,5) + 2*C(n-3,2)*t3 + 2*t5

The E[fwd^4] formula involves clusters of at most 4 adjacent X_i's.
The cluster of 4 adjacent positions involves 5 vertices,
and contributes to E[X_i X_{i+1} X_{i+2} X_{i+3}].

E[X_i X_{i+1} X_{i+2} X_{i+3}] = #fwd5path / (n(n-1)(n-2)(n-3)(n-4))
= [C(n,5) + 2*C(n-3,2)*t3 + 2*t5] / P(n,5)

where P(n,5) = n!/(n-5)! = n(n-1)(n-2)(n-3)(n-4).

The t5 coefficient in this expectation is 2/P(n,5).

In E[fwd^4], the contribution from 4-clusters is:
2*(n-4) * [E[X_i X_{i+1} X_{i+2} X_{i+3}] - lower-order terms]

Wait, this is getting complicated. But the key insight:
the t5 dependence in E[fwd^4] comes ONLY from 4-clusters,
and each 4-cluster samples a 5-vertex sub-tournament.

The coefficient 2/C(n,4) = 2*4!/(n(n-1)(n-2)(n-3)) = 48/P(n,4).
And 2/P(n,5) = 2/(n(n-1)(n-2)(n-3)(n-4)).

The ratio is (n-4):1, which is the number of 4-clusters in n-1 positions.
Actually, there are (n-4) consecutive quadruples (i, i+1, i+2, i+3) for i=0..n-5.

Hmm, let me compute this more carefully.

Author: opus-2026-03-07-S46d
"""
from fractions import Fraction
from math import comb, factorial, perm

# Check the claim about fwd5path
# #fwd5path = sum over C(V,5) of H(T[S])
# For the transitive tournament: H(T[S]) = 1 for all S (transitive restriction is transitive)
# Wait no, H(transitive on 5) = 1? No! H(transitive on n) = 1 only when the permutation is identity.
# Actually no: H(transitive T5) = # Hamiltonian paths = sum A(5,k) [just the total perms that are HPs]
# For transitive tournament, EVERY permutation is a HP if all edges go forward, but
# a permutation sigma is an HP iff sigma[0]->sigma[1]->...->sigma[n-1] in T.
# For transitive T, adj[i][j]=1 iff i<j. So sigma is HP iff sigma is the identity.
# Wait that's wrong too. Actually H(T) counts directed Hamiltonian PATHS, not cycles.
# An HP is a permutation sigma such that T[sigma[i]][sigma[i+1]] = 1 for all i.
# For transitive T5 with adj[i][j]=1 iff i<j, the only HP is the identity: 0->1->2->3->4.
# H(transitive T5) = 1. But we said F(trans,x) = Eulerian poly, so sum F(trans,k) = n! for all m.
# Hmm, H(T) = F(T,1)? No, F(T,1) = sum F_k = n! always.
# F_k = #{sigma : fwd(sigma) = k} and sum F_k = n!.
# For the transitive tournament, fwd(sigma) = des(sigma), so F_k = A(n,k).
# A(5,0)=1, A(5,1)=26, A(5,2)=66, A(5,3)=26, A(5,4)=1.
# H(T) = F_{n-1} = A(n,n-1) = 1 for the transitive tournament. OK so H(trans)=1.

# But for a non-transitive 5-vertex tournament, H can be > 1.
# And H(T) = F_{n-1}(T) = #{sigma : all consecutive edges forward} = #{HPs}.

# OK so #fwd5path = #{ordered 5-tuples a->b->c->d->e all forward}
# This is the same as the number of HPs in 5-vertex sub-tournaments, summed:
# #fwd5path = sum_{S in C(V,5)} H(T[S])

# But wait, "a->b->c->d->e all forward" means a,b,c,d,e are 5 DISTINCT vertices
# with T[a][b]=T[b][c]=T[c][d]=T[d][e]=1. This is a directed path of length 4
# through exactly 5 vertices. And for each 5-element subset S, H(T[S]) counts
# exactly the number of such paths within S.
# YES, so #fwd5path = sum_S H(T[S]).

# Now H(T[S]) at 5 vertices = I(Omega(T[S]), 2) = 1 + 2*t3(S) + 2*t5(S)
# (since alpha_2 = 0 at n=5)

# So #fwd5path = sum_S [1 + 2*t3(S) + 2*t5(S)]
#             = C(n,5) + 2*[sum_S t3(S)] + 2*[sum_S t5(S)]

# Each directed 3-cycle on vertices {a,b,c} is contained in C(n-3, 2) 5-vertex subsets.
# So sum_S t3(S) = t3 * C(n-3, 2).

# Each directed 5-cycle on vertices {a,b,c,d,e} is in exactly 1 five-vertex subset.
# Actually t5(S) counts directed 5-cycles within S. If S has 5 vertices, there can be
# 0 or more directed 5-cycles. Each directed 5-cycle uses exactly 5 vertices.
# So a given 5-cycle appears in exactly C(n-5, 0) = 1 subset (its own vertex set).
# Therefore sum_S t5(S) = t5.

# #fwd5path = C(n,5) + 2*C(n-3,2)*t3 + 2*t5

print("FORWARD 5-PATH COUNT")
print("=" * 60)
for n in [5, 6, 7]:
    formula = f"C({n},5) + 2*C({n-3},2)*t3 + 2*t5"
    c5 = comb(n, 5)
    cn32 = comb(n-3, 2)
    print(f"n={n}: #fwd5path = {c5} + {2*cn32}*t3 + 2*t5")

# Check against known: at n=5, #fwd4path = C(5,4) + 2*3*t3 = 5 + 6*t3
# Wait, that's 4-paths (4 edges, 5 vertices). Let me clarify.
# fwd-r-path uses r consecutive edges, hence r+1 vertices.
# fwd-4-path: 4 edges, 5 vertices.
# From THM-092: #fwd4path = C(n,4) + 2(n-2)*t3.
# At n=5: 5 + 6*t3.
# But above I computed #fwd5path which is a 5-VERTEX path = 4-edge path.
# So they should be the same! Let me check:
# C(n,5) + 2*C(n-3,2)*t3 + 2*t5 at n=5:
# C(5,5)=1 + 2*C(2,2)*t3 + 2*t5 = 1 + 2*t3 + 2*t5.
# But THM-092 says #fwd4path(n=5) = C(5,4) + 2*3*t3 = 5 + 6*t3.
# These don't match! Let me re-examine.

# Ah, the confusion: "fwd4path" in THM-092 means a 4-step forward path,
# not a path using 4 vertices. Let me re-read THM-092.

# From THM-092:
# #fwd3path = C(n,3) + 2*t3
# #fwd4path = C(n,4) + 2(n-2)*t3
# These count ordered directed paths through r+1 vertices of each (r+1)-subset.

# Wait, "fwd3path" counts paths using 3 consecutive edges? Or 3 vertices?
# From the theorem: #fwd2path = C(n,2), #fwd3path = C(n,3) + 2*t3.
# If fwd2path = C(n,2) = #{ordered pairs (a,b) : a->b}, that's just the number of edges.
# For tournament: = C(n,2). ✓ (every pair has exactly one directed edge)

# fwd3path = #{ordered triples (a,b,c) distinct : a->b->c}
# = sum over C(V,3) of #{directed 2-paths in T[S]}
# In a 3-vertex tournament T[S]: either 0 or 2 directed 2-paths.
# #{2-paths} = 2 if T[S] is transitive, 0 if T[S] is a 3-cycle.
# Wait no: in a transitive 3-tournament (say 0->1->2, 0->2), the 2-paths are
# 0->1->2 only. And the paths going "forward" are 0->1->2, 2->0->1 is NOT forward (2->0 not in T).
# Actually wait, #{ordered triples (a,b,c) : a->b->c} counts ALL directed 2-paths.
# In transitive {0,1,2} with 0->1, 1->2, 0->2:
#   2-paths: 0->1->2. That's it. Oh wait, we need a->b->c with a,b,c all distinct.
#   0->1->2: yes. 0->2->?: 2 beats nobody (if 0->1, 0->2, 1->2, this is T[S]=transitive)
#   Actually 2 has out-degree 0, so no 2-path starting at 2.
#   1->2->?: 2 has out-degree 0.
#   So only 1 directed 2-path. Hmm.

# Actually for a 3-vertex transitive tournament 0->1, 0->2, 1->2:
# Out-degrees: out(0)=2, out(1)=1, out(2)=0
# 2-paths a->b->c: need b to have out-degree >= 1 and in-degree >= 1 (from a).
# a=0: 0->1->2 (1 path), 0->2->? (2 has out=0, 0 paths). Total from a=0: 1.
# a=1: 1->2->? (2 has out=0). Total: 0.
# a=2: 2->? (2 has out=0). Total: 0.
# Grand total: 1 directed 2-path.

# For a 3-cycle 0->1->2->0 (i.e., 0->1, 1->2, 2->0):
# 2-paths: 0->1->2, 1->2->0, 2->0->1. Total: 3 directed 2-paths.

# So #2paths = sum_S #{2-paths in T[S]}:
# Transitive triple: 1 path. Cyclic triple: 3 paths.
# #2paths = (C(n,3) - t3)*1 + t3*3 = C(n,3) + 2*t3. ✓!

# So "fwd(r+1)path" counts #{ordered (r+1)-tuples distinct forming a directed r-path}

# For r=4 (4-edge path, 5 vertices):
# #4paths = sum_{S in C(V,5)} #{directed 4-paths in T[S]}
# where #{directed 4-paths in T[S]} = H(T[S]) = 1 + 2*t3(S) + 2*t5(S)

# Hmm but H(T[S]) = #{Hamiltonian paths of T[S]}
# A Hamiltonian path in T[S] is an ordering v1->v2->v3->v4->v5 of all 5 vertices
# such that all edges go forward.
# So H(T[S]) = #{4-edge directed paths using all 5 vertices}.
# But directed 4-paths can also use FEWER than 5 vertices? No, a 4-path visits 5 distinct vertices.
# So #4paths = sum_S H(T[S]). ✓

# Now at n=5, sum_S H(T[S]) has only S=V, so #4paths = H(T).
# H(T) at n=5 = 1 + 2*t3 + 2*t5 (from OCF).
# But THM-092 says #fwd4path = C(5,4) + 2*3*t3 = 5 + 6*t3.
# These don't match! At n=5 with t3=0, t5=0: H=1 but C(5,4)=5.

# AH, the issue: THM-092's "fwd4path" must be counting SUBSETS of 4 vertices,
# not subsets of 5. Let me re-read.

# From THM-092:
# #fwd2path = C(n,2) -- this is 2 vertices, 1 edge
# #fwd3path = C(n,3) + 2*t3 -- this is 3 vertices, 2 edges
# #fwd4path = C(n,4) + 2(n-2)*t3 -- this is 4 vertices, 3 edges

# So fwd-r-path uses r vertices and r-1 edges!

# fwd4path = #{ordered 4-tuples (a,b,c,d) distinct : a->b->c->d}
# = sum_{S in C(V,4)} #{directed 3-paths in T[S]}
# Wait, that's #{3-edge directed paths through 4 distinct vertices}.
# In a 4-vertex tournament T[S], #{3-edge directed HP} = H(T[S]).
# H(T[S]) at n=4: 1 + 2*t3(S). (alpha_2=0 for n=4, and t5 doesn't exist)

# OK so:
# #fwd(r)path (r vertices, r-1 edges) = sum_{S in C(V,r)} H(T[S])

# For r=4: sum_S H(T[S]) = sum_S [1 + 2*t3(S)]
# Each 3-cycle {a,b,c} is in C(n-3, 1) = n-3 four-vertex subsets.
# sum_S t3(S) = t3 * (n-3).
# #fwd4path = C(n,4) + 2*(n-3)*t3.

# At n=5: C(5,4) + 2*2*t3 = 5 + 4*t3.
# But THM-092 says C(5,4) + 2*3*t3 = 5 + 6*t3.
# Hmm, (n-3) = 2 at n=5, but THM-092 has coefficient 2*(n-2) = 6.

# Let me recheck. At n=5, there are C(5,4)=5 four-vertex subsets.
# Each 3-cycle on 3 vertices is contained in (n-3) = 2 four-vertex subsets?
# Wait: to form a 4-subset containing a given 3-element subset,
# we add 1 vertex from the remaining n-3 = 2 vertices.
# So yes, n-3 = 2 at n=5.
# #fwd4path = 5 + 2*2*t3 = 5 + 4*t3.

# But THM-092 from the previous session claims #fwd4path = C(n,4) + 2(n-2)*t3.
# At n=5 that's 5 + 6*t3. Let me verify numerically.

from itertools import permutations as perms, combinations as combs

def t_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

n = 5
print(f"\nNUMERICAL CHECK AT n={n}")
for bits in [0]:  # transitive (all 0-bits = all j beats i)
    adj = t_from_bits(n, bits)
    # Also construct actual transitive
    adj = [[1 if i < j else 0 for j in range(n)] for i in range(n)]

    # Count directed 3-edge paths (4 distinct vertices, a->b->c->d)
    count_3path = 0
    for a,b,c,d in perms(range(n), 4):
        if adj[a][b] and adj[b][c] and adj[c][d]:
            count_3path += 1

    # Count directed 4-edge paths (5 distinct vertices)
    count_4path = 0
    for p in perms(range(n)):
        if all(adj[p[i]][p[i+1]] for i in range(4)):
            count_4path += 1

    t3 = sum(1 for i,j,k in combs(range(n), 3)
             if (adj[i][j] and adj[j][k] and adj[k][i]) or
                (adj[i][k] and adj[k][j] and adj[j][i]))

    print(f"Transitive T{n}: t3={t3}")
    print(f"  #3-edge-paths (4 vertices) = {count_3path}")
    print(f"  Expected C({n},4) + 2*(n-3)*t3 = {comb(n,4)} + {2*(n-3)*t3} = {comb(n,4) + 2*(n-3)*t3}")
    print(f"  Expected C({n},4) + 2*(n-2)*t3 = {comb(n,4)} + {2*(n-2)*t3} = {comb(n,4) + 2*(n-2)*t3}")
    print(f"  #4-edge-paths (5 vertices) = {count_4path}")
    print(f"  H(T{n}) = {count_4path}")

# Check with a non-transitive tournament
print()
# Tournament with t3=1: flip one edge
adj = [[1 if i < j else 0 for j in range(n)] for i in range(n)]
adj[0][1] = 0; adj[1][0] = 1  # flip 0->1 to 1->0, creating cycle 1->0->2->1? No...
# Let me create one with a known 3-cycle: 0->1, 1->2, 2->0 (plus 0->3,0->4,1->3,1->4,2->3,2->4,3->4)
adj2 = [[0]*n for _ in range(n)]
# 3-cycle on {0,1,2}: 0->1, 1->2, 2->0
adj2[0][1] = 1; adj2[1][2] = 1; adj2[2][0] = 1
# 3,4 beat by 0,1,2: 0->3, 0->4, 1->3, 1->4, 2->3, 2->4
adj2[0][3] = 1; adj2[0][4] = 1
adj2[1][3] = 1; adj2[1][4] = 1
adj2[2][3] = 1; adj2[2][4] = 1
# 3->4
adj2[3][4] = 1

t3_2 = sum(1 for i,j,k in combs(range(n), 3)
           if (adj2[i][j] and adj2[j][k] and adj2[k][i]) or
              (adj2[i][k] and adj2[k][j] and adj2[j][i]))

count_3path_2 = sum(1 for a,b,c,d in perms(range(n), 4)
                    if adj2[a][b] and adj2[b][c] and adj2[c][d])

print(f"Tournament with t3={t3_2}:")
print(f"  #3-edge-paths = {count_3path_2}")
print(f"  C(5,4) + 2*2*t3 = {comb(5,4) + 4*t3_2}")
print(f"  C(5,4) + 2*3*t3 = {comb(5,4) + 6*t3_2}")

# The n-2 vs n-3 question:
# In a 4-vertex tournament T[S], H(T[S]) = #{3-edge paths in T[S]}.
# For a transitive 4-vertex tournament: how many directed 3-paths?
# 0->1->2->3 is the only HP. So H = 1.
# For 4-vertex tournament with one 3-cycle:
# H can be > 1.

# But sum_S H(T[S]) counts the total over all C(n,4) subsets.
# The t3 coefficient comes from: each 3-cycle {a,b,c} adds 2 to H(T[S])
# for each S containing {a,b,c} (since a 3-cycle adds 2 HPs to the 4-vertex tournament).
# Wait, does it? In T[{a,b,c,d}], if {a,b,c} form a 3-cycle, how many extra HPs?

# 3-cycle a->b->c->a, plus vertex d.
# d can be source, sink, or mixed relative to {a,b,c}.
# If d is sink (all of a,b,c beat d): HPs are a->b->c->d, b->c->a->d, c->a->b->d = 3.
# Transitive 4-tournament has 1 HP. So extra = 2.
# If d is source: same by reversal, 3 HPs.
# If d has mixed relation: fewer HPs?

# Let me count explicitly for a specific case.
print("\n\n4-VERTEX TOURNAMENT HP ANALYSIS")
n4 = 4
for bits in range(1 << (n4*(n4-1)//2)):
    adj4 = t_from_bits(n4, bits)
    t3_4 = sum(1 for i,j,k in combs(range(n4), 3)
               if (adj4[i][j] and adj4[j][k] and adj4[k][i]) or
                  (adj4[i][k] and adj4[k][j] and adj4[j][i]))
    h = sum(1 for p in perms(range(n4))
            if all(adj4[p[i]][p[i+1]] for i in range(3)))
    # Also count 3-edge paths (= H for 4 vertices)
    count = sum(1 for a,b,c,d in perms(range(n4), 4)
                if adj4[a][b] and adj4[b][c] and adj4[c][d])
    # 3-edge path count should equal H since it's a full permutation
    # Wait, perms(range(4), 4) gives all 4! orderings, same as HP count.
    # Actually no, for 4 vertices, a 3-edge path IS a Hamiltonian path.
    # So count = H(T[4]). ✓

print("4-vertex: H(T[S]) as function of t3(S):")
from collections import Counter
t3_to_h = Counter()
t3_count = Counter()
for bits in range(1 << (n4*(n4-1)//2)):
    adj4 = t_from_bits(n4, bits)
    t3_4 = sum(1 for i,j,k in combs(range(n4), 3)
               if (adj4[i][j] and adj4[j][k] and adj4[k][i]) or
                  (adj4[i][k] and adj4[k][j] and adj4[j][i]))
    h = sum(1 for p in perms(range(n4))
            if all(adj4[p[i]][p[i+1]] for i in range(3)))
    t3_to_h[(t3_4, h)] += 1
    t3_count[t3_4] += 1

for key in sorted(t3_to_h):
    print(f"  t3={key[0]}, H={key[1]}: {t3_to_h[key]} tournaments")

# H(4-tournament) = 1 + 2*t3(4-tournament)
# Since at n=4: OCF = 1 + 2*alpha_1 = 1 + 2*t3 (no alpha_2 possible)
# At n=4, t3 can be 0 or 1: transitive (t3=0, H=1) or cyclic (t3=1, H=3).

# So #fwd4path = sum_S H(T[S]) = sum_S [1 + 2*t3(S)]
# = C(n,4) + 2*sum_S t3(S)
# And sum_S t3(S) = t3 * C(n-3, 4-3) = t3 * C(n-3, 1) = t3*(n-3).
# #fwd4path = C(n,4) + 2*(n-3)*t3.

# At n=5: 5 + 2*2*1 = 5 + 4*t3 (for t3=1: 9).
# Let me verify this numerically for the tournament with t3=1 above.
count_check = sum(1 for a,b,c,d in perms(range(5), 4)
                  if adj2[a][b] and adj2[b][c] and adj2[c][d])
print(f"\nVerification: t3={t3_2}, #3-edge-paths = {count_check}")
print(f"  C(5,4) + 2*(5-3)*{t3_2} = {comb(5,4) + 2*2*t3_2}")
print(f"  Match: {count_check == comb(5,4) + 2*2*t3_2}")

# So THM-092's claim of 2(n-2)*t3 was WRONG for fwd4path!
# The correct coefficient is 2*(n-3) not 2*(n-2).
# Let me verify at n=6 and n=7 to be sure.

for n in [6, 7]:
    print(f"\nn={n}: sampling to check 3-edge-path formula")
    import random
    random.seed(42)
    me = n*(n-1)//2
    for trial in range(3):
        bits = random.getrandbits(me)
        adj_n = t_from_bits(n, bits)
        t3_n = sum(1 for i,j,k in combs(range(n), 3)
                   if (adj_n[i][j] and adj_n[j][k] and adj_n[k][i]) or
                      (adj_n[i][k] and adj_n[k][j] and adj_n[j][i]))
        count_n = sum(1 for a,b,c,d in perms(range(n), 4)
                      if adj_n[a][b] and adj_n[b][c] and adj_n[c][d])
        expected = comb(n,4) + 2*(n-3)*t3_n
        print(f"  t3={t3_n}: actual={count_n}, C(n,4)+2(n-3)t3={expected}, match={count_n==expected}")
