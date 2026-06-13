#!/usr/bin/env python3
"""
WHY is Cov(H_sym, H_anti) = 0 exactly?

H_sym(x) = (H(x) + H(~x))/2
H_anti(x) = (H(x) - H(~x))/2

Cov(H_sym, H_anti) = E[H_sym * H_anti] - E[H_sym]*E[H_anti]

Note: E[H_anti] = 0 always (sum over all x of H(x)-H(~x) = 0 by pairing).

So Cov = E[H_sym * H_anti] = E[(H(x)+H(~x))(H(x)-H(~x))/4]
       = E[(H(x)^2 - H(~x)^2)/4]
       = (E[H^2] - E[H^2])/4 = 0

WAIT. This is TRIVIALLY zero! Because H(x)^2 and H(~x)^2 have the same
distribution (~ is just a relabeling of the tiling space).

So the orthogonality is NOT a deep fact about tournaments — it's a trivial
consequence of the GS flip being an involution on a uniform measure space!

Any involution sigma on a finite set X gives:
  f_sym(x) = (f(x) + f(sigma(x)))/2
  f_anti(x) = (f(x) - f(sigma(x)))/2

And Cov(f_sym, f_anti) = 0 always, for ANY f and ANY involution sigma.

Proof: E[f_sym * f_anti] = (1/N) sum_x [(f(x)+f(sigma(x)))(f(x)-f(sigma(x)))]/4
     = (1/4N) sum_x [f(x)^2 - f(sigma(x))^2]
     = (1/4N) [sum_x f(x)^2 - sum_x f(sigma(x))^2]
     = 0 (since sigma is a bijection)

E[f_anti] = 0 similarly. So Cov = 0.

OK so this is trivial. The INTERESTING question is different:
What is the ratio Var(H_sym)/Var(H_anti)?

For n=4: Var(H_sym)/Var(H_anti) = 0.50/1.50 = 1/3
For n=5: 6.17/6.69 ≈ 0.92
For n=6: 67.25/37.31 ≈ 1.80

The ratio shifts from anti-dominated to sym-dominated as n grows.
This IS a non-trivial fact about H.

Actually, the MORE interesting question is: what does the GS flip DO
to the FORWARD-EDGE POLYNOMIAL F(T,x)?

If GS flip reverses all non-backbone edges, and P has k non-backbone
forward edges, then in T_flip, those k edges become backward and the
(steps_nbb - k) backward ones become forward.

Let's be precise. For a permutation P, the step P_i -> P_{i+1} is either:
- A backbone step (|P_i - P_{i+1}| = 1): unchanged by GS flip
- A non-backbone step (|P_i - P_{i+1}| >= 2): forward/backward flipped

So if P has b backbone-forward steps and f non-backbone-forward steps,
and s non-backbone steps total, then:
  fwd(P, T) = b + f
  fwd(P, T_flip) = b + (s - f) = b + s - f

So fwd(P, T) + fwd(P, T_flip) = 2b + s = constant for each P!

This means: FOR EACH PERMUTATION P INDIVIDUALLY,
  fwd(P,T) + fwd(P,T_flip) = 2*bb_fwd(P) + nbb_steps(P)

Let's verify this and see what it implies.
"""
from itertools import permutations
import numpy as np

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = []
    for a in range(n):
        for b in range(a):
            if a - b >= 2:
                tiles.append((a, b))
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A

print("=== Per-permutation invariant under GS flip ===")
for n in [4, 5]:
    tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
    m = len(tiles)

    # For a few tilings, check the invariant
    for bits in [0, 7, (1 << m) - 1]:
        A = tiling_to_tournament(bits, n)
        flip = bits ^ ((1 << m) - 1)
        Af = tiling_to_tournament(flip, n)

        invariant_vals = set()
        for P in permutations(range(n)):
            fwd_T = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
            fwd_Tf = sum(1 for i in range(n-1) if Af[P[i]][P[i+1]])

            bb_fwd = sum(1 for i in range(n-1) if abs(P[i]-P[i+1])==1 and A[P[i]][P[i+1]])
            nbb_steps = sum(1 for i in range(n-1) if abs(P[i]-P[i+1])>=2)

            predicted = 2*bb_fwd + nbb_steps
            actual = fwd_T + fwd_Tf
            if predicted != actual:
                print(f"  MISMATCH at n={n}, bits={bits}, P={P}")
            invariant_vals.add((fwd_T + fwd_Tf, nbb_steps, bb_fwd))

        unique_sums = set(v[0] for v in invariant_vals)
        print(f"  n={n}, bits={bits}: fwd(P,T)+fwd(P,T_flip) values = {sorted(unique_sums)}")

# The key insight: fwd(P,T) + fwd(P,T_flip) depends only on P, not on T!
# (for fixed backbone structure)
#
# More precisely: it depends on which STEPS of P are backbone vs non-backbone.
# The backbone edges i->i-1 are always present (A[i][i-1]=1), so:
# bb_fwd(P) = number of steps P_i -> P_{i+1} with P_{i+1} = P_i - 1

print("\n\n=== Consequence for H ===")
# H(T) = F(T, 0) = number of perms with fwd=0 (all backward)
# Wait, that's wrong. H = total HP count = sum_P indicator(P is a ham path)
#
# Actually: H(T) = sum_P prod_{i=0}^{n-2} A[P_i][P_{i+1}]
# and each A[P_i][P_{i+1}} is either 0 or 1.
# H counts exactly those permutations where ALL edges in the path exist in T.
# This is NOT the same as F(T,0).
#
# F(T,x) = sum_P x^{fwd(P)} where fwd counts ONLY the forward edges among
# the edges that DO exist. Wait, no — F counts forward STEPS among ALL n-1
# steps, where a step P_i->P_{i+1} is "forward" if A[P_i][P_{i+1}]=1.
#
# So: H(T) = sum_{P: all edges present} 1 = F(T,0)? No!
# F(T,0) = sum_P 0^{fwd(P)} = number of perms with fwd=0.
# But perms with fwd=0 means ALL edges backward, i.e., A[P_i][P_{i+1}]=0
# for all i. That's NOT a Hamiltonian path in T!
#
# Actually F(T,x) = sum_P x^{fwd(P)} where EVERY permutation is summed over,
# not just Hamiltonian paths. And fwd(P) counts how many edges P_i->P_{i+1}
# exist in T (direction matters).
#
# So F(T,1) = sum_P 1 = n! (all perms equally weighted)
# F(T,0) = #{P : fwd(P)=0} = #{P : no forward edges}
#
# H(T) = #{P : all edges present in T} = F(T, ???)
#
# Actually no. Let me re-read. The standard forward-edge polynomial:
# For a tournament T and permutation P = (v_0,...,v_{n-1}):
#   fwd(P) = |{i : T has edge v_i -> v_{i+1}}| = |{i : A[v_i][v_{i+1}]=1}|
#
# So fwd(P) + bwd(P) = n-1 for ALL P (since T is a tournament, exactly one
# direction exists for each edge).
#
# H(T) = sum_P [fwd(P) == n-1] = coefficient of x^{n-1} in F(T,x) times
# (normalizing...) Actually:
# F(T,x) = sum_P x^{fwd(P)} = sum_{k=0}^{n-1} D_k * x^k
# where D_k = #{P : fwd(P) = k}
#
# H(T) = D_{n-1} = #{P : all n-1 edges forward} = #{Hamiltonian paths in T}
#
# YES! H = D_{n-1} = F(T, x) evaluated... no. D_{n-1} is the coefficient of
# x^{n-1}, which equals F(T,1) - sum_{k<n-1} D_k... that's circular.
#
# But by palindromicity: D_k = D_{n-1-k}. So D_{n-1} = D_0.
# D_0 = #{P : fwd(P) = 0} = #{P : all edges backward}
# = #{Hamiltonian paths in T^op (reversed tournament)}
#
# So H(T) = D_{n-1} = D_0 = H(T^op). This is the palindrome consequence.
#
# Now: what happens to D_k under GS flip?
# Under GS flip, for each P:
#   fwd(P, T_flip) = bb_fwd(P) + nbb_steps(P) - nbb_fwd(P, T)
#                   = bb_fwd(P) + nbb_steps(P) - (fwd(P,T) - bb_fwd(P))
#                   = 2*bb_fwd(P) + nbb_steps(P) - fwd(P,T)
#
# Let c(P) = 2*bb_fwd(P) + nbb_steps(P) (constant for P, independent of T)
# Then fwd(P, T_flip) = c(P) - fwd(P, T)
#
# So the FORWARD-EDGE POLYNOMIAL transforms as:
# F(T_flip, x) = sum_P x^{c(P) - fwd(P,T)}
#              = sum_P x^{c(P)} * x^{-fwd(P,T)}
#              = sum_P x^{c(P)} * (1/x)^{fwd(P,T)}
#
# This is NOT a simple substitution because c(P) varies with P!

# Let's compute c(P) distribution
print("Distribution of c(P) = 2*bb_fwd(P) + nbb_steps(P):")
for n in [4, 5, 6]:
    c_vals = []
    bb_vals = []
    nbb_vals = []
    for P in permutations(range(n)):
        # bb_fwd = steps where P_{i+1} = P_i - 1 (backbone edge, always forward)
        bb_fwd = sum(1 for i in range(n-1) if P[i+1] == P[i] - 1)
        nbb_steps = sum(1 for i in range(n-1) if abs(P[i]-P[i+1]) >= 2)
        c = 2*bb_fwd + nbb_steps
        c_vals.append(c)
        bb_vals.append(bb_fwd)
        nbb_vals.append(nbb_steps)

    from collections import Counter
    c_dist = Counter(c_vals)
    print(f"\n  n={n}: c(P) distribution: {dict(sorted(c_dist.items()))}")
    print(f"    bb_fwd range: [{min(bb_vals)}, {max(bb_vals)}]")
    print(f"    nbb_steps range: [{min(nbb_vals)}, {max(nbb_vals)}]")

    # Is c(P) always = n-1?
    if len(c_dist) == 1:
        print(f"    c(P) is CONSTANT = {list(c_dist.keys())[0]}!")
    else:
        print(f"    c(P) is NOT constant")

# If c(P) = n-1 for all P, then fwd(T_flip) = (n-1) - fwd(T)
# and F(T_flip, x) = x^{n-1} * F(T, 1/x)
# which, combined with palindromicity, would give F(T_flip, x) = F(T, x)!
# That would mean GS flip preserves the forward-edge polynomial!
# And hence H(T) = H(T_flip), meaning H_anti = 0 always.
# But we KNOW H_anti ≠ 0, so c(P) can't be constant.
