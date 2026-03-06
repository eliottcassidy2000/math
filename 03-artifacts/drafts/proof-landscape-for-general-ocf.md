# Proof Landscape for General OCF: H(T) = I(Omega(T), 2)

**Status:** PROVED for n <= 7 (exhaustive), VERIFIED for n <= 10 (sampling)
**Author:** opus-2026-03-05-S3
**Goal:** Prove for ALL n

---

## The Identity

For any tournament T on n vertices:
  H(T) = I(Omega(T), 2)

where H(T) = #{Hamiltonian paths} and I(Omega(T), x) = independence polynomial of the
odd-cycle conflict graph, evaluated at x=2.

## Equivalent Formulations

1. **Arc-flip induction (THM-015):** For any arc i->j in T, with T' = flip(T, i, j):
   adj(i,j) - adj'(j,i) = Delta_I
   (Suffices because any T is reachable from the transitive tournament by arc flips.)

2. **Vertex-deletion recursion (Claim A):** For any vertex v:
   H(T) - H(T-v) = 2 * sum_{C containing v} mu(C)

3. **Swap involution (THM-014):**
   #U_T - #U_T' = Delta_I
   where U_T = unmatched T-paths under the swap map.

## What's Known

- **Claim B is PROVED:** I(Omega(T),2) - I(Omega(T-v),2) = 2*sum mu(C).
- **Formulation 1 is PROVED at n<=7** by exhaustive symbolic verification.
- **Delta_I has a closed form (THM-013):** Delta_I = sum_{k>=1} 2^k Delta(alpha_k).
  Simplified: Delta_I = -2*sum_x s_x*H(B_x) + 2*(D5-C5) + higher corrections.

## Approaches Tried and Their Status

### A. Per-vertex decomposition (DEAD END)
Attempted to show that for each vertex x, the contribution to U_T matches s_x*H(B_x).
**Result:** The ratio nadj(x,i,j)/H(B_x) is NOT constant. Varies from 2 to 7+.
The identity holds only globally, not vertex by vertex.

### B. Subset convolution framework (PROMISING but incomplete)
adj(i,j) = sum_{S subset V_0} f_i(S) * g_j(V_0\S)
adj'(j,i) = sum_{S} f_j(S) * g_i(V_0\S)
where f_i, g_j depend on {T[x][i]} and {T[j][x]} (s-value structure).
**Status:** Correct framework. Both sides are multilinear polynomials.
The key is showing the convolution difference equals Delta_I.
**Obstacle:** No known way to simplify the convolution into cycle-counting terms.

### C. Transfer matrix / permanent (EXPLORED, no result)
H(T) is NOT the permanent of the adjacency matrix. No clean matrix expression found.
The DP-based computation (transfer through subsets) doesn't simplify.

### D. Bijective proof (SUGGESTIVE but unclear)
I(Omega,2) counts "2-colored VD odd cycle collections."
Each such object should biject to a Ham path.
At n=3: all 3 paths have the cycle consecutive, so the coloring doesn't obviously
distinguish them. The bijection, if it exists, is non-trivial.

### E. Induction on n via vertex deletion (STANDARD approach, stuck)
Inductive step: H(T) - H(T-v) = 2*sum mu(C).
By induction: H(T-v) = I(Omega(T-v), 2). And Claim B gives I(Omega(T),2) - I(Omega(T-v),2).
Need: H(T) - H(T-v) = I(Omega(T),2) - I(Omega(T-v),2), i.e., #{TypeII} = sum mu(C).
At n<=5: trivial (mu=1 always). At n>=6: requires understanding how Type II positions
relate to cycle mu-weights. No per-path formula works.

### F. Arc-flip induction + strong induction (CURRENT BEST)
Assume OCF for < n vertices. Prove Delta_H = Delta_I for n-vertex tournaments.
THM-013 gives Delta_I. THM-014 gives Delta_H = U_T - U_T'.
The complement terms in Delta_I use OCF for (n-3), (n-5), ... vertex sub-tournaments
(by induction). So the inductive hypothesis is fully used.
**Obstacle:** Need to show the swap involution's unmatched count difference
equals the cycle-counting formula. This is the CORE DIFFICULTY.

## Key Structural Facts for a General Proof

1. s_x = 1 - T[x][i] - T[j][x] classifies vertices:
   - s=-1: blocks T-paths only (via pred or succ position)
   - s=+1: blocks T'-paths only
   - s=0: blocks neither (swap always valid)

2. At most one affected cycle in any independent set (A-clique property).

3. Affected cycles all contain {i,j}; their complements are unchanged by the flip.

4. The identity holds as a MULTILINEAR polynomial in arc variables.
   Verifying on all {0,1} assignments proves it for fixed n.
   But 2^{C(n,2)} grows too fast for exhaustive checking at large n.

## Most Promising Path Forward

The subset convolution identity:
  sum_S [f_i(S)*g_j(R) - f_j(S)*g_i(R)] = sum_{k>=1} 2^k Delta(alpha_k)

Both sides are multilinear. The RHS has a recursive structure (alpha_k terms).
A proof likely needs to:
1. Express f_i, g_j, f_j, g_i in terms of sub-tournament Ham path counts
2. Use the inductive hypothesis (OCF for < n) to convert to independence polynomials
3. Show the resulting expression telescopes to the cycle-counting formula

Alternatively: find a SIGN-REVERSING INVOLUTION on the LHS that cancels everything
except the cycle terms on the RHS.

## Computational Targets

- n=8 proof: 2^27 = 134M cases. Feasible with C/Cython optimization (~hours).
- n=9 proof: 2^35 = 34B cases. Needs HPC or clever partitioning.
- n>=10: impractical by exhaustive enumeration.
