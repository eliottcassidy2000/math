# Signed Adjacency Identity (new proof angle)

**Discovered by:** opus-2026-03-05-S4
**Status:** Verified as polynomial identity at n=3,4,5; equivalent to even-odd split and OCF

---

## Statement

For any tournament T on n vertices with arc i→j, define:

- **Position generating functions:**
  F(x) = sum_k adj_k(i→j, T) * x^k  (T-paths using i→j, grouped by position)
  G(x) = sum_k adj_k(j→i, T') * x^k (T'-paths using j→i, grouped by position)

  where adj_k counts paths with exactly k W-vertices before the {i,j} segment.

- **Signed adjacency counts:**
  F(-1) = "signed T-adjacency count"
  G(-1) = "signed T'-adjacency count"

**Identity:** F(-1) = G(-1) for all tournaments T and all arcs i→j.

Equivalently: D(x) = F(x) - G(x) has x = -1 as a root, so D(x) = (1+x)*Q(x).

**This is a POLYNOMIAL IDENTITY** in the arc variables (verified with real-valued arcs at n=3,4,5).

---

## Equivalences

The signed adjacency identity is EQUIVALENT to:

1. **Even-odd split lemma:** sum_{|S| even} Delta(S,R) = sum_{|S| odd} Delta(S,R)
2. **Alternating sum vanishing:** sum_S (-1)^|S| Delta(S, W\S) = 0
3. **Alternating subset convolution symmetry:** B(L_i, R_j) = B(L_j, R_i)
   where B(f,g) = sum_S (-1)^|S| f(S) g(W\S)
4. **Combined with delta_H = delta_I (THM-015), proves OCF for all n where both hold.**

---

## Key Properties

1. **x = -1 is the UNIQUE universal root** of D(x). At all other x values, D(x) != 0 generically.

2. **Proof for n=3 (one line):**
   D(x) = -(s_w)(1+x) where s_w = 1 - T[w][i] - T[j][w].
   Proof: h_end({v,w}, v) + h_start({v,w}, v) = T[w][v] + T[v][w] = 1 for all v,w.

3. **Proof for n=4 (algebraic):**
   D(x) = (1+x)*Q(x) where Q = [-s_{W[pi(0)]}, -s_{W[pi(1)]}] and pi depends on internal arc.
   Proved by hand in algebraic_proof_n4.py.

4. **Polynomial identity nature:**
   Since D(-1) is multilinear in the arc variables, and verified to vanish at all {0,1}^k inputs for n<=8,
   it is identically zero as a polynomial for each n<=8.

---

## Structural Analysis

**Swap involution interpretation:**
- Matched paths (swap works): contribute equally to F(-1) and G(-1). Cancel in D(-1).
- Unmatched T-paths: blocked by type C vertices (s_w = -1, i.e., w→i and j→w)
- Unmatched T'-paths: blocked by type D vertices (s_w = +1, i.e., i→w and w→j)
- The identity reduces to: signed count of C-blocked paths = signed count of D-blocked paths.

**What DOESN'T work for a proof:**
- Per-vertex decomposition: C(u,v) + C(v,u) != 0 in general
- Type-signature determination: B_signed depends on internal arcs, not just (p_w, q_w)
- Vertex deletion reduction: D_L and D_R don't cleanly relate to (n-1)-case

**What MIGHT work:**
- Algebraic proof exploiting multilinearity and the complementary structure of i↔j interface arcs
- Transfer matrix / generating function approach
- The fact that T[w][v] + T[v][w] = 1 is the KEY algebraic identity (proved n=3; at n>=4, internal arcs complicate)

---

## Key Algebraic Reduction (NEW)

The identity B(Li, Rj) = B(Lj, Ri) is equivalent to the SIGMA-INVARIANCE of B:

Define sigma: p_w -> 1-q_w, q_w -> 1-p_w (keeping internal arcs fixed).
This negates s_w and preserves d_w = p_w - q_w.

**CONFIRMED:** B(Li, Rj) is an EVEN FUNCTION of the s-variables.
Reparametrize p_w = (1-s_w+d_w)/2, q_w = (1-s_w-d_w)/2.
Then B(+s, d, t) = B(-s, d, t) for all (s, d, internal arcs).

Since each monomial in B uses at most 2 interface arcs (one from the left boundary p_{w_k},
one from the right boundary q_{u_1}), the maximum s-degree is 2. Therefore:

**The evenness condition reduces to: all s-degree-1 coefficients vanish.**

WARNING: This does NOT mean C_w + D_w = 0 everywhere. C_w + D_w is a function of
all arc variables including other s-values. The correct statement is:

**alpha_w := dB/ds_w |_{s=0} = 0 for each w in W.**

At n=4: C_a + D_a = s_b (zero only at s=0). The s-polynomial is B = B_0(d,t) - sa*sb/2.
At n=5: all s-linear coefficients verified zero (< 1e-11). Quadratic terms are only
CROSS-terms s_w*s_u (w != u); diagonal s_w^2 terms are zero.

Verified computationally at n=3,4,5,6,7: B(Li,Rj) = B(Lj,Ri) with real-valued arcs in [-2,3].

**Subsidiary identity:** B(Li, Lj) = (-1)^m * B(Ri, Rj), equivalently B(Lj, Li) = B(Ri, Rj).

---

## Proof Structure (INDUCTIVE)

The identity B(Li,Rj) = B(Lj,Ri) can be proved by INDUCTION on |W| = m:

**Base cases:** m=0 (trivial), m=1 (complementarity: T[a][b]+T[b][a]=1), m=2 (hand algebra).

**Inductive step:** Assume the identity for |W| < m. The difference D = B(Li,Rj)-B(Lj,Ri)
is LINEAR in the s-variables (s_w = 1-p_w-q_w). The coefficient of s_v decomposes into:

1. **d-dependent part:** sum_{u!=v} (d_u/2)(left(u)-right(u)). By the inductive hypothesis
   applied to sub-tournaments on W\{u,v} (size m-2), left(u) = right(u), so this VANISHES.

2. **d-independent part:** boundary + d_indep = 0, which reduces to:

**CLAIM (B):** For any tournament on m vertices W with vertex v:
  sum_{S ⊆ W\{v}} (-1)^|S| H(S) h_start(W\S, v) = (-1)^{m+1} h_end(W, v)
where H(S) = total Hamiltonian path weight on S (H(∅) = 1).

**STATUS of Claim (B):**
- Proved algebraically at m=1,2,3
- Verified computationally at m=1,...,8 (polynomial identity in arc variables)
- This is a STANDALONE identity about tournament Hamiltonian paths (no interface variables)
- Does NOT involve vertices i,j at all — purely about internal tournament structure

If Claim (B) is proved for all m, the full induction closes, proving OCF for all n.

---

## Connection to OCF Proof

The signed adjacency identity IS the even-odd split, which IS equivalent to OCF.
Proving Claim (B) for all m proves the main conjecture (Claim A).

The proof chain: Claim (B) for all m → B(Li,Rj)=B(Lj,Ri) for all n → OCF for all n → Claim A.
