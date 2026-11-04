# One-child modulo-four separation of quotient converses (THM-4864)

kind-pasteur-2026-11-04-S9.

**Status.**
- PROVED: Lemmas 1–3, Proposition 4, Theorems A and B, and the Corollary, given the cited reviewed inputs.
- FINITE-EXACT: the pre-registered confirmations in `PREREGISTRATION.md` (outcomes appended there).
- An independent audit is recorded in `audit/REVIEW.md` once it is complete.

## 0. Setting and the question answered

All names are physical. For a named tournament T the ports are attached with x->y fixed, and a_v=[v->x], c_v=[v->y].
F_T(a,c) is the Hamilton count of the augmented tournament.

A core is W = Q[A_1,...,A_m]. It is compared with W' = Q^op[A'_1,...,A'_m] on the same child sets, where:
- the A'_i are arbitrary tournaments with |A'_i| = |A_i|;
- no profile equality is needed for anything below.

**The question.** Session189 (`session189_response_frontier_20261104/collision/PROOF.md`, reviewed) proves the
following. In a full-response collision, a transitive child of order >= 5 at a highest changed prime node can only be
**unchanged or wholly reversed**. The remaining binary cube then has to be separated.

**The answer.** Below, both cases are separated modulo four by explicit actual queries, at every order and in every
fixed ambient. The same mechanism also separates every quotient converse with a singleton child, from coarse queries
alone.

## 1. Odd-cycle parity and the weighted port-cycle polynomial

**Odd-cycle expansion.** Modulo four the classical odd-cycle Hamilton expansion reads
H(R) = 1 + 2·(number of directed odd cycles of R) (mod 4).
- Source: Grinberg–Stanley, arXiv:2307.05569, Theorem 1.39.
- It is used in exactly this form in `session188_structural_bridges_20261104/parity_profile_transfer/PROOF.md` and
  in `odd_pair_quotient/PROOF.md`.

**The polynomial g_Z.** Let Z be a tournament whose vertices carry weights w in F2; the ports weigh 1.
g_Z(a,c) in F2 is the number, modulo two, of directed cycles of the augmented tournament Z(a,c) that pass through at
least one port and have odd total weight. It is treated as its unique multilinear polynomial. Enumerating port
positions gives:

    g_Z = sum_{P: w(P) even} [(1+a_s)a_e + (1+c_s)c_e]  +  sum_{P: w(P) odd} (1+c_s) a_e
        + sum_{(P,P') disjoint, w(P)+w(P') odd} (1+a_s) c_e (1+c_s') a_e'              (1)

Here P runs over directed paths of Z (singletons included), with start s and end e, and P' has start s' and end e'.
The cycles behind each term are:
- x→P→x, for the first a-term;
- y→P→y, for the first c-term;
- x→y→P→x, for the second sum;
- x→P→y→P'→x, for the double sum.

These are all of them, since y→x is impossible.

**Lemma 2 (coefficients).** Extracting coefficients from (1), for u≠v:

    [a_u c_v] g_Z     = N_odd(v⇝u) + Pi(→v,→u),
    [a_u a_v c_v] g_Z = #{P ⊆ Z−v ending at u : w(P)+w(v) odd},

all modulo two. The terms are:
- N_odd(v⇝u) counts odd-weight paths from v to u.
- Pi(→v,→u) counts ordered pairs (P, P') of disjoint paths, with P ending at v, P' ending at u, and odd total weight.

**Converse duality.** Reversing every cycle and exchanging x with y gives g_{Z^op}(a,c) = g_Z(1+c, 1+a). So every
coefficient of g_{Z^op} is the same formula with paths reversed.

**Lemma 1 (contraction modulo four).** Let A_1,...,A_r be disjoint modules of W, and let a query be uniform on every
A_i. Then the number modulo two of odd port cycles of the augmented W equals g_Z. Here Z contracts each A_i to one
vertex of weight [|A_i| = 1], and every other vertex keeps weight 1. Hence
(F_W − F_W')/2 ≡ const + g_Z − g_Z' (mod 2), with a constant independent of the query, when W' has the same module
sets.

*Proof.*
1. **Several segments cancel.** Suppose a port cycle meets A_i in r >= 2 maximal segments. Permuting the segment
   contents among their slots preserves every boundary incidence, by the module property and uniform ports. The cycle
   has a distinguished port, and its segments are disjoint and nonempty, so S_r acts freely. Such cycles therefore
   cancel modulo two.
2. **One segment.** A single segment through A_i is a Hamilton path of some A_i[S] with S nonempty. Its number of
   realizations of a given length parity is the sum of H(A_i[S]) over the S of that parity. Since H(A_i[S]) is odd
   (Rédei), this is congruent modulo two to the number of nonempty subsets of that parity.
   - For |A_i| >= 2: 2^{|A_i|−1}−1 ≡ 1 for even |S|, and 2^{|A_i|−1} ≡ 0 for odd |S|.
   - For |A_i| = 1: only odd |S|.
   So a contracted vertex behaves as weight 0 or 1 respectively.
3. **Port-free cycles** contribute a query-independent constant. ∎

This recovers the reviewed all-nontrivial coarse law (2) of `parity_profile_transfer`.

**Lemma 3 (endpoint parity).** For every tournament R and every v:
- the number of all paths starting at v, singleton included, is congruent to E_v(R) modulo two;
- the number of all paths ending at v is congruent to S_v(R).

Here S_v and E_v count Hamilton starts and ends. *Proof:* this is `parity_profile_transfer` §2. Induct on |R|: the
number of ending paths at v is 1 + sum_{w→v} (number of ending paths at w in R−v), which is ≡ 1 + sum_{w→v}
S_w(R−v) = sum_{v→w} S_w(R−v) = S_v(R), because H(R−v) is odd. Converse gives the other identity.

Consequence: with eps_v(S) = S_v(R[S]) + E_v(R[S]), we have sum_{v ∈ S ⊊ R} eps_v(S) ≡ 0 for every R and v.

## 2. Weighted E1

**Proposition 4.** Let Z be weighted on n >= 3 vertices, with delta_Z = g_Z − g_{Z^op} (same weights). For every v,

    L_v := sum_{u≠v} [a_u c_v] delta_Z + sum_{v→u} [a_u a_v c_v] delta_Z ≡ w(v)  (mod 2).

With w ≡ 1 this is the incoming E1 identity (`two_port_parity_separation_20261104/PROOF.md`).

*Proof.* Use Lemma 2 for Z and for Z^op.

1. **The quadratic part.**
   - The N-terms give O_v: the number of odd-weight paths with >= 2 vertices that start at v, plus the number that end
     at v.
   - By Rédei, the Pi-terms give sum_{S∋v} eps_v(S)·f(S), where f(S) = #{nonempty S' ⊆ V∖S : w(S') ≡ 1+w(S)} mod 2.
2. **The cubic part.**
   - Prepending v to the paths that start at u (with v→u) gives the odd-weight paths starting at v.
   - For the paths that end at u, write sum_{v→u} = sum_{all u} − sum_{u→v}.
     - Appending v turns the second sum into the odd-weight paths ending at v.
     - By Rédei the first sum is g(v) = #{nonempty S ⊆ V−v : w(S) ≡ 1+w(v)} mod 2.
   - Together the cubic part is O_v + g(v).
3. **Combining.** Hence L_v ≡ g(v) + sum_{S∋v} eps_v(S) f(S).
4. **The value of g(v).** g(v) ≡ w(v) for n >= 3: subsets of V−v are equidistributed by weight parity when V−v
   contains a weight-1 vertex, and otherwise all have weight 0.
5. **The correction vanishes.** f(V) = 0. For S ≠ V, f(S) ≡ w(S) + omega(S)·[|V∖S| = 1], where omega(S) is the
   indicator that V∖S contains a weight-1 vertex. By Lemma 3, for x ≠ v,
   sum_{S∋v,x; S≠V} eps_v(S) ≡ eps_v(V∖x), while sum_{v∈S⊊V} eps_v(S) ≡ 0. Therefore both
   sum_{S≠V} eps_v(S)·w(S) and sum_S eps_v(S)·omega(S)[|V∖S|=1] equal sum_{x≠v} w(x)·eps_v(V∖x), and they cancel. ∎

The failure boundary is n = 2: there g(v) can differ from w(v).

## 3. Singleton children

**Theorem A.** Let m >= 3. If some child A_z of W is a singleton, then some child-uniform actual query has
F_W − F_W' ≡ 2 (mod 4). This holds for every tournament Q and all same-size target interiors.

*Proof.*
1. On child-uniform queries, Lemma 1 contracts W to Q with weights [|A_i|=1], and W' to Q^op with the same weights.
2. Modulo two and up to a constant, the half-contrast is delta_Q.
3. Proposition 4 at z gives L_z = 1. Since L_z is a sum of nonconstant coefficients, some child-uniform query has an
   odd half-contrast. ∎

L_z needs at most 4(m+1) distinct coarse queries.

## 4. Transitive children: unchanged or reversed

**Theorem B.** Let A_j be transitive with order k >= 2 and physical order b_1→...→b_k. Let u be b_1 or b_k. Suppose u
is also the source or the sink of A'_j, which A'_j may be any transitive tournament on the same names satisfying.
Then some actual query of the following form has F_W − F_W' ≡ 2 (mod 4):
- uniform on every other child;
- uniform on beta = A_j − u;
- with arbitrary bits at u.

In particular **A'_j = A_j (unchanged) and A'_j = A_j^op (reversed) are both separated.**

*Proof.*
1. **Kernel.** beta is a module of A_j and of A'_j, because u is an end of both orders. Lemma 1 contracts the other
   children, giving them weight [size 1], and contracts beta, giving it weight w_beta = [k=2]. Take u = b_1; the case
   u = b_k is the mirror image. The kernel is Z = Q[u→beta] with w(u) = 1.
2. **Reversed type (u becomes the sink of A'_j).** The target kernel is Q^op[beta→u] = Z^op. Proposition 4 at u gives
   L_u = 1, since n = m+1 >= 3.
3. **Unchanged type (u stays the source).** The target kernel is Z_u = Q^op[u→beta].
   - **(a) Swap-orbit equality, for every k and every weight.**
     - Let Y = Q[beta→u] = Z with the u–beta arc reversed. By converse duality, g_{Z_u} = g_Y∘J, where
       J(a,c) = (1+c, 1+a). Hence the half-contrast is delta_Z + (g_Z − g_Y)∘J.
     - Every port cycle of Z through the arc u→beta corresponds to the cycle of Y with u and beta swapped. Interior
       incidences agree by the module property. The vertex set, and hence the weight, is unchanged.
     - Hence g_Z − g_Y = sum (h − h∘tau), where tau exchanges the query variables of u and beta. The same holds after
       composing with J, since J commutes with tau. So the half-contrast and delta_Z have **equal sums over every
       tau-orbit of monomials**.
   - **Some other child is a singleton z.** The monomial set of L_z is a union of tau-orbits, because beta and u have
     the same incidence to z. So L_z(half-contrast) = L_z(delta_Z) = w(z) = 1, by (a) and Proposition 4.
   - **Otherwise (hard case).** Every other weight is 0.
     - (b) **k >= 3.**
       - The monomial set of L_u, together with that of L_beta and the single monomial a_u a_beta c_beta, is a union of
         tau-orbits. By (a), its sum on the half-contrast equals its sum on delta_Z, which is
         w(u) + w(beta) + h_ub = 1 + 0 + h_ub by Proposition 4.
       - By Lemma 2, h_ub counts the paths of Z−beta through u at either end. Every such path has odd weight, so by
         Lemma 3 h_ub ≡ S_j(Q) + E_j(Q) = p_j(Q).
       - Hence if p_j(Q) = 0, an anchored query separates.
       - If p_j(Q) = 1, the coarse linear coefficient [alpha_j] equals p_j(Q) by Lemma 3, so a coarse query separates.
         This is also `parity_profile_transfer` (1).
     - (c) **k = 2 (w_beta = 1).**
       - The half-contrast is delta_Z + (g_{Z^op} − g_{Z_u}). Z^op and Z_u differ only in the u–beta arc.
       - A port cycle using that arc has weight equal to its number of ports modulo two, because u and beta weigh 1 each
         and the other vertices weigh 0. So it is odd only with exactly one port, and contributes pure-a or pure-c
         monomials.
       - L_u consists of mixed monomials only, so L_u(half-contrast) = L_u(delta_Z) = 1. ∎

**Battery.** At most 9(m+1) distinct anchored queries realize L_u, L_beta and the extra cubic. Each kernel bit expands
uniformly onto its child, and onto beta.

## 5. Ambient lift and the original-module frontier

**Lifting.** Every core query above with residue 2 modulo 4 lifts to every fixed ambient when both ports copy the core's
exterior incidences. Source: `parity_profile_transfer` §3 (reviewed). The universal cover residues q_2, q_3 and
q_{>=4} modulo 4 depend only on the augmented order, and the outer Hamilton coefficient is odd.

**Corollary.** Suppose T differs from an original-module XOR target U but F_T = F_U. Then at every highest changed
original strong node M, with prime quotient Q of order >= 7 (S188 prime_seven_frontier):
1. no child of M is a singleton (Theorem A);
2. p(Q) = 0 (`parity_profile_transfer` §3);
3. no child is transitive of order 2 or 3;
4. no child is transitive of order >= 5;
5. every transitive child of order 4 is mapped in U to an order whose two ends are exactly its two middle names.

So every child of every highest changed prime node is **non-transitive, or a TT4 with that crossed reordering**.

*Proof.*
1. **Highest-node reduction and restoration.** Use `highest_quotient_reduction` and `prime_seven_frontier` §3, both
   reviewed. Compare T with V, where V changes only the interior of M. The query's ports copy M's exterior incidences,
   and every other highest node is restored by its equal complete profile.
2. **Transitive children stay transitive.** Hereditary equality (`module_decomposition`) keeps U[A] transitive when
   T[A] is.
3. **Order 2.** A two-name transitive U[A] is A or A^op.
4. **Order 3.** Every reordering of TT3 has b_1 or b_3 at an end.
5. **Order >= 5.** Session189 §5 forces U[A] ∈ {A, A^op}.
6. Theorem B, with the lift, separates each of these cases, a contradiction. ∎

## 6. Scope

- **Not claimed.**
  - No unrestricted full-response injectivity.
  - No separation for non-transitive children, or for TT4 crossed reorderings. The FINITE-EXACT evidence for those is
    in `PREREGISTRATION.md`, with 0 blind cases.
  - No collision at any order.
  - No statement for linear nodes beyond the reviewed S188 results.
- **Integer gaps.** The modulus 4 is essential to the proofs above. Integer gaps are not addressed.
- **Sharpness.** The S188 prime-eight endpoint-obstruction quotient is blind to all coarse queries modulo 4, yet it is
  separated by anchored queries at every TT_k child: P4. So the anchor, not the coarse layer, carries the separation.
