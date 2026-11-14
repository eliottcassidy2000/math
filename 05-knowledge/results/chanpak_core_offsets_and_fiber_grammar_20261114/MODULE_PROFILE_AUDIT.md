# Independent audit: finite responses for the modular path-cover state

**Unpublished.** This is a mathematical verification, not a claim of literature priority. The objects counted here are tournament path covers and Hamiltonian paths, not tournament isotypes or rooted metagraph edges.

## Setting and known inputs

For a nonempty tournament T of order n, let pc_T(r) count its spanning covers by r unordered nonempty directed paths. Write F_T(r)=r! pc_T(r), P_T(q)=sum_r pc_T(r) (q)_r, with (q)_r the falling factorial. In particular H(T)=F_T(1)=P_T(1). Use the known odd-cycle formula
`P_T(q)=q^n+sum_{j=1}^d w_j q^{n-2j}` with `d=floor((n-1)/2)`; each w_j is the sum of 2^|K| over all vertex-disjoint directed odd-cycle collections K of total excess `sum_{C in K} (|C|-1)=2j`. Thus every w_j is an even nonnegative integer.

The modular walk-content kernel in [THM-3121](https://github.com/eliottcassidy2000/math/blob/8f3e642bd29de5bfbb7b0bc1fa558fbfe18b6b90/01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md) implies that full pc equality persists under any tournament substitution in which the differing block occupies one quotient slot and all other blocks are identical. Its C3 special case gives quotient walk multiplicity 3 when all three positive visit counts are equal, 1 when their range is 1, and 0 otherwise. [THM-4900](https://github.com/eliottcassidy2000/math/blob/8f3e642bd29de5bfbb7b0bc1fa558fbfe18b6b90/01-canon/theorems/THM-4900-gate-and-doubling-in-the-marked-resolvent-golden-towers-and-tower-ocf.md) credits and uses the known transformation for the single-vertex gate Γ=C3[T,1,1]: `P_(ΓT)=LP_T`, where `(LP)(q)=q(P(q+1)+qP(q)-P(q-1))`.

## Finite response theorem: pass

At fixed n with d≥1, the first d gate responses H(Γ^r T), r=0,...,d−1, determine P_T, hence the full path-cover profile. Indeed, expanding `(L^rP_T)(1)` permits steps +1,0,−1 in the integer evaluation argument. Every branch that encounters q=0 before its last step vanishes from the factor q, and q=0 at the last step vanishes because P_T(0)=0. The remaining terms involve only the samples P_T(1),...,P_T(r+1). Reaching r+1 requires r consecutive +1 steps, so its coefficient is r!. The positive samples up to d therefore follow by triangular back-substitution. Then `(P_T(q)-q^n)/q^(n-2d)` is a polynomial of degree at most d−1 in q², and P_T(1),...,P_T(d) give its values at distinct squares. The falling factorial polynomials form a triangular basis, so P_T determines pc_T. If d=0 (n=1 or 2), the order already determines P_T.

Together with the sufficiency of THM-3121, this characterizes the coarsest **unmarked, fixed-order module-substitution** equivalence that preserves H in all contexts. Any nested context used here is itself a valid one-hole substitution context, because the original block remains a module after each enclosing substitution. The theorem says nothing about contexts that attach to individual vertices, preserve the automorphism group, or preserve arc-reversal neighborhoods.

Independent non-nested proof: set `R_0=H(T)` and `R_k=H(C3[T,TT_k,TT_k])` for k≥1. Only visit counts r≤k+1 in T are possible. At r=k+1 the other two counts must both be k, giving exactly coefficient (k!)² on F_T(k+1) because pc_(TT_k)(k)=1. Thus R_0,...,R_(d−1) determine F_T(1),...,F_T(d). These in turn determine P_T(1),...,P_T(d), since (q)_r=0 for r>q. Use the same interpolation. This corroborates the conclusion without assuming the gate operator.

## First obstruction to scalar H

Number the edges (u,v), 0≤u<v≤4, lexicographically starting at bit 0; set bits orient v→u, clear bits orient u→v. At order 5 the masks 69 and 76 give H=15 in both cases. Counting all 120 permutations by bad adjacent arcs gives:

| mask | reversed pairs | zero bad | exactly one bad | F_T(2) | H(ΓT) |
|---:|---|---:|---:|---:|---:|
|69|(0,1), (0,3), (1,4)|15|18|78|123|
|76|(0,3), (0,4), (1,4)|15|30|90|135|

An ordered 2-cover corresponds precisely to a permutation and a cut at one gap so that neither piece has a bad adjacent arc. Therefore `F_T(2)=4·(#zero bad)+(#one bad)`. The C3 kernel gives `H(ΓT)=3H(T)+F_T(2)`. The resulting power basis polynomials are respectively q⁵+8q³+6q and q⁵+10q³+4q. This is an obstruction at the smallest possible order for fixed-order scalar H: orders ≤4 have d≤1, so n and H determine pc.

## Fixed-H bound and truncated tests

If H(T)=h, then sum_j w_j=h−1. In fact, if w_j>0, take a contributing packing K with pairwise disjoint vertex sets V_i and cycle lengths 2ℓ_i+1, where sum ℓ_i=j. Inside T[V_i] there are at least 2ℓ_i+1 Hamiltonian paths, from opening the chosen cycle at any edge. By positivity of the odd-cycle formula, the choices of one arbitrary odd-cycle packing per V_i form a subfamily of all packings of T. Their weights multiply and their union gives an injective description of that subfamily. Consequently `h ≥ prod_i H(T[V_i]) ≥ prod_i (1+2ℓ_i) ≥ 1+2j`.

For admissible odd h and n, set `e=min(d,(h-1)/2)`. All coefficients w_j beyond e vanish. If e≥1, the profile count of the exact-h fiber is at most `binom((h-1)/2+e-1,e-1)`, by applying stars and bars to w_j/2; this bound counts profiles, not classes, and may substantially overcount attainable profiles. Either triangular argument then needs only e probes including the already known H. If e=0, the only possible nonempty fiber has h=1, and it has the unique profile P_T(q)=q^n. Even h or h<1 are impossible.

## Optional: one larger strong context already separates all profiles

Let N=n!·2^n. Choose a positive integer k satisfying `(n/(n+1))^k < 1/(100N)`. Such k is O(n² log(n+1)); for example `ceil((n+1) log(100N))+1` works, by log(1−x)≤−x. For this same k, the single response `R_k(T)=H(C3[T,TT_k,TT_k])` distinguishes every pair of length-n path-cover profiles.

To see this, write `R_k(T)=sum_r A_r F_T(r)`, where A_r is the sum of products `B(b)B(c)` weighted by the C3 walk multiplicity at (r,b,c), and `B(j)=F_(TT_k)(j)`. In a transitive tournament there is exactly one directed path on any nonempty vertex subset, so B(j) is the number of onto maps from a k-element set to j labelled bins. In particular
`j^k−j(j−1)^k ≤ B(j) ≤ j^k`.
The summand at (b,c)=(r+1,r+1) has multiplicity one. For r≤n, `(r+1)(r/(r+1))^k < (n+1)/(100N) < 1/2`, giving `A_r ≥ (r+1)^(2k)/4`. For s<r, every allowed b,c≤s+1≤r. There are at most seven pairs b,c (within one of s and with range at most one including s); the central pair has multiplicity three and others have multiplicity one. Hence `A_s ≤ 9r^(2k)`. Since `(r+1)^(2k) > (100N)^2 r^(2k)`, we have `A_r >9N r^(2k)`.

For a tournament T, an ordered r-cover injects into (permutation, r−1 chosen cuts), by concatenating its ordered paths and retaining their boundaries. Thus `sum_r F_T(r) ≤ n! 2^(n−1)` and `sum_r |F_T(r)−F_U(r)| ≤ N`. Let r be the largest index at which two profiles differ. Its nonzero integer difference times A_r has absolute value greater than 9N r^(2k). All lower-index differences combined contribute at most 9N r^(2k), so cancellation is impossible. This is a separation statement, not a claim that such a large oracle query is efficient, nor that its output separates isotypes. The smaller triangular probes are preferable as direct linear reconstruction.

## Prior-search scope

Reviewed the inherited canon statements THM-1975, THM-3121, and THM-4900 at repository head 8f3e642bd29de5bfbb7b0bc1fa558fbfe18b6b90. They establish the raw path-cover substitution, gate-operator, and odd-cycle ingredients. This note establishes the finite recovery by explicit interpolation and the alternative C3 collar estimates; limited textual checks do not rule out overlap elsewhere in the repository or literature.
