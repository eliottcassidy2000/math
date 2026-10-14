# Euclid avoidance, rigid unit-distance hinges, and boundary-stable critical growth

**Date:** 2026-10-14. **Status:** PROVED constructions and conditional algebraic theorem below; FINITE-EXACT checks; two explicitly OPEN hypotheses. **Novelty unconfirmed.** The classical spindle, valuation-colouring method, Euclid–Mullin graph, and the inherited clutch resultant are credited rather than presented as new. Neither Hadwiger–Nelson nor planar JC is resolved.

## 1. What the connection actually does

Two operations transfer Euclid's multiplicative avoidance into geometry.

1. A rigid triangular-lattice rod of length `n sqrt(3)` can be hinged into an exactly four-chromatic induced unit-distance graph. Every injective realization forces `sqrt(12n^2-1)` into its real coordinate field. Making `n` divisible by the old odd prime factors forces genuinely fresh ramification. The odd-exponent part, not just prime support, is decisive.
2. In the repository's clutch construction, `B_N=1+rad(V)^N` avoids every old boundary component and preserves arbitrarily many boundary jets. Nevertheless its interior critical resultant grows with slope `5 deg(rad(V))`. This excludes these high-exponent candidates as Keller components. Avoidance creates new interior obstructions instead of removing the old obstruction budget.

A third calculation prevents a false bridge: boundary-avoiding cyclic covers merely multiply the escape defect of an étale curve map by the covering degree. New ramification alone supplies no bounded-complexity contradiction.

### Inheritance, board, and scope

The live board is: multiplicative avoidance; prime-exponent parity; rigid port relations; quadratic local fields; boundary jets; critical resultants; escape multiplicity.

- Closest proved geometry mechanism: the classical Moser spindle; the explicit `n=1` construction below is exactly that mechanism. The earlier `planar-recovery` six-colour list-core bound remains a separate, fixed-port result. No claim here upgrades that bound.
- Closest algebraic mechanism: [THM-3225, affine Jacobian clutch resultant and two-boundary no-escape](../../01-canon/theorems/THM-3225-affine-jacobian-clutch-resultant-and-two-boundary-no-escape.md). Its universal resultant is reused; its old affine proof is not re-audited. The new operation uses unbounded-degree clutches with frozen boundary jets.
- Canonical hostile: a fresh prime divisor need not have odd valuation (`3+1=4` gives no quadratic extension). A flexible chain does not force its intended endpoint distance. A high-chromatic residue quotient does not lift a high-chromatic real graph.
- Corrected near miss: fresh quadratic extensions can grow without bound while the whole coordinate-field plane stays four-colourable. The local certificate below proves this, rather than merely warning about it.
- Least-used sidecar: injective triangle propagation, and the parity of the local squareclass at 2.
- No tournament is used: neither unit distance nor divisor coprimality supplies an intrinsic antisymmetric binary relation.

## 2. Literature boundary: which Mullin question?

For the first Euclid–Mullin sequence, `p_(r+1)` is the **least** prime divisor of `1+P_r`, where `P_r=p_1...p_r` and `P_0=1`. This gives `2,3,7,43,13,53,5,6221671,38709183810571,...`. The assertion that every prime eventually occurs remains open in Booker–Simon's January 2026 account. The second sequence selects the largest prime divisor and has different known omission results.

Booker–Irvine's graph has vertices `m` and edges `m -> mp` whenever `p | m+1`; the product forgets the path ordering and the graph is not a tree. See [Booker–Irvine, The Euclid–Mullin graph](https://arxiv.org/abs/1508.03039), Definition preceding Theorem 1.1 and Lemma 2.2.

Polynomial and cyclotomic Euclid–Mullin generalizations are already established territory: [Booker–Simon, A generalisation of the Euclid–Mullin sequences](https://arxiv.org/abs/2601.21901), especially Definition 1.1 and Proposition 1.5. The constructions here do not assert a new general definition of an Euclid sequence. Their proposed contribution is an explicit carrier into rigidity/colouring and into the inherited JC resultant.

## 3. A rigid prime-producing unit-distance hinge

Put `e_1=(1,0)`, `e_2=(1/2,sqrt(3)/2)` and

`P_n={i e_1+j e_2: 0<=i,j<=n}`, for `n>=1`.

Take every unit edge of this patch. Its endpoints are `o=0` and `t=n(e_1+e_2)`, with squared separation `3n^2`.

### Lemma 3.1: rigidity and a real colour relation

The patch has `(n+1)^2` vertices and `3n^2+2n` edges. Every injective planar unit-distance realization of this abstract graph is congruent to the displayed patch. In every proper three-colouring its endpoints have the same colour.

**Proof.** The unit triangles form an edge-connected triangulation of the rhombus. Once one triangle is fixed, the third vertex of an adjacent triangle must be the other equilateral choice across the shared edge: using the same choice would identify two distinct vertices. Triangle propagation therefore determines all vertices. It similarly determines all colours from the first triangle. The explicit colouring is `i-j mod 3`, which assigns the same colour to `(0,0)` and `(n,n)`. Counting horizontal, vertical and diagonal unit edges gives `2n(n+1)+n^2`. ∎

Injectivity is load-bearing. A chain of diamonds sharing only endpoints is a flexible equality gadget, so its drawn length cannot support this field-obstruction theorem.

Set

`Delta_n=12n^2-1`,

`c_n=1-1/(6n^2)`, `s_n=sqrt(Delta_n)/(6n^2)`.

These satisfy `c_n^2+s_n^2=1`. Let `R_n` be the corresponding rotation. Define `H_n` by the points `P_n union R_n P_n`, identifying their common origin and retaining **every** unit edge.

### Theorem 3.2: exact graph and intrinsic field obstruction

For every `n>=1`:

- the two patches intersect only at their origin;
- the only new edge between their nonzero vertices joins the two tips;
- `|V(H_n)|=2(n+1)^2-1`, `|E(H_n)|=6n^2+4n+1`, and `chi(H_n)=4`;
- if `H_n` has an injective planar unit-distance realization with both coordinates in a real field `K`, then `sqrt(3),sqrt(Delta_n) in K`;
- conversely, the displayed realization lies in `Q(sqrt(3),sqrt(Delta_n))^2`.

**Proof of the seam claim.** `Delta_n` is not a square, since it is `2 mod 3`, and it is not three times a rational square, since its 3-adic valuation is zero. Hence `sqrt(Delta_n)` is not in `Q(sqrt(3))`. For old-field vectors `u,v`,

`|u-R_n v|^2=|u|^2+|v|^2-2c_n(u dot v)+2s_n det(u,v)`

(up to the harmless sign convention for the determinant). A unit cross-distance therefore requires `det(u,v)=0`. The same separation shows that a nonzero point cannot lie in both patches.

Collinear nonzero lattice vectors in the positive rhombus have the form `u=r w`, `v=s w`, where `r,s` are positive integers and `w=a e_1+b e_2` is primitive. Write `q=a^2+ab+b^2>=1`. The cross-distance becomes

`q((r-s)^2+rs/(3n^2))`.

If `r!=s`, this is strictly greater than one. If `r=s`, it equals one only when `|u|^2=3n^2`. The norm `i^2+ij+j^2` has this maximum on `0<=i,j<=n` only at `(n,n)`. Thus only the tips are joined.

Both tips must match the common origin in three colours, contradicting their edge. To four-colour, use `i-j mod 3` on both patches and recolour just the second tip with a fourth colour. This remains proper because the only cross-edge is the tip edge.

For the field obstruction, one unit triangle has determinant `+/-sqrt(3)/2` after subtraction of its first vertex, so `sqrt(3) in K`. Patch rigidity forces both tip vectors `u,v` to have squared norm `3n^2`, and their connecting edge forces `u dot v=3n^2-1/2`. Consequently

`det(u,v)^2=(3n^2)^2-(3n^2-1/2)^2=Delta_n/4`.

Their coordinates lie in `K`, so `sqrt(Delta_n) in K`. No assumed global rotation of the coordinate field is needed. ∎

This theorem concerns the entire abstract graph under injective unit realizations, not just the displayed placement. It does not assert that every four-chromatic graph needs fresh ramification.

## 4. A Euclid tower with unbounded field degree and exactly four colours

Let `n_1=1`, `D_j=12n_j^2-1`, `n_(j+1)=n_j D_j`, and

`K_r=Q(sqrt(3),sqrt(D_1),...,sqrt(D_r))`.

Every later `D_j` is `-1` modulo every earlier `D_i`, so the `D_i` are pairwise coprime. Each is nonsquare (`2 mod 3`). Thus each has an odd-valuation prime absent from all earlier radicands and from 3. The squareclasses `3,D_1,...,D_r` are independent, giving

`[K_r:Q]=2^(r+1)`.

This is the exact role of prime powers: even valuations disappear from squareclasses; at least one fresh odd valuation is required. The first radicands are `11,1451,3057042251`.

Nevertheless, for every `r>=1`,

`chi(K_r^2)=4`.

Here is a self-contained upper-bound proof. All `n_j` are odd, so `D_j=3 mod 8`, and `D_j/3` is a square in `Q_2`. Hence `K_r` embeds in `L=Q_2(sqrt(3))`, whose valuation ring has residue field `F_2`. Change coordinates by

`alpha=x+y/sqrt(3)`, `beta=2y/sqrt(3)`.

Then `x^2+y^2=alpha^2-alpha beta+beta^2`. This binary form is anisotropic modulo the maximal ideal: over `F_2` it is one at every nonzero pair. Therefore every unit vector has integral `alpha,beta`, and its residue pair is nonzero. Partition `L^2` into additive cosets of the integral lattice, choose a representative in each coset, and colour each point by the residue pair after subtracting its representative. This gives four colours; a unit edge stays in one coset and changes the residue pair. Pull back to `K_r^2`. The lower bound is `H_1`.

The valuation/coset method is classical; see [Madore, The Hadwiger–Nelson problem over certain fields](https://arxiv.org/abs/1509.07023), Proposition 3.2 and the quadratic-form discussion. The new application here is the explicit fresh-ramification tower and its rigid graph carriers. No literature priority is claimed for that application.

For **even** `n`, `Delta_n=7 mod 8`. Locally adjoining its square root adds `sqrt(-1)` to `Q_2(sqrt(3))`, and the specific anisotropy certificate above fails. Failure of this certificate is not proof of five-chromaticity. It is a precise arithmetic filter for choosing a different search sector.

## 5. A reusable conditional colour escalator

Let `G` be a finite `k`-colourable unit-distance graph with terminals `a,b` that are equal-coloured in every proper `k`-colouring. Suppose a realization has terminal separation `R>=1/2`. Two copies can be rotated about their common `a` so that their `b` terminals are unit-distance, using

`cos(theta)=1-1/(2R^2)`, `sin(theta)^2=(4R^2-1)/(4R^4)`.

If the combined placement is injective except for the intended common vertex, its complete unit-distance graph is not `k`-colourable. Extra unit edges do not invalidate this lower bound. Rigidity is needed only if one also claims an intrinsic field obstruction or a prescribed distance in all realizations.

For `k=5` or `6`, producing the required equality relation is a genuine new obligation. Deleting an edge from an already `(k+1)`-critical unit graph can yield an abstract equality gadget but recycles the stronger obstruction; its removed geometric unit edge is still present in the full unit-distance graph. That is not an independent escalation mechanism.

## 6. Euclid clutches force new interior critical points

Let the ground field be algebraically closed of characteristic zero. Fix nonconstant `V(x)` and a polynomial `A(x)` satisfying

`2 V A' - A V' = 2 V`.

For any polynomial `B`, set

`P_B(x,z)=(V z^2+B z)^2+A z+x`, `J_B=2 V B'-B V'`.

On `V!=0`, write `y=Vz`, `L=y^2+B y`. The inherited exact gradient reduction is

`R_1=2L(2y+B)+VA`, `R_2=V^3+V^2 y+J_B yL`.

The universal resultant from THM-3225 is `Res_y(R_1,R_2)=V^3 K_B`, where

`K_B = -A^3J^3+12A^2J^2V^2-4AB^3J^2V+4AB^2JV^2`

`       +24ABJV^3-48AJV^4-16AV^4-8B^4JV^2+8B^3JV^3`

`       +32B^2V^4-96BV^5+64V^6`, with `J=J_B`.

Let `g=rad(V)` be monic, `B_N=1+g^N`, and let `K_*` denote `K_B` at `B=1`. Assume `K_*` is nonzero. Put

`v=deg V`, `a=deg A`, `s=deg g`, `M=sum_(V(alpha)=0) ord_alpha(K_*)`.

### Theorem 6.1: boundary-stable linear critical growth

For all sufficiently large positive integers `N`,

`ord_alpha(K_(B_N))=ord_alpha(K_*)` for every root `alpha` of `V`,

`deg K_(B_N)=5Ns+3v+a-2`,

and there are exactly `5Ns+3v+a-2-M` units of resultant-root multiplicity away from `V=0`. In particular, `P_(B_N)` has a finite critical point and cannot have a polynomial constant-nonzero-Jacobian mate.

**Proof.** Since `B_N-1` is divisible by `g^N`, and

`J_(B_N)-J_1=2V(g^N)' - g^N V'`,

this difference also has order at least `N` at every root of `V`. Substitution in the displayed polynomial gives `K_(B_N)-K_*` divisible by `g^N`. For `N>M`, every old boundary multiplicity is therefore unchanged. This is stronger than mere disjointness: arbitrarily many jets are frozen.

Write `b=Ns`. The response identity implies `a>=1`. If `a>1`, its leading terms force `v=2a`. For `b` sufficiently large, the unique highest-degree term of `K_B` is `-4AB^3J^2V`, of degree `5b+3v+a-2`; its leading coefficient is nonzero since `J` has leading coefficient `(2b-v) lc(V) lc(B)`.

If `a=1`, the terms `-4AB^3J^2V` and `-8B^4JV^2` share that degree. The response identity gives `lc(A)=2/(2-v)`, so their combined leading coefficient is

`-4 lc(B)^5 lc(V)^3 (2b-v) (lc(A)(2b-v)+2)`.

It is nonzero for `b>v`. All other terms have lower degree once, for example, `b>2(v+a+1)`. Thus the degree formula holds. Remove the fixed boundary contribution `M`; the remaining polynomial degree is positive for large `N`.

Every resultant root with `V!=0` yields an affine common zero of `R_1,R_2`, because the leading coefficient of the cubic `R_1` in `y` is the constant 4. The inherited gradient equivalence makes it a critical point of `P_B`. ∎

For the actual THM-3225 accessory families, its nonzero constant-clutch resultant supplies the input `K_*!=0`. This theorem therefore rules out an unbounded-degree family of new clutches there. It does not quantify the threshold for those accessory coefficients, and it says nothing about arbitrary polynomial maps or all possible nonlinear repairs.

### Exact toy family

Take `V=x`, `A=2x`, `g=x`:

`P_N=(xz^2+(1+x^N)z)^2+2xz+x`.

The independent exact resultant checks for `N=1,...,8` give `deg K_N=5N+2`, `ord_0 K_N=3`, hence interior multiplicity `5N-1`: `4,9,14,19,24,29,34,39`. For `N>=2` in this range the leading coefficient is `-16N(2N-1)`. These multiplicities are not counts of distinct Morse points. The general theorem, not the finite range alone, supplies the unbounded conclusion.

## 7. Why a fresh cyclic cover does not prove properness

Let `g:C->A^1` be an étale morphism from a smooth connected affine complex curve, extending to a degree-`d` map of its smooth projective completion. Let its finite punctures be `b`, with local degrees `e_b`. Riemann–Hurwitz gives the exact escape mass

`delta(g)=d-chi(C)=sum_(finite punctures b) e_b`.

Choose a squarefree polynomial `H(t)=1+Q(t)^m`, where `Q` vanishes at every finite escape value. Then `H=1` at those values. Its zeros are away from all escape values. Let `B'` be the affine curve `z^ell=H(t)` and `C'` the corresponding pullback of `g`, normalized if needed. With `r=deg H`, the finite-branch Euler characteristic identities are

`chi(B')=ell-(ell-1)r`,

`chi(C')=ell chi(C)-(ell-1)dr`.

Thus `d chi(B')-chi(C')=ell delta(g)`. The old escape punctures split into `ell` copies, carrying the same local degrees. Nothing forces the normalized defect to increase. Growing genus in this family is allowed; a contradiction would require an additional justified bounded-complexity invariant.

For a planar Keller map one may apply this bookkeeping to a smooth generic first-coordinate fibre and its étale second-coordinate map, with generic connectedness/degree and all punctures retained. It is not a proof that an étale map is proper.

## 8. Two new, explicit hypotheses to try next

Both are **OPEN / novelty unconfirmed**. They are search targets, not conclusions.

### H-G: mixed-parity hinge coupling can exceed four colours

**REFUTED later in the same session.** The historical proposal below is superseded by the [global independent-hinge seam theorem and explicit four-colouring](independent_hinge_seam_four_colouring_20261014.md), which proves `chi(G_m)=4` for every `m>=1`. The obstruction survives overlapping seam norms; the next operation uses affine attachment.

Let `R_1` have cosine `5/6`, sine `sqrt(11)/6`, and `R_2` have cosine `23/24`, sine `sqrt(47)/24`. For `m>=1`, let `G_m` be the complete unit-distance graph on the distinct points

`{R_1^a R_2^b (i e_1+j e_2): |a|,|b|,|i|,|j|<=m}`.

**Hypothesis:** some `G_m` is not four-colourable.

This is a specific mechanism claim inside `Q(sqrt(3),sqrt(11),sqrt(47))^2`, rather than a claim that fresh radicals automatically amplify chromatic number. Odd-only Euclid towers above give a proved hostile control; the second rotation deliberately enters the missing local squareclass.

**Prediction and decisive tests:** enumerate every exact edge for `m=1,2,...`, retain an independently checkable four-colouring or an UNSAT certificate, and inspect the first coupled boundary relation that cannot be realized. The candidate size is at most `(2m+1)^4` before point identifications. The cheapest tests are `m=1` and `m=2`; these mixed-rotation tests have not been run in this checkpoint. A compatible four-colouring of the entire union disproves the hypothesis. Four-colouring any one finite stage only rejects that stage. The known individual hinged patches, even those with even `n`, are exactly four-colourable and are not positive evidence by themselves.

### H-J: all-degree unit clutches retain an interior obstruction

For each of the fixed accessory pairs `V,A` in THM-3225, and every polynomial `B` with `gcd(B,V)=1`, the saturated polynomial

`K_B / product_(V(alpha)=0) (x-alpha)^(ord_alpha K_B)`

has positive degree (with `K_B=0` treated separately as an even stronger common-root case).

**Hypothesis:** this holds for every such `B`.

It would exclude every polynomial clutch in that fixed source architecture. The inherited affine theorem and the high-exponent Euclid theorem cover different nontrivial sectors; neither proves the whole statement.

**Prediction and decisive tests:** start with `B=b_0+b_1x+b_2x^2`, then cubic `B`. Solve the coefficient conditions saying the entire resultant is supported on `V=0`, while saturating by `Res(B,V)` so forbidden boundary collisions are not counted as successes. A single admissible `B` with no critical point refutes the hypothesis. Keep high-multiplicity boundary jets and leading-coefficient cancellation in the test; degree growth without those coordinates is insufficient.

## 9. Reproduction and finite scope

Run:

`python3 04-computation/experiments/euclid_rigid_hinges_and_clutch_growth_20261014.py`

`python3 -O 04-computation/experiments/euclid_rigid_hinges_and_clutch_growth_20261014.py`

Both produce the same JSON byte-for-byte. The checks use exact integers and SymPy rational polynomial resultants, no floating-point distance thresholds. They verify all seams and explicit four-colourings for `n=1,...,8`, three fresh odd-tower stages, and an independent resultant for every toy exponent `N=0,...,8`. General rigidity, field necessity, the local four-colouring, and the asymptotic clutch theorem are mathematical proofs above, not inferred from the finite census.

Files: this note, the experiment script, and `euclid_rigid_hinges_and_clutch_growth_20261014.json`. No theorem number is reserved, no old result is recanonized, and no live checkout or remote was modified by the producing worker.
