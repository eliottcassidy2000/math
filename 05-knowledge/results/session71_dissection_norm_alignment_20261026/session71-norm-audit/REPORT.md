# Session 71 independent norm / geometry / boundary-rank audit

Status: **PROVED elementary interface results + CITED dissection minimum + FINITE-EXACT controls**. No novelty or formal-proof claim. No shared-repository edits or commits.

Audit basis: repository commit `c2039c32f9726b09ed0e09d5f5c126fc72975897`; source hashes and precise paths in `SOURCE_PINS.json`. The assignment's user transcript is treated as trusted history, not as a mathematical proof. “44 and 5” as historical references remains unresolved in this lane; the integers below are explicitly hostile controls, not an asserted identification of that history.

## 1. Findings first

1. Fermat's theorem genuinely governs rational/integer-coordinate squares and Gaussian squaring genuinely governs Pythagorean triples. Neither scalar norm membership nor area determines a polygonal dissection's minimum number of pieces.
2. The Dudeney connection is an exact **labelled linear compatibility gate** for each proposed matching template. A rank-only summary that forgets the positions of relations fails.
3. Root/Penrose's common-coordinate obstruction is correct, and its lower bound can be strengthened: **every number field** containing both equal-area target polygons has degree at least eight, even without a Galois hypothesis. If degree eight is Galois, its group is necessarily dihedral of order eight.
4. None of the arithmetic or field obstructions below proves the four-piece lower bound. They concern coordinate fields or a proposed template's linear relaxation. The cited geometric proof retains the missing matching, angle, topology, and no-flip constraints.

## 2. Inheritance and scope

The closest proved repository map is **THM-3333**, `01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md`: for `u=(m,n)`,

`Phi(u)=(m²−n²,2mn,m²+n²)`, with `⟨Phi(u),Phi(v)⟩_L=2 det(u,v)²`.

It gives an actual graph isomorphism from primitive projective lattice rays with determinant-one adjacency to the **canonical raw image** with Lorentz-pairing-two adjacency. It does not identify arbitrary integral null triples, normalized triples, or dissections. Its explicit hostile is `X=(4,3,5), Y=(6,8,10)`: pairing two, but X has no integral raw spinor because `C+A=9` is not twice a square.

The least-used relevant sidecar is **THM-3336**, `01-canon/theorems/THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md`: multiplication by `s=a+ib`, followed by primitive division, requires `d_s(u)=gcd(G_su)`. In particular the determinant transforms by `N(s)/(d_s(u)d_s(v))`, not just by the scalar norm. That theorem explicitly separates raw similarity, projective action, and primitive normalization.

Concept board: Gaussian norms / canonical Pythagorean image / labelled matching-length systems / golden-length modules / simultaneous coordinate fields. The operative methods are “type every analogy,” “compute the repair quotient,” and “audit saturation and basis covariance,” from `META-PATTERNS.md`.

The external paper is Demaine–Kamata–Uehara, **arXiv:2412.03865v3**. Its Theorem 1 rules out at most three polygonal pieces **without flips**; the four-piece construction establishes optimality in that category. Its §4 chooses triangle side `τ=2`, square side `σ=3^(1/4)`, equal area `√3`. Observation 3 and Lemma 9 yield equality along edge paths and alternating length relations. In its worked case, the contradiction needs merely `σ≠τ`, not Fermat's theorem or the full irrationality of `σ/τ`. The paper leaves flips and curved pieces outside the established theorem. Source-version observation: arXiv's v3 landing page points to v4 dated 2025-08-01, so v3 must not be called latest.

## 3. Norms: the real bridge and its exact boundary

For integer `n>0`, the classical two-square criterion says

`n=N(a+ib)=a²+b² (a,b∈Z)` iff every prime `p≡3 mod4` has even valuation in n.

The identical criterion characterizes `n=x²+y²` for rational x,y: clearing a common denominator multiplies n by a square, which does not alter those valuation parities; the converse follows from the integer representation. More generally it characterizes positive rational Gaussian norms using all integer-valued valuations, including negative ones.

Consequently `5=1²+2²` is a norm, whereas `44=2²·11` is not. This is not a statement about the square shape in real geometry. **A square of integer area n can have all its vertices in Q² iff n is a Gaussian integer norm.** Indeed a rational edge vector `(x,y)` gives area `x²+y²`; conversely an integer representation gives the four vertices `(0,0),(a,b),(a−b,a+b),(−b,a)`.

No nondegenerate equilateral triangle can have all vertices in Q². Its squared side length q is positive rational, and its shoelace area is rational, whereas its area is `q√3/4`, a contradiction. Thus a proposed rational-Cartesian modelling of both Dudeney targets already fails at the target vertices. Arbitrary real rotations and cuts are permitted by the cited dissection problem.

The norm map is multiplicative but not injective; at norm65, `(1,8)` and `(4,7)` are neither unit associates nor conjugate-unit associates. Squaring retains additional directional information:

`(m+in)²=(m²−n²)+2mn i`, `N(m+in)=m²+n²`.

This gives a right triangle after suitable signs/nondegeneracy conditions. Primitive positive triples require `m>n>0`, `gcd(m,n)=1`, and opposite parity. Two-square representability alone does not imply a nondegenerate integral right triangle with that integer as hypotenuse: n=2 is a minimal counterexample. Nor should n=44 being a nonnorm be confused with an exclusion of rational right triangles of hypotenuse44: `(132/5,176/5,44)` is one.

### Decisive minimum-piece hostiles

Let `k(P,Q)` be the minimum in the no-flip polygonal category. A common similarity of positive scale λ transports dissections in both directions, so `k(λP,λQ)=k(P,Q)`.

For **every real A>0**, an equal-area equilateral-triangle/square pair of area A has k=4 by similarity and the cited theorem. Taking A=5 and A=44 changes integer norm membership but not k. Conversely a square paired with itself has k=1 at either area. Hence area, norm-membership of area, and their combination do not determine k. The first failed implication is forgetting target shape and admissible piece-matching structure.

## 4. Exact labelled boundary-rank theorem

**Theorem (linear-template compatibility).** Let

`A∈Q^(q×n), B∈Q^(q×m), L=Cβ∈R^m`,

where `C∈Q^(m×r)` and the entries of `β∈R^r` are linearly independent over Q. Then

`∃x∈R^n: Ax=BL` iff `rank_Q[A | BC]=rank_Q A`.

Equivalently, for every `u∈ker_Q(A^T)`, `u^TBC=0`. In the consistent case a solution can be chosen with every x-coordinate in `span_Q{β_1,…,β_r}`.

**Proof.** A rational basis of the left kernel of A is also a real basis of its real left kernel. If Ax=BCβ, multiplying by such u gives `(u^TBC)β=0`. Rational independence of β forces every coefficient of `u^TBC` to vanish. Thus every column of BC lies in the rational column space of A, exactly the displayed rank condition. Conversely solve `AY=BC` over Q and take `x=Yβ`. This also proves the final statement. ∎

**Boundary-relation formulation.** Define

`W={B^T u : u∈ker_Q(A^T)}⊂Q^m`,
`Ann_Q(L)={v∈Q^m : v^T L=0}=ker_Q(C^T)`.

Compatibility is exactly `W⊂Ann_Q(L)`. If r is the Q-span dimension of the boundary lengths, `dim Ann_Q(L)=m−r`; the dimensional bound `dim W≤m−r` is necessary but not sufficient. The inclusion, with labels, is the actual predicate.

**Interpretation.** Once a proposed dissection matching template supplies equality/segmentation equations, x can represent its unknown segment lengths and L its named target-side lengths. This theorem eliminates x exactly for the **unconstrained linear equalities only**. It does not generate or enumerate valid templates. It does not impose positivity, angles, closure, embedded nonoverlapping polygons, cyclic orders, orientation, or a piece budget. These remain sidecars. There is no global scissors-congruence obstruction from boundary Q-rank alone.

For Dudeney, `(τ,σ)=(2,3^(1/4))` are Q-linearly independent: `X⁴−3` is Eisenstein at3. Therefore a forced rational relation `aτ+bσ=0` requires a=b=0. But some paper cases need only inequalities or `σ≠τ`; do not attribute every case to irrationality.

For golden lengths `(1,φ,φ²)`, with `φ²=φ+1`, use

`C=[[1,0],[0,1],[1,1]], β=(1,φ)`.

The exact annihilator is generated by `(-1,-1,1)`. Replacing it by the rank-one space generated by `(1,-1,0)` preserves dimension but destroys feasibility. This is a labelled-module bridge to Penrose-type golden length arithmetic, **not an isomorphism of tilings, matching diagrams, or dissection minimizers**. A separate hostile is `x=−1`: its equality-system rank test passes but positive segment length is impossible.

## 5. Independent audit of the common-coordinate field lemma

Credit: proposed by root/Penrose lane; this section independently verifies it and strengthens the degree scope.

**Claim.** Suppose two nondegenerate polygons, one equilateral triangle and one square of the same area, have all their complex vertex coordinates in a number field `K⊂C`.

1. K contains `Q(i,√3)`. A square's successive nonzero edge quotient is ±i. An equilateral triangle's two incident-edge quotient is `(1±i√3)/2`. These two quotients give i and √3, independent of placement, scale, and orientation.
2. No such K can be a finite abelian extension of Q. More strongly, if K is Galois, complex conjugation cannot be central in its Galois group.
3. `[K:Q]≥8`, without assuming K Galois.
4. The bound8 is attained by the target vertices in `K0=Q(i,α)`, `α=3^(1/4)`. K0 is the splitting field of `X⁴−3`, with dihedral group of order8.
5. Every degree-eight **Galois** such K has that abstract Galois group.

**Proof of 2.** Let u be a triangle side and z a square side. Since K is Galois, complex conjugation c preserves it. Set w=z/u. Area equality gives `w c(w)=√3/4`. If c were central, then for every `σ∈Gal(K/Q)`,

`σ(w c(w))=σ(w)c(σ(w))=|σ(w)|²≥0`.

But some σ sends √3 to −√3, because K is normal and contains this quadratic field. Contradiction. Abelian groups have central c. ∎

**Proof of 3.** By1, degree is a multiple of4. If degree4, K equals `Q(i,√3)`, which is an abelian Galois field, excluded by2. Thus degree is at least8. This proof uses no normality assumption on the initial K. ∎

**Proof of 4.** The triangle vertices `0,2,1+i√3` and square vertices `0,α,α+iα,iα` lie in K0 and both areas equal √3. The polynomial `X⁴−3` is Eisenstein, Q(α) has degree4 and is real, so adjoining i doubles degree. All four roots `±α,±iα` lie in K0. The automorphisms `r:α↦iα,i↦i` and `c:α↦α,i↦−i` satisfy `r⁴=c²=1`, `crc=r⁻¹`; hence Gal(K0/Q)=D4 of order8. ∎

**Proof of 5.** An order-eight abelian group is excluded by2. The other possibilities are D4 and Q8. Quaternion Q8 has only one involution, necessarily central; complex conjugation is a nontrivial involution because i∈K, so Q8 is excluded. ∎

The noncentrality is visible at `w=α/2`: `r(wc(w))=−√3/4`, whereas `r(w)c(r(w))=+√3/4`. These differ precisely because r and c do not commute.

### Real Cartesian coordinates are a different field convention

If instead every real x- and y-coordinate lies in a real number field `F⊂R`, the sharp degree is **four**, not eight. For triangle edge `(a,b)` and square edge `(c,d)`, area equality gives

`√3=4(c²+d²)/(a²+b²)∈F`.

Thus `[F:Q]` is even. Degree two would force `F=Q(√3)`, and applying its negative embedding to the displayed ratio is impossible: numerator and denominator stay positive sums of squares. More generally the same proof excludes every totally real F, since that negative embedding extends to a real embedding of F. The field `F0=Q(α)` realizes the displayed target coordinates and has degree four; it is real but not totally real. Its associated **complex** coordinate field is `F0(i)=Q(α,i)`, of degree eight. An arbitrary minimal F need not have Galois F(i); the D4 classification above is explicitly conditional on K being Galois. This comparison agrees with the producer's Cartesian `K²` theorem, not a contradiction to it.

**Do not overread:** this theorem constructs only the two targets in K0. It does not assert that Dudeney's cut vertices lie in K0, does not bound algebraic degrees of cut coordinates, and is not a minimum-piece proof. The D4 here is a field automorphism group, not a physical symmetry group or one of the repository's unrelated Keller D4 obstructions.

## 6. Verification and residuals

Run:

```sh
python3 check_norm_rank.py > check-normal.json
python3 -O check_norm_rank.py > check-optimized.json
cmp check-normal.json check-optimized.json
```

Both runs PASS with **38,815 checks**, byte-identical, using SymPy1.14.0. No load-bearing Python `assert` is used. Seven source snapshots are hash-verified. The exhaustive matrix universe consists of A of size2×1 and B of size2×m, all entries in `{−1,0,1}`, for the two specified C matrices:729 Dudeney-type systems (73 equality-feasible),6561 golden-type systems (513 equality-feasible). Those counts are abstract linear systems, not counts of dissections or tilings. Fraction elimination is cross-checked against direct left-kernel formulas on every system and against SymPy on a specified sparse subuniverse. Integer norm membership is checked for1≤n≤1000; norm and Lorentz identities on coordinate box[−5,5]; parity content on primitive parameters in[−20,20]. The all-order statements rest on the proofs above, not these bounded controls.

Completed: independent mathematical audit, exact controls, field-scope sharpening, source pins, and producer norm-lane cross-review (`CROSS_REVIEW.md`). No task-essential action was denied. Ancillary blocked lookups: the repository's route `CORE-PAPERS-PYTHAGOREAN.md` is absent at the pinned commit; AFP theory-page browser opens failed, so only its successfully read entry page is cited and no AFP build is claimed. Processor message peek cannot run because this assigned read-only lane has no `.machine-id`; no identity was invented. No shared checkout was modified.
