# The sphaleron winding group, Steinitz numbers, and the missing completion

**Date:** 2026-10-20 UTC.  
**Status:** CITED physical input; PROVED elementary reconstruction below; FINITE-EXACT checks; no new physical equivalence claimed.  
**Question:** “The SU(2)L gauge-Higgs sphaleron path is isomorphic with Steinitz greatest supernatural number.”

## 1. Verdict and strongest defensible statement

The literal sentence does not yet specify two objects in one category or a map between them. A sphaleron is a field configuration; a sphaleron path or loop is a continuous family of configurations; a supernatural number is a divisibility invariant. They are different types of objects.

There is nevertheless an exact, useful repair:

> An oriented degree-one sphaleron loop generates an infinite cyclic subgroup of the relevant configuration-space fundamental group. Its **profinite completion**, formed by retaining every finite cyclic winding quotient, is isomorphic to the absolute Galois group of a finite field. This completed procyclic group has the greatest supernatural order, Ω = ∏p p^∞.

For characteristic two the maps can be written explicitly:

\[
\langle[\gamma]\rangle\cong\mathbb Z
\ \lhook\joinrel\longrightarrow\ \widehat{\mathbb Z}
\ \xrightarrow{\sim}\ \operatorname{Gal}(\overline{\mathbb F}_2/\mathbb F_2),
\qquad [\gamma]^k\longmapsto k\longmapsto\operatorname{Frob}_2^k.
\]

The first arrow is a dense embedding, not an isomorphism. The second is an isomorphism of topological groups. Ω is the **supernatural order of the target group**, not that group itself and not the winding of a single loop. This bridge applies equally to any chosen infinite cyclic winding group; it is not evidence for a special electroweak–arithmetic dynamical equivalence.

## 2. Precisely which physical topology is available?

Work first with the SU(2) gauge field and a complex Higgs doublet, omitting hypercharge. This is the zero-weak-mixing-angle model. In the radial-gauge construction, finite-energy boundary fields determine maps from the spatial sphere at infinity S² to SU(2) ≅ S³. A loop parameter supplies S¹. The based boundary conditions collapse S¹∨S², giving

\[
U:S^1\wedge S^2\cong S^3\longrightarrow SU(2)\cong S^3.
\]

Its homotopy class lies in π₃(S³) ≅ Z. The standard loop has degree ±1, depending on orientation. These are statements about the **loop's boundary family**, not a nonzero monopole charge of the static sphaleron. [Klinkhamer–Rupp, §§III–IV.A](https://arxiv.org/pdf/hep-th/0304167)

Let B = Map₊(S²,S³), with basepoints fixed. Adjunction gives π₁(B) ≅ π₃(S³) ≅ Z. In a specified gauge/framing model C admitting the boundary map b:C→B, set

\[
w=b_*:\pi_1(C)\longrightarrow\mathbb Z.
\]

For the standard loop, choose orientation so w([γ])=1. Then H=⟨[γ]⟩ is infinite cyclic: if [γ]^k=1, applying w gives k=0. The restriction w|H is an isomorphism onto Z. **This argument needs no assertion that all of π₁(C) is Z.** Manton's mapping-space description supplies the corresponding physical topology. [Manton 2019, §4, equation (5)](https://arxiv.org/pdf/1903.11573)

### Gauge-quotient and existence qualifications

Another familiar description fixes asymptotic data and uses based gauge transformations G₊=Map₊(S³,SU(2)). Its components are π₀(G₊)=π₃(SU(2))=Z. If an ambient field space E is contractible and the chosen completion makes E→E/G₊ a principal gauge bundle, the homotopy exact sequence yields π₁(E/G₊)≅Z. Those are explicit functional-analytic/model assumptions, not a license to identify arbitrary gauge quotients.

Quotienting only by the identity component G₊⁰ retains the large-gauge winding labels. Quotienting by the full based group identifies the endpoint vacua and turns the connecting path into a loop. Changing framing or boundary compactification can change this description. A representative path before the endpoint identification is not automatically a based loop.

The topology motivates a minimax construction; it does not, by itself, supply compactness or a solution of the field equations. Manton's original paper proves the noncontractible-loop input and proposes the unstable solution. [Manton 1983](https://doi.org/10.1103/PhysRevD.28.2019) The subsequent original construction is [Klinkhamer–Manton 1984, *A saddle-point solution in the Weinberg–Salam theory*](https://doi.org/10.1103/PhysRevD.30.2212).

In the standard unwrapped description adjacent vacua differ by ΔN_CS=1, and the sphaleron sits at N_CS=1/2 modulo integers. This coordinate is gauge-dependent up to integer shifts. It is not a prime-divisibility coordinate. [Klinkhamer–Rupp, §VI.A, equations (6.9)–(6.10), figure 11](https://arxiv.org/pdf/hep-th/0304167)

## 3. Supernatural numbers and the exact completion theorem

A supernatural number is a formal product S=∏p p^{e_p}, with e_p∈{0,1,2,…,∞}; infinitely many nonzero exponents are allowed. Divisibility is coordinatewise comparison. Its greatest element is

\[
\Omega=\prod_{p\ \mathrm{prime}}p^\infty.
\]

“Greatest” means greatest in this divisibility order. It is not a largest ordinary integer, real number, ordinal, or cardinal. For a profinite group G, its supernatural order is the least common multiple of the orders of its finite continuous quotients. This is the convention used here. [Ribes, §4, definitions 4.1–4.3](https://www.tau.ac.il/~haran/Teaching/ProfiniteGroups/ribes.pdf)

**Proposition 1 (winding completion; PROVED).** Suppose w:G→Z is a homomorphism and g∈G satisfies w(g)=1. Then H=⟨g⟩≅Z, and

\[
\widehat H=\varprojlim_{n\mid m}H/nH
\cong\varprojlim_{n\mid m}\mathbb Z/n\mathbb Z
=\widehat{\mathbb Z}.
\]

Here the notation nH uses H's additive identification, and the transition maps are reduction modulo n. The supernatural order of the inverse limit is Ω.

**Proof.** The preceding winding calculation shows H≅Z. Its finite-index subgroups are exactly nZ. Each projection from the inverse limit to Z/nZ is surjective. Hence every p^a occurs among its finite quotient orders, for every prime p and every a≥1. The exponent of every prime in its supernatural order is therefore ∞. Conversely every finite quotient order is a finite integer and so divides Ω. ∎

The canonical map j:Z→Ẑ sends k to (k mod n)ₙ. It is injective because an integer divisible by every n is zero. It is dense because any compatible finite list of congruence coordinates is represented by an integer. It is not surjective: Z is countable while Ẑ contains an uncountable p-adic factor. Equivalently, Ẑ≅∏p Z_p by the Chinese remainder theorem.

The infinite exponents were proved from the **chosen family of all finite quotients**. They were not inferred from the mere fact that the winding group is infinite, or that a physical trajectory crosses a barrier.

### The tower is part of the data

For a divisibility-directed family D of positive integers, the restricted completion

\[
G_D=\varprojlim_{n\in D}\mathbb Z/n\mathbb Z
\]

has supernatural order lcm(D). In particular:

- D={2^a:a≥0} gives Z₂ and order 2^∞.
- D={odd positive integers} gives ∏_{p odd}Z_p and order ∏_{p odd}p^∞.
- D={N!:N≥1}, cofinal among all positive integers, gives Ẑ and Ω.

No finite numerical experiment proves cofinality. The definitions and the all-prime argument establish it. If only a bounded or restricted physical measurement family is supplied, replacing it by all n is an additional modeling choice.

## 4. The finite-field realization, with its actual map

Inside a fixed algebraic closure of F₂, put F_{2^n}={a:a^{2^n}=a}. These fields satisfy F_{2^n}⊆F_{2^m} exactly when n|m, and their union is F₂bar. Each finite Galois group is generated by arithmetic Frobenius a↦a² and has order n. Normal bases exist. [Milne, *Fields and Galois Theory*, finite fields, infinite Galois theory, and theorem 5.18](https://www.jmilne.org/math/CourseNotes/FT.pdf)

For a coherent residue family a=(a_n mod n)ₙ∈Ẑ define

\[
\Phi_a(x)=\operatorname{Frob}_2^{a_n}(x),\qquad x\in\mathbb F_{2^n}.
\]

This is well-defined: two field choices embed in their common extension of degree lcm(n,m), and the coherence conditions give the same restrictions. It is an automorphism because each finite restriction is one; its inverse is Φ_{−a}. Every automorphism has a unique such residue family. Composition is Φ_aΦ_b=Φ_{a+b}. Fixing F_{2^n} corresponds exactly to the open kernel a_n=0. Thus Φ is a topological group isomorphism.

The corresponding extension-degree invariant is also Ω:

\[
[\overline{\mathbb F}_2:\mathbb F_2]_{\mathrm{supernatural}}
=\operatorname{lcm}\{[K:\mathbb F_2]:K\subset\overline{\mathbb F}_2,\ [K:\mathbb F_2]<\infty\}
=\operatorname{lcm}_{n\ge1}n=\Omega.
\]

This is a degree invariant, not the cardinality of F₂bar; the latter is countable. More generally S=∏p p^{e_p} specifies the subfield union of F_{2^n} over finite n|S. Closure under lcm ensures this union is a field. Conversely finite subfield degrees of any algebraic extension of F₂ inside F₂bar form a divisor-closed, lcm-closed family, recovered from their supernatural lcm. This gives the finite-field extension interpretation associated with Steinitz. The original bibliographic source is [Steinitz 1910, *Algebraische Theorie der Körper*](https://doi.org/10.1515/crll.1910.137.167); the argument here is self-contained rather than an attempted transcription of that paper.

No property of the sphaleron selects characteristic two. The same construction works over every finite field F_q; binary registers explain why F₂ is a convenient separate realization.

## 5. The useful LFSR bridge and its prime-two hostile

There are two distinct linear actions on a finite field; confusing them breaks the claimed bridge.

### A. Primitive-element multiplication

For a primitive element α∈F_{2^n}*, the linear map x↦αx fixes zero and cycles the 2^n−1 nonzero states. Its order is 2^n−1, always odd. Over all n,

\[
\operatorname{lcm}_{n\ge1}(2^n-1)=\prod_{p\text{ odd}}p^\infty=\Omega_{\mathrm{odd}}.
\]

Proof: 2 never divides any term. For each odd prime power p^a, 2 is a unit modulo p^a, so some positive multiplicative order d satisfies 2^d≡1 mod p^a. Thus p^a divides one of the displayed terms. This proves all odd exponents are unbounded and the exponent at 2 is zero. Primitive maximal-length binary LFSR cycles alone therefore miss Ω.

### B. Frobenius / cyclic-shift action

Choose a normal basis (β,β²,…,β^{2^{n−1}}) of F_{2^n} over F₂. Squaring permutes these n basis vectors in a cycle. Its matrix is the n-position cyclic shift, with characteristic and minimal polynomial x^n−1 and operator order exactly n. In characteristic two, if n=2^a m with m odd,

\[
x^n-1=(x^m-1)^{2^a}.
\]

Thus even Frobenius orders arise through repeated primary factors of the linear operator, even though the finite-field extension itself is separable. The distinction is between the polynomial defining the field and the minimal polynomial of Frobenius as an F₂-linear operator. Since n ranges over all positive integers, this coherent action produces Ω. Normal bases at different degrees need not be compatible; the field-inclusion restrictions supply the coherence before any basis is chosen.

The broader invertible-register formula is also exact. If f=∏ᵢpᵢ^{eᵢ} with pᵢ irreducible and f(0)=1, let rᵢ be the multiplicative order of a root of pᵢ. Then multiplication by x on F₂[x]/(f) has order

\[
\operatorname{lcm}_i\bigl(r_i\,2^{\lceil\log_2 e_i\rceil}\bigr).
\]

Indeed x^{rᵢ}−1 has a simple pᵢ factor because rᵢ is odd; raising it to 2^b multiplies that valuation by 2^b. Chinese remaindering then gives the lcm. Squarefree registers have odd order; unbounded thickness adds 2^∞. Operator order does not imply that every state has that period: zero is always fixed. The companion LFSR investigation develops the state-by-state refinements; this note uses only the operator-clock statement.

## 6. Hostiles that rule out stronger literal readings

1. **Connectedness.** The image of any continuous map from an interval or circle into a profinite group is a singleton. Proof: every finite discrete quotient coordinate is constant, and those coordinates separate points. Consequently a nonconstant sphaleron path cannot be continuously identified with Ẑ or its adding-one motion.

2. **The midpoint.** No additive map from Z[1/2] to Ẑ sends 1 to 1: its image h of 1/2 would satisfy 2h=1, impossible modulo 2. Thus the standard N_CS=1/2 midpoint cannot be inserted into the integer-residue map while preserving this arithmetic. This does not prevent unrelated set-theoretic encodings, which carry none of the requested structure.

3. **Generic universality.** A degree-k loop with k≠0 still generates an abstract copy of Z and still has completed supernatural order Ω. Its image in the ambient completed winding group is kẐ, of index |k|. Ω alone loses this index, the winding orientation, and whether the chosen loop was primitive.

4. **Prime support.** A dyadic tower has order 2^∞; a primitive binary multiplication tower has Ω_odd. Neither equals Ω. “Infinite” does not specify a supernatural number.

5. **Energy and dynamics.** Winding and its finite quotients contain no radial field profiles, Higgs/gauge couplings, energy functional, metric, negative fluctuation mode, physical time, or transition rate. A degree-one abstract circle construction would give the same completion without any sphaleron.

6. **Continuous characters.** The winding character k↦exp(2πiαk) for irrational α has no continuous extension to Ẑ. A continuous homomorphism Ẑ→U(1) has finite image: continuity supplies an open subgroup mapped inside a sufficiently small arc, and no nontrivial subgroup of U(1) stays in that arc. Thus the finite-monodromy completion does not retain every character of the discrete winding group. This is an abstract loss statement; no claim about an observable electroweak theta parameter is needed.

## 7. What precisely is preserved?

| Source | Target and map | Preserved | Lost / additional choice |
|---|---|---|---|
| Framed loop classes in C | w:π₁(C)→Z from boundary degree | Concatenation, orientation-sensitive integer degree | Field profiles, energy, possible kernel of w; gauge and boundary model required |
| H=⟨[γ]⟩ with w([γ])=1 | j:H→Ẑ, [γ]^k↦(k mod n) | Every finite cyclic quotient, exact integer distinction | New completion points; not a continuous encoding of the path parameter |
| Ẑ | a↦Φ_a in Gal(F₂bar/F₂) | Topological group law, all finite restriction maps, distinguished generator | Choice of finite field and arithmetic Frobenius; no physical field interpretation |
| Procyclic group | Its supernatural order | Classification within the procyclic category | Generator, inclusion index, dynamics outside that category |
| Frobenius on F_{2^n} | Normal-basis n-bit cyclic shift | F₂-linear operator conjugacy, period n | Basis is chosen and need not be compatible across n |

A genuinely categorical strengthening is available at the finite-monodromy level: finite covers of the oriented circle and finite étale schemes over F₂ are both described by finite sets equipped with a permutation, equivalently finite continuous Ẑ-sets. Choose the circle generator to correspond to arithmetic Frobenius. Connected n-fold covers correspond to Spec(F_{2^n}); the equivalence preserves the geometric fiber and its monodromy. This is an equivalence of **covering categories**, not a homeomorphism S¹≅Spec(F₂), still less an isomorphism of physical configuration spaces. The finite étale/Galois-action equivalence is documented in [Stacks, lemma 58.2.2, tag 03QR](https://stacks.math.columbia.edu/tag/03QR).

## 8. Universal UHF algebra: same invariant, another category

The factorial matrix system

\[
M_{1!}(\mathbb C)\longrightarrow M_{2!}(\mathbb C)\longrightarrow\cdots,
\qquad a\longmapsto a\otimes I_{n+1},
\]

has norm-completed inductive limit the universal UHF algebra Q, with supernatural invariant Ω and ordered scaled K₀ group (Q,Q₊,1). [Ghasemi–Kubiś, §8.1](https://arxiv.org/pdf/1903.10392)

This does not identify that C*-algebra with Ẑ or with a sphaleron path. Glimm's classification concerns UHF algebras **within their category**. The procyclic-group classification concerns procyclic groups within a different category. A common complete invariant in two separate categories permits a dictionary of isomorphism classes, not a cross-category isomorphism. A physical UHF interpretation would need a specified observable algebra, finite-stage embeddings, involution, norm, and dynamics; none is supplied by winding alone.

## 9. Exact checks, provenance, and stopping boundary

The accompanying `verify_bridge.py` computes the boundary map's degree symbolically, verifies Frobenius/multiplication distinctions over F_{2^n} for 1≤n≤8, verifies normal-basis cyclic shifts, and records finite prime-power divisibility witnesses. Run:

```bash
python3 session9-sphaleron/verify_bridge.py
python3 -O session9-sphaleron/verify_bridge.py
```

The output is deterministic JSON. All checks use explicit exceptions, not assert statements; optimization does not remove validation. The exact degree calculation uses the determinant of (y,∂μy,∂θy,∂φy), producing (1−cosθ)sin²μ sinθ and normalized integral 1. The computed Frobenius cycles are independently checked against the binary-necklace formula and the fixed-point count 2^{gcd(n,k)} for the kth iterate. The finite-field census is a consistency check, not a proof of the infinite completion theorem.

Repository inheritance was read-only at cached commit `5a362f37b195`: no sphaleron, supernatural-number, or LFSR result was found by those exact terms in canon/results. The closest verified neighboring mechanism is `01-canon/theorems/THM-3512-rule30-van-der-put-haar-cocycle-and-profinite-automaton-boundary.md`, whose clock is specifically Z₂ and whose sidecars retain actual orbit information. It is an instructive counterexample to silently replacing a dyadic clock by Ẑ; no Rule 30 statement is imported into sphaleron physics. `THM-3778-rule30-odd-period-finite-scale-cycle-projective-profile-no-go.md` reinforces the existing algebraic-carrier versus physical-realization boundary. General frontier problems remain unchanged.

**Research boundary.** The abstract completion/Galois/finite-monodromy bridge is established. A physical dynamical isomorphism remains undefined and unsupported. To go further, specify an observable or finite-cover family naturally selected by the gauge-Higgs theory, show its maps are gauge-compatible, and supply an intertwining map preserving a requested physical structure. Without that additional input, searching for arbitrary matching prime factors would add no evidence.

**Compact replacement for the original sentence:** “The sphaleron loop's integer winding group has a universal procyclic completion of greatest supernatural order; this completion is isomorphic to Gal(F₂bar/F₂), and Frobenius realizes its finite quotients as binary cyclic-shift clocks.”
