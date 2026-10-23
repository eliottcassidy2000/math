# Session 35: 57, the intermediate Lie algebra, and the information lost by dimension

Status: CITED structural inputs; PROVED elementary deductions below; FINITE-EXACT root census and arithmetic; OPEN interpretations explicitly identified. No claim of new Lie-theoretic or VOA classification. All work is outside the repository; no Git mutation.

## Inheritance and concept board

Closest proved mechanism: forgetting graded/quantum structure retains a finite integer but can destroy reconstruction. Canonical hostile: the session-34 vertex note records distinct Virasoro models with identical C2-algebras. Corrected near miss: THM-4671 says the **adjusted** boundary syzygy lift vanishes, not the unadjusted vector. Least-used sidecar: the first quantum-transfer coefficient in THM-4671 and THM-4695.

Read-only recovery used AGENTS.md, bounded startup, START-HERE, CURRENT-FRONTIER, ACTIVE-GUARDRAILS, RESEARCH-PROTOCOL, CORE-PAPERS, and targeted Git object reads. Relevant recovered paths:

- `01-canon/theorems/THM-4671-virasoro-minimal-model-first-quantum-correction-to-the-arc-space.md`.
- `05-knowledge/results/session34_gauge_representations_20261023/vertex/vertex_simplices_associator_bridge.md` (which gives the corrected inherited scope and pointers).

Live board: (1) arithmetic57; (2) contact grading and Heisenberg bracket; (3) invariant bilinear form; (4) VOA rationality; (5) dimension-forgetting and first quantum transfer. The anchor is an exact E7½ obstruction, niche is the N=57 Virasoro collision, wildcard is whether Blum/Leyland data can select a bracket. No actual LRC implication is claimed.

## 1. What 57 really measures here

CITED: Landsberg–Manivel construct the intermediate algebra from a five-step grading of a simple algebra. For E8, its semisimple degree-zero part is E7 and its positive part is a Heisenberg radical. Thus

    g = e7 ⋉ (V56 ⊕ Cz),   [u,v] = ω(u,v)z,
    dim(g) = 133 + 56 + 1 = 190.

Here V56 is the symplectic fundamental module. We write H57 for its **total dimension**; literature sometimes writes H56 or H28 using other conventions. This is not a new simple Lie algebra of fractional rank. Reference: [Landsberg–Manivel, The sextonions and E7½](https://people.tamu.edu/~jml/LMsexpub.pdf), introduction and intermediate-algebra construction.

PROVED arithmetic: 57=3·19, both primes 3 modulo4, hence 57 is a Blum integer; also 57=2^5+5^2 is a Leyland number. These facts describe the integer, not a construction of the bracket.

FINITE-EXACT independent root check: generate the 240 E8 roots as ±e_i±e_j and half-sign vectors with even minus parity. For θ=e1+e2, grades (α,θ) have counts 1,56,126,56,1 at grades −2,−1,0,1,2. The grade-one roots form 28 pairs {α,θ−α}. Whenever two grade-one roots add to a root, that sum is θ. The standard root-space bracket therefore supplies the Heisenberg pairing (root-space nonvanishing is standard structure theory, not inferred merely from a count).

The zero-grade Cartan has dimension8, so the parabolic has dimension191. Removing its grading direction leaves E7½ of dimension190. This prevents a frequent one-dimensional bookkeeping error. The computation uses an integral Chevalley root datum, so the same dimensions work for the split real realization e7(7)⋉H57(R) inside split e8(8). Complexification yields the complex intermediate algebra. “Split real,” “semidirect split extension,” and “simple Lie algebra” are different predicates. The Heisenberg central extension itself is nonsplit.

## 2. Stronger theorem: the 57-dimensional radical is invisible to every invariant symmetric form

PROVED over C (and by complexification for the specified split real form): every invariant symmetric bilinear form B on g=e7⋉H57 is pulled back from a scalar multiple of the Killing form on e7. Consequently a nonzero such B has radical exactly H57, and g admits no nondegenerate invariant symmetric form.

Proof. Write s=e7 and H=V⊕Cz. The module V is irreducible and nontrivial, so [s,V]=V. Also [s,s]=s and [V,V]=Cz. Thus [g,g]=g, whereas z is central. Invariance gives B(z,[a,b])=B([z,a],b)=0, so B(z,g)=0. The pairing between V and s is an s-equivariant map V→s*, which is zero because V and the adjoint module are nonisomorphic irreducibles (dimensions56 and133). Finally, for u,v∈V and a∈s,

    B(u,[v,a]) = B([u,v],a) = ω(u,v)B(z,a) = 0.

Since [V,s]=V, B(V,V)=0. Hence B(H,g)=0. Invariant symmetric forms on the simple algebra s are scalar multiples of its Killing form, and all their pullbacks are invariant on g. QED.

Immediate hostile: dimension190 alone also allows an abelian Lie algebra, which has many nondegenerate invariant forms. Same dimension, opposite quadratic behavior. Even recording 57 as a radical dimension is insufficient without its bracket and Levi-module action.

Consequence: the ordinary inverse-metric Sugawara construction cannot simply be applied to g as if it were a simple metric Lie algebra. This does not prohibit every possible vertex algebra or conformal construction associated with the intermediate object.

## 3. Rational number is not rational VOA

CITED: Dong–Mason, Theorem1, proves that V1 is reductive for a strongly rational VOA. Here “strongly rational” includes their strong-CFT-type, rationality and C2-cofiniteness hypotheses. [Primary paper](https://arxiv.org/pdf/math/0201318).

NEW APPLICATION (elementary corollary, no novelty claim): there is no strongly rational VOA whose full weight-one Lie algebra is isomorphic to this g, because its Heisenberg radical is nonabelian. This is a condition on the full weight-one Lie algebra, not on the mere presence of an embedded nonreductive subalgebra.

CITED with status discipline: Lee–Sun–Wang associate formal dimension190 and c=38/5 to the level-one intermediate series, and explicitly state their positive-level rationality range as a **conjecture**. Their paper uses intermediate constructions and proposed extensions of familiar affine formulas. The identity 190/(24+1)=38/5 is exact arithmetic, not by itself a valid ordinary Sugawara derivation for g. [Primary paper, abstract and Conjecture1.1](https://arxiv.org/html/2306.09230v1).

CITED: Kawasetsu constructs intermediate vertex subalgebras between E7 and E8 lattice vertex algebras and proves modular differential equations for suitable module characters. [Primary paper](https://arxiv.org/abs/1305.6463). Character modularity is not a replacement for all the strongly-rational hypotheses above. The obstruction does not refute a differently defined “rational intermediate” object; it identifies hypotheses any such interpretation must check rather than silently assume.

## 4. A genuinely sharp dimension57/Virasoro bridge

Apply inherited THM-4671, with its stated minimal-model inputs, rather than guessing from numerology. For Vir_(p,q), coprime 2≤p<q, the C2-algebra is C[x]/(x^N), N=(p−1)(q−1)/2, and c=1−6(p−q)^2/(pq).

PROVED finite classification at N=57: (p−1)(q−1)=114, so the factor pairs give exactly

    (p,q) = (2,115), (3,58), (4,39), (7,20).

All pairs are coprime. All four share the same **graded C2-algebra** C[x]/(x^57), with x of conformal weight2, not merely the same dimension. Their central charges are, respectively,

    −38192/115, −2996/29, −1199/26, −437/70.

The inherited first quantum-transfer coefficient is

    A = 57·115·(c−c_(2,115))/8.

Its values are 0, 10873149/58, 48741099/208, 29904309/112. Thus the three nonboundary models have a one-dimensional first arc-to-associated-graded kernel at conformal weight117, while the boundary model has no such first defect (the all-weight classical-freeness input is separately cited in the inherited canon).

This is an actual source→target map with a tested failure boundary:

- Source: Virasoro minimal-model vertex algebras.
- Target: their C2-algebra, then its integer dimension57.
- Preserved: N and the degree-zero nilpotent relation x^57=0.
- Lost: central charge and quantum syzygy-transfer coefficient.
- Sidecar: c, or equivalently the first transfer amplitude A at fixed N.
- Cheapest decisive test: factor114, evaluate four rational central charges and A.

The four models are not proved to have pairwise different singular supports by the first coefficient alone: the three nonzero coefficients generate the same first kernel line after rescaling. The decisive boundary distinction is the zero/nonzero split. Their differing central charges nevertheless distinguish them as conformal VOAs.

## 5. Where a new hypothesis could genuinely live

OPEN: construct a functorial assignment from arithmetic data (not just the integer57) to a symplectic module with a specified Levi action, and ask whether Blum or Leyland structure constrains that assignment. No such canonical assignment is established here. Every odd dimension2n+1 already supports a Heisenberg Lie algebra, including non-Blum/non-Leyland dimensions, so arithmetic57 alone cannot explain the Heisenberg structure.

The surviving research principle is rigorous: dimension is an integer-valued decategorification; it preserves direct-sum addition and vector-space tensor-product multiplication but discards brackets, real form, invariant metric, and filtered quantum corrections. It is useful as an arithmetic filter, not an inverse classification map. No canonical Lie-algebra tensor multiplication follows merely from vector-space dimension multiplication.

Reproduce: `python3 session35-lie/exact_checks.py`. Standard library only. Output: `exact_checks.json`. Universe: all240 E8 roots; all factorizations114; positive controls θ-complement pairing and arithmetic identities; hostile controls parabolic191 versus intermediate190 and four same-C2-algebra models. Finite checks do not substitute for the cited representation-theoretic or inherited VOA theorems.
