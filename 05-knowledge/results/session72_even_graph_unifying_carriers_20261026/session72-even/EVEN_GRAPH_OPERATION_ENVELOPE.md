# The tournament/even-graph identity as an affine harmonic correspondence

**Status:** PROVED elementary operation-envelope theorem and switching-gap corollary below; INHERITED / CITED counting, Walsh, codec and obstruction results; FINITE-EXACT implementation controls as specified below. No historical novelty claim. This packet proposes a larger organizing theorem, not a new proof of the classical counting identity and not a natural bijection of labelled objects.

Base: `4950589d056a64f8807774cf7da893dcafcd63a0`. `SOURCE_PINS.json` identifies exact repository source paths, frozen copies and SHA256 values. Missing sparse-worktree files were read with `git show BASE:path`, not treated as absent mathematics.

## 1. Recover the actual identity before generalizing it

There are **three different parities**, not one interchangeable class of even graphs:

1. **Sign-even / automorphism-even / Royle-even.** For a simple graph X on [n],
   `epsilon_X(g)=(-1)^{#{ij in E(X):i<j, g(i)>g(j)}}`, for g in Aut(X).
   X is sign-even iff this character is trivial. The parity is independent of the reference orientation. It is not the sign of the permutation induced on the edge set.
2. **Degree-even / Eulerian.** Every vertex has even degree. These edge sets form the binary cycle space Z_n.
3. **Even number of edges.** This is the character controlling global tournament reversal in the Fourier basis; neither of the first two conditions implies it.

The three-vertex path is sign-even but not degree-even. The triangle is degree-even but sign-odd. A three-leaf star is sign-even with an odd number of edges. The empty graph is sign-even in every order.

The classical identity is

    # {tournaments on n vertices}/isomorphism
      = # {sign-even simple graphs on n vertices}/isomorphism.

It has the SAME n on both sides and counts UNLABELLED types. For example n=3 gives two types on each side but eight labelled tournaments versus four labelled sign-even graphs. At n=4 the labelled counts are 64 versus 20, and both type counts are four. No n-to-n+1 shift is present. The primary paper's Theorem 1.1 and Sections 1–4 verify precisely this scope: Royle, Praeger, Glasby, Freedman, Devillers, [*Tournaments and Even Graphs are Equinumerous*, arXiv:2204.01947v2](https://arxiv.org/html/2204.01947v2), published in *J. Algebraic Combinatorics* 57 (2023), 515–524. Primary text checked 2026-10-26.

For m=binom(n,2), both type counts equal

    (1/n!) sum_{g in S_n, ord(g) odd} 2^{c_2(g)},

where c_2(g) is the number of orbits of g on unordered pairs. For n=0,1 equality is trivial. The sequence begins 1,1,1,2,4,12,56,456,6880. In particular the equality is not the assertion that degree-even graph types count tournaments: at n=4 the degree-even types number three, not four.

### Existing unification, not a new claim here

THM-4531, `twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes`, already proves for every S_n-invariant U in V=F_2^{E(K_n)}:

    #(Tour/U)/S_n = #{sign-even X in U^perp}/S_n.

Tour is an affine V-torsor, and U acts by reversing its edge sets. U=0 is the classical identity; U=Cut gives tournament switching classes versus **sign-even AND Eulerian** graphs. THM-4721, `twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis`, supplies the signed Walsh basis and prime-local refinements. THM-4724, `affine-fixed-point-criterion-restricted-twisted-burnside-for-every-submodule-even-n-switching-half-and-the-dc-obstruction-to-uniform-circulant-rules`, gives the odd-P fixed-point extension for every U.

By contrast THM-4069, `even-graph-basis-dependence-and-canonical-cycle-envelope`, and THM-4078, `even-graph-triangle-quotient-spectrum-and-boolean-noncommutation`, use the **linear** cycle space Z_n. Its Fourier dual is V/Cut with the ordinary signed-graph action, not the affine tournament action. Their dual modes are ordinary signed-graph switching classes, not tournament switching classes. This distinction is the missing bridge coordinate.

## 2. Proposed overarching theorem: affine Fourier descent with operations

Let M be a finite vector space over F_2, let A be a torsor for M, and let a finite group G act affinely on A with linear part rho:G→GL(M). Choose a reference a_0 solely to write

    g(a_0+x)=a_0+rho_g x+c_g,
    c_{gh}=c_g+rho_g c_h.

Let M^vee be the character group, with chi(x) in {+1,-1} and `(g chi)(x)=chi(rho_g^{-1}x)`. For chi let G_chi be its stabilizer and define

    epsilon_chi(g)=chi(c_g),   g in G_chi.

This is a character of G_chi. Its triviality is independent of a_0. Call the G-orbit of chi **admissible** if epsilon_chi is trivial.

### Theorem A — counting, joint spectrum, and exact operation visibility

Over Q (or any characteristic-zero field containing these character values):

**A1.** `Q^A` has its usual Walsh decomposition into character lines. Every G-orbit Omega in M^vee contributes one invariant line exactly when Omega is admissible, and zero otherwise. Consequently

    dim (Q^A)^G = |A/G| = # admissible G-orbits in M^vee.

The line is intrinsic; its displayed generator is determined only up to sign. A chosen reference and one representative in Omega give a signed orbit sum `Phi_Omega`. This is a canonical line decomposition, not a canonical choice of a signed basis vector and not a set bijection.

**A2.** For a **source-independent, linear-G-invariant** kernel k:M→Q, define

    (C_k f)(a)=sum_{d in M} k(d) f(a+d).

This preserves G-invariant functions and

    C_k Phi_Omega = k_hat(chi) Phi_Omega,
    k_hat(chi)=sum_d k(d) chi(d),   chi in Omega.

Thus all such weighted operations commute after orbit descent, with multiplicities and loops retained.

**A3.** As k ranges over all linear-G-invariant rational kernels, their restrictions to `(Q^A)^G` form the **entire diagonal algebra** on the admissible lines. The exact equality criterion is

    C_k|_(Q^A)^G = C_l|_(Q^A)^G
      iff k_hat(chi)=l_hat(chi) for every admissible chi.

In particular the operation kernel is the space of invariant kernels whose Fourier transform vanishes on all admissible dual orbits. Different labelled operation families can therefore become indistinguishable on unlabelled affine objects.

The full algebra includes signed rational kernels; it does not assert that every diagonal projection is itself a nonnegative Markov move. Nonnegative k with sum k=1 does give a Markov operator, and its quotient has stationary mass proportional to labelled orbit size.

### Proof of A

Write `f_chi(a_0+x)=chi(x)`. The pullback action `g f(a)=f(g^{-1}a)` sends f_chi to a signed multiple of f_(g chi); explicitly the multiplier is `(g chi)(c_g)`. On the stabilizer it is epsilon_chi. A nontrivial stabilizer sign forces an invariant coefficient to be its negative; a trivial stabilizer sign permits a freely chosen coefficient to propagate consistently around the orbit. Different dual orbits have disjoint Walsh supports. This proves A1, including completeness and the reference-choice boundary.

Translation by d acts on f_chi by chi(d), so C_k acts by k_hat(chi). Linear G-invariance makes that multiplier constant on Omega and makes C_k commute with the affine action. This proves A2.

Ordinary Fourier inversion identifies the space of linear-G-invariant k with all G-invariant functions on M^vee. For each dual orbit Omega, the explicit kernel

    e_Omega(d)=|M|^{-1} sum_{chi in Omega} chi(d)

has Fourier transform equal to the indicator of Omega. It projects onto the corresponding Walsh orbit space; on affine invariants it is the single admissible-line projection or zero. Hence the available restrictions are all diagonal operators and A3 follows. Uniform measure on A is preserved by each translation; since d=-d, each C_k is self-adjoint for this measure. Orbit descent yields the stated reversible weight. QED.

### Theorem B — coarsening is dual restriction, naturally with every operation

Let N≤M be G-invariant and let q:A→A/N. The pullback q* embeds functions on A/N into functions on A. Its dual support is exactly N^perp. Under the same reference convention:

* admissible orbits in N^perp are unchanged, because stabilizers and the restricted affine cocycle are unchanged there;
* normalized fibre averaging

      (E_N f)(a+N)=|N|^{-1} sum_{u in N} f(a+u)

  fixes the surviving character modes and kills the others;
* for the pushed kernel `k_N(d+N)=sum_{u in N} k(d+u)` — a **sum, not an average** — both diagrams commute:

      C_k q* = q* C_(k_N),
      E_N C_k = C_(k_N) E_N.

These statements restrict to the G-invariant function spaces. For tournament submodules U⊆W they apply to M=V/U and N=W/U: coarsening Tour/U→Tour/W deletes precisely the signed graph modes outside W^perp.

**Proof.** A character is constant on N-cosets iff it annihilates N. The finite character sum over N is |N| in that case and zero otherwise. The two intertwining equations follow by collecting the d-summands by d+N, and interchange of the two finite sums proves the averaging equation. The stabilizer sign is unchanged because characters in N^perp kill the cocycle change modulo N. QED.

### Sharp operation scope

A labelled operator written `(Lf)(a)=sum_d k_a(d)f(a+d)` commutes with **every translation of A** iff k_a(d) is independent of a. This follows by testing point indicators; translation conjugation replaces k_a by k_(a+u). Thus source-independent kernels are an exact class, not an aesthetic preference. A G-equivariant source-dependent operator can still descend to G-orbits, but its Fourier matrix need not be diagonal and it need not commute with coarsening by N. No assertion is made that every source-dependent kernel fails after every quotient.

## 3. Concrete payoff: the switching-class arc-flip gap is exactly eight

Let C_n be the isomorphism types of tournament switching classes on [n]. Retain the operator

    (A_n f)([T])=sum_{e in E(K_n)} f([T reversed at e]),
    L_n=m I-A_n,   m=binom(n,2).

Cuts commute with arc reversal, and relabellings permute the named edges, so this is a well-defined multiplicity-weighted quotient operation. Distinct flips leading to the same class still contribute separately, and a flip staying in the same class contributes a loop.

**Corollary C (PROVED).** The entire Laplacian spectrum is

    Spec(L_n) = { 2|E(X)| : X a sign-even Eulerian graph type on n vertices },

one entry per type, including repeated numeric eigenvalues. For every n≥4,

    gap(L_n)=8, with a one-dimensional eigenspace.

For n≤3 there is just one switching-class type; no positive gap is asserted. For n≥4 the nonlazy **algebraic** stochastic gap is 8/m; for `(I+A_n/m)/2` the gap is 4/m and relaxation time m/4.

**Proof.** Apply A to M=V/Cut and A=Tour/Cut. Its dual is the Eulerian space Z_n. Admissible dual orbits are precisely the sign-even Eulerian graph types, by the inherited cocycle identification. The edge-flip kernel has Fourier multiplier `sum_e (-1)^[e in X]=m-2|E(X)|`, even when several edge increments coincide modulo cuts; their multiplicity is retained. This gives the spectrum. A nonempty degree-even simple graph has at least three edges. The only three-edge possibility is C3 plus isolates; a transposition of triangle vertices reverses an odd number of edges, so it is not sign-even. The only four-edge possibility is C4 plus isolates; orient the cycle cyclically, and its dihedral automorphisms preserve that orientation or reverse all four edges, so it is sign-even. Isolated-vertex permutations contribute no edge reversals. Thus four is the unique minimum edge count among nonempty admissible types. QED.

A concrete eigenfunction is

    F(T)=tr(K(T)^4)-n(n-1)(2n-3),

where K is the **skew** sign matrix of a tournament. Its repeated-vertex closed 4-walks contribute exactly the subtracted constant; every remaining monomial is a four-cycle, and the coefficient of each signed four-cycle is ±8, nonzero. Hence F is a pure nonzero Walsh-degree-four function, switching/relabel invariant, and L_n F=8F. This is not the symmetric Rademacher quadratic-form statement of THM-2501: for skew K, x^T K x is identically zero. The spectral proof does not use such a quadratic form.

**Inheritance boundary.** THM-584/587/588 already prove the ordinary tournament arc-flip spectrum, its signed cycle-index multiplicities, global reversal trace, and gap four. Reinterpreting those multiplicities as edge counts of sign-even graph types is a reconciliation with THM-4721, not a new ordinary-tournament spectrum. The simultaneous operation algebra and coarsening diagram are the enlarged organizing statement here. Corollary C is a further elementary quotient consequence; targeted searches found no earlier exact switching-gap-eight statement in the frozen sources, but no literature-wide priority is claimed.

The two gaps refer to distinct quotients: ordinary tournament cyclicity survives at Walsh degree two, whereas switching erases it because a two-edge path is not Eulerian. The first surviving mode is the signed four-cycle. This is a concrete explanation of why a quotient changes the slowest observable, not a heuristic analogy.

## 4. The actual bijections and the limits of naturality

### A. SCC/Lyndon codec: all orders, type-valued, rank-dependent

`05-knowledge/results/tournament_even_bijection_20261011c.md` provides an effective inverse pair conditional on the inherited all-order count theorem:

    tournament type → ordered SCC-type word → nonincreasing Lyndon factors
       → equal-rank connected sign-even graph atoms → disjoint union.

The inverse decodes the atom ranks, sorts factors in nonincreasing order and concatenates. The alphabet order is `(vertex count, canonical tournament mask)`; connected graph atoms are sorted by canonical edge mask. It preserves vertex count and identifies the multiset of Lyndon-factor weights with the graph's component sizes. This is not a rank-free proof of equinumerosity, not label-equivariant, and not ordinary-join to disjoint-union multiplicativity.

`05-knowledge/results/tournament_converse_context_20261011d.md` transports ordered join to an associative **conditional fusion** product: adjacent Lyndon factors L<M fuse to LM, and every legal fusion order has the same normal form. It transports converse to an anti-involution D. With a the singleton and b=C3,

    Phi(a^r b)=K_(1,r+2),   Phi(ba^r)=P3 disjoint-union r K1,
    D(K_(1,r+2))=P3 disjoint-union r K1.

This supplies a genuine operation-compatible type correspondence, but the target operation is conditional fusion, not ordinary disjoint union. The literal reversible labelled fan operation retains the centre and a two-element leaf subset. The unmarked completion has exactly binom(n-1,2) preimages per star, so that many sidecar values are necessary and sufficient for this completion. An equivariant unmarked leaf-pair choice does not exist. There is no general chromatic-number preservation: the inherited n=6 masks 110 and 123 have chromatic numbers two and three and are interchanged by D.

### B. Prime circulants: explicit local routing, with correctly typed fibres

`05-knowledge/results/circulant_necklace_bijection_20261020/CIRCULANT_NECKLACE_BIJECTION.md` identifies prime-order circulants with negacyclic and cyclic binary necklaces, then uses differentiation, integration and a specified branch cut in the finite derivative functional graph. It is a type bijection for every odd prime p. Its forward realization is an AGL(1,p)-equivariant **surjection** on translation-invariant labelled objects. It is not a labelled bijection: already p=3 has two labelled source circulants and one target graph. The map preserves odd multiplier subgroups; its labelled fibres are powers of two. The primitive-root and off-cycle branch choices remain part of the recipe. The constant-step derivative D_c cannot be universally substituted for this source-aware stopping rule.

### C. What symmetry refinements do and do not say

The repaired inherited statements are:

* Equality of type counts with a specified odd p-group subgroup predicate is proved, separately for each odd prime p (`prime_local_walsh_bijection_20261020/SYLOW_ACTION_ADDENDUM.md`).
* Equality for **every odd-order mixed-prime H** is false: H=C15 acting as a 3-cycle, a disjoint 5-cycle and two fixed vertices at n=10 gives 32 tournament types versus 31 sign-even graph types (`c15_symmetry_counterexample_20261020/C15_ORDER_TEN_COUNTEREXAMPLE.md`).
* For every n≥5 no S_n-equivariant map from labelled tournaments to labelled graphs induces a bijection onto sign-even graph types (THM-4722 and the corrected THM-4531). The reverse-direction obstruction is already the empty graph: it is S_n-fixed, no tournament is when n≥2.
* Equal single-prime marginals do not supply one joint all-prime preserving bijection. The empty graph admits both 3- and 5-symmetries; no tournament can simultaneously contain the relevant Sylow actions.

Therefore the strongest natural common object is currently the signed harmonic decomposition with operations, not an identification of unmarked objects or their automorphism groups.

## 5. Exact controls and decisive hostiles

Run from this directory:

    python3 -B verify_affine_fourier.py > normal.json
    python3 -B -O verify_affine_fourier.py > optimized.json
    cmp normal.json optimized.json

The script has only standard-library imports, no assert statements, no floating-point constants, explicit finite universes, and identical normal/-O outputs. It checks **81,177 gates** across 25 models: n=0,...,5 with U=0, Cut, span(all edges), and V, plus n=6 with U=Cut. Admissibility is computed by transport consistency under adjacent transpositions, independently checked by every automorphism permutation for ordinary graphs through n=5. All orbit-of-difference kernels are checked against every invariant eigenline; exact rational elimination verifies both the signed basis and the full diagonal operation algebra. Coarsening kernels are checked as literal integer counts, not inferred from ranks.

Switching results:

| n | labelled switching classes | unlabelled switching classes | Laplacian eigenvalues with multiplicity |
|---:|---:|---:|---|
| 0,1,2 | 1 | 1 | 0 |
| 3 | 2 | 1 | 0 |
| 4 | 8 | 2 | 0,8 |
| 5 | 64 | 2 | 0,8 |
| 6 | 1024 | 6 | 0,8,12,16,16,24 |

The n=6 weighted matrix and representatives are in `normal.json`. Finitely many controls are not the proof of Corollary C; the minimum-support argument is the all-order proof.

Hostile controls, with exact receipts:

1. **Omit the affine cocycle:** n=4 yields 11 ordinary graph orbits instead of four tournament orbits.
2. **Call sign-even a cycle space:** P3 is sign-even, C3 is not; sign-even sets are not closed under symmetric difference. Admissible graph types label Fourier lines; they do not inherit the addition of V.
3. **Assume distinct kernels survive:** on n=3, one-edge and two-edge reversal kernels are distinct labelled difference orbits but both descend to `[[2,1],[3,0]]`. The admissible characters all have even support, so global reversal is invisible.
4. **Erase loops/multiplicity:** at n=4, the commuting one-edge and two-edge Hamming-layer operators have loop-deleted Boolean supports B1,B2 with `(B1 B2)[0,1]=2`, `(B2 B1)[0,1]=1` in the stored type order `[0,2,4,5]`. Thus Booleanization leaves the commuting algebra.
5. **Remove the source-independent hypothesis:** labelled directed-triangle reversal fixes state 0, but sends cyclic state 2 to 5. Translation by 2 does not commute with it. This example is explicitly labelled; after the two-type quotient the same map happens to be identity, so it is not falsely promoted to an orbit-level failure.
6. **Mistake a positive gap for an absolute nonlazy mixing rate:** at n=4 switching the adjacency eigenvalues are 6 and -2. The algebraic stochastic gap is 4/3, whereas the lazy gap is 2/3. Laziness is part of the stated relaxation conclusion.

## 6. Scope ledger and useful larger questions

| Claim | Status | Preserved / lost |
|---|---|---|
| Same-n unlabelled tournament/sign-even graph count | CITED / INHERITED PROVED | Cardinality, not labels or automorphism group |
| All-U twisted count and signed Walsh basis | INHERITED PROVED | Invariant-function space, with stabilizer phase retained |
| Theorems A–B, operation visibility and coarsening | PROVED elementary here | All invariant source-independent weighted operations; no Booleanization claim |
| Switching-class full arc spectrum and gap eight | PROVED corollary here; FINITE-EXACT n≤6 | Correct multiplicities, loops, graph parity and quotient |
| SCC/Lyndon inverse pair and conditional fusion | INHERITED PROVED relative to count theorem | Weight and decoded word; ranking choices retained |
| Prime derivative-path map | INHERITED PROVED, scoped type bijection | Odd multiplier predicates; labelled fibres not discarded |
| Universal mixed-prime symmetry matching | REFUTED | C15 n=10, 32≠31 |
| One natural S_n-equivariant type map, n≥5 | REFUTED in the precise inherited scope | Arbitrary type bijections not excluded |
| A rank-free, operation-selected atom pairing for all n | OPEN | Needs a declared operation and equal refined fibres first |
| Positive graph-native realization of the whole harmonic operation correspondence | OPEN / candidate program | Fourier change of basis is signed; count equality does not transport positivity or guards |

The proposed grander organizing principle is therefore precise: **count equality is the dimension shadow of a phase-filtered harmonic decomposition; operation transport is available when the source difference kernel survives the quotient, and the exact kernel of transport is visible on the surviving dual modes.** The two independent failure modes are (i) losing stabilizer phase and (ii) replacing a faithful weighted/source-aware operation by an unsupported unweighted or context-free one. The codec and marked fan show a separate repair: retain the operation's actual context and accept a new target operation, rather than pretending a count bijection preserves the old one.
