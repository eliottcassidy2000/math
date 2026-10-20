---
id: HYP-9340
title: "(a) REFUTED for general odd-order H: C15 on ten vertices gives 32 tournament types versus 31 even-graph types; PROVED for every odd p-group H and the stated odd-pronormal-base quotient extensions. (b) Prime-local Walsh and all-prime circulant constructions; simultaneous preservation of all odd-H predicates is impossible already at n=5; labelled affine surjection versus type bijection. (c) Sequence conjectures: odd graphs 0, 1, 2, 7, 22, 100, 588, 5466 (b-file to n = 60), labelled even graphs 1, 1, 4, 20, 396, 16432, 1312056, 200582824, connected odd graphs 0, 1, 1, 4, 13, 69, 455, 4704 are not in OEIS (2026-10-19); A334335 = connected even graphs extended to n = 128"
status: >
  (a) REFUTED (PROVED + FINITE-EXACT): H=C15 with vertex-orbit sizes (3,5,1,1) at n=10 has
  tau_T(H)=32 and tau_E(H)=31. The earlier n<=8 and selected n=9 counts remain valid finite evidence.
  PROVED survivor: equality for every odd p-group H<=S_n, via a matching preserving its Sylow p action.
  PROVED extension: odd pronormal K, K <= H <= N(K), and H/K a q-group for one fixed odd prime q;
  one matching on the sector admitting K preserves all these extension predicates simultaneously.
  PROVED classification plus FINITE-EXACT counts now cover every odd mixed-prime subgroup class through n=12.
  The least failure degree is ten; the only failure classes through twelve are C15 on 5+3+1^f, f=2,3,4.
  THM-4721's ambient Sylow theorem and multiplicity-weighted normalizer identity remain valid.
  (b) PROVED: for each fixed odd prime p, a Walsh bijection preserves all p-group predicates simultaneously;
  no single bijection preserves all odd-H predicates (already C3 and C5 at n=5). An explicit all-odd-prime
  circulant construction is affine-equivariant and surjective on labelled objects, bijective on types.
  An odd-order Gamma-equivariant labelled bijection requires a chosen Gamma-invariant tournament T0;
  it is not a unique canonical choice. Abstract multiplicative-group naturality has a sharp dichotomy:
  a generator-free rule for p=3 mod4 and the p=5 exception; no natural type injection for p=1 mod4, p>5.
  Separately, Boolean-affine multiplier-equivariant type bijections exist iff p=3 mod4 or p=5.
  Necessary Boolean degree is unbounded over primes; exact minima are two at p=13 and three at p=17.
  Nonlinear field-additive algorithms are not excluded. The displayed parity
  projection itself is not uniformly first-order definable; definability of some other prime-field rule remains OPEN.
  (c) FINITE-EXACT (census n <= 8, Burnside formula for odd graphs to n = 60, inverse Euler transform of the repo's
  A000568 b-file to n = 128); OEIS lookups 2026-10-19 (JSON API, generic UA): odd graphs, labelled even graphs,
  connected odd graphs: no match; A334335 found with 19 terms and the comment "it appears ... connected even graphs
  (verified n <= 10)", now PROVED by THM-4721 Theorem A.
source: mac-mini-2026-10-19-evengraphs
related:
  - 01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md
  - 01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md
  - 05-knowledge/hypotheses/HYP-9172-no-natural-bijection-for-twisted-mallows-sloane-at-any-n-ge-5.md
  - 05-knowledge/hypotheses/HYP-9339-the-two-three-axis-numbers-fixed-points-and-relations.md
note: 05-knowledge/results/even_graphs_20261019_sylow_refinement_components_walsh_duality.md
---

# HYP-9340 — refined equinumerosity for every odd-order symmetry, and what an exact bijection can be

## (a) The refined conjecture, now refuted in general

For `H <= S_n` of odd order let `tau_T(H)` (resp. `tau_E(H)`) be the number of isomorphism types of tournaments (resp. sign-even graphs) whose automorphism group contains a conjugate of `H`.

**REFUTED, 2026-10-20.** The historical conjecture asserted `tau_T(H)=tau_E(H)` for every odd-order `H`. For `H=<(0 1 2)(3 4 5 6 7)> <= S_10`, fixing 8 and 9, the exact values are **32 and 31**. The [complete proof and 512-by-512 invariant census](../results/c15_symmetry_counterexample_20261020/C15_ORDER_TEN_COUNTEREXAMPLE.md), with independent review in that packet, identify the single obstruction as `K_{5,5}`. Two fixed vertices join the 3-orbit's independent twin class; the new equal-sized part swap is odd but does not normalize the chosen `H`. Thus normalizer-evenness does not imply full evenness.

**PROVED survivor.** For every fixed odd prime `p` and every `p`-group `H<=S_n`, equality holds. The [Sylow-action Walsh matching](../results/prime_local_walsh_bijection_20261020/SYLOW_ACTION_ADDENDUM.md) preserves the `S_n`-conjugacy class of the full Sylow `p` action, hence all `p`-group symmetry predicates simultaneously. This extends THM-4721's ambient-Sylow result; it does not combine different primes into one matching.

**PROVED mixed-prime survivor (2026-10-20).** The [pronormal quotient theorem](../results/pronormal_quotient_bijections_20261020/PRONORMAL_QUOTIENT_THEOREM.md) applies to an odd subgroup K pronormal in S_n. Pronormality identifies K-fixed normalizer orbits with full isomorphism types and makes normalizer-evenness equivalent to full sign-evenness. For each fixed odd prime q, one matching preserves the complete Sylow-q action in N(K)/K and hence every predicate for `K <= H <= N(K)` with `H/K` a q-group. H may have mixed prime divisors and need not itself be pronormal. This is a sufficient ambient descent condition, not a proof for arbitrary odd H. Abstract normality of K in H is not enough. The inherited all-p-group theorem already includes C3 on six vertices; the exact single-3-cycle and synchronized-two-3-cycle counts are respectively 11 and 3 on each side.

The earlier exact evidence is retained: every tested odd-order subgroup class through `n=8`, and fifteen classes at `n=9`, including `(3,1^6): 2954/2954`, four `C3 x C3` actions, `C3 wr C3`, `C9` and `C15`. In particular `5+3` at order eight gives `2/2`, and `5+3+1` at order nine gives `8/8`. Those selected finite counts did not imply the universal claim. The later [complete mixed-prime census through ten](../results/odd_mixed_prime_census_20261020/README.md) proves that ten is the least failure degree, with a unique failure class there. Its [complete extension through twelve](../results/mixed_prime_census_extension_20261020/README.md) has 8 and 10 classes at n=11 and 12; the only failures remain C15 on 5+3 with fixed vertices, with counts 172/168 and 1376/1360. This does not classify n>=13.

The later [double-triple descent boundary](../results/double_triple_symmetry_boundary_20261020/README.md) structurally proves equality for the synchronized C5 times C3 and independent C5 times C3^2 actions on 5+3+3 with zero or one fixed vertex. With two fixed vertices, both have exact deficit six (selected finite totals 768/762 and 256/250). This is a theorem for two specified degree-thirteen embeddings, not a complete census there.

THM-4721 Proposition Q remains the correct multiplicity-weighted statement:
`sum_T m_T(H) = sum_{X even} m_X(H) + sum_{X odd, N(H)-even} m_X(H)`, where `m_x(H)` counts `Aut(x)`-classes of conjugates of `H` inside `Aut(x)`. Its extra normalizer-even term cannot simply be erased. The per-permutation identity also remains false, as recorded in THM-4721 section 6.

## (b) What an exact bijection can be

1. **The canonical linear object survives.** The signed even-subgraph counts `Phi_X` and type indicators are `Z[1/2]`-bases of the invariant space, with Walsh orthogonality and determinant a signed power of two (THM-4721 Theorem W). A linear change of basis is not itself a set bijection.
2. **Keep labelled and orbit maps separate.** For odd-order `Gamma`, choose a `Gamma`-invariant tournament `T0`; reversing `T0` along each graph gives a `Gamma`-equivariant labelled bijection. Different choices of `T0` need not give the same map. For prime circulants the [explicit derivative-path construction](../results/circulant_necklace_bijection_20261020/CIRCULANT_NECKLACE_BIJECTION.md) gives an `AGL(1,p)`-equivariant surjection from labelled circulant tournaments onto labelled sign-even circulant graphs, inducing a bijection of types, for every odd prime. It is not a labelled bijection: at `p=3` the labelled counts are two versus one. The [uniform obstruction](../results/uniform_equivariance_obstruction_20261020/UNIFORM_EQUIVARIANCE_OBSTRUCTION.md) now proves HYP-9172(a) for every n>=5, precisely for S_n-equivariant maps inducing a type bijection; arbitrary unlabelled bijections remain possible.
3. **Separate one-prime preservation from a joint coupling.** A single matching preserves every `p`-group predicate for one chosen odd `p`. Simultaneous preservation of all odd-`H` predicates is impossible already at `n=5`: the empty graph admits both `C3` and `C5`, while no five-vertex tournament admits both. Equal individual counts do not imply equal joint intersections. The [Sylow-action addendum, section 5](../results/prime_local_walsh_bijection_20261020/SYLOW_ACTION_ADDENDUM.md) records this obstruction; the order-ten counterexample also disproves an individual mixed-prime equality.
4. **The constructive circulant target is solved in its stated algorithmic sense.** Fixed multiplier derivatives still fail at the documented examples, but a bounded all-prime derivative-path rule exists. Write `p-1=2^r q`, `q` odd. The `K_max`-fixed tournament sector has `2^(2^(r-1))/2^r` multiplier types. Only for `r<=2` is its unique type forced to the empty graph. At `p=17`, `K_max` is trivial and all sixteen types lie in that sector. The [generator-free naturality theorem](../results/generator_free_circulants_20261020/GENERATOR_FREE_NATURALITY.md) now separates those questions. For p=3 mod4, the parity projection on nonzero squares is affine-equivariant, exactly two-to-one on labelled circulants and bijective on types. The p=5 type map is unique. For p=1 mod4, p>5, no type injection can respect all automorphisms of the abstract multiplicative cyclic group; this does not rule out using field addition. A cited AC0 parity bound excludes a uniform first-order ring-plus-input formula for this particular projection, not for every possible bijection. Over even-degree finite-field extensions, a Frobenius obstruction concerns the explicitly defined multiplier/necklace target, not automatically full sign-even graph types or cyclic rings Z/p^r.

**Later Boolean-degree boundary (2026-10-20).** The [affine ideal theorem](../results/affine_circulant_ideal_obstruction_20261020/AFFINE_CIRCULANT_IDEAL_OBSTRUCTION.md) classifies every Boolean-affine multiplier-equivariant rule: a type bijection exists exactly for p=3 mod4 or p=5. This is a different obstruction from abstract multiplicative-group naturality. The [degree theorem](../results/circulant_boolean_degree_20261020/BOOLEAN_DEGREE_BOUNDARY.md) proves degree >= ceil((2^a-1)/a) when the dyadic part of (p-1)/2 is 2^a>=4, unbounded necessity across primes, and exact minima two at p=13 and three at p=17. The all-prime nonlinear construction survives. The [square-field family](../results/affine_square_field_census_20261020/AFFINE_SQUARE_FIELD_CENSUS.md) also gives tau_T=tau_E=1 for its explicitly embedded F_q-additive semidirect groups; it is not a census of all prime-power-order objects.

## (c) Sequences (candidates for OEIS submission, with the THM-4720 and S48 lists)

| sequence | terms (n = 1, 2, ...) | status |
|---|---|---|
| odd graphs = A000088 − A000568 | 0, 1, 2, 7, 22, 100, 588, 5466, 83132, 2272112, 115244616, 10982861424 | not in OEIS; b-file to n = 60 (`odd_graphs_bfile.txt`) |
| labelled even graphs | 1, 1, 4, 20, 396, 16432, 1312056, 200582824 | not in OEIS; census n <= 8 |
| connected even graphs = A334335 | 1, 0, 1, 2, 8, 43, 398, 6413, 184596, 9540998, ... | interpretation now PROVED; b-file to n = 128 (`connected_even_graphs_bfile.txt`) |
| connected odd graphs = A001349 − A334335 | 0, 1, 1, 4, 13, 69, 455, 4704 | not in OEIS; census n <= 8 |
| strong tournaments | 1, 0, 1, 1, 6, 35, 353, 6008 | A051337; `EULER(A334335) = INVERT(A051337) = A000568` |
| even 4-, 6-coloured complete graphs | = THM-4720's `T_4`, `T_6` | THM-4721 Theorem K |
