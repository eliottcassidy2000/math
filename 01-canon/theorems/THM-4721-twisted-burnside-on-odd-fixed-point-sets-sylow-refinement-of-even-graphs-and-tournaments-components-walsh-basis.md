---
id: THM-4721
title: "The twisted Burnside principle on the fixed-point sets of odd-order groups, and its consequences for even graphs versus tournaments: for every odd-order P <= S_n and every Gamma <= N(P), the Gamma-orbits of P-invariant tournaments are equinumerous with the Gamma-orbits of Gamma-even P-invariant graphs (parity lemma: the F_2-sum of the P-translates of the antipodal graph of an even cycle of h is <P,h>-invariant and h-odd). Hence (i) the Sylow refinement of Royle-Praeger-Glasby-Freedman-Devillers: for every odd prime p and Sylow p-subgroup P of S_n, the tournaments with a P-symmetry and the even graphs with a P-symmetry are equinumerous; at n = p: circulant tournaments and even circulant graphs of prime order are equinumerous, t_p = (1/(p-1)) sum_{d | p-1, d odd} phi(d) 2^{(p-1)/(2d)}, with an AGL(1,p)-equivariant surjection inducing a bijection on types; (ii) the Gamma-refined identity for every Gamma <= S_n (natural equivariant bijections exist for odd |Gamma|) and the coloured identity: k-arc-coloured tournaments = even 2k-coloured complete graphs (HYP-9339 (3)); (iii) a graph is even iff all its components are, so A334335 = connected even graphs (von Broemssen's 'it appears'); (iv) the signed even-subgraph counts Phi_X are a Z[1/2]-basis of the tournament invariants, with two orthogonality relations, det = +-2^k, and prod_T |Aut T| / prod_{X even} |Aut X| is a power of 2"
status: >
  PROVED (pre-audit):
  - Lemma P (parity lemma) and Theorem R (restricted twisted Burnside principle, every odd-order P and Gamma <= N(P));
  - Theorem S (Sylow refinement: tau_T(P) = tau_E(P) for every Sylow p-subgroup P of S_n, p odd) with Corollary C
    (prime circulants: t_p = e_p, the formula, the refinement per odd multiplier subgroup K, existence of an
    AGL(1,p)-equivariant surjection inducing a type bijection by THM-4531's Lemma N1 and Hall on upsets);
  - Theorem G' (P = 1: the Gamma-refined identity for every Gamma <= S_n; Gamma-equivariant bijections after choosing a Gamma-invariant reference tournament for
    odd |Gamma|; none from Gamma-even graphs to tournaments for even |Gamma|) and Theorem K (coloured identity);
  - Lemma E (components) and Theorem A (A334335 = connected even graphs; Euler(connected even) = INVERT(strong));
  - Theorem W (Walsh basis over Z[1/2], orthogonality relations, det M = +-2^k, odd parts of the Aut-products agree);
  - Proposition Q (class identity: sum over tournament types of the number of Aut-classes of elements of odd type
    lambda = the same over even graphs + the centraliser-even classes over odd graphs).
  FINITE-EXACT: tau_T(H) = tau_E(H) for every odd-order subgroup class H of S_n, n <= 8, and for fifteen classes at
  n = 9 including the Sylow 3-subgroup Z3 wr Z3 of order 81, Z3^3, Z9, Z15, Z7:Z3 and four actions of Z3xZ3
  (subgroups_fix.out, subgroups_fix_n8.out, subgroups_fix_n9.out; e.g. 2954 = 2954 for Z3 of type (3,1^6) at n = 9);
  the Theorem S closed form reproduces the census counts for the Sylow 3-, 5-, 7-subgroups of S_6, S_7, S_8
  (sylow_formula_check.out); parity lemma and Theorem R for 440 pairs (P,h) and Gamma in
  {C(P), N(P)}, n <= 7 (parity_lemma_check.out); t_p = e_p for all 45 primes p <= 200 by the formulas, by brute force
  for p <= 31 with the refined K-counts, and the evenness criterion against networkx for p <= 13
  (circulant_prime.out); the coloured identity for 2, 4, 6 colours and n <= 4 (coloured_even.out); the census n <= 8
  (census_n7.out, census_n8.out: even = tournaments, labelled even 1, 4, 20, 396, 16432, 1312056, 200582824,
  connected even 1, 0, 1, 2, 8, 43, 398, 6413 = A334335, self-converse sums = A002785, Hall deficits 1, 2, 5, 12 as in
  THM-4531); orthogonality and det for n <= 5 (walsh_check.out), odd parts of the Aut-products for n <= 8.
  CITED: Royle, Praeger, Glasby, Freedman, Devillers, J. Algebraic Combin. 57 (2023) 515-524, arXiv:2204.01947
  (definition of even graphs, Lemma 4.3: a self-paired cycle of g_A carries an odd number of inversions, Theorem 4.4);
  THM-4531 (Theorem G, Lemma N1, Theorem O3; "natural" = S_n-equivariant); THM-4524 (E1); Sylow and Frattini;
  Turner 1967 / Elspas-Turner 1970 (isomorphic circulants of prime order are multiplier-equivalent; re-derived here
  from Sylow); Alspach 1973 (Aut of a prime circulant is Z_p x| Stab(S); used only in the brute-force evenness
  criterion, not in the proofs); OEIS A334335 (von Broemssen 2020, comment "it appears that a(n) is the number of
  connected even graphs", verified n <= 10) and A000568.
  CURRENT REPAIR (2026-10-20): HYP-9340(a) is REFUTED for general odd H (C15 at n=10: 32 versus 31),
  but PROVED for all odd p-groups by the linked Sylow-action Walsh addendum. The linked derivative-path
  construction solves the all-odd-prime circulant type-bijection target; fixed multiplier rules still fail.
  FURTHER EXTENSION (2026-10-20): odd-pronormal-base/one-prime normalizer-quotient matching proves the
  stated mixed-prime survivor. Abstract cyclic-group naturality has the linked p=3 mod4 / p=5 /
  p=1 mod4,p>5 dichotomy. HYP-9172(a) is now proved for every n>=5 in its precise equivariant scope.
  Notation: sgn_X(g) = (-1)^{# edges {u<v} of X with g(u) > g(v)} for g in Aut X (Royle et al.'s sign, THM-4524's
  eps); X is even iff sgn_X = 1 on Aut X; Fix_T(P), Fix_G(P) = P-invariant labelled tournaments / graphs on [n];
  tau_T(H) = number of isomorphism types of tournaments whose automorphism group contains a conjugate of H, tau_E(H)
  the same for even graphs; c_2(G) = number of orbits of G on unordered pairs; r(h, X) = parity of the number of
  edges of X reversed by h (defined for every permutation h and every graph X, a cocycle).
source: mac-mini-2026-10-19-evengraphs (owner's prompt: an exact bijection between tournaments and even graphs; past sequence work; the two attached blueprints; Minkowski's ?)
depends_on:
  - 01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md
  - 01-canon/theorems/THM-4720-k-valued-tournaments-by-burnside-the-parity-dichotomy-of-the-tie-state-and-three-hypertournaments.md
related:
  - 05-knowledge/hypotheses/HYP-9340-refined-equinumerosity-for-every-odd-order-symmetry-and-what-a-tournaments-even-graphs-bijection-can-be.md
  - 05-knowledge/hypotheses/HYP-9339-the-two-three-axis-numbers-fixed-points-and-relations.md
  - 05-knowledge/hypotheses/HYP-9172-no-natural-bijection-for-twisted-mallows-sloane-at-any-n-ge-5.md
  - 01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md
note: 05-knowledge/results/even_graphs_20261019_sylow_refinement_components_walsh_duality.md
scripts: 04-computation/experiments/even_graphs_20261019/{census,classes_identity,subgroups_fix,parity_lemma_check,circulant_prime,coloured_even,walsh_check,odd_graphs_bfile,sylow_formula_check}.py
outputs: 04-computation/experiments/even_graphs_20261019/*.out, connected_even_graphs_bfile.txt (A334335 to n = 128), odd_graphs_bfile.txt (n <= 60)
---

# THM-4721 — the twisted Burnside principle on odd fixed-point sets; the Sylow refinement of "tournaments and even graphs are equinumerous"

**PROVED (pre-audit)**, with the finite checks listed in the status. Full note:
[even_graphs_20261019_sylow_refinement_components_walsh_duality](../../05-knowledge/results/even_graphs_20261019_sylow_refinement_components_walsh_duality.md).

**Current continuation, 2026-10-20.** [Pronormal descent and quotient refinement](../../05-knowledge/results/pronormal_quotient_bijections_20261020/PRONORMAL_QUOTIENT_THEOREM.md) extend the equality to odd pronormal bases and their one-prime normalizer-quotient extensions. The unrestricted odd-H assertion remains refuted. [Generator-free circulant naturality](../../05-knowledge/results/generator_free_circulants_20261020/GENERATOR_FREE_NATURALITY.md) gives an exact abstract-group naturality dichotomy, while [the uniform obstruction](../../05-knowledge/results/uniform_equivariance_obstruction_20261020/UNIFORM_EQUIVARIANCE_OBSTRUCTION.md) resolves HYP-9172(a) for all n>=5 without excluding arbitrary unlabelled type bijections. These links supplement, rather than alter, the proofs and old finite evidence below.

## 0. Setting

Graphs and tournaments live on `[n]`. For a permutation `h` and a graph `X` let `r(h, X)` be the parity of the number of edges `{u < v}` of `X` with `h(u) > h(v)`. This is a cocycle: `r(st, X) = r(t, X) + r(s, tX)`, and it is additive under symmetric difference. For `h` in `Aut X` it is Royle et al.'s sign: `sgn_X(h) = (-1)^{r(h,X)}`, a homomorphism `Aut X -> {+-1}` whose kernel has index 1 (`X` even) or 2 (`X` odd). Royle et al.'s Lemma 4.3 gives `r(h, X) = (number of h-orbits on the edge set of X that are self-paired, i.e. some power of h swaps the two ends) mod 2`. For a group `G` of permutations, `X_O` denotes the graph whose edge set is a `G`-orbit `O` on pairs, and `chi_O(h) = sgn_{X_O}(h)` for `h` in `G`; each `chi_O` is a character of `G` and `sgn_X(h) = prod_{O subset X} chi_O(h)` for `X` in `Fix_G(G)`.

The tool behind everything is the **signed count over an invariant family**: for a group `G` and `h` in `G`,

    sum_{X in Fix_G(G)} sgn_X(h) = prod_{O} (1 + chi_O(h)) = 2^{c_2(G)} if chi_O(h) = 1 for every G-orbit O, and 0 otherwise.     (0.1)

The tournament side is `|Fix_T(G)| = 2^{c_2(G)}` if `|G|` is odd and `0` otherwise (Royle et al., Corollary 3.2).

## 1. The parity lemma

**Lemma P.** Let `P <= S_n` have odd order, let `h` normalise `P` and have even order, and let `G = <P, h> = P<h>`. Then some `G`-invariant graph `X` has `r(h, X) = 1`. Explicitly: let `C` be an even cycle of `h`, `Y_C` the graph of the `|C|/2` antipodal pairs `{u, h^{|C|/2} u}` of `C` (one self-paired `h`-orbit, so `r(h, Y_C) = 1`), and

    X = sum over k in P of k(Y_C)   (symmetric difference).

Then `X` is `G`-invariant and `r(h, X) = 1`.

*Proof.* `G`-invariance: `P` permutes the summands; `h k Y_C = (h k h^{-1}) h Y_C = k' Y_C` with `k' in P`, since `Y_C` is `<h>`-invariant. The sign: by additivity `r(h, X) = sum_k r(h, kY_C)`. By the cocycle, `r(hk, Y_C) = r(k, Y_C) + r(h, kY_C)`, and writing `hk = k'h`, `r(k'h, Y_C) = r(h, Y_C) + r(k', hY_C) = r(h, Y_C) + r(k', Y_C)`. So `r(h, kY_C) = r(h, Y_C) + r(k', Y_C) - r(k, Y_C)`. Summing over `k in P`, the map `k -> k' = hkh^{-1}` is a bijection of `P`, so the last two sums cancel and `r(h, X) = |P| r(h, Y_C) = 1 (mod 2)`. In particular `X` is not empty. ∎

Consequently, for such `G` and `h`, the functional `r(h, .)` on the `F_2`-space `Fix_G(G)` is nonzero, so (0.1) vanishes: `sum_{X in Fix_G(G)} sgn_X(h) = 0 = |Fix_T(G)|`.

## 2. The restricted twisted Burnside principle

**Theorem R.** Let `P <= S_n` have odd order and let `Gamma <= N_{S_n}(P)`. Call a `P`-invariant graph `X` *Gamma-even* if `sgn_X(h) = 1` for every `h` in `Gamma ∩ Aut X`. Then

    #(Gamma-orbits on Fix_T(P)) = #(Gamma-orbits on the Gamma-even members of Fix_G(P)).

*Proof.* Cauchy-Frobenius for `Gamma` on `Fix_T(P)` gives `(1/|Gamma|) sum_{h in Gamma} |Fix_T(<P,h>)|`. The signed Cauchy-Frobenius for `Gamma` on `Fix_G(P)` (each orbit contributes `(1/|Stab|) sum_{h in Stab} sgn_X(h) = [X is Gamma-even]`) gives `(1/|Gamma|) sum_{h in Gamma} sum_{X in Fix_G(<P,h>)} sgn_X(h)`. The terms agree: if `h` has odd order, `<P,h> = P<h>` has odd order, `h` reverses nothing, and both terms are `2^{c_2(<P,h>)}`; if `h` has even order, both are `0` by Lemma P. ∎

`P = 1` is **Theorem G'**: for every `Gamma <= S_n`, `#(Gamma-orbits on tournaments) = #(Gamma-orbits on Gamma-even graphs)`. `Gamma = S_n` is Royle et al.; `Gamma = S_{n-1}` is the rooted form; `Gamma` trivial is the labelled triviality `2^{C(n,2)} = 2^{C(n,2)}`. If `|Gamma|` is odd every graph is `Gamma`-even and the map `X -> T_0 reversed on X` (`T_0` any `Gamma`-invariant tournament, e.g. the rotational tournament for `Gamma = Z_n`, `n` odd) is a `Gamma`-equivariant bijection. If `|Gamma|` is even, no `Gamma`-equivariant map from the `Gamma`-even graphs to tournaments is bijective on orbits (the empty graph has stabiliser `Gamma`, no tournament does), and in the other direction THM-4531's Theorem O3 applies to `Gamma = S_n`.

**Theorem K (coloured).** Give each pair a colour in `Z/2k` and call an automorphism `g` of a colouring *odd* if the number of pairs reversed by `g` that carry an odd colour is odd. Then the even `2k`-colourings up to isomorphism are equinumerous with the `k`-arc-coloured tournaments of THM-4720: both are `sum_{lambda odd} (2k)^{c_2(lambda)} / z_lambda`. *Proof.* The sign `(-1)^{(sum of colours of reversed pairs)}` is a character of the automorphism group (the mod-2 reduction of the colouring makes it the sign of the graph of odd-coloured pairs), and in the signed count a self-paired orbit of `g` contributes `sum_{a in Z/2k} (-1)^a = 0` while a non-self-paired orbit contributes `2k`. So the signed Burnside sum is `2k^{c_2(g)}` over odd-order `g` and `0` otherwise, which is THM-4720's `T_{2k}(n)` (there `o_+ = c_2` for odd cycle types). ∎ This settles HYP-9339 (3) in its correct form (even `2k`-colourings, not `k`-colourings). Checked exactly for `2k = 2, 4, 6` and `n <= 4` (`coloured_even.out`: 4, 176, 1956 at `n = 4`).

## 3. The Sylow refinement

**Theorem S.** Let `p` be an odd prime and `P` a Sylow `p`-subgroup of `S_n`. Then

    tau_T(P) = tau_E(P):

the number of isomorphism types of tournaments admitting an automorphism group containing a conjugate of `P` equals the number of isomorphism types of even graphs admitting one. Both equal `(1/|N(P)|) sum_{h in N(P), |<P,h>| odd} 2^{c_2(<P,h>)}`.

*Proof.* (a) *Types are `N(P)`-orbits.* If `x, x' in Fix(P)` (tournaments or graphs) are isomorphic via `s`, then `sPs^{-1}` and `P` are Sylow `p`-subgroups of `Aut(x')`, hence conjugate by some `a` in `Aut(x')`; then `as in N(P)` and `(as)x = x'`. So the types meeting `Fix(P)` are exactly the `N(P)`-orbits on `Fix(P)`. (b) *Frattini.* If `X in Fix_G(P)` is odd, let `A = Aut X` and `A^+ = ker sgn_X` (index 2). `P <= A^+` is a Sylow `p`-subgroup of the normal subgroup `A^+`, so `A = A^+ N_A(P)`; as `A != A^+`, `N_A(P)` contains an odd automorphism. Hence on `Fix_G(P)`, `N(P)`-even means even. (c) Theorem R with `Gamma = N(P)`. ∎

**Corollary C (prime order).** For `n = p` prime, `P = Z_p`, `N(P) = AGL(1,p)`: the circulant tournaments of order `p` up to isomorphism are equinumerous with the even circulant graphs of order `p` up to isomorphism, and

    t_p = e_p = (1/(p-1)) sum_{d | p-1, d odd} phi(d) 2^{(p-1)/(2d)}   (1, 1, 2, 4, 6, 16, 30, 94, 586, 1096, 7286, ... for p = 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37).

(The `N(P)`-orbit description is Turner's theorem for prime order, re-derived by (a).) Refined: for every odd subgroup `K <= Z_p^*`, apply Theorem S' = Theorem R to `P_K = Z_p x| K` (odd order, normal in `AGL(1,p)`; (a) and (b) go through with the Sylow `Z_p`): the circulant tournament types with multiplier stabiliser containing `K` are equinumerous with the even circulant types with stabiliser containing `K`. The odd subgroups of `Z_p^*` form a divisor lattice, so inclusion-exclusion over joins gives equality on every upset, Hall's condition holds for THM-4531's stabiliser-containment graph, and by Lemma N1 an **`AGL(1,p)`-equivariant surjection from labelled circulant tournaments onto labelled even circulant graphs, inducing a bijection on types, exists for every odd prime `p`** (for `p = 7`: `QR_7 -> empty`, the rigid circulant `-> complement of C_7`). A labelled bijection is impossible already at `p=3` (two tournaments, one even graph); stabilizer containment gives an orbit-bijective map, not stabilizer equality. The [audited derivative-path construction](../../05-knowledge/results/circulant_necklace_bijection_20261020/CIRCULANT_NECKLACE_BIJECTION.md) now supplies an explicit all-odd-prime rule. The original fixed-multiplier tests remain: `S -> {x : [x in S] != [cx in S]}` is a bijection for some `p` and not others for every fixed `c in {+-2, +-3, +-4}` (`circulant_prime.out`).

The identity `tau_T(H) = tau_E(H)` holds in the census for **every** odd-order subgroup class `H` of `S_n`, `n <= 8`, and for the fifteen classes tested at `n = 9` (not only Sylow ones; `subgroups_fix*.out`): e.g. `Z_3` of type `(3,1^5)`: 280 = 280; `Z_3 x Z_3` in `S_8`: 16 = 16; `Z_15`: 2 = 2. For non-Sylow `H` the two ingredients (a), (b) fail separately (the tournament with `Aut = Z_3 x Z_3` on 6 vertices has two classes of 3-cycle modules; `K_{3,3}` is odd but `N(Z_3)`-even) and cancel in those tested cases; the universal HYP-9340(a) is now [refuted at n=10](../../05-knowledge/results/c15_symmetry_counterexample_20261020/C15_ORDER_TEN_COUNTEREXAMPLE.md), while [all odd p-groups are proved](../../05-knowledge/results/prime_local_walsh_bijection_20261020/SYLOW_ACTION_ADDENDUM.md). What is proved for every odd `P` is the class form, **Proposition Q**: `sum_T m_T(P) = sum_{X even} m_X(P) + sum_{X odd, N(P)-even} m_X(P)`, with `m_x(P)` the number of `Aut(x)`-classes of `S_n`-conjugates of `P` inside `Aut(x)` (Theorem R with `Gamma = N(P)`, orbits counted with multiplicity), and its element form for odd cycle types (`classes_identity.out`, n <= 7, all OK).

## 4. Components and A334335

Repository precedence: the component-even lemma and the all-order SCC/Lyndon type codec were already proved in the [11 October compositional note](../../05-knowledge/results/tournament_even_bijection_20261011c.md), sections 2–4. The argument and sequence extension below are retained.

**Lemma E.** A graph is even iff every connected component is even. *Proof.* An odd automorphism of a component extends by the identity to an odd automorphism of `X`. Conversely let `g` be an odd automorphism of `X`, so the number of cycles of `g` whose antipodal pairs are edges ("antipodal-edge cycles") is odd. `g` permutes the components; group the components into `g`-cycles of components, of lengths `k_1, ..., k_s`. If `C` is an antipodal-edge cycle of `g` of length `L`, its antipodal pair `{u, g^{L/2} u}` is an edge, so both ends lie in one component, so the component-cycle through `u`'s component has length `k` dividing `L/2`, and `C` visits each of its `k` components in `L/k` (even) vertices. For a component `H` of that component-cycle, `g^k` restricts to an automorphism of `H`, `C ∩ H` is one `g^k`-cycle of even length `L/k`, and its antipodal pairs are the antipodal pairs of `C` lying in `H` (again `{v, g^{L/2} v}`), hence edges; conversely every antipodal-edge cycle of `g^k|_H` arises this way from the `g`-cycle through it. So the number of antipodal-edge cycles of `g^k|_H` equals the number `a_i` of antipodal-edge cycles of `g` attached to the component-cycle `i` of `H`. As `sum_i a_i` is odd, some `a_i` is odd, and `g^{k_i}` restricted to any component of that component-cycle is an odd automorphism of that component. ∎

**Theorem A.** The generating function of even graphs is the Euler transform of that of connected even graphs; hence the inverse Euler transform of A000568, OEIS A334335 (von Broemssen 2020, "it appears that a(n) is the number of connected even graphs"), counts the connected even graphs: 1, 0, 1, 2, 8, 43, 398, 6413, 184596, ... (census agreement to `n = 8`; `connected_even_graphs_bfile.txt` extends A334335's 19 terms to `n = 128` from the repo's A000568 b-file). Since tournaments are also sequences of strong tournaments, `EULER(connected even) = INVERT(strong tournaments, A051337)`: multisets of connected even graphs and sequences of strong tournaments are equinumerous.

## 5. The Walsh basis: the exact linear bijection

For an even labelled graph `X` and a tournament `T` let `phi_X(T) = prod_{e in X} T(e)` (`T(e) = +1` if `T` orients `e` upward in the label order), and `Phi_X(T) = sum_{X' in the S_n-orbit of X} eps(X') phi_{X'}(T)`, the signed count of the copies of `X` in `K_n` by the parity class of the orientation `T` induces on them. `Phi_X` is well defined and nonzero exactly when `X` is even (the cocycle sign on the orbit is single-valued iff `sgn_X` is trivial on `Aut X`); for odd `X` the orbit sum vanishes. Examples: `Phi_empty = 1`; `Phi_{P_3}(T) = #(directed 2-paths) - #(2-stars)`.

**Theorem W.** (i) The `Phi_X`, `X` even, form a basis of the isomorphism-invariant functions on tournaments, over `Z[1/2]`. (ii) With `m = C(n,2)`, `|orb T| = n!/|Aut T|`, `|orb X| = n!/|Aut X|`:

    sum_T |orb T| Phi_X(T) Phi_{X'}(T) = 2^m |orb X| delta_{X X'},      sum_{X even} |Aut X| Phi_X(T) Phi_X(T') = 2^m |Aut T| delta_{T T'}.

(iii) The matrix `M = (Phi_X(T))` (types x even graphs) has `det M = +-2^k`, `k = mN/2 + (1/2) log_2 (prod_T |Aut T| / prod_X |Aut X|)`, `N` = number of types; in particular `prod_T |Aut T|` and `prod_{X even} |Aut X|` have the same odd part (n <= 8: `3^328 5^12 7^4` at n = 8).

*Proof.* Over `R = Z[1/2]` the Walsh characters `phi_X` are an `R`-basis of `R[tournaments]` and `S_n` acts monomially, `g phi_X = eps(g,X) phi_{gX}` (THM-4531, Theorem G); the invariants are spanned by the orbit sums, which are `R Phi_X` for `X` even and `0` for `X` odd, while the orbit indicators are another `R`-basis of the invariants; so `M` is invertible over `R`. (ii) is Parseval for the Walsh basis summed over orbits; (iii) takes determinants in `M^T D_T M = 2^m D_X`. ∎ (`walsh_check.out`: n = 3, 4, 5: det = -2^2, -2^8, 2^49.)

So the exact object relating tournament types and even graphs is this change of basis: **every isomorphism invariant of tournaments is a unique `Z[1/2]`-combination of signed even-subgraph counts**. It is canonical and natural; it is not a bijection of sets, and THM-4531 shows no `S_n`-equivariant set bijection exists for `5 <= n <= 15`.

## 6. Two negative facts

- The per-permutation identity behind Theorem S is global, not termwise: at `n = 4`, type `(3,1)`, the labelled tournaments with a 3-cycle automorphism number 16, the graphs 10, and the difference is cancelled by the three double transpositions (`+6 - 3 x 2`).
- No bijection preserves `|Aut|` (THM-4531), and the `|Aut|` multisets differ already in their odd parts (n = 7: even graphs have `|Aut|` with odd parts 15, 45, 315); only the products agree (Theorem W (iii)).

## 7. Loose links to the attached blueprints

Both attached proofs turn on a parity class preserved by a symmetry: the strip colouring `floor(y/h) mod 2` that makes equilateral triangles the only exception in the Erdős-173 blueprint, and the sign classes of countable colourings of `C`-linear graphs. The dichotomy here is of the same shape: a graph is even exactly when its automorphisms preserve the two orientation-parity classes, and the signed Burnside sums are the bookkeeping of that preservation. No theorem transfers; Minkowski's question mark was looked at (Walsh characters are the dyadic Fourier basis) and nothing honest was found.
