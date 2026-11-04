---
id: THM-4871
title: "The C_pq (p, q, 1^k) defect family of HYP-9340 beyond p-groups. For coprime odd prime powers p != q and H = <c_P c_Q> <= S_n of orbit type (p, q, 1^k): (1) the normaliser N = Hol(Z_p) x Hol(Z_q) x Sym(F) has exactly 2 t_p t_q b_k(4) orbits on H-fixed tournaments (circulant multiplier classes on the two orbits, the cross orientation, a 4-coloured tournament on the k fixed points) and an H-fixed graph is N-even iff its two circulant parts are multiplier-even and its 4-coloured fixed-point graph is colour-sign-even; (2) in a permutation group of odd order the supports of two cycles of prime lengths are equal or disjoint, and disjoint for distinct primes; (3) for primes p < q every conjugate of H inside the automorphism group of an H-fixed tournament has p-orbit P or a p-subset of the fixed points and q-orbit Q or a q-subset, so the tournament multiplicity excess E_T vanishes for k < p and equals 4 t_p^2 t_q at k = p; (4) K_{q,q} on Q against P plus q - p fixed points is an N-even sign-odd H-fixed graph, so tau_T(H) - tau_E(H) >= 1 at k = q - p whenever p < q < 2p: HYP-9340 (a) fails for every such prime pair at degree 2q; (5) census (3,5) k <= 6, (3,7) k <= 6, (5,7) k <= 6, (3,11) k <= 4, (3,13) k <= 3, (7,13) k <= 3, (7,11) k <= 4, (11,13) k <= 2, (9,5) k <= 4, (9,7) k <= 3, (9,11) k <= 3, (9,13) k <= 4, (3,25) k <= 1 (62 rows, thirteen families): HYP-9344 (2)'s pre-registered defect(5) = 88 is refuted (84; defect(6) = 672), the first defect sits at k = |p - q| with value 1 in every family, and tau_T >= tau_E in every row; (7) for primes p < q < 2p the defect is 0 for k < q - p and exactly 1 at k = q - p (below the first defect every N-even H-fixed graph is sign-even with one class of conjugates of H; at k = q - p the only sign-odd N-even type is K_{q,q}), so the least degree at which (p, q) violates HYP-9340 (a) is 2q"
status: >
  PROVED: statements 1-4 and 7 (elementary and classical group theory: THM-4721 Theorem R for the orbit identity;
  the sign homomorphism; Jordan's prime-cycle theorem, Burnside's prime-degree theorem, Galois's theorem on
  solvable primitive groups with Feit-Thompson (or Jones 2014) for the lemma; Gallai's modular decomposition for
  statement 7; CRT). FINITE-EXACT: statements 5-6 (cpq_defect_census.py with pynauty 2.8.8.1 canonical forms and
  automorphism generators; every row verifies THM-4721 Theorem R and THM-4728's identity
  tau_T - tau_E = W + E_even + E_wit - E_T; the (3,5) rows k <= 4 reproduce THM-4728 / HYP-9344 exactly; the
  brute-force mode over all 2^c H-fixed objects with the full normaliser action agrees at (3,5) k <= 3
  (cpq_3_5_k4_brute3.out), (9,5) k <= 1 and (3,25) k = 0, and audit A's own definition-level brute force adds
  (3,5) k = 4, (5,7) k <= 3, (3,7) k = 3). REFUTED: HYP-9344 (2)'s prediction defect(5) = 88 (and the alternative 64)
  for C15 (3,5,1^5): the value is 84, confirmed by audit A's independent implementation. CITED: THM-4721
  (Theorem R, the sign cocycle), THM-4728 (defect identity, census through n = 13), THM-4723, HYP-9340,
  HYP-9344, the Codex packets of 2026-10-20; Royle-Praeger-Glasby-Freedman-Devillers for the sign convention;
  Alspach-Parsons / Babai (Z_p is a CI-group); Palfy / Muzychuk (Z_9, Z_25 are not); Gallai 1967. NOT PROVED
  (HYP-9463): the first-defect law for q > 2p and for prime powers (statement 7 covers primes p < q < 2p), the
  formula defect(k) = b_(k-d)(4) - E_T^F(k - d), and the sign tau_T >= tau_E. Audited: audit A (2026-11-04, PASS
  WITH CORRECTIONS, applied: statement 3(c)'s constant 4 t_p t_q -> 4 t_p^2 t_q, the brute-force and row-count
  claims, the Jones citation; MISTAKE-710); statement 7 was added after audit A and is covered by audit B.
session: mac-mini-2026-11-04-oddprimes
source: 05-knowledge/results/cpq_defect_family_mirror_folds_20261104.md
scripts:
  - 04-computation/experiments/cpq_defect_20261104/cpq_defect_census.py (census; --brute cross-check), witness_examples.py, lemma_check.py (the 838-cycle sanity check), cpq_defect_k7.py (streaming version); outputs cpq_<p>_<q>_k<k>.out, cpq_3_5_k4_brute3.out
  - 04-computation/experiments/cpq_defect_20261104/audit_A/ (audit A's independent implementations and re-runs); 05-knowledge/results/cpq_defect_family_audit_A_20261104.md
related:
  - 01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md
  - 01-canon/theorems/THM-4728-hyp-9340-beyond-p-groups-defect-identity-pronormal-frattini-lemma-and-the-smallest-failure-c15-at-n-10.md
  - 01-canon/theorems/THM-4723-sylow-profile-form-of-hyp-9340-and-the-defect-reduction-for-the-three-cycle-type.md
  - 05-knowledge/hypotheses/HYP-9340-refined-equinumerosity-for-every-odd-order-symmetry-and-what-a-tournaments-even-graphs-bijection-can-be.md
  - 05-knowledge/hypotheses/HYP-9344-after-the-sylow-day-defect-growth-of-c15-tau-t-at-least-tau-e-affine-and-first-order-rules-for-prime-circulants-and-prime-power-circulants.md
  - 05-knowledge/hypotheses/HYP-9463-cpq-defect-laws-first-defect-at-the-prime-distance-the-b-k-minus-d-formula-the-sign-and-first-order-circulant-rules.md
  - 05-knowledge/results/c15_symmetry_counterexample_20261020/C15_ORDER_TEN_COUNTEREXAMPLE.md, odd_mixed_prime_census_20261020, mixed_prime_census_extension_20261020 (Codex, 2026-10-20)
  - C. Jordan, Bull. SMF 1 (1873) 40-71 (prime cycles in primitive groups; Wielandt, Finite permutation groups, Thm 13.9); W. Burnside, Theory of groups of finite order, 2nd ed., sect. 169 (transitive groups of prime degree); G. A. Jones, Primitive permutation groups containing a cycle, Bull. Aust. Math. Soc. 89 (2014) 159-165 (arXiv:1209.5169); T. Gallai, Acta Math. Acad. Sci. Hungar. 18 (1967) 25-66 (modular decomposition); B. Alspach, T. D. Parsons, Discrete Math. 25 (1979); P. P. Palfy, European J. Combin. 8 (1987) 35-43; M. Muzychuk, Discrete Math. 167/168 (1997) 497-510
---

# THM-4871 — the (p, q, 1^k) defect family

## Setting

`p != q` are coprime odd prime powers (in statements 2-4 they are primes). `H = <h>`, `h = c_P c_Q`, with
`c_P = (0 1 ... p-1)` on `P = {0, ..., p-1}`, `c_Q = (p ... p+q-1)` on `Q`, and `F = {p+q, ..., n-1}` the `k`
fixed points, `n = p + q + k`. `H` has order `pq` and orbit type `(p, q, 1^k)`. `N = N_(S_n)(H)`.

For `m` an odd prime power, a *skew set* is `S` with `S` and `-S` partitioning `Z_m \ {0}` (a circulant tournament,
arcs `i -> i + s`), a *symmetric set* is `C = -C` not containing `0` (a circulant graph); the unit group `(Z/m)^*`
acts on both by dilation, with `t_m` and `g_m` orbits. `Hol(Z_m) = Z_m x| (Z/m)^*` is the holomorph.

Sign convention (Royle-Praeger-Glasby-Freedman-Devillers; THM-4721, the S45 census): for a graph `X` on `{0, ..., n-1}`
and `g in Aut X`, `sgn_X(g) = (-1)^#{edges {u < v} of X with g(u) > g(v)}`. For automorphisms this is a homomorphism
`Aut X -> {+-1}` (the inversion indicator satisfies `I_(gg')(e) = I_g(g' e) + I_(g')(e)` and `g'` permutes the edge set)
and does not depend on the vertex order (an adjacent transposition of the order changes the count by `[{a,b} in E] +
[{g^(-1)a, g^(-1)b} in E] = 0 mod 2`). `X` is *sign-even* if `sgn_X = 1` on `Aut X` and *`N`-even* if `sgn_X = 1` on
`N cap Aut X`. A circulant graph `(Z_m, C)` is *multiplier-even* if every dilation in its stabiliser `M(C) = {a : aC = C}`
has even sign on it; `e_m` is the number of dilation classes of multiplier-even symmetric sets.

`b_k(4)` is the number of isomorphism types of 4-coloured tournaments on `k` points (Burnside over permutations with odd
cycles only: `1, 4, 16, 88, 704, 8912, 188288, 6847200, 437069312` for `k = 0..8`; not in OEIS 2026-11-04);
`ec_k(4)` the number of types of 4-coloured graphs on `k` points all of whose colour-preserving automorphisms have even
sign (`k` points, graph edges only).

Types with symmetry (THM-4721 Theorem R, THM-4728): `tau_T(H)` (resp. `tau_E(H)`) is the number of isomorphism types of
tournaments (resp. sign-even graphs) whose automorphism group contains an `S_n`-conjugate of `H`; for a type `X`,
`m_X(H)` is the number of `N`-orbits of `H`-fixed objects isomorphic to `X` (= the number of `Aut X`-classes of
`S_n`-conjugates of `H` inside `Aut X`); `E_T = sum_T (m_T - 1)` over tournament types, `E_even` over sign-even graph
types, `W` the number of sign-odd graph types with an `N`-even `H`-fixed representative (witnesses), `E_wit` their
excess counted over `N`-even representatives; the *defect* is `tau_T(H) - tau_E(H) = W + E_even + E_wit - E_T`
(THM-4728 statement 1). All quantities are functions of `k` for fixed `(p, q)`.

## Statements

1. **Normaliser-orbit parametrisation (PROVED).**
   (a) `N = Hol(Z_p) x Hol(Z_q) x Sym(F)`, of order `p phi(p) q phi(q) k!`. The `H`-orbits on pairs are the
   `(p - 1)/2` distance classes `{i, i + d}` inside `P` (size `p` each), the `(q - 1)/2` inside `Q`, the single orbit
   `P x Q` (size `pq`), for each `f in F` the orbits `{f} x P` and `{f} x Q`, and the `C(k,2)` pairs inside `F`; so there
   are `2^c` `H`-fixed tournaments and `2^c` `H`-fixed graphs, `c = (p + q)/2 + 2k + C(k, 2)`.
   (b) An `H`-fixed tournament is a quadruple `(S_P, S_Q, cross, kappa)`: skew sets on `Z_p`, `Z_q`, the orientation
   `P -> Q` or `Q -> P`, and a 4-coloured tournament `kappa` on `F` whose colour at `f` is `(x_f, y_f) = ([P -> f], [Q -> f])`.
   `N` acts by dilating `S_P`, `S_Q` and permuting `kappa`; its orbits are (dilation class of `S_P`) x (of `S_Q`) x
   `cross` x (isomorphism type of `kappa`), so `R(k) := #(N`-orbits on `H`-fixed tournaments`) = 2 t_p t_q b_k(4)`.
   (c) An `H`-fixed graph is `(C_P, C_Q, cross, gamma)` with `gamma` a 4-coloured graph on `F` (`x_f = [f ~ P]`,
   `y_f = [f ~ Q]`), and `N cap Aut X = (Z_p x| M(C_P)) x (Z_q x| M(C_Q)) x Aut(gamma)` with
   `sgn_X(alpha beta pi) = sgn_(C_P)(alpha) sgn_(C_Q)(beta) sgn_gamma(pi)`, each factor computed inside its own block.
   Hence `X` is `N`-even iff `C_P` and `C_Q` are multiplier-even and `gamma` is colour-sign-even, the `N`-even
   `N`-orbits number `2 e_p e_q ec_k(4)`, and Theorem R gives `e_p e_q ec_k(4) = t_p t_q b_k(4)`; in fact `e_m = t_m`
   (Theorem R for `H = <c_m>`, `k = 0`) and `ec_k(4) = b_k(4)` (Theorem R for the trivial group and each Young subgroup).
2. **Disjoint-orbit lemma (PROVED).** Let `G <= Sym(Omega)` have odd order and let `g, g' in G` be cycles of prime
   lengths `r, s` with supports `O, O'`. Then either `O = O'` (and then `r = s`) or `O cap O' = empty`.
3. **Tournament multiplicities for primes `p < q` (PROVED).**
   (a) If `T` is an `H`-fixed tournament and `H' <= Aut T` is `S_n`-conjugate to `H`, then the `p`-orbit of `H'` is `P`
   or a `p`-subset of `F`, and its `q`-orbit is `Q` or a `q`-subset of `F`; if both equal `P`, `Q` then `H' = H`.
   (b) `E_T(k) = 0` for `k < p` (every tournament type with `H`-symmetry has `m_T(H) = 1`).
   (c) `E_T(p) = 4 t_p^2 t_q`: exactly the tournaments whose fixed-point tournament `kappa` is circulant on `F` with a
   constant colour carry a second class of `H`-conjugates (the `p`-orbit `F`), each such type has `m_T = 2`, and
   there are `8 t_p^2 t_q` such `N`-orbits. (`t_3 = t_5 = 1`, `t_7 = 2`, `t_11 = 4`, `t_13 = 6`; `Z_p` is a CI-group.
   The census sees only `p = 3, 5`, where `t_p = 1`; audit A computed the subfamily at `(7,11)` and `(7,13)`,
   `k = 7`: `E_T(7) = 64 = 4 . 4 . 4` and `96 = 4 . 4 . 6`.)
4. **Witnesses and the first defect (PROVED).**
   (a) For `k = q - p` (`q > p`), the complete bipartite graph `K_(q,q)` with parts `Q` and `P cup F` is `H`-fixed,
   `N`-even and sign-odd, so `W(q - p) >= 1`. If `p < q < 2p` are primes, then `E_T(q - p) = 0` by 3(b) and
   `tau_T(H) - tau_E(H) = W + E_even + E_wit - E_T >= 1`: the refined equinumerosity HYP-9340 (a) fails for `H` of
   type `(p, q, 1^(q-p))` at degree `n = 2q` (`(3,5)` at `10`, the known C15 failure; `(5,7)` at `14`, `(7,11)` at `22`,
   `(11,13)` at `26`, `(13,17)` at `34`, ... an infinite family of mixed-prime failures with an explicit obstruction).
   (b) For `p = 3`, `q >= 7`, `k = 3`: the `2 t_q` graphs with `P cup F` carrying `K_(3,3)` between `P` and `F`, `Q`
   carrying a multiplier-even circulant `C_Q`, and `Q` joined to all of `P cup F` or to none, are pairwise non-isomorphic
   witnesses, so `W(3) >= 2 t_q`, while `E_T(3) = 4 t_q` by 3(c) (`t_3 = 1`).
5. **Census (FINITE-EXACT; 62 rows, thirteen families).** `tau_T, tau_E, W, E_T, E_even, E_wit` (tournament side first; `R = 2 t_p t_q b_k(4)`):

   | family | k = 0, 1, 2, ... | defect | W | E_T | E_even | E_wit |
   |---|---|---|---|---|---|---|
   | (3,5), t = 1,1 | tau_T 2, 8, 32, 172, 1376, 17564, 373804; tau_E 2, 8, 31, 168, 1360, 17480, 373132 | 0, 0, 1, 4, 16, **84**, 672 | 0, 0, 1, 4, 14, 60, 334 | 0, 0, 0, 4, 32, 260, 2772 | 0, 0, 0, 4, 34, 282, 3102 | 0, 0, 0, 0, 0, 2, 8 |
   | (3,7), t = 1,2 | tau_T 4, 16, 64, 344, 2752, 35136, 747672; tau_E 4, 16, 64, 344, 2751, 35132, 747656 | 0, 0, 0, 0, 1, 4, 16 | 0, 0, 0, 4, 17, 82, 482 | 0, 0, 0, 8, 64, 512, 5480 | 0, 0, 0, 4, 48, 434, 5014 | 0 |
   | (5,7), t = 1,2 | tau_T 4, 16, 64, 352, 2816, 35640, 753088; tau_E 4, 16, 63, 348, 2800, 35552, 752384 | 0, 0, 1, 4, 16, 88, 704 | 0, 0, 1, 2, 6, 24, 112 | 0, 0, 0, 0, 0, 8, 64 | 0, 0, 0, 2, 10, 72, 656 | 0 |
   | (3,11), t = 1,4 | tau_T 8, 32, 128, 688, 5504 | 0, 0, 0, 0, 0 | 0, 0, 0, 8, 32 | 0, 0, 0, 16, 128 | 0, 0, 0, 8, 96 | 0 |
   | (3,13), t = 1,6 | tau_T 12, 48, 192, 1032 | 0, 0, 0, 0 | 0, 0, 0, 12 | 0, 0, 0, 24 | 0, 0, 0, 12 | 0 |
   | (7,13), t = 2,6 | tau_T 24, 96, 384, 2112 | 0, 0, 0, 0 | 0 | 0 | 0 | 0 |
   | (7,11), t = 2,4 | tau_T 16, 64, 256, 1408, 11264; tau_E 16, 64, 256, 1408, 11263 | 0, 0, 0, 0, 1 | 0, 0, 0, 0, 1 | 0 | 0 | 0 |
   | (11,13), t = 4,6 | tau_T 48, 192, 768; tau_E 48, 192, 767 | 0, 0, 1 | 0, 0, 1 | 0 | 0 | 0 |
   | (9,5), t = 4,1 | tau_T 6, 24, 96, 528, 4224; tau_E 6, 24, 96, 528, 4223 | 0, 0, 0, 0, 1 | 2, 8, 32, 176, 1409 | 2, 8, 32, 176, 1408 | 0 | 0 |
   | (9,7), t = 4,2 | tau_T 12, 48, 192, 1056; tau_E 12, 48, 191, 1052 | 0, 0, 1, 4 | 4, 16, 65, 354 | 4, 16, 64, 352 | 0, 0, 0, 2 | 0 |
   | (9,11), t = 4,4 | tau_T 24, 96, 384, 2112; tau_E 24, 96, 383, 2108 | 0, 0, 1, 4 | 8, 32, 129, 706 | 8, 32, 128, 704 | 0, 0, 0, 2 | 0 |
   | (9,13), t = 4,6 | tau_T 36, 144, 576, 3168, 25344; tau_E 36, 144, 576, 3168, 25343 | 0, 0, 0, 0, 1 | 12, 48, 192, 1056, 8449 | 12, 48, 192, 1056, 8448 | 0 | 0 |
   | (3,25), t = 1,208 | tau_T 410, 1640 = tau_E | 0, 0 | 4, 16 | 6, 24 | 0 | 2, 8 |

   HYP-9344 (2) pre-registered `defect(5) = 88` (against the alternative `64`) for C15 `(3,5,1^5)` at `n = 13`: the
   value is `84`, with `W = 60`, `E_T = 260`, `E_even = 282` and the first nonzero `E_wit = 2`; `defect(6) = 672`.
   In every family the first nonzero defect is at `k = |p - q|` and equals `1`, followed by `4` and then `16`
   where computed; `tau_T >= tau_E` in all 62 rows; THM-4721 Theorem R and THM-4728's identity hold in all 62 rows.
   The rows `(7,11) k = 4` and `(11,13) k = 2` are the next two members of the family of statement 4(a) and show
   exactly one witness and no excess, as statement 7 proves.
6. **Prime powers (FINITE-EXACT).** For `p = 9` the four dilation classes of circulant 9-tournaments fall into three
   isomorphism types (`Z_9` is not a CI-group), which gives `E_T(k) = 2 t_q b_k(4)` already at `k = 0` and an equal
   number of witnesses `W(k) = E_T(k)` (sign-odd circulant graphs on `Z_9` that are multiplier-even), cancelling exactly
   for `k < |p - q|`; the first defect is again `1` at `k = |p - q|` (`k = 4` for `(9,5)`, `(9,13)`; `k = 2` for `(9,7)`,
   `(9,11)`), and `(3,25)` shows the same non-CI cancellation with witness multiplicities `E_wit = 2, 8`.
7. **First-defect law for primes `p < q < 2p` (PROVED; added after audit A, audited separately as audit B).**
   Let `p < q` be primes and `H` of type `(p, q, 1^k)`.
   (a) If `k < q - p`, every `N`-even `H`-fixed graph is sign-even and carries a single `Aut`-class of `S_n`-conjugates of
   `H`; so `W = E_even = E_wit = 0`, and with 3(b) (`k < q - p < p`) also `E_T = 0`: `tau_T(H) = tau_E(H)`.
   (b) If `k = q - p`, the only sign-odd `N`-even `H`-fixed graph type is `K_(q,q)` (parts `Q` and `P cup F`), every
   `N`-even type carries a single `Aut`-class of conjugates of `H`, and `tau_T(H) - tau_E(H) = 1`.
   Hence for `p < q < 2p` the defect vanishes for `k < q - p` and equals `1` at `k = q - p`: HYP-9463 (1) holds in the
   range `q < 2p`, and the least degree at which the pair `(p, q)` violates the refined equinumerosity HYP-9340 (a) is
   exactly `n = 2q` (for `(3, 5)` this is the Codex "least failure degree ten").

## Proofs

**1.** (a) An element of `S_n` normalising `<h>` permutes the `<h>`-orbits preserving their sizes; since `p != q` it
fixes `P`, `Q`, `F` setwise, and its restrictions normalise `<c_P>` and `<c_Q>`, i.e. lie in the holomorphs.
Conversely for `(alpha, beta, pi)` with `alpha c_P alpha^(-1) = c_P^a`, `beta c_Q beta^(-1) = c_Q^b`, CRT gives `j` with
`j = a (mod p)`, `j = b (mod q)` and then `(alpha, beta, pi) h (alpha, beta, pi)^(-1) = h^j`. The pair-orbit list is direct:
inside `P` the translate of `{i, i+d}` by `j` equals it only for `j = 0` (`p` odd), so each distance class is one orbit of
size `p`; a pair `{u, v}` with `u in P`, `v in Q` has orbit `{(c_P^j u, c_Q^j v)}` of size `lcm(p, q) = pq`; `{f, u}` has
orbit size `p` or `q`; pairs inside `F` are fixed.
(b) An `H`-fixed tournament orients each pair-orbit consistently: a distance class `{i, i+d}` as `i -> i+d` for all `i`
(`d in S_P`) or the reverse (`-d in S_P`), the cross orbit as `P -> Q` or `Q -> P`, `{f} x P` as `P -> f` or `f -> P`,
`{f} x Q` likewise; pairs inside `F` freely. Translations of `Z_p`, `Z_q` act trivially on this datum; a dilation by `a`
replaces `S_P` by `aS_P`; `pi in Sym(F)` permutes `kappa`; the colour of `f` is preserved by `alpha, beta` because they
preserve `P` and `Q`. So the orbit set is the stated product and `R(k) = t_p t_q 2 b_k(4)`.
(c) Likewise for graphs. `(alpha, beta, pi) in Aut X` iff `a C_P = C_P`, `b C_Q = C_Q`, `pi in Aut(gamma)`. For the sign
use the vertex order `P < Q < F` (the statement is order-free). The automorphism `alpha` (identity outside `P`) can
reverse only edges inside `P`: an edge `{u, w}` with `u in P`, `w notin P` has `u < w` and `alpha(u) < w`. Similarly
`beta` reverses only edges inside `Q` (`w in P` precedes both `u` and `beta(u)`; `w in F` follows both), and `pi` only
edges inside `F`. Translations have odd order, hence sign `+1`. Since `sgn_X` is a homomorphism the sign of a product
is the product of the three block signs, and the `N`-even criterion follows. The orbit count on `N`-even graphs is the
product count with the two circulant factors restricted to multiplier-even classes and `gamma` to colour-sign-even
types; THM-4721 Theorem R (`Gamma = N`) equates it with `R(k)`. The two separate identities are Theorem R for the
cyclic group `<c_m>` on `m` points (`N = Hol(Z_m)`: `N`-orbits on circulant tournaments = `N`-orbits on `N`-even
circulant graphs, i.e. `t_m = e_m`) and Theorem R for the trivial group with `Gamma` the Young subgroup
`S_(k_1) x ... x S_(k_4)` (summing over compositions `k_1 + ... + k_4 = k` gives `b_k(4) = ec_k(4)`).

**2.** Suppose `O cap O'` is nonempty and `O != O'`; put `Delta = O cup O'`, `m = |Delta| <= r + s - 1`, and
`K = <g, g'> <= G` (odd order). `K` fixes `Omega \ Delta` pointwise and is transitive on `Delta` (the two orbits `O`, `O'`
of `g`, `g'` meet). *Primitivity.* Let `B` be a system of blocks of size `b`, `1 < b < m`. For the `r`-cycle `g`: a block
moved by `g` lies in a `g`-orbit of `r` disjoint blocks, each meeting `O` (a block disjoint from `O` is fixed pointwise by
`g`), so each contains exactly one point of `O`; a block fixed by `g` meets `O` in `empty` or `O` (`g` is transitive on
`O`, `r` prime). Hence either (i) `O` lies inside one `g`-invariant block, or (ii) every block meets `O` in at most one
point. The same dichotomy holds for `g'` and `O'`. If (ii) holds for both, a block containing a point of `O cap O'` is a
singleton, contradicting `b > 1`. If (i) holds for `g` with `O <= B_0` and (ii) for `g'`, then `B_0 cap O'` is a single
point, so `|O cap O'| = 1` and every other block meeting `O'` is disjoint from `O`, hence a singleton: contradiction.
If (i) holds for both, the blocks containing `O` and `O'` meet (in `O cap O'`), so coincide and contain `Delta`:
`b = m`, contradiction. So `K` is primitive on `Delta`.
*Sizes.* If `O' <= O` then `Delta = O`, `m = r` is prime; by Burnside's theorem (Wielandt 11.7) a transitive group of
prime degree is 2-transitive or a subgroup of `AGL(1, r)`; `K` has odd order, so it is not 2-transitive (`r(r-1)` would
divide `|K|`), hence `K <= AGL(1, r)`, in which every non-identity element fixes at most one point; but `g'` fixes the
`r - s >= 2` points of `O \ O'` (both lengths are odd primes). Symmetrically `O` is not inside `O'`. So `t := m - r >= 1`
and `m - s >= 1`. If `t >= 3`, Jordan's theorem (a primitive group of degree `m` containing a cycle of prime length
`r <= m - 3` contains `A_m`; Wielandt 13.9) gives `A_m <= K`, of even order: impossible. If `t = 1`, `m = r + 1` is
even and `m` divides `|K|` (transitivity): impossible. If `t = 2`, `K` is a primitive group of odd, hence (Feit-Thompson)
solvable, order; by Galois's theorem it has a regular normal elementary abelian subgroup, `m = l^e` for a prime `l`
(odd, as `m = r + 2` is odd), and `K <= AGL(e, l)`; after conjugation `g` fixes `0`, so `g in GL(e, l)` and its fixed-point
set is the subspace `ker(g - 1)`, of size `l^j`; but `g` fixes exactly `2` points of `Delta` and `l^j = 2` is impossible for
odd `l`. (Alternatively: by Jones's classification of primitive groups containing a cycle, every such group other than
`C_l <= K <= AGL(1, l)` has even order, and in `AGL(1, l)` non-identity elements fix at most one point.) Hence
`O cap O' = empty` or `O = O'`, and `O = O'` forces `r = |O| = |O'| = s`. ∎

**3.** (a) Let `H' = <h'> <= Aut T` with `h' = c'_P c'_Q`, `c'_P` a `p`-cycle with support `P'`, `c'_Q` a `q`-cycle with
support `Q'`. Both are powers of `h'` (`c'_P = (h')^(q u)` with `q u = 1 (mod p)`), so they lie in `Aut T`, as do `c_P, c_Q`.
`Aut T` has odd order (an involution would contain a transposition reversing an arc). By statement 2, `P'` equals or is
disjoint from `P`, and is disjoint from `Q` (`p != q`); so `P' = P` or `P' <= F`; likewise `Q' = Q` or `Q' <= F`. If
`P' = P` and `Q' = Q`: the `p`-cycles of `Aut T` with support `P` generate, restricted to `P`, a transitive group of
prime degree and odd order, hence a subgroup of `AGL(1, p)` whose unique Sylow `p`-subgroup is `<c_P|_P>`; so
`c'_P in <c_P>`, similarly `c'_Q in <c_Q>`, and `h' = c_P^a c_Q^b = h^j` with `j = a (mod p)`, `j = b (mod q)`: `H' = H`.
(b) For `k < p` neither `P' <= F` nor `Q' <= F` is possible, so every conjugate of `H` in `Aut T` equals `H`: `m_T(H) = 1`.
(c) For `k = p`, the only alternative is `P' = F` (and `Q' = Q`). As in (a), all `p`-cycles of `Aut T` with support `F`
generate one cyclic group `<c_F>`, so there is at most one further class, `H' = <c_F c_Q>`, and `m_T(H) <= 2`.
`c_F in Aut T` iff `T[F]` is a circulant tournament for `c_F` and every vertex outside `F` beats all of `F` or loses to all
of `F`, i.e. iff `kappa` is circulant with constant colour `(x, y)`. In that case `H` and `H'` are not `Aut T`-conjugate:
a conjugating automorphism would map `P` onto `F` while all arcs between `P` and `F` point the same way (constant `x`).
The datum of `T` relative to `H'` is `(cross' = y, colour' = (1 - x, cross), S_Q, kappa' = T[P])`, which differs from the
datum relative to `H` in the `x`-coordinate, so the two `N`-orbits are distinct and each such type has exactly two.
The `N`-orbits with circulant constant-colour `kappa` number `t_p (S_P) x t_q (S_Q) x 2 (cross) x 4 (colour) x t_p`
(circulant tournaments on `p` points up to isomorphism; `Z_p` is a CI-group, so these are the dilation classes), i.e.
`8 t_p^2 t_q`, giving `4 t_p^2 t_q` types with `m_T = 2` and `E_T(p) = 4 t_p^2 t_q`.

**4.** (a) `X = K_(q,q)` with parts `A = P cup F` (`|A| = q`) and `B = Q` is `H`-fixed (`H` preserves `P`, `Q`, `F`). Every
element of `N` preserves `A` and `B`, so `N <= Aut X`; for `g in N` no edge `{a, b}` (`a in A`, `b in B`) can be swapped
by a power of `g`, so by the orbit form of the sign (the parity of order reversals over a `<g>`-orbit of pairs is the
indicator that some power of `g` swaps the pair) `sgn_X(g) = 1`: `X` is `N`-even. The side swap `sigma` (any involution
exchanging `A` and `B`) is an automorphism whose `<sigma>`-orbits on edges are the `q` swapped edges `{a, sigma a}` and
unflipped pairs, so `sgn_X(sigma) = (-1)^q = -1`: `X` is sign-odd, hence a witness. With `p < q < 2p` we have
`q - p < p`, so 3(b) gives `E_T(q - p) = 0` and THM-4728's identity gives the inequality.
(b) The configurations are `H`-fixed and `N`-even by 1(c) (`C_P` empty, `gamma` = the three fixed points independent
with constant colour `(1, y)`, `y = cross`; its colour-preserving automorphisms are permutations of three pairwise
non-adjacent twins, hence even). The swap `P <-> F` is an automorphism of odd sign: in the orbit form its three
self-paired edges `{i, sigma i}` are the only flipped orbits (the direct count in the order `P < Q < F` is `9` or
`9 + 6q` reversed edges, odd either way). The type determines `cross` (whether `Q` has neighbours outside `Q`) and the
isomorphism type of `X[Q] = (Z_q, C_Q)` (a CI-class, `q` prime), so the `2 t_q` graphs are pairwise non-isomorphic. ∎

**5, 6.** Computation (`cpq_defect_census.py`): the `N`-orbit representatives of statement 1 are built explicitly
(circulant classes by dilation orbits, coloured types by nauty canonical forms with the colour-class shape, 4-coloured
objects on `F` up to isomorphism), each full object is canonically labelled by nauty, and for graphs the automorphism
generators give the sign; `tau_T, tau_E, W, E_T, E_even, E_wit` are read off the type multiplicities. The brute-force
mode enumerates all `2^c` `H`-fixed objects, computes types directly and the `N`-action by closure. Audit A re-derived
the rows with two independent implementations (a definition-level brute force and an own orbit method with VF2
cross-checks of the nauty coincidences). ∎

**7.** Throughout, `X` is an `N`-even `H`-fixed graph with datum `(C_P, C_Q, cross, gamma)` and `k <= q - p < p`;
"module" means a homogeneous set of `X` (every outside vertex adjacent to all or to none of it), and an *independent
module* is a module inducing the empty graph, i.e. a set of pairwise non-adjacent twins; `Sym(I) <= Aut X` for every
independent module `I`, and the union of two overlapping independent modules is an independent module (a vertex of
one outside the other is non-adjacent to the intersection, hence to all of the other). Three facts are used.
*(i) Prime circulants.* A circulant graph on `Z_r`, `r` prime, other than the empty and complete graphs, is connected
and co-connected and vertex-transitive of prime order, so its automorphism group is a transitive group of prime
degree that is not 2-transitive, hence lies in `AGL(1, r)` (Burnside), and the graph is prime: the maximal strong
modules of its modular decomposition are permuted transitively, so they are singletons, and a prime root with
singleton children leaves only the trivial modules (Gallai). In particular no vertex of such a circulant is a twin of
any other vertex of `X` (its neighbourhood inside the circulant is a nonempty proper subset, while an outside vertex is
uniformly attached). The complete circulant `K_r` is never multiplier-even (the dilation by a primitive root has cycle
type `(1, r-1)`, an odd permutation, and on `K_r` the reversal sign is the permutation sign), so `C_P` and `C_Q` are
empty or prime circulants.
*(ii) Two overlapping cycles.* If `c, c' in Aut X` are cycles of the same prime length `r` with distinct overlapping
supports `O, O'`, then `K = <c, c'>` is transitive and primitive on `Delta = O cup O'` (the block argument in the proof
of statement 2 uses no parity) and 2-transitive: for `|Delta| >= r + 3` by Jordan (`A_Delta <= K`); for `|Delta| = r + 1`
because a point stabiliser contains the `r`-cycle; for `|Delta| = r + 2`, with `c` fixing `a, b`: `K_a >= <c>` is
transitive on `Delta \ {a, b}`, and if `K_a` fixed `b`, then `K_a <= K_b` with `|K_a| = |K_b|` gives `K_a = K_b` and
`Fix(K_a) = {a, b}`, and an element carrying `a` to `b` normalises `K_a` and must carry `b` into `Fix(K_a)`, so `{a, b}`
would be a block, against primitivity; hence `K_a` is transitive on `Delta \ {a}`. Since `K` fixes `V \ Delta`
pointwise, `Delta` is then a module inducing the empty or the complete graph.
*(iii) Signs across invariant sets.* If `g in Aut X` preserves each part of a partition of `V`, then `sgn_X(g)` is the
product of the signs of the restrictions of `g` to the induced subgraphs on the parts: an edge between two parts is
never swapped by a power of `g`, so by the orbit form of the sign it contributes nothing.

*Step 1: `Q` lies in an `Aut X`-invariant module `M` that is independent or equals `Q`.* For `g in Aut X`,
`|g(Q)| = q > p + k = |V \ Q|`, so `g(Q)` meets `Q`. If `g(Q) != Q`, (ii) for `c_Q` and `g c_Q g^(-1)` makes
`Q cup g(Q)` a module inducing the empty or complete graph; complete is impossible by (i) (`X[Q]` would be complete),
so `Q cup g(Q)` is an independent module, which by (i) cannot exist when `X[Q]` is a prime circulant. Hence: if `X[Q]`
is a prime circulant, `g(Q) = Q` for every `g`, and we put `M = Q`; if `X[Q]` is empty, let `M` be the union of all
independent modules containing `Q`, an independent module, and `g(M) = M` for every `g` (`g(M)` is an independent
module meeting `M`, since `g(Q)` meets `Q`, so `g(M) cup M <= M`).
*Step 2: `P` lies in an `Aut X`-invariant module `D`, either inside `M` or disjoint from it.* If `C_P` is a prime
circulant, put `D = P`: by (i) no vertex of `P` is a twin of any vertex, so `P cap M = empty` when `M` is
independent (and trivially when `M = Q`), and for `g in Aut X`, `g(P) <= V \ M <= P cup F` with `|F| < p` gives
`g(P) cap P != empty`; `g(P) != P` would make `P cup g(P)` a module inducing the empty or complete graph by (ii),
impossible since it contains the prime circulant `X[P]`; so `g(P) = P`. If `C_P` is empty, all vertices of `P` are twins
of one another (common neighbourhood `(Q if cross) cup {f : x_f = 1}`), so either `P <= M` (when some, hence every,
vertex of `P` is a twin of the vertices of `Q`, which requires `X[Q]` empty) or `P cap M = empty`; in the latter case
let `D` be the union of all independent modules containing `P`: `D cap M = empty` (otherwise `D cup M` is an
independent module containing `P` and `Q`, forcing `P <= M`), and `g(D) = D` for all `g` by the argument of Step 1
(`g(P)` meets `P` because `g(P) <= V \ M` and `|F| < p`; then `g(P) = P` or `P cup g(P)` is an independent module, so
`g(P) <= D` and `g(D) cup D <= D`). In the case `P <= M` put `D = M`.
*Step 3: every automorphism is even.* Put `F' = F \ (M cup D)`. Each `g in Aut X` preserves `M`, `D` and `F'`, so by
(iii) `sgn_X(g) = sgn_(X[M])(g|M) sgn_(X[D])(g|D) sgn_(X[F'])(g|F')` (with the middle factor omitted when `D = M`).
The first two factors are `+1`: an independent set has sign `+1` under every permutation, and when `M = Q` with `X[Q]`
a prime circulant (or `D = P` with `X[P]` one), `g|Q in Aut X[Q] <= Hol(Z_q)` by (i) is a translation (odd order, sign
`+1`) times a dilation in `M(C_Q)` (sign `+1` by multiplier-evenness). For the last factor, extend `g|F'` by the identity
on `F cap (M cup D)` to a permutation `g^` of `F`. It is a colour-preserving automorphism of `(F, gamma)`: for
`f in F'`, `x_f = [f ~ D]` and `y_f = [f ~ M]` (an outside vertex is uniformly attached to a module), both preserved
because `g` preserves `D` and `M`; a vertex of `F cap M` is a twin of the vertices of `Q`, so its neighbours in `F'` are
the `f` with `y_f = 1`, a set `g` preserves, and its colour is the constant `(cross, 0)`; a vertex of `F cap D` (when
`D != M`) is a twin of the vertices of `P`, with neighbours `{f in F' : x_f = 1}` and colour `(0, cross)`; the edges
between `F cap M` and `F cap D` are all present or all absent according to `cross`; there are no edges inside
`F cap M` or inside `F cap D`; and `g` preserves the edges inside `F'`. By 1(c), `gamma` is colour-sign-even, so
`sgn_gamma(g^) = +1`, and by (iii) applied to `F = (F cap M) cup (F cap D) cup F'`, `sgn_gamma(g^) = sgn_(X[F'])(g|F')`.
Hence `sgn_X(g) = +1` for every automorphism: `X` is sign-even, so `W = E_wit = 0`.
*Step 4: one class of conjugates.* Let `H' = <c' c''> <= Aut X` be `S_n`-conjugate to `H`, with `p`-orbit `P'` and
`q`-orbit `Q'`; we find `a in Aut X` with `a H' a^(-1) = H`. As in Step 1, `Q'` meets `Q`. If `X[Q]` is a prime
circulant, then `Q' = Q` (otherwise `Q cup Q'` would be an independent module), and `c''|Q in Aut X[Q] <= AGL(1, q)`
lies in its unique Sylow `q`-subgroup `<c_Q|Q>`: `c'' in <c_Q>`. If `X[Q]` is empty, `Q cup Q' <= M` and some element
of `Sym(M) <= Aut X` carries `Q'` to `Q`; conjugating by it we may assume `Q' = Q`, and then `Sym(Q) <= Aut X` (fixing
`V \ Q` pointwise, hence `P'`) conjugates `c''` into `<c_Q>`. Now `P' <= V \ Q = P cup F` and `|F| < p`, so `P'` meets
`P`. If `C_P` is a prime circulant, `P' = P` (as in Step 2) and `c'|P in <c_P|P>` by the Sylow argument in `AGL(1, p)`.
If `C_P` is empty and `P cap M = empty`: `P' = P`, or `P cup P'` is an independent module by (ii), contained in `D`;
`Sym(D) <= Aut X` fixes `Q` pointwise (`D cap M = empty`) and carries `P'` to `P`, and then `Sym(P)` conjugates `c'`
into `<c_P>`. If `C_P` is empty and `P <= M`: `c'` preserves `M` (Step 1) and `P'` meets `M` (`|P'| = p > |F \ M|`),
so the `<c'>`-orbit `P'` of a point of `M` lies in `M`; `Sym(M) <= Aut X` carries `(P', Q)` to `(P, Q)` (both pairs of
disjoint subsets of `M` of sizes `p, q`) and then `Sym(P) x Sym(Q) <= Sym(M)` conjugates the two cycles into
`<c_P>`, `<c_Q>`. In every case `a H' a^(-1) = <c_P^i c_Q^j> = H` (CRT): `m_X(H) = 1`, so `E_even = 0`. With
`E_T = 0` (3(b)), THM-4728's identity gives (a).
*Step 5: `k = q - p`.* Now `|V \ Q| = q`, and `g(Q) cap Q = empty` is possible only with `g(Q) = P cup F`. Then
`X[P cup F] ~= X[Q]` is a circulant on `q` vertices containing the module `P` of size `p`, `1 < p < q`, so by (i) `X[Q]`
is empty and `P cup F` is independent; `Q` and `P cup F` are modules (the latter as the orbit of the `q`-cycle
`g c_Q g^(-1)`), so between them `X` is complete or empty: `X = K_(q,q)` or `X` is empty. For every other `N`-even `X`
each automorphism satisfies `g(Q) cap Q != empty`, and Steps 1-3 apply verbatim (they used `k < p` and that
overlap): `X` is sign-even. The empty graph is sign-even, so the sign-odd `N`-even types are exactly `K_(q,q)`
(`W = 1`; sign-odd and `N`-even by 4(a)). Step 4 applies to every `N`-even `X` other than `K_(q,q)` and the empty
graph (a conjugate with `Q' cap Q = empty` has `Q' = P cup F` and forces, by the argument just given with `c''` in place
of `g c_Q g^(-1)`, `X in {K_(q,q), empty}`), the empty graph has `Sym(V) <= Aut X`, and in `K_(q,q)` a `q`-cycle
preserves the bipartition (an odd power of a part-swapping automorphism swaps the parts), so `Q'` is a part, `P'` lies in
the other, and `Sym(part) x Sym(part) <= Aut X` conjugates `H'` to `H`: `m_X(H) = 1` for all `N`-even `X`, i.e.
`E_even = E_wit = 0`, while `E_T = 0` by 3(b). Hence `tau_T - tau_E = W = 1`. ∎

## Scope and boundary

* Statements 3 and 4(a) need `p`, `q` prime (the lemma is about prime-length cycles; a 9-cycle's cube is a product of
  three 3-cycles). For prime powers the census shows additional coincidences of a different origin (non-CI pairs) that
  cancel against additional witnesses; no proof is offered there.
* Statement 7 settles the first-defect law for primes `p < q < 2p`: `0` below `k = q - p` and exactly `1` there. For
  `q > 2p` (`(3,7)`, `(3,11)`, `(3,13)`, `(7,13)`) and for prime powers (`(9,q)`, `(3,25)`) the census shows the same law,
  but there the first defect `1` is a net value of many witnesses and coincidences (`(3,7) k = 4`: `17 + 48 - 64`;
  `(9,5) k = 4`: `1409 - 1408`) and nothing is proved beyond `E_T = 0` for `k < p` (primes); the complete cancellation
  at `k = p = 3` (`W(3) = E_even(3) = 2 t_q` against `E_T(3) = 4 t_q`) is proved as lower bounds plus census.
* Statement 7's Step 2 uses `|F| < p` and Step 1 uses `|V \ Q| < q`, i.e. `k < min(p, q - p) = q - p`; the proof does
  not extend to `q - p >= p`, where `g(Q)` can miss `Q` without `X` being `K_(q,q)` or empty.
* The lemma uses Feit-Thompson only in the case `m = r + 2`; Jones 2014 gives the same case without it. The sanity
  check over all tournaments of orders 3-8 (838 prime-length cycles in automorphism groups; `lemma_check.py`,
  reproduced by audit A's `audit_b_lemma_search.py`) found no overlapping supports at all.
* Nothing here concerns `H` with two orbits of the same prime size, or non-cyclic `H`.

## Reading

The owner asked to prove HYP-9340 beyond Sylow in terms of odd primes and their powers. For a cyclic symmetry with
two prime orbits the obstruction is the *distance* `|p - q|` between the primes: the larger orbit `Q` can be twinned
with the smaller orbit padded by `q - p` fixed points, and the resulting `K_(q,q)` has the odd side swap that the
normaliser cannot see. Twinning the smaller orbit with `p` fixed points (`K_(p,p)`) also produces witnesses, but these
are exactly cancelled by the tournament and even-graph coincidences that the same `p` fixed points create (a second
`p`-cycle in the automorphism group), which the disjoint-orbit lemma locates precisely. The p-group theorem of
HYP-9340 (no defect for a single prime) is the degenerate case of the same bookkeeping, and prime powers add
non-CI coincidence/witness pairs that cancel in the same way. For primes `p < q < 2p` this is now a theorem (statement 7): nothing happens before `k = q - p`, and at `k = q - p`
exactly `K_(q,q)` happens. The value of the first defect, `1`, then `4`, `16`, ..., follows the 4-coloured tournament
counts `b_j(4)` corrected by the tournament coincidences `E_T(j)`: HYP-9463.
