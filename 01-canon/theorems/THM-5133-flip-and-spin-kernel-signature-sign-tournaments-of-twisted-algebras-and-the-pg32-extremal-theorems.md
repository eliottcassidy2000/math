---
id: THM-5133
title: "Flip and spin for Hamiltonian path counts. (A) Visit matrix: for a tournament K with two ports, H(K[A, B]) = sum N_ij P_i(A) P_j(B) with N_ij the number of walks of K visiting the ports i and j times and every other vertex once; with at most two visits the kernel is a 2 x 2 form G in (H, P_2), N_22 H = L(A) L'(B) + det(G) H(A) H(B), and it is elliptic, parabolic or hyperbolic by the sign of det G. (B) For a closing of THM-5130 with the two ports joined, det G = (offset)(slope): seeds are elliptic, anti-seeds hyperbolic, the sign flip of THM-5130 is a change of signature; first hyperbolic kernel at order 5, H = 3 Pi_2(A) Pi_2(B) - H(A) H(B); Gaussian kernel at order 4, H = Pi_2(A) Pi_2(B) + H(A) H(B). (C) Census of all two-port kernels of order at most 9; plus carriers d = 3, 21, 35 (known) and one minus carrier d = 31: 31 2^(a+b+1) - 2^(a+4) - 1 is a path count for all a, b >= 1. (D) Sign tournaments of twisted group algebras on Z_2^n: a line is cyclic iff it is a quaternion algebra and transitive, with the unit of square -1 in the middle, iff it is split; (e_a + e_b)(e_c +- e_d) = 0 iff a + b = c + d and the 4-set contains an even number of cyclic triangles. (E) Fano plane: the tournaments with every line complement odd are one orbit of 1024 under switchings and line coboundaries, 16 with cyclic lines (the octonions); the number of even line complements is never 1. (F) PG(3,2): no tournament has all 105 zero-sum 4-sets odd; at most 8 of the 15 planes are odd, in exactly 135 ways (15 + 120); the least number of even 4-sets is exactly 21; with all lines cyclic it is 42, attained by the sedenions. (G) The Cayley-Dickson sign tournament on 2^n - 1 vertices has all lines cyclic and is doubly regular for every n; H = 198,335,025 at n = 4. (H) Fermat modules: Pi_2 = (-1)^n H (mod 3) restricts them, and they exist at orders at most 9 only for levels 1, 2, 3, 7."
status: >
  PRE-AUDIT DRAFT (independent audit pending; statuses below are the author's).
  (A) PROVED (path surgery). The principle is KNOWN: Session193 profile_substitution/PROOF.md Theorem 1 (multi-port
  transfer), THM-3121 (walk-content substitution). Specific here: the explicit two-port visit matrix, the
  factorisation with det G and the trichotomy. Bilinear laws for single kernels are KNOWN: Session260
  QUASIPERFECT_NORMS.md eq. (1) (C_3, x x' + 2 h h'), Session261 BILINEAR_KERNEL.md (37 H = x x' + 329 H H'), s91
  README (one-arc joins). The four-term form for transitive modules is the module family of the S21/S22 packets.
  (B) PROVED from the closing law (THM-5130 (A), row h is Session216's) and multiplicativity of H and Pi_2 under
  the ordered join. The two identities are recounted on 120 random module pairs.
  (C) FINITE-EXACT (kernels.c orders 3-7, kernels_sat.c orders 8-9, validated against each other at orders 5-7; the
  order-8 signature counts agree with the independent table of modfam.c). Each carrier identity is PROVED by (A)
  once its visit matrix is computed. The plus carriers d = 3, 21 are KNOWN (multiplicand_saturation_chanpak_20261111,
  OFFSET_TRANSPORT.md: (16, 24), (4, 42)). New: the minus carrier d = 31.
  (D) PROVED (elementary). The algebraic criterion for two-term zero divisors is classical in the literature on
  Cayley-Dickson zero divisors (Moreno; de Marrais); the readings in terms of tournaments are specific here. The
  line rule is PROVED for line-associative algebras and FINITE-EXACT for all generalized Cayley-Dickson algebras of
  dimension 8, 16, 32.
  (E) FINITE-EXACT (exact linear algebra over F_2 and enumeration of all 2^21 tournaments). KNOWN: the 16
  alternative orientations and the H values 189 / 171 (collatz_mod6_20260922_paley_fano_octonion_design.md,
  S14-S16).
  (F) FINITE-EXACT (exhaustive exact linear algebra; five programs). KNOWN / CITED: the sedenions have 84 zero
  products (Moreno; de Marrais, arXiv math/0011260) and 8 octonion + 7 quasi-octonion subalgebras (Cawagas;
  Wilmot, arXiv:2505.11747). Specific here: the impossibility for every sign pattern, the bound 8 with its 135
  extremal sets, the exact minimum 21, the minimum 42 for cyclic lines.
  (G) PROVED by induction from the doubling formula (lines; Hadamard). No novelty claimed for the Hadamard
  property (elementary; not located in a brief web search on 2026-11-11). FINITE-EXACT: the values at n = 4, 5, 6;
  the counts 84 and 1,260 are CITED (Wilmot).
  (H) KNOWN congruence (THM-4756, step 2); FINITE-EXACT census (198,948 classes of order at most 9; one-apex words
  of length at most 22). Pre-registered predictions F1 and F2 FAILED.
  Chan-Pak remains OPEN. Pre-registration: PREREGISTRATION.md of the packet, pushed before the runs (b91157e81d);
  outcomes F1 FAILED, F2 FAILED, K1 CONFIRMED, K2 CONFIRMED, C1 CONFIRMED, C2 CONFIRMED (the count in C2 was
  already FINITE-EXACT in canon).
session: kind-pasteur-2026-11-11-S23
source: 05-knowledge/results/kps_s23_flip_spin_kernels_20261111/README.md
scripts:
  - 05-knowledge/results/kps_s23_flip_spin_kernels_20261111/scripts/kernels.c, kernels_sat.c, ker_stats.py, ker_census.py, ker_verify.py
  - 05-knowledge/results/kps_s23_flip_spin_kernels_20261111/scripts/cd_signs.py, cd_more.py, cd_lines_aut.py, cd_orbit.py, cd_extremal.py, cd_min.py, cd_f5_orbits.py
  - 05-knowledge/results/kps_s23_flip_spin_kernels_20261111/scripts/pairspec.c, classgen.c, e1_summary.py, words_f2.c
related:
  - 01-canon/theorems/THM-5130-sign-flip-towers-gadget-transfer-law-and-seeds-of-both-signs.md
  - 01-canon/theorems/THM-4756-forward-step-polynomial-mod-3-is-h-times-eulerian.md
  - 05-knowledge/hypotheses/HYP-9664-signatures-carriers-of-both-signs-and-extremal-sign-patterns.md
---

Notation. `H(T)` = number of Hamiltonian paths. `P_i(T)` = number of ordered `i`-tuples of disjoint nonempty paths
covering `T`; `P_1 = H`, `P_2 = Pi_2 - 2H`, `Pi_2` = ordered pairs of disjoint possibly empty paths covering `T`.
`K[A, B]` = the tournament `K` with the vertex `x` replaced by `A` and `y` by `B`. Pair `(B_0, O)`, state
`(r, M, ...)`, closing `C(B_0, O; X)`: THM-5130. Proofs: sections 2-4 of the source note.

## (A) Visit matrix and signature

For a tournament `K` of order `m` and ports `x != y` let `N_ij` be the number of walks in `K` that visit every vertex
other than `x, y` exactly once, `x` exactly `i` times and `y` exactly `j` times. Then `N_ij = 0` unless
`1 <= i, j <= m - 1`, and for all tournaments `A, B`

    H(K[A, B]) = sum_{i,j} N_ij P_i(A) P_j(B).

If no walk visits a port three times, `G = [[N_11, N_12], [N_21, N_22]]` and `h = H`, `pi = Pi_2`:

    H(K[A, B]) = al h h' + be pi h' + ga h pi' + de pi pi',
    de = N_22, be = N_21 - 2 N_22, ga = N_12 - 2 N_22, al = N_11 - 2 N_12 - 2 N_21 + 4 N_22, al de - be ga = det G,
    N_22 H(K[A, B]) = (N_22 P_2(A) + N_12 H(A)) (N_22 P_2(B) + N_21 H(B)) + det(G) H(A) H(B),
    H(K[TT_a, TT_b]) = al + be 2^a + ga 2^b + de 2^(a+b)      (a, b >= 1).

The kernel is called elliptic, parabolic, hyperbolic as `det G > 0, = 0, < 0`.

## (B) Signature and offset

For a pair `(B_0, O)` with state `(r, M)` the kernel `K = C(B_0, O; x => y)` has at most two visits,
`(al, be, ga, de) = (r, 0, 0, M)` and `det G = r M`:

    H(K[A, B]) = r H(A) H(B) + M Pi_2(A) Pi_2(B).

It is elliptic iff `r > 0` and hyperbolic iff `r < 0`. Order 4, `(r, M) = (1, 1)`:
`H = Pi_2(A) Pi_2(B) + H(A) H(B)`. Order 5, `(B_0, O) = (TT_2, both)`, `(r, M) = (-1, 3)`:
`H = 3 Pi_2(A) Pi_2(B) - H(A) H(B)`, `G = [[11, 6], [6, 3]]`.

## (C) Census and carriers

Ordered port pairs with at most two visits (join-like means `N_12 = N_21 = N_22 = 0`):

| order | all | at most two visits | join-like | elliptic | parabolic | hyperbolic |
|---|---|---|---|---|---|---|
| 3 | 12 | 12 | 6 | 6 | 0 | 0 |
| 4 | 48 | 46 | 12 | 22 | 12 | 0 |
| 5 | 240 | 160 | 26 | 76 | 54 | 4 |
| 6 | 1,680 | 630 | 60 | 290 | 218 | 62 |
| 7 | 19,152 | 2,892 | 162 | 1,242 | 946 | 542 |
| 8 | 385,280 | 16,962 | 602 | 6,732 | 5,008 | 4,620 |
| 9 | 13,790,592 | 148,596 | 3,804 | 53,804 | 39,306 | 51,682 |

Every odd `N > 1` other than `2^s + 1` is `N+(d, r, s) = 1 + 2^s + d 2^(r+s)` and every odd `N >= 1` is
`N-(d, r, s) = d 2^(r+s) - 2^s - 1`, each for unique `d` odd and `r, s >= 1`. A kernel with
`(al, be, ga, de) = (1, 2^u, 0, d 2^w)` gives `N+(d, r, s)` for all `s >= u + 1`, `r >= max(1, w - u + 1)`; with
`(-1, -2^u, 0, d 2^w)` it gives `N-(d, r, s)` for the same ranges. Among the kernels of order at most 8: plus
carriers exactly for `d = 3, 21` (KNOWN: the carriers (16, 24), (4, 42) of OFFSET_TRANSPORT.md) and a minus carrier exactly for `d = 31` (order 8, code 17082497, ports
(0, 7), `G = [[215, 124], [108, 62]]`):

    H(K[TT_a, TT_b]) = 31 2^(a+b+1) - 2^(a+4) - 1      for all a, b >= 1.

At order 9 the only further carrier is the plus carrier `d = 35` (KNOWN, the carrier (32, 35) of OFFSET_TRANSPORT.md);
so through order 9 the plus carriers are `d = 3, 21, 35`, the only minus carrier is `d = 31`, and there is none for
`d = 9` on either side.

## (D) Sign tournaments

Let `V` be the set of nonzero vectors of `Z_2^n`. A sign pattern `g(a, b) = -g(b, a) = +-1` on pairs of distinct
points is a tournament (`a -> b` iff `g(a, b) = 1`) and the table of the imaginary units of any algebra with basis
`1, e_a` and `e_a e_b = g(a, b) e_(a+b)` (`a != b`).

1. If `1, e_a, e_b, e_(a+b)` span an associative subalgebra, the line `{a, b, a+b}` is a cyclic triangle iff the
   three units have square `-1`; otherwise exactly one has square `-1`, the line is a transitive triangle and that
   unit is its middle vertex.
2. For distinct `a, b, c, d`: `(e_a + e_b)(e_c + s e_d) = 0` for some sign `s` iff `a + b = c + d` and the number of
   arcs from `{a, b}` to `{c, d}` is even. A zero-sum 4-set `Q` (`a + b + c + d = 0`) admits such a product iff it
   contains an even number of cyclic triangles (`Q` transitive or strong: `Q` is *even*), and then for exactly two
   of its three pairings; never if it contains exactly one (`Q` is *odd*).
3. The parities of all zero-sum 4-sets are invariant under switching at a vertex and under line coboundaries
   (reverse the arcs `{v, w}` with `v + w = f`); for `n = 3, 4` these moves span the kernel of the parity map
   (dimensions 10 and 25).

## (E) Fano plane

The tournaments on `PG(2, 2)` whose seven line complements are all odd are exactly 1024, one orbit of the Paley
tournament under the moves of (D3); 16 of them have all lines cyclic (the octonion tables), and the other 112
tournaments with all lines cyclic have exactly 6 even line complements. Over all `2^21` tournaments the number of
even line complements is 0, 2, 3, 4, 5, 6, 7 for 1024, 21504, 143360, 322560, 688128, 637952, 282624 of them;
it is never 1.

## (F) PG(3, 2)

15 points, 15 planes, 105 zero-sum 4-sets (seven per plane). A plane is odd if its seven 4-sets are odd.

1. No tournament on `PG(3, 2)` has all 105 zero-sum 4-sets odd.
2. At most 8 planes are odd in one tournament; the sets of 8 planes that are odd in some tournament are exactly 135:
   the 15 sets of all planes off a point, and the 120 sets formed by a plane and the seven planes through a point
   not on it.
3. In the second type the other seven planes have 6 even 4-sets each (42 in all); in the first type the total is
   at least 21, with equality for 64 violation patterns (3 in each of the seven planes).
4. The least number of even zero-sum 4-sets of a tournament on `PG(3, 2)` is 21.
5. If all 35 lines are cyclic, the number of even 4-sets is `6 (15 - number of odd planes) >= 42`, with equality iff
   8 planes are odd; the sedenion sign tournament attains it.

## (G) Cayley-Dickson sign tournaments

With the doubling `(a, b)(c, d) = (ac - conj(d) b, d a + b conj(c))` the sign table of order `2^(n+1)` is
`[[G, G^T], [G D, -G^T D]]`, `D = diag(1, -1, ..., -1)`. For every `n >= 2` the sign tournament `CD_n` on `2^n - 1`
vertices has all lines cyclic and is doubly regular. `CD_3` is the Paley tournament on 7 points; `CD_4` has
`H = 198,335,025`, 21 automorphisms, 8 odd planes (Paley) and 7 planes with 6 even 4-sets (`H = 171`); the numbers
of even 4-sets are 42, 630, 6510 for `n = 4, 5, 6` (`42 [n-1 choose 3]_2`). Split forms: the split quaternions give
the transitive triangle, the split octonions a tournament with `H = 45`.

## (H) Fermat modules

`Pi_2(X) = (-1)^n H(X) (mod 3)` (THM-4756). If `H(X) = 2^s + 1` and the odd part of `Pi_2(X)` divides 9, then
`Pi_2(X)` is a power of two for `s` even and has odd part 3 or 9 for `s` odd. Among all classes of order at most 9
such `X` exist only for `s = 1, 2, 3, 7`; among one-apex words of chain length at most 22 only for `s = 1, 2`.
For every tournament `X`, every pair with state `(r, M)` and every `t >= 0`, `r H(X) + M 2^t Pi_2(X)` is a path count.
