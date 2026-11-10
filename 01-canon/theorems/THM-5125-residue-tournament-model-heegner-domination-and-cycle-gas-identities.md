---
id: THM-5125
title: "The quadratic-residue tournament model of the Hamiltonian-path spectrum. (1) On the primes congruent to 3 mod 4, p -> q iff p is a quadratic residue mod q defines a tournament R; it has the extension property, so it is the countable generic tournament and every finite tournament is induced on some set of primes (also inside the class 3 mod 8); H of the tournament on the prime factors of a squarefree n counts the residue chains of n, and Chan-Pak is the statement that these counts are exactly the odd numbers other than 7 and 21. (2) Heegner domination: a prime q = 3 mod 4 of class number one beats, in R, every prime p = 3 mod 4 below (q+1)/4. (3) The tournament on the seven Heegner primes 3, 7, 11, 19, 43, 67, 163 has 71 Hamiltonian paths; its prefixes have 1, 1, 3, 5, 9, 33, 71; 231 = 3*7*11, 4389 = 3*7*11*19 and 188727 = 3*7*11*19*43 are the least n with 3, 5 and 9 residue chains. (4) Cycle-gas identities: per(I+A), det(I+A), H and H mod 3 are four evaluations of the cycle-packing enumerator, and the cubic conditions nullity(A^2+A+I) > 0 and nullity(M^2+M+I) > 0 (Monsky convention) are parity balances of cycle packings graded mod 3."
status: >
  PROVED (opus-2026-11-10-S90): statement 1 (Dirichlet's theorem, the Chinese remainder theorem and quadratic
  reciprocity; the finite realisation inside the class 3 mod 8 is KNOWN in the repository,
  bsd_selmer_growth_20261013a.md section 3; the countable statement is the tournament analogue of the classical
  quadratic-residue model of the Rado graph, no priority claimed), statement 2 (from the classical norm-form
  argument for class number one), statement 4 (linear algebra over F_2 and F_4 plus the odd-cycle formula; the
  Monsky formula itself is CITED from the repository note, not reproved). FINITE-EXACT: statement 3 and the table
  of least realisers (redei_model.py), the identities of statement 4 on every labelled tournament of orders 3-6
  (cycle_gas.py). CITED: the list of imaginary quadratic fields of class number one (Heegner, Baker, Stark);
  Rabinowitsch's criterion. UNVERIFIED and not used: whether Redei's 1934 parity theorem was motivated by his work
  on class groups of quadratic fields (the sources consulted on 2026-11-10 do not say). The Chan-Pak conjecture
  remains OPEN; statement 1 is a reformulation, not progress on it.
session: opus-2026-11-10-S90
source: 05-knowledge/results/s90_rigidity_redei_20261110/README.md
scripts:
  - 04-computation/experiments/s90_rigidity_redei_20261110/redei_model.py (+ redei_model.out)
  - 04-computation/experiments/s90_rigidity_redei_20261110/cycle_gas.py (+ cycle_gas.out)
related:
  - 05-knowledge/results/bsd_selmer_growth_20261013a.md, bsd_component_recursion_20261013c.md (Monsky matrix, reciprocity tournament)
  - 01-canon/theorems/THM-4666-* (class number as Paley ranking excess), THM-4751-*, 05-knowledge/results/ramanujan_heegner7_trivial_cycle_20261006.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md, THM-4094-*, THM-4940-*, THM-5071-*
  - 05-knowledge/hypotheses/HYP-9657-* (mass profile and refuted bridges), 01-canon/theorems/THM-5124-*
  - Grinberg-Stanley, arXiv:2307.05569; L. Redei, Ein kombinatorischer Satz, Acta Sci. Math. (Szeged) 7 (1934) 39-43
---

# THM-5125 - the quadratic-residue tournament model

## Setting

For distinct primes `p, q` congruent to `3` modulo `4`, quadratic reciprocity gives `(p/q) = -(q/p)`. Write `p -> q`
when `p` is a quadratic residue modulo `q`. This is a tournament `R` on the set `P3` of primes `= 3 (mod 4)`. For a
squarefree `n` with all prime factors in `P3`, `T_n` is the tournament induced on its prime factors, and a *residue
chain* of `n` is an ordering `q_1, ..., q_k` of its prime factors with `q_i` a quadratic residue modulo `q_(i+1)` for
every `i`; the number of residue chains is `H(T_n)`, odd by Redei's theorem.

## Statements

1. **Universality.**
   (a) (*Extension property.*) For disjoint finite `U, V <= P3` there are infinitely many primes `r` in `P3`, which
   may be taken `= 3 (mod 8)`, with `r -> u` for all `u` in `U` and `v -> r` for all `v` in `V`.
   (b) Hence every finite tournament is isomorphic to `T_n` for infinitely many `n` (KNOWN inside the class
   `3 mod 8`: bsd_selmer_growth_20261013a.md), and `R` is isomorphic to the countable generic (random) tournament.
   (c) The Chan-Pak conjecture is equivalent to: the number of residue chains of a squarefree `n` with all prime
   factors `= 3 (mod 4)` takes every odd value except `7` and `21` (it never takes those two: THM-1370).
2. **Heegner domination.** Let `q` be a prime, `q = 3 (mod 4)`, such that `Q(sqrt(-q))` has class number one
   (`q` in `{3, 7, 11, 19, 43, 67, 163}`, CITED). Then `q -> p` for every prime `p = 3 (mod 4)` with `p < (q + 1)/4`.
   Thus in the Heegner tournament the larger prime wins unless the smaller one is at least `(q + 1)/4`.
3. **The Heegner tournament (FINITE-EXACT).** Arcs (`p -> q` iff `p` is a residue mod `q`):

       3 -> 11;  7 -> 3, 19;  11 -> 7, 19, 43;  19 -> 3, 67;  43 -> 3, 7, 19, 163;  67 -> 3, 7, 11, 43;
       163 -> 3, 7, 11, 19, 67.

   Scores `(1, 2, 3, 2, 4, 4, 5)`, eight triangles, `H = 71`. Prefixes `{3}, {3,7}, ..., {3,...,163}` have
   `H = 1, 1, 3, 5, 9, 33, 71`: `{3, 7, 11}` is cyclic and `{3, 7, 11, 19}` is the strong 4-tournament. The six
   Heegner primes `= 3 (mod 8)` (`n = 294441081`) have `H = 29` and Monsky nullity `4`. Its sub-tournament counts are
   `{1, 3, 5, 9, 13, 15, 17, 23, 29, 33, 71}`.
   Least realisers (certified least over all squarefree `n` with prime factors `= 3 mod 4`): `a(1) = 3`,
   `a(3) = 231`, `a(5) = 4389`, `a(9) = 188727 = 3*7*11*19*43`, `a(11) = 100947`, `a(13) = 206283`,
   `a(15) = 311619`, `a(17) = 48894657`, `a(19) = 8378601`, `a(23) = 4340721`; the table through `N = 45` is in
   `redei_model.out` (the entry for `N = 35` is least only inside the search box).
4. **Cycle-gas identities.** Let `A` be the adjacency matrix of a tournament `T` (`a_ij = 1` iff `i -> j`). A cycle
   packing is a set of pairwise disjoint directed cycles; write `ev(P)` for its number of even cycles and `|V(P)|`
   for the number of covered vertices.
   (a) `per(I + A) = #{cycle packings}`, `det(I + A) = sum_P (-1)^ev(P)`,
       `H(T) = sum over odd packings of 2^|P|`, and `H(T) = sum over odd packings of (-1)^|P|  (mod 3)`.
   (b) Over `F_2`, `det(x I + A) = sum_k pi_k x^(n-k)` with `pi_k` the parity of the number of cycle packings
       covering exactly `k` vertices. Consequently `A` has the eigenvalue `omega` in `F_4` (equivalently
       `nullity_F2(A^2 + A + I) > 0`) iff `P_0 = P_1 = P_2 (mod 2)`, where `P_j` is the number of cycle packings, the
       empty one included, with `|V(P)| = j (mod 3)`.
   (c) For `M = A + diag(row parities)` (the Monsky convention), `det(x I + M) = sum_P x^e(P) (x + 1)^o(P)` over
       `F_2`, where `e(P)`, `o(P)` count the vertices not covered by `P` with even, resp. odd, out-degree; and
       `nullity_F2(M^2 + M + I) > 0` iff `Q_0 = Q_1 = Q_2 (mod 2)`, `Q_j` the number of cycle packings with
       `e(P) - o(P) = j (mod 3)`. By the cited Monsky formula this nullity is the pure 2-Selmer dimension of
       `y^2 = x^3 - n^2 x` when `T = T_n` and all `p_i = 3 (mod 8)`. (The nullity does not depend on which of the two
       converse orientations is used: passing to the converse replaces `M` by `M^T` or `M^T + I`, and
       `x^2 + x + 1` is invariant under `x -> x + 1` in characteristic `2`.)

## Proofs

**1.** (a) For each `u` in `U` choose a nonzero residue class `c_u` mod `u` that is a square, and for each `v` in `V`
a class `c_v` mod `v` that is a non-square (possible since `v >= 3`). By the Chinese remainder theorem the conditions
`r = c_w (mod w)` for `w` in `U u V` and `r = 3 (mod 8)` define one residue class modulo `8 prod w`, coprime to the
modulus; by Dirichlet's theorem it contains infinitely many primes `r`. Then `(r/u) = 1`, i.e. `r -> u`, and
`(r/v) = -1`, so `(v/r) = 1` by reciprocity, i.e. `v -> r`.
(b) Realise a tournament on `[k]` vertex by vertex using (a). A countable tournament with the extension property
is the Fraisse limit of the finite tournaments, unique up to isomorphism by back-and-forth.
(c) By (b) the set of values `H(T_n)` is the set of all Hamiltonian-path counts of tournaments.

**2.** Let `K = Q(sqrt(-q))`, `q = 3 (mod 4)`, with ring of integers `Z[w]`, `w = (1 + sqrt(-q))/2`, and norm form
`N(x + y w) = x^2 + x y + ((q + 1)/4) y^2 = (x + y/2)^2 + q y^2/4`. If `y != 0` the norm is at least `q/4`, hence at
least `(q + 1)/4`. Let `p < (q + 1)/4` be a prime. If `p` were not inert there would be a prime ideal of norm `p`,
principal because the class number is one, generated by some `x + y w` of norm `p`; `y = 0` gives `p = x^2`,
impossible, and `y != 0` gives `p >= (q + 1)/4`, a contradiction. So `p` is inert: `(-q/p) = -1`. For `p = 3 (mod 4)`,
`(-1/p) = -1`, hence `(q/p) = 1`: `q` is a quadratic residue modulo `p`, that is `q -> p`.

**3.** Computation (`redei_model.py`): Legendre symbols by Euler's criterion, path counts by subset dynamic
programming; least realisers by exhaustive search over products below `3 * 10^9` of at most seven primes
`= 3 (mod 4)` below `600`, certified when the count first occurs at the order found (exhaustive order spectra,
orders at most `6`) and `n < 600 * (product of the k - 1 smallest primes = 3 mod 4)`, so that every smaller candidate
lies in the box.

**4.** (a) Expanding `per(I + A)` over permutations, a permutation contributes `1` iff each of its nontrivial cycles
is a directed cycle of `T` (there are no 2-cycles), i.e. iff it is a cycle packing; in `det(I + A)` it carries the
sign `prod over cycles (-1)^(length - 1) = (-1)^ev`. The third identity is the odd-cycle formula; the fourth is its
reduction modulo `3`, since `2 = -1`.
(b) Over `F_2` the coefficient of `x^(n-k)` in `det(x I + A)` is the sum of the `k x k` principal minors, each equal
modulo `2` to the permanent of the principal submatrix, the number of cycle covers of `T[S]`. Then
`det(omega I + A) = sum_k pi_k omega^(n-k)`; multiplying by `omega^(-n)` and using `omega^3 = 1`,
`omega^2 = omega + 1`, the sum `P_0 + P_1 omega^2 + P_2 omega = (P_0 + P_1) + (P_1 + P_2) omega` vanishes iff
`P_0 = P_1 = P_2`. Since `x^2 + x + 1 = (x + omega)(x + omega^2)` over `F_4` with conjugate kernels,
`nullity_F2(A^2 + A + I) = 2 dim_F4 ker(A + omega I)`.
(c) The same expansion with diagonal entries `x + d_v`, `d_v` the out-degree parity, gives the product
`x^e (x + 1)^o` over the uncovered vertices; at `x = omega`, `omega + 1 = omega^2`, so the term is
`omega^(e + 2o) = omega^(e - o)`, and the balance criterion follows as in (b). QED.

## Remarks

1. *What statement 1 buys.* Nothing toward completeness: by (b) the residue tournaments are all tournaments. Its
   use is a dictionary. Redei's parity theorem becomes: the number of residue chains of `n` is odd. The Monsky
   matrix of the same tournament computes a 2-Selmer rank (CITED), and statement 4 shows that this rank and `H` are
   two characters of one cycle gas: `H` grades cycles by parity with weight `2` (reversal acts on a `k`-cycle by
   `(-1)^(k-1)`, so even cycles cancel and odd cycles double in the Redei-Berge function), while the Selmer
   condition grades packings modulo `3` with coefficients in `F_2` (the cubic `x^2 + x + 1`, i.e. the action of an
   element of order `3` on the 2-torsion of the curve). No congruence between the two beyond these identities
   survives small orders (HYP-9657 (C)).
2. *What statement 2 buys.* It explains the shape of the Heegner tournament (10 of its 21 arcs are forced by
   domination: `19 -> 3`; `43 -> 3, 7`; `67 -> 3, 7, 11`; `163 -> 3, 7, 11, 19`; the other 11 are computed) and
   nothing more; `71`, `29` and the Monsky nullity `4` are data. The same fact (class number one
   forces all small primes to be inert, so `q` and `(q + 1)/4` are prime) is what makes Heegner numbers light for
   sum-of-two-products families (HYP-9657 (B)).
3. *Typed analogies* (no map claimed): the exclusion of `7` and `21` is a statement of local-global type with the
   only local condition at `2` (oddness); the countable object `R` plays the role of an adelic ambient space in
   which every finite configuration occurs; see the note for the ledger of analogies and for those refuted.
