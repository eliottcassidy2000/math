# Fermat tournaments, {-3}-self-duality, the few-cycle table and 31 | 32 | 33: which coincidences have a mechanism

**Session:** kind-pasteur-2026-11-04-S10 (third owner prompt of the session). Hypothesis: HYP-9464 (pre-registered in
3fc61b533c, confirmed at order 10).

Canon from the follow-up: THM-4874, the uniquely Hamiltonian deletion law. It proves HYP-9464 (B) and its hc = 5
extension, and was confirmed to order 15 by HYP-9467.

**Predecessors, credited throughout.**
- THM-4868 (this session): non-critical runs of uniquely Hamiltonian tournaments, and the rooted count
  R(m) = F(2m - 5) of THM-4868 (F). Garey 1972 counts the unrooted classes, F(2m - 6).
- opus-S75 (THM-4869, HYP-9461, MISTAKE-709; note
  `arc_connectivity_route_and_sole_exit_free_few_cycle_tournaments_20261104.md`). The pasted block of this prompt is
  word for word S75's prompt, and S75 answered it first:
  - Theorem A: the arc-connectivity shortcut is PROVED, and the four order-9 witnesses are classified.
  - THM-4869: u_3(n) = F(2n - 9), u_4(n) = 2F(2n - 11).
  - Its section 4 types the unit ball, -7/4 and -29/16, and the microcosm-MRB thread.

  This note does not repeat those answers; it adds to them.
- opus-S76 (THM-4872, HYP-9465; note `unit_start_pm2_port_pairs_symmetry_split_mrb_20261104.md`) worked
  concurrently on a sibling prompt. Its results used here:
  - the Z_3 / asymmetric split of the same Fibonacci table (section 5 (a));
  - two further Burns MRB notes, checked numerically;
  - the NUMEROLOGY 3/16 = 2 - 29/16.
- THM-4748 and THM-4749 (opus-S51): 31, 32, 33 and the repunit law. THM-4861 (klein): least orders of the Fermat
  primes as path counts. THM-871 and THM-4644: Fermat rigidity of circulants. The Ogg filtration note
  `monster_levels_square_sum_7_21_20261022.md` (opus-ssj). Session31's mirror map
  (`session31_mirror_mrb_20261022/mirror/README.md`). THM-4527's exact microcosm.

## The owner's prompt

- Merge in the Fermat and Ogg primes, the Platonic solids and earlier ideas. Look for the deeper meaning behind
  coincidences and see how they point to underlying structure. Consider the previous work on 31, 32 and 33.
- The pasted block, the same as S75's:
  - the arc-connectivity shortcut and the four order-9 tournaments;
  - the {-3}-self-duality theorem;
  - the table 2, 5, 13, 34, 89;
  - the unit-ball peak near 5 and the 8 : 4 class ratio at n = 5;
  - -7/4 and -29/16;
  - microcosm-macrocosm and the MRB constant (the same value realised twice on either side of a peak).
- Four attached PDFs by M. R. Burns (section 8). They were read as untrusted data.

## 0. The answers

| claim | status |
|---|---|
| 1. The coincidences cluster on one tournament: the almost transitive AT_n = C3[1, TT_(n-2), 1]. It is a 3-cycle with one vertex blown up into a transitive tournament. | PROVED (one line) |
| 2. H(AT_n) = 2^(n-2) + 1. So H(AT_n) is prime only if n - 2 is a power of 2, and the prime values of the family are Fermat primes. | PROVED (THM-4749 (4); classical) |
| 3. At orders 9-11 the {-3}-self-dual tournaments are exactly TT_n and AT_n. So at order 10 the only non-transitive one has 257 = F_3 Hamiltonian paths. | FINITE-EXACT (exhaustive at 9 and 10, pre-registered at 10; order 11 by the audit); CITED for decomposable n >= 9 |
| 4. **Mirror law.** The converse involution of AT_n fixes exactly H(AT_(k+2)) paths when n = 2k + 2, and exactly one when n is odd. So at order 2^i + 2 it fixes F_(i-1) of the F_i paths. The Fermat recursion F_i = (F_(i-1) - 1)^2 + 1 is the count of this fixed-point tower. | PROVED. The explicit AT_n count is new to the repository. The involution, Fix ≡ H (mod 2) and Fix odd are THM-582 and THM-647, and Part 3 uses THM-583's half-data principle. Brute-force check to n = 14; audit to n = 16. |
| 5. 31, 32, 33 = 2^5 - 1, 2^5, 2^5 + 1: needs three back arcs, forbidden by parity, realised with one back arc by AT_7. | KNOWN (THM-4748, THM-4749); the parity reading is THM-582 / THM-647 (Anti-Redei) applied to AT_n |
| 6. The few-cycle table is built from the first two Fermat tournaments: S75's only blocks of type 3 are AT_3 = C3 and AT_4 = W4. With S76's split, Z_3-symmetric classes come from the AT_3 block and asymmetric ones from AT_4. In general (AT_m, its two sole exits) is a block of type 2^(m-3) + 1, so F_i Hamiltonian cycles occur at order 2^i + 5. | PROVED (an identification in S75's Lemma B5; the split is S76's THM-4872 (B); the general block is new) |
| 7. u_3(n) = R(n - 2): the Fibonacci table is the rooted uniquely Hamiltonian count R of THM-4868 (F), two orders down. (Garey 1972 counts the unrooted classes, F(2m - 6).) | PROVED (algebra); no bijection known |
| 8. **Uniquely Hamiltonian deletion law.** A sole-exit-free tournament with exactly 3 Hamiltonian cycles has exactly 3 uniquely Hamiltonian deletions; with exactly 4 cycles (n >= 7), exactly 2. | PROVED in its range (n >= 6 for hc = 3, n >= 7 for hc = 4) by THM-4874: the core block decides hc(T - v) for every v. FINITE-EXACT at orders 6-15 (10-15 pre-registered as HYP-9467) |
| 9. Fermat ∩ Ogg = {3, 5, 17}, and the Platonic face sizes {3, 4, 5} contain the Fermat primes 3, 5. | NUMEROLOGY: both are forced by size |
| 10. The "three 9s" (the self-duality threshold; the first order with route-(b) witnesses; hc = |Aut| = 9). | ANALOGY (S75), sharpened by a module-type dichotomy: transitive modules only (FINITE-EXACT at orders 7 and 8, CITED from 9); no C3 module (PROVED at 7 and 8, CITED from 9) |
| 11. The MRB mirror and the microcosm. The x^(1/x) mirror pairs values across the peak at e. The AT_n mirror pairs the binomial row across its peak, and its fixed part is the half-size AT. | the two laws PROVED; the parallel is ANALOGY |

## 1. The hub: the almost transitive tournament AT_n

AT_n has vertices 1 < 2 < ... < n. Every arc points forward except n -> 1 (Busch's "nearly transitive", Moon's NT_n,
THM-4748/4749's NT_n). It is **not** the bineutral A_n of THM-4868; the two agree only at n = 4 (S75, MISTAKE-709).

**Facts (PROVED).**
1. **Substitution form.** The interior I = {2, ..., n-1} is a transitive module: 1 beats all of it and all of it
   beats n. The quotient on {1, I, n} is 1 -> I -> n -> 1, a 3-cycle. So AT_n = C3[1, TT_(n-2), 1]. In particular
   AT_3 = C3 and AT_4 = W4, the unique strong 4-tournament.
2. **Hereditary self-converse.** An induced subtournament of AT_n is AT_k if it contains both 1 and n, and TT_k
   otherwise. Both are self-converse (via i -> k + 1 - i). So AT_n is {-k}-self-dual for every k, and so is TT_n.

   This is KNOWN. Reid and Thomassen (1976, Monatsh. Math. 81, 291-304) also prove the converse: a strongly self-dual
   tournament (every induced subtournament self-dual) with at least 8 vertices is transitive or almost transitive.
   The citation is as quoted in arXiv:1204.2513 section 1.2, where the audit read it; the original was not accessed.
3. **Path count.** A Hamiltonian path either avoids n -> 1, and is then 1 2 ... n, or it is A n 1 B, with A and B the
   increasing listings of a split I = A ⊔ B. So H(AT_n) = 2^(n-2) + 1.

   More generally, reversing one arc of span d >= 2 in TT_n gives 2^(d-1) + 1 (Moon 1972; THM-250). Span 1 gives
   TT_n again, with H = 1. The single-back-arc values are 3, 5, 9, 17, 33, 65, 129, 257, ...
4. **Uniquely Hamiltonian** (n >= 4). AT_n has exactly one Hamiltonian cycle, a non-critical run of length n - 2,
   3 cut arcs, and two sole exits n - 1, n (THM-4868 (C), (D); THM-4868 assumes n >= 4). AT_3 = C3 has three sole
   exits and no non-critical vertex.
5. **Fermat.** 2^m + 1 is prime only if m is a power of 2 (classical). So H(AT_n) is prime iff n = 2^i + 2 and F_i is
   prime: n = 3, 4, 6, 10, 18 give 3, 5, 17, 257, 65537. F_5 to F_32 are known to be composite, so no other n <= 2^32 + 2
   gives a prime.

**What is new here: the {-3}-self-duality side (HYP-9464 (A)).**
- The cited theorem (Achour-Boudabbous-Boussairi, arXiv:1204.2513, Theorem 4) says: a decomposable tournament on at
  least 9 vertices is {-3}-self-dual iff it is transitive or almost transitive.
- Exhaustively, the {-3}-self-dual classes number 56, 10, 18, 2, 2 at orders 6, 7, 8, 9, 10. At orders 9 and 10 they
  are exactly TT_n and AT_n, with no indecomposable exception.
- Order 10 was pre-registered and run afterwards. All 9,733,056 classes were checked, and networkx identified the two
  survivors.
- So **at order 10 the only non-transitive {-3}-self-dual tournament has F_3 = 257 Hamiltonian paths**. By the
  cited theorem, at order 18 the only decomposable non-transitive one has F_4 = 65537.
- The audit's complete pruned search at order 11 finds only TT_11 and AT_11. The search was validated by recovering
  10 and 18 classes at orders 7 and 8. So no indecomposable {-3}-self-dual tournament exists at orders 9-11; from
  order 12 on this is OPEN.
- Compare THM-4861 (klein): 257 is first reached at order 8 (thirty classes), and 65537 at order 11. The self-dual
  realisation sits further out, at 10 and 18. THM-4861 (c) also shows that symmetry of the other kind (vertex
  transitivity) never gives a prime count except for C3.

**Mechanism.** The cited theorem forces the skeleton (TT or AT). Fact 3 then forces the count 2^(n-2) + 1, and Fact 5
the Fermat condition. Nothing beyond these three steps is needed. The Fermat prime at order 10 is a consequence, not a
coincidence.

## 2. The mirror law (PROVED; an explicit AT_n instance of THM-582, THM-583 and THM-647)

The involution below, the congruence H ≡ Fix (mod 2) and the oddness of Fix are canon:
- THM-582 (mac-mini): palindromic Hamiltonian paths of a self-converse tournament;
- THM-647 (monad): the Anti-Redei theorem;
- THM-024: every self-converse tournament has an involutive anti-automorphism.

Determining a palindromic path by half its data is THM-583's principle. THM-583 states it for an anti-automorphism with
one fixed vertex; for even n, sigma has none. What is new here is the exact fixed-point count for AT_n and its
identification with the half-size AT.

sigma(i) = n + 1 - i is an isomorphism from AT_n onto its converse. It acts on Hamiltonian paths as an involution:

    iota(P) = reverse(sigma(P)),    iota(A n 1 B) = sigma(B) n 1 sigma(A),    iota(1 2 ... n) = 1 2 ... n.

**Theorem (mirror law).** Put m = n - 2.
1. iota maps the paths with |A| = j onto the paths with |A| = m - j. The back-arc paths form the binomial row
   C(m, j), and iota is its mirror about the peak j = m/2.
2. Fix(iota) = 1 + 2^(m/2) if m is even, and 1 if m is odd.
3. Hence Fix(iota on AT_(2k+2)) = H(AT_(k+2)). At n = 2^i + 2 (i >= 1) the involution fixes F_(i-1) of the F_i
   paths. The Fermat recursion F_i - 1 = (F_(i-1) - 1)^2 is "all subsets of I = pairs of half-subsets", and the fixed
   paths are the diagonal.

*Proof.*
- **Part 1.** The second display: sigma reverses the order of I, so sigma(B) listed increasingly is reverse(sigma(B)).
- **Part 2.** iota fixes A n 1 B iff A = sigma(I - A), i.e. A contains exactly one vertex of each pair {j, n + 1 - j}
  in I.
  - For odd m the middle vertex is its own pair, which is impossible, so only the transitive path is fixed.
  - For even m there are 2^(m/2) such A.
- **Part 3.** A fixed A is determined by A ∩ {2, ..., k + 1}, an arbitrary subset of a k-set. That is the same data as
  a back-arc path of AT_(k+2), and the transitive path matches the transitive path. QED.

`at_mirror_fixed_paths_20261104.py` enumerates every path by DFS, without the formula, and checks parts 1-3 for
3 <= n <= 14 (`.out`): Fix = 1, 3, 1, 5, 1, 9, 1, 17, 1, 33, 1, 65 for n = 3..14.

**Readings.**
- **Parity.** iota is an involution, so H(AT_n) ≡ Fix(iota) (mod 2), and Fix is odd. This is THM-582 (and THM-647,
  Anti-Redei) for the anti-automorphism sigma. Here Fix is computed exactly, so its oddness is visible directly. For
  odd m the only fixed path is the transitive one: 33 = H(AT_7) = 1 + 2 x 16 mirror pairs.
- **Exact microcosm.** The fixed part of the order-(2k + 2) object is a full copy of the order-(k + 2) object. This is
  the same shape as THM-4527's exact microcosm S(4A + 1) = S(A) for the Syracuse map. That parallel is ANALOGY; both
  laws are PROVED.
- **MRB.** The owner's "same value realised twice on either side of the peak" is literal here. Each binomial count
  C(m, j) is realised at j and at m - j, the converse symmetry exchanges them, and the peak carries the fixed set. The
  x^(1/x) mirror (session31; the PDFs' r_n) exchanges n > e with r_n < e. Its rational integer points are only
  2 <-> 4, and its asymmetry ("small changes near 1 are large changes near infinity") has no counterpart in the
  symmetric binomial row. The parallel is ANALOGY.

  One exact contact: for x > 0 with x != 1, (x, x^2) is an x^(1/x) mirror pair iff x^(x^2) = (x^2)^x iff x = 2
  (x = 1 gives the degenerate pair (1, 1)). So the integer mirror pair
  2 <-> 4 is the first squaring step of the tower H(AT_(2^i + 2)) - 1 = 2, 4, 16, 256, 65536. This is PROVED, but it
  explains nothing further: NUMEROLOGY.

## 3. 31 | 32 | 33

The repository's previous work (opus-S51, THM-4748 and THM-4749; audited, MISTAKE-657):
- **33 = 2^5 + 1 = H(AT_7) = H(NT_7).** It needs one back arc. N(n, 33) has degree w(33) = 2 in n, from 33 = 3 x 11
  (two strong factors).
- **32** is impossible: H is odd (Redei).
- **31 = 2^5 - 1** needs at least 3 back arcs (THM-4748 (3), the seven-mod-eight lines). The minimum is attained by
  the "Fibonacci fan": an apex losing to a, b, a + b. The name refers to a + b, not to Fibonacci numbers. 31 is prime
  and a strong value, so w(31) = 1 and N(n, 31)/N(n, 33) -> 0. Exact counts (THM-4749 (6)): 10/23 at n = 8, 18/49 at
  n = 10, 34/125 at n = 14.

What this session adds:
- **The Fermat side is the self-dual side.** 2^k + 1 is realised by AT_(k+2): hereditarily self-converse, uniquely
  Hamiltonian, one back arc.
- **The Mersenne side.** 2^k - 1 needs at least three back arcs for k >= 5 (THM-4748 (3)). Three suffice, via the
  Fibonacci fan, which is an apex construction (THM-4748 (4)).
- **Why.** The "+1" is the transitive path, which a single reversed arc can only add to. Statement 2 of THM-4748 (the
  seven-mod-eight lines) is the precise form of "two reversed arcs cannot reach 2^k - 1 for k >= 5". Two arcs do
  reach 15 = 2^4 - 1 and 23 = 3 x 2^3 - 1.
- **The mirror view of the parity.** H(AT_7) = 33 = 1 + 2 x 16: the single fixed path plus 16 mirror pairs, so 32
  would need an even fixed count. This restates canon:
  - every self-converse tournament has an involutive anti-automorphism (THM-024);
  - H ≡ Fix (mod 2) for each (THM-582);
  - Fix is odd (THM-647).

  For all tournaments the oddness is Redei's theorem.
- **33 recurs** as the mirror-fixed count of AT_12 (part 3: Fix(AT_12) = H(AT_7)). It sits one level down the
  microcosm tower. This is PROVED, and illustrates the mechanism rather than adding a new one.
- The Fibonacci fan and the few-cycle Fibonacci table (section 5) share only a name. NUMEROLOGY.

## 4. Fermat primes, Ogg primes, Platonic solids

| link | mechanism | status |
|---|---|---|
| Fermat primes = prime path counts of single-back-arc tournaments | 2^m + 1 prime ⟹ m = 2^i, with section 1, Fact 3 | PROVED |
| Fermat primes = primes p all of whose circulant tournaments have Aut = Z_p exactly (rotation-rigid) | (Z/p)^× is a 2-group, and tournament automorphism groups have odd order | PROVED (THM-871, THM-4644, Akin 2018) |
| A vertex-transitive tournament has a prime count only if it is C3 | |Aut| divides H (THM-4747) | PROVED (THM-4861 (c)) |
| The two Fermat appearances (2^(n-2) + 1 paths; (Z/p)^× a 2-group) are the same condition 2^m + 1 prime | same arithmetic, different objects. They meet only at C3 = AT_3, the circulant on Z_3 with F_0 = 3 paths. For n >= 4, AT_n has trivial automorphism group and is not vertex-transitive. Its score classes are {1, 2}, {n - 1, n} and singletons, and an automorphism cannot swap the two ends of an arc. No circulant has a Fermat prime >= 5 as its count (THM-4861 (c)). | PROVED as stated; no shared tournament mechanism beyond C3 |
| Fermat ∩ Ogg's supersingular primes = {3, 5, 17}. On the level filtration M, B, Fi22 of `monster_levels_square_sum_7_21_20261022.md` (lists from Aricheta 2019, CITED), 3 and 5 are at level 6 (Fi22) and 17 at level 2 (B); 257 and 65537 lie outside | every prime <= 31 is supersingular, and 257 > 71 | NUMEROLOGY (size) |
| Platonic face sizes {3, 4, 5}; the odd primes among them, 3 and 5, are Fermat | 1/p + 1/q > 1/2 forces p <= 5, and every odd prime below 7 is Fermat. The classical link is Gauss-Wantzel: each face is a constructible polygon, as is any p-gon with p <= 6. | NUMEROLOGY (small numbers) |
| Planar groups {Z_n, D_n, A_4, S_4, A_5} and McKay | treated in this session's THM-4870 note (the five 5's and {7, 21}). Tournament automorphism groups have odd order, so none of A_4, S_4, A_5 is one. | see THM-4870 |

**The structural statement behind the table.** "Fermat" enters tournament theory through a prime being one more than
a power of 2 (p - 1 = 2^m).
- In the circulant case that power of 2 is the number p - 1 of multipliers.
- In the AT case it is the number 2^(n-2) = H(AT_n) - 1 of subsets of the interior. Both are PROVED. The Ogg and Platonic overlaps are explained by size alone.

## 5. The few-cycle table 2, 5, 13, 34, 89

The counts are u_3(n) for sole-exit-free tournaments with exactly three Hamiltonian cycles. S75 proved
u_3(n) = F(2n - 9) (THM-4869).

**(a) Built from the first two Fermat tournaments (PROVED).** In S75's structure theory a sole-exit-free tournament
with hc = 3 has:
- a terminal core: a transitive chain followed by a block of type 3 (Lemma B4);
- hc(T) = s_l(T - D) x e_u(D) with e_u(D) >= 2 (Lemmas B2, B3). So e_u(D) = 3 and s_l = 1: either T - D is the
  single vertex l, or the contraction Q = T/D is uniquely Hamiltonian.

Lemma B5 lists the type-3 blocks as (C3, all three vertices) and (W4, its two sole exits). By section 1, Fact 1,
these are AT_3 and AT_4, the Fermat tournaments for F_0 = 3 and F_1 = 5.

**With opus-S76's symmetry split (THM-4872 (B), proved concurrently, credited).** The C3-block classes are exactly the
classes with Aut = Z_3, and number F(2n - 10). The W4-block classes are asymmetric, and number F(2n - 11). In the
language of this note, the Fibonacci step F(2n - 9) = F(2n - 10) + F(2n - 11) sorts the table by its Fermat atom:
- AT_3: a module, rotated by Z_3, the larger part;
- AT_4: rigid, the smaller part.

The counts are also the partial-sum relation of S75's core counts. The C3 cores (2^(d-5) at order d >= 5) are the
partial sums of the W4 cores (1 at order 5, then 2^(d-6)), and sum_(j <= m) F(2j - 1) = F(2m). The proof is S76's;
only the AT reading is added here.

**More generally (PROVED).** For every m >= 4, (AT_m, {m - 1, m}) is a block of type 2^(m-3) + 1 = H(AT_(m-1)).
- **Type.** end_m(AT_m) = 1 (only the transitive path ends at m), and end_(m-1)(AT_m) = 2^(m-3) (the paths A m 1 B
  with m - 1 ∈ B).
- **The other block conditions.** Each vertex of {1, ..., m - 2} has out-degree >= 2 in AT_m. Let Z be non-empty.
  - If Z meets A, it would need no leaving arc, so Z = W by strongness; but then |Z ∩ A| = 2.
  - If Z misses A, every z ∈ Z beats m - 1 and every z ≠ 1 beats m, so at least |Z| arcs leave Z, and {1} alone
    sends m - 2 >= 2 arcs out.

  So no Z has exactly one leaving arc plus |Z ∩ A| = 1.
- **The smallest tournament on this block.** Take D = AT_m + u, with u beaten by m - 1 and m and beating the rest.
  Add l with only u -> l leaving D, and l beating all of AT_m. This T has order m + 2 and no sole exit.
- **Its cycle count.** Every Hamiltonian cycle enters l from u, leaves it into AT_m and returns to u from m - 1 or m.
  So hc(T) = end_(m-1) + end_m = 2^(m-3) + 1.
- **Fermat values.** At m = 2^i + 3 the count is F_i: 3, 5, 17, 257 Hamiltonian cycles at orders 6, 7, 9, 13. Path
  counts F_i at orders 2^i + 2 thus have cycle-count twins F_i at orders 2^i + 5.

`at_block_fermat_cycles_20261104.py` (+ `.out`) checks the block conditions over all subsets Z, the type, the absence
of a sole exit, and hc(T) by DP for 4 <= m <= 12.

**(b) Two orders down (PROVED).** u_3(n) = F(2n - 9) = F(2(n - 2) - 5) = R(n - 2).
- R(m) = F(2m - 5) counts the rooted pairs (T, sole exit) of uniquely Hamiltonian tournaments of order m (THM-4868
  (F)). Garey 1972 counts the unrooted classes, F(2m - 6). Equivalently, R(m) is the number of uniquely Hamiltonian
  classes of order m + 1 with two sole exits (confirmed by the audit at orders up to 14).
- So u_3(n) also equals the number of uniquely Hamiltonian classes of order n - 1 with two sole exits.
- In generating functions, U_3(x) = x^2 (R(x) - x^2 - x^3), with R(2) = F(-1) = 1 (a convention) and
  R(3) = F(1) = 1 (C3 rooted at a vertex).

The equality of counts is algebra. A bijection by deletion was probed (`u3_deletion_bijection_probe_20261104.py`,
`.out`): each T has exactly three uniquely Hamiltonian deletions T - v, classified by their numbers of sole exits and
sole entries. Neither "the deletion with two sole exits" nor "the deletion with two sole entries" is defined for
every T. At n = 10, 22 of the 89 classes have three deletions with one sole exit and one sole entry each. A
bijection is OPEN.

**(c) The uniquely Hamiltonian deletion law (HYP-9464 (B)).** The table counts sole-exit-free classes with hc = 3 and
hc = 4 by the number of vertices v such that T - v is uniquely Hamiltonian.

| hc | n = 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|
| 3 | 2 classes, all 3 | 5, all 3 | 13, all 3 | 34, all 3 | **89, all 3** |
| 4 | 3 classes: 2, 3, 4 | 4, all 2 | 10, all 2 | 26, all 2 | **68, all 2** |
| 5 | 4: 0, 2, 2, 6 | 4: 2 or 3 | 9: 2 or 3 | 23: 2 or 3 | not run |

Order 10 (bold) was pre-registered in 3fc61b533c and run afterwards. The independent audit
(`fermat_self_duality_20261104_audit_S10/aud_fewcycle.c`, every strong class with hc <= 5) reproduces orders 6-10 and
extends the law:

| hc | n = 11 | 12 | 13 | 14 |
|---|---|---|---|---|
| 3 | 233, all 3 | 610, all 3 | 1,597, all 3 | 4,181, all 3 |
| 4 | 178, all 2 | 466, all 2 | 1,220, all 2 | 3,194, all 2 |

**PROVED (THM-4874).** The deletion count is decided by the block (W, A) of the terminal core:
- T - v is uniquely Hamiltonian iff v ∈ W and sum over x in A - v of end_x(W - v) = 1;
- the C3 block gives its three vertices, and the W4 = AT_4 blocks give a, d, plus c for type 3.

The proof is the one-arc product law hc(T - v) = e_u(D - v) s_l(R - v), and the law is confirmed to order 15
(HYP-9467).

**Audit observation (not pre-registered).** For hc = 5 at orders 7-14, L(2n - 12) classes have exactly 2 uniquely
Hamiltonian deletions (3, 7, 18, 47, ..., 2,207), and F(2n - 13) have exactly 3 (1, 2, 5, 13, ..., 610). This matches
S75's u_5 = 4F(2n - 12) - 3F(2n - 14) = L(2n - 12) + F(2n - 13) term by term. It is now PROVED (THM-4874): the
three-deletion classes are exactly those whose core block is (AT_5, its two sole exits).

**(d) Disjointness.** For orders 7-10 no {-3}-self-dual tournament has hc = 3 or 4. Their hc values are:
- order 7: 0, 1, 8, 10, 13, 14, 17;
- order 8: 0, 1, 16, 22, 29, 32, 34, 37, 48, 49, 52;
- orders 9 and 10: 0 and 1. The two halves of the prompt meet only through AT_3 and AT_4 as
atoms, not as whole tournaments.

## 6. The three 9s

S75 typed the shared order 9 as ANALOGY: "order 9 is where a decomposable structure is first forced into an (almost)
transitive skeleton". The module data sharpen it into a dichotomy:
- **Self-dual side.** Every module of every decomposable {-3}-self-dual class is transitive.
  - At orders 7 and 8 this is an exhaustive check (FINITE-EXACT; `selfdual3_modules_20261104.py` over the census
    record, reproduced by the audit).
  - For every n >= 9 such a class is TT_n or AT_n (CITED Theorem 4). Their modules are intervals: TT_2, ..., TT_(n-1),
    respectively TT_2, ..., TT_(n-2).
- **Witness side.** Each of the four route-(b) witnesses has two C3 modules (v => C3 => x => y => C3' => v), so
  hc = 3 x 1 x 1 x 3 = 9 and |Aut| = 3 x 3 = 9 (S75 section 1; S74 audit; re-checked by the audit's
  `aud_witnesses.py`).

So the two 9s sit on opposite sides of one line: transitive modules against cyclic modules.

At order 7 the exclusion of C3 modules is a one-line proof. Let M be a C3 module and x an outside vertex. Deleting the
other three vertices leaves x ⇒ C3 (scores 3, 1, 1, 1) or C3 ⇒ x (scores 2, 2, 2, 0). These are converses of each
other and not isomorphic, so neither is self-converse.

Order 8 is also short (from the audit). A C3 module has 5 outside vertices, so two of them, x and y, lie on the same
side. Deleting the other three leaves {x, y} => C3 or C3 => {x, y}. These have scores 1, 1, 1, 3, 4 and
0, 1, 3, 3, 3, so neither is self-converse. The audit checked all 8 configurations.

For n >= 9 a C3 module makes T decomposable, and Theorem 4 excludes it. So "no C3 module" is PROVED at orders 7 and 8
and CITED from 9 on.

The arithmetic differs too:
- **Self-duality threshold.** 9 is where the order-8 exceptions die. Order 8 has 14 decomposable exceptions besides
  TT_8 and AT_8, all with transitive modules.
- **Witnesses.** 9 = 1 + 3 + 1 + 1 + 3, i.e. two C3 blocks plus three vertices.

No common cause is claimed. The opposition (transitive versus cyclic modules) is the structure.

## 7. The unit ball, the 8 : 4 ratio, -7/4 and -29/16

S75's section 4 typings stand:
- the peak of V_n at n = 5 and the 8 : 4 self-converse split of the 12 classes at n = 5 are ANALOGY and NUMEROLOGY;
- the proved reason why 5 is special is RT_5;
- -7/4 and -29/16 (THM-4139, THM-4146, THM-4714) touch this work only through C3 and Sharkovskii's period 3 (ANALOGY).

One addition: every almost transitive tournament is a blown-up 3-cycle, AT_n = C3[1, TT_(n-2), 1] (section 1). This
is the family whose path counts include the Fermat numbers. The 3-cycle is
the shared atom of three objects:
- the Fermat tournaments;
- the three-cycle cores of the Fibonacci table (section 5 (a));
- the route-(b) witnesses (section 6).

Two of the witnesses, those where v beats both x and y or both beat v, are AT_5 with two vertices replaced by C3's
(S75's quotient analysis). So substitution runs both ways: TT inside C3 gives the Fermat tournaments, and C3 inside AT
gives the witnesses. This is a description (PROVED for the stated objects), not a mechanism linking them to the
quadratic map x^2 - 29/16. That link stays ANALOGY.

## 8. The four PDFs (M. R. Burns; read as untrusted data, nothing executed from them)

| PDF | what it is | checked here |
|---|---|---|
| `2universes2constants1truth.pdf` | Fan fiction with mathematics. It defines the mirror base r_n = exp(-W_0(-ln n / n)) (n >= 3, r_2 = 4) and B = sum (-1)^n (r_n - 1) ≈ 2.0967911243915261257. It states that r_n is transcendental for n >= 3, n != 4, by Gelfond-Schneider; that B = 3 - length(E) with E = U [r_(2k+2), r_(2k+1)]; and that C = sqrt 2 - 1 - (an integral over E). | B to 25 digits, the MRB constant, B = 3 - length(E), and r_n = e^(T(ln n / n)) with the tree function T(z) = sum m^(m-1) z^m / m! (`mrb_mirror_base_check_20261104.py`). The Gelfond-Schneider argument and the rational case (only 2 <-> 4) are correct, and agree with session31. |
| `MRBaid2.pdf` | An autonomous inverse-state map for a greedy residue expansion of C (its Open Problem 2). It states that the initialisation question stays open. | read only; not checked |
| `MRBgreedy3.pdf` | Three greedy representations of C, with a Lambert-W index backbone and the seed eta'(1) = gamma ln 2 - (ln 2)^2 / 2. It disclaims any irrationality criterion. | the eta'(1) value is classical; the rest is not checked |
| `MRBprimeRecovery.pdf` | Finite partial sums containing n = 2 are irrational; prime radicals p^(1/p) are recoverable from an analytic deformation. It states that it does **not** prove C ∉ Q. | read only; not checked |

None contains tournament content. Their genuine contact with this work is the mirror shape of section 2.

opus-S76 checked two further Burns notes, "The MRB Constant via Lambert W" and "The Rational Term 3/16 in the MRB
Constant". It typed C_MRB = 3/16 + R with 3/16 = 2 - 29/16 as NUMEROLOGY, which touches the -29/16 thread of
section 7. Neither is re-checked here. One numerical
caution: mpmath's default `nsum` on length(E) is wrong from the third or fourth decimal, depending on precision.
- With it, 3 - length(E) comes out as 2.09732 (dps 15), 2.09693 (dps 30) or 2.09688 (dps 40).
- The true value is B = 2.0967911243... (audit `aud_mrb.out`).
- Euler-Maclaurin on the real-argument mirror gives B to 21 digits; the script uses it and notes the pitfall in a
  comment.
- The audit's own Euler-Boole and Cohen-Villegas-Zagier code gives B = 2.0967911243915261257250306750849. All 28
  digits printed in the PDF are correct.

## 9. Typing summary

- **PROVED:**
  - AT_n = C3[1, TT_(n-2), 1];
  - AT_n is hereditarily self-converse (KNOWN: Reid-Thomassen 1976);
  - H(AT_n) = 2^(n-2) + 1, and the Fermat prime criterion;
  - the mirror law (section 2), an explicit AT_n instance of THM-582, THM-583 and THM-647. The fixed-point count and
    its identification with H(AT_(k+2)) are the new part;
  - u_3(n) = R(n - 2);
  - the type-3 blocks are AT_3 and AT_4;
  - (AT_m, {m - 1, m}) is a block of type 2^(m-3) + 1, giving sole-exit-free tournaments with F_i Hamiltonian cycles
    at orders 2^i + 5;
  - no {-3}-self-dual tournament of order 7 or 8 has a C3 module;
  - for x > 0, x != 1, (x, x^2) is an x^(1/x) mirror pair iff x = 2.
- **FINITE-EXACT:**
  - the {-3}-self-dual classes at orders 6-11, with orders 9-11 exactly {TT_n, AT_n} (order 11 by the audit);
  - the uniquely Hamiltonian deletion law at orders 6-15. The law itself is now PROVED (THM-4874).
  - transitive modules only, at orders 7 and 8;
  - the PDF numerics.
- **CITED:**
  - Achour-Boudabbous-Boussairi Theorem 4;
  - Reid-Thomassen 1976 (as quoted there);
  - Moon 1972; Akin 2018;
  - Gauss-Wantzel; Gelfond-Schneider;
  - the compositeness of F_5 to F_32.
- **KNOWN:**
  - THM-4748, THM-4749, THM-4861, THM-871, THM-4644, THM-4869;
  - THM-4872 (B) (S76's symmetry split);
  - THM-024, THM-582, THM-583, THM-647.
- **OPEN:**
  - indecomposable {-3}-self-dual tournaments at orders >= 12;
  - a bijection for u_3(n) = R(n - 2).
- **ANALOGY / NUMEROLOGY:** sections 4, 6 and 7 as marked, and the MRB parallel.

## 10. Reproduction

All paths are in `04-computation/experiments/`.
- `selfdual3_fewcycle_20261104.c`: S73's census engine plus the modes `selfdual3`, `hcdump` and `autspec`. Build with
  `gcc -O3 -fopenmp`, then run `selfdual3 10` (38 s), `hcdump 10 3 2` and `hcdump 10 4 2`. Regenerate the `hcdump`
  files with these commands before running the readers below, or use the committed copies in
  `selfdual3_fewcycle_20261104_records/`.
- `uh_deletion_count_20261104.py`: reads the `hcdump` files.
- `selfdual3_modules_20261104.py` (+ `.out`): module types of the {-3}-self-dual classes (section 6).
- `selfdual3_networkx_check_20261104.py` (+ `.out`): independent networkx code. It checks that AT_n and TT_n are
  {-3}-self-dual for n = 9-11, that random tournaments of order 10 are not, and the identity of the two order-10
  classes.
- Records: `selfdual3_fewcycle_20261104_explore.out` (before the pre-registration) and `_confirm.out` (after). The raw
  outputs are in `selfdual3_fewcycle_20261104_records/`.
- `at_mirror_fixed_paths_20261104.py` (+ `.out`): the mirror law, n <= 14.
- `u3_deletion_bijection_probe_20261104.py` (+ `.out`): the deletion-bijection probe of section 5 (b), run on the
  `hcdump` records.
- `at_block_fermat_cycles_20261104.py` (+ `.out`): the AT_m blocks of section 5 (a), 4 <= m <= 12.
- `mrb_mirror_base_check_20261104.py` (+ `.out`): the PDF numerics.
- The audit: `fermat_self_duality_20261104_audit_S10/` (AUDIT.md and six programs with outputs).

## 11. Audit

Independent adversarial audit: **PASS WITH CORRECTIONS**, 1 major and 11 minor findings, all applied (MISTAKE-712).
The report is `04-computation/experiments/fermat_self_duality_20261104_audit_S10/AUDIT.md`.
- **The major finding.** The parity reading of the mirror (section 2) restates canon (THM-024, THM-582, THM-583,
  THM-647) and had been presented as new. It is now credited; only the AT_n fixed-point count is new.
- **The minor findings:**
  - the Fermat condition was stated backwards in section 4;
  - the x = 1 case of the mirror pair;
  - the module dichotomy's range and typing;
  - the order-8 proof;
  - four credits (Garey, Reid-Thomassen, S75/S74, the Ogg note);
  - an overstatement of THM-4748;
  - the nsum numbers;
  - uncommitted records;
  - domain qualifiers.

The audit reproduced every number with its own code: a labelled-tournament census to order 8, a pruned search to order
11, a few-cycle census to order 14, and the AT-family, witness and MRB checks. It extended (A) to order 11 and (B) to
order 14.
