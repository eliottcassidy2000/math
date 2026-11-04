# Independent audit: THM-4869 / HYP-9461 (opus-2026-11-04-S75)

**Audited.** Worktree `C:/tmp/wt-s75`, branch `opus-2026-11-04-S75`, checkpoint `60620d918c`.
- Note: `05-knowledge/results/arc_connectivity_route_and_sole_exit_free_few_cycle_tournaments_20261104.md` (read fully).
- Canon: `01-canon/theorems/THM-4869-...md`.
- Pre-registration: `05-knowledge/hypotheses/HYP-9461-...md`.
- Programs: `fewcycle_sole_exit_free_20261104.c`, `fewcycle_core_lineage_20261104.c`, `fewcycle_blocks_20261104.c`, and their
  outputs. These were read, not reused. Every number below comes from the audit's own code in this folder.

## Verdict: PASS WITH CORRECTIONS

Theorems A, B1-B5 and C are correct. The arguments survive an adversarial reading, and the claimed order of proof has
no circularity. Every number in the note was reproduced by independent code:
- all labelled tournaments of orders 5-8, counted into classes by Burnside's lemma (no canonical form);
- an exhaustive class enumeration through order 10 (all 9,733,056 classes; no Moon theorem);
- the audit's own Moon-extension generator for orders 11-14;
- direct enumeration of cores and blocks through order 14.

There are no MAJOR findings. There are 13 MINOR findings:
- one false identification in the ANALOGY paragraph (the almost transitive tournament is not THM-4868's A_n);
- a repeat of MISTAKE-707 (the word "strong" dropped);
- a missing hypothesis in the converse of B3;
- a same-day credit (S74 had the Fraisse-Thomassen cut-vertex remark);
- two lineage errors;
- misquotes, paths and wording.

## 1. Claim-by-claim verdicts

| claim | verdict | evidence |
|---|---|---|
| A: hc = S+ * prod H(D_i) * S-, S+, S- >= 3, hc >= 9 | **correct** | proof checked; formula and both bounds hold on every (class, cut vertex) pair of every 2-arc-strong class with a cut vertex, at orders 7-10 (2,928 pairs at order 10) and at order 11 with hc <= 40 (3,195 pairs) |
| A: minimum 9 at orders 7-11 | **reproduced** | 1, 4, 9, 18, 36 classes with hc = 9 at orders 7-11; nothing below 9 |
| A': UH implies lambda = 1 (Fraisse-Thomassen) | **correct** | FT statement confirmed in two secondary sources (section 5); the cut-vertex half needs only S+, S- >= 2 |
| forced-arc structure; order-9 witnesses | **correct**; credit (m4) | exactly 4 classes at order 9, none at orders <= 8; all are v => C3 => x => y => C3' => v with forced arc x -> y, automorphism group of order 9, converse pattern as stated. At order 10 (new): 41 classes and 49 forced arcs, all with a cut vertex, of types junction (48) or D_1/D_m (1). There are no middle or at-v arcs, as the note's classification says |
| B1 uniqueness/disjointness | **correct** | submodularity and the D_1 u D_2 = V case checked; exactly one minimal set in every u-class at orders 6-14 |
| B2 product law | **correct** | 0 failures, orders 6-14 |
| B3 core properties, e_u >= 2 | **correct**; converse needs "D strong" (m3) | |
| B4 chain + block, 2^(d-w-2) | **correct** | the direct core census equals the formula at every d = 5-14 |
| B5 block list (k <= 5) | **correct** | direct enumeration of (W, A), hc(W) <= 2, orders 3-14 |
| B5 hc(W) = 2 Redei argument | **correct** | see section 3 |
| C induction P(n)/S(n), bound n <= 12 | **correct** | see section 3 |
| base P(6)..P(12) | **reproduced** | exhaustive through 10 and Moon extension at 11-14: 2-strong classes with hc <= 5 exist only at order 5 (RT_5) and order 6 (hc 4: 1 class, hc 5: 4 classes) |
| census completeness (Moon + hc(T - v) <= hc(T)) | **correct** | the session's K = 40 census has 1, 1, 6, 35, 353, 5,524, 35,355, 148,185 classes at orders 3-10, exactly the Moon-free exhaustive counts of strong classes with hc <= 40 |
| (C3), (C4) laws and generating functions | **correct** | algebra checked to n = 39 by exact series; census u_3, u_4, u_5 at n = 6-14 equal the laws, plus the 2-strong classes at n = 6 |
| min hc of 2-strong: 2, 4, 8, 13, 19, 30 (n = 5-10); > 40 at n = 11 | **reproduced** | |
| HYP-9461 P1-P10 | **reproduced**; blind range is 11-14 (m10) | |
| lineage (section 3.1) | **two errors** (m5, m6) | |
| section 4 inspiration facts | mostly correct; m1, m11, m12 | section 6 |

## 2. Findings

**MAJOR: none.**

**m1. The almost transitive tournament is not THM-4868's A_n.**
- Note lines 120-121 say: "an almost transitive tournament is a transitive one with the arc between its two extreme
  vertices reversed. It is the uniquely Hamiltonian A_n, the run-length-2 class of THM-4868."
- This is false for n >= 5.
  - The almost transitive AT_n (Achour-Boudabbous-Boussairi's definition, also Reid-Thomassen; Busch 2006 calls it
    "nearly transitive") is uniquely Hamiltonian. It has **n - 2** non-critical vertices (u_2 ... u_(n-1)) and 3 cut
    arcs.
  - THM-4868's run-length-2 class is Burzio-Demaria's bineutral A_n (Definition 2: x_i -> x_j iff j < i - 1 or
    j = i + 1; Las Vergnas). It has 2 non-critical vertices and n - 1 cut arcs.
  - The two coincide only at n = 4, where both are W4 (`aud_small_checks.out`, n = 4-10).
- Line 125 also calls the quotient "the almost transitive A_5". The almost-transitive part is right: the quotient is
  AT_5 exactly when v beats both x and y or both beat v. But the witness with x -> v, v -> y has the **bineutral**
  A_5 as its quotient (`aud_quotient5.out`), so the label A_5 is doubly misleading.
- The parallel stays ANALOGY and no theorem is affected.
- **Fix:** "an almost transitive tournament ... reversed (Busch's 'nearly transitive'). It is uniquely Hamiltonian with
  non-critical run of length n - 2 and 3 cut arcs. For n >= 5 it is not the bineutral A_n (run length 2) of THM-4868."
  At line 125, write AT_5.

**m2. MISTAKE-707 repeated: "strong" dropped (note lines 32-33).**
- The novelty sentence reads: "the bound 'minimum out-degree >= 2 and n >= 6 implies at least 3 Hamiltonian cycles'".
- This is false without "strong": a source beating RT_5 has n = 6, minimum out-degree 2 and hc = 0.
- (C3) at line 274 and the THM-4869 title do say "strong".
- **Fix:** "every strong tournament with minimum out-degree >= 2 and n >= 6 has at least 3 Hamiltonian cycles".

**m3. The converse of B3 needs "D strong".**
- Location: note lines 166-173; THM-4869 lines 82-91.
- The three listed properties (no out-degree-1 vertex in D - u; no Z with exactly one arc into D - Z; e_u(D) >= 2) do
  not imply that D is strong.
- Counterexample: D = W + u with W strong and W => u. It satisfies all three (e_u = H(W) >= 3). But in
  Q[z <- (D, u)], u has out-degree 1, a sole exit, so the converse "has no sole exit and has D as its terminal core"
  fails.
- HYP-9461's definition of a core includes "D is strong", and B4's configurations are strong (u -> w_1, or u beats
  W - A when p = 0), so no count is affected.
- **Fix:** add "D strong" to the definition in the note and in THM-4869 B3.

**m4. Same-day credit for "a forced arc needs a cut vertex" (note line 97; THM-4869 lines 70-71).** The
Fraisse-Thomassen deduction was recorded by 12:51 the same day, in S74's note
(`uniquely_hamiltonian_one_arc_cut_shelor_20261104.md` lines 131-133) and in the S74 audit
(`uh_one_arc_cut_20261104_audit_S74/AUDIT.md` lines 208-211 and 391). **Fix:** credit both. New in S75: the product
formula, hc >= 9, the forced-arc classification, and the explicit v => C3 => x => y => C3 => v description (S74 gave
only out-rows).

**m5. Lineage: "exactly half the cores of each order are attached" (note lines 344-345).** False at small orders,
in the session's own `fewcycle_core_lineage_20261104.out`: d = 6 gives 2 of 3 (k = 3), 2 of 2 (k = 4), 0 of 4 (k = 5);
d = 7, k = 5 gives 4 of 5. It holds once every core has p >= 2 (d >= 7 for k = 3, 4; d >= 8 for k = 5), because a
p = 1 core is attached iff A != W and p = 0 cores are orphans. **Fix:** state the range and the reason; adjust the
typing row at line 403.

**m6. Lineage: wrong exponent for k = 4 (note line 347).** The orphan polynomial is **2z^2 (1 + z)^(d-7)**, not
(d-6). The session's own output at d = 12 is 2, 10, 20, 20, 10, 2 = z^2 * 2(1 + z)^5; the note's formula would give
2, 12, 30, 40, 30, 12, 2. Lemma B4 gives d - 7 (two blocks of order 4, p = d - 5, p - 2 free arcs;
`aud_gf_check.out`). The k = 3 polynomial z(1 + z)^(d-7)(1 + z + z^2) is right (checked for d = 7-15).

**m7. Wording at note line 235.** "So A = {z, q}" should read "So A contains {z, q}". The bound sum >= H(Q_1) + 2 uses
only containment.

**m8. Misplaced colon in THM-4869 line 57.** "The bound is attained by v => C3 => C3 => v:" is followed by the case
analysis proving S+ >= 3, not by an attainment argument. **Fix:** "S+ >= 3, by cases: (C3 ...; W4 ...; hc(D_1) >= 2
...; uniquely Hamiltonian of order >= 5 ...). The bound is attained by v => C3 => C3 => v."

**m9. Wording at note line 341.** "no 2-strong or cut-vertex classes with hc <= 5" is wrong as written: a strong class
that is not 2-strong has a cut vertex by definition (every uniquely Hamiltonian class does), so almost all classes with
hc <= 5 have one. **Fix:** "no 2-strong classes, and no 2-arc-strong classes with a cut vertex".

**m10. Pre-registration scope (note line 330; THM-4869 line 10; HYP-9461 status).** At order 10, P3 and the P5-P7
values were klein's exhaustive numbers (HYP-9459, 11:59; P7 even quotes "klein's order-10 value is 60"), and the u_5
law was fitted on data including that value. The run outputs first appear in `60620d918c` (14:13), after
`c355721ce6` (13:48), as required; the generator changed in between only cosmetically (MAXK macro, int hc in the
hash). **Fix:** "confirmed blind at orders 11-14; order 10 reproduces klein's table" (HYP-9461 line 19 already says so).

**m11. Session31 misquote (note line 382).** "2 <-> 4 the only rational pair" is not what the cited README proves.
README lines 18 and 50 say that the rational mirrors of **integers** n >= 2 are exactly r_2 = 4 and r_4 = 2.
Non-integer rational mirror pairs exist: 9/4 <-> 27/8, and in general ((1 + 1/k)^k, (1 + 1/k)^(k+1)).
**Fix:** "the only rational mirrors of integer inputs".

**m12. Path (note line 354).** The path `gsma_norm_collapse_20261031/volume/THEORY.md` should be
`session127_gsma_norm_collapse_20261031/volume/THEORY.md`. VOL-1 there is correctly quoted, including
V_6/V_5 = 5 pi/16 < 1 at line 60.

**m13. "Not read" literature can be read (note lines 33-35 and 321-322).**
- Busch 2006 (EJC 13, #N3, four pages, open access) was read for this audit. It states:
  - Thomassen 1980: every 2-connected tournament of order n has at least alpha^(n/32 - 1) Hamiltonian cycles,
    alpha = 6^(1/4) ~ 1.565;
  - Busch's Corollary 2 improves this to beta^(n/32 - 1), beta = 5^(1/3) ~ 1.710.
- Both bounds are vacuous below n ~ 32. Kuehn-Osthus (arXiv:1006.0590, section 4.2) says the minimum for strongly
  2-connected tournaments "is still open".
- **Fix:** quote these exact statements. They support the note's novelty claim and do not contradict it.

**Optional (not errors).**
- (C3) closes the gap that `session189_response_frontier_20261104/twohc_comparison/COMPARISON.md` section 3 flagged:
  every two-cycle tournament of order >= 6 has a sole exit.
  - So klein's digit-2 ready counts (THM-4862 (g); HYP-9454 (5); A377857(n - 2)) count all two-cycle tournaments.
  - The census's P2 already gives t_2 = 1,775, 5,257, 15,360, 44,394 at orders 11-14 over all strong tournaments.
    This audit reproduced those values independently.
  - Worth a line in the note and a cross-reference in HYP-9454.
- HYP-9459 is still filed OPEN, and THM-4868 (B) still calls the two-cycle half open. After this audit, update them,
  or notify klein and kind-pasteur, citing THM-4869 (C3)-(C4).

## 3. Proof review (adversarial)

**Theorem A.**
- Only v enters D_1. 2-arc-strong then gives |N+(v) & D_1| >= 2, so |D_1| >= 3.
- Every Hamiltonian cycle runs v -> D_1 -> ... -> D_m -> v, block by block, with free junctions. This gives the
  product formula, with an empty middle product when m = 2.
- S+ >= 3 holds case by case:
  - hc(D_1) >= 2: each start >= hc;
  - C3: in-degree 1, so v beats all three;
  - W4: start_a + start_b = 1 + 2, recomputed by hand;
  - uniquely Hamiltonian of order m >= 5: the non-critical sole entry y of THM-4868 (C) has start_y >= m - 2.
- Corollary A' needs only S+, S- >= 2 together with Fraisse-Thomassen ((k + 1)-strong minus any k arcs is
  Hamiltonian, Graphs Combin. 3(3) (1987) 239-250).
- The forced-arc analysis is right: rotations of a middle cycle miss every arc, and |N+(v) & D_1| >= 2 frees the arcs
  at v.

**B1-B4.**
- B1: submodularity forces d+(D_1 & D_2) = 1. In the case D_1 u D_2 = V, the core C must send one arc into A and one
  into B, which leaves A and B unjoined. Disjoint sets satisfy |D_1||D_2| <= 2. A minimal set of order <= 2 forces a
  sink or a sole exit.
- B3, e_u >= 2: e_u = 1 forces D uniquely Hamiltonian with u its only sole exit. THM-4868 (C) then makes u
  non-critical with e_u >= |D| - 2 >= 2.
- B4: Z = W_r gives |A| >= 2. For Z meeting the chain the condition holds automatically (checked: arcs >= 2 in every
  case), so the block conditions are exactly minimality restricted to W. The counts 2^(d-w-2) come from the free arcs
  u -- w_2..w_p. Configurations are distinguished by the canonical chain.

**B5 with hc(W) = 2 (the Redei step).**
- A sole exit z exists for m >= 6 by (C3).
- If z is non-critical: e_z >= m - 2 >= 4, plus another end >= 2, gives a sum >= 6.
- If z is critical: z' lies in Q_1 (only z enters Q_1), so hc(W) = start_(z')(Q_1) * prod_(i>=2) H(Q_i) = 2.
  - Redei (H odd) forces start = 2 and singleton Q_i, so the last singleton q is a sole exit and lies in A.
  - Every Hamiltonian path of W - z ends at q, so e_z = H(Q_1). H(Q_1) is odd and >= 2, so >= 3.
  - Then sum <= 5 forces H(Q_1) = 3, so Q_1 = C3, whose starts are 1. Contradiction.
- Orders <= 5 are finite: hc = 2 occurs only at order 5 (two classes), and neither gives a block.

**Circularity.**
- The k <= 3 list and "no type-2 blocks" use no part of Theorem C. For k = 2: hc(W) >= 2 gives a sum >= 4; C3 gives 3;
  W4 gives >= 3; a uniquely Hamiltonian W of order >= 5 gives >= 4.
- P(n)/S(n) use only those facts, Theorem A, B1-B4 and the census.
- (C3) then holds for all n. Only after that do the k = 4, 5 lists use it.
- The note's "Order of the proof" paragraph (lines 244-246) is accurate.

**Theorem C.**
- S(n) from P(n): T is not 2-arc-strong (Theorem A), so it has a unique core with e_u >= 3, which forces s_l = 1 and
  Q uniquely Hamiltonian. A sole entry of Q other than z survives: for |Q| = 3 all three vertices are sole entries;
  for |Q| >= 4, z has in-degree >= 2.
- P(n) for n >= 13: a sole exit of T - v has out-degree 2 in T and beats v; dually for a sole entry. So
  V is contained in N+(O_2) u N-(I_2).
  - C(o, 2) <= 2o gives o <= 5, and o = 5 makes O_2 closed, so o <= 4.
  - |N+(O_2)| <= 3o - C(o, 2) takes the values 0, 3, 5, 6, 6 for o = 0..4, so is at most 6.
  - Hence n <= 12. Correct.
- The base, orders 6-12, is reproduced below.

## 4. Reproduction (audit code; 12 threads)

**(a) Labelled brute force with Burnside class counts.**
`aud_labelled_burnside.c`, output `aud_labelled_burnside_5_8.out`, 3 min 24 s. It covers all 2^28 labelled 8-tournaments
plus fixed points of odd-type permutations. Every Burnside remainder vanishes. Classes 12, 56, 456, 6,880 and strong
6, 35, 353, 6,008 match A000568 and A051337.

| n | t_1..t_5 | u_1..u_5 | min hc, 2-strong | min hc, 2-arc-strong with cut vertex |
|---|---|---|---|---|
| 5 | 3 2 1 0 0 | 0 1 0 0 0 | 2 | none |
| 6 | 8 5 11 6 5 | 0 0 2 3 4 | 4 (1 class; 4 classes at 5) | none |
| 7 | 21 18 38 25 38 | 0 0 5 4 4 | 8 | 9 (1 class) |
| 8 | 55 60 124 91 138 | 0 0 13 10 9 | 13 | 9 (4 classes) |

This equals klein's table (HYP-9459; sweep README section 8) entry by entry for k <= 10 at n = 5-8.

**(b) Exhaustive class enumeration (no Moon theorem).**
`aud_census.c exh 10`, output `aud_census_exh_10.out`, about 1.5 min.
- Method: every class of order m + 1 is a class of order m plus a vertex with out-set O, over all 2^m choices of O.
- Canonical form: the audit's own individualization-refinement (`aud_common.h`).
- Classes 191,536 and 9,733,056; strong 178,133 and 9,355,949 (A000568, A051337 exact).
- Order 9: t_1..t_5 = 144 191 390 315 466 and u = 0 0 34 26 23.
- Order 10: t = 377 589 1,193 1,049 1,500 and u = 0 0 89 68 60. This is klein's full k <= 10 row; for example
  k = 6-10 give u = 227 117 191 372 190.
- Minimum hc of 2-strong classes: 19 at order 9 and 30 at order 10.
- 2-arc-strong classes with a cut vertex: minimum 9 at both orders.
- Theorem A formula: 154 and 2,928 (class, cut vertex) pairs, 0 failures, min(S+, S-) = 3.
- Forced arcs: 4 classes at order 9, all matched to the four v => C3 => x => y => C3' => v patterns
  (`aud_small_checks.out`); 41 classes and 49 arcs at order 10 (`aud_forced10.out`).

**(b') Moon extension from the exhaustive order-10 list.**
`aud_census.c ext ... 5 14`, output `aud_census_ext_K5_14.out`, 2.9 min.
- Method: the child hc comes from a lowest-bit recursion over the out-set O. Every kept child was re-counted by a
  direct DP (3.70M at the last step). A sample of rejected children was re-counted too (1.05M), with 0 mismatches.

| n | t_1..t_5 | u_3, u_4, u_5 | sole-exit-free classes with hc <= 5 checked |
|---|---|---|---|
| 11 | 987 1,775 3,574 3,390 4,664 | 233 178 157 | 568 |
| 12 | 2,584 5,257 10,537 10,705 14,142 | 610 466 411 | 1,487 |
| 13 | 6,765 15,360 30,676 33,194 42,082 | 1,597 1,220 1,076 | 3,893 |
| 14 | 17,711 44,394 88,400 101,420 123,414 | 4,181 3,194 2,817 | 10,192 |

At every order 6-14:
- no class has lambda >= 2, and each has exactly one minimal one-arc-leaving set D, with |D| >= 3 and D strong;
- product-law failures 0; s_l = 1 always;
- Q = T/D is uniquely Hamiltonian, and z is a sole exit of Q;
- every class has a sole entry;
- the class count by core order equals c_k(d) R(n - d + 1) for every k and d;
- R(m) = F(2m - 5) (rooted canonical forms), m = 3-14;
- no 2-strong class and no 2-arc-strong class with a cut vertex has hc <= 5 at orders 7-14.

`ext ... 40 11` (`aud_census_ext_K40_11.out`, 13 s) gives 545,489 strong classes with hc <= 40 at order 11:
- no 2-strong class (so the minimum is > 40);
- 2-arc-strong classes with a cut vertex have minimum 9 (36 classes);
- Theorem A formula holds on all 3,195 pairs.

**(c) Direct core and block enumeration.**
Every strong D with hc <= 5 and every root u was tested against the core definition (with D strong). Results:
- c_3(d) = 2, 3, 6, 12, 24, 48, 96, 192, 384, 768 for d = 5-14;
- c_4(d) = 2, 2, 4, ..., 512;
- c_5(d) = 0, 4, 5, 10, ..., 640;
- c_2 = 0 throughout;
- every value equals the formula.

Blocks (W, A), with W strong and hc(W) <= 2, canonical up to isomorphism of (W, A):
- order 3: (C3, all);
- order 4: k = 3: 1 block, k = 4: 2 blocks, k = 5: 1 block;
- order 5: three blocks with k = 5, |A| = 2 (two of them on the same W);
- orders 6-14: none.
This is exactly Lemma B5.

**(d) Algebra (`aud_gf_check.py`).**
- R, C_3, C_4, C_5 and U_3, U_4, U_5 equal their closed forms as exact series to x^39.
- U_k equals the convolution sum_d c_k(d) R(n - d + 1).
- u_5 = 4F(2n - 12) - 3F(2n - 14) = L(2n - 12) + F(2n - 13).
- At n = 6 the convolution gives u_4 = 2, u_5 = 0; the census has 3 and 4 (plus 1 and 4 2-strong classes).

## 5. Prior art and same-day overlap

**Bang-Jensen, Gutin, Yeo, CPC 6 (1997) 255-261** (abstract, Cambridge page).
- It is a vertex-connectivity strengthening of Fraisse-Thomassen: for a partition X_1..X_l of V, T minus all arcs
  inside the parts is Hamiltonian when k >= sum_(i<l) floor(|X_i|/2) + |X_l|, and this is sharp.
- **It does not treat arc-strong connectivity or forced arcs. It does not characterise the 2-arc-strong exceptions.**

**Bang-Jensen, Havet, Yeo, arXiv:1905.11019** (introduction read). Its Theorem 5 states Fraisse-Thomassen exactly as
used, and it records Thomassen's 2-strong tournaments with an arc on no Hamiltonian cycle. It says no characterisation
of the arcs on Hamiltonian cycles is known; its arc-strong analogues concern spanning eulerian subdigraphs.
Bang-Jensen, Depres, Yeo (arXiv:1907.00853, Theorem 4) gives the same FT statement.

**Thomassen 1980, Busch 2006, Kuehn-Osthus.** These give exponential bounds that are vacuous at small n, and the
general minimum is open (m13). Nothing found gives "minimum out-degree 2, strong, n >= 6 implies hc >= 3", the
sole-exit-free counts, or the cores.

**Italian school.** Burzio-Demaria 1988 (Acta Univ. Carolinae 29, read) classifies Hamiltonian tournaments by
non-coned cycles and defines the bineutral A_n (used in m1). Demaria-Kiihl (UNICAMP reports) enumerate Douglas
(uniquely Hamiltonian) tournaments only. Nothing on k >= 2 Hamiltonian cycles was found.

**OEIS (searched 2026-11-04).** No entry for 1, 11, 38, 124, 390, 1193; 6, 25, 91, 315, 1049; 5, 38, 138, 466, 1500;
3, 4, 10, 26, 68, 178; 4, 9, 23, 60, 157; or 2, 4, 8, 13, 19, 30. The values 2, 5, 13, 34, 89 are A001519 and
5, 18, 60, 191, ... is A377857; neither entry mentions tournaments.

**Repository.**
- I read the full same-day log (`gitlog_since_20261104.txt`, 193 commits on origin/main since 00:00, plus later
  fetches through 14:44). The search covered HYP-9459, "sole-exit-free" and "minimum out-degree".
- No other lane proved HYP-9459 (1)-(3) or the 2-strong bound on 2026-11-04:
  - klein's later commits are digit spectra;
  - kind-pasteur's THM-4870 (14:04) and mac-mini's THM-4871 are unrelated;
  - Session189's two-cycle comparison (14:26) treats the gap as open.
- The only overlap is m4 (S74's Fraisse-Thomassen remark).
- The inputs THM-4868 (C), (F) and Garey/Shelor are credited.

## 6. Section 4 (inspiration) checks

- **VOL-1:** correct (m12 path). V_5 is the integer maximum.
- **Self-converse share:** correct. 1, 1/2, 2/3, 12/56, 88/456 at n = 3-7, with the unique interior local maximum at
  n = 5. Sources: S74 audit m7 (lines 270-272) and `two_selector_boolean_20261102/INSPIRATION.md` line 136
  (8 of 12).
- **Session31 mirror map:** a decreasing involution fixing e, with mu e^(-mu) = lambda e^(-lambda). Correct apart from
  the rationality misquote (m11).
- **THM-062:** palindromy a_k = a_(n-1-k), equivalently W(-r) = (-1)^(n-1) W(r); reciprocal root pairs follow. Correct.
- **THM-4139/THM-4146** (order-six lift, eight quarter-integer preperiodic points, fibre firewall) and **THM-4714**
  (-7/4 is the airplane root, angle 3/7): correct. The 3-cycle -1/4 -> -7/4 -> 5/4 was recomputed by hand. The
  reflection x -> -x - 1/2 is at `arithmetic_seams_20260921_synthesis.md` line 145; F(a) = (a^2 - 29)/4 and the mod-290
  descent are in `session146_carrier_primegame_20261103/prime290`; the Zsigmondy triad c in {0, -1, -2} plus -7/4 is
  at `collatz_mod6_20260917_zsigmondy_triad.md` line 159.
- **Lucas readings:** -29/16 = -L(7)/2^4 and -7/4 = -L(4)/2^2 are arithmetically right and typed NUMEROLOGY.
- **Achour-Boudabbous-Boussairi, arXiv:1204.2513:** the abstract and Theorem 4 were read. The owner's sentence is
  Theorem 4 verbatim, and the definitions of {-k}-self dual and almost transitive are quoted correctly. The
  identification with A_n is wrong (m1).
- **The analogies themselves** (Sharkovskii vs pancyclicity; microcosm-macrocosm via the self-reciprocal
  1 - 3x + x^2) are correctly typed. The pole-cancellation algebra is right.

## 7. Files (this folder)

- `aud_common.h`: canonical form (individualization-refinement), DPs, connectivity.
- `aud_labelled_burnside.c` -> `aud_labelled_burnside_5_8.out`.
- `aud_census.c`:
  - `exh 10` -> `aud_census_exh_10.out`; it writes `aud_strong_hcle40_10.bin`, deleted after use and regenerated by
    this command;
  - `ext aud_strong_hcle40_10.bin 5 14` -> `aud_census_ext_K5_14.out`;
  - `ext aud_strong_hcle40_10.bin 40 11` -> `aud_census_ext_K40_11.out`.
- `aud_gf_check.py`, `aud_small_checks.py`, `aud_forced10.py`, `aud_quotient5.py` and their `.out` files.
- `gitlog_since_20261104.txt`.

Build: `gcc -O3 -march=native -fopenmp -o aud_census aud_census.c`, and likewise for `aud_labelled_burnside.c`.
