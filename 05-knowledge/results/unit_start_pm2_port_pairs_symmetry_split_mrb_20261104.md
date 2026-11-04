# The +-2 port pair of order 7, checked; sporadic unit-start +-2 pairs through order 14; the Z_3 / trivial split of the three-cycle Fibonacci counts; the two MRB notes; and where the rest of the prompt already stands (THM-4870)

**Session:** opus-2026-11-04-S76. Worktree `C:/tmp/wt-s76`. Canon: THM-4872. Pre-registration: HYP-9465.

## The owner's prompt

The prompt, verbatim, minus the attachments:

> merge in {6,4,8,9,12} and the 8 thomassen graphs and prior work on 105. consider that every planar triangulation
> contains linearly many pairwise non-touching reducible configurations or obstructing cycles of length at most 5.
> tournament vs. even-automorphism graph equinumerosity and attempts at a bijection between the two categories and the
> concepts we were working on in those sessions. consider how planar hypohamiltonian and hypotraceable graphs relate
> with prior work on the sum-square problem and adjacent ideas that were brough up in those sessions. consider planar
> groups and Maschke and how the 5 {Z_n,Z_n,A_4,S_4,A_5} relate to {2,5,13,34,89} and other 5's we have studied in
> this repo and how the trivial group and the cyclic group of order 2 relate to {7,21} and other pairs we have
> studied. two strong seven-vertex components have matching port response but path counts differing by two, and a
> four-arc route between them keeps the intermediate wrappers below the starting count. If the checks hold, this gives
> a reusable local ±2 adjustment across an entire clocked family, though it still does not supply the ternary parent
> needed for full spectrum coverage. two strongly connected seven-vertex components differ by two Hamiltonian paths
> while matching the selected root statistics. An independent search found a four-edge reversal route between them
> through strong tournaments; every shortest route must temporarily leave the unit root-start condition.

**Attachments** (read in full):
- M. R. Burns, "The MRB Constant via Lambert W" (Nov 2026);
- M. R. Burns, "The Rational Term 3/16 in the MRB Constant" (Nov 2026);
- C. Grogan, "Prime Pairs at the Prime Number Theorem Threshold" (20 Sep 2026). This is the paper already considered
  in S73; nothing tournament-related is used here.

**What was already done.** Earlier the same day kind-pasteur-2026-11-04-S10 answered the first half of this prompt,
nearly verbatim, in THM-4870 and its audited note `thomassen_merge_hypohamiltonian_groups_20261104.md`:
- the 8 Thomassen graphs and {6,4,8,9,12} (the distinct crossing numbers of the six non-planar ones);
- 105 (Thomassen's planar hypohamiltonian graph);
- the 2026 four-colour theorem (Inoue-Kawarabayashi-Miyashita-Mohar-Thomassen-Thorup) against r + k = n + 1;
- equinumerosity;
- hypohamiltonian graphs against the square-sum graphs;
- Maschke and the fives;
- the trivial group and C_2 against {3, 7, 21}.

kind-pasteur's HYP-9464 (confirmed at 15:34, 9aa6e0836f) continues with the {-3}-self-dual classes at order 10.

**Session190 also answered this first half** before this session's first commit (S76 audit, M1):
- CP1 (3f218bec77, 15:08): the four-colour, hypohamiltonian and Thomassen synthesis, with a boundary-state patch
  transfer;
- CP2 (a8108cfdb2, 15:22), `session190_constructive_links_20261104/`:
  - `planar/`: planar interfaces and 124 square-sum deletion cards;
  - `bijection/`: a tournament/even-graph Walsh coupling, with no natural bijection;
  - `groups/session190-groups/PROOF.md`: SO(3) against Maschke's planar groups; P_7 with F_21 and AGL(1,7); Theorem 3,
    the signed metagraph-edge stabiliser is 1 or C_2.

Session190 CP3 (15:34), after this session's checkpoint, independently found the same order-7 pair and route,
credits S76, and adds a mod-3 proof that the pair cannot give the ternary parent.

This session does not redo any of that. It adds sections 1-3 and records the cross-references in section 4.

## Status

| item | typing |
|---|---|
| 1. The pasted claims on the order-7 pair | CHECKED, all hold; the pair is the unique one at order 7 (FINITE-EXACT, exhaustive over all 2^21 labelled tournaments) |
| 1. Session189's wrapper formula H(D_m) = h + 2^m F | re-verified directly (220 cases) |
| 1. Unit-start +-2 pairs at orders 7-14: 1, 2, 5, 1, 6, 10, 4, 8 | FINITE-EXACT (orders <= 9 exhaustive; 10-14 by the proved FUNNELS F1/F2 grammar); P3 pre-registered |
| 2. F(2n-9) = F(2n-10) (Aut = Z_3) + F(2n-11) (Aut = 1); hc = 4, 5 classes rigid | PROVED (from THM-4869 and THM-4868 (E)); FINITE-EXACT at orders 6-14; P1 pre-registered at 13-14 |
| 3. The two MRB notes | CHECKED numerically (30-40 digits); readings typed ANALOGY or NUMEROLOGY |
| 4. Thomassen, four-colour theorem, equinumerosity, hypohamiltonian, Maschke | see THM-4870 (kind-pasteur, audited) |

## 1. The +-2 port pair of order 7

**Definitions** (`adaptive_metagraph_reformulation_20261103/ADAPTIVE.md` section 2). A rooted tournament (T, v) has the
profile (h, s, q, F):
- h = H(T);
- s = the number of Hamiltonian paths starting at v (the "root start"; *unit* root start means s = 1);
- q = H(T - v);
- F = the number of ordered spanning two-slot covers (P, Q) with P non-empty ending at v.

**The wrapper.** Session189's dyadic wrapper (`session189_response_frontier_20261104/returns/PROOF.md` (1)) adds a
transitive block Z of order m and a root w, with w -> Z -> U and r -> w -> U - {r}. It has H(D_m(U, r)) = h + 2^m F.
The formula holds for every rooted U and every m >= 1. If w comes first, Z follows, giving h paths. Otherwise r -> w,
and Z splits freely between the two slots, giving 2^m F paths. At m = 0 the count is h + F - s.

This is the zero-ready clock of `triangle_return_170807_20261103/ZERO_READY.md`, which Session189 inherits. This
program re-verified it directly on all 55 strong rooted classes of order 7 with s = 1, for m = 1..4 (220 cases,
0 failures). The S76 audit re-checked it on the literal wrappers for m = 1..5.

**The pair** (`plus_minus_two_port_routes_20261104.out`; exhaustive over all 2^21 labelled tournaments on {0..6} with
root 0). There are 3,040 rooted classes, of which 2,395 are strong.
- Pairs of strong rooted classes with equal F and h' = h + 2: there are 689. Exactly one has s = s' = 1:

  | | out-rows (root 0) | (h, s, q, F) |
  |---|---|---|
  | A | 8 73 3 4 15 31 61 | (29, 1, 15, 157) |
  | B | 2 84 1 7 13 31 61 | (31, 1, 9, 157) |

- Requiring also s = s' and q = q' leaves 47 pairs (91 with q = q' alone); none has s = s' = 1.
- So the "selected root statistics" that match are s = 1 and F = 157. The other statistics differ: q (15 against 9),
  the end count (10 against 9), P(2) (580 against 644) and S_w(2) (34 against 31).
- Other attachments break the +-2: a sole exit x -> 0 gives 53/53, and x beating N+(0) gives 49/73 (S76 audit).
- So the reusable +-2 is specific to the D_m wrappers. It also holds at m = 0, where h + F - s gives 185/187.

**The routes.** Breadth-first search over strong labelled tournaments, root fixed, one arc reversal per step.
- The shortest distance between the two rooted classes is **4**, in both directions. There are 38 shortest routes, the
  orderings of 5 four-arc sets.
- **None of the 38 keeps s = 1 throughout.** If s = 1 is required at every step, the distance becomes 7. So every
  shortest route must temporarily leave the unit root-start condition.
- **Exactly 2 of the 38** in each direction, two orderings of one arc set, keep every intermediate wrapper below the
  start at every clock
  (F_i <= F_0 and h_i + 2F_i < h_0 + 2F_0). For example:

      A = (29,1,15,157) -> (41,2,15,134) -> (29,1,5,89) -> (49,2,9,134) -> (31,1,9,157) = B,

  reversing the arcs {0,1}, {1,2}, {3,4}, {0,3}.

**What it gives.** For every clock m >= 1, H(D_m(B)) - H(D_m(A)) = 2. The values are 343/345, 657/659, 1285/1287, ...,
at orders m + 8.
- This is a reusable local +-2 step inside the clocked family.
- It sidesteps Session189's spacing obstruction, which concerns wrappers of a *fixed* component (no two differ by 2),
  because it switches components.

**Prior art and scope** (S76 audit, M2).
- Session176 section 1 (`session176_tournament_information_20261104/spectrum/PROOF.md`, 4fde8ff005, 2026-11-03)
  already proved a persistent single-edge +-2 gate across a dyadic clock: H = 6 * 2^r -+ 1, both strong at every r.
  Its section 2 explains why the gate does not recurse into ternary closure.
- Session185 section 5 (`session185_proof_frontiers_20261104/alternative/PROOF.md`) has a three-state version.
- Without the unit-start condition, D_m-wrapper +-2 pairs exist already at order 5: (h, s, F) = (13, 3, 15) and
  (15, 3, 15), 3 arc reversals apart, and 27 profile pairs at order 6. Unit start is exactly what makes order 7 the
  first.
- B itself is already in the repository: it is the order-7 nested core of the named H = 170807 source
  (`zero_ready_component_returns_20261103/component/COMPONENT.md` line 141).
- New here: A, the pair, the routes, the uniqueness at order 7 and the counts at orders 8-14.

**What it does not give.** The values h + 2^m F form one exponentially sparse family per pair. The +-2 step moves inside
that family. Nothing in it produces a parent near N/3, the 3H + 2 / 3H + 4 step of the C2 births (THM-4862), which
full spectrum coverage needs. So the prompt's caveat is correct.

A one-line check: 29 + 157 * 2^m is congruent to 2 + (-1)^m (mod 3), so the lower wrapper is never 2 mod 3. Session190's
spectrum packet (CP3, section 5.1) proves that the pair cannot give the ternary parent.

**Beyond order 7.**
- Unique-start rooted strong tournaments are counted by F(2n - 4) (FUNNELS F2, `polarized_path_calculus_20261103`).
  This is re-confirmed here exhaustively at orders 3-9 (1, 3, 8, 21, 55, 144, 377) and by the F1 grammar at 10-14.
- The +-2 port pairs (s = s' = 1, F = F', h' = h + 2) number:

  | order | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
  |---|---|---|---|---|---|---|---|---|
  | pairs | 1 | 2 | 5 | 1 | 6 | 10 | 4 | 8 |

- The first three looked like odd Fibonacci numbers; order 10 refutes that. The pairs are sporadic coincidences of F,
  present at every order 7-14.
- Pairs that also match q first appear at order 14 (2 pairs).
- The smallest-h pairs are (29, 31), (59, 61), (97, 99), (241, 243), (379, 381), (493, 495), (1267, 1269),
  (1959, 1961).

**NUMEROLOGY.** At order 9, the pair listed first in F order has h = 105, 107 (F = 573). 105 is the order of
Thomassen's planar hypohamiltonian graph. There is no mechanism.

## 2. The symmetry split of the three-cycle Fibonacci counts

**Theorem (PROVED, given THM-4869 and THM-4868 (E)).** For n >= 6, consider the F(2n - 9) sole-exit-free tournaments with
exactly 3 Hamiltonian cycles:
- exactly F(2n - 10) have automorphism group Z_3;
- the other F(2n - 11) are asymmetric.

For n >= 7, every sole-exit-free tournament with 4 or 5 Hamiltonian cycles is asymmetric.

*Proof.*
- **Aut(T) fixes the core.** By THM-4869 B1 the terminal core D is unique, so Aut(T) fixes D, its tail u and the head
  l.
- **It fixes the prefix pointwise.** The prefix Q = T/D is uniquely Hamiltonian with the contracted vertex z as a sole
  exit, so its automorphisms fixing z are trivial: THM-4868 (E) for |Q| >= 4, and the stabiliser of a vertex of C3
  for |Q| = 3. (The degenerate prefix is the single vertex l.)
- **It fixes the chain.** Inside D, the chain vertices are the trivial strong components of D - u, in order.
- **So Aut(T) = Aut(W, A).** Aut(T) is the group of automorphisms of the block W preserving A, since every vertex
  outside W relates to W uniformly or through A.
  - For (C3, all three vertices) this is Z_3.
  - For the W4 blocks and the three uniquely Hamiltonian 5-blocks it is trivial, because W is rigid.
- **Counting.** The C3-block cores number 2^(d-5) at each order d >= 5 (THM-4869 B4), with generating function
  x^5/(1 - 2x). Multiplying by the prefix factor x(1 - 2x)/(1 - 3x + x^2) gives x^6/(1 - 3x + x^2), i.e. F(2n - 10).
  QED.

**FINITE-EXACT** (`fewcycle_symmetry_and_unit_pairs_20261104*.out`; orders 13-14 pre-registered as HYP-9465 P1):

| n | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|
| asymmetric | 1 | 2 | 5 | 13 | 34 | 89 | 233 | 610 | 1,597 |
| Z_3 | 1 | 3 | 8 | 21 | 55 | 144 | 377 | 987 | 2,584 |

**For the prompt** (INTERPRETATION).
- The five values {2, 5, 13, 34, 89} are u_3(6..10): klein's HYP-9459 (2), proved as THM-4869, also recorded in HYP-9464
  (C) and in mac-mini's `cpq_defect_family_mirror_folds_20261104.md` section 0. By symmetry they split into consecutive
  Fibonacci pairs:
  1 + 1, 2 + 3, 5 + 8, 13 + 21, 34 + 55. For n >= 7 the trivial group takes the smaller part and Z_3 the larger (at
  n = 6 the parts are 1 + 1).
- So the Fibonacci recursion F(2n - 9) = F(2n - 10) + F(2n - 11) is this symmetry split.
- Z_3 is the rotation group of the triangle, a cyclic planar group in Maschke's list. This is consistent with
  kind-pasteur's observation that tournament symmetry is planar iff it is cyclic. C_2 cannot act on a tournament, so on
  the tournament side the pair of groups is {1, Z_3}; kind-pasteur found {1, C_2} only through the converse and the
  sign character.

## 3. The two MRB notes

**"The Rational Term 3/16".** Checked with mpmath at 30-40 digits:
- the Abel-Plana integral C_MRB = int_0^inf Im((1 + it)^(1/(1 + it))) / sinh(pi t) dt;
- M_0 = 1/4 and M_1 = 1/8;
- h(t) = t - t^3/2 - (3/4) t^5 + (9/8) t^7 + ...;
- C_MRB = 3/16 + R, with R = 0.000359642462067120248517934..., equal to the remainder integral.

Theorem 1 holds. The no-go theorem for fixed damping and the transcendence of M_k(2) are standard and were not
re-derived here.

**"The MRB Constant via Lambert W".** Checked:
- the switching points: N(a_n) = n via W_(-1) (n = 3, 8);
- the 200-term Euler-transformed interval series reproduces C_MRB to 40 digits, as in the note's table.

**A cautionary echo of the note's warning.** One mpmath extrapolation of the *paired* series returned 0.19164 here,
which is wrong. The S76 audit did not reproduce that exact value in 36 set-ups. On the paired series the alternating,
Levin and Richardson methods give 0.192-0.197; on the unpaired series Euler-Maclaurin gives 0.171 or 0.180. The correct
value comes from the alternating method on the unpaired series and from Euler-Maclaurin on the paired series.

The reproducible MRB checks are the audit's `aud_mrb.py`. Two observations for the owner:
- the attachment's table is labelled M = 51..201, but its code prints M = 50..200;
- the 200-term value is right to 67 digits; "40 digits" was a precision limit.

**The two branches are the previous prompt's mirror** (ANALOGY with an exact identity).
- On t = n^(1/n) - 1, the W_(-1) branch returns x = n (> e).
- The W_0 branch returns the mirror r_n < e with r_n^(1/r_n) = n^(1/n): r_3 = 2.478, r_4 = 2, r_5 = 1.765,
  r_8 = 1.4625. This is session31's mirror map. The same reading is in mac-mini's same-day note (section 5, 14:44)
  and in the S75 note.
- So the "same value realized twice on either side of the peak" is literally the pair of real Lambert-W branches.
- The note's parity selector lives on the far branch (switches accumulating at t = 0, x -> infinity). The mirrors
  accumulate at x -> 1.

**NUMEROLOGY.** 3/16 = 2 - 29/16. The MRB rational term equals the distance from c = -29/16 (the rational 3-cycle
-1/4 -> -7/4 -> 5/4 of x^2 + c; THM-4139) to the Mandelbrot tip c = -2.

## 4. The rest of the prompt (cross-references; nothing redone)

| prompt item | where | one line |
|---|---|---|
| {6,4,8,9,12}, the 8 Thomassen graphs, 105 | THM-4870 note, section 0 | crossing numbers of the six non-planar Thomassen graphs; 105 = Thomassen's planar hypohamiltonian graph; the LRC co-occurrence is forced |
| linearly many reducible configurations (arXiv:2603.24880) | THM-4870 (A) | the tournament analogue r + k = n + 1 (INTERPRETATION), mean fraction phi/sqrt(5) |
| equinumerosity, bijection attempts | THM-4870 section 5; THM-4720/4721, HYP-9340, THM-4531/4722/4725 | no bijection; uniquely Hamiltonian rigidity obstructs Aut-preserving maps. Section 2 adds that the sole-exit-free families are rigid or Z_3 |
| planar hypohamiltonian/hypotraceable against square-sum | THM-4870 (C) | no Q_N with N >= 3 is hypo-anything; vertex-deletion data |
| Maschke, the fives, {2,5,13,34,89} | THM-4870 section 4; section 2 here | tournament symmetry is planar iff cyclic; the five values split as {1, Z_3} Fibonacci pairs |
| trivial group, C_2, {7, 21} | THM-4870 (B) | the asymmetric spectrum misses exactly {3, 7, 21} through H <= 1125 |
| (not in this prompt) {-3}-self duality at order 10 | HYP-9464 (kind-pasteur, confirmed at 15:34, 9aa6e0836f) | TT_10 and AT_10 (H = 257 = F_3) |
| the same first half, independently | Session190 CP1/CP2 (`session190_constructive_links_20261104/`) | planar and square-sum deletion cards; Walsh coupling with no natural bijection; SO(3)/Maschke, P_7 with F_21; signed edge stabilisers 1 or C_2 |

## 5. Typing

| claim | type |
|---|---|
| order-7 pair: s = 1, F = 157, h = 29, 31; unique at order 7 | FINITE-EXACT (exhaustive) |
| route distance 4; 38 shortest routes; none keeps s = 1 (s = 1 distance 7); 2 keep wrappers below the start | FINITE-EXACT (exhaustive BFS) |
| H(D_m(B)) - H(D_m(A)) = 2 for every m | PROVED (Session189 formula (1), re-verified) |
| unit-start pairs 1, 2, 5, 1, 6, 10, 4, 8 (orders 7-14) | FINITE-EXACT |
| symmetry split F(2n-10) + F(2n-11); hc = 4, 5 rigid | PROVED; FINITE-EXACT 6-14 |
| MRB notes | numerically CHECKED |
| 105, 3/16 = 2 - 29/16 | NUMEROLOGY |
| Maschke and mirror readings | INTERPRETATION / ANALOGY |

## 6. Pre-registration record

HYP-9465 (commit 9036e4dcb9) pre-registered P1-P3 for orders 13-14 after the exploratory runs at orders <= 13 (pairs)
and <= 12 (symmetry). All hold:
- P1: 987/610 and 2,584/1,597, with hc = 4, 5 rigid;
- P2: 46,368 classes, all with s = 1;
- P3: 8 pairs at order 14.

**Typing of the predictions** (S76 audit, m10). P2 is the PROVED F2 count, and P1 follows from the proof sketch. Both
are consistency checks. P3 is a low-risk existence prediction. The order 10-13 counts quoted in HYP-9465 had no
committed output at 9036e4dcb9.

## 7. Audit

Independent adversarial audit:
[`AUDIT.md`](../../04-computation/experiments/unit_start_pm2_pairs_20261104_audit_S76/AUDIT.md). Verdict:
**PASS WITH CORRECTIONS.**

**Reproduced with different methods.**
- Burnside counting at order 7;
- a different unit-start generator, checked by brute force through order 8;
- its own census, canonical form and automorphism counter for the symmetry split (orders 6-14);
- the literal wrappers, m = 1..5;
- the MRB checks to 40-60 digits.

The Theorem B proof is correct.

**Majors (provenance only), applied.**
- **M1.** Session190's CP1 and CP2 had answered the same first half before this session's first commit.
- **M2.** The +-2 clocked adjustment has prior art: Session176's gate H = 6 * 2^r -+ 1, Session185, non-unit pairs at
  order 5, and B as the core of the 170807 source.

**Minors, all applied.**
- The 47/91 counts.
- Order ranges in the title.
- Proof dependencies.
- The general wrapper proof, credited to the zero-ready clock.
- Arc sets of the routes.
- The scope of "port response".
- The source of {2, 5, 13, 34, 89}.
- MRB wording and credits.
- The section 4 rows.
- Typing of the predictions.
- Stale program comments.
- The n = 6 wording and the mod-3 check.
- The INDEX lines.

**Logged as MISTAKE-711.**

## 8. Reproduction

    gcc -O3 -march=native -fopenmp -o pm2 plus_minus_two_port_routes_20261104.c && ./pm2            # order 7, exhaustive
    gcc -O3 -march=native -o unit unit_start_pm2_pairs_20261104.c && ./unit 9                       # all classes to 9
    gcc -O3 -march=native -fopenmp -o ugram unit_start_pm2_pairs_grammar_20261104.c && ./ugram 14   # F1/F2 grammar
    gcc -O3 -march=native -fopenmp -o sym fewcycle_symmetry_and_unit_pairs_20261104.c && ./sym 14   # symmetry split
