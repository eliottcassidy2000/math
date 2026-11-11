# Independent audit: `character_of_two_records_20261111`

Scope: the packet at commits `1ca9c5f83b` (pre-audit checkpoint), `054cd956e2` (Run I, conjecture (bK), R12) and
`f2e03e2935` (section 5.4, literature credits, final outputs of H3 and I4); canon candidates THM-5136, THM-5137,
HYP-9667; `PREREGISTRATION.md`. Line numbers refer to the files at `f2e03e2935`. Everything below was derived by
hand or recomputed with programs written for this audit (`audit/aud_*`); no module of the author was imported.
The working tree was being edited while the audit ran (after `f2e03e2935`: README and `cycle_parities.py` modified,
HYP-9667 renamed, new untracked files `descendants23.py`, `switch_quadratic.py` and outputs); those changes are
outside the audited scope.

## Verdict

**PASS WITH CORRECTIONS for the theorems, FAIL for the conjecture as stated.**

- THM-5136 (Q), (E), (P), THM-5137 (D), (R), Proposition H, the reduction of section 5.3 and the elasticity remark:
  every proof was re-derived and every number recomputed; **no mathematical error in any statement typed PROVED**.
- **HYP-9667 (bS) and (bK) are false.** Explicit counterexamples at `N = 24` (M1). The forms restricted to doubly
  regular tournaments, (b4), (b8) and the normalised form of (bK), survive, and the audit extends them to **all 37
  doubly regular tournaments of order 23**.
- Five results or observations typed new or "not found" were on record (M2, M3, M4, M6, M7), one registration (R8)
  and half of another (R1) were decidable from the repository and are not admitted (M2), and one identification is
  false (M5).

The packet is not canonisable as submitted; it is after M1-M7 and the minor points are applied.

---

## MAJOR findings

### M1. HYP-9667 (bS) and (bK) are refuted at `N = 24`

**Where.** README summary item 4 (l. 42-54), section 5.3 (statements l. 412-418, tables l. 452-481, "The law is
global" l. 486), section 8 (row "Hadamard matrices" l. 684, remark "Small numbers again" l. 703-708: "one (the skew
Hadamard form) that has survived every test to order 28"), section 9 (R11, R12), section 10; HYP-9667 title, status
(l. 13, 15-22), statements, evidence tables, falsifiers; PREREGISTRATION "Conjecture to be recorded" of blocks 5 and
6 (l. 173, 210); commit message of `054cd956e2`.

**Counterexamples** (`aud_counterexample.py`, `aud_counterexample.out`, `aud_counterexample_single.out`; adjacency
and out-mask files `counterexample_N24*.txt`). Let `D23_0` be the first doubly regular tournament of order 23 of
`04-computation/experiments/sierpinski_ocf_20261022/drt23_codes_all.txt` (class `D23_0` of `drt23_classes.json`,
`H = 15759953146384333`), vertex 0 a new source, vertices 1..23 those of `D23_0`.

| tournament `U` | `S S^T = 23 I` | `H(U)` | `H mod 4` ((bS) predicts 1) | parities of `(c_3, ..., c_23)` ((bK) predicts all 0) |
|---|---|---|---|---|
| `(source + D23_0)` switched at the single vertex 1 | yes | `39150169277169599 (mod 2^64)` | **3** | `00000000001` (`c_23` odd) |
| `(source + D23_0)` switched at `X = 10110111011011111111001` | yes | `74048555293286923 (mod 2^64)` | **3** | `00000000001` |
| the same class, `X = 11010000100001010110010` | yes | `1 (mod 4)` | 1 | `00000000110` (`c_19`, `c_21` odd): (bK) fails although (bS) holds |

Each value was obtained three ways with separate code: subset dynamic programming modulo `2^32` (`aud_ham.c`, also on
a relabelled converse), inclusion-exclusion over vertex sets modulo `2^64` (`aud_ham_ie.c`), and cycle counts by
length (`aud_cyc.c`) with `1 + 2 c_odd = H (mod 4)`. As a cross-check only, the author's own `cycle_lengths.c` run on
`counterexample_N24*_outmasks.txt` returns `c_23 = 7` and `173 (mod 256)`: odd.

**Extent** (`aud_refute.py`, `aud_refute_*.out`, `aud_n20.py`, `aud_class.c`).

| order | switching class of `(source) +` | members examined | violating (bK) | violating (bS) |
|---|---|---|---|---|
| 4, 8, 12 | `P_3`, `P_7`, `P_11` | **all** 8, 128, 2048 | 0 | 0 |
| 16 | `CD_4` | **all** 32768 | 0 | 0 |
| 20 | `P_19`; `D_19` | 19 single-vertex switchings + 150 random, each | 0 | 0 |
| 24 | `P_23`; double of `P_11` | 6 random each (+ 1 and 23 single-vertex switchings) | 0 | 0 |
| 24 | each of the ten annealed classes (codes 0-9) | 6 random each | **45 of 60** (every class) | **30 of 60** (nine classes) |
| 24 | `D23_0` | the 23 single-vertex switchings | **18 of 23** | **12 of 23** |

In all 63 violating members found only the lengths 19, 21, 23 have the wrong parity, and `c_19`, `c_21` have equal
parity (an observation, not a claim). By the classification (one doubly regular tournament of orders 3, 7, 11; two of
orders 15 and 19) the rows for `N <= 16` cover **every** skew Hadamard tournament of these orders (the class of the
converse of `CD_4` consists of the converses), so (bS) and (bK) are true for `N <= 16`, unrefuted at `N = 20` (two
classes, all there are), false at `N = 24`: a fourth "small-number law" in the sense of the note.

**Cause.** The class-wide evidence at `N = 24` (R11, run I2) sampled only the classes of `P_23^+` and `D(P_11)^+`, the
two with large automorphism groups (253 and 55); the ten classes with trivial group were in the same data file and
were used in I1 as normalised members only. At `N = 28` only normalised members were tested.

**What survives.** For doubly regular tournaments (normalised members): (b4), (b8) and
`c_(2k+1) = [n = 3 mod 8] C((n-1)/2, k) (mod 2)` (equivalently the path form). The 24 descendants of the ten annealed
classes and of the double are doubly regular with `H = 5 (mod 8)` and all odd cycle counts even
(`aud_refute_desc*.out`). By nauty (`aud_order23_census.py`) the fourteen recorded classes, their converses and all
their descendants are **37 isomorphism classes, the full number of doubly regular tournaments of order 23** (Spence;
McKay's catalogue): 19 up to converse, only `P_23` self-converse, automorphism groups `1` (16 times), `5`, `55`,
`253`; one representative of each has `H = 5 (mod 8)` and all odd cycle counts even. So these three statements hold
for every doubly regular tournament of order at most 23; the first open order is 27 (40 of 722 classes, author's
runs).

**Correction.** Retype (bS), (bK) as REFUTED at `N = 24` (true for `N <= 16`, exhaustively). Restate the conjecture for
doubly regular tournaments only (keep (b4), (b8), the path form) with evidence "complete through order 23". Rewrite
summary item 4, the tables of 5.3, "The law is global", the merge row, the remark of section 8, the outcome text of
R11/R12 (confirmed on the registered objects; the conjecture they were meant to support is false), and the HYP file.
The lemma "no member is self-converse" and the reduction stay; the reduction now has an explicit failing instance.

### M2. `H(P_q) mod 8` is HYP-9273 (B); R8 and half of R1 were decidable from the repository

**Where.** README l. 158-161 ("Not found ...: `H(P_q) mod 4` or `mod 8`"), l. 314-318 (table: `q = 31` "not computed";
"`H(P_27) mod 256 = 95` is new"), l. 323-327 ("Observation, not proved (HYP-9667 (b8))"), section 9 (R1, R8, "Two
registrations were decidable"); PREREGISTRATION l. 32-34, 85 ("`H(P_27)` was not found in the repository"), 98;
THM-5136 status l. 22-23, 27; HYP-9667 (b8), (c).

**Prior art.** `05-knowledge/hypotheses/HYP-9273-...md` (opus-2026-10-14-S33; source
`05-knowledge/results/number_tournament_bridges_20261014.md`, sections 5, 5b, 5c):

> "(B) ... For every prime power `q ≡ 3 (mod 4)`, the Paley tournament `P_q` ... has `H(P_q) / (q(q−1)/2) ≡ 1 (mod 8)`."
> "Since `|Aff|` is odd, (B) says `H(P_q) ≡ q(q−1)/2 (mod 8)`."
> "OPEN; seven cases (`q = 3, 7, 11, 19, 23, 27, 31`). The `q = 31` case was named in advance ... `H(P_31)/465 ≡ 41 (mod 256)`."

and `number_tournament_bridges_20261014.out:55`: "`H(P_27) mod 8424 = 3159` ... under the affine group ... `H/351 = 1 mod 8`".

**Consequences.**
1. The mod-8 observation of section 4, the conjectures "to be recorded" after R1 and R8, and the Paley part of (b8)
   are HYP-9273 (B), with one more case (`q = 31`: `H = 121 (mod 256)`, so the table entry "not computed" is on record).
2. Theorem (Paley) (3)-(4) is new as a **proof of HYP-9273 (B) modulo 4** (`K_q = 1 (mod 4)`); it should be typed that
   way and HYP-9273 should receive the cross-reference. `H(P_27)` was known modulo 8424; modulo 256 it is new.
3. **R8 was decidable from the repository** (`3159 = 7 mod 8`), and so was the `P_27` half of R1 (`H = 3 mod 4` gives `F`
   odd by Session29). Neither is admitted; the sentence of the third block is wrong (HYP-9273 is a hypotheses file).
   Out-of-sample in R1 was `P_43` only.
4. HYP-9273 records "B is not a property ... of doubly regular tournaments" (control `T_4`). That control is the
   orbit form `H/|Aut| = 1 (mod 8)` (`9444525 = 5`); in the form `H = n(n-1)/2 (mod 8)` the same tournament conforms
   (`198335025 = 1 = 105`). The two hypothesis files must reference each other and say this.

### M3. Section 7 (elasticity 5/2) is on record; "because of 63" is wrong as stated

**Where.** README l. 161 ("Not found ...: factorisation lengths of path counts"), l. 169 ("Session266 ... None
treats ..."), section 7 (title l. 647, remark l. 653-663), section 8 rows l. 680 and "strong components as primes";
THM-5137 remark l. 94-98.

**Prior art.**
- `05-knowledge/results/session160_ordered_programs_20261103/axioms/NOTE.md`, section 4 (2026-11-03): "The abstract
  monoid M has elasticity exactly 5/2 ... Equality occurs at `63·63 = 49·3·3·3·3 = 3969`, with lengths 2 and 5 ...
  the actual spectrum has arithmetic elasticity at least 5/2 unconditionally; equality is conditional on the
  conjecture" (proof by the additive weight `Ω(n/7^v) + v/2`; the atom list is re-proved in the same note).
- `05-knowledge/results/session266_fiber_radix_20261111/smooth/PROOF.md`, section 6 "Exact factorization elasticity:
  the maximum length ratio is 5/2" (commit `03576e3f93`, 11:33, before the note's first content commit 12:02).
- THM-4749: `w(h)` = "the longest factorization of `h` into strong values".

New in section 7 are only the explicit length set, "exactly at the even powers of 63" and the smallest element 315
(all three confirmed by brute force below 130,000).

**Wrong causal claim.** The title "not half-factorial, because of 63" and the merge row "63 as the atom that breaks
unique factorisation length": `1225 = 35 · 35 = 5 · 5 · 49` has lengths 2 and 3 and does not involve 63; 502 of the
1525 elements below 130,000 with two lengths are not divisible by 63; `7^6 = 49^3 = 343^2`. Without 63 the monoid is
still not half-factorial (elasticity 3/2 on the elements not divisible by 63); 63 raises the elasticity from 3/2 to
5/2 and gives the smallest example (315).

**Correction.** Type section 7 KNOWN (Session160 section 4; Session266 section 6; THM-4749), keep the three additions
as a remark, and reword the title and the merge row.

### M4. Theorem D (1) is a corollary of THM-155 / HYP-668; the `F_2` facts of D (2) are on record

**Where.** README section 5.1 typing (l. 364-368: "Not found in the repository: (1) as a statement"), THM-5137 status
l. 20-22 ("NEW here: (D) 1 as a statement in this repository"), merge row "5".

**Prior art.** `01-canon/theorems/THM-155-disjoint-3cycle-identity.md` (kind-pasteur, 2026-03-12, PROVED): for every
regular tournament `c_5 − 2 ov_1 − 2 ov_2 = −3n(n²−1)(n²−9)/320`; HYP-668 (historical index l. 3936): "`c5 + 2*ov2 =
n(n^2-1)(n^2-9)/160` for ALL regular tournaments (PROVED)". Since `ov_2 = c_4`, this is the identity now cited from
the literature, and `n(n²−1)(n²−9)/160 = 12 n C((n+3)/2, 4)/5` is divisible by 4, so `c_5` is even. The Savchenko
references themselves are in `00-navigation/INVESTIGATION-BACKLOG.md` (INV-053, March 2026: "c5, c6 formulas (2016)
... c7 formula (2017) ... c8(DRT_n) is INDEPENDENT of which DRT is chosen"), not "found during the session".
For D (2): `lenses_tm_ca_fractran_tournaments_20261101.md` item 4 ("The dichotomy holds for every doubly regular
tournament (`A^(-1) = A^T` over GF(2) when `n = 3 mod 8`; `A^2 = A` when `n = 7 mod 8`)") and
`oai3_two_orbits_twos_and_threes_20261007.md` section 5 ("Every doubly regular tournament of order `n ≡ 7 (mod 8)`
gives ... a self-orthogonal `F_2`-code of dimension `(n−1)/2`").

**Correction.** D (1): KNOWN (THM-155 / HYP-668; Rowlinson, Savchenko), with THM-5132 (1) as a second proof. D (2):
credit the two notes for the linear algebra; new are the weights modulo 4 of the kernel and the value of `G_35`.

### M5. `CD_k` is not the doubling tower of THM-448 for `k >= 5`

**Where.** README section 2 (l. 182: "it is the doubling tower of THM-448 (THM-5133 (G))"), PREREGISTRATION l. 8-9,
section 5.1 proof (3) and typing ("KNOWN: THM-448, THM-5133 (G)"), R6 ("the first part was in canon (THM-448,
THM-451)"), merge rows "Cayley-Dickson".

**Fact.** The audited THM-5133 (G) says: "`CD_4` is `T_15` of THM-448 ...; `CD_5` is not `T_31` of THM-448". The audit
confirms with nauty on its own constructions (`aud_nauty.py`): `CD_4 ≅ D(CD_3)`; `CD_5` is isomorphic neither to
`D(CD_4)` nor to its converse; the same for `CD_6` and `D(CD_5)`; `CD_5`, `CD_6` are not self-converse (all groups
of order 21).

**Consequences.** Theorem D (3) stands (double regularity of `CD_k` is THM-5133 (G), the Hadamard property), but not
through THM-448. The non-existence of a converse symmetry was in canon for `CD_4` only; for `CD_5` it was a
legitimate registration (R6 is over-admitted there). "Theorems Q, E do not apply (no converse symmetry)" is
verified for `k = 4, 5, 6` and unsupported beyond.

### M6. The reach-set reading of A001187 is in the OEIS

**Where.** README section 6.3 ("the reading by reach sets was not found"), THM-5137 status l. 22-24 ("NEW here: ...
the reach-set reading `l_n = A001187(n)`").

**Prior art.** OEIS A001187, comment: "a(n) is also the number of tournaments on {1,...,n} for which 1 is reachable
from every vertex. - Don Knuth, Aug 06 2020" (the converse statement with another label). The sequence is cited by
number in the note. **Correction.** R3 (a): KNOWN (OEIS, Knuth); what remains is the refinement by max-ordered type
(R3 (c)), an instance of Proposition 6.1 of the cited paper.

### M7. A parity formula for `F` is in Session130 (the cycle-cut bridge)

**Where.** README summary item 1 (l. 25-28), l. 158 ("Not found ...: the parity of `F`"), section 3 typing (l. 268-271:
"Not found in the repository: the tournament `Z`, the parity formula, the corollaries"); THM-5136 status l. 22.

**Prior art.** `05-knowledge/results/session130_blue_black_theorems_20261031/hamiltonian/PROOFS.md`, Theorem 2
(2026-10-31; summarised in the packet README section 4 and `00-navigation/S130.md`). With `z` the fixed vertex,
`S` its out-neighbours, `B = T - z`, `W_uv` the number of ordered pairs of non-empty directed paths partitioning
`V(B)`, the first ending at `u`, the second starting at `v`, and `b_uv = (W_uv + W_vu)/2`:

> "`q = sum_(u in S,u->rho(u)) b_(u,rho(u)) (mod 2)`. (2.3) Thus the parity of *all fixed odd cycles* in the centered
> extension is detectable from the opposite-pair two-path-cover weights of the smaller tournament."

In the notation of the note `q = F`, `u = a_i^+`, `rho(u) = a_i^-`, and the sum runs over the pairs of `D_1`. So the
shape of Theorem Q (3), `F = sum over D_1 of a weight (mod 2)`, and Corollary (a) for the total `F` (empty sum, hence
`H = 1 mod 4`) were on record. The audit checked both formulas on every tournament with the fixed anti-involution
of orders 3, 5, 7 (`aud_session130.py`): no violation of either, and the weights agree **termwise**,
`b_(a_i^+, a_i^-) = #{paths of Z ending at i} (mod 2)`, for all 6210 pairs of `D_1` met.

**New in Theorem Q:** the fold `Z` (a tournament of half the size whose path counts give the weights, independent of
`D_1`), the refinement by cycle length (statement (2)), Corollaries (b), (c), and with them Theorem E for
`n = 3 (mod 8)` and the Paley theorem. **Correction.** Credit Session130 Theorem 2 in section 1 ("Converse symmetry"),
in the typing of section 3 and in THM-5136; replace "the parity of `F`: not found" by this comparison.

---

## Minor findings

1. **Trivial orders.** (bS), "all scores have one parity" and the lemma "no member is self-converse" need
   `N = 0 (mod 4)`: for `N = 2` the single arc has `S^2 = -I`, `H = 1` against `1 + N/2 = 2`, scores 1 and 0, and is
   self-converse. Add the hypothesis (README l. 412, 496; HYP statements).
2. **"Equivalently".** HYP-9667 title ("(bS) ... equivalently `H mod 4` is constant on the switching class") and
   l. 83 ("So (bS) is equivalent to ...") are inaccurate: constancy does not give the value; (bS) = constancy + (b4).
   The README sentence was corrected at `054cd956e2`, the hypothesis file was not. (Moot after M1, but the reduction
   sentence should stay in the README form.)
3. **Count.** "the fourteen doubly regular classes of order 23, twelve of them without converse symmetry" (README
   l. 478; HYP l. 21): thirteen (twelve annealed and the double of `P_11`; only `P_23` is self-converse).
4. **Evidence counts.** "86 of 86" contains 7 instances of Theorem E (`P_23`; `P_19`, `D_19`; `P_27` and three
   self-converse classes) and 2 converses of other entries (`D23_10`, `D23_11`; converse tournaments have the same
   cycle counts and the same `H`). Likewise "fourteen classes" of order 23 are twelve independent values of `H`.
5. **Section 7, proof.** "use `w <= b`: `L >= c + a - b`" does not follow from the displayed bound (that gives
   `a - b + 2c/3`); it follows from `L = c + a + v + x - w >= c + a - w`. Reword.
6. **Score-parity form.** "true for `|R| <= 5`" (l. 705, 764; HYP "true for `N <= 6`") holds for odd `|R|` only:
   `<s, g>` is odd for 2 of the 4 classes of order 4 and 32 of the 56 of order 6 (`aud_switch.out`).
7. **Orthogonal row.** "true for every tournament tried with `N <= 8`" can be stated exhaustively: a tournament has
   row `v` of `S` orthogonal to the others iff it is a switching of `(source) + R'` with `R'` regular, and all pairs
   for `N = 4, 6, 8` agree (4, 16, 192 pairs over all regular classes); at `N = 10` exactly 1152 of 3840 pairs differ
   (9 of the 15 regular classes of order 9).
8. **`det S`.** The order-4 case of `det S = 4H - 3 (mod 16)` is the recorded `|Pf| = H (mod 4)` on all 64 labelled
   4-tournaments (`collatz_procgen_20260925_incoming_synthesis.md`, F3). Credit.
9. **Registrations.** R4, R5 follow from THM-5132 (3) and the recorded `F_2` facts (M4) by the two-line argument of
   Theorem D; R7 follows from the classification (both doubly regular tournaments of order 15 have
   `H = 198335025`). They are consistency checks, not evidence; say so in the outcome table.
10. **Canonical lift.** `E[U_T^nc]` is `n!` times the canonical lift of the sprout sequence `E[U_T]/n!` (seed
    `(2+t)/(2-t)`, `b_k = 2^(1-k)` for odd `k`, `0` for even `k`); README l. 608 and THM-5137 (R) 2 drop the `n!`.
11. **THM-5136 status l. 14.** "against the known `H` for `q <= 27`": exact for `q <= 23`; `q = 27` modulo 256 only.
12. **Literature for section 5.4.** Verified: the references (Rowlinson 1986; Savchenko 2016, 2017; arXiv:2403.07629),
    the formula `c_8(DR_n) = (n+1)n(n-1)(n-3)^2(n-7)(n^2-7n+1)/2048`, Conjecture A there ("for each `m <= 9`", the
    restriction "essential"), and the identity `c_5 + 2c_4` (stated there without attribution). A paper "On the
    number of 9-cycles in a doubly-regular tournament" is announced as in preparation (arXiv:2403.07629,
    arXiv:2403.19555): cite it next to the remark that a formula for `c_9` would settle `k = 4`. Not verified: the
    text of the 2017 formula for `c_7` (the audit's derivation and exact counts confirm the note's formula) and
    that Rowlinson 1986 contains the identity.
13. **Section 5.3, short cycles.** Confirmed and strengthened: `c_3, ..., c_9` are identical (exact values) for all
    37 doubly regular tournaments of order 23 (`c_9 = 97774380`) and both of order 19 (`c_9 = 12156390`); `c_10`
    takes 11 values among the 14 recorded classes.
14. **HYP-9667 falsifiers.** "The smallest untested objects: the remaining classes of order 23" is settled by the
    census of M1; the smallest untested objects are the remaining classes of order 27.
15. **Follow-up.** For any statement on `H` as a polynomial in the switching variables (the untracked
    `switch_quadratic.py`), check THM-4524 first: "In the loop spins, `H` is a polynomial of degree at most
    `2 floor(N/3)` ... for `N <= 5`, `H` on a switching class is an Ising Hamiltonian".

---

## A. Proofs (re-derived by hand)

| item | result |
|---|---|
| Theorem Q (1)-(3) | Correct. `y_ij = y_ji` and the four transition rules follow from `psi`; `M = zI + g 1^T (mod 2)` checked in both cases; a cut term is `y(first cut) prod z`; for `r >= 3` the lists `(Q_0, Q_1, Q_sigma(2), ..., Q_sigma(r))` are distinct (disjoint non-empty segments) and have the same weight; for `r = 1, 2`, `B(u,w) = B(w,u) (mod 2)` by Redei on the two orientations (holds for `|V'| = 2` and for segments of length 1). The half-path form needs that the reflection induced on an odd cycle fixes a vertex; correct. |
| Corollaries (a)-(d), remark | Correct; in (c) the arcs inside the pairs are free under `psi`. Conjugacy of anti-involutions: correct (Session29). |
| Theorem E | Correct. The counts `t`, `t + 1` and the invariance of `E` checked; for `a` with `a -> v_0` the set contains `v_0` and the same conclusion follows. |
| Theorem (Paley) | Correct, including prime powers (`e` odd, `chi(2) = (2/p)`, `q = p (mod 8)`), the isomorphism `a -> a^2`, and `K_q = 1 (mod 4)`. |
| Theorem D | Correct given THM-5132 (1), (3): rank by Sylvester and `ker A ∩ ker A^T = 0`; `n_00 >= 1` gives "least possible". |
| R1, R2, R3, dictionary | Correct (Pieri gives each hook twice; Lemmas 3.4, 4.1, Definition 5.4, Proposition 6.1, Example 11 of the paper say what is used; `Psi_k -> l_k/(k-1)!`). |
| Section 7 | Atom list, length formula, bound, equality case correct (wording: minor 5); prior art: M3. |
| Section 5.3 | Reduction, near-Hamiltonian form (via `E_(r-1) = sum_(i even) S_i` for odd `r`), endpoint lemma, score parities (for `N = 0 mod 4`), `c_3 = (n+1)/4`, `c_5 = 0`, the double in tournament form: correct. |
| (bK) and the path form | Subset identity quoted correctly; unitriangular in the odd `j`; both formulas right (`C(2a, b) = [b even] C(a, b/2)` and `C(m,k) C(m-k, i-k) = C(m,i) C(i,k)`); (bK) implies (bS); `k = 1, 2` from THM-5132 (1) for `N = 0 (mod 4)` and scores of one parity. The statement (bK) itself is false (M1). |
| No self-converse member | Correct for `N = 0 (mod 4)`. |
| Proposition H | Correct: at most three repeated vertices, each twice, alternating; `sum t_u f_u = N_1 + 2N_2 + 3N_3`; types `(1,2,1,3)` and `(1,2,2,2)` with `2t^2(t+1)` and `t^2(t+1)` per arc; `N_3 = 3 t c_3`; `(A^3)_wu = 2t(t+1)`; `s_7 = 7q^3 - 14q^2 + 7q - 1`; both closed forms; `8 | c_7`. |

## B. Recomputation

| program (in `audit/`) | what | result | agrees |
|---|---|---|---|
| `aud_thmq.c` | every labelled tournament with the fixed anti-involution, orders 3, 5, 7, 9 (4, 64, 4096, 1048576), cycles from the definition | no violation of (2) for any `k`, of `c = F (mod 2)`, of both forms of `H mod 4`, of (a), (b), (c) | yes |
| `aud_paley.py`, `aud_ham.c` | `q = 3..31`: `D_1`, fold, `H(Z) = 1, 3, 15, 3357, 92411, 3703011, 196922055`, `F_k mod 2`, `F mod 4 = 1, 0, 3, 1, 0, 3, 0`; exact `H(P_q)`, `q <= 23`; `K_q = 1, 9, 1729, 6857869865, 62293308207033` | as in the note | yes |
| `aud_ham8.c`, `aud_paley27.py` | paths of `P_27` starting `0 -> 1` | `1 (mod 256)`, `H = 95 (mod 256)` | yes |
| `aud_drt.py` | `CD_4` (`H = 198335025`, no anti-automorphism); `D_19` (`H = 1172732961771`, one anti-involution, 9 of 9 triples cyclic); the 15 codes of order 23 and the double of `P_11`, exact `H` | all on record, all `5 (mod 8)`; modulo 16 both 5 and 13 occur | yes |
| `aud_cube.c`, `aud_theoremD.py` | `G_35` by direct count (`P_3`, `P_7`, `CD_3`, `P_11`), by walks over all extensions (`CD_4`, `P_19`, `P_23`, `P_27`, `P_31`, `CD_5`), kernel data to `P_47`, `CD_7`; `c_5` of 2524 regular tournaments | as Theorem D; `c_5` always even | yes |
| `aud_records.py`, `aud_ln.c` | all labelled tournaments, `n <= 5` (R1, R2), `n <= 6` (R3), `l_n` for `n <= 8` by complete enumeration, identity for `n <= 16`, `q_n` against hooks | `l_n = 1, 1, 4, 38, 728, 26704, 1866256, 251548592`; all identities hold | yes |
| `aud_elasticity.py` | below 130,000 | atoms as listed; ratio `5/2` only at 3969; first two-length element 315; 502 two-length elements not divisible by 63 | yes (M3) |
| `aud_switch.py` | classes at `N = 4..24`, reduction, three local forms, doubles | as the note, plus minors 6, 7; annealed classes fail (M1) | yes, except (bS) |
| `aud_cyc.c`, `aud_cycles.py` | `CD_4` cycle counts (two methods), per-vertex counts (13-cycles odd at 8 of 15 vertices), order-23 classes modulo 4, path form at orders 3-23 | as the note | yes |
| `aud_class.c` | whole classes at `N = 4, 8, 12, 16` | (bK), (bS) hold for all members; `H mod 8`: 3; 5; 7; 1 (24576) and 5 (8192) | yes |
| `aud_refute.py`, `aud_n20.py`, `aud_order23_census.py` | M1 | counterexamples; all 37 classes of order 23 conform to (b4), (b8) | **refutes (bS), (bK)** |
| `aud_heptagon.py` | Proposition H symbolically, identities by brute force, exact `c_7`, `c_8`, short cycles | all hold | yes |
| `aud_nauty.py` | `CD_k` against `D(CD_(k-1))` | M5 | **disagrees with the note** |
| `aud_session130.py` | Session130 (2.3) against Theorem Q (3), orders 3, 5, 7 | both hold; weights agree termwise | M7 |

## C. Pre-registration

- **Order of commits.** Each of the six blocks (`c98a09856a` 11:24, `953e580672` 11:38, `1432c0ac0d` 11:46,
  `6e387cf531` 12:02, `88d90c6a80` 12:25, `dcad3e3237` 12:56) precedes the first commit of the outputs of its runs
  (A-F in `6e387cf531`, G in `88d90c6a80` and `1ca9c5f83b`, H1/H2 in `1ca9c5f83b`, I in `054cd956e2`, H3 and the
  extension of I4 in `f2e03e2935`). Git shows the order of commits, not when a run was started.
- **Verdicts.** Every row of the outcome table matches its output file, including H3 (nine classes, all `7 (mod 8)`,
  folding agrees on the three self-converse ones) and I4 (10 of 10; extension 40 of 40, declared unregistered).
- **Decidable from the repository.** Admitted: R6 (correct for `CD_4` only, M5), R9/R10 at orders 19 and 23
  (confirmed). **Not admitted: R8 and the `P_27` half of R1 (M2).** Derivable: R4, R5, R7 (minor 9).
- **R11, R12.** Confirmed on the registered objects; the sampling of I2 and H1 could not detect the failure (M1).

## D. Prior art (summary)

On record and missed: HYP-9273 (B) (M2); Session160 and Session266 elasticity, THM-4749 (M3); THM-155 / HYP-668,
INV-053, the two `F_2` notes (M4); THM-5133 (G) as audited (M5); OEIS A001187 (M6); Session130 Theorem 2 (M7); the
`|Pf|` remark (minor 8).
Searched and not found: the fold tournament `Z`, the parities `F_k` by length and Corollaries (b), (c) (Session29
half-path and converse-congruences reports read in full; THM-552, THM-791, THM-4743, THM-4747, THM-4930, THM-582,
THM-4730 scanned); `H mod 4` of doubly regular tournaments or under switching; cycle parities by length; the path
form; the cube sums `G_35`; the mean Redei-Berge function as `q_n` (the scalar forms are credited correctly). Literature (generic queries): nothing
found on `H` of Paley or self-converse tournaments modulo 4; numbers of doubly regular tournaments 1, 1, 1, 2, 2,
37, 722 at orders 3 to 27 (McKay's catalogue, after Spence).

## E. Typing and wording

- Chan-Pak is stated OPEN throughout and nothing suggests progress on it.
- Merge table: the "no bridge" rows are honest; the fold ladder is typed REFUTED; rows to change: "7, 21, 63" and
  "strong components as primes" (M3), "Cayley-Dickson" (M5), "Hadamard matrices" and the remark "Small numbers
  again" (M1), "5" (M4).
- README, HYP-9667 and THM-5137 agree on "lengths 3, 5, 7 proved for doubly regular tournaments; (b4) is equivalent
  to an even number of odd cycles of length at least 9"; with M1 these refer to the doubly regular form only.

## Not checked

- The order-27 classes (H3, I4 and its extension) were not rebuilt; only `P_27` was recomputed.
- `F mod 4` at `q = 43, 47` (HYP-9667 (c)).
- The texts of Savchenko 2017 and Rowlinson 1986.
- Whether the two exact values of `H` in M1 are below `2^64` was not proved (only residues are used).

## Reproduction

    cc -O2 -o B/aud_thmq aud_thmq.c; cc -O2 -o B/aud_ham aud_ham.c; cc -O2 -o B/aud_ham_ie aud_ham_ie.c
    cc -O2 -o B/aud_ham8 aud_ham8.c; cc -O2 -o B/aud_cube aud_cube.c; cc -O2 -o B/aud_ln aud_ln.c
    cc -O2 -o B/aud_cyc aud_cyc.c; cc -O2 -o B/aud_class aud_class.c
    for m in 1 2 3 4; do B/aud_thmq $m; done
    python3 aud_paley.py B/aud_ham; python3 aud_paley27.py B/aud_ham8; python3 aud_drt.py B/aud_ham ROOT all23
    python3 aud_theoremD.py B/aud_cube big; python3 aud_records.py; for n in 1 2 3 4 5 6 7 8; do B/aud_ln $n; done
    python3 aud_elasticity.py 130000; python3 aud_switch.py B/aud_ham n24; python3 aud_cycles.py B/aud_cyc ROOT
    B/aud_class < (n and rows of P_3, P_7, P_11, CD_4)                      # CD_4: 8 minutes
    python3 aud_counterexample.py B/aud_ham B/aud_ham_ie B/aud_cyc ROOT write
    python3 aud_counterexample.py B/aud_ham B/aud_ham_ie B/aud_cyc ROOT write 10000000000000000000000 N24_single
    python3 aud_refute.py B/aud_ham B/aud_cyc ROOT sample 6; python3 aud_refute.py B/aud_ham B/aud_cyc ROOT single 0
    python3 aud_refute.py B/aud_ham B/aud_cyc ROOT desc 0 1 2 3 4 5 6 7 8 9 D
    python3 aud_n20.py B/aud_cyc ROOT; python3 aud_order23_census.py B/aud_ham B/aud_cyc ROOT
    python3 aud_heptagon.py B/aud_cyc ROOT; python3 aud_nauty.py; python3 aud_cd_tower.py; python3 aud_session130.py

`B` is a directory for binaries, `ROOT` the repository root. Largest memory use: 0.9 GB (`aud_ham` at order 24).
