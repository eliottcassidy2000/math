# Audit K1 of THM-4930 and HYP-9520 (eternal colours of odd Hamiltonian-path counts)

2026-11-08. Independent audit of the packet at commit `6abb620db7` (branch `mac-mini/session-20261108-firstrealizers`,
pushed to `origin/main`):
- `01-canon/theorems/THM-4930-eternal-colours-which-odd-path-counts-are-ever-self-converse.md`;
- `05-knowledge/hypotheses/HYP-9520-the-eternally-black-values.md`;
- `05-knowledge/results/eternal_colours_self_converse_values_20261108.md` (the note);
- the THM-4930 / HYP-9520 line of `05-knowledge/hypotheses/INDEX.md`;
- `04-computation/experiments/eternal_colours_20261108/` (programs and outputs).

**Independence.** All of the auditor's code (`04-computation/experiments/eternal_colours_20261108/audit_K1/`) was written
from scratch. It differs from the packet's code in four ways:
- another labelling of the anti-involution (`σ(v) = v XOR 1`, not the packet's `(i, i+5)` pairing);
- arc orbits computed generically, not by hand-coded x/y rules;
- a pull-style subset DP for H;
- Landau's score criterion for strongness (the packet uses breadth-first search).

The packet's programs were only read, and run for their own validation builds. The one shared input is THM-4922's
census generator `classgen.c`, and every field of its output was re-verified class by class (§1.1).

## Verdict: PASS WITH CORRECTIONS

- **No mathematical error was found.** Every number in THM-4930, HYP-9520, the note and the INDEX line was reproduced
  from the auditor's own data. That includes:
  - both exact lists;
  - the birth-colour comparison (65/44; 27 + 17; 12/10 and 15/7);
  - the delays (0/1/2/3/5: 135/95/5/13/1);
  - the composite-only pair {121, 275} and the order-13 realisation of 275;
  - the share figures and the 36 unresolved values;
  - the refinement table.
- **The exactness argument is sound.** THM-791(a) is proved in the repository for all even `n`: its proof section
  proves it for every tournament with a fixed-point-free involutive anti-automorphism, not only for small `n`. The parity
  steps of statement 2 therefore stand. The completeness of the order-10 and order-11 inventories was confirmed by
  independent enumerations. The order-11 enumeration uses no insertion monotonicity at all.
- **The corrections concern credit and typing.** The main credit point: the membership test with Busch decidability,
  the anti-involution enumeration of all self-converse tournaments through `n = 11`, and the order-10 self-converse
  value set are already in the repository (Session24, THM-4743 (6) / audit S49, Session29, and THM-4922's audit). They
  must not be listed as "new here". The remaining items are typing and cosmetic.
- **Supplementary results (the auditor's, beyond the packet).**
  - **377.** All self-converse 12-tournaments with `H < 625` were enumerated, and none has `H = 377`. Since 377 is
    squarefree and `B(377) = 12`, this makes **377 eternally black**, so every odd value below 625 is now decided. Below
    625 there are 62 eternally black values (7, 21 and 60 attainable); 33 of them lie in `[225, 625)`, all born at
    order 8.
  - **HYP-9520's 36.** The complete order-11 census (all H) shows that two of the 36 "unresolved" values in
    `[625, 2001)`, 1551 and 1627, are self-converse values at order 11. At most 34 remain.
  - Both are new data for the packet to adopt or not. They correct nothing: the packet marked 377 OPEN, and its 36 was
    an upper bound.

| Item | Verdict |
|---|---|
| THM-4930 statement 1 (criterion, palindrome, converse construction) | PASS WITH CORRECTIONS (credit) |
| THM-4930 statement 2 (Busch decidability, parity) | PASS WITH CORRECTIONS (credit; proof dependency) |
| THM-4930 statement 3 (below 225) | PASS (cosmetic wording) |
| THM-4930 statement 4 (225 to 625; 377 open; delays; composite-only) | PASS |
| THM-4930 proofs | PASS WITH CORRECTIONS (typing: the scope of the `-DCHECK` test; citation) |
| THM-4930 remarks | PASS WITH CORRECTIONS (citation) |
| THM-4930 status block ("new here") | PASS WITH CORRECTIONS (overclaim / credit) |
| HYP-9520 | PASS WITH CORRECTIONS (typing in the title; cosmetic) |
| Note §0-§2 (credits, definitions, data) | PASS WITH CORRECTIONS (credit) |
| Note §3-§4 (results) | PASS (cosmetic wording) |
| Note §5 (refinement table) | PASS (reproduced byte-identically, and independently) |
| INDEX line | PASS WITH CORRECTIONS (audit status, credit, 5.2%) |
| Programs and outputs | PASS (two cosmetic defects, no claim affected) |

## 1. Independent recomputation

### 1.1 The census at orders 1-9 (`aud_census.c`)

THM-4922's census (`classgen 9 cls`) was regenerated. For every class, the auditor's code recomputed five fields
independently and compared them with the census:
- H, by a pull subset DP;
- strongness, by Landau's criterion;
- c3, from the scores;
- self-converseness, by a backtracking search for an anti-automorphism;
- |Aut|, by backtracking.

| order | classes | self-converse | mismatches | sum of n!/\|Aut\| |
|---|---|---|---|---|
| 1-9 | 1, 1, 2, 4, 12, 56, 456, 6880, 191536 | 1, 1, 2, 2, 8, 12, 88, 176, 2752 | 0 | = 2^C(n,2) at every n |

The class counts are A000568 and the self-converse counts are A002785. The orbit-stabilizer sums are exact. The least
strong H per order is 3, 5, 9, 15, 25, 45, 75 at orders 3-9, which is Busch's f, and max H(9) = 3357.

### 1.2 Orders 1-9 in the anti-invariant model (`aud_sc.c`, `aud_compare_le9.py`)

All `σ`-anti-invariant labelled tournaments were enumerated with `σ(v) = v XOR 1`, plus a fixed vertex `n−1` at odd `n`.
Two modes were run:
- **full**: every arc-orbit bit;
- **reduced**: every within-pair arc oriented `2k → 2k+1`. This is legitimate because conjugating by the transposition
  `(2k 2k+1)` commutes with `σ` and reverses exactly that arc.

Every generated tournament was asserted to be a tournament and `σ`-anti-invariant. At every order 1-9 and in both modes,
the set of self-converse H values, with strong flags, is identical to the census's. At order 9 there are 1,215 values,
1,173 of them strong. This independently validates the odd-order (fixed-vertex) model that the packet's order-11 run
relies on. The packet itself validated only even orders.

### 1.3 Order 10, three ways (`aud_compare_o10.py`)

| source | self-converse values | strong |
|---|---|---|
| auditor, reduced model (`2^20` codes) | 3114 | 3039 |
| packet `sc_even10.out` (all `2^25` codes) | 3114 | 3039 |
| THM-4922 audit class census `aud_o10_byH.txt` (canonical augmentation; no strong flag) | 3114 | n/a |

The three agree exactly, including every strong flag where one is recorded. All 3,114 values are 1 (mod 4), and the
least strong value is 125 = f(10). Rebuilding the packet's own validations, `sc_even10.c -DN=6` and `-DN=8` reproduce
the census's value sets and strong flags exactly.

### 1.4 Order 11: all self-converse 11-tournaments, no pruning

`aud_sc 11 red` enumerated all `2^25` reduced anti-invariant codes of order 11, three workers on a partition of the code
range. This is every self-converse class of order 11, with no H bound, and it does not use insertion monotonicity.

- **Totals.** The labelled counts sum to exactly `2^25`. There are 41,912 distinct self-converse values at order 11,
  from 1 to 95,095 (the maximum of H at order 11), and 41,323 of them have a strong realisation. The least strong value
  is 225 = f(11).
- **Below 625.** There are 212 values, 48 of them strong. **This is identical to the packet's `sc_odd11.out`, including
  every strong flag.**
- **Checks.** Every self-converse value of orders 1, 3, 5, 7 and 9 recurs (padding). Neither 275 nor 377 occurs.

So the packet's order-11 inventory below 625 is confirmed by a method that uses none of the three things it was built
with: insertion monotonicity, the split-at-`v0` evaluation, and the four-part merge.

### 1.5 The packet's `-DCHECK` self-test, and a stronger one

- `sc_odd11_fast.c -DCHECK` passes.
- Its scope is narrower than THM-4930 says. The test applies the same `H(T') < LIM = 625` filter, so it compares only
  the extensions of the 7,469 parents with `H(T') < 625` among the codes up to `2^17`, not "every extension of the first
  `2^17` codes".
- A copy of the program, with the exit at `2^12` and `LIM = 65000`, compared all 32 extensions of each of the first
  4,097 parents of any H. All agree.
- One robustness defect, not triggered by any packet run: with `LIM > 65536` the arrays `seen11` and `strong11` (`2^16`
  entries) overflow and corrupt the DP tables. A first attempt with `LIM = 10^6` produced a spurious mismatch for this
  reason.

### 1.6 Supplementary: order 12 below 625 (`aud_sc12.c`, validated by `aud_scpair.c`)

**Method.** Every self-converse 12-tournament is isomorphic to a reduced `σ`-anti-invariant one. Write `T'' = T − {10, 11}`,
a reduced order-10 code. Then `H(T'') ≤ H(T − 11) ≤ H(T)`, by insertion monotonicity twice. The program therefore:
- takes the 27,784 reduced order-10 codes with `H < 625`;
- extends each by the `2^10` choices of vertex 10's arcs, evaluated by splitting paths at vertex 10 (Gray-code
  incremental) and keeping those with `H < 625` (1,410,552 order-11 intermediates);
- adds vertex `11 = σ(10)` with its forced arcs and computes `H(T)` by a direct DP.

The run checked the split value against direct DP 8,036 times and checked anti-invariance of every order-12 tournament.

**Validation.** The same code with the order as a parameter (`aud_scpair.c`) was run with no H bound at `n = 4, 6, 8, 10`.
It reproduces the direct reduced enumeration count for count (h, labelled count, strong count). The two enumerate the
same labelled set through an exact bijection. With `LIM = 625` at `n = 10`, it reproduces the 116 values below 625
exactly, which tests the pruning.

**Result.** There are 117 self-converse values below 625 at order 12, all ≡ 1 (mod 4) (THM-791, observed). The 12
strong ones are 405, 425, 441, 481, 513, 549, 553, 561, 565, 577, 585 and 601. Every self-converse value below 625 of
orders 2, 4, 6, 8 and 10 recurs, as padding requires. **377 does not occur**, so 377 is not a strong self-converse value
at any order: orders ≤ 11 are complete below 625 (§1.4), order 12 is now complete below 625, and `f(13) = 625`. Since
377 = 13·29 is squarefree, **377 is eternally black**. 275 does not occur either, as THM-791 predicts. The parent
session's concurrent `sc_pair12 377` run should be compared with this.

### 1.7 The lists (`aud_lists.py`, own data only)

**Inputs.** `aud_lists.py` reads only the auditor's data:
- the re-verified census values;
- the auditor's orders 10 and 11;
- optionally, order 12 below 625.

It applies the criterion, its own Busch function, and THM-791 only for orders that were not enumerated.
`out/aud_lists.out` is the run without order 12, which is the packet's setting. It reproduces every list and count of
THM-4930, HYP-9520 and the note:
- **Eternally black below 625.** Exactly the packet's 61 values. 29 are below 225 and 32 in `[225, 625)`. Mod 4 they
  split 49 / 12, and by birth order {6: 3, 7: 11, 8: 45} (plus 7 and 21). The only open value is 377. Without order-12
  data, every value other than 377 is decided, so "exact below 375, and for `h ≡ 3 (mod 4)` below 625" holds.
- **Below 225.** 109 attainable values; 65 blue at birth and 44 black. The 44 split into 27 eternally black and 17 later
  blue. Parity-forced: 22, split 12 / 10, with exactly the note's lists. Not forced: 15 / 7, also exactly the note's
  lists. The late-blue orders are 25, 27, 31 at 7; 57, 73, 89 at 8; ten values at 9; 121 at 10.
- **Below 625.** 309 attainable values; 135 blue at birth and 174 black. The 174 split into 59 eternally black, 114 later
  blue and 377. The delays are 0 (135), 1 (95), 2 (5), 3 (13), and none for 275.
- **Composite-only blue.** Exactly 121 and 275, both decided, and the claim is exact without order-12 data:
  - every other eternally blue value below 625 has an explicit strong self-converse realisation at order ≤ 11;
  - `121 ∉ Σ_sc^str` because `B(121) = 9` and the census is complete;
  - `275 ∉ Σ_sc^str` because `B(275) = 11` and orders ≤ 11 are complete below 625;
  - below 225, 121 is the only one.
- **Order 9.** 1,215 self-converse values among the 1,679 odd values up to 3357.
- **`[625, 2001)`.** The packet's data (witnesses, `d²` closure, orders ≤ 10, and order 11 below 625) leave its 36
  values. **With the auditor's complete order 11, two of the 36, 1551 and 1627, are self-converse values at order 11.**
  At most 34 remain unresolved (28 of them ≡ 3 (mod 4)), at most 34/688 = 4.9%. HYP-9520's "at most 36" stays true but
  can be tightened.
  - Explicit strong self-converse 11-tournaments, found by `aud_find11.c` and re-checked in Python (`out/aud_find11.out`,
    `out/aud_find11_check.out`), have these census-convention upper-triangle codes:
    - H = 1551: `26516080775662589`;
    - H = 1627: `26520624255470589`.

`out/aud_lists_with12.out` adds order 12 below 625. The only change is that 377 becomes eternally black. The black list
below 625 is then 62 values: 33 in `[225, 625)`, all born at order 8, with mod 4 split 49 / 13. The birth split becomes
174 = 60 eternally black + 114 later blue.

### 1.8 Witnesses (`aud_witness.c`)

Every one of the 857 lines of `sc_witness.out` was rebuilt as a tournament: start from `TT_n`, then reverse each listed
arc and its mirror. For each:
- it is a tournament;
- `μ(v) = n−1−v` is an anti-automorphism, which proves self-converseness. An independent backtracking anti-isomorphism
  search confirmed this for 694 witnesses: all those with `n ≤ 9`, and every 25th line above;
- H equals the stated value in every case.

There is one cosmetic defect. The line `W 1089 12 1 rev 0-4 1-5` is flagged strong, but the listed tournament is
`{0..5} ⇒ {6..11}` (a `33²` palindrome), which is not strong. `sc_witness.c` sets `best_s[h] = 1` when a later strong
tournament of the same order appears, but keeps the first arc list. No claim uses the witness strong flags.

### 1.9 Session24 and other inputs

- Session24's `low_H_formulas.json` has `self_converse_gf = 0` for exactly
  `{7, 19, 21, 23, 35, 39, 43, 47, 53, 55, 59, 63}` (H ≤ 73). This is the packet's list below 73. Its first-occurrence
  exponents (25, 27, 31 at `z^7`; 57, 73 at `z^8`) agree with the packet's late-blue orders.
- The birth orders and birth colours below 625 agree with THM-4922's audit fiber table (`aud_fibers_le10.txt`) on all
  1,520 values born by order 9. Only 7 and 21 are unattained below 625.
- **275.** Its only decompositions are `d = 1, e = 275` and `d = 5, e = 11`.
  - The only strong tournament with H = 5 is the order-4 class (`B(5) = 4`).
  - The only strong tournament with H = 11 is the unique order-5 class with H = 11 (`B(11) = 5`), which is
    self-converse.
  - So every self-converse realisation of 275 is a sandwich around that atom, of order at least `4 + 5 + 4 = 13`.
  - `aud_275.py` builds `B₅ ⇒ C₁₁ ⇒ B₅^op` (order 13, H = 275, self-converse) and `B₁₁ ⇒ B₁₁^op` (order 10, H = 121,
    self-converse).

### 1.10 The refinement table (note §5)

- `refine_dichotomy.py`, rerun on the regenerated census, reproduces `refine_dichotomy.out` byte for byte.
- `aud_refine.py` is an independent reimplementation. It re-checks every two-class converse pair by an isomorphism
  search rather than the census's converse-code field. It gives the same counts: 68, 208, 74, 476, 282, 475 and 616 of
  1,519.

## 2. Mathematical review

**A. The criterion (statement 1).** Correct.
- **The palindrome.** An anti-isomorphism maps strong components to strong components and reverses the condensation
  order. So `T[C_(k+1−i)] ≅ T[C_i]^op`, and the centre is strong and self-converse.
- **Multiplicativity.** H is multiplicative over `⇒`, and `H(C^op) = H(C)`.
- **The converse.** `B ⇒ Z ⇒ B^op` is self-converse for any `B` and any self-converse `Z`.
- **d ranges over Σ.** `d = ∏ H(B_i)` ranges over all of Σ, because every tournament's H is the product over its
  condensation.
- **A non-strong centre** is itself a palindrome and is absorbed.

**B. Decidability and parity (statement 2).**
- **Busch.** `f(m)` matches the census minima at orders 3-9. Busch's theorem itself is CITED (THM-4744 cites Busch,
  EJC 13 (2006) N3, Cor. 1).
- **THM-791(a) is PROVED for all even n.** Its proof section states and proves the general theorem: a tournament with a
  fixed-point-free involutive anti-automorphism `ρ` has `c_odd` even and `H ≡ 1 (mod 4)`. The argument was checked:
  - a cycle `C` fixed by `C ↦ reverse(ρC)` has a `ρ`-invariant support;
  - a `ρ`-invariant set is a union of 2-element orbits, so it has even size;
  - hence no odd cycle is fixed, the odd cycles pair up, and `c_odd` is even;
  - by OCF (THM-466, Grinberg-Stanley), `H ≡ 1 + 2c_odd ≡ 1 (mod 4)`.
  - A side remark, outside this packet: THM-791's title still calls (a) a "VERIFIED LAW (even n)". Its status line and
    its proof section say PROVED, and the proof is general.
- **The bridge from "self-converse" to "has such a ρ"** is the anti-involution lemma. `|Aut|` is odd, so if `α` is an
  anti-automorphism of order `2s`, then `s` is odd and `α^s` is an anti-involution. It fixes at most one vertex, and none
  at even order. THM-4930 states this lemma only in the proof of 3-4 and in a remark, not where statement 2 uses it
  (correction 2).
- **Independent proofs exist.** Session29's converse-congruences report (§§1, 3) proves the same parity, with a mod-8
  refinement `H ≡ 1 + 4I (mod 8)` at even order. The THM-4922 audit observed it on all 8,784 self-converse classes of
  order 10.
- **The parity of the centre** (THM-4866(E)): the centre has order `n − 2Σ|B_i| ≡ n (mod 2)`. Correct.
- **The exactness logic of statement 4** was re-derived:
  - for `h ≡ 3 (mod 4)`, every decomposition `h = d²e` has `e ≡ 3 (mod 4)`, because `d² ≡ 1 (mod 8)`;
  - a strong centre with `375 ≤ e < 625` would need order 12, which parity excludes;
  - every `e` with `d ≥ 3` is below 70.

**C. The enumeration arguments.**
- Every self-converse tournament has an anti-involution, with no fixed vertex at even order and exactly one at odd
  order. Involutions of one cycle type are conjugate in `S_n`. So a single fixed `σ` covers every class.
- The packet's x/y arc rules were checked by hand against `u → v ⟺ σ(v) → σ(u)`. They are also confirmed by the exact
  agreement of `sc_even10.out` with the auditor's generically built orbits.
- Insertion monotonicity holds: insert `v` at its first feasible position, an injection (THM-1370 (2)).
- The 32 extensions at the fixed vertex are correct: `v0 → a_i ⟺ b_i → v0`.
- The 4-way split `[k·2^23, (k+1)·2^23)` covers `[0, 2^25)` exactly. The awk merge ORs the strong flags correctly (tested).

**D. HYP-9520.**
- The conjecture is well posed. Natural density among the odd numbers is defined, and `E` is a set of odd integers.
- "(b) implies a weak form of Chan-Pak" is stated correctly. `E` contains every odd non-path-count, so a finite `E`
  means only finitely many odd non-path-counts. The repository consistently records Chan-Pak's Conjecture 6.2 (TCS 1015
  (2024), arXiv:2308.10214) as "the only odd non-path-counts are 7 and 21" (THM-115, THM-4876, THM-4887).
- The heuristic numbers are right:
  - `B(e) ≈ (3/ln 5) ln e = 1.864 ln e`;
  - `n(e) ≈ ln e / ln ln e`, from `ln(n!/2^(n−1)) ~ n ln n`;
  - at order 9, 1,215 of the 1,679 odd values up to `max H(9) = 3357` are self-converse values.
- The share figures are 27/111 = 24.3%, 32/200 = 16% and 36/688 = 5.23%.

## 3. Priority and credits

The packet credits Session24 eq. (2), Session224 §5, THM-4866(E), THM-791, THM-1370, THM-4744, THM-4875 and THM-4922.
Those credits are correct. Four earlier items are missing, and three of them make the packet's "new here" list
overclaim.

1. **Busch decidability is in Session24** (`actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md`).
   - §1 bounds every strong component by Busch: `|C| ≤ 1 + 3 log_5 q`.
   - Its "Completeness in H" paragraph declares the strong and strong-self-converse inventories globally complete for
     every odd `H ≤ 73`.
   - Its paragraph on merged multiplicities already notes that 19, 35, 39 and 63 "have no self-converse classes at any
     order".
   - So "the criterion stated as a membership test with Busch decidability" is KNOWN. What is new is the exact `B(e)` and
     the explicit lists above 73.
2. **The anti-involution enumeration of all self-converse tournaments through n = 11 is THM-4743 (6)**, from its
   independent audit S49 (`04-computation/experiments/seven_point_spheres_20261021/audit_S49/tourn_sc_census.py`).
   - It uses the same completeness argument and ran all `2^30` labellings at `n = 11`.
   - It recorded automorphism and exceptional-subset data, not H.
   - The existence and conjugacy of anti-involutions are also proved in Session29's half-path report §1 and
     converse-congruences report §1.
3. **The order-10 self-converse value set was already computed.** THM-4922's audit table
   `first_creators_genealogy_20261107/audit/aud_o10_byH.txt` gives the self-converse class count for every H at order 10.
   New in the packet are the strong flags at order 10 and the order-11 self-converse values below 625.
4. **Even-order parity.** Session29 has many sub-reports, so the remark's "Session29 section 4" is ambiguous.
   - The matching text is the converse-congruences report: §3 proves `H ≡ 1 (mod 4)` at even order, and §4 shows that
     `a_n(h)` is even for even `n` and `h ≡ 3 (mod 4)`.
   - The citation should name that report.
   - Session29's half-path report §2 has the same "no fixed odd cycle at even order" step.
   - THM-4743 (1b) also states that an anti-automorphism maps `C_i` onto `C_(m+1−i)`.

## 4. Corrections

Severity: **error**, **overclaim**, **typing**, **credit**, **cosmetic**. No item is an error.

### THM-4930

1. **(overclaim / credit) Status block, "New here".**
   - Old:
     ```
       New here:
       - the eternal colour as the stable form of the owner's two kinds;
       - the criterion stated as a membership test with Busch decidability;
       - the exact lists above 73;
       - the comparison of birth colours with eternal colours;
       - the order-10 and order-11 self-converse enumerations.
     ```
   - New:
     ```
       KNOWN, restated here:
       - the membership test and its Busch decidability (Session24 PROOF section 1: strong components have
         |C| <= 1 + 3 log_5 q; its "Completeness in H" paragraph);
       - the enumeration of all self-converse tournaments through one fixed anti-involution (THM-4743 (6), audit S49,
         all 2^30 labellings at n = 11; existence and conjugacy of anti-involutions also Session29 half-path and
         converse-congruences reports, section 1);
       - the order-10 self-converse value set (THM-4922 audit, aud_o10_byH.txt).
       New here:
       - the eternal colour as the stable form of the owner's two kinds (J_h != 0 in Session24's notation);
       - the exact lists above 73;
       - the comparison of birth colours with eternal colours;
       - the strong flags of the order-10 self-converse values and the order-11 self-converse values below 625.
     ```
2. **(typing: proof dependency) Statement 2, second bullet.**
   - Old: ``Every self-converse tournament of even order has `H ≡ 1 (mod 4)` (THM-791).``
   - New: ``Every self-converse tournament of even order has `H ≡ 1 (mod 4)`: it has an anti-involution (`|Aut|` is odd,
     so an odd power of an anti-automorphism is one), with no fixed vertex at even order, and THM-791(a) applies (also
     Session29 converse-congruences section 3).``
   - Add the same sentence as a first bullet of **Proof 2**, since the lemma is currently given only in the proof of 3-4.
3. **(typing) Proof 3-4, order-11 bullet.**
   - Old: "It was checked against direct DP on every extension of the first `2^17` codes".
   - New: "It was checked against direct DP on all 32 extensions of each of the 7,469 parents with `H(T') < 625` among
     the codes up to `2^17` (audit K1 extended the check to all extensions of the first 4,097 codes, of any H)".
4. **(credit: citation) Remark "THM-791's parity, proved".**
   - Old: "(THM-791(a); Session29 section 4)".
   - New: "(THM-791(a); Session29 converse-congruences report, sections 1 and 3:
     `05-knowledge/results/session29_framed_profiles_20261022/session29-converse-congruences/REPORT.md`)".
5. **(credit) Proof 1, first bullet.**
   - Old: "(THM-4866(E); Session24 PROOF eq. (2))".
   - New: "(THM-4866(E); Session24 PROOF eq. (2); THM-4743 (1b))".
6. **(cosmetic) Statement 3, last bullet.**
   - Old: "22 of the 44 are black at birth only because they are born at an even order with `h ≡ 3 (mod 4)` (THM-791)."
   - New: "22 of the 44 are forced black at birth by parity: they are born at an even order and `h ≡ 3 (mod 4)`
     (THM-791)."
7. **(optional update, not a correction) 377.**
   - If the packet adopts §1.6, then in the title, status, statement 4 and the INDEX line:
     - "except 377, which needs order 12" → "and 377 (decided by audit K1's order-12 enumeration below 625)";
     - "59 attainable values" → "60 attainable values";
     - "32 in [225, 625)" → "33 in [225, 625)";
     - "61 eternally black values ... 49 are ≡ 3 (mod 4)" → "62 ... 49 are ≡ 3 (mod 4)";
     - "The 174 split into 59 eternally black, 114 later blue, and 377" → "60 eternally black and 114 later blue".
   - Statement 4's list gains 377 between 359 and 415.

### HYP-9520

8. **(typing) Title.**
   - Old: "... a self-converse tournament (THM-4930) have natural density zero, and conjecturally ..."
   - New: "... a self-converse tournament (THM-4930) conjecturally have natural density zero, and (EXPLORATORY) only
     finitely many ..."
   - The title currently asserts the open conjecture (a) as a fact.
9. **(cosmetic) Title, data table and INDEX line.**
   - Old: "at most 5.2%".
   - New: "at most 5.3%" (36/688 = 5.23%).
10. **(cosmetic) Data table, row `[3, 225)`.**
    - Old: "24% of 111".
    - New: "24% of 111 (27 of the 109 attainable values, 25%)".
11. **(cosmetic) Tests.**
    - Old: "Settle the `h ≡ 1 (mod 4)` values in `[375, 625)` that remain open."
    - New: "Settle `h = 377`, the only value below 625 left open."
    - Or, if §1.6 is adopted: "Done (audit K1): 377 is eternally black; all values below 625 are decided."
    - With §1.6 adopted, the counts change: "Below 625 there are exactly 61 (with 7 and 21; 377 undecided)" → "62 (with
      7 and 21)"; the `[225, 625)` row becomes "33 | 16.5% of 200"; "49 of the 61" → "49 of the 62".
    - **(optional update)** With the auditor's complete order 11 (§1.7):
      - "at most 36 ... at most 5.2% of 688" → "at most 34 ... at most 4.9% of 688";
      - "652 of the 688" → "654 of the 688";
      - remove 1551 and 1627 from the list of values without a witness (both are self-converse values at order 11;
        `out/o11_full.txt`);
      - "(30 of them ≡ 3 mod 4)" → "(28 of them ≡ 3 mod 4)".

### The note

12. **(credit) §0, "Earlier".** Add three items:
    - "THM-4743 (6) and its audit S49: all self-converse tournaments through n = 11 via one anti-involution (method of
      §2)";
    - "Session29 half-path and converse-congruences reports: existence and conjugacy of anti-involutions; even-order
      `H ≡ 1 (mod 4)`";
    - "THM-4922 audit: the self-converse class count of every H at order 10 (`aud_o10_byH.txt`)".
13. **(credit) §1.**
    - Old: "**Decidability** (PROVED)."
    - New: "**Decidability** (KNOWN: Session24 PROOF section 1 and its 'Completeness in H' paragraph; restated with the
      exact `B(e)`)."
14. **(cosmetic) §0.**
    - Old: "the blue share near 0.83".
    - New: "the blue share outside the even-order black cell, 82% in aggregate (0.825 at order 10)".
15. **(cosmetic) §3.**
    - Old: "Half of the black-at-birth values are an artefact of THM-791's even-order parity."
    - New: "Half of the black-at-birth values are forced black by THM-791's even-order parity."
16. **(optional) §4 and §6.** The remark that order 12 "needs numba" is superseded: a pair extension of the 27,784
    order-10 parents settles it in minutes (§1.6).

### INDEX line

17. **(typing / credit)**
    - "audit PENDING" → "audit K1 PASS WITH CORRECTIONS".
    - "<= 5.2%" → "<= 5.3%".
    - After "order 10 all 2^25 anti-invariant labellings", add "(method: THM-4743 (6) / Session29)".

### Programs (packet files; cosmetic; no claim affected)

18. **`sc_witness.c` / `sc_witness.out`.** When `best_s[h]` flips to 1, update `wit[h]`, or document the flag as "some
    enumerated witness of this order is strong". The line for 1089 prints a non-strong witness under flag 1.
19. **`sc_odd11_fast.c` header.**
    - Old: "verified against direct DP at order 9 via the -DCHECK self-test".
    - New: "... at order 11 ...".
    - Also guard `LIM ≤ 65536` in `sc_odd11.c` and `sc_odd11_fast.c`, whose `seen11` and `strong11` arrays have `2^16`
      entries.

## 5. Genuine mathematical errors (for a MISTAKE entry)

None. The MISTAKE entry, if any, should record the credit overclaims as a process lesson:
- the "new here" list omitted repository priors for Busch decidability (Session24), for the anti-involution enumeration
  through `n = 11` (THM-4743 (6) / S49), and for the order-10 self-converse value set (THM-4922 audit);
- HYP-9520's title asserted an open conjecture as fact.

## 6. Not checked

- Busch's theorem itself (CITED via THM-4744). Only its values against the census minima at orders 3-9, and the least
  strong self-converse values at orders 10-12, were checked.
- THM-4875 beyond what is used. Only `d ≤ 45` matters here, and attainability of every odd `d ≤ 45` except 7 and 21 was
  read off the census.
- The text of Chan-Pak's paper (no download). The repository's consistent statement of Conjecture 6.2 was used.
- The note's statement that Codex Sessions 225-228 do not overlap.
- The packet's full order-10 and order-11 runs were not rerun with the packet's programs. They were replaced by the
  auditor's own complete enumerations (§§1.3-1.4), and the packet's validation builds were rerun.

## 7. Files

`04-computation/experiments/eternal_colours_20261108/audit_K1/`:
- **Programs:**
  - `aud_census.c`: census re-verification;
  - `aud_sc.c`: anti-invariant census, any order, full or reduced;
  - `aud_sc12.c`: order 12, `H < 625`;
  - `aud_scpair.c`: its validation copy, with the order as a parameter;
  - `aud_witness.c`;
  - `aud_find11.c`: explicit order-11 witnesses for 1551 and 1627;
  - `aud_compare_le9.py`, `aud_compare_o10.py`, `aud_lists.py`, `aud_275.py`, `aud_refine.py`, `aud_merge.py`;
  - `aud_run.sh`: reproduction.
- **Outputs:** `out/` (the small result files, listed in `out/README.txt`). Among them is `o11_full.txt`: every
  self-converse H value of order 11, with labelled counts and strong counts.
- The census files `cls_n.txt` (30 MB) are regenerated by `aud_run.sh` and are not committed.
- **Machine.** The 8-core mac-mini was shared, with load average about 98 during the run. The auditor used at most four
  processes at a time, each well under 100 MB. The order-11 census took 35 min of wall time on three workers; order 12
  took about 9 min on one.
