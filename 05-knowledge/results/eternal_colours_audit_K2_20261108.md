# Audit K2 of THM-4930 statement 5 and the post-K1 additions (eternal colours to 1131)

2026-11-08. Independent audit of the additions made after audit K1, at commit `2b2dec3735` (branch
`mac-mini/session-20261108-firstrealizers`, `origin/main`). It also covers the uncommitted files the coordinator
announced: `sc_even12_377.out`, `sc_walk.out`, `check_walk.out` and the regenerated `sc_holes.out`.
- **Objects:**
  - THM-4930 statement 5, its proof bullets and two new remarks;
  - HYP-9520's new data;
  - note §§2, 4b, 6, 7;
  - the INDEX line and MISTAKE-738;
  - the programs `sc_odd13.c`, `sc_even_pair.c`, `sc_walk.c`, `check_witness_codes.py`, `check_witness_rev.py` and the
    rewritten `sc_holes.py`.
- **The coordinator's later text edits are audited too** (uncommitted in the worktree; `git diff` against `2b2dec3735`):
  - HYP-9520's HEURISTIC bullet "Evidence from the walks" (§1.6, correction 5);
  - THM-4930's order-12 bullet: `sc_even_pair.c` "uses no reduction: all `2^25` labelled order-10 bases ... 6,000
    order-12 values checked by direct DP ... finds exactly the audit's 66 values up to 377, and so does `sc_odd13.c`".
    This is accurate. The program enumerates all `2^(M + 2C(M,2)) = 2^25` base codes, within-pair bits included, and
    leaves the `a`-`b` arc free. Its 66 values up to 377 equal audit K1's and the auditor's order-12 data (§1.4).

**Independence.** The auditor's code is in `04-computation/experiments/eternal_colours_20261108/audit_K2/`, written from
scratch, and extends audit K1's own programs. `aud_sc13.c` uses:
- K1's labelling (`σ(v) = v XOR 1` on `{0..11}`, fixed vertex 12), not the packet's `(i, i+M)` pairing;
- K1's generic orbit construction;
- at level 3 (order 12), a full pull DP, where the packet splits at `b`;
- at level 4, the backward table of the order-12 tournament read off its forward table by anti-invariance,
  `G[S][w] = F[σS][σw]`. The packet builds separate tables;
- Landau's score criterion for strongness.

The order of extension (base 10, then the pair, then the fixed vertex) is the packet's, as the coordinator suggested.
The validation is independent: it matches the auditor's own complete censuses count for count.

## Verdict: PASS WITH CORRECTIONS. Statement 5 is confirmed independently; no mathematical error

- **Orders 12 and 13, H ≤ 1124.**
  - The auditor's enumeration gives 210 self-converse values at order 12 (41 strong) and 433 at order 13 (26 strong).
  - This is **identical to `sc_odd13_1124.out`, value for value and in every strong flag.**
  - Below 625 the order-12 data equal K1's `aud_sc12.c` in every labelled count.
- **Every number of statement 5, HYP-9520's new table and note §4b was recomputed from the auditor's data alone:**
  - the exact range `h < 1131` and the 21 black values in `[625, 1131)`, all born at order 9 (17 / 4 mod 4);
  - 83 black values in total below 1131, 66 of them `≡ 3 (mod 4)`;
  - the birth split 189 / 64 = 21 + 43, and the delays 27 / 13 / 3, with 711, 855 and 1011 at delay 4;
  - the composite-only set {121, 275, 711, 855};
  - the 26 order-13 strong values, and 41 at order 12 from 405;
  - `f_sc(6) = 17` and `f_sc(12) = 405`, with `f_sc = f` at the other orders 3-13;
  - 1125, 1127 and 1129 strong self-converse at order 9;
  - the 13 values in `[1131, 2001)` without a realisation;
  - the shares 8.30% and 2.99%.
- **The packet's validations reproduce.** `sc_odd13.c -DM=2, 3` equals the census at orders 6-9, and
  `sc_full1011_le2001.out` equals the auditor's complete orders 10 and 11 up to 2001. `sc_holes.out` is reproduced byte
  for byte by the packet's `sc_holes.py` on the regenerated census. `sc_even12_377.out` equals K1's order-12 data up to
  377 (66 values, flags identical, 377 absent).
- **Witnesses.** All 169 lines of `sc_walk.out` and all 857 lines of the corrected `sc_witness.out` re-check: H,
  mirror anti-invariance and the strong flag. The K1 defect at 1089 is gone.
- **K1 corrections.** All 19 landed.
- **The corrections below are typing, credit and cosmetic:**
  - K1's order-11 run is described as "all `2^30` labellings"; it was the `2^25` reduced codes;
  - HYP-9520's walk evidence cannot be reproduced from the committed `sc_walk.out`;
  - the f_sc remark's proof covers only `m = 3q`, and HYP-9356 should be credited;
  - the wording on 1011, and two validation claims with no committed record.

| Item | Verdict |
|---|---|
| THM-4930 statement 5 (list, exactness, spectra, birth split, delays, composite-only) | PASS (cosmetic wording on 1011) |
| THM-4930 title and status additions | PASS |
| Proof bullet "Orders 12 and 13, H < 1125" | PASS (cosmetic: two validation claims have no committed record) |
| Proof bullet "Order 11" (description of audit K1) | PASS WITH CORRECTIONS (typing: `2^30` → `2^25` reduced codes) |
| Order-12 bullet (`sc_even_pair.c`, 377) | PASS |
| Decision bullets | PASS |
| Remark: delays | PASS |
| Remark: `f_sc(m) = f(m)` iff `m ≢ 0 (mod 6)` | PASS WITH CORRECTIONS (proof completeness; credit HYP-9356) |
| HYP-9520 (table, lists, birth orders, mod 4, tests) | PASS WITH CORRECTIONS (typing: walk evidence) |
| Note §2 | PASS WITH CORRECTIONS (typing: `2^30`; cosmetic: the order-9 walk claim) |
| Note §4b, §6, §7 | PASS |
| INDEX line | PASS (update the audit status) |
| MISTAKE-738 (and 737 taken by kind-pasteur S17: confirmed) | PASS |
| `sc_odd13.c`, `sc_even_pair.c`, `sc_walk.c`, checkers, `sc_holes.py` | PASS |
| The 19 K1 corrections | all landed |

## 1. Independent recomputation

### 1.1 `aud_sc13.c` and its validation

**Method.** Every self-converse 13-tournament is isomorphic to a reduced `σ`-anti-invariant one, by K1's conjugation
argument. Write `T_12 = T − v0` and `T_10 = T_12 − {10, 11}`; both are anti-invariant. By insertion monotonicity,
`H(T_10) ≤ H(T_12 − 11) ≤ H(T_12) ≤ H(T)`. So pruning every level at `H < 1125` loses nothing. The levels are:
1. the reduced order-10 bases;
2. vertex 10: 1024 out-masks, Gray-code order, evaluated by splitting paths at the new vertex;
3. vertex 11: forced arcs and `10 → 11`, with a full DP that also gives the forward table `F` of `T_12`;
4. the fixed vertex: 64 orbit masks, Gray-code order, a split at `v0` with `G[S][w] = F[σS][σw]`.

The `CHECK` build recomputes sampled level-2 and level-4 values by direct DP. It also verifies that every order-12 and
every checked order-13 tournament is a tournament and anti-invariant.

**Validation.** The program is parametrized by the final odd order `NO`.
- With no H bound, at `NO = 5, 7, 9, 11`, the E lines (order `NO − 1`) and O lines (order `NO`) equal the auditor's
  complete reduced censuses count for count, including labelled counts and strong counts. These are K1's
  `anti_n_red.txt` and, at order 11, K1's complete `o11_full.txt`: 41,912 values.
- The `NO = 11` run made 702 level-2 and 159,040 level-4 direct-DP checks.
- With `LIM = 1125` at `NO = 11`, the output equals the censuses restricted to `H < 1125`. This tests the pruning.

### 1.2 Production: orders 12 and 13 with `H ≤ 1124`

The run used two niced workers over interleaved quarters of the base codes, about 11 min of wall time.

| quantity | value |
|---|---|
| reduced order-10 parents with `H < 1125` | 42,704 (= the count in K1's order-10 census) |
| order-11 intermediates kept | 3,157,912 |
| order-12 tournaments with `H < 1125` (reduced, labelled) | 317,728: 210 values, 41 strong |
| order-13 tournaments with `H < 1125` (reduced, labelled) | 745,552: 433 values, 26 strong |
| direct-DP checks (level 2 / level 4) | 29,016 / 96,512, no failure |

- **The packet's order 12 and 13 agree exactly.** Both value sets and strong flags of `sc_odd13_1124.out` agree.
- **The order-12 data agree with K1.** Below 625 they equal K1's `aud_sc12.c` output in every labelled count and strong
  count.
- **THM-791 is observed at order 12.** All 210 order-12 values are `≡ 1 (mod 4)`.
- **The order-13 strong values below 1125** are exactly the theorem's 26: 625, 675, 729, 775, 819, 825, 841, 867, 873,
  881, 937, 945, 947, 961, 999, 1011, 1017, 1025, 1031, 1035, 1049, 1053, 1083, 1089, 1111 and 1123.
- **`f_sc`.** The least strong self-converse value per order `m = 3..13` is 3, 5, 9, 17, 25, 45, 75, 125, 225, 405, 625.
  So `f_sc = f` except at `m = 6` and `m = 12`.

### 1.3 The lists (`aud_lists2.py`)

The script uses only the auditor's data:
- the re-verified census for orders up to 9;
- K1's complete orders 10 and 11, with every H;
- the run above for orders 12 and 13.

The decision rule is K1's, with a per-order completeness table: a centre `e` is decided once every order `m ≤ B(e)` is
complete at `e`, apart from even `m` when `e ≡ 3 (mod 4)`. The witness files are used only above the exact range, after
the line-by-line re-checks of §1.5. Results (`out/aud_lists2.out`):

- **Exact range.** Every odd `h < 1131` is decided; 1131 is the first open value. The exactness argument works as follows:
  - for `625 ≤ e < 1125`, `B(e) = 13`, and orders 12 and 13 are complete to 1124;
  - 1125, 1127 and 1129 are strong self-converse values at order 9. Each is realised by exactly one self-converse class,
    which is strong, among 241, 187 and 225 classes with that H (re-verified census);
  - 1131 is `≡ 3 (mod 4)`, so order 14 is excluded, but order 13 is complete only to 1124.
- **Black in `[625, 1131)`.** These are 671, 683, 707, 719, 723, 755, 785, 795, 815, 851, 869, 887, 907, 931, 949, 983,
  1003, 1009, 1067, 1087, 1115.
  - There are 21, all born at order 9; 17 are `≡ 3 (mod 4)` and 4 are `≡ 1 (mod 4)`.
  - Below 1131 there are 83 in all, 66 of them `≡ 3 (mod 4)`. Below 625 there are exactly K1's 62.
- **Birth against eternal colour on `[625, 1131)`.**
  - All 253 values are born by order 9: 15 at order 8 and 238 at order 9.
  - 189 are blue at birth, and 64 black = 21 eternally black + 43 later blue.
  - The delays are 1 (27), 2 (13) and 4 (3): 711, 855 and 1011, each born at 9 and first self-converse at 13.
- **Composite-only blue below 1131:** 121, 275, 711, 855.
  - `711 = 3²·79` and `855 = 3²·95` have only `d ∈ {1, 3}`.
  - With `d = 1`, `e = h` is not a strong self-converse value at any order `≤ B(h) = 13`, a decided "no".
  - With `d = 3`, the centres are 79 and 95, strong self-converse at order 7 only (one class each). The smallest
    sandwich is therefore `C₃ ⇒ Z ⇒ C₃^op` at order `3 + 7 + 3 = 13`. `aud_sandwich.py` builds both and checks H = 711
    and 855 and self-converseness.
  - 1011 = 3·337 is squarefree and first self-converse at order 13, as a strong atom. `sc_walk.out` and
    `sc_witness.out` both carry a re-checked strong order-13 witness.
- **`[1131, 2001)`.** 422 of the 435 values are realised. The unresolved 13 are 1131, 1145, 1231, 1245, 1263, 1391,
  1459, 1475, 1535, 1567, 1583, 1699, 1891.
  - 11 of them are `≡ 3 (mod 4)`.
  - All are born at order 9.
  - The tests' requirements are right: order 13 above 1124 (and order 15 for 1891) for the values `≡ 3 (mod 4)`, and
    orders 12-14 for 1145 and 1245.
  - 1475 = 5²·59 has `d = 5` with the black centre 59, so it too needs a strong atom.

### 1.4 The packet's validations, rerun

- `sc_odd13.c -DM=2` and `-DM=3` (with no bound) reproduce the census's self-converse values and strong flags at orders
  6, 7, 8 and 9.
- `sc_full1011_le2001.out` equals the auditor's complete orders 10 and 11 on `H ≤ 2001` (356 and 844 values).
- `sc_even12_377.out` equals K1's `o12_lt625.txt` up to 377: 66 values, flags identical, 377 absent.
- `sc_even_pair.c` was read. Its a/b arc rules satisfy `u → v ⟺ σ(v) → σ(u)` for `σ = (i, i+M)(a b)`. The `a`-`b` arc is
  left free (unreduced), and its tables are initialised correctly.
- `sc_holes.py` was read. Its per-order completeness table, `exact_e` and consistency block are correct. Run on the
  regenerated census and the packet's files, it reproduces the current `sc_holes.out` byte for byte.
  - One cosmetic point: `first_sc` does not use the parity exclusion for incomplete even orders. It never matters below
    1131.

### 1.5 Witnesses

- `sc_walk.out`: all 169 lines rebuilt from `TT_n` with their mirror orbits reversed, at orders 12-13. Every one is a
  tournament, mirror anti-invariant, with the stated H and strong flag (`aud_witness16.c`; `out/aud_walk_check.out`).
- `sc_witness.out`, corrected: all 857 lines pass. The 1089 line is now `rev 0-3 0-5 0-11` at order 12 and strong.
- `check_witness_rev.py` and `check_witness_codes.py` were read. They are sound independent checkers.

### 1.6 HYP-9520's walk evidence cannot be read off `sc_walk.out`

`sc_walk.out` holds one witness per value, at the least order where it was found. Against the exact inventories:
- **Order 12, `H < 1125`.** The file's order-12 lines are exactly the 41 strong values.
- **Order 13, `H < 1125`.** Its order-13 lines cover only 22 of the 26. The other four, 841, 945, 1025 and 1089, appear
  only as order-12 lines; they are strong at both orders.
- **Above 1124.** The file has 91 order-13 lines, where HYP-9520 says 95.
- **Order 14.** The file has no order-14 lines at all, where HYP-9520 says "the order-14 walk ... finds only 1125 and
  1225". Those two are recorded at orders 12 and 13.

The per-order statements may well be true of the individual walk runs, but the committed file cannot show them. On the
merged reading the evidence is still clean: every strong value of the order-12 and order-13 inventories below 1125
occurs in the file, and no line below 1125 lies outside them.

## 2. Mathematical review

- **Exactness of statement 5.** Correct.
  - For odd `h < 1125`, every decomposition `h = d²e` with `d ≥ 3` has `e < 125`, which orders ≤ 9 decide.
  - With `d = 1`, the centre has order `≤ B(h) ≤ 13`. Orders ≤ 11 are complete for every H. Orders 12 and 13 are
    complete to 1124.
  - Order 12 is excluded by THM-791 when `h ≡ 3 (mod 4)`.
  - `[1125, 1131)` contains only blue values.
- **The order-12/13 method is sound.**
  - Conjugating by `(a_p b_p)` is legitimate; it commutes with `σ`.
  - The two deletions keep anti-invariance: `v0` is fixed, and `{a, b}` is a `σ`-orbit.
  - Insertion monotonicity holds.
  - The forced arcs of `b`: `w → b ⟺ a → σ(w)`.
  - The fixed vertex's arcs: `v0 → a_p ⟺ b_p → v0`.
- **Delays and composite-only.** Correct, as in §1.3. The minimal-order arguments for 711 and 855 need the smallest
  tournament with `H = 3` (`C₃`, order 3) and the smallest strong self-converse atoms with `H = 79` and `H = 95`
  (order 7, since `max H(6) = 45`).
- **`f_sc(m) = f(m)` iff `m ≢ 0 (mod 6)`.** Correct, given Moon-Yang (THM-4745 (6): the Busch-extremal strong
  tournaments are the optimal chains, self-converse iff palindromic).
  - At `m = 3q` the optimal multiset is uniquely `{2, 3^(q−1)}`: the product `3·5^(q−1)` has a single factor 3. A
    palindrome needs the 2 in the middle, so `q` must be odd.
  - At `m = 3q+1` all parts are 3.
  - At `m = 3q+2` the arrangement `(2, 3, …, 3, 2)` is a palindrome.
  - The remark proves only the `3q` case (correction 3).
  - The data agree at `m = 3..13`, and `f(6k) = 3·5^(2k−1) ≡ 3 (mod 4)` explains the parity hole.

## 3. Priority and credit

The search covered both statements and methods: "self-converse" together with "order 12/13", "least strong
self-converse", "palindromic chain/composition", "mod 6", the new list values, and the enumeration technique.
- **Self-converse inventories of orders 12 and 13.** Nothing earlier in the repository, apart from audit K1's order 12
  below 625, which is credited.
  - THM-4744 (3) and `strong_ears_625.py` give strong-class inventories to order 13 below 625, without self-converse
    flags. They also lie below 625, not 1125.
  - THM-4743 (6) / S49 enumerated self-converse tournaments only through `n = 11`, and recorded automorphism data.
  - Session29's eight-pair probe concerns maxima of fixed-path counts.
- **Eternal colours above 625.** Nothing earlier.
  - The concurrent kind-pasteur S17 lane (`first_creator_generation_20261108/`) calls THM-4930 "a different question".
  - A grep hit in its `lineage_le10.txt` is a value ledger.
- **`f_sc` per order.** The substance is in the repository, and the packet types it correctly as KNOWN in substance:
  - THM-4744 (4)'s table has a "self-converse" row of extremal class counts, 1, 1, 2, 0, 1, 1, 1, 1, 3, 0, 1 at
    `m = 3..13`, with zeros exactly at 6 and 12;
  - HYP-9356 says "the self-converse ones are the palindromes";
  - THM-4745 (6) gives "self-converse iff palindromic".
  - The explicit `m ≢ 0 (mod 6)` criterion is not stated earlier, but follows at once. HYP-9356 should be added to the
    credit (correction 4).
- **Methods.** The pair extension with insertion monotonicity is K1's `aud_sc12.c` (credited), and the reduction is
  K1's (credited). The split-at-new-vertex evaluation is the packet's own `sc_odd11_fast.c` idea. Gray-code
  incremental updates are a standard device.

## 4. The 19 K1 corrections

All 19 were applied as intended.
- **THM-4930:** the "KNOWN, restated" and "New here" lists, the anti-involution step in statement 2 and Proof 2, the
  `-DCHECK` scope (7,469 parents), the remark's citation, the THM-4743 (1b) credit, and "forced black at birth by
  parity". The 377 adoption is consistent throughout: 60 + 114, 62 values, and 49 / 13 mod 4.
- **HYP-9520:** "conjecturally" in the title, the share wording, and the tests.
- **The note:** the §0 credits, "Decidability (KNOWN ...)", the blue share at 82%, and "forced black".
- **INDEX:** the audit status and the method credit.
- **Programs:**
  - the `sc_witness.c` flag fix (verified: 857/857 lines consistent);
  - the `sc_odd11_fast.c` header ("order 11");
  - the `LIM ≤ 65536` guards in both order-11 programs.

## 5. Corrections

Severity: **error**, **overclaim**, **typing**, **credit**, **cosmetic**. No item is an error.

1. **(typing) THM-4930, Proofs 3-5, order-11 bullet, and note §2.**
   - Old (THM): "Audit K1 enumerated all `2^30` anti-invariant labellings at order 11 with no pruning
     (`audit_K1/aud_sc.c`); its values below 625 agree."
   - New: "Audit K1 enumerated all `2^25` reduced anti-invariant codes at order 11 (within-pair arcs fixed by
     conjugation, so every class occurs), with no H bound and no pruning (`audit_K1/aud_sc.c`); its values below 625
     agree, and audit K2 reproduced its complete order 11 count for count."
   - Old (note §2): "(all `2^30` labellings, no pruning)".
   - New: "(all `2^20` and `2^25` reduced codes, every class, no pruning)".
   - The `2^30` run is THM-4743 (6)'s S49 audit, and is correctly cited as such elsewhere.
2. **(cosmetic) THM-4930 statement 5, birth bullet.**
   - Old: "1011 is a strong atom born at order 9."
   - New: "1011 (= 3·337, squarefree) is born at order 9 and is first self-converse at order 13, as a strong atom."
3. **(typing: proof completeness) THM-4930 remark on `f_sc`.**
   - Old: "at `m = 3q` the single part 2 must sit in the middle, which needs `q` odd."
   - New: "at `m = 3q` the optimal multiset is `{2, 3^(q−1)}`, and the single 2 must sit in the middle, which needs `q`
     odd; at `m = 3q+1` the composition `(3, …, 3)` and at `m = 3q+2` the composition `(2, 3, …, 3, 2)` are
     palindromes."
4. **(credit) THM-4930 remark on `f_sc`, and `depends_on`.**
   - Add HYP-9356 ("the self-converse extremal chains are the palindromes") to THM-4744 (4) / THM-4745 (6).
   - Add THM-4745 to `depends_on`.
5. **(typing) HYP-9520, "Evidence from the walks".** The per-order statements cannot be reproduced from the committed
   `sc_walk.out` (§1.6):
   - "26 values" at order 13;
   - "95 further" values;
   - "the order-14 walk ... finds only 1125 and 1225".

   Either commit the per-order walk outputs (or their summaries), or restate: "`sc_walk.out` (one witness per value, at
   the least order found) contains every strong self-converse value of the exact order-12 and order-13 inventories
   below 1125 (four of the order-13 values are recorded at order 12), and none outside them; above 1125 it has 91
   order-13 and 15 order-12 strong values, none among the 13."
6. **(cosmetic) Note §2, the walk at order 9** ("recovers all 1,173 strong self-converse values of the census"). No
   output for that run is committed. Commit it or drop the claim.
7. **(cosmetic) THM-4930, Proofs 3-5, orders 12 and 13 bullet.** Two claims have no committed record:
   - "audit K1's complete unpruned order 11 (41,912 values)" (the file keeps `H ≤ 2001` only);
   - "18,000 sampled ... values were recomputed by direct DP".

   Add "(compared in full during the session; confirmed independently by audit K2)", or commit a checksum and the
   check log. Audit K2 confirms the first claim independently (§1.1).
8. **(update, not a correction) THM-4930 status and INDEX line.**
   - "statement 5 (added after audit K1; PRE-AUDIT, audit K2 pending)" → "statement 5 (audit K2: PASS WITH
     CORRECTIONS, confirmed by an independent enumeration of orders 12 and 13 below 1125)".
   - The same in HYP-9520's status and note §4b's heading.

## 6. Genuine mathematical errors (for a MISTAKE entry)

None. Items 1 and 5 are typing: a misdescribed method, and evidence not recoverable from the committed file. If a
MISTAKE entry is wanted, the lesson is this. A HEURISTIC claim about per-run statistics needs the per-run output
committed, not a merged file that keeps only the least order per value. When describing another lane's run, copy its
own size statement (K1 §1.4: "all `2^25` reduced codes").

## 7. Not checked

- The walks themselves (`sc_walk.c`) were not rerun. Only their committed witnesses were checked; the walks are
  heuristic and not used for exactness.
- Edits made to the worktree after this report's `git diff` (about 21:40) were not seen.
- Busch's theorem and the Moon-Yang theorem are cited (THM-4744, THM-4745). They were checked only against the data:
  `f` against the census minima at orders 3-9, and `f_sc` at orders 3-13.

## 8. Files

`04-computation/experiments/eternal_colours_20261108/audit_K2/`:
- **Programs:**
  - `aud_sc13.c`: orders `NO − 1` and `NO`, bounded by `LIM`, with validation and `CHECK` modes;
  - `aud_lists2.py`: lists, comparison and the walk check;
  - `aud_witness16.c`: witness checker to order 16;
  - `aud_sandwich.py`: 711, 855;
  - `aud_run2.sh`: reproduction. It needs audit K1's work directory.
- **Outputs:** `out/`, listed in `out/README.txt`. They include:
  - `o1213_lt1125.txt`: every self-converse value of orders 12 and 13 with `H ≤ 1124`, with labelled and strong counts;
  - the validation and production logs.
- **Machine.** At most two compute processes, niced, each under 10 MB. The load average was about 7-30 during the run.
