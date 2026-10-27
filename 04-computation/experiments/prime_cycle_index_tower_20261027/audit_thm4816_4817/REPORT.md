# Independent audit of THM-4816 and THM-4817 (prime cycle-index tower; order-12, 13, 14 ledgers) and MISTAKE-683

**Verdict: PASS WITH CORRECTIONS** (2026-10-27). No major finding. The tower principle (THM-4816 A) is correctly
proved, and every number in both theorems was reproduced exactly with independent code: a non-rigid census through order
13, an exact-rational ledger and an independent Burnside count. MISTAKE-683's numbers are confirmed. There are eleven
minor findings: wording, overstated evidence, incomplete propagation of MISTAKE-683, a conditional that does not apply, a
latent array bound in the producer's code, and a credit nit. The auditor wrote 18 exact edits (`required_edits.json`;
each old string occurs once, the edits apply in sequence, and the edited front matter parses), and the producer applied
all 18. The harness blocked the auditor's own report file, so this is its report as delivered, saved by the producer.
No repository file was edited and nothing was committed by the auditor.

## 1. Proofs (THM-4816 A): correct

- **A1 (cycle lemma).** Correct. Tournament automorphisms have odd order, and a lone p-cycle with at least one fixed
  point is a module, since each fixed vertex relates uniformly to every cycle.
  - Completeness needs only prime-order elements. Every non-rigid group contains one (Cauchy), and nauty recomputes the
    full group of each class found. Elements of order 9, 15, 21, … are therefore never needed as generators: the
    order-9 elements at orders 12 and 13 are found through their cubes, of types `1^3 3^3` and `1^4 3^3`.
  - The allowed type lists are complete: order 11 `1^5 3^2, 1^2 3^3, 1^1 5^2, 11`; order 12 `1^6 3^2, 1^3 3^3, 3^4,
    1^2 5^2`; order 13 `1^7 3^2, 1^4 3^3, 1^1 3^4, 1^3 5^2, 13`.
  - The side claims hold: no element of order 7 at orders 8-13, none of order 9 at order 11, and none of order 15 at
    orders ≤ 15.
- **A2 (invariant parametrisation).** Correct, including the bit count and the totals 657,440 / 15,736,832 / 373,358,656.
  Taking `T[F]` up to isomorphism is legitimate: the permutations of `F` that fix every cycle vertex commute with the
  representative `g0`.
- **A3 (rigid self-converse anti-automorphism type), A4 (ledger closure).** Correct; A4 re-derived.
- **A5 (labelled identity).** Correct, with two remarks (findings 4 and 5).

## 2. Independent census (own C and nauty code, written before reading `nonrigid_primes.c`)

| order | prime invariant tournaments | non-rigid classes | agreement with the producer |
|---|---:|---:|---|
| 11 | 331,808 | 12,155 | identical |
| 12 | 10,394,448 | 285,705 | identical (33 s) |
| 13 | 240,616,480 | 11,256,145 | identical (16.4 min CPU, 269 MB) |

- The per-case counts (enumerated, prime, new) are identical to the producer's logs at orders 12 and 13.
- Class by class at orders 3-12: every line of the producer's `nrclasses{n}.txt` was rebuilt and re-canonicalised. The
  class sets are identical, as are all stated `|Aut|` and self-converse flags. At order 12 the SHA-256 of the sorted keys
  is `c5a8bb9d…09ac135` for both.
- The full type tables (this census plus independent rigid rows) reproduce klein's tables for orders 3-10, the
  producer's `t11`, `t12`, `t13`, and klein's order-11 table in THM-4818.
- The primality test was cross-checked against a brute force over all subsets on 6,000 random order-9 tournaments: no
  mismatch.

## 3. Independent ledger (exact fractions)

- `A(n)` and `SC(n)` by Burnside (constraint propagation on each cycle type) equal OEIS A000568 and A002785 through
  `n = 16`.
- All 56 entries of THM-4817's table are reproduced exactly (`P(12)` to `P(14)`, `SCP`, `c_s`, `c`, `s_lab`, `P_lab`,
  `l`), as are the ratios and the parity statements. The THM-4803 values of `c(n)` and `l(n)` for `n ≤ 10` come out too.
- The labelled identity holds at every order 3-13, and the Burnside defects 1337728/165, 8571902/45, 4389979016/585 and
  3574801896/7 are reproduced.
- OEIS: none of the new sequences (self-converse primes, labelled primes, `c`, `c_s`, `l`) is found; the control search
  finds A000568. Aside: A003505 has the same name as A280772 but gives 203 and 3,785 at orders 7 and 8; the census
  supports A280772.

## 4. Circulants, histograms, narrative

- Circulants by direct canonical forms: 4 classes at order 11 (`|Aut|` 11 ×3, 55 ×1) and 6 at order 13 (13 ×5, 39 ×1);
  Burnside agrees.
- The `|Aut|` histograms at orders 11-13 are as stated, including `Z_9` = 8 and `Z_3^2` = 24 at order 12 and `Z_9` = 16
  at order 13.
- Every row obeys the cycle lemma and the 2 mod 4 law for anti-automorphisms (THM-4818 (A)).

## 5. MISTAKE-683: confirmed

- The printed `P_lab(12)` is exactly the nearest double to the true value: 1,024 above it, a multiple of 8,192. `l(12)`
  is off by +512.
- Mechanism, in klein's `species12b.py`: `Plab[4]` is the integer 0 of an empty sum, and `Plab[4] / factorial(4)` is the
  float 0.0. That makes the whole non-prime sum a float from `n = 5` on. The non-prime count itself happens to be exactly
  representable; the subtraction from `s_lab(12)` is not.
- `S_4 ⇒ S_4`: brute force over `8!` permutations gives `|Aut| = 1` and exactly one anti-automorphism, so it is rigid and
  self-converse; klein's `an8.txt` lists it with `sc=1`.
- The replacement examples hold: the order-6 class is rigid, not self-converse, with `L = 1`; both `S_4 ⇒ S_5'` are
  rigid, not self-converse, with `delta0 = 1` and no eligible vertex.

## 6. Order-14 partial census (new; cheap)

- Type `7^2`: the 8,192 invariant tournaments give 98 prime classes, 93 with `|Aut| = 7` and 5 with `|Aut| = 21`. So the
  conditional "≡ 6 (mod 7)" of THM-4817 and HYP-9423 (6) does not apply. The correct 7-adic test,
  `93 · 6/7 + 5 · 20/21 ≡ 1/7 ≡ 3574801896/7 (mod Z_(7))`, passes.
- Type `1^4 5^2`: 4,200 classes, all with `|Aut| = 5`; `4,200 ≡ 0 (mod 5)`, so the 5-part is integral, as the
  denominator 7 requires.
- The order-14 non-rigid count is at most 766,028,109.
- The producer's own `out/nr14_p5.txt` and `out/nr14_p7.txt`, committed at 10:32 during this audit, are identical to the
  audit's.

## 7. Credits

Fair overall. The producer's reservation (09:14) predates klein's commit (10:06:44), klein's packet acknowledges this and
takes no credit for the ledger, and THM-4816 / THM-4817 credit klein's THM-4814 lemmas and the concurrent THM-4818 (B)
accurately (one nit: finding 11).

## Minor findings (all corrected by the 18 edits)

1. THM-4816 E: "cosets of the cubes `{1, 3, 9}`" is wrong. `{1, 3, 9}` is the subgroup of order 3 of `Z_13^*` (the
   fourth powers); the cubes are `{1, 5, 8, 12}`.
2. Note section 6: "the first `Z_9` … groups" at order 12 is false, since `Aut ≅ Z_9` already occurs at order 9 (two
   circulants). Only `Z_3^2` is new at 12.
3. THM-4817 status and title: the agreements `S(n) = A051337(n)` and `s_lab(n) = A054946(n)` follow from A000568 and
   `2^C(n,2)` alone. They check the inputs, not the tower; the same holds for klein's "independent check"
   A000568(12) − A051337(12).
4. THM-4816 A5 scope: the labelled identity tests only `P(n)` and the single sum `Σ (1 − 1/|Aut|)` over non-rigid classes.
   It is blind to cycle types, self-converse flags and anti-automorphism types, which this audit now confirms through
   order 13.
5. Notation (THM-4816 A5, THM-4817): "1 + T = 1/(1 − S)" is inconsistent with "(T − 1)^k". It should be `T = 1/(1 − S)`,
   with `T` including the empty tournament.
6. THM-4817 and HYP-9423 (6): the "≡ 6 (mod 7)" conditional does not apply (section 6). The verified 7-adic and 5-adic
   parts should be recorded.
7. MISTAKE-683, THM-4817, THM-4818 wording: only `P_lab(12)` is a double rounding; `l(12)` was computed exactly from it
   and inherits the error. The mechanism of section 5 should be stated.
8. MISTAKE-683 propagation: the rounded values were still uncorrected in `05-knowledge/results/INDEX.md` line 3,
   HYP-9422 line 87 (klein's update paragraph), and klein's packet README lines 92 and 95. The Repair line should say
   whether the packet is left as filed.
9. THM-4816 F / HYP-9423 (5): the expression called an "invariant count" is a class estimate. With the labelled prime
   fraction (0.3845, 0.4638, 0.5393, 0.6029 at orders 10-13) it reproduces the counts 700 / 10,808 / 268,116 / 10,960,216
   to within 0.003%; with the unweighted fraction, only to within 1.4-3%.
10. `nonrigid_primes.c`: `build_frep` holds at most 2,048 representatives, with no guard, while `N ≤ 16` is accepted.
    Since `t(8) = 6,880`, an order-14 run of type `1^8 3^2` would overflow. No filed run is affected.
11. HYP-9422 status line: it should also credit klein's concurrent THM-4818 (B) resolution.

## Files in this directory

- `required_edits.json` (18 edits), `reproduce.sh`;
- census: `audit_census.c`, `audit_census_p.c` (single-prime variant, used for order 14), `my_nr{3..13}.txt`,
  `my_nr14_p{5,7}.txt`, `my_stats*.txt`;
- class comparison: `audit_reclass.c`, `compare_classes.out`, `class_digests_n12.txt`;
- ledger: `audit_burnside.py`, `audit_ledger.py`, `audit_main.py`, `audit_main.out`, `audit_main_prod13.out`,
  `audit_ledger_values.txt`;
- checks: `narrative_checks.py`, `circulants.py`, `primality_crosscheck.py`, `check_s4s4.py`,
  `check_mistake683_examples.py`, `prime_fraction.c` and their outputs;
- `.gitignore` for binaries and the 12-14 MB order-12 class lists, which `reproduce.sh` regenerates in about a minute.
