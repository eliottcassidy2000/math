# Independent audit: THM-5156 (the Monsky nullity zero sector below 125)

Independent adversarial audit, 2026-11-13, of:
- `01-canon/theorems/THM-5156-the-monsky-nullity-zero-sector-attains-every-odd-count-below-125-except-3-7-13-19-21.md`;
- the README section "Selmer-refined spectrum (THM-5156)";
- `scripts/monsky.py`.

Every repository file was read-only except this report. Nothing was stashed, committed or pushed. **Chan–Pak and BSD
remain OPEN.** The THM claims neither.

## Verdict: PASS WITH CORRECTIONS

| object | verdict |
|---|---|
| (1) additivity of nu; the nu = 0 class is the ordinal sums of strong nu = 0 classes; its spectrum is the multiplicative closure | **PASS.** The cited proof is correct (re-derived below). The block identity was re-checked on all 14,410 non-strong classes of orders 2-9. |
| (2) Busch reduction to strong classes of order ≤ 9 | **PASS.** |
| (3) the gaps below 125 are exactly 3, 7, 13, 19, 21 | **PASS.** Reproduced with independent code on independently regenerated class lists. The same five gaps hold up to 2657 (M1). |
| (4) nu distributions, n = 3..9 | **PASS**, exact. They count isomorphism classes, not merged classes (m4). |
| (5) Monsky/Selmer dictionary | **PASS WITH CORRECTIONS.** Accurate in substance and not overclaimed. The section pointer is wrong, nu is said to be defined by the census, and the first in-repo record is uncredited (m1). One caveat should be added (m2). |
| "The next undecided value is 127" (README); "127 is undecided" (THM status) | **FALSE (M1).** Six strong nu = 0 classes of order 8 have H = 127. |
| converse invariance; the diagonal convention | **PASS.** Phi(M(T^op)) = Phi(M(T)) holds as a labelled identity. All four parity-diagonal conventions give the same nu; only the zero diagonal differs. |
| credits | **MINOR corrections** (m3). |
| replay | **PASS WITH CORRECTIONS.** `monsky.py` reproduces every README number. The 1,007-class additivity check is not shipped, and classgen is not located (m4). |

**In one paragraph.**
- **Inputs.** The class lists were checked with an own canonical form. They were also regenerated independently, as all
  one-vertex extensions of the previous order, which reproduce the author's canonical sets exactly at every order
  2-9.
- **Recomputation.** H, strongness, |Aut| and the converse were recomputed for all 198,948 classes, and nu by two
  unrelated methods. Every number in (1)-(4) and every `monsky.py` output reproduces.
- **Proofs.** The additivity proof and the Busch reduction are correct.
- **The one false statement** is the closing remark about 127. In fact every odd value from 125 to 2655 is a nu = 0
  count, so (3) holds verbatim with 2657 in place of 125.
- **Dictionary.** It is used correctly. It was checked against the literal 2t × 2t Monsky matrix built from Legendre
  symbols, and on explicit prime realizations of the counts 3, 13 and 19.

---

## MAJOR finding

### M1. 127 is decided: every odd value from 125 to 2655 is attained, and the first undecided value is 2657

**Exact sentences.**
- README, "Result": "The next undecided value is 127."
- THM-5156 status: "NOT CLAIMED: values from 125 up (127 is undecided)".

**The problem.** `monsky.py` looks only below 125, so 127 was never examined. The author's own `cls_7..9.txt` contain
these strong classes with H = 127:

| H = 127 | order 7 | order 8 | order 9 |
|---|---:|---:|---:|
| nu = 0 | 0 | **6** | 3 |
| nu = 2 | 1 | 2 | 0 |

One realizer is the order-8 class with code 134554753 (re-verified in pure Python: strong, H = 127, nu = 0).

More is true.
- Every odd value in [125, 2655] is the count of an ordinal sum of strong nu = 0 classes of order ≤ 9. For example, 125
  is the count of the strong order-7 class with code 1259590.
- The first odd value that orders ≤ 9 leave undecided is **2657**. It is prime, and its 14 strong classes of order 9
  all have nu = 2. Strong classes of orders 10-15 can still reach it, since f(16) = 3125 > 2657.
- The next undecided values are 2699, 2869, 2877, 2887, 2897, ...

**Consequence: a free strengthening.** The odd numbers below 2657 that are the Hamiltonian-path count of no nu = 0
tournament, at any order, are exactly 3, 7, 13, 19, 21.
- Below 125 this is the THM's own argument.
- On [125, 2655] every value has an explicit realizer, so no Busch argument is needed there.

**Fix.**
- README: replace the sentence with "Every odd value from 125 to 2655 is attained by an ordinal sum of strong nu = 0
  classes of order ≤ 9; 127 is attained by six strong classes of order 8. The first value that orders ≤ 9 do not
  decide is 2657, a prime whose 14 strong order-9 classes all have nu = 2."
- THM status: "NOT CLAIMED: values from 2657 up (2657 is undecided)".
- Optionally restate (3) as "below 2657", which is FINITE-EXACT on the same inputs, and adjust the title. Otherwise keep
  125 in the title and add the extension as a remark.

(1)-(4) are unaffected.

---

## MINOR findings

### m1. Dictionary: the section, where nu is defined, and the first in-repo record

**Exact sentences.**
- THM status: "the nu-to-Selmer dictionary is Monsky's formula as recorded in 05-knowledge/results/integer32_twoedge_20261113
  (bridges/RESULT.md B2), whose census also defines nu".
- README: "It is the BSD bridge of the integer-32 census (integer32_twoedge_20261113, bridges B2)".

**The problems.**
- **B1, not B2.** B1 contains both the definition of nu (`M(T) = A + diag(outdegree mod 2)`, Phi = X² + X + I) and the
  dictionary:
  - orientation i → j iff p_j is a nonresidue mod p_i;
  - full Monsky matrix [[M+I, I], [I, M]], with kernel map x ↦ (x, (M+I)x);
  - "nu is the pure 2-Selmer dimension, dim Sel_2 − dim E_N(Q)[2]".

  B2 is the full-Q/scores hostile. It restates the ordinal-join identity and gives the T8/T10 prime realizations.
- **The census does not define nu.** `integer32_twoedge_20261113/census/` mentions "marked Selmer nullity" once, and only
  to point at the bridges.
- **First in-repo record.** The first record of the dictionary is `bsd_selmer_growth_20261013a.md` §§2-3 (2026-10-13).
  It is restated in `bsd_component_recursion_20261013c.md` §2, which the THM cites for additivity only, and in THM-5125
  (4)(c).
- **Literature source.** It is Monsky's matrix (his appendix to Heath-Brown's 1994 congruent-number Selmer paper), which
  B1 cites through Watkins et al. §4.1.
- **Attribution.** "Monsky's dictionary" attributes the tournament reading to Monsky. His result is the matrix formula.
  The Legendre-tournament reading is the repository's: the canon name is "Legendre tournament" (THM-5125 (1),
  MISTAKE-752).

**Fix.**
- Write: "Monsky's formula, read on the Legendre tournament (cited via bridges/RESULT.md **B1**; first recorded in
  bsd_selmer_growth_20261013a §§2-3)".
- Cite B1 for the definition and dictionary, and B2 for the join identity.
- Drop "whose census also defines nu".
- Add `bsd_selmer_growth_20261013a.md` to `related`.

**Checked.**
- **Specialization.** For an odd product of primes ≡ 3 mod 8, Monsky's matrix is [[A + D_2, D_2], [D_2, A + D_{−2}]],
  with D_2 = I and D_{−2} = 0. Its kernel is {(x, (A+I)x) : Phi(A)x = 0}.
- **Numerical check.** The auditor built the full 2t × 2t matrix from literal symbols on 3,000 random sets of 1-8
  distinct primes ≡ 3 mod 8. The symbols came from Euler's criterion, with D_2 and D_{−2} computed rather than assumed.
  The nullity equalled nu of the Legendre tournament every time, and reciprocity made every orientation a tournament.

### m2. "Force nontrivial pure 2-Selmer": state the exact dimension and the Kummer caveat

The sentence is correct and not overclaimed. Two sharpenings come free.

- **Exactly 2.** Every tournament with H ∈ {3, 13, 19} has nu = 2 exactly.
  - These counts are prime, and Busch confines their strong carriers: C3 for 3; the order-5 class 305 for 13; the
    order-6 classes 17057 and 6658 for 19.
  - Each carrier has nu = 2, and every other component is a singleton.
- **Selmer, not rank.** B1 records that nu = rank E_N(Q) + dim Sha(E_N)[2], so these counts force Selmer, not rank.
  - Concrete case: the H = 3 realization N = 627 = 3·11·19 has nu = 2.
  - Its Tunnell counts are A = 48 and B = 40, so A ≠ 2B and L(E_627, 1) ≠ 0.
  - Hence the rank is 0 (Coates–Wiles), and both dimensions are Sha[2]. This agrees with
    bsd_tunnell_terminals_20261013c.
- **Not vacuous.** The audit realized all four carrier classes on primes ≡ 3 mod 8 ("Dictionary data" below). Each has
  full Monsky nullity 2, so the statement has content for every one of the three counts.

**Fix.** "... the counts 3, 13, 19 force pure 2-Selmer dimension exactly 2 (that is, rank + dim Sha[2]; it does not
force positive rank, as N = 627 shows)".

### m3. Credits

- **Precursor of the 13 case.** opus-2026-11-10-S90 recorded that the order-5 tournament with out-masks [18, 8, 3, 5, 14]
  has H = 13 and positive Monsky nullity. It was the counterexample to "nu > 0 ⇒ 3 | H" (HYP-9657 (C)2;
  s90_rigidity_redei_20261110/README.md §4). The audit checked that it is class 305, with nu = 2. Credit it.
- **Unsourced credit line.** "This session's sweep pre-checked 3, 13, 19 through order 6 with sampling at order 7"
  stands under CREDITED but names no session or packet. Name it, or move it out of CREDITED.
- **Multiplicativity of H.** (1) uses that H is multiplicative over ordinal sums but does not cite it. THM-5152 cites
  THM-3166 for this.
- **Converse invariance.** If nu is ever read on the packet's merged-class tables, cite THM-5125 (4)(c) and B1.
- **Novelty.** The canon and results were grepped for Monsky, Selmer and nullity together with spectrum and path-count
  terms, and for the counts 13 and 19 with nu. No prior nu = 0 spectrum was found. "The decision to 125 is new here"
  stands.

### m4. Replay and labelling

- **Replay works.** `monsky.py` was run on the author's `cls_1..9.txt`. It prints exactly the README distributions, the
  strong values below 125 and the gap list (auditor `replay_monsky.out`).
- **Unshipped check.**
  - "Re-verified on all 1,007 non-strong classes of orders 2-8" is produced by no shipped script; `monsky.py` computes
    no strong components.
  - The count is right: 1 + 1 + 3 + 6 + 21 + 103 + 872 = 1,007.
  - Ship a check, or cite this audit's.
- **Inputs.** classgen's location is still not stated, as the THM-5155 audit (m5) also noted. It is
  `05-knowledge/results/first_creators_genealogy_20261107/classgen.c`. The THM's replay line should name it too.
- **Which classes.** The nu table counts isomorphism classes (A000568: 2, 4, 12, 56, 456, 6,880, 191,536), while the
  packet's other tables count merged classes. Label the row "isomorphism classes". The merged distribution is given
  below; nu is converse-invariant, so it is well defined.

---

## Proof check

**(1) Additivity.**
- **Set-up.** Let T = T1 ⇒ T2 with b = |T2|, and let B, C be the row-zero matrices of T1 and T2. Every vertex of T1
  gains b out-arcs, so

      M(T) = [ B + bI   J ]
             [   0      C ].

- **The off-diagonal block vanishes.**
  - BJ = 0, because the rows of B sum to 0.
  - Column v of C sums to (b − 1 − s_v) + (s_v mod 2) ≡ b − 1, so JC = (b − 1)J.
  - The upper-right block of Phi(M(T)) is therefore (B + bI)J + JC + J = (b + (b − 1) + 1)J = 0.
- **The diagonal blocks are unchanged.** They are Phi(B + bI) = Phi(B), because Phi(X + I) = Phi(X) over F2, and
  Phi(C).
- **Conclusion.** Phi(M(T)) = Phi(B) ⊕ Phi(C) as matrices (bsd_component_recursion (2); integer32 B2), so nu adds.
  - Iterating along the strong-component chain gives nu(T) = Σ nu(C_i) ≥ 0, so nu(T) = 0 iff every component has
    nu = 0.
  - With H(T) = Π H(C_i), the nu = 0 counts are exactly the multiplicative monoid generated by the strong nu = 0 counts,
    including 1 = H(TT1).

No gap.

**(2) Busch.**
- THM-4744 gives f(3q) = 3·5^(q−1), f(3q+1) = 5^q and f(3q+2) = 9·5^(q−1).
- Since 3 ≥ 5^(2/3) and 9 ≥ 5^(4/3), f(n) ≥ 5^((n−1)/3) ≥ 125 for every n ≥ 10.
- So every strong component of a nu = 0 tournament with H < 125 has order ≤ 9, H < 125 and nu = 0.
- The data agree: the least strong H at orders 3..9 is 3, 5, 9, 15, 25, 45, 75 = f(n).

No gap.

**Converse invariance.**
- M(T^op) = M(T)^t + (n−1)I, and M^t = M + J + I.
- With MJ = 0, JM = (n−1)J and J² = nJ, Phi(M + J) = Phi(M) + ((n−1) + n + 1)J = Phi(M).
- Hence Phi(M(T^op)) = Phi(M(T)) on the same labels (B1; THM-5125 (4)(c) for nu). Merged-class statements are
  therefore fine.
- Checked as a labelled matrix identity on all 198,948 classes.

**Edge cases of the convention.**
- The diagonal is forced: M_ii = s_i mod 2.
- Three other conventions give the same Phi up to transpose, hence the same nu. Checked on every class.
  - "Rows sum to 1" gives M + I.
  - In-adjacency with rows summing to 0 gives M(T^op).
  - "Columns sum to 0" gives M + (n−1)I.
- Only the zero diagonal (pure adjacency) is a different invariant. It differs from nu on 3, 19, 194, 3,154 and 86,361
  classes at n = 5..9. HYP-9657 (C)2 records that its behaviour differs. The THM uses the right matrix.
- Small cases: n = 1 gives Phi = [1], so nu(TT1) = 0. n = 2 gives Phi = I.
- nu is always even: ker Phi(M) is a module over F2[x]/(x² + x + 1) = F4. Every distribution agrees.

---

## Recomputed independently

All auditor code is new and shares nothing with the packet or classgen. It lives in
`C:\Users\Eliott\AppData\Local\Temp\claude\2f1de76a73e81b8a\a7c1f371-6162-4db2-99bd-68121f2beaaf\scratchpad\audit30b\`
and was not committed:
- `nu_audit.c` (`run9.out`);
- `closure.py` and `beyond.py`, a divisor recursion with a least-order DP rather than the packet's fixpoint loop;
- `witness_check.py`, all pure Python: dense matrices, own H DP, explicit ordinal sums, literal Monsky matrices and
  Tunnell counts;
- `nu_sets.py`, `merged_dist.py` and `labelled_check.py`.

| check | result |
|---|---|
| class lists, n = 1..9 | Own canonical forms (colour refinement seeded by scores, then the exhaustive minimum code within cells; \|Aut\| = the number of minimal leaves) are pairwise distinct. Counts are 1, 1, 2, 4, 12, 56, 456, 6,880, 191,536. **Independent generation:** the one-vertex extensions of the (n−1)-list (1,761,280 at n = 9) give exactly the canonical set of the n-list for n = 2..9. H (subset DP), strongness (reachability), \|Aut\|, the SC flag and the converse code have 0 mismatches over 198,948 classes. |
| nu by two methods | (A) bitmask rows with basis insertion; (B) integer matrices, an explicit square and a brute-force kernel count over all 2^n vectors. They agree on every class. A third, dense pure-Python implementation re-checked the 1,321 strong representatives used for the 2657 extension. |
| (4) | All-class distributions exactly as in the README. Strong: n = 5: 0:3, 2:2, 4:1. n = 6: 0:20, 2:14, 4:1. n = 7: 0:192, 2:135, 4:24, 6:2. n = 8: 0:3,401, 2:2,270, 4:327, 6:10. n = 9: 0:100,695, 2:67,628, 4:9,451, 6:353, 8:6. The labelled sums Σ n!/\|Aut\| reproduce the older bsd_component_recursion census (n = 5: 624 + 360t² + 40t⁴, strong 264/240/40; n = 6 strong 13,920/7,680/720). |
| additivity and block identity | On 14,410 non-strong classes (1,007 at orders 2-8 and 13,403 at order 9), Phi(M(T)) is zero off the component blocks and equals Phi(M(C)) on each block. 0 failures. |
| converse | nu(T^op) = nu(T), and the labelled Phi matrices are equal, on all classes. |
| strong nu = 0 values below 125 | 5, 9, 11, 15, 17 and every odd number from 23 to 123: 56 values, as in `monsky.py`. |
| gaps below 125 | {3, 7, 13, 19, 21}. 3: the only strong class is C3, with nu = 2. 13: one strong class (order 5, code 305), nu = 2. 19: two strong classes (order 6, codes 17057 and 6658), both nu = 2. 7 and 21 are not path counts (THM-1370). |
| least orders | Every attained value below 125 is attained at order ≤ 9. The only composite least-order witness is 75 = 5·15, and 75 is also a strong order-9 value. |
| beyond 125 | See M1. 117 explicit ordinal sums were built and their H and nu checked directly, with 0 failures; these were every value below 125 and a random sample above. |
| packet replay | `monsky.py` output is identical to the README. |

### Dictionary data (m1, m2)

**Realizations of the gap counts on distinct primes ≡ 3 mod 8.** Each N realizes the carrier class as its Legendre
tournament (verified), and each has full Monsky nullity 2.

| H | class | primes | N |
|---|---|---|---:|
| 3 | C3 | 3, 11, 19 | 627 |
| 13 | 305 | 3, 11, 19, 43, 67 | 1,806,387 |
| 19 | 17057 | 3, 11, 19, 67, 283, 619 | 7,359,010,593 |
| 19 | 6658 | 3, 19, 11, 43, 163, 2251 | 9,892,341,393 |

**BSD-conditional sanity check (Tunnell).** It found no inconsistency.
- All 66 two-prime N, all of which have nu = 0, have L(E_N, 1) ≠ 0.
- The same holds for all 69 three-prime N ≤ 150,000 with nu = 0.
- Of the 34 three-prime N with nu = 2, 18 have L = 0 and 16 have L ≠ 0, matching B1's remark that nu = 2 does not say
  which term supplies the dimensions.

### New data offered to the author

1. **The extension to 2657** (M1).
2. **The exact set of nu for each odd H < 125, at any order** (`nu_sets.out`). It is valid because Busch confines every
   component to order ≤ 9, H multiplies and nu adds.
   - **Forced nu = 0:** H ∈ {1, 5, 11, 23, 35, 43, 47, 63}. Under the dictionary, any N of this family whose Legendre
     tournament has one of these counts has trivial pure 2-Selmer, hence rank 0 and Sha[2] = 0.
   - **Forced nu = 2:** H ∈ {3, 13, 19}.
   - **Largest value:** nu = 8 at H = 81, e.g. by C3 ⇒ C3 ⇒ C3 ⇒ C3.
3. **Merged-class distribution** (nu is converse-invariant):

   | n | merged classes by nu |
   |---|---|
   | 5 | 0:5, 2:4, 4:1 |
   | 6 | 0:19, 2:12, 4:3 |
   | 7 | 0:142, 2:107, 4:21, 6:2 |
   | 8 | 0:1,996, 2:1,316, 4:208, 6:8 |
   | 9 | 0:54,692, 2:36,922, 4:5,309, 6:216, 8:5 |

   The totals are 10, 34, 272, 3,528 and 97,144, the packet's merged counts.
