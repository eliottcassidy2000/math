# Audit F report: HYP-9244 and the coalescence phase-diagram note, sections 0, 1, 2, 5

Independent adversarial audit, run 2026-10-08 as a session subagent. The harness did not let the subagent write this file, so the main session saved its final message here verbatim. All scripts and outputs are in this folder.

**Target.** The pre-audit checkpoint `8de8989efb`. The author's session was editing the note in the working tree during the audit: §1.4 was rewritten around "rivers", and THM-4606's numbers were updated after audit E. The findings mark where that in-progress revision already fixes something.

**Bottom line.** Every number I re-derived independently reproduces, mostly exactly, otherwise within about 2σ. The conjecture as stated is false, because of congruence obstructions that the "non-degenerate" clause does not exclude. The sheet-blindness reading (P1) is false for the same reason. With an accessibility hypothesis added, all my data are consistent with the rank trichotomy, and the "no conspiracy" premise holds in every case I measured.

## Verdicts

| Claim | Verdict | Severity |
|---|---|---|
| HYP Setting: pair-chain recursion, coupling rule `i ≡ Mj + e`, Haar digits, mean-zero debt walk | CONFIRMED | — |
| Ranks, signs of Λ, `d \| m_i i + r_i`, branches onto `Z_d` (all 10 maps) | CONFIRMED | — |
| HYP Statement and mechanism "coalescence iff contraction and recurrence" | **WRONG as stated** (F1) | major |
| HYP evidence table, q(T) values | CONFIRMED WITH FIXES (F4, N1) | minor |
| HYP and §2 "debt visits" column | CONFIRMED for ranks 2, 3 and expanding rows; WRONG for contracting rank-1 rows (F3) | minor |
| HYP status, "PROVED parts" | CONFIRMED WITH FIXES (F5) | minor |
| Peres–Popov–Sousi (PPS) caveat | CONFIRMED (fair) WITH FIXES (F6) | minor |
| Reading: "Collatz is rank one in every clock" | CONFIRMED WITH FIXES (N11) | nit |
| Reading and §5: orphan-law 1/2 "is" the rank-one first-return exponent | NOT CONFIRMED; heuristic stated as fact (F7) | minor |
| Reading P1, §5 P1, §0 "sheet-blind" | **WRONG** (F2) | major |
| Reading and §5 P12 ("third unit → rank two") | NOT CONFIRMED; loose (N8) | nit |
| §1.1 rivers = level sets of the odd count `o`; product identity; offset; THM-4605 citation | CONFIRMED (finite-exact; data file reproduced byte-for-byte); citation fair | nits N4–N6, N12 |
| §1.2 span table | Numbers CONFIRMED; "O(√K)" NOT CONFIRMED (F8) | minor |
| §1.3 jumps, births, orphan rates, braid widths | CONFIRMED exactly | nits N3, N9 |
| §1.4 (pre-audit) landing integers and bound `N_D < (3/2)^D` | Numbers CONFIRMED exactly (different method); per-depth Kesten–Goldie reading NOT CONFIRMED (F9; fixed in the working tree) | minor |
| §1.5 anchor census | CONFIRMED exactly; the proof sketch has a wrong step count (F10) | minor |
| §2 phase diagram | CONFIRMED WITH FIXES (F3, F4) | minor |

## Findings

### F1 (major): HYP-9244 is false as written because of congruence obstructions

**Lemma.** Let p be a prime not dividing `d·Πm_i`, and suppose some `c ∈ Z_p` satisfies `(m_i − d)c + r_i ≡ 0 (mod p)` for every i. Then along the pair chain from `(1, e_0)`:

    e_n − (1 − M_n)c ≡ P_n e_0 (mod p),   P_n = Π_k m_{i_k}/d   (a p-adic unit).

The one-step identity is `e' − (1 − M')c = (m_i/d)(e − (1 − M)c)` mod p. So whenever `M_n = 1` we have `e_n ≢ 0`, the state (1, 0) is never reached, and a meeting with `M ≠ 1` is Haar-null. Hence if p ∤ e_0, then y and y + e_0 never meet.

**Counterexamples.** All are contracting, and all satisfy the HYP's non-degeneracy clause.

| Map | Rank | Obstruction |
|---|---|---|
| **3x+5** on Z_2 | 1 | p = 5, c = 0 |
| Z_3, m = (1,1,5), r = (0,2,2) | 1 | parity |
| Z_3, m = (1,2,5), r = (0,7,14) (the HYP's rank-2 multipliers) | 2 | p = 7 |
| Z_3, m = (1,2,5), r = (9,1,5) | 2 | affine, c = 1 |
| Z_5, m = (1,2,3,7,1), r ≡ 0 mod 11 | 3 | p = 11 |

**Checks.**
* The exact identity held with 0 violations in 60,000 chain steps per map.
* Direct integer orbits gave 0 meetings in 400 samples to T = 4096 for each map.
* Unobstructed controls with the same multipliers merge: 3x+1 q(4096) = 0.172; Z_3 (1,1,5) with r = (0,−1,−1) gives 0.028.
* 3x+5 with offset 5 behaves exactly like 3x+1 with offset 1 (√T q = 10.27 against 10.35 at T = 4096).
* For 3x+r the clean statement is: y and y+e merge almost surely iff `r/3^{v_3(r)}` divides e (THM-4581 3′ plus the lemma).

**Related problem.** For a start e_0 prime to d, the very first coupling `j ↦ j + e_0` is a d-cycle whose steps already span Γ⊗R. So the clause "coupling laws that occur span Γ⊗R" excludes nothing in rank ≥ 1 (see F11).

**Fix.** Add "(1, 0) is accessible from every state reachable from the start", or at least "no obstruction mod p^k as above with v_p(e_0) < k". Restrict starts to the admissible module (for px + r: e ∈ rZ[1/p]). Restate the mechanism as contraction × recurrence × accessibility. The rank-0 "reaches 0" clause becomes a special case.
* The HYP's own refutation criterion ("rank-1 contracting map with a tail exponent other than 1/2") is met by 3x+5, where q ≡ 1.
* The HYP's ten maps have no such obstruction for any p < 2000 (`obstruction_scan.out`), so its evidence is unaffected.

### F2 (major): sheet-blindness is false

The claim "Every statement here depends only on the multipliers; the r_i enter only through rank-0 connectivity", and the §0/§5 claim "sheet-blind (px + r is conjugate to px + 1)", both fail.
* The conjugacy x ↦ rx rescales the offset e ↦ re. THM-4606 statement 6 says this correctly; the note's summary drops it.
* Same multipliers, different sheets, offset 1: 3x+1 merges almost surely, 3x+5 never.
* Z_3 (1,2,5) at T = 16384, 1000–2000 samples each:

  | r | q(16384) |
  |---|---|
  | (0,1,2) | 0.421 |
  | (0,1,−1) | 0.499 |
  | (0,4,2) | 0.634 |
  | (0,−2,5) | 0.678 |
  | (0,7,14) | 1 |

* The narrower claim survives: 3x+1 and 3x−1 have identical coalescence laws (via x ↦ −x, offsets ±1).

**Fix.** Say that coalescence laws are invariant under affine conjugacy, which rescales offsets, and depend on the sheet through the offset's arithmetic and accessibility.

### F3 (minor): the visits diagnostic for contracting rank-1 rows

"~√T" for 3x+1 and Z_3 (1,2,4) is not what the author's own diagnostic gives. That diagnostic is a mean over all samples, frozen at absorption, and it saturates:

| Map | T | Mean visits |
|---|---|---|
| 3x+1 | 1024 / 4096 / 16384 / 65536 | 9.0 / 11.9 / 14.1 / 15.2 |
| Z_3 (1,2,4) | 1024 / 4096 / 16384 | 1.44 / 1.54 / 1.61 |

No .out files exist for these two rows. "Pólya predicts √T, log T or bounded" applies only to unabsorbed walks, which the expanding rows confirm:
* 5x+1: ≈ 1.08√T.
* Z_3 (1,4,16): 0.544√T.
* Z_3 (1,5,7): about log T, 2.86 at T = 4096.

**Fix.** Use "visits in (T/4, T] per chain alive at T/4":
* 3x+1 (rank 1): about 8–13.
* Z3_125 (rank 2): about 0.21–0.26, flat.
* Z5_12371 (rank 3): 0.117 → 0.001, roughly T^(−1/2).

Also: "in rank 3 it stops at 0.54" is the all-sample mean. Among unmerged chains the cumulative mean is 0.216.

### F4 (minor): the Z_3 (1,2,5) headline numbers

"0.598 → 0.338" comes from the lowest of the author's three runs (500 samples). The other runs give 0.617 / 0.658 → 0.416 / 0.430 at T = 16384. My 2000-sample run gives:

| T | 16 | 64 | 256 | 1024 | 4096 | 16384 | 65536 |
|---|---|---|---|---|---|---|---|
| q (±0.011) | 0.643 | 0.582 | 0.530 | 0.487 | 0.451 | 0.421 | 0.391 |

The qualitative claim holds strongly: 1/q rises by 0.164, 0.169, 0.164, 0.167, 0.155, 0.183 per factor of 4. Quote pooled values with standard errors.

### F5 (minor): scope of the "PROVED parts"

* "Expanding case for d = 2" is proved only for px+1 (m_0 = 1, r = (0,1)) with integer offsets, and px + r with offsets in rZ.
* "Rank zero with d = 2" is proved only for x+1. Under x+5, y and y+1 never meet.

### F6 (minor): the PPS citation

The paraphrase is accurate. I checked the arXiv PDF:
* Definition 1.1: "d-dimensional" means an invertible covariance matrix.
* Theorem 1.2: two d-dimensional mean-zero laws with 2+β moments in d ≥ 3 give transience under any adapted rule.
* Theorem 1.5: d fully supported laws chosen by a position-dependent rule (the maximal coordinate) give recurrence.

Fixes:
* Add the moment hypothesis.
* Give the journal title: "On recurrence and transience of self-interacting random walks", Bull. Braz. Math. Soc. (N.S.) 44(4) (2013) 841–867. The arXiv title is "Self-interacting random walks".
* Note that the debt walks use many laws, some degenerate. Examples: the zero law at the identity coupling, and the j ↦ 4j coupling for Z5_12371, whose covariance has rank 2. Such laws violate Theorem 1.3's trace condition, so no PPS result applies directly.
* The rank-2 claim is not a corollary either: zero-drift walks in d = 2 can be transient (Georgiou–Menshikov–Mijatović–Wade, Adv. Appl. Probab. 48A (2016) 99–118). Cite this.

### F7 (minor): the orphan-law exponent

In HYP-9242, 1/2 is the exponent of a sketch-level upper bound. The measured exponents are:
* 0.42 for K = 2, 3;
* 0.54–0.63 for K = 9–33 ("which excludes 1/2");
* 0.59 on the Mersenne line (95% CI [0.42, 0.75]).

The note's own §1.3 gives births ~X^0.41, i.e. orphan density ~X^(−0.59). Label the identification as a heuristic.

### F8 (minor): river spans

Mean span/√K by window is 1.85 ± 0.45, 2.11 ± 0.24, 2.38 ± 0.21, 2.68 ± 0.25. The rise is about 2σ, the last window is right-censored at 12800, and the mean span's local exponent is about 0.6. The data support "span/√K ≤ 5.9 for K ≤ 12800", not O(√K).

### F9 (minor): §1.4 in the pre-audit version

My reproduction uses a different method: 2-adic truncation plus the parity-vector formula. The pre-audit numbers reproduce exactly:
* all N_D positive, max 880, median 4, s_0/D ∈ [1.500, 2.882];
* n·P(N_D > n) = 1.510, 1.737, 2.789, 3.522, 3.634, 3.586, 2.497, 2.817, 3.842;
* N_D < (3/2)^D for every D ≤ 4499.

But the 3998 depths are not independent samples, so per-depth n·P is not a tail estimate:
* There are only 40 distinct N_D values.
* 3545 of 3997 consecutive depths share N_D.
* There are 186 distinct landings (N_D, s_0).
* The 30 depths with N_D > 512 come from just two rivers (880 × 8, 732 × 22).

The working-tree revision fixes this, and its numbers reproduce exactly:
* 195 classes for 3 ≤ D < 4500;
* 3929 of 4496 consecutive pairs share a landing;
* largest classes 144, 126, 124, 114, 114;
* the stated records;
* per-river P(N > 4, 16, 64, 256) = 0.467, 0.221, 0.067, 0.026.

Its "consistent with 2.87/x" is generous at 64–256: 5 rivers above 256 against 2.2 expected. The note's "excess at most 9" also checks out: the maximum of odd(N) − σ_T(N)/2 for N ≤ 880 is 8.5, at N = 871.

### F10 (minor): §1.5 proof sketch

"Two common odd steps, a common even step, then u even and v odd" is wrong. The pattern is (odd, odd), (even, even), (u even, v odd): three steps, landing on (3c+1)/8. I checked all 2048 test residues.

### F11 (minor): the non-degeneracy clause is imprecise and vacuous

* "Laws that occur" is undefined: once, infinitely often, or with positive frequency?
* As written it is automatic for starts prime to d (F1).

## Nits

* **N1.** The table has ten maps, not "nine".
* **N2.** The widest river by span is [10367, 10972] (605 exponents, 360 members). [7763, 8298] is the widest by span/√K and the largest by membership.
* **N3.** R(X)/X^0.41 declines from 3.09 to 2.82; the local exponents are 0.39, 0.39, 0.36. "≈ 2.9 K^0.41" is rough.
* **N4.** "ε ∈ (0, 0.33)" is empirical: the maximum over n ≤ 10^6 is 0.3256 at n = 993, and over M_K it is 0.3217. "Rivers = level sets of o" is finite-exact for K ≤ 12800, not exact in general.
* **N5.** "Up to 2^(−K)" should read "up to 2^(−K)/ln 2".
* **N6.** In "jump tuning is automatic", only the bound ‖Δo·log₂3‖ ≤ max ε − min ε (about 0.32) is automatic. The median 0.063 reflects the concentration of ε (interquartile range 0.16–0.27).
* **N7.** See F3: "stops at 0.54" mixes the all-sample and unmerged-chain means.
* **N8.** The debt rank is at most d − 1 and does not grow with the number of primes; d = 2 maps have rank ≤ 1 even for 15x+1. So "a third unit makes the debt rank two" is loose.
* **N9.** The orphan rate by v₂(K−1) has χ² = 7.4 on 4 degrees of freedom (p ≈ 0.11): flat, but borderline.
* **N10.** Several cited numbers have no saved .out: the 3x+1 phase and visits runs, and the Z3_124 visits run.
* **N11.** "Every clock" is fine for the Terras clock, the odd-step clock (debt group 2^Z, checked algebraically; that map has infinitely many branches, outside the HYP's class) and any linear clock. "Always T^(−1/2)" needs e ∈ Z[1/3]: y and y + 1/5 never meet under 3x+1.
* **N12.** The river key is a function of o for all K ≥ 2, not only K ≥ 5.

## Details by task

**A. Independent simulation.** I iterated actual integer orbits of y and y+1, with y uniform mod d^(T+64). The orbits advance in blocks using the composed affine map, validated bit-for-bit against naive one-step iteration. The debt is an exponent vector built from the two branch sequences.
* The pair-chain recursion matches the integer orbits with 0 mismatches over about 200k steps (10 maps, offsets 1 and 7). `B_n e_n ∈ Z` always holds, so e is an integer at every return to M = 1.
* There were no meetings with nonzero debt in any run.

| Map | Samples | Mine | Author |
|---|---|---|---|
| 3x+1 | 4000 | √T q = 2.48, 4.13, 6.44, 8.60, 10.35 (T = 4096), 10.69, 10.75 (T = 65536) | 2.4 → 9.6 |
| Z_3 (1,2,4) | 2000 | √T q 1.59 → 2.50 (T = 4096), 2.11 ± 0.36 (T = 16384) | 1.53 → 2.69 |
| Z_3 (1,2,5) | 2000 | see F4 | see F4 |
| Z_5 (1,2,3,1,1) | 1500 | 0.396 → 0.210 (T = 16384) | 0.423 → 0.243 |
| Z_5 (1,2,4,7,1) | 1500 | 0.699 → 0.487 | 0.692 → 0.478 |
| Z_5 (1,2,3,7,1) | 2000 | 0.715 (T = 16), then 0.689, 0.6875, 0.686, 0.6855 for T = 1024 … 65536; visits 0.516 | 0.674 → 0.654; 0.54 |
| 5x+1 | 5000 | 0.9930 ± 0.0012 | 0.9917 |
| Z_3 (1,4,16) | 1500 | 0.803 ± 0.010 | 0.77 |
| Z_3 (1,5,7) | 1000 | 1.000, visits 2.86 | 1.000, 2.70 |
| x+1 | 2000 | 0 by T = 16, visits 2.00 | 0 by T = 16, 2.02 |

**B. Map checks.** All ranks, Λ signs, divisibility conditions and onto-ness (exhaustive mod d^K) pass. Λ values: x+1 −0.693, 3x+1 −0.144, (1,2,4) −0.405, (1,2,5) −0.331, (1,2,3,1,1) −1.251, (1,2,4,7,1) −0.804, (1,2,3,7,1) −0.862, 5x+1 +0.112, (1,4,16) +0.288, (1,5,7) +0.087.

**C. Counterexample hunt.** Besides F1, no map I tried breaks the rank rule.
* **Rank 2, 12 maps, none transient.** The maps were Z₃ (1,2,7), (1,4,5) [Λ = −0.10] and (1,2,11) [Λ = −0.068]; Z₅ (1,3,9,2,1) and (1,6,11,1,1) [translation-only couplings]; Z₄ (1,3,5,1); Z₇ (1,2,3,4,6,1,1); a product-like Z₆ map m_i = α(i mod 2)·β(i mod 3); and three other sheets of Z₃ (1,2,5). In every case the per-survivor visit window stays flat at 0.11–0.66 and q decreases. The Z₃ (1,2,7) and (1,2,11) runs and the extra (1,2,5) sheets sit at the slow end of 1/q growth.
  * Weakly contracting rank-2 maps have q almost flat over accessible T: Z₃ (1,2,11) goes 0.757 → 0.735 from T = 1024 to 65536, while its visit window holds at about 0.55–0.66. So "q freezes" alone is not a rank-3 signature.
* **Rank 3, 6 maps, none with q → 0.** The maps were Z₅ (1,2,3,7,1), Z₄ (1,3,5,7), Z₅ (1,8,3,7,12) [Λ = −0.088], Z₅ (1,2,3,7,6), Z₇ (1,2,3,5,1,1,1), and the translation-only Z₅ (1,6,11,16,1). All have q flat by T ≈ 256–1024, at levels 0.45–0.98, with the visit window going to 0.
* **Rank 4.** Z₅ (1,2,3,7,11): q flat at 0.949.
* **Tied multipliers.** Z₅ (1,2,2,2,2) is rank 1 and behaves like it (√T q ≈ 0.9–1.4).
* **Anisotropy test (supports the "no conspiracy" premise).** For Z3_125, Z5_12311, Z5_12371, Z₄ (1,3,5,7) and Z₃ (1,2,11), the effective dimension `d_eff = E|ξ_w|²/E[(ξ_w·u)²]` equals the rank to within 0.003, overall and in each half-space (M > 1 vs M < 1, and x₀ ≷ 0). The coupling states are exactly uniform: identity share 1/6, 1/20 and 1/8 for d = 3, 5, 4. Equidistribution of (M mod d, e mod d), independent of the debt, looks like a provable route to the rank criterion once accessibility is assumed.

**D. PPS.** See F6.

**E. Battery.**
* My 2¹⁶-ary Terras recomputation of (o, σ_T) for all K ≤ 12800 matches the data file with 0 mismatches.
* The product identity holds at K = 101, 1001, 5001 to 1e-12.
* All §1.1–1.3 numbers reproduce exactly. Side observation: all 106 river births at K ≥ 200 occur at odd K.
* §1.4: see F9.
* §1.5: 82 cycles, 154 anchors (60 with first letter 2, 94 with ≥ 3), outcomes DRIFT 1450 / ANCHOR 192 / SHIFT 4222 / ABSORB 296. All 94 anchors with first letter ≥ 3 merge at D = 1 at time 3. 43 of the 60 are barriers. 7/23 merges at D = 5, 6; 7/503 at D = 15–18, 25, 26, 35, 36.
* "Never merge" in §1.5 means never absorb: the 192 ANCHOR outcomes are numeric meetings with unpaid debt.

**F. Overclaims.** See F2, F7, N8 and N11.

## Files in this folder

`maps_check.py`, `bigint_pairs.py`, `chain_vs_orbits.py`, `obstruction_check.py`, `obstruction_scan.py`, `anisotropy.py`, `battery_check.py`, `anchor_check.py` and `run_*.sh`, with outputs `maps_check.out`, `partA_bigint.out`, `chain_vs_orbits.out`, `obstruction_scan.out`, `partC_maps.out`, `partC4_extra.out`, `partC3_anisotropy.out`, `partE_battery.out`, `partE_landing.out`, `partE_landing4000.out` and `partE_anchors.out`.

## Not verified

* I did not attempt proofs of the conjecture or of THM-4606.
* T beyond 65536 was not run.
* I have not checked whether the absence of local obstructions implies accessibility.
* I did not rigorously check the working-tree revision's new "density-one" claim, which is outside scope.

Sources: arXiv:1203.3459 (https://arxiv.org/abs/1203.3459); journal BibTeX for Bull. Braz. Math. Soc. 44(4) (https://www.cmup.pt/publications-node/9763/export); GMMW, arXiv:1506.08541 (https://arxiv.org/pdf/1506.08541)
