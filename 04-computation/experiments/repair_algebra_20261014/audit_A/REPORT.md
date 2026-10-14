# Audit A — mac-mini-2026-10-14-wilson (THM-4640, THM-4641, HYP-9270, results note)

Independent audit agent, 2026-10-14. Audited commit `611ad440b3` (checkpoint 3), including the coordinator's three scope additions. The audit edited no audited file and committed nothing. The agent could not write this file itself, so the coordinator saved its final report here verbatim. Its scripts and their `.out` files are in this folder.

## Verdicts

| file | verdict |
|---|---|
| THM-4640 | **CORRECT WITH FIXES.** Statements 1–7 hold, and I reproduced every number with my own code. Statement 8's numbers are exact, but its definition doesn't describe them, and its growth reading is wrong. The Reading has three false sentences (M1, M3, M5), and 5(d) has one false example sentence (m1). |
| THM-4641 | **CORRECT WITH FIXES.** Every numerical row checks out. Row 2 and the title call the carries "the generators", which contradicts THM-4640 (5d) (M4). The row-8 chance baseline is about 80× too small (m7). |
| HYP-9270 | **INCORRECT as stated.** Part (b) is refuted by an elementary bound, `c(E) ≤ 3/8` for every E. Its open question is answered by THM-4581: `c(E) → 3/8`, and → 1/2 once the root switch is counted. Part (a)'s constant and exponent are not supported by the census to E = 29. |
| Results note | **CORRECT WITH FIXES.** It inherits M1–M5. §3 has a mis-citation, a mis-attribution of the cycle bound and an outdated JC2 range. The new "anchor obstruction" paragraph is correct. |

## Major findings

**M1. `c(E)` silently leaves out the root collision.**
- **Where:** THM-4640 (8), its title ("4.9% of odd 2-adic sources"), its Reading, note §1 and HYP-9270.
- **What is wrong:** `c(E)` is defined as coverage "from a primitive collision with exponents ≤ E, root twins included". The root `(a) ~ (2, a−2)` is such a collision; THM-4640 (2) itself says "root included". Its switch is THM-4555's reset switch, `m = (n−1)/2`, which covers exactly 1/4 of odd sources. `switch_coverage.py` skips `k == −1`, so the printed numbers are coverage by non-root collisions only.
- **The canonical part adds nothing.** Every canonical family starts with a letter ≥ 3 and needs `r ≥ D ≥ 2`, so it lies inside the reset switch's domain.
- **Evidence (`a7_coverage`).** An independent exact computation reproduces the authors' `c(E)` to all 6 printed digits for E = 8..22, and splits it:

  | E | c(E) as printed | canonical part | twin part (reset-2 sources) | twin share of reset-2 mass | with root (1/4 + twin) |
  |---|---|---|---|---|---|
  | 16 | 0.037497 | 0.012499 | 0.024998 | 10.0% | 0.274998 |
  | 22 | 0.048872 | 0.016291 | 0.032582 | 13.0% | 0.282582 |
  | 24 | 0.052095 | 0.017365 | 0.034730 | 13.9% | 0.284730 |

- **Exact identity.** twin = 2 × canonical, so `c = 3 × canonical`. The map `s ↦ (1)·s` keeps head totals, prefix relations and `cmin`, and lowers the switch distance by one.
- **Monte Carlo on integers (`a7b`).** On 200,000 random odd 400-bit n, the observed frequencies match the exact values:

  | part | observed | exact |
  |---|---|---|
  | canonical | 0.01214 ± 0.00024 | 0.012499 |
  | twin | 0.02555 ± 0.00035 | 0.024998 |
  | root | 0.2482 ± 0.0010 | 0.25 |
  | any | 0.2737 ± 0.0010 | 0.274998 |

  The merge was confirmed on all 54,747 covered integers.
- **Fix:** report three numbers: the reset switch (1/4), the twin coverage of reset-2 sources (13.0% of that mass at E = 22), and the total (0.283 at E = 22).

**M2. HYP-9270 (b) is false, and its open question is already answered.**
- **The bound.** Canonical pairs have D even (4b) and non-zero (m5), so D ≥ 2. Their cylinders lie in {r ≥ 2, first letter ≥ 3}, mass 1/8. Twin cylinders lie in {r ≥ 1, first letter = 2}, mass 1/4. So `c(E) ≤ 3/8` for every E.
- **Why (b) fails.** Increments `~ κ/E` would make `c(E)` unbounded. Sources with r = 1 and first letter ≥ 3 (mass 1/8) can only be switched by the excluded root.
- **The limit.** THM-4581 6(e) (PROVED, audited) says almost every first-reset-2 source merges with `(n−1)/2` at equal odd-step time. THM-4581 (3) at the start (2, 8) gives the same for `(n+1)/4 − 1`. Every such absorption is exactly a THM-4555 collision: n makes D more halvings, so (*) gives `Q_w = Q_w'`. By parity it is a twin (D = 1) or canonical (D = 2) family.
- **Consequence.** `c(E) → 3/8`, and including the root the coverage → 1/2, i.e. all trailing-ones sources.
- **The HYP's own test is already drifting.** `E·Δ` at the step midpoints, using `c(24)`:

  | midpoint E | 9 | 11 | 13 | 15 | 17 | 19 | 21 | 23 |
  |---|---|---|---|---|---|---|---|---|
  | E·Δ | 0.0791 | 0.0745 | 0.0713 | 0.0710 | 0.0708 | 0.0713 | 0.0726 | 0.0741 |

- **Fix:** retire (b) and state the bound and limit instead. Keep the log fit only as a description up to E ≈ 24.

**M3. "The −1 mechanism alone cannot descend most trailing-ones sources" is false.**
- **Where:** THM-4640 Reading, note §1, and note §3 ("union grows like log (this note)").
- **Why:** with the reset switch, coverage at E = 22 is 0.2826, which is 56.5% of the trailing-ones mass. It tends to 100% (M2).

**M4. "Generators" contradicts 5(d).**
- **Where:** THM-4640 5(a) ("So the generators are indexed by {2,4,8}"), THM-4641 title, row 2 and Reading, note §0 and §1.
- **Evidence (`a5_carries`).** No pair is derivable by carries: 0/133 with 0–3 extra leading 1s on both sides, exponent floor 0 or −6, and no BFS state cap (largest component has 85 states).
- **Fix:** call them "the two-position relations, plus the padding relation `(1)·s ~ s`".

**M5. "Higher, non-degenerate {2,3}-unit relations" is false for about 40% of primitive collisions.**
- **Where:** THM-4640 Reading; note §3 says "no generating pattern".
- **Counts (`a3_reps`, `a3b_degenerate`).**
  - Exponents ≤ 18: 737 of 1,749 primitive pairs are degenerate.
  - Exponents ≤ 16: 203 of 474 are degenerate, and 193 of those contain the 3-term Gersonides relation `4 − 3 − 1 = 0` as a vanishing sub-sum.
- **First example:** 259, `(2,8) ~ (2,3,5,6)`, splits as `18·(3+1−4) + 4·(64+3−27−24−16) = 0`. 125 is non-degenerate.
- Statement 6 is worded correctly. 5(d) remains true only in the rewriting sense.

## Minor findings

- **m1.** 5(d) says "(1,1,7) admits no carry at all". This is false: the Catalan carry with b = 1 sends `(1,1,7)` to `(4,1,1)` (val 152 for both). The authors' own `carry_derivation.out` shows a component of size 2.
- **m2.** The "6-term relation `27 = 3 + 8 + 16`" has 4 terms. The 6-term relation is `2^7 + 27 = 3 + 2^5 + 3·2^4 + 9·2^3`.
- **m3.** Statement 7's "533 primitive pairs" means 4 × 133 + 1: the 4 twin variants of each canonical pair, plus the root. 5(d) says 133. State 133 canonical pairs.
- **m4.** My independent census matches through E = 24 and extends to 29:

  | E | N(E) | N√E/2^E | local exponent |
  |---|---|---|---|
  | 25 | 172,073 | 0.025641 | 0.507 |
  | 26 | 337,478 | 0.025642 | 0.499 |
  | 27 | 662,823 | 0.025661 | 0.481 |
  | 28 | 1,303,189 | 0.025689 | 0.470 |
  | 29 | 2,563,623 | 0.025715 | 0.471 |

  C ≈ 0.0255 is below every observed value, and the exponent 1/2 is not established. Drop the law from THM-4640's title.
- **m5.** Two lemmas are missing from the proofs:
  - (i) Distinct representations of the same K have different lengths: compare the 2-adic valuation at the first differing index. Statement 3 needs this for D ≥ 1, and the canonical pairs need it for D ≥ 2. I checked it on all 524,288 canonical s with exponents ≤ 20.
  - (ii) An empty head pairs only with `(2)` at i = −1. Both the parity proof and the converse range `c ≥ max(1, 1−2i)` need this.
- **m6.** Statement 1's attribution is incomplete.
  - The formula and the collision criterion are already in THM-4555's Remarks, and in its 2026-10-07 Update as the gap lemma from openai/math #004.
  - The classical source is Böhm–Sontacchi (1978) and Lagarias (1985), with `K = (B_w − 3^k)/2`.
  - Related: `val(s)` is a 3-smooth representation in the sense of Blecksmith–McCallum–Selfridge (Amer. Math. Monthly 105 (1998)).
- **m7.** THM-4641 row 8: the "0.1% (0.104³)" chance applies a global average to small primes. Using size-dependent rates for 5, 13 and 563 gives about 0.09: Hardy–Littlewood rates 0.82, 0.51, 0.21; empirically 3/4 for p ≤ 13 and 0.16 near 563. The COINCIDENCE typing is unaffected. The Carmichael numbers below 10^6 with c + 2 prime are 561, 2465 and 656601.
- **m8.** "The next eight primes of that form are not Wilson": there are nine with a ≤ 24, and ten with a ≤ 44 (adding `2^29 − 3`).
- **m9.** Row 2's "anchors −1, 1, 1/5 for {2,4,8}" is numerology. The Gersonides object attached to 8 is the 2-cycle {−5, −7} of row 2b, not `c_3 = 1/5`.
- **m10.** Note §3 has two citation errors.
  - THM-4606 and THM-4607 are non-coalescence results, so they don't support "coalescence holds almost surely".
  - The `k ≥ 72,057,431,991` cycle bound comes from Eliahou's continued-fraction budget plus Barina's verification (THM-4633), not from Baker's theory.
- **m11.** The planar JC range has moved beyond Moh's 100: 104 (Nguyen, arXiv:1902.05923) and 108 (Guccione et al., arXiv:2204.14178).
- **m12.** Thakur's infinitude holds for every q; only the `P'' = 0` criterion is restricted to odd characteristic. Because that criterion is algebraic, the function-field result does not support a `log log x` count. For example, over F_p it rules out every degree d with 2 ≤ d < p.

**Nits:**
- "Factors 1.91–1.96" holds only from E = 18 on.
- HYP-9270's E·Δ list was computed from rounded values.
- Statement 6 should read `±2^a 3^b`.
- Row 9's "exponents < 64" undersells the result: positivity makes it complete.
- The anchor row contains the anchors' reciprocals, not the anchors.
- Cite Miller (1903) for row 3.
- The R2b census does not assert that every orbit closed; my census does, and confirms it.

## What I verified and how (all independent, exact arithmetic)

- **THM-4640 statements 1–3** (`a1_words`): Fraction iteration of all 524,287 reduced words of total ≤ 20. This covers the numerator form, cancellation, ladders, the converse with c ≤ 25, the root twins, and the 300/15,073/14 counts.
- **Positivity lemma and anchor rigidity** (`a3_reps`): the lemma holds on all canonical s with exponents ≤ 20. Searching to a = 40 finds second representations only at a = 7, 13, 17, 21, 28, 31.
- **Statements 4–5:** the mod-8 rule holds. Brute force over lengths 2–5 and exponents 0–12 finds only the unit and Catalan carries. `3^D(2^x − 1) = 2^y − 1` has no other solutions for D ≤ 200, x ≤ 1200. I also checked the Zsigmondy proof by hand, and it is correct.
- **Census and gaps** match exactly (`a2*`, `a3_reps`).
- **Switches on integers** (`a6`): 7,182 sources across 2,394 families, including the twin switches at distance D−1 and r > D. The authors had not checked twin switches on integers.
- **THM-4641** (`a8_wilson`): every row reproduced. `T(p)` was checked for p ≤ 61 with `A_p` summed directly over partitions; the authors used a residual for p < 48. Row 2b's cycles, the −17 cycle with all 7 rotation numerators divisible by 139, and THM-307's table all match.
- **Coordinator additions:** the positivity proof is correct. The reworded "no trailing one / THM-4601" bullet is correct. "Cannot descend most trailing-ones sources" is false (M3). The anchor obstruction is correct. The Hadwiger–Nelson citations don't contradict anything I know, and I didn't re-verify them, per the coordinator's instruction.

## Prior art

- I found no source for the general collision switch or the ladder classification at −1.
- The ingredients are classical: Böhm–Sontacchi, Lagarias, Zsigmondy, Gersonides, and the S-unit theorem (cited correctly).
- Costa–Gerbicz–Harvey's 2·10^13 is confirmed, and I found no newer search.

## Not checked

- Full texts of the cited papers (no downloads).
- THM-4581's proof itself: M2's limit relies on its audited PROVED status.
- N(E) beyond E = 29 (memory).

## Sources

- Thakur, ANT 7 (2013): https://msp.org/ant/2013/7-8/p03.xhtml
- Sauerberg–Shu–Thakur–Todd, Acta Arith.: https://impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/157/1/83829/infinitude-of-wilson-primes-for-mathbb-f-q-t
- Guccione et al., arXiv:2204.14178: https://arxiv.org/abs/2204.14178
- Nguyen, arXiv:1902.05923: https://arxiv.org/html/1902.05923v5
- Costa–Gerbicz–Harvey, arXiv:1209.3436: https://arxiv.org/abs/1209.3436
- Math. Comp. version: https://www.ams.org/mcom/2014-83-290/S0025-5718-2014-02800-7/
- Goodwin, arXiv:1504.03040: https://arxiv.org/abs/1504.03040
- Burson, arXiv:2005.09456: https://arxiv.org/abs/2005.09456
- Williams, arXiv:2607.01718: https://arxiv.org/abs/2607.01718
- OEIS A226383 (BMS citation): https://oeis.org/A226383
