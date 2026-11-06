# Pre-registration addendum P6: the residue-resolved Poisson law

kind-pasteur-2026-11-06-S13. Committed and pushed **before** the P6 runs, together with the frozen program
`residue_law.c` (sha256 10db375b...) and seed list `p6_seeds.txt` (sha256 a57b20c5...). The theorem number in
`PREREGISTRATION.md` (THM-4887) became THM-4888 after a reservation race (commit 48ba576e1f); nothing else changes.

## The law

For a seed s, a window W of odd integers and an odd modulus M, let, for each class c (mod M):
- n_c = #{odd N in W : N = c (mod M)};
- p_c = #{m >= 0 : H_s(m) in W, H_s(m) = c (mod M)}, parameters counted with multiplicity;
- v_c = number of distinct such values; miss_c = 1 - v_c/n_c; lambda_c = p_c/n_c.

**Prediction:** miss_c = exp(-kappa lambda_c), with kappa = 1/2 for the root t (its converse involution pairs the
parameters) and kappa = 1 for every other seed. It refines the local-factor model of `PREREGISTRATION.md`: confinement
by a prime p | gcd(d, R) is the case lambda_c = 0 off the allowed classes.

## Exploration behind it (in-sample, not evidence)

- Root, [2^20, 2^21) and [2^24, 2^25), M = 3, 5, 7, 9, 16, 21: every class within 0.0023 at 2^24. The mod-7 miss
  rates span 0.44-0.58 at 2^20, 0.46-0.56 at 2^22 and 0.46-0.55 at 2^24; THM-4878 (5) recorded this bias of the
  one-apex image without an explanation. The pattern changes between windows, and the law tracks it.
- 33 in-sample generic seeds of bank A (d > 1, R != 2d, R != h, no prime >= 5 in gcd(d, R)), [2^22, 2^23),
  M = 63: all 1890 live pairs within 0.0023 (95th percentile 0.0007). Mean error 0.00024 against 0.0080 for the
  class-blind prediction; per-seed spreads of the class miss rates 0.015-0.085.
- Root, [2^22, 2^23), M = 63: max 0.0046 (class size 33k, noise-dominated).
- Not a prediction, recorded honestly: bank-level independence fails. Bank A's union misses 0.0219 of the odd N in
  [2^22, 2^23); the product of its 85 members' miss rates per class mod 3 is 0.0036. Removing the 35 balanced seeds
  nested by JZ = OJ gives 0.0157 (union unchanged), and residue or position resolution does not close the gap.

## Confirmatory predictions

- **P6a (root).** In each of [2^25, 2^26) and [2^26, 2^27), for M = 7, 9, 21 and 63 (aggregated from the mod-63
  counts) and every class c: |miss_c - exp(-lambda_c/2)| <= 0.006. Non-triviality: the mod-7 miss rates span at least
  0.05 in each window.
- **P6b (generic seeds).** The 102 seeds of `p6_seeds.txt` are the P1 seeds without an odd prime >= 5 dividing
  gcd(d, R). Window [2^24, 2^25), M = 63. Over all (seed, class) pairs with lambda_c > 0,
  |miss_c - exp(-lambda_c)| <= 0.004 for at least 95% of the pairs and <= 0.010 for all of them. Every class with
  p_c = 0 has v_c = 0 (automatic). Reported, not predicted: the mean errors of the resolved and the class-blind law.

## Commands

```bash
gcc -O2 -o residue_law residue_law.c
./residue_law root.txt 33554432 67108864 63 > p6_root_w25.txt
./residue_law root.txt 67108864 134217728 63 > p6_root_w26.txt
./residue_law p6_seeds.txt 16777216 33554432 63 > p6_seeds_w24.txt
```

A failed prediction is recorded as REFUTED in HYP-9480, not repaired after the fact.
