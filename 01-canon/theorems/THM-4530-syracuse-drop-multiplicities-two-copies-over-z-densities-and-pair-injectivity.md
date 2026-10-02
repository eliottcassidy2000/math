---
id: THM-4530
title: "Syracuse drop multiplicities. Over all odd integers, every d in Z is a drop A - S(A) = 2d exactly 2 + N(6d+1) times: two unit preimages 8d+1 and -4d-1, plus one for each v >= 3 with (2^v - 3) | 6d+1. This is the exact form of the owner's 'two copies'. The densities of the multiplicities exist and are computable (delta_1 = 0.696174903440..., delta_2 = 0.266238036677..., delta_3 = 0.035389558839...). The maximal order lies between sqrt(2 log2 X) - 2.5 and 1 + log2(6X+1)/2. Two consecutive drops determine the point. Among odd q, only q = 3 has a drop map that hits every integer. Drop statistics cannot constrain cycles"
status: >
  PROVED + INDEPENDENTLY AUDITED:
  - Prop. 2.6, two copies over Z;
  - Theorem 2.1, the profinite law of m(d) and the existence of the densities, with Lemma 2.2 (absolute convergence);
  - Theorem 2.5, bounds on the maximal order;
  - Theorem 3.1, the summatory function of K, and Theorem 3.2, the Jacobsthal sum;
  - Theorem 4.1, the k-step identity;
  - Prop. 4.4, two copies at k = 2;
  - Theorem 4.5, pair-injectivity on both sheets (proof read line by line);
  - Theorem 5.2, unit branches: all d >= 0 iff q = 2^a - 1, all d < 0 iff q = 2^a + 1, all of Z only for q = 3;
  - Section 6, the mod-3^k identities.
  FINITE-EXACT + AUDITED: the densities to 18 digits (exact rationals for v <= 64, error < 1.1e-19; recomputed by the
  orchestrator's own implementation); the 12 multiplicity records up to 13 copies; the k = 2, 3 multiplicity laws.
  CITED: Terras/Everett (the joint law of drops along orbits); Rhin (the irrationality measure of log 3/log 2, used for the
  k-step mean multiplicities); P. Borwein 1991 (irrationality of sum 1/(2^v - 3)).
  CONJECTURE, now HYP-9171: max_{d <= X} m(d) = (1+o(1)) sqrt(2 log2 X), and pair-injectivity for every odd q.
  NO consequence for Collatz: the 3x-1 sheet has the same drop statistics and the same pair-injectivity, yet it has three
  cycles (Theorem 4.7).
  Builds on THM-4527 (opus S15, twentieth note), which proved 6K + 1 = (2^v - 3)S and m(d) = #{v >= 2 : (2^v - 3) | 6d+1}
  on the positive integers. This lane re-verified it independently, and the orchestrator re-verified it again.
  Notation: S(A) = oddpart(3A+1), v = v_2(3A+1), K(A) = (A - S(A))/2, M_v = 2^v - 3, N(n) = #{v >= 3 : M_v | n}, and
  m(d) = 1 + N(6d+1) for d >= 0.
source: collatz-procgen-20260922 session, cdiff lane (2026-10-01), answering the owner's "2N-1 -> 2N-1-K_N ... exactly 2 copies of each difference ... please investigate the structure yourself" and opus S15's directions D74/D75; audited and promoted by the session orchestrator 2026-10-01
depends_on:
  - 01-canon/theorems/THM-4527-syracuse-in-owner-labels-microcosm-and-difference-spectrum.md (opus S15: the label map F, 6K + 1 = (2^v - 3)S, the multiplicity formula on N)
related:
  - 01-canon/theorems/THM-4484-free-and-sporadic-cycles.md (the unit words are exactly the free cycles)
  - 05-knowledge/hypotheses/HYP-9171-drop-multiplicity-maximal-order-and-pair-injectivity.md
  - 05-knowledge/results/collatz_label_map_paley_snarks_20261001.md
note: 05-knowledge/results/procgen_cdiff_20261001_drop_multiplicities.md
scripts: 04-computation/experiments/procgen_cdiff_20261001_run.py
script_audit: 04-computation/experiments/procgen_cdiff_20261001_orchestrator_check.py
output: 05-knowledge/results/procgen_cdiff_20261001.out
output_sha256: f84405a1c0e8b2f0043e1b4b5d4927ed74c2d2bbcd5c9bf3cad70f88c8686daa
output_audit: 05-knowledge/results/procgen_cdiff_20261001_orchestrator_check.out
output_audit_sha256: 7ef2bf128e691395d358dfa922dab0dd8629b8521567c6ac9523c3d5e34d0d87
hash_basis: raw bytes
audit: >
  The orchestrator read and found sound the proofs of Prop. 2.6, Theorem 2.1 and Lemma 2.2, and Theorem 4.5.
  Theorem 4.5 was checked step by step: the elimination giving n_1 E = M_v M_v' (2^delta - 1) and E > 0; the divisibility
  E | g(2^delta - 1); the bound forcing u = 2 and e in {1, 2}; the e = 2 contradiction; the interval argument leaving only
  delta = 1, v = 3, where n_1 = 65 is not 1 mod 6.
  Independent code (procgen_cdiff_20261001_orchestrator_check.py; the lane's code was not read), 19 checks in 15 s:
  - two copies over Z for |d| <= 1500, by brute force over odd A in [-12050, 12050], and the positive-A version;
  - the densities delta_1..3 recomputed by the orchestrator's own implementation (factor M_v for v <= 64, condition on the
    valuations at the 10 shared primes 5, 11, 13, 19, 23, 29, 37, 47, 71, 431, sum independent private events): agreement
    with the lane to 3e-23, and a sieve over [0, 2e7] agreeing to 1e-6;
  - m(d_k) = k + 1 for all 12 record values, and minimality for k <= 5 by sieve;
  - both summatory identities;
  - the k-step identity on 3000 random cases;
  - the two 2-step preimages of every D in [-3000, -1], and mu_2 = 2 + [5 | D];
  - P(mu_2 = 0) and P(mu_3 = 0) on finite ranges within 1e-3 of the exact laws;
  - pair-injectivity for A < 10^6 on both sheets;
  - Theorem 5.2 for all odd q <= 33, and the 59.27% miss rate of 5x+1;
  - the mod-3^k identities.
  The lane's runner was re-run (160 s, 95 checks, ALL CHECKS PASSED); identical up to timing lines.
  Scope notes:
  - The density computation uses the same natural method in both implementations (conditioning on the shared primes).
    The agreement is between independent codes, not independent methods; the sieve is the method-independent check.
  - The minimality of the records for k >= 6 (exhaustive below 2^100) and the Rhin-based k-step asymptotics are lane-level.
---

# THM-4530 — Syracuse drop multiplicities

**PROVED / FINITE-EXACT + INDEPENDENTLY AUDITED.** Full note: [procgen_cdiff_20261001_drop_multiplicities](../../05-knowledge/results/procgen_cdiff_20261001_drop_multiplicities.md).

## 1. The owner's "two copies of each difference"

The owner wrote odd numbers as `A = 2N - 1` and their Syracuse images as `A - K`, and guessed that every difference appears exactly twice.

- **Over the positive integers (S15, THM-4527).** The multiplicity of a drop `d >= 0` is the number of `v >= 2` with `2^v - 3 | 6d + 1`. Its mean is `1.3437`, not 2.
- **Over all odd integers (this theorem).** The owner's picture is exact:
  - every `d` has exactly two **unit preimages**, `8d + 1` (branch `2^2 - 3 = 1`) and `-4d - 1` (branch `2^1 - 3 = -1`);
  - every other preimage comes from a divisor `2^v - 3 = 5, 13, 29, 61, 125, ...` of `6d + 1`.
- **On the positive integers.** One of the two unit copies is an ascent (`A = 3 mod 4`) and the other a shallow descent. The extra copies sit in S15's "microcosm" labels.

## 2. How often each multiplicity occurs

`m(d) - 1` counts the moduli `2^v - 3` that divide a random integer. These moduli are almost coprime: among `v <= 64` they share only ten primes, which makes the law exactly computable.
- **Densities.** 69.6% of drops occur once, 26.6% twice, 3.5% three times.
- **Maximal order.** The multiplicity grows at most like `log X` and at least like `sqrt(log X)`. The data and HYP-9171 point to `sqrt(2 log2 X)`.

## 3. What drops know and do not know

- **Two consecutive drops pin down the starting number.** This holds on both sheets.
- **No constraint on cycles.**
  - The k-step drop equation at `D = 0` is exactly the cycle equation `(2^p - 3^k) x = c`.
  - The `3x-1` sheet has identical drop statistics, yet it has extra cycles. So drop statistics cannot rule out cycles.
- **What is special to `q = 3`.** It is the only odd multiplier with a unit branch on both sides (`3 = 2 + 1 = 4 - 1`). Hence it is the only one whose drop map hits every integer.

**FORMALIZED 2026-10-01 (Lean round 2; package [`04-computation/lean/ProcgenSelfieEdim/`](../../04-computation/lean/ProcgenSelfieEdim/README.md); orchestrator rebuilt from scratch, verify.py PASS).**
- Prop. 2.6, two copies over Z: `dropZ_unit_pos`, `dropZ_unit_neg`, `dropZ_fibre`, `fibreZ_length` (exactly 2 + N(6d+1) preimages).
- Theorem 4.5, two consecutive drops determine the point, for 3x+1 and 3x-1: `drop_pair_injective`, `dropM_pair_injective`.
- Theorem 5.2, both directions, over positive odd A: `units_nonneg_iff`, `units_neg_iff`, `units_all_iff` (every integer is a drop iff q = 3).
