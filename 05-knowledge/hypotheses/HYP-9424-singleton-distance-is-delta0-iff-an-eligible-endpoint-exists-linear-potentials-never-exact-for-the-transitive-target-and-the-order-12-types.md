---
id: HYP-9424
title: "After THM-4818: (1) the native distance from any tournament to the singleton locus is delta^0 iff some delta^0-attaining vertex is an eligible path endpoint (small score in the initial strong component or large score in the terminal one), and delta^0 + 1 otherwise; (2) for every n >= 6 the fractional dipole LP, hence every linear negation-invariant potential, fails to be exact for the transitive target on some rigid non-self-converse class; (3) the symmetric primes of order 12 come from the automorphism types 3^2 1^6, 3^3 1^3, 3^4, 5^2 1^2 and the anti-types 2^6, 6 2^3, 6^2, 10 2, and number a few million; (4) the distance to the transitive class is attained by a word whose scores decrease L at every step except the first"
status: >
  OPEN (klein-2026-10-27). (1) "no eligible vertex => > delta^0" PROVED (THM-4818 (C2)); the converse and the value
  delta^0 + 1 are FINITE-EXACT at orders <= 8 (0, 1, 3, 16 failures, all with delta^0 = 1), so the "+1" rests on
  thin evidence; the strong case is THM-4811 (D2) / HYP-9421 (6). (2) FINITE-EXACT at n = 5 (exact except the three
  |Aut| = 3 classes) and n = 6 (gaps on rigid classes); the LP at n = 7 was not run over all labellings (5,040 per
  class). (3) a task: the candidate counts are 56*2^17, 2^24, 2^22, 2^14 and, for the anti-types with the 2-cycle
  normalisation, 56*2^21, 2^18, 2^15, 2^8 (the 2 mod 4 law of THM-4818 (A) excludes every type with a 4- or
  8-cycle); THM-4817's "type tower beyond" is reserved by mac-mini-2026-10-27-dipole. (4) a guess from the L-descent
  counts of the census (a class at distance L + k has geodesic words whose L-profile is not monotone); untested.
session: klein-2026-10-27 (codex/session-overnight-20261027)
source: 05-knowledge/results/order11_prime_types_n12_ledger_dipole_census_20261027/README.md (§3)
related:
  - 01-canon/theorems/THM-4818-anti-automorphism-cycles-are-2-mod-4-order-11-prime-types-by-symmetry-order-12-ledger-dipole-census-condensation-criterion-and-p7-class-distance-six.md
  - 05-knowledge/hypotheses/HYP-9421-black-descent-exceptions-tc-spectrum-and-the-native-forcing-tail.md
  - 05-knowledge/hypotheses/HYP-9422-regular-tournaments-exactly-m-black-moves-from-a-singleton-plus-strong-core-nuclear-norm-maximizers-off-the-hadamard-orders-and-the-order-11-types.md
---

# HYP-9424 — eligibility, the LP gap, the order-12 types

1. **Eligibility decides the singleton distance.** With `δ⁰ = min_v min(s(v), n−1−s(v))`, call `w` eligible when
   `s(w) = δ⁰` and `w` is in the initial strong component, or `n−1−s(w) = δ⁰` and `w` is in the terminal one.
   Prediction: distance `δ⁰` iff an eligible vertex exists, and `δ⁰ + 1` otherwise (one move to relocate the
   components, then THM-4811 (D)'s descent). What would refute it: a non-strong class at order 9 or 10 without an
   eligible vertex at distance `≥ δ⁰ + 2`, or one with an eligible vertex at distance `> δ⁰` (the strong version is
   HYP-9421 (6)).
2. **Linear potentials are never exact for the transitive target beyond order 5.** The fractional dipole LP of
   THM-4818 (C4) is the best bound any linear negation-invariant potential can give. Prediction: for every `n ≥ 6`
   some rigid non-self-converse class has `d(T, TT_n) > LP(T)`; the gap classes are exactly those whose LP-optimal
   decompositions use paths that are not directed paths of the intermediate tournaments. Test: run `lp_dipole.py 7`
   on the 182 `L`-failures (5,040 labellings each) and compare with `an7.txt`.
3. **Order-12 symmetric primes.** The same symmetry-first enumeration at order 12 (automorphism types with at least
   two `p`-cycles: `3^2 1^6`, `3^3 1^3`, `3^4`, `5^2 1^2`; a lone 11-cycle is a module; anti-types with no fixed
   point and all cycles `≡ 2 (mod 4)`: `2^6`, `6 2^3`, `6^2`, `10 2`) has about 1.4·10^8 candidates after the two
   reductions, an hour in C. Prediction from the Burnside defect (`P(n) − P_lab(n)/n!` grows by a factor ≈ 13.6 from
   10 to 11): about 1.6·10^5 non-rigid primes and, by the self-converse decomposition at order 12, 1,170,949
   self-converse primes (already fixed by A002785(12)); the type tower would then give the order-13 ledger with
   A000568(13) = 48,542,114,686,912 and A002785(13) = 95,458,560.
4. **Geodesics to the transitive class.** In the census every class at distance `L + k` (`k ≥ 1`) from `TT_n` has a
   geodesic along which `L` is not monotone; guess: one non-descending step suffices (`k ≤ 1` in a suitably
   generalised potential). Untested beyond the counts in `an5.txt`–`an8.txt`.
