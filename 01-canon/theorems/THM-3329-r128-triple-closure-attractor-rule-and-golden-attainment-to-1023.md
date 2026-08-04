---
id: THM-3329
title: "AMM 12592: R=128 closes three independent ways, a deterministic attractor rule with explicit master equation closes epochs through R=512 with o(R) slack, and the golden rate is attained for all n <= 1023"
status: >
  PROVED + VERIFIED-EXACT / AWAITING CROSS-SESSION HOSTILE AUDIT.
  (W) WITNESSES. The gamma* floor profile d_i = floor(gamma*(R+i)), D0 = 0,
  closes at R = 128 by THREE independently produced witnesses: (W1) direct
  l1deg beam (width 1000, two new exact prunes: recursive coefficient caps
  cap_i = B_{i+1} + shift(cap_{i+1}) and the evaluation bound
  |sigma_i(1)| <= R-1-i) finished by an exact 8-row banded attractor repair;
  (W2) LP-clamp + nearest-parity rounding with an evaluation-at-1 ballot-cone
  cut (see (E)); (W3) the parameter-free attractor rule (see (R)). All three,
  plus a chip-fire-slimmed R=64 witness, PASS a fresh-implementation referee
  (no solver code imported; exact Fibonacci/Lucas integer floor tests
  5^d < phi^(2m) < 5^(d+1); 12/12 negative controls behave). Effective rate
  58/97 < gamma* at every one. With THM-3302 this attains
  T(n) = n + 1 + floor(gamma* n) for all n <= 255.
  (R) THE RULE. Delta_i := AdmClamp_{d_i}(trunc_{d_i}(sigma_{i-1} -
  x E_{R-2-i})), E_m = -1 + x + ... + x^m, AdmClamp = Bernstein-cellwise box
  clamp + parity fix toward zero: deterministic, search-free, closes exact-
  floor epochs R = 8..128 (6 s at 128) and, with additive slack, R = 256 at
  D0 = 1 (103 s) and R = 512 at D0 = 8 (2328 s), each witness re-verified
  exactly. Since D0(R) = o(R) per epoch suffices for C* <= 1 + gamma*
  (confirmed proof), these closures attain C = 1 + gamma* = log_5(5 phi^2)
  with additive slack <= 9: T(n) <= n + 1 + floor(gamma* n) + 8 for ALL
  critical values n <= 1023.
  (S) STRUCTURE (referee-confirmed exhaustively to R = 1024). The reduced
  coordinate has the witness-independent closed form C = q^{R-1} - x E_{R-2}
  solving qC = p^R + q^R - p(p-q); the ballot normal form Delta_i = (p-q) +
  c_i (i <= R-2), Delta_{R-1} = -1 reduces (*) to the MASTER EQUATION
  sum_{i<=R-2} x^i c_i = q^{R-1} - E_{R-1}, via the identity
  (p-q)(1+...+x^{R-2}) = E_{R-2} + 2x^{R-1}; and q^{R-1} - E_{R-1} is even
  coefficientwise IFF R = 2^t (dyadic parity theorem) — so at dyadic epochs
  parity vanishes entirely and the all-R problem is pure box-capacity
  transportation of the explicit target G_R = (q^{R-1} - E_{R-1})/2. The
  normal form characterizes rule/slim witnesses, NOT all witnesses (the fat
  doubling, direct-beam and lp R=128 witnesses violate it yet verify).
  (E) THE BALLOT CUT. delta_{i,0} = Delta_i(1) = +-1 is forced every row, so
  |sigma_i(1)| <= R-1-i with sigma_i(1) == i+1 mod 2 — a scalar cut invisible
  to per-cell clamping, STRUCTURALLY NECESSARY at R = 128 (a no-cut hunt
  reached row 125 with every state provably dead at sigma(1) = 16 vs budget
  2), slack-covered at R <= 64. Any all-R rounding proof must schedule the
  sign word as an explicit ballot path; the winning LP trajectory rides the
  cone boundary at exactly 1/row over the last 14 rows.
  (T) TRANSPORTATION FACTS. In f-coordinates (f = (binom - delta)/2) parity
  is free and the real relaxation at R = 128 is feasible (exact certificate
  committed); 21.24% of cells saturate including every forced top cell, so
  LP margin can never certify integer feasibility — a rounding theorem is
  the correct target. An exact S_L endgame decision procedure (L <= 2
  complete at e = 0, sandwich at e = 1; selftests 60/60, 30/30, full-scale
  positive controls 12/12 and 8/8) decides absorbability of endgame
  residuals.
  SCOPE: upper-bound attainments matching the balanced-block floor; the
  general-class floor stays OPEN (MISTAKE-361). OPEN: R = 1024 (rule ladder
  D0 = 16/32/64 running), the D0(R) growth law (data: 0 through R = 128,
  1 at 256, <= 8 at 512 — sublinear so far), and the all-R transient
  theorem, now reduced to: represent the explicit alternating-binomial
  target G_R by even boxed digits with a scheduled ballot path.
source: boxeph-2026-08-03-multifront-r128
depends_on:
  - THM-3302
  - THM-3002
  - THM-3026
  - THM-2966
related:
  - THM-3029
  - THM-3009
  - HYP-9061
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
script: 04-computation/amm12592_r128_independent_referee_boxeph.py
witnesses:
  - 04-computation/amm12592_floor_witness_R128_direct.json
  - 04-computation/amm12592_floor_witness_R128_lp.json
  - 04-computation/amm12592_floor_witness_R128_rule_boxeph.json
  - 04-computation/amm12592_witness_R256_ruleA_D0_1_boxeph.json
  - 04-computation/amm12592_witness_R512_ruleA_D0_8_boxeph.json
output: 05-knowledge/results/amm12592_r128_independent_referee_boxeph.out
---

# THM-3329 -- R=128 triple closure, the attractor rule, and attainment to n = 1023

## 0. Statement

With `gamma* = log(phi)/log(sqrt 5)` and the epoch representation problem
(*) of THM-3002:

1. `(*)` at the exact `gamma*` floor profile has solutions at `R = 128`
   (three independent witnesses, referee-verified), extending THM-3302's
   attainment to all `n <= 255` at `T(n) = n + 1 + floor(gamma* n)`.
2. The deterministic attractor rule closes `R = 256` at `D0 = 1` and
   `R = 512` at `D0 = 8`; assembling all closed epochs,
   `T(n) <= n + 1 + floor(gamma* n) + 8` for all `n <= 1023`, i.e. the
   golden deadline slope `C = log_5(5 phi^2)` is attained with additive
   slack at most `9` through `n = 1023`.
3. The all-R problem is equivalent (master equation) to representing the
   explicit target `G_R = (q^{R-1} - E_{R-1})/2` by even boxed digits, with
   parity vanishing identically at dyadic `R`.

Bodies of the proofs, the rule, the ballot cut, the S_L procedure, and the
transportation anatomy are in the committed scripts/notes listed in the
frontmatter and in `05-knowledge/results/amm12592-allR-family-hunt-boxeph.md`,
`amm12592-r128-endgame-algebra-boxeph.md`, and the referee output. All
computations exact (int/Fraction); no floats in any decision.

## 1. Why this matters for (Q)

The bracket for `C*` is now: upper — golden attained to `n = 1023` with
bounded slack and a search-free rule whose only obstruction is a
quantified transient; lower — balanced-block floor at golden
(audit-hardened), general class OPEN (MISTAKE-361). The remaining
obligations are therefore exactly:

- **all-R:** a transient estimate for the rule (or a rounding theorem for
  the `G_R` transportation problem with a scheduled ballot path), giving
  `C* <= log_5(5 phi^2)` unconditionally;
- **general floor:** a deadline-bounded routing window, giving
  `C* >= log_5(5 phi^2)`.

Together they would settle `C* = log_5(5 phi^2)` for AMM 12592's part-(c)
frontier question.
