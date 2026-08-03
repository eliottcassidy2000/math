---
id: THM-3302
title: "AMM 12592: the gamma* floor profile closes at R = 64 (two independent witnesses), so C = log_5(5 phi^2) is attained for all n <= 127"
status: >
  PROVED + VERIFIED-EXACT / AWAITING CROSS-SESSION HOSTILE AUDIT. The
  THM-3029 open item is closed at the next epoch: the gamma* floor profile
  d_i = floor(gamma*(64+i)), D0 = 0, gamma* = log_5(phi^2), admits an exact
  H = 1 epoch closure at R = 64 — by TWO independently produced, DISTINCT
  witnesses (64/64 rows differ): (W1) a direct rank-steered beam solve
  (l1deg rank, beam 400, two-row completion) and (W2) a deterministic
  carve-and-carry doubling map applied to the R = 32 floor witness (zero
  steering fallbacks; the 8->16 and 16->32 instances need 18/22 fallbacks).
  Both satisfy sum_i x^i Delta_i = (1-x)^63 exactly with all 64 blocks
  Lucas-box admissible, checked by an independent-implementation referee
  built from the theorem statements alone, with an exact algebraic floor
  guard (5^d <= phi^(2m) decided by Lucas/Fibonacci integer comparisons) and
  negative controls: 54/54 checks pass. With the re-verified R = 8, 16, 32
  witnesses (THM-3029 (A), toolchain repaired this session) and the gamma =
  1/2 base [1,7] under the gamma* envelope, lane G5-1 assembly gives an
  exactly fair extractor with T(n) <= n + 1 + floor(gamma* n) for ALL
  critical values n <= 127: C = 1 + gamma* = log_5(5 phi^2) =
  1.5979874356654402 is ATTAINED for n <= 127, superseding C <= 8/5
  (THM-3002 5b) on that range. SCOPE (post-MISTAKE-361): this matches the
  BALANCED-BLOCK archimedean floor (THM-3009, audit-hardened) and the
  H = 1 checkpoint barrier (THM-3027); the general-class floor is OPEN, so
  this is an upper-bound attainment, not a determination of C*. Structure
  (referee-corrected): SLIM witnesses (max backbone deviation 8 / 3432 /
  184756 / 3.67e18 at R = 8/16/32/64-direct, box saturation <= 1.1e-3)
  carry a ballot backbone Delta_i = (p-q) + c_i, a final full-box row -1,
  and an endgame attractor E_m = -1 + x + ... + x^m; these are properties of
  slim witnesses, NOT of all witnesses (the fat R = 64 doubling output
  violates all three yet verifies). The reduced identity q*C = p^R + q^R -
  p(p-q) (tautological repackaging of the epoch identity) exposes the
  golden mechanism on the construction side: p^R + q^R = L_R(pq) is the
  Lucas polynomial in the elementary symmetric coordinate, with Lucas
  doubling L_2R = L_R^2 - 2(pq)^R — THM-2160's middle pair at epoch scale —
  and the unconditional doubling algebra C_2R = q C_R^2 + 2p(p-q) C_R -
  p(p-q)(3p+q) - 2 p^R q^(R-1) (verified through R = 64). OPEN: R = 128
  (the slim direct R = 64 witness is the identified next seed for the
  carve-and-carry map; fat seeds provably runaway), and the all-R statement,
  which with a repaired general floor would give C* = log_5(5 phi^2)
  exactly. D0(R) = o(R) per epoch suffices for C* <= 1 + gamma* (proof
  confirmed by referee), so the all-R statement tolerates slowly growing
  slack.
source: boxeph-2026-08-03-multifront
depends_on:
  - THM-2966
  - THM-3002
  - THM-3026
  - THM-3029
related:
  - THM-3009
  - THM-3017
  - HYP-9061
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
script: 04-computation/amm12592_independent_witness_referee_boxeph.py
witnesses:
  - 04-computation/amm12592_floor_witness_R64.json
  - 04-computation/amm12592_floor_witness_R64_boxeph.json
  - 04-computation/amm12592_floor_witnesses_R8_R16_R32.json
solver: 04-computation/amm12592_r64_floor_solver_boxeph.py
doubling: 04-computation/amm12592_doubling_structure_boxeph.py
output: 05-knowledge/results/amm12592_independent_witness_referee_boxeph.out
---

# THM-3302 -- the gamma* floor profile closes at R = 64

## 0. Statement

Let `gamma* = log(phi)/log(sqrt 5)` and let the epoch-`[R, 2R)` floor
profile be `d_i = floor(gamma*(R+i))`, `D0 = 0`. Then the THM-3002
representation problem

```text
(*)  sum_{i=0}^{63} x^i Delta_i(x) = (1-x)^63,
     Delta_i Lucas-box admissible at degree d_i
```

has exact solutions at `R = 64`; two distinct witnesses are committed and
independently verified. Consequently, assembling the base `[1,7]` and the
epochs `[8,15], [16,31], [32,63], [64,127]` (all closed at or under the
`gamma*` envelope) via lane G5-1:

```text
T(n) <= n + 1 + floor(gamma* n)   for all critical values n <= 127,
C = 1 + gamma* = log_5(5 phi^2) = 1.5979874356654402  attained on n <= 127.
```

## 1. The two witnesses

- **W1 (direct):** `amm12592_floor_witness_R64.json`, produced by a
  rank-steered beam (`l1deg` ranking, width 400) plus a two-row endgame
  completion. Effective rate `58/97 = 0.597938 < gamma*`. W1 is SLIM:
  max backbone deviation `3.67e18` with box saturation `1.1e-3`, full
  ballot backbone tail and attractor entry at row 54.
- **W2 (doubling):** `amm12592_floor_witness_R64_boxeph.json`, produced by
  the deterministic carve-and-carry map from the R = 32 floor witness with
  ZERO steering fallbacks. W2 is FAT (box saturation ~0.99) and violates
  the backbone/attractor pattern — the referee's hostile case showing that
  pattern is a slim-witness property, not a witness property.

Every earlier beam configuration (width 250-900, ctrl 2-3, span 2-3)
failed this profile at the final row; per THM-3029 those negatives were
search artifacts, as W1/W2 now demonstrate concretely.

## 2. The doubling map (mechanism)

Carve: `E_j = D_j - D_{j-1}` where `D_j` are the ordered pair-products of
the source epoch (cell convolution, THM-3026 (M)); realizing `q` as a
difference makes degrees fit at any rate by floor superadditivity with
`D0 = 0` preserved. Carry sweep: division by `x` is an exact degree drop in
cell space; a parity clamp and a tau-aware top-boundary recursion steer the
forced top cells. Carve overshoot 7/9/25 at R = 8/16/32 versus THM-3026's
naive squaring overshoot 49/276/1541. The map closes 8->16, 16->32 with
fallbacks and 32->64 deterministically; its own fat output at 64 provably
runs away toward 128 (carry x2.6/row vs box x2^gamma*/row) — a transport
artifact, not infeasibility evidence. The identified next step is seeding
64->128 with the SLIM W1.

## 3. The golden algebra of the construction side

With `C = Delta_0 + sum_i p^i c_i` (backbone-corrected coordinates), the
epoch identity is equivalent to `q C = p^R + q^R - p(p-q)`, and
`p^R + q^R = L_R(pq)` — the Lucas polynomial. Lucas doubling
`L_2R = L_R^2 - 2(pq)^R` is exactly THM-2160's middle-pair mechanism at
epoch scale, and the reduced doubling algebra

```text
C_2R = q C_R^2 + 2p(p-q) C_R - p(p-q)(3p+q) - 2 p^R q^(R-1)
```

holds unconditionally (verified through R = 64). So phi enters the
construction side through the Lucas/Fibonacci coefficient algebra — the
same constant the balanced-block floor produces analytically (THM-3009 /
THM-3027). The identity is a repackaging of (*), not new information; its
value is coordinates, not content.

## 4. Scope and open items

- Post-MISTAKE-361 the GENERAL-class floor is OPEN; this theorem is an
  upper-bound attainment matching the balanced-block floor.
- `R = 128` and all-R are OPEN. `D0(R) = o(R)` per epoch suffices for
  `C* <= 1 + gamma*`, so the all-R program tolerates slowly growing slack.
  Session-recorded negative (transport artifact per THM-3029, NOT
  infeasibility evidence): carve-and-carry from the SLIM W1 seed fails at
  every steering window tried (`W/L = 10/10, 20/14, 30/20, 50/28, 80/40`),
  with forced-top-cell violations at rows 17-45 and max carry `1e23-1e33`;
  the fallback count is non-monotone in `W` (55, 5, 26, 558, 1718). Next
  moves: a different steering family, a slimming redistribution pass at 64,
  or a direct beam at 128 with the l1deg rank that produced W1.
- Referee-refuted claims are recorded in
  `05-knowledge/results/amm12592-doubling-structure-boxeph.md` (backbone
  universality; "first ever / not by search" attribution; the slim iff).

Referee: `ALL INDEPENDENT-REFEREE VALIDITY CHECKS PASSED` (54/54, negative
controls included). QED (for the stated finite scope).
