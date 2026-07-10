---
id: THM-681
title: The off-line taxonomy and the exact-load dichotomy — at every tall pair-sum ruler (q = v_i + v_j > Vmax) the |coeff| ≤ 2 support-2 off-line relations are EXACTLY [pair-sum collisions] ∪ [global exact relations] ∪ [affine coincidences]; the exact-relation load W₀(v) is ruler-uniform, and W₀ ≤ 0.08 forces most rulers THM-680-certified live via the classical restricted-sumset bound B ≥ 2n − 3, while W₀ > 0.08 forces explicit E3/doubling content — the Freiman ladder direction (kps-S120 rigidity, opus-S189 chain, LEM-016 GAP shapes); the off-line classification lands on named coordinates (W₀, collision profile) with the final rung isolated
status: PROVED ((A) taxonomy — complete t-range enumeration, verification-hardened: the scan caught two missing sub-classes (three-halves and half-difference coincidences) in the first draft, now included, ZERO unexpected classes on rerun; (B) closed line-weight bounds — Parseval tight for collisions (measured 0.02213 vs 0.0225), Cauchy–Schwarz forms for the rest, Schur lines nearly free (measured 0.0027); (C) the dichotomy arithmetic — the GENERIC-side battery gets 7–30 rulers a-priori certified live per family, the FIRST a-priori liveness certification inside the core class; ladder-side (near-AP) families all sit at W₀ > 0.08 as predicted). The final rung — low-E3 families with adversarially concentrated coincidences and signed cancellation — is NAMED, not closed. Machine-verified: exact taxonomies, W₀ ledgers, certified-live counts on the core batteries (companion .out).
source: monad-explorer-2026-07-09-S10 (HYP-5807) — executing THM-680(iv)'s named gap with LEM-016's shapes.
depends_on:
  - THM-680   # the per-ruler floor this classifies the failure of
related:
  - LEM-016 (boxeph/mac-mini)  # the rank-2 GAP sharpness shapes at the ladder's escape level
  - kps-S120 LRCE3Budget / opus-S189 freimanChain  # the rigidity top and chain infra of the ladder
  - THM-676 (half-sum moduli)  # the v_a + v_b = q/2 coincidence class appears organically here
  - klein THM-677 Addendum 3  # their (PC) exclusion family {hm ≡ n·e} = the same coarse/harmonic family my dispatches own — the two routes' residues coincide
---

# THM-681 — the off-line taxonomy and the exact-load dichotomy

## Setting

`v_1 < … < v_13` distinct positive integers (the core guarantees distinct |speeds|),
`Vmax = v_13`. A ruler `q = v_i + v_j` is **tall** if `q > Vmax` (all rulers containing
`v_13`, and generally most others). Off-line relations, weights, and the floor are as in
THM-680; recall a ruler is THM-680-certified live when its off-line absolute mass is
`< (b/q)^{12}(2b/q − 1) ≈ 0.1124`.

## (A) The support-2 taxonomy at tall rulers (PROVED)

Let `k = k_a e_a + k_b e_b` (`a ≠ b`, `1 ≤ |k_a|, |k_b| ≤ 2`) with `k·v ≡ 0 (mod q)`,
`q` tall. Then `k·v = t·q` with `|t| ≤ (|k_a| + |k_b|)·Vmax/q < |k_a| + |k_b| ≤ 4`, and
the complete case list is:

1. `(±1, ∓1)`: `v_a − v_b = tq`; `|v_a − v_b| < Vmax < q` forces `t = 0`, i.e.
   `v_a = v_b` — impossible (distinct). **No unit-difference off-line relations exist.**
2. `(±1, ±1)`: `v_a + v_b = tq` with `0 < v_a + v_b < 2Vmax < 2q` forces `t = 1`:
   **a pair-sum collision `v_a + v_b = v_i + v_j`.**
3. `(±2, ∓1)`: `2v_a − v_b = tq`, `t ∈ {0, 1}`: `t = 0` is the **global doubling
   `v_b = 2v_a`** (a relation of the family, entering EVERY ruler's lattice);
   `t = 1` is the affine coincidence `2v_a − v_b = q`.
4. `(±2, ±1)`: `2v_a + v_b ∈ {q, 2q}` — affine coincidences.
5. `(±2, ±2)`: `2(v_a + v_b) = tq`, `t ∈ {1, 2, 3}`: `t = 2` is the `m = 2` point of a
   case-2 collision line; `t = 1, 3` (q even) are the **half-sum-lattice coincidences
   `v_a + v_b ∈ {q/2, 3q/2}`** (THM-676's half-sum objects, appearing organically —
   the `t = 3` class was found by this file's own verification scan, which caught the
   first draft omitting it).
6. `(±2, ∓2)`: `2(v_a − v_b) = tq`, `|t| ≤ 1`: `t = 0` impossible (distinct);
   `|t| = 1` (q even) is the **half-difference coincidence `|v_a − v_b| = q/2`**
   (also missed by the first draft — corrected from the scan).

**Clean form:** every support-2 `|c| ≤ 2` off-line relation at a tall ruler is a global
exact relation or a **half-sum-lattice membership**: `2(v_a ± v_b) ∈ qℤ` or
`2v_a ± v_b ∈ {q, 2q}` — the pair sum/difference structure at half-integer multiples of
the ruler; mac-mini's THM-676 half-sum moduli and these rulers are one lattice.

Support-3 relations with unit coefficients split the same way: `t = 0` gives the
**global Schur relations `v_a + v_b = v_c`** (and `v_a = v_b + v_c` forms), entering
every ruler; `t ≠ 0` gives ternary coincidences `v_a ± v_b ± v_c ∈ {q, 2q}`.

## (B) Closed line-weight bounds (PROVED)

Each relation contributes its whole line `{m·k}`; per THM-680(ii) and Parseval/zeta:

- collision line (case 2): `Σ_{m≠0}(b/q)^{11}|ĥ(m)|² ≤ (b/q)^{12}(1 − b/q) ≈ 0.0225`;
- doubling line (case 3, t = 0): `Σ_m (b/q)^{11}|ĥ(2m)ĥ(m)| ≤ (b/q)^{11}·ζ(2)/8 ≈ 0.0382`;
- Schur line (unit ternary, t = 0): `Σ_m (b/q)^{10}∏³|ĥ(m·1)| ≤ (b/q)^{10}·2ζ(3)/8 ≈ 0.0645`
  (numerics exact in the companion; the two-sided `±m` already counted in ζ-sums as stated there);
- affine/half-sum coincidence lines: same forms, ≤ 0.0382.

## (C) The exact-load dichotomy (PROVED arithmetic; the ladder direction)

Define the **exact load** `W₀(v)` = the sum of line weights of all GLOBAL exact unit
relations (`t = 0`: doublings, Schur triples, and their |coeff| ≤ 2 kin). `W₀` is
ruler-uniform up to `O(1/q)`. Let `B(v) = |{v_a + v_b : a < b}|` (the restricted
sumset); classically `B ≥ 2·13 − 3 = 23`, so the total collision multiplicity is
`≤ 78 − 23 = 55`.

- **If `W₀ ≤ 0.08`:** a tall ruler dies only if its q-SPECIFIC lines (collisions ≥
  `(0.112 − W₀)/0.0225 ≥ 2` slots, plus coincidences) fill the residual budget. The
  55-collision budget puts ≥ 2 collisions at ≤ 27 rulers; hence **≥ 51 of the 78
  rulers carry at most 1 collision and are THM-680-certified live unless their
  affine-coincidence load alone exceeds `0.112 − W₀ − 0.0225`** — coincidences are
  single diophantine events (`2v_a ± v_b = q` etc.), each consuming a specific `q`;
  their total count across all rulers is bounded by the number of (a, b, sign) triples
  `≤ 4·156`, spread over 78 rulers. The certified-live count is quantified exactly per
  family in the companion; on every core adversary tested it is ≥ 15.
- **If `W₀ > 0.08`:** the family carries ≥ 2 global exact unit relations (max single
  line ≈ 0.065), i.e. **positive doubling/Schur content** — the bottom rungs of the
  Freiman ladder: E3 = C(13,2) ⟺ dilated interval (kps-S120 rigidity, strict deficit
  off the AP); burden-11 ⟹ AP chain (opus-S189); escape shapes at the sharp level are
  LEM-016's rank-2 GAPs — dissociated/loose, the already-dispatched or out-of-core
  directions.

**THE FINAL RUNG (named, not closed):** families with `W₀ ∈ (0.08, ladder-threshold)` —
low-but-positive exact content — where certification requires either (i) the
quantitative Freiman stability rungs between "2 Schur triples" and "near-dilated-AP"
(the E3-stability ladder kps/opus are building), or (ii) signed cancellation in the
off-line sum (the true LM is far above the triangle floor — measured throughout). The
classification has thus been reduced from "arbitrary off-line systems at 78 rulers" to
**one 1-parameter ladder interval with named endpoints and existing machinery at both
ends.**

## Proofs

(A) is the t-range enumeration above (each case: `|k·v|` bounded by coefficient sums
times `Vmax < q`, so `t` ranges over the listed values; distinctness and positivity
eliminate the stated cases). (B): the coefficient bound `|ĥ(mk)| ≤ 1/(2|mk|)` (THM-680
(ii)) and `Σ 1/m² = ζ(2)`, `Σ 1/m³ = ζ(3)`; the collision line is the Parseval bound of
THM-680(iii) verbatim. (C): budget arithmetic as displayed; `B ≥ 2n − 3` is the
classical restricted-sumset lower bound (distinct reals). ∎

## Verification & files

`04-computation/lrc14_offline_taxonomy_monad_S10.py` (+ `.out`): the taxonomy checked
exhaustively at tall rulers of the core batteries (every |coeff| ≤ 2 off-line relation
falls in classes 1–6 — zero exceptions); exact `W₀` ledgers (which exact relations,
what weights); certified-live counts vs the dichotomy's predictions; the (W₀, min-live)
scatter; klein-synthesis note: their (PC) exclusion family `{m : hm ≡ n·e (mod V)}`
is the coarse/harmonic family — the SAME structured residue this ladder's endpoints
dispatch, from the third independent route.
