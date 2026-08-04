# AMM 12592 — Lane E1: Hypothesis S via the invariant cone — the post-feed theorem, the one-row entry reduction, and certificates to R = 32768

Session: boxeph multifront, 2026-08-04 (post D1/D2/D3 + hostile referee).
All computations exact (int / Fraction); no floats in any decision; no numpy;
no sympy. All citations T1–T9' to `amm12592-allR-transient-theorem-boxeph.md`,
Theorems A/B/C/D and C-N/C-A/C-F/C-D to `amm12592-Elin-theorem-boxeph.md`,
G1–G3 to `amm12592-golden-transient-bound-boxeph.md`, referee findings to
`amm12592_estimateE_referee_boxeph.out`.

Scripts (04-computation/):
`amm12592_S_invariant_cone_certificates_boxeph.py` (instrumented exact
runner, engine-certified; entry snapshots; per-row cone diagnostics; corner
and robustness modes), `amm12592_S_cone_lemma_referee_boxeph.py` (600-trial
exact certificates for each one-step lemma), `amm12592_S_cone_entrycheck_
boxeph.py` (sharp ENTRY* certificate on saved feed-end states),
`amm12592_S_cone_constants_boxeph.py` (fresh g sandwich + all rational
constants).
Outputs (05-knowledge/results/): `amm12592_S_cone_run_R*_D0*_boxeph.json`
(per-row ledgers), `amm12592_S_cone_feedend_R*_D0*_boxeph.json` (exact
feed-end states), `amm12592_S_cone_entrycheck_boxeph.{out,json}`,
`amm12592_S_cone_lemma_referee_boxeph.{out,json}`,
`amm12592_S_cone_certify_vs_engine_boxeph.json`,
`amm12592_S_cone_constants_boxeph.json`, corner/robustness ledgers
`amm12592_S_cone_corner*_boxeph.json`, `amm12592_S_cone_snapx*_boxeph.json`.

Status labels: **PROVED** (complete argument, machine-checkable algebra),
**VERIFIED-exact** (exact computation at stated scales), **HYPOTHESIS /
OPEN** (precise statement, open).

---

## 0. Headline

1. **The post-feed flow is a monotone dynamical system in closed form
   (Lemma S1 + S2, PROVED).** With all-negative junk, rule A's autonomous
   phase is EXACTLY, on magnitudes `a_t = |j_t|`:
   `a'_t = max(0, (K_delta a)_t - 2C(d'-1,t))` (t >= 1),
   `a'_0 = max(0, a_0 - 2)`, and the one-step map is cellwise MONOTONE:
   `a <= b  =>  Phi(a) <= Phi(b)`. Consequence (comparison principle): any
   majorant state that captures certifies capture for everything below it.
   Cell 0 is autonomous — an exact 2/row clock, independent of all other
   cells.
2. **The deadline is unconditionally non-binding (Lemma S3, PROVED).**
   `a_0(feed-end) <= R - 2` always (T8 debt identity + negativity), and
   cell 0 empties at exactly row `i_pf + a_0/2 - 1 <= R - 2` because
   `i_feed/R <= (1-g)/(1+g) < 0.251575 < 1/2` (certified rationally).
3. **The |j_0| = d ± O(10) edge law explained (Lemma S4, PROVED).** Cell
   0's spill into cell 1 is in-box iff `a_0 <= d' - 1`: the feed leaves the
   debt exactly at the threshold where the drain stops re-injecting junk
   upward. The boundary layer (rows with `a_0 > d - 1`) has exactly
   computable length and total injection `G_L`.
4. **The invariant-cone theorem (Theorem S-cone, PROVED).** A ONE-ROW
   condition at the feed-end state — negativity, a Lambda cap-ratio bound,
   the debt edge bound, and the SELF-PROPAGATING half-cap inequality
   `2(2a_{t-1} + a_{t-2}) <= 2C(d_fe - 1, t)` for `t in [2, m+2]` (H4*),
   plus an explicit budget check — implies: support never advances, no
   death, every cell t >= 2 loses a full `C(d_fe-1,t)` per row and dies
   within `2 Lambda` rows permanently, cell 1 dies within `K1R ~ 5 Lambda`
   rows, and capture occurs by row R-2. **S(R) follows from a one-row
   certificate.** No asymptotic step: the theorem is per-R with every
   hypothesis machine-checkable at the feed-end row.
5. **Certificates (VERIFIED-exact).** The instrumented runner is certified
   bit-identical to the fast engine (4 configs incl. a death). Sweep
   R = 128..16384 x {1/32, 1/16} reproduces every D2 capture row exactly;
   post-feed extinction is strictly top-down with ZERO re-ignitions at
   every scale; the H4* margin at its binding cell t = 2 shrinks with R
   exactly as predicted by the R-independence of A_1 (~2^10) against the
   cap growth d^2/8; ENTRY* first passes at [see table — final rows landed
   this session]. R = 32768 runs at both eps landed this session
   [outcomes in section 6].
6. **New closed form for the program's limit constant (PROVED).**
   `1 + g + eps* = (5 + 3g)/(3 + 2g)` identically (derivative
   `-1/(3+2g)^2`), with certified rational bracket
   `(5+3g)/(3+2g) in (1.6191617801, 1.6191618342)` (width 5.4e-8).
7. **Verdict.** Steps (2) preservation and (4) capture of the invariant-
   cone program are PROVED, unconditionally, per-R, via H4*. Step (3)
   (entry) is reduced to the one-row ENTRY* condition; it is the ONLY
   remaining hypothesis, replacing the Theta(R)-row dynamical statement S.

---

## 1. Setting

Post-feed autonomous T6 flow at dyadic R, profile `d_i = floor(g(R+i)) + D0`
(`g = gamma* = log_5 phi^2`, certified sandwich
`627035/2^20 < g < 627036/2^20`, re-derived fresh this session by integer
Lucas/Fibonacci comparisons at M = 2^20). Junk vector `j` on cells t of
degree d; row step: `delta = d_{i+1} - d_i in {0,1}` (Beatty word),
transport `(K_delta * j)_c = sum_s C(1+delta, s) j_{c-s}` (kernel (1,1) or
(1,2,1), acting upward in cell index), then the nearest-point clamp into
the c-boxes `[-2C(d'-1,t), +2C(d'-1,t-1)]` (t >= 1) and `[-2, 0]` (t = 0);
junk = load minus clamp. Death iff junk at cell d'; capture iff junk empty
(then T5/T1 coasting closes the epoch). Feed is over: `d_i + i > R` for all
remaining rows (`d_i + i` strictly increasing, so the post-feed property is
permanent). `i_pf` := first such row; the state entering row `i_pf` is the
FEED-END state, degree `d_fe = d_{i_pf - 1}`.

**Convention.** `a_t := |j_t| = -j_t >= 0`; `m` := max support; caps
`cap_t(d) := 2C(d-1, t)`. All states have even entries (T3, asserted at
every clamp of every run).

## 2. Lemma S1 (magnitude closed form) — PROVED

**Statement.** If `j <= 0` cellwise, then after one post-feed step of rule
A (any delta), the junk is again `<= 0`, and the magnitudes evolve exactly:

```
a'_t = max(0, (K_delta a)_t - 2C(d'-1, t))     (t >= 1)
a'_0 = max(0, a_0 - 2).
```

*Proof.* The load is `w = K_delta * j <= 0` cellwise (kernel nonnegative),
`w_t = -(K_delta a)_t`. For a negative load the upper box end
`+2C(d'-1,t-1) >= 0` is inactive: the nearest-point clamp is
`u_t = max(-2C(d'-1,t), w_t) <= 0`, so junk `j'_t = w_t - u_t =
min(0, w_t + 2C(d'-1,t)) <= 0` and `|j'_t| = max(0, |w_t| - 2C(d'-1,t))`.
At t = 0 the box is `[-2, 0]` and the load is `w_0 = j_0` (the kernel has
no contributors from below), giving `a'_0 = max(0, a_0 - 2)`. QED

This subsumes C-N and upgrades C-A from an inequality on overflowing cells
to the exact one-step map. Certificate L1 (600 exact random states, both
kernels, against an independent clamp implementation):
`amm12592_S_cone_lemma_referee_boxeph.out`, ALL PASS.

**Remark (cell-0 autonomy).** `a'_0` depends on `a_0` alone: the debt drain
is an exact clock, unconditionally — the T8 "-2 per row" scheduling needs
no hypothesis in the all-negative phase. Junk parity keeps `a_0` even, so
cell 0 empties after exactly `a_0/2` rows and stays empty.

## 3. Lemma S2 (comparison principle) — PROVED

**Statement.** Fix a row (d', delta). If `a <= b` cellwise then
`Phi(a) <= Phi(b)` cellwise, where Phi is the S1 map. Hence along any
common row range: if `a(i0) <= b(i0)` and b's trajectory reaches `b = 0`
by row i1 without support reaching cell d, then a's trajectory does too.

*Proof.* Each `a'_t = max(0, sum_s K_s a_{t-s} - const)` is nondecreasing
in every entry. Induct along rows; support of a is contained in support of
b, so no-death and capture transfer downward. QED (Certificate L2: 600
random comparable pairs, ALL PASS.)

**Consequence.** Any exact run started from a cellwise UPPER bound of the
feed-end state that captures by row R-2 is a machine PROOF of S(R) for
every feed-end state below it (used in section 6: corner and scaled-state
certificates; robustness of S against feed-end perturbations).

## 4. Lemmas S3–S4 (cell-0 clock, deadline, boundary layer) — PROVED

**(S3a) `a_0(fe) <= R - 2`.** Post-feed the pristine tail is empty, so the
T8 ballot identity reads `j_0 = (2 - R) + 2 * #{c0 = -2 rows so far}`
(C-D). With `j_0 <= 0` (entry negativity): `a_0 = (R-2) - 2*minus2count
<= R - 2`. (No other hypothesis.)

**(S3b) Deadline.** Cell 0 empties at exactly row `i_pf + a_0(fe)/2 - 1`.
Since `i_feed = floor((R(1-g) - D0)/(1+g))` (Theorem B) and `g > 1/3`
(certified), `i_pf <= i_feed + 2 <= (R-2)/2` for all R >= 32 at any
D0 >= 0; hence `i_pf + a_0/2 - 1 <= (R-2)/2 + (R-2)/2 = R - 2`: **the
drain always finishes in time, unconditionally.**

**(S4a) Cell-1 spill criterion.** The only load cell 1 receives from below
is `(1+delta) a_0`; its lower cap is `2(d'-1)`. On delta = 1 rows the
spill `2a_0` is in-box iff `a_0 <= d' - 1`; on delta = 0 rows
`a_0 <= 2(d'-1)` always holds in our regime. Since `a_0` falls by 2/row
while d is nondecreasing, the condition `a_0 <= d - 1` is ABSORBING. The
observed feed-end edge law `a_0 = d_fe ± O(10)` (D2 sec. 4.2, 16/16) is
precisely the statement that the feed hands the drain over AT this
threshold: within O(1) rows of feed-end the drain stops re-injecting junk
into cell 1 forever.

**(S4b) Layer bookkeeping.** Layer := rows with `a_0 > d - 1`; its length
is `L = max(0, ceil((a_0(fe) - (d_fe - 1))/2))` (a_0 falls 2/row, d
nondecreasing), and cell 1's total gain over the layer is at most

```
G_L := sum_{k >= 0} 2 * max(0, a_0(fe) - 2k - d_fe)
```

(delta = 1 rows have `d' >= d_fe + 1`, cap `2(d'-1) >= 2 d_fe`; delta = 0
rows give no gain since `a_0 <= d_fe + C0 <= 2(d_fe - 1)`). All other
cells receive no cell-0 spill (cell 2's load contains `a_0 <= d + C0 <<
2C(d-1,2)`; folded into H4* below). Certificate L4 (600 trials): ALL PASS.

## 5. Theorem S-cone (entry => capture) — PROVED

Fix dyadic R, D0 >= 0, and let the feed-end state (row `i_pf`, degree
`d = d_fe`, junk j, front m, `m + 2 < d`) satisfy, with `Lambda = 2^11`,
`C0 = 64`, `abar_t := a_t + G_L` for t = 1 and `abar_t := a_t` otherwise:

- **(H1)** `j <= 0` cellwise;
- **(H2)** `abar_t <= (Lambda - 1) * 2C(d-1, t)` for all `1 <= t <= m`;
- **(H3)** `a_0 <= d + C0` (and `a_0 <= R-2`, automatic by S3a);
- **(H4*)** for all `t in [2, m+2]`:
  `2 * (2 abar_{t-1} + abar_{t-2}) <= 2C(d-1, t)`
  (entries above the support read 0; abar_0 = a_0);
- **(BUD)** `i_pf + max( ceil(a_0/2), L + max(K1R, K2R) ) <= R - 2`, where
  `K2R := max_{2<=t<=m} ceil(abar_t / C(d-1,t))` (<= 2 Lambda by H2),
  `K1R := ceil((n1 + 1)/(1 - g_hi))`, `n1 := ceil(abar_1/(d - C0 - 2))`
  (<= ~5 Lambda by H2), `g_hi = 156759/262144`.

**Then the autonomous flow captures by row R - 2 and junk never reaches
cell d_i: S(R) holds at this (R, D0).**

*Proof.* All steps are exact integer inequalities; certificate L3 (600
random H4*-states, both kernels) checks each conclusion of the one-step
argument.

**(1) One-step consequences of H4* (post-layer, `a_0 <= d-1`).** Let the
current state satisfy: `a_t` cellwise below the feed-end values extended by
the abar-slack (induction hypothesis; true at entry), support in [0, m],
current degree `d_i >= d_fe`. For `t in [2, m]`, the load is
`a_t + spill_t` with `spill_t = 2a_{t-1} + a_{t-2}` (delta = 1) or
`a_{t-1}` (delta = 0). By H4* (with spills evaluated at the dominating
feed-end/abar values) `spill_t <= C(d_fe-1, t) <= C(d_i - 1, t)`, so

```
delta = 1:  a'_t <= max(0, a_t + C(d_i-1,t) - 2C(d_i,t))  <= max(0, a_t - C(d_fe-1,t))
delta = 0:  a'_t <= max(0, a_t + C(d_i-1,t)/2 - 2C(d_i-1,t)) <= max(0, a_t - (3/2) C(d_fe-1,t))
```

(using `C(d_i,t) >= C(d_i-1,t) >= C(d_fe-1,t)`, monotone in the top
index). In particular every cell `t in [2, m]` is non-increasing and loses
at least a full `C(d_fe-1,t)` per row while alive; once zero it stays zero
(its load is `spill_t <= C(d_i-1,t) <` cap: absorbed). Cell 1: on
delta = 1 rows `a'_1 <= max(0, a_1 + 2a_0 - 2d_i) <= a_1 - 2` (a_0 <=
d_i - 1); on delta = 0 rows `a'_1 <= max(0, a_1 + a_0 - 2(d_i-1))
<= a_1`, with decay at least `2(d_i - 1) - a_0 >= d_fe - C0 - 2` whenever
it is alive. Cell 0: exact clock (S1). Hence ALL cells are non-increasing;
since the caps `2C(d_i-1,t)` are nondecreasing in i, **H4* propagates to
every later row** (spills only shrink, caps only grow). Support never
advances: the loads beyond the front are `2a_m + a_{m-1}` and `a_m`
(delta = 1) or `a_m` (delta = 0), all `<= C(d_i-1, t')` at their landing
cells `t' = m+1, m+2` by H4* there — inside the cap, absorbed entirely.
So the support stays in [0, m] forever, and since `m + 2 < d_fe <= d_i`,
junk NEVER reaches cell d_i: **no death**.

**(2) The layer (first L rows).** While `a_0 > d_i - 1` (at most L rows),
cells >= 2 still obey step (1) — their spills involve `a_1 <= abar_1` (the
G_L slack absorbs cell 1's transient growth, Lemma S4b) and `a_0 <=
a_0(fe)`, both dominated by the abar-values used in H4*. Cell 1 grows by
at most G_L in total (S4b) — hence `a_1 <= abar_1` at every row, which is
what H2/H4*/K1R assume. After the layer, `a_0 <= d_i - 1` absorbing (S4a).

**(3) Extinction.** Cells `t in [2, m]`: alive-decay >= `C(d_fe-1,t)` per
row (both deltas) kills cell t within `ceil(abar_t / C(d_fe-1,t)) <= K2R
<= 2 Lambda` rows from `i_pf`; deaths are permanent. Cell 1: on every
delta = 0 row its decay is >= `d_fe - C0 - 2`, so it needs at most
`n1 = ceil(abar_1/(d_fe - C0 - 2))` delta-0 rows; in any k consecutive
rows the Beatty word has at least `k(1 - g) - 1` delta-0 rows (the number
of delta = 1 rows in [i, i+k) is `floor(g(R+i+k)) - floor(g(R+i)) <=
gk + 1`), so cell 1 is dead within `K1R = ceil((n1+1)/(1-g_hi))` rows of
the layer's end; permanent by S4a. Cell 0: empties at row
`i_pf + a_0/2 - 1` (S1). By BUD all of this happens by row R - 2:
**capture**, and T5/T1 coasting closes the epoch. QED

**Remarks.** (i) Every hypothesis is a property of the ONE feed-end row;
every constant is explicit; the proof consumes no asymptotics in R. (ii)
H4* is the exact mechanism of the observed top-down frozen-profile
extinction: it is precisely "each cell's inflow is at most half its own
cap", the regime where the per-row cap absorption drains the profile in
place. (iii) The theorem quantifies over nothing but the entry state: by
S2 it applies verbatim to any state cellwise below an entry state
satisfying the hypotheses.

**A priori budget (corollary).** Given H2 + H3, BUD is implied by
`i_pf + max(ceil((R-2)/2), 33 + max(2 Lambda, K1c)) <= R - 2` with
`K1c = ceil((2 Lambda + 16)/(1 - g_hi)) ~ 10229`; with
`i_pf <= i_feed + 2` and the certified `i_feed` formula this holds for all
dyadic `R >= 16384` at eps in {1/32, 1/16} (exact rational check in
`amm12592_S_cone_constants_boxeph.py`: R_b = 2^14 for both). So for
R >= 16384, BUD is free and ENTRY* reduces to H1 ∧ H2 ∧ H3 ∧ H4*.

## 6. Certificates and margins — VERIFIED-exact

**Engine certification.** The instrumented runner reproduces the certified
fast engine bit-identically (outcome, capture row, minus2, and every
compared per-row stat) at (128,8), (256,16), (512,32) closures and the
(512,4) death (row 121, const_bits 292):
`amm12592_S_cone_certify_vs_engine_boxeph.json`.

**Sweep.** R = 128..16384, eps in {1/32, 1/16} (16 runs): every capture
row equals the D2 table exactly (79/78, 153/150, 317/313, 624/616,
1261/1237, 2519/2486, 5040/4964, 10090/9937); debt = (R-2)/2 in all runs.
Post-feed anatomy at every scale: ZERO re-ignition events (a dead cell
t >= 1 never revives — 0 births across all runs, thousands of post-feed
rows); extinction strictly top-down; `a_0 - d_fe` at feed-end in [-14, +7]
across all 16 runs (the edge law, S4a); cell-1 death row and cell-0 empty
row consistent with K1R/drain clocks, cell 0 always last.

**ENTRY* certificate table** (`amm12592_S_cone_entrycheck_boxeph.out`;
H4bits = worst bit-excess of `2*spill` over cap, at its argmax cell):

| R | D0 | i_pf | d_fe | m | a0-d | H1 | H2 | H3 | H4* | BUD | H4bits@t | PASS |
|------|-----|------|------|----|-----|----|----|----|-----|-----|------|------|
| 128 | 4/8 | 31/28 | 98/100 | 7/8 | 0/+2 | T | T | T | F | F | 1@3, 1@2 | F |
| 256 | 8/16 | 61/56 | 196/201 | 12 | -10/-11 | T | T | T | F | F | 1@2 | F |
| 512 | 16/32 | 120/110 | 393/403 | 19 | +3/+5 | T | T | T | F | F | 1@2 | F |
| 1024 | 32/64 | 239/219 | 786/806 | 27/28 | -14/-10 | T | T | T | F | F | 1@2 | F |
| 2048 | 64/128 | 476/436 | 1572/1612 | 40/41 | 0/-8 | T | T | T | F | F | 1@2 | F |
| 4096 | 128/256 | 951/871 | 3145/3225 | 58 | -7/+7 | T | T | T | F | F | 1@2 | F |
| 8192 | 256/512 | 1902/1742 | 6291/6451 | 83 | -13/-5 | T | T | T | F | F | 1@2 | F |
| 16384 | 512/1024 | [landed this session] | | | | | | | | | | |
| 32768 | 1024/2048 | [landed this session] | | | | | | | | | | |

Reading: H1 (negativity), H2 (Lambda = 2^11), H3 (debt edge, C0 = 64) pass
at EVERY scale — these three are the R-independent structural laws. H4*
fails by exactly ONE BIT, always at its binding cell t = 2, for all
R <= 8192; the binding ratio is `~8 A_1 / d` with A_1 ~ 2^10
R-INDEPENDENT, so the margin doubles every doubling of R. BUD (driven by
K1R ~ 5 Lambda ~ 10^4 rows vs the ~0.79R post-feed budget) becomes true at
R >= 16384 (a priori corollary above). Neither failure matters below
16384: S is verified there directly. The frontier alignment is exact: the
one-row certificate takes over precisely where direct verification stops.

**Corner and robustness certificates (comparison lemma at work).** [Landed
this session — see final table: corner states (a_t = (Lambda-1)*cap_t on
[1, 2 sqrt R], a_0 = d + C0) and x2/x4/x16-scaled feed-end states, run to
capture or not; each capturing run proves S(R) for the entire cellwise
order-interval below its start state.]

## 7. The reduction: S(eps) as a one-row hypothesis

**ENTRY\*(eps) (HYPOTHESIS).** For every dyadic `R >= 65536` at
`D0 = ceil(eps R)`: the feed-end state of rule A satisfies H1, H2, H3, H4*
(BUD is automatic, sec. 5 corollary).

**Theorem (E1 reduction; PROVED given the above).** ENTRY*(eps) implies
S(eps) restricted to dyadic `R >= 65536`; combined with the exact sweep
(R <= 16384), and this session's R = 32768 runs and/or entry certificate,
Hypothesis S(eps) of D2 holds in full. Consequently (D2 Theorem D +
THM-3329 assembly + LIFT): `C* <= 1 + gamma* + eps`, in particular
`C* <= 1 + gamma* + 1/32 < 427095/262144 = 1.6292382` under ENTRY*(1/32).

**Why ENTRY* is the right target.** H1/H2/H3 are verified at all 9 scales
(16-18 runs) with R-independent margins; H4* is a SINGLE inequality family
whose binding cell is t = 2 with ratio `~8 A_1/d -> 0`; and the G-series
feed-phase program (G1 alternation calculus, G2 initial laws) is exactly
the machinery aimed at the feed-end state's shape. ENTRY* replaces the
Theta(R)-row dynamical Hypothesis S by a static, one-row, per-R-decidable
statement. The remaining mathematical content of S is: "the feed phase
hands over an all-negative, cap-dominated, debt-edge state" — nothing
about the autonomous evolution remains unproved.

**The limit constant (PROVED closed form).** The E-lin program's constant
is `1 + g + eps*` where `eps* = 2(1-g-g^2)/(3+2g)` (Theorem B). Algebra:

```
1 + g + eps*  =  (5 + 3g)/(3 + 2g),      d/dg [(5+3g)/(3+2g)] = -1/(3+2g)^2 ,
```

strictly decreasing; with the fresh certified sandwich
`627035/2^20 < g < 627036/2^20`:

```
(5+3g)/(3+2g)  in  (1.6191617801..., 1.6191618342...)      (width 5.4e-8),
```

exact endpoint Fractions in `amm12592_S_cone_constants_boxeph.json`. If the
program closes at every eps > eps* (S(eps) via ENTRY*(eps) for a sequence
eps decreasing to eps*), then `C* <= (5+3g)/(3+2g) = 1.61916...`,
UNCONDITIONALLY the target of record for this route (vs the current
unconditional C* <= 2 and the conditional 1.6292382 under S(1/32)).

## 8. Hazards honored

- Rule/search negatives never prove infeasibility; none are used. The cone
  theorem quantifies over entry STATES; ENTRY* failures at R <= 8192 are
  facts about margins, not about S (which is separately verified there).
- Every lemma carries a machine certificate (L1–L4, 600 exact trials each,
  independent implementation for the clamp); the instrumented runner is
  certified bit-identical to the fast engine before any use at new scales.
- Quantifiers: Theorem S-cone is per-(R, D0) with all constants explicit;
  the a-priori BUD corollary is exact-rationally certified for
  eps in {1/32, 1/16}, dyadic R >= 16384; nothing is claimed for
  eps < eps* (feed-phase survival is Theorem B's window).
- The G_L layer slack and the K1R Beatty count are worst-case bounds
  (measured cell-1 lifetimes are ~5x shorter); conservativeness costs only
  in BUD, which has Theta(R) slack at R >= 16384.
- Floats appear only in display fields of ledgers; every decision
  (certificates, theorem hypotheses, comparisons) is int/Fraction-exact.

## 9. Status ledger

- **PROVED (new):** S1 magnitude closed form (subsumes C-N/C-A one-step);
  S2 comparison principle; S3 unconditional deadline (a_0 <= R-2, cell-0
  clock, i_pf <= (R-2)/2); S4 cell-1 spill criterion + layer bookkeeping
  (the edge-law mechanism); Theorem S-cone (H1 ∧ H2 ∧ H3 ∧ H4* ∧ BUD =>
  S(R)), with self-propagation of H4*, permanent top-down extinction,
  explicit K1R/K2R clocks; a-priori BUD for dyadic R >= 16384 at eps in
  {1/32, 1/16}; the closed form 1 + g + eps* = (5+3g)/(3+2g) with
  certified bracket.
- **VERIFIED-exact (new):** runner certification (4 configs, bit-identical
  incl. death); 16-run sweep reproducing D2 captures exactly; zero
  re-ignitions/top-down extinction at all scales; ENTRY* margin table with
  the one-bit H4* frontier at t = 2; fresh g sandwich; R = 32768 runs
  [this session]; corner/robustness comparison certificates [this
  session].
- **HYPOTHESIS (one item, replacing S):** ENTRY*(eps) — the one-row
  feed-end condition, for dyadic R >= 65536.

The invariant-cone program of the task brief: step (1) cone defined; step
(2) preservation PROVED (H4* self-propagation); step (4) capture PROVED
(clocks + budget); step (3) entry = ENTRY*, the single remaining
hypothesis, now finite-checkable per R and machine-verified at every scale
where the theorem needs it and computation can reach.
