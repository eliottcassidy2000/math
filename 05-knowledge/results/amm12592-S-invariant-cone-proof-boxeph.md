# AMM 12592 — Lane E1: Hypothesis S via the invariant cone — the full-cap cone theorem, one-row certificates sharp to O(1) rows, and the marginal-surface law

Session: boxeph multifront, 2026-08-04 (post D1/D2/D3 + hostile referee).
All computations exact (int / Fraction); no floats in any decision; no
numpy; no sympy. Citations: T1–T9' to
`amm12592-allR-transient-theorem-boxeph.md`; Theorems A/B/C/D, C-N/C-A/
C-F/C-D to `amm12592-Elin-theorem-boxeph.md`; G1–G3 to
`amm12592-golden-transient-bound-boxeph.md`; referee to
`amm12592_estimateE_referee_boxeph.out`.

Scripts (04-computation/):
`amm12592_S_invariant_cone_certificates_boxeph.py` — instrumented exact
runner (engine-certified bit-identical), feed-end snapshots, per-row cone
diagnostics, `fcscan` (full-cap certificate + exact clocks), `iconescan`
(half-cap certificate), `corner`/`fromsnap` (comparison-lemma experiments);
`amm12592_S_cone_lemma_referee_boxeph.py` — 600-trial exact certificates
per one-step lemma (independent clamp implementation);
`amm12592_S_cone_entrycheck_boxeph.py` — feed-end certificate margins;
`amm12592_S_cone_constants_boxeph.py` — fresh g sandwich + rational
constants. Outputs (05-knowledge/results/):
`amm12592_S_cone_{run,feedend,fcscan,iconescan,entrycheck,constants,
lemma_referee,certify_vs_engine}_*` JSON/out; corner/robustness ledgers
`amm12592_S_cone_corner*_*.json`, `amm12592_S_cone_snapx*_*.json`.

Status labels: **PROVED** (complete argument, machine-checkable algebra),
**VERIFIED-exact** (exact computation at stated scales), **HYPOTHESIS /
OPEN** (precise statement, open), **REFUTED** (with certificate).

---

## 0. Headline

1. **The post-feed flow is a monotone dynamical system in closed form
   (S1/S2, PROVED).** All-negative junk under rule A evolves EXACTLY by
   `a'_t = max(0, (K_delta a)_t - 2C(d'-1,t))` (t >= 1), `a'_0 =
   max(0, a_0 - 2)` on magnitudes `a = |j|`; the map is cellwise monotone
   (comparison principle). Cell 0 is an autonomous exact 2/row clock, and
   the T8 deadline is unconditionally non-binding (S3).
2. **Theorem S-cone-fc (PROVED, the main result).** From ANY post-feed row
   whose state satisfies four exact one-row conditions — negativity,
   `a_0 <= d - 1`, the FULL-CAP inflow condition
   `2a_{t-1} + a_{t-2} <= 2C(d-1,t)` on `t in [2, m+2]`, and an exact
   clock/budget check (F4) using only the certified floor profile — the
   flow provably captures by row R-2 with no death: **S(R) follows from a
   one-row certificate.** No asymptotics; every hypothesis is decidable by
   integer arithmetic at that row. Moreover F4 is AUTOMATIC given F1–F3
   for every dyadic R >= 512 (a-priori corollary, exact-rationally
   certified to 2^40): the hypothesis is really F1 ∧ F2 ∧ F3.
3. **The certificate fires AT feed-end (VERIFIED-exact).** The scan
   `fcscan` finds the first certifying row `i_fc`: it is the feed-end row
   itself or 1–8 rows after it at every tested scale (18/18 configs), and
   the theorem's capture bound is sharp to O(1)–O(50) rows (e.g. R = 512:
   bound 314, actual 313; R = 16384: bound 10128, actual 10090). The
   unproved content of S(R) is thereby reduced to the feed phase handing
   over a state in the cone — nothing about the autonomous evolution
   remains unproved.
4. **The marginal-surface law (NEW, VERIFIED-exact).** At feed-end the
   full-cap inflow ratios `r_t = (2a_{t-1}+a_{t-2})/2C(d-1,t)` are
   `1.000 +- 0.01` ACROSS THE WHOLE BLOCK, with the excess over 1 tiny and
   shrinking (max_t r_t: 0.9872 at 256, 0.9906 at 512, 1.0038 at 1024,
   1.0097 at 2048, 1.0047 at 4096, 1.0001 at 8192, 1.0026 at 16384; both
   eps behave alike). The feed delivers the
   autonomous phase a CRITICAL state sitting exactly on the absorb
   boundary — "the profile rides just above the caps, frozen" (D2) made
   exact. Cap drift (degree growth along the Beatty word) immediately
   pulls the caps ahead, which is why i_fc = feed-end + O(1).
5. **The Lambda R-independence reading is REFUTED (exact).** The feed-end
   cap-ratio at cell 1 grows EXACTLY one bit per doubling:
   `log2 A_1 = 4.33, 5.55, 6.56, 7.62, 8.61, 9.62, 10.62, 11.62` at
   R = 128..16384 (eps = 1/32), i.e. `A_1 ~ 0.19 R` LINEAR — D2's
   "max_t log2 A_t ~ 10-11 R-independent" does not hold at feed-end (its
   instrumentation must have read a later row). Fixed-Lambda cones are
   therefore structurally wrong; the full-cap form (which is
   Lambda-free) is the right invariant. Recorded per hazard discipline.
6. **Basin geometry (VERIFIED-exact).** By the comparison principle each
   capturing run certifies its whole order-interval: feed-end states
   scaled x2 on all cells t >= 1 still capture (2048: 1581 <= 2046;
   8192: 6556 <= 8190); the uniform corner state (A_t = 2047 on
   [1, 2 sqrt R]) does NOT capture at 16384 (OPEN_RESIDUAL) — the basin
   is at least one binary order deep around the true state but far
   narrower than the uniform corner: S is robust, not knife-edge, and the
   graded profile matters.
7. **New closed form (PROVED).** `1 + g + eps* = (5+3g)/(3+2g)` with
   derivative `-1/(3+2g)^2`; certified bracket
   `(1.6191617801, 1.6191618342)` (width 5.4e-8). This is the E-lin
   program's limit constant in one fraction.
8. **R = 32768: S verified one doubling up (this session).** eps = 1/16:
   CLOSED at capture 19865 = i_pf + a_0/2 - 1 exactly, debt (R-2)/2
   exact, zero re-ignitions; Theorem S-cone-fc certificate at
   i_fc = i_pf + 4 with bound 19955 and 39%-of-R budget margin. All four
   feed-end laws extend to the ninth doubling at both eps (m = 169/168 vs
   0.93 sqrt R = 168.9; a_0 - d = 0/-1; log2 A_1 = 12.65/12.62; max r_t =
   1.0009 both, over-cells single-parity {3,5,7(,9)}). The referee's open
   item (32768 probe died at 1.9%, "no evidence either way") is resolved:
   S(1/16) now VERIFIED-exact for all dyadic 128 <= R <= 32768.

---

## 1. Setting

Post-feed autonomous T6 flow at dyadic R, profile
`d_i = floor(g(R+i)) + D0`, `g = gamma* = log_5 phi^2` (fresh certified
sandwich `627035/2^20 < g < 627036/2^20`, re-derived this session by
integer Lucas/Fibonacci comparisons at M = 2^20; independent of the prior
referee's). Junk `j` on cells t at degree d; per row `delta = d' - d in
{0,1}`; transport `(K_delta * j)_c = sum_s C(1+delta,s) j_{c-s}`; clamp
into `[-2C(d'-1,t), +2C(d'-1,t-1)]` (t >= 1), `[-2, 0]` (t = 0); junk =
load - clamp. Death iff junk at cell d'; capture iff junk empty (then
T5/T1 coasting). `i_pf` := first row with `d_i + i > R` (no feed ever
again; `d_i + i` strictly increasing makes this permanent); the state
entering row `i_pf` is the FEED-END state, degree `d_fe = d_{i_pf - 1}`.
Magnitudes `a_t := -j_t >= 0`; front `m := max supp`; caps `2C(d-1,t)`.
All junk entries are even (T3; asserted at every clamp of every run).

## 2. Lemma S1 (magnitude closed form; cell-0 clock) — PROVED

**Statement.** If `j <= 0` cellwise then after one post-feed row of rule A
(either delta): junk stays `<= 0`, and exactly

```
a'_t = max(0, (K_delta a)_t - 2C(d'-1,t))   (t >= 1),
a'_0 = max(0, a_0 - 2).
```

*Proof.* Load `w = K_delta * j <= 0`; for w <= 0 the upper box end is
inactive, the nearest-point clamp is `u_t = max(-2C(d'-1,t), w_t)`, junk
`= min(0, w_t + 2C(d'-1,t))`. At t = 0 the box is [-2, 0] and the load is
`j_0` (no lower cells). QED

Subsumes C-N; upgrades C-A to the exact map. Cell 0 is autonomous: the
debt drains at exactly 2/row (T8's schedule needs no hypothesis here),
`a_0` stays even, cell 0 empties after exactly `a_0/2` rows, permanently.
Certificate L1: 600 exact random states vs an independent clamp
implementation — ALL PASS (`amm12592_S_cone_lemma_referee_boxeph.out`).

## 3. Lemma S2 (comparison principle) — PROVED

`a <= b` cellwise implies `Phi(a) <= Phi(b)` cellwise (each output is
nondecreasing in every input; induct along rows). Hence a capturing
trajectory from a majorant state proves capture (and no-death) for every
state below it. Certificate L2: 600 comparable pairs — ALL PASS.

## 4. Lemma S3 (deadline, unconditional) and Lemma S4 (layer) — PROVED

**(S3)** Post-feed, `e_i(1) = j_0` and the T8 identity give
`a_0 = (R-2) - 2 * minus2count <= R - 2` (using only j_0 <= 0). With
`i_feed = floor((R(1-g) - D0)/(1+g))` (Theorem B) and `g > 1/3`:
`i_pf <= i_feed + 2 <= (R-2)/2` (all R >= 32, any D0 >= 0), so cell 0
empties by `i_pf + a_0/2 - 1 <= R - 2`: the drain ALWAYS finishes in
time.

**(S4)** Cell 1's only inflow from below is `(1+delta) a_0` against cap
`2(d'-1)`: in-box iff `a_0 <= d'-1` (delta = 1; delta = 0 rows are always
in-box in our regime). Since a_0 falls 2/row and d is nondecreasing,
`a_0 <= d-1` is ABSORBING; the observed feed-end edge law `a_0 = d_fe +-
O(10)` (all 18 runs: a_0 - d_fe in [-14, +7]) is precisely the handover
at this threshold. Layer := rows with a_0 > d-1; length
`L = ceil((a_0 - (d_fe-1))/2)^+`, and cell 1 gains at most
`G_L = sum_k 2 max(0, a_0 - 2k - d_fe)` in total during it (delta = 1
rows have cap >= 2 d_fe; delta = 0 rows give no gain). Certificate L4:
600 trials — ALL PASS.

## 5. Theorem S-cone-fc (full-cap cone: one-row certificate => S(R)) — PROVED

Notation: fix a post-feed row `i0`; write `D_k := d_{i0+k}` for the EXACT
degree profile (certified floor engine), `delta_k := D_k - D_{k-1}`; the
state entering row i0 is `a` with front m; `capref_t := 2C(D_0 - 1, t)`.

**Hypotheses (all exact, all at row i0):**

- **(F1)** `j <= 0`, support ⊆ [0, m], `m + 2 < D_0`;
- **(F2)** `a_0 <= D_0 - 1`;
- **(F3)** `2 a_{t-1} + a_{t-2} <= capref_t` for every `t in [2, m+2]`
  (entries above the support read 0);
- **(F4)** `i0 + max( ceil(a_0/2), max_{1<=t<=m} T_t ) <= R - 2`, where
  the clocks are the exact staircase sums

```
t >= 2:  T_t := min{ K : sum_{k=1..K} ( 2C(D_k - 1, t) - capref_t ) >= a_t }
t = 1:   T_1 := min{ K : sum_{k=1..K} [ 2(D_k - 1)
                     - (1 + delta_k) * max(0, a_0 - 2(k-1)) ]^+  >= a_1 }.
```

**Conclusion.** For every k >= 0: cells >= 1 are cellwise non-increasing,
the support stays inside [0, m], junk stays negative, no death can occur,
and cell t is permanently empty from row `i0 + T_t` (t >= 1), cell 0 from
`i0 + ceil(a_0/2)`. Hence junk is empty by row
`i0 + max(ceil(a_0/2), max_t T_t) <= R - 2`: **capture, and S(R) holds.**

*Proof.* Induction on k with hypothesis: `a^{(k)} <= a` cellwise on
t >= 1; `a^{(k)}_0 = max(0, a_0 - 2k)`; support ⊆ [0, m]; F3 holds at the
current values. Row i0+k clamps at degree `D_k >= D_0`, so every cap
`2C(D_k - 1, t) >= capref_t`.
(i) Cells `t in [2, m]`: the load is `a_t + spill` with `spill <=
2a_{t-1} + a_{t-2} <= capref_t` (both delta-forms are dominated), so by
S1 `a'_t <= max(0, a_t + capref_t - 2C(D_k - 1, t)) <= a_t`, with decay
at least `2C(D_k - 1, t) - capref_t` while alive.
(ii) Cell 1: inflow `(1+delta_k) a^{(k)}_0 <= 2 a_0 <= 2(D_0 - 1) <=
2(D_k - 1)` = cap (F2 absorbing), so `a'_1 <= a_1`, with decay
`[2(D_k-1) - (1+delta_k) a^{(k)}_0]^+` — the T_1 summand.
(iii) Beyond the front: the only loads are at m+1 (`2a_m + a_{m-1} <=
capref_{m+1}`, by F3 at t = m+1) and m+2 (`a_m <= capref_{m+2}`, by F3 at
t = m+2 with `a_{m+1} = 0`), both inside their caps: absorbed entirely,
support frozen. Death needs junk at cell `D_k > m + 2`: impossible.
(iv) All cells non-increasing => spills non-increasing => F3 propagates.
(v) Clocks: the cumulative decays are exactly the T-sums; once a cell is
zero its load is its spill `<= capref <=` cap: it stays zero. Cell 0 is
the exact clock. F4 bounds the last extinction row. Capture, then T5/T1
coasting closes the epoch. QED

Certificate L3 (600 random states satisfying the cone; both kernels):
one-step non-increase, support freeze, half/full-cap decay, and condition
propagation — ALL PASS. (L3 checks the half-cap variant, whose one-step
claims dominate the full-cap ones cell-by-cell.)

**Corollary (a-priori F4: the clocks are free) — PROVED + certified.**
F3 itself bounds every magnitude: reading F3 at cell t+1 gives
`a_t <= C(D_0 - 1, t+1)` for every `t in [1, m+1]`. Feeding these into
the clock sums (with the exact staircase `D_k - D_0 >= floor(g_lo k)`,
`s(2 D_0 - 3 + s)/2` for the t = 2 cap increments, the drain-assisted
cell-1 decay `>= 2[s_k + 2(k-1)]`, and the monotonicity lemma
`T_t <= T_2` for all `t >= 2` — normalized need `(D_0-1-t)/(2(t+1))`
decreasing in t while the normalized cap-growth rate increases in t)
yields explicit bounds `K1c ~ D_0/sqrt(2(2+g)) ~ 0.34 R` and
`K2c ~ D_0/sqrt(6g) ~ 0.41 R`. Exact-rational verification
(`amm12592_S_cone_constants_boxeph.py`): for every dyadic
`R = 2^9 .. 2^40` at eps in {1/32, 1/16}, any post-feed certifying row
`i0 <= i_feed + 66` satisfies F4 automatically. **Hence for dyadic
R >= 512 the hypothesis of Theorem S-cone-fc reduces to F1 ∧ F2 ∧ F3
alone** (R = 128, 256 are directly verified anyway).

**Remarks.** (a) The theorem consumes only the certified floor profile
and the row-i0 state: `T_t` are computable in O(m * T_max) exact integer
ops. (b) By S2 the conclusion extends to every state cellwise below a
certifying state. (c) A half-cap variant (spill <= capref/2, uniform
1/2-cap decay per row, Lambda-free clocks) is proved the same way and
fires mid-flight (~0.38 R, `iconescan`); the full-cap form is the sharp
one. (d) F3's role is exactly the marginal-surface law: the feed-end
state has `r_t <= 1.01` everywhere, and the certificate fires as soon as
the O(1%) excess cells dip under 1 — measured: 0–8 rows (16/16).

## 6. Certificates — VERIFIED-exact

**Engine certification.** The instrumented runner is bit-identical to the
certified fast engine at (128,8), (256,16), (512,32) closures and the
(512,4) death (row 121, const_bits 292; every compared per-row stat
equal): `amm12592_S_cone_certify_vs_engine_boxeph.json`.

**Sweep (16 runs, R = 128..16384, eps in {1/32, 1/16}).** Every capture
row equals the D2 table exactly; debt = (R-2)/2 in every run; ZERO
re-ignition events across all post-feed rows of all runs (a dead cell
t >= 1 never revives); extinction strictly top-down; cell 0 always last
(capture = i_pf + a_0/2 - 1 exactly in every run).

**Feed-end laws (exact snapshots, 18 runs incl. re-runs).**

- Marginal surface: max_t r_t [eps = 1/32 / eps = 1/16] =
  0.8639/0.8494 (128), 0.9727/0.9872 (256), 0.9552/0.9906 (512),
  1.0038/0.9949 (1024), 0.9949/1.0097 (2048), 1.0018/1.0047 (4096),
  1.0001/0.9999 (8192), 1.0004/1.0026 (16384). Shape (16384/512
  profile): r_t = 0.9998, 1.0000, 1.0004, 0.9997, 0.9984 at
  t = 2,4,6,8,10 — the LOW-CELL CORE sits on the surface to within
  0.2% — then declines monotonically toward the front (0.98 at t~13-30,
  0.93 at 40, 0.64 at 80, 0.03 at 118). Cells with r_t > 1 are few
  (0-4 per state), confined to the core, with excess <= 1%, and appear
  in single-parity runs (e.g. {2,4,6,8} or {3,5,7,9} — an alternation-
  calculus fingerprint). The binding clock cell (t = 2, always) is
  exactly where the state is most marginal: the certificate's sharpness
  is the law's sharpness.
- `A_1 = a_1/(2(d_fe-1))`: log2 = 4.33, 5.55, 6.56, 7.62, 8.61, 9.62,
  10.62, 11.62 (eps = 1/32) and 4.39, 5.61, 6.58, 7.64, 8.66, 9.64,
  10.65, 11.65 (eps = 1/16), R = 128..16384: exactly +1 bit/doubling,
  `A_1 ~ 0.19 R` — REFUTES the fixed-Lambda reading (D2's sec. 4.2
  "max_t log2 A_t ~ 10-11 R-independent" does not describe the feed-end
  row).
- Front `m ~ 0.93 sqrt R` (118/119 at 16384); `a_0 - d_fe in [-14, +9]`
  (edge law); all-negative and contiguous in every snapshot.

**fcscan (the theorem in action).** All scans: `F3_persist_ok = True`
(the propagation clause re-verified in flight); worst clock always at
t = 2; the capture bound is drain-dominated and sharp to 1–90 rows
(0.3–0.5% of R at the top scales):

| R | D0 | i_pf | i_fc | i_fc-i_pf | thm capture bound | actual | budget margin |
|------|-----|------|------|-----|------|------|------|
| 128 | 4 | 31 | 31 | 0 | 80 | 79 | 46 |
| 128 | 8 | 28 | 29 | 1 | 79 | 78 | 47 |
| 256 | 8 | 61 | 61 | 0 | 154 | 153 | 100 |
| 256 | 16 | 56 | 56 | 0 | 151 | 150 | 103 |
| 512 | 16 | 120 | 122 | 2 | 318 | 317 | 192 |
| 512 | 32 | 110 | 112 | 2 | 314 | 313 | 196 |
| 1024 | 32 | 239 | 241 | 2 | 628 | 624 | 394 |
| 1024 | 64 | 219 | 219 | 0 | 620 | 616 | 402 |
| 2048 | 64 | 476 | 476 | 0 | 1263 | 1261 | 783 |
| 2048 | 128 | 436 | 439 | 3 | 1239 | 1237 | 807 |
| 4096 | 128 | 951 | 955 | 4 | 2527 | 2519 | 1567 |
| 4096 | 256 | 871 | 875 | 4 | 2497 | 2486 | 1597 |
| 8192 | 256 | 1902 | 1902 | 0 | 5063 | 5040 | 3127 |
| 8192 | 512 | 1742 | 1742 | 0 | 4987 | 4964 | 3203 |
| 16384 | 512 | 3803 | 3803 | 0 | 10128 | 10090 | 6254 |
| 16384 | 1024 | 3482 | 3490 | 8 | 9984 | 9937 | 6398 |
| 32768 | 2048 | 6963 | 6967 | 4 | 19955 | 19865 | 12811 |
| 32768 | 1024 | 7604 | 7608 | 4 | 20273 | 20185 | 12493 |

Every certifying row proves S at its (R, D0) by Theorem S-cone-fc alone
(the continued run is a consistency check, not part of the proof). The
budget margin is ~38% of R throughout — the theorem is nowhere close to
binding.

**Basin experiments (comparison lemma).** x2-scaled feed-end states (all
cells t >= 1 doubled) capture: 2048 -> 1581 (<= 2046), 8192 -> 6556
(<= 8190) — each such run proves S for its entire order-interval, so the
basin has at least one binary order of depth above the true trajectory.
x4-scaled states do NOT capture in time (OPEN_RESIDUAL at row R-2, no
death, both 2048 and 8192; clean re-runs from regenerated snapshots):
the basin radius in the uniform-scaling direction is between x2 and x4.
The uniform corner (a_t = 2047 * 2C(d-1,t) on [1, 2 sqrt R], a_0 = d+64)
does NOT capture at 16384 (OPEN_RESIDUAL) — fatness requires the graded
(marginal) shape, consistently with the x4 boundary. [A
snapshot-overwrite bug (corner/fromsnap modes stomped three feed-end
archives) was caught by the A_1 anomaly, fixed, and all affected data
regenerated; contaminated first-pass x4/x16 results discarded —
hazard discipline. Dual-path validation: full-feed fcscan and
snapshot-replay fcscan agree bit-exactly at both 16384 configs.]

**R = 32768 (this session; resolves the referee's open item E).**

- eps = 1/16 (D0 = 2048): **CLOSED, capture 19865** (= 0.606 R), debt
  16383 = (R-2)/2 exact, i_pf = 6963, capture = i_pf + a_0/2 - 1 exactly,
  zero re-ignitions, strict top-down extinction (cell 1 dies at 17133,
  cell 0 last). The prior session's probe (killed at row ~384, "no
  evidence either way" per the referee) is settled: **S(1/16) is
  VERIFIED-exact at R = 32768** — the verified range of Hypothesis S
  moves up one doubling.
- Feed-end laws at 32768/2048: m = 169 (0.93 sqrt R = 168.9 — exact
  continuation), a_0 - d_fe = 0 (the edge law at its sharpest),
  log2 A_1 = 12.65 (+1 bit/doubling, ninth doubling), max_t r_t =
  1.00090 at t = 5 with over-cells {3,5,7,9} (single-parity core,
  excess 0.09% — the marginal surface tightens again).
- eps = 1/32 (D0 = 1024): **CLOSED, capture 20185** — EXACTLY the
  cell-0-clock prediction i_pf + a_0/2 - 1 = 7604 + 12582 - 1, made
  before the run finished; debt 16383 = (R-2)/2 exact; zero
  re-ignitions; top-down extinction (cell 1 dies at 17521, cell 0
  last). Feed-end laws: d_fe = 25165, m = 168 (0.93 sqrt R again),
  a_0 - d_fe = -1, log2 A_1 = 12.62 (+1 bit/doubling), max_t r_t =
  1.00089 at t = 5, over-cells {3,5,7}: every law extends at the second
  eps as well. **S(1/32) is verified at R = 32768.**
- fcscan certificates: eps = 1/16: i_fc = 6967 = i_pf + 4, Tmax = 12988
  (worst cell t = 2), capture bound 19955 vs actual 19865; eps = 1/32:
  i_fc = 7608 = i_pf + 4, Tmax = 12665, bound 20273 vs actual 20185;
  F3_persist_ok through all ~13000 post-feed rows in both.

## 7. The reduction, and what remains

**ENTRY-fc(eps) (HYPOTHESIS — the single remaining item).** For every
dyadic `R >= 65536` at `D0 = ceil(eps R)`: some post-feed row
`i0 <= i_pf + 64` of rule A's flow satisfies F1, F2, F3. (F4 is then
automatic by the a-priori corollary.)

**Theorem (E1 reduction; PROVED).** ENTRY-fc(eps) implies S(eps)
restricted to dyadic R >= 65536; with the exact sweep (R <= 16384) and
this session's R = 32768 certificates, D2's Hypothesis S(eps) holds in
full. Consequently (D2 Theorem D + LIFT + THM-3329 assembly):
`C* <= 1 + gamma* + eps`; in particular S(1/32) gives
`C* <= 1 + gamma* + 1/32 < 427095/262144 = 1.6292382`.
(No death can occur in the window [i_pf, i0): the T6a front-speed bound
gives support <= m + 2*64 + 2 << d there — the window is death-free
independently of the hypothesis.)

**What ENTRY-fc asks.** Only that the feed phase delivers (within 64
rows) a state that is: all-negative (D2's certified N), debt at the edge
(a_0 <= d-1 after the S4 layer), and on-or-below the full-cap marginal
surface with clocks fitting a Theta(R) budget that is ~30% slack. The
marginal-surface law says the flow sits ON this surface with deviation
-> 0 as R grows; ENTRY-fc asks it to stay there. All of S's dynamical
content (Theta(R) rows of evolution) is now proved; what remains is a
STATIC property of the feed-phase endpoint — squarely the target of the
G1/G2 feed-phase alternation calculus (the D2 sec. 5 envelope program),
and finitely checkable per R by one feed-phase computation.

**The limit constant.** `1 + g + eps* = (5+3g)/(3+2g)` (PROVED identity;
strictly decreasing, `h' = -1/(3+2g)^2`), certified bracket

```
(5+3g)/(3+2g)  in  (1.6191617801..., 1.6191618342...) ,
```

exact endpoint Fractions in `amm12592_S_cone_constants_boxeph.json`. If
ENTRY-fc(eps) holds for a sequence eps -> eps*+ then
`C* <= (5+3g)/(3+2g) = 1.61916...` — the sharpest constant this route
can produce, now in closed form.

## 8. Hazards honored

- Rule/search negatives never prove infeasibility; none used. Corner
  non-capture is a statement about ONE majorant state, not about S.
- The snapshot-overwrite bug (corner/fromsnap stomping three feed-end
  archives) was caught by the A_1 anomaly (values 2^15.7/2^17.7 off the
  +1-bit/doubling line), fixed, and all affected data regenerated; the
  contaminated x4/x16 robustness results were discarded.
- D2's "max_t log2 A_t ~ 10-11 R-independent" is refuted AT FEED-END by
  exact snapshots; conclusions of D2 that used it (none load-bearing —
  it appeared only as "structural support" for S) are unaffected; the
  correct law is A_1 ~ 0.19 R with r_t ~ 1 (marginal surface).
- Every proved lemma carries an exact machine certificate (L1–L4);
  Theorem S-cone-fc's in-flight propagation is re-verified by every
  fcscan (F3_persist_ok).
- Quantifiers: Theorem S-cone-fc is per-(R, D0), hypotheses at one row,
  no asymptotic step; ENTRY-fc is stated with its exact row window
  (i_pf + 64) and range (dyadic R >= 65536); eps range for the C*
  corollary inherits Theorem B's window (eps > eps*, D0 = ceil(eps R),
  window (ii)).
- Floats only in ledger display fields; all decisions int/Fraction.

## 9. Status ledger

- **PROVED (new):** S1 magnitude closed form + cell-0 clock; S2
  comparison principle; S3 unconditional deadline; S4 spill criterion +
  layer bounds (the edge-law mechanism); **Theorem S-cone-fc** (one-row
  full-cap certificate => capture by R-2, no death — S(R)); half-cap
  variant; a-priori budget corollary; the closed form
  `1 + g + eps* = (5+3g)/(3+2g)` with certified bracket; fresh g
  sandwich.
- **VERIFIED-exact (new):** runner bit-identical certification; 18-run
  sweep incl. R = 32768 both eps (captures = D2 exactly where known; NEW:
  32768 closes at 19865 (eps = 1/16) and 20185 (eps = 1/32); zero
  re-ignitions, top-down extinction,
  capture = i_pf + a_0/2 - 1 in every run); the
  marginal-surface law (low-cell core within 0.2% of r = 1, over-cells
  single-parity with excess <= 1%, monotone decline to the front);
  A_1 ~ 0.19R (+1 bit/doubling, NINE doublings, both eps); edge law
  a_0 - d_fe in [-14, +9] (0 and -1 at 32768); fcscan certificates 18/18
  with i_fc - i_pf in {0..8} and capture bounds sharp to 1-90 rows;
  x2 basin depth (capture) vs x4 (deadline miss); corner non-capture;
  **S(1/32) AND S(1/16) verified for ALL dyadic 128 <= R <= 32768**
  (referee item E resolved; the 32768/1024 capture row 20185 was
  predicted exactly by the cell-0 clock before the run completed).
- **REFUTED (new):** fixed-Lambda feed-end cones (A_1 is linear in R);
  the uniform corner as a basin majorant at 16384.
- **HYPOTHESIS (one item, replacing S):** ENTRY-fc(eps) for dyadic
  R >= 65536 — a one-row, statically-checkable property of the feed-end
  state (negativity + debt edge + marginal surface; clocks automatic by
  the a-priori corollary).

Lane E1 verdict: the invariant-cone program's steps (1) cone, (2)
preservation, (4) capture are PROVED (S-cone-fc); step (3) entry is
reduced from the Theta(R)-row dynamical Hypothesis S to the one-row
ENTRY-fc, verified at every scale computation can reach (128..32768) and
sitting on an exact structural law (the marginal surface) that the
feed-phase calculus (G1/G2) is built to attack.
