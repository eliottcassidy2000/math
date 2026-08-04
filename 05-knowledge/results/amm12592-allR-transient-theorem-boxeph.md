# AMM 12592 — Angle B2: the transient theorem for attractor rule A (exact skeleton)

Session: boxeph multifront, 2026-08-03. Post-THM-3329. All computations exact
(int; no floats in any decision; sympy used only to evaluate one analytic
constant, flagged where used).

Scripts (04-computation/):
`amm12592_transient_error_dynamics_boxeph.py` (dynamics + conjugacy + traces),
`amm12592_transient_initial_decode_boxeph.py` (closed-form initial decode +
margin windows), `amm12592_transient_lemma_certificates_boxeph.py` (T1/T6/T7/
master-equation finite certificates). Outputs (05-knowledge/results/):
`amm12592_transient_conjugacy_{small,R256,R512}_boxeph.json`,
`amm12592_transient_trace_R{128,256,512,1024}_D0*_boxeph.json`,
`amm12592_transient_initial_decode_boxeph.json`, `amm12592_transient_*.out`.
New witness: `04-computation/amm12592_witness_R512_ruleA_D0_5_boxeph.json`.

---

## 0. Headline

1. **Rule A is exactly conjugate to a clamped-Pascal junk flow** with explicit
   initial data, explicit absorbing boxes, and boundary feed on two cells
   (Lemmas T2, T6). All structure of the transient (S5's measured avalanche)
   becomes provable inequalities in these coordinates.
2. **Closed form for the initial cell decode** (Lemma T4): the row-0 load
   profile of 2G_R is three binomials; the initially safe window is the top
   `~32.77%` of cells (constant `1-tau*`, an explicit entropy-equation root).
3. **Death-delay theorem** (T6b, proved): the rule cannot die before row
   `d_0 - t_lo + 1 ~ 0.196 R`. Measured deaths (R=256, 512) occur within 5–11
   rows of this bound — the bound is near sharp, and every death occurs in the
   feed phase `i <= i_feed ~ 0.2516 R`.
4. **Ballot-debt theorem** (T8, proved): capture cannot occur before row
   `(R-2)/2 - 1`; the c0-word is the scheduled ballot path, and the endgame is a
   forced cone-boundary ride draining the residual debt at exactly 2/row
   (drain Lemma T7). Verified: the closures pay the debt exactly —
   `#{-2 rows} = (R-2)/2` in every winning trace.
5. **New record slack at R = 512: D0 = 5 closes** (deterministic rule A;
   witness saved + verified; D0 = 0..4 die at rows 107/110/113/116/121).
   With THM-3329's assembly this improves the n <= 1023 attainment to
   `T(n) <= n + 1 + floor(gamma* n) + 5`.
6. The all-R theorem is reduced to ONE precisely stated estimate
   (**Estimate E**, Section 9: front stall + junk annihilation), with exact
   finite-margin data at R = 128, 256, 512 quantifying how far the plain rule
   is from satisfying it at each slack level.

Status labels used below: **PROVED** (complete argument, machine-checkable
algebra), **VERIFIED-finite** (exact computation at stated R), **CONJECTURED**
(precise statement, open).

---

## 1. Setup (conventions fixed once)

`p = x`, `q = 1-x`, `E_m = -1 + x + ... + x^m` (`E_{-1} = 0`),
`gamma* = log(phi)/log(sqrt 5)`. Epoch problem (*) at slack D0:
find blocks `Delta_i` admissible at `d_i = floor(gamma*(R+i)) + D0`
(Bernstein cells `Delta = sum_t delta_t x^{d-t} q^t`, `|delta_t| <= C(d,t)`,
`delta_t == C(d,t) mod 2`) with `sum_{i<R} x^i Delta_i = q^{R-1}`.

Rule A (plain, tozero): `sigma_{-1} = q^{R-1}`;
`ideal_i = sigma_{i-1} - x E_{R-2-i}`;
`Delta_i = AdmClamp_{d_i}(trunc_{d_i}(ideal_i))`;
die unless `x | sigma_{i-1} - Delta_i`; `sigma_i = (sigma_{i-1}-Delta_i)/x`;
closure iff `sigma_{R-1} = 0`.

Cell decode of `P = sum_j a_j x^j` at degree d (toolbox `poly_to_bern`):
`w_t = sum_j a_j C(d-j, t)`; encode/decode are mutually inverse on deg <= d.

## 2. Lemma T1 (coasting; PROVED, known)

If `sigma_{i-1} = E_{R-1-i}` then `ideal_i = 2x-1`, whose degree-d decode is
the ballot vector `b_t = C(d-1,t) - C(d-1,t-1)`; `|b_t| <= C(d,t)` with
correct parity, so no clamping, `Delta_i = 2x-1`, `sigma_i = E_{R-2-i}`.
Final row: `sigma_{R-2} = E_0 = -1` gives `Delta_{R-1} = -1` (decode
`-C(d,t)`, the box corner) and `sigma_{R-1} = 0`. So `E`-track states coast
to closure. (Identities: `x E_m + (2x-1) = E_{m+1}`; `b` is the decode of
`2x-1` because `2C(d-1,t) - C(d,t) = C(d-1,t) - C(d-1,t-1)`.)

## 3. Lemma T2 (error conjugacy; PROVED + VERIFIED-finite)

Define `e_i := sigma_i - E_{R-2-i}` (`E_{-1}:=0`), so
`e_{-1} = q^{R-1} - E_{R-1} = 2 G_R` and, coefficientwise,
`[x^0] 2G_R = 2`, `[x^j] 2G_R = (-1)^j C(R-1,j) - 1` for `1 <= j <= R-1`
(re-derived and machine-verified here, not taken from the state note).

**Claim.** For `0 <= i <= R-2` rule A is exactly the map
```
w   = decode_{d_i}(trunc_{d_i}(e_{i-1}))
c_t = clamp of w_t into the asymmetric even box [-2C(d_i-1,t), +2C(d_i-1,t-1)]
      (+ tozero parity fix relative to the ballot corner; see T3)
die unless [x^0](e_{i-1} - c) = 0   (equivalently e_{i-1}[0] in {0,2})
e_i = (e_{i-1} - c)/x ,      with Delta_i = (2x-1) + c
```
and for `i = R-1` the same with corner `-C(d,t)` and box `[0, +2C(d,t)]`.
Closure iff `e_{R-1} = 0`; summing the recursion recovers the master equation
`sum_{i<=R-2} x^i c_i = 2 G_R` (THM-3329 S2).

*Proof.* `trunc` and `decode` are linear and `trunc_d(e + (2x-1)) =
trunc_d(e) + (2x-1)` for `d >= 1`; the decode of `2x-1` is `b`; clamping
`b_t + w_t` into `[-C(d,t), C(d,t)]` equals `b_t +` (clamp of `w_t` into
`[-C(d,t)-b_t, C(d,t)-b_t] = [-2C(d-1,t), +2C(d-1,t-1)]`). The recursion
`e_i = (e_{i-1}-c_i)/x` is `x E_m + (2x-1) = E_{m+1}` again. The die
condition is `[x^0]`-matching; `[x^0]Delta = delta_d in {+-1}` is the forced
unit, so survival iff `e_{i-1}[0] in {0,2}`. QED

**Certificates.** Block-for-block equality of the two implementations
(independent code paths) at R = 8, 16, 32, 64, 128 (D0=0), plus full
admissibility + epoch identity of the error-dynamics output
(`amm12592_transient_conjugacy_small_boxeph.json`: all True); equality against
the STORED rule-A witnesses at R = 256 (D0=1) and R = 512 (D0=8)
(`..._R256_...json`, `..._R512_...json`: `witness_equal: true`).

## 4. Lemma T3 (dyadic parity inertness; PROVED + VERIFIED-finite)

At `R = 2^t`, `2 G_R` is even coefficientwise (dyadic parity theorem,
re-checked). Evenness of `e` is invariant: `w` is an integer combination of
even coefficients, the box ends `-2C(d-1,t), +2C(d-1,t-1)` are even, so the
min/max clamp of an even value is even and `c` (hence `e_i = (e-c)/x`) stays
even. Moreover `b_t - C(d,t) = -2C(d-1,t-1)` is even, so after clamping the
parity condition `delta_t == C(d,t) (mod 2)` holds automatically: **the
parity-fix branch of AdmClamp never executes at dyadic R**, and the three
variants tozero/toward/away coincide. Certificate: `parity_fires = 0` in
every run at R = 8 … 512, all D0 probed.

Consequence: at dyadic epochs the rule is a pure box clamp; the transient
problem is parity-free (matches THM-3329 S3).

## 5. Lemma T4 (initial decode in closed form; PROVED + VERIFIED-finite)

For `0 <= t <= d <= R-2`:
```
w_t(R,d) := decode_d(trunc_d(2G_R))_t
          = (-1)^{d-t} C(R-2-t, d-t)  -  C(d+1, t+1)  +  2 C(d, t)
```
*Proof.* Vandermonde: `sum_j (-1)^j C(R-1,j) C(d-j,t) = [x^{d-t}](1-x)^{R-2-t}
= (-1)^{d-t} C(R-2-t, d-t)`; hockey stick: `sum_{j<=d} C(d-j,t) = C(d+1,t+1)`;
`a_0 = 2` contributes `+2C(d,t)` beyond the uniform formula. QED
Certificate: equality with the toolbox decode for
R in {6, 8, 12, 16, 25, 32, 64, 100, 128, 256} at d in {2, 3, d0-2, d0, d0+3,
R-2} (also non-dyadic R: the formula is general).

**Margin window.** With caps `[-2C(d-1,t), +2C(d-1,t-1)]`, define
`t_lo(R,D0)` = least t such that cells `t..d` are all in-box at
`d = floor(gamma* R) + D0`. Exact table
(`amm12592_transient_initial_decode_boxeph.json`), excerpt at D0 = 0:

| R | d_0 | t_lo | in-box fraction | worst overflow (bits) |
|------|------|------|------|------|
| 128 | 76 | 51 | 33.7% | 119 |
| 256 | 153 | 102 | 33.7% | 242 |
| 512 | 306 | 205 | 33.2% | 491 |
| 1024 | 612 | 411 | 32.9% | 988 |
| 8192 | 4898 | 3293 | 32.7% | 7955 |

**Limit constant** (first-order; CONJECTURED as a limit, exact at the level
of the entropy comparison): `t_lo/d -> tau*`, the root of
`(1 - g tau) H(g(1-tau)/(1-g tau)) = g H(tau)` with `g = gamma*`,
H = binary entropy: `tau* = 0.67227593818...`; in-box fraction
`1 - tau* = 0.32772406...` (measured 0.67232 at R = 8192). The death-delay
coefficient below is `c_1 = g(1-tau*) = 0.19597487...` (no nicer closed form
found; `H(g)/5 = 0.19442...` is ruled out by digit 3).

## 6. Lemma T5 (peel, one-shot absorption, capture; PROVED)

(a) If at row i no cell clamps, then `c = trunc_{d_i}(e_{i-1})` exactly
(encode/decode inverse), hence `e_i[j] = e_{i-1}[j+1]` for `j >= d_i` and
`= 0` below: a pure peel; junk polynomial `J_i := trunc(e_{i-1}) - c = 0`.
(b) If moreover `deg e_{i-1} <= d_i`, then `e_i = 0`.
(c) If `e_i = 0`, the rule coasts to closure (T1). So "capture" = first row
with `e_i = 0`; after it, zero clamping and `Delta = 2x-1` to the end.
Deciding whether a given residual is absorbable in L rows is exactly the S_L
endgame procedure (`amm12592_r128_endgame_algebra_boxeph.py`), which remains
the decision tool for Phase II below.

## 7. Lemma T6 (junk transport = clamped Pascal flow; PROVED)

Write `H_i` = part of `e_{i-1}` above degree `d_i` (the pristine tail) and
`J_i` = `trunc_{d_i}(e_{i-1}) - c_i` (junk created at row i; supported on
cells where w overflows). Then `e_i = (H_i + J_i)/x` and, with
`delta_i := d_{i+1} - d_i in {0,1}` (Beatty word of gamma*):

**(T6-kernel)** the next row's load is
`w_{i+1} = K * J_i + feed_{i+1}` where `K` is the Pascal kernel
`(1,1)` if `delta_i = 0`, `(1,2,1)` if `delta_i = 1` acting on cell index
(`x^{d-t-1} q^t = sum_s C(1+delta_i, s) x^{d_{i+1}-(t+s)} q^{t+s} / ...`
via `x + q = 1` applied `1 + delta_i` times), and

**(T6-feed)** `feed_{i+1}` = decode of the freshly exposed band of `H_i/x`,
supported on cells `{0, 1}` only (band width `1 + delta_i`), with values that
are *original coefficients of 2G_R*: the tail above the truncation line is
never modified by clamping (junk has degree <= d_i), so the coefficient
entering at row i is literally `[x^{j}] 2G_R` at `j = d + i + 1 - ...`
(position bookkeeping: coefficient j of e_i equals coefficient j+i+1 of 2G_R
whenever j >= (junk ceiling)). Feed is nonzero only while
`i <= i_feed := max{ i : d_i + i + 1 <= R-1 }`, and
`i_feed = (1-gamma*)R/(1+gamma*) + O(D0)`; exactly: i_feed(128,0) = 31,
i_feed(256,0) = 63, i_feed(512,D0=0..8) = 128..123, i_feed(1024,0) = 256.

**Corollary T6a (avalanche law, proves S5's measured bounds).** The junk
front `F(i) := max{t : J_{i,t} != 0}` satisfies
`F(i+1) <= F(i) + 1 + delta_i`, and cell amplitudes multiply by at most
`sum K = 2^{1+delta_i} <= 4` per row. (Measured: x2–4/row. Now an exact
inequality.)

**Corollary T6b (death-delay theorem).** Death at row i requires junk at
cells `>= d_{i-1} - 1` (only cells d, d-1 of `J` touch `[x^0]e`, via
`[x^1](x^{d-t}q^t)`). Since `F(0) <= t_lo - 1` (T4) and F advances at most
`1 + delta` per row while the bottom recedes `delta` per row:
```
                 no death before row   d_0 - t_lo + 1  ,
```
which by the T4 window is `>= (c_1 - o(1)) R`, `c_1 = gamma*(1-tau*) =
0.19597...`.

**Sharpness (VERIFIED-finite).** Deaths (all detected by `e[0]` leaving
{0,2}, the only failure mode ever observed):

| R | D0 | bound `d_0-t_lo+1` | actual death row | margin |
|-----|----|------|------|-----|
| 256 | 0 | 52 | 61 | 9 |
| 512 | 0 | 102 | 107 | 5 |
| 512 | 1 | 104 | 110 | 6 |
| 512 | 2 | 106 | 113 | 7 |
| 512 | 3 | 108 | 116 | 8 |
| 512 | 4 | 110 | 121 | 11 |

In every death trace the gap `d_i - tmax(i)` closes at exactly 1/row after a
short dawdle: the front marches at the maximal speed permitted by T6a.
All deaths occur in the feed phase (`death row <= i_feed`).

## 8. Lemmas T7–T8 (endgame drain and ballot debt; PROVED + VERIFIED-finite)

**T8 (ballot debt).** `c_i(1) = c_{i,0} in {-2, 0}` always (box at t = 0),
so `s_i := e_i(1)` starts at `2 - R`, is even, and moves by `0` or `+2` per
row: the c0-word IS the scheduled ballot path of THM-3329 (E). Capture
(`e_i = 0`, hence `s_i = 0`) therefore requires `(R-2)/2` of the first
`i_c + 1` rows to be `-2`-rows:
```
                 capture row  i_c  >=  (R-2)/2 - 1 .
```
The plain rule's transient is provably Theta(R) rows long — harmless for
closure (only closure by row R-1 matters), decisive for understanding: the
mid-phase `c0 = -2` lock observed in every winning trace is forced
scheduling, not an artifact. Verified exactly: `#{-2 rows} = (R-2)/2` in all
three winning traces (63 at R=128; 255 at R=512, D0 = 5 and 8), all paid by
capture; capture rows 81 (0.63R), 312 (0.61R), 319 (0.62R).

**T7 (drain lemma).** If `e_{i-1} = v x^{d_i - 1}` with v even and
`-d_i <= v <= -2`, then (`delta = 0` case) `w = (v, v, 0, ...)`, or
(`delta = 1` case) `w = (v, 2v, v, 0, ...)`; cells >= 1 absorb, cell 0 emits
`c_0 = -2`, and `e_i = (v+2) x^{d_{i+1}-1}` of the same shape: the singleton
drains at exactly `+2` per row, riding the ballot cone boundary, and
vanishes in `|v|/2` rows. (Direct computation, both `delta_i` cases;
certificate `amm12592_transient_lemma_certificates_boxeph.out`. The
hypothesis `|v| <= d_i` is necessary: at a d-increment row cell 1's box
`[-2d, +2]` clamps `2v` when `v < -d` — caught by the certificate run, which
falsified the naive range `|v| <= 2(d-1)`.) This is precisely the observed last-phase of every winning trace
(final 12–24 rows) and the mechanism behind THM-3329's "1/row over the last
14 rows" LP observation. Closure demands the residual debt at the singleton
stage satisfy `|v|/2 <=` remaining rows — a *deadline* on Phase II.

## 9. The remaining analytic estimate (CONJECTURED, precisely stated)

The proved skeleton reduces the all-R transient theorem to:

**Estimate E (front stall + annihilation).** There exists `f(R) = o(R)`
(target: O(log R) or O(1)) such that for every dyadic R, at slack
`D0 = f(R)`, the clamped-Pascal flow of T6 with initial data T4 satisfies:

- **(E1, no-march)** for all rows `i <= i_feed`: for every cell
  `t' > F(i)`, `|(K * J_i)_{t'}| <= 2 C(d_{i+1}-1, t')` — i.e. junk spread
  beyond the current front is absorbed, so `F(i)` never reaches `d_i - 1`
  before the feed stops; and
- **(E2, annihilation deadline)** after `i_feed`, the junk L1 contracts to a
  cell-0 singleton with `|v| <= 2(R - i) - O(1)` in time for the T7 drain,
  i.e. capture occurs by row `R - 1 - |v|/2`.

**Then** (by T2+T3+T5+T6+T7) rule A closes epoch R at slack f(R), and since
o(R) per-epoch slack suffices for the assembly (THM-3329), we get
`C* <= 1 + gamma* = log_5(5 phi^2)` unconditionally.

**What the exact data says about E** (finite margins, no asymptotic claim):

- The dichotomy is razor sharp in D0: at R = 512, D0 <= 4 marches to death,
  D0 >= 5 stalls (front recedes monotonically after a <= 16-row advance) and
  closes. At R = 256 the flip is between 0 and 1; at R <= 128 D0 = 0 stalls.
- Local mechanism at the front (from T4's closed form): behind the crossover
  the load grows per cell (going down-cell) by
  `(1 - g tau*)/(g(1-tau*)) ≈ 3.05`, while caps grow by
  `tau*/(1-tau*) ≈ 2.05`; a marching front needs the crossover-cell excess
  ratio to exceed `~0.1` of cap — a Diophantine property of where the floor
  profile cuts the smooth crossover, which D0 shifts multiplicatively.
- **Known D0 thresholds (exact): D0(<=128) = 0, D0(256) = 1, D0(512) = 5.**
  Growth 0, 1, 5 across one doubling each: the data does not yet separate
  `f(R) = O(polylog)` from `f(R) = Omega(R^{1+eps})`; the running R = 1024
  ladder (control: D0 = 16/32/64) will discriminate.
- **Pre-registered predictions (R = 1024, probes launched this session):**
  exact T6b bounds `d_0 - t_lo + 1` are **202** at D0 = 0 (d = 612,
  t_lo = 411) and **218** at D0 = 8 (d = 620, t_lo = 403). If these slacks
  are below the 1024 stall threshold, deaths should land within ~5–15 rows
  above the bounds (margins were 5–11 at R = 256/512); a death BELOW either
  bound falsifies T6b; a stall at D0 = 8 (front receding by row ~300) would
  put D0(1024) <= 8 and kill the superlinear-growth scenario.

**Hazard notes.** (i) Rule-A deaths at small D0 are NOT evidence against
epoch feasibility at those profiles (search/rule negatives never are): the
ballot-scheduled repair (choose the c0-word globally instead of greedily)
is exactly what E1's failure localizes — the estimate pinpoints WHERE a
scheduled rule must deviate from greedy (the <= 11-row death margin says the
plain schedule loses by a bounded dawdle, not by capacity). (ii) The
`tau*` limit statement is first-order entropy comparison; all decisions in
this note use the exact integer window `t_lo(R, D0)`.

## 10. Corollary of the new witness (VERIFIED-exact, pending hostile audit)

Rule A closes R = 512 at `d_i = floor(gamma*(512+i)) + 5` (blocks
admissible, epoch identity `sum x^i Delta_i = q^{511}` exact; conjugate
implementation, block-equal to plain rule A at all verified scales). With
THM-3329's assembly: `T(n) <= n + 1 + floor(gamma* n) + 5` for all critical
n <= 1023 (slack constant improved from 8; D0 = 6, 7, 8 also close —
monotone in D0 on the tested range, so the D0scan gap 5..7 is now closed and
`D0(512) = 5` exactly, since D0 = 4 dies at row 121).

## 11. Status ledger

- PROVED: T1 coasting; T2 conjugacy + master-equation recovery; T3 parity
  inertness at dyadic R (variants coincide); T4 closed form `w_t`; T5
  peel/absorption/capture; T6 transport kernel + feed support + i_feed
  formula; T6a avalanche bounds; T6b death-delay `>= d_0 - t_lo + 1`;
  T7 drain; T8 ballot debt + capture `>= (R-2)/2 - 1`.
- VERIFIED-finite: conjugacy R = 8..512 (incl. stored-witness equality);
  parity zero-fires everywhere; T4 tables to R = 8192; death table + march
  speed; winning-trace anatomy (stall, capture rows, debt bookkeeping,
  singleton drain); D0(512) = 5 with D0 <= 4 deaths.
- CONJECTURED: Estimate E (E1 no-march + E2 annihilation deadline) with
  `f(R) = o(R)`; `tau*` limit identity; D0(R) growth law (open — decisive
  data point: R = 1024 ladder).

The all-R theorem = Estimate E. Everything else in the chain
`(*) -> C* <= log_5(5 phi^2)` is now proved or exactly verified.

---
---

# PART II — second sitting (same day): the fast T6 flow, the exact edge lemma, and the fate of Estimate E

New scripts: `amm12592_transient_fast_junkflow_boxeph.py` (T6 clamped-Pascal
flow as an O(R x band) engine + certification battery),
`amm12592_transient_T4b_edge_lemma_boxeph.py` (edge lemma exact grid
certificate). New outputs: `amm12592_transient_fastflow_cert_boxeph.json`,
`amm12592_fastflow_growthlaw_boxeph.json`,
`amm12592_fastflow_trace_R{1024,2048,4096,8192}_D0*_boxeph.json`,
`amm12592_transient_T4b_certificate_boxeph.json`, witness
`04-computation/amm12592_witness_R1024_ruleA_D0_15_fastflow_boxeph.json`.

## 12. The fast engine (T6 made executable; PROVED conjugate + certified)

Lemma T6 says the entire dynamics is the junk vector alone: w_i = K_delta *
j_{i-1} + feed_i (feed on cells {0,1}, values = original 2G_R coefficients),
j_i = overflow of w_i, death iff junk reaches cell d_i, capture iff j = {} after
feed stops. Implementing exactly this (sparse dict, incremental exact
binomials) gives a ~10^4 speedup over the polynomial implementation (R = 512:
0.1 s vs ~10^3 s) with zero change in semantics. Certification
(`..._fastflow_cert_boxeph.json`, all exact):

- FULL block-for-block equality vs the slow error dynamics at R = 8..128,
  plus independent admissibility + epoch-identity checks of the reconstructed
  witnesses: all True.
- Per-row trace equality (d, nclamped, tmin, tmax, junkL1_bits, c0, e_in0)
  against every STORED slow trace: R = 256 D0 = 0/1, R = 512 D0 = 0/4/5/8.
  Zero mismatches on all compared rows; deaths bit-identical (rows 61/107/121,
  same const sizes); the only differences are rows after capture, which the
  fast engine provably need not simulate (T1 coasting; stored traces confirm
  nclamped = 0 on every post-capture row).
- Evenness asserted at every clamp (T3 parity certificate, now to R = 8192).

## 13. Lemma T4b (exact edge law; PROVED skeleton + exact grid certificate)

**Statement.** For d in the window R/2 < d < 2R/3 (with the explicit boundary
inequality (ii)), the row-0 load w_t (T4 closed form) violates the c-box
exactly for t <= R-2-d and lies inside it for all t >= R-1-d. Hence:

- initial junk front F(0) = R-2-d_0 EXACTLY; t_lo = R-1-d_0 EXACTLY;
- T6b death-delay bound = d_0 - F(0) = **2 d_0 - R + 2 = (2 gamma* - 1) R
  + 2 D0 + O(1)** in closed form;
- the T4 "entropy-root" constants collapse: tau* = (1-gamma*)/gamma* and
  c1 = 2 gamma* - 1 = 0.19597487133...  (the entropy equation is satisfied
  IDENTICALLY at tau = (1-g)/g via H(p) = H(1-p); the numeric root of the
  first sitting was this point; "no nicer closed form" is superseded).

**Proof skeleton** (each step an exact integer inequality, certified on a
hostile grid incl. all dyadic R = 128..8192 at many D0 and non-dyadic R):
(I) for t <= R-3-d: A_t := C(R-2-t, d-t) >= C(d+1, d-t) = C(d+1, t+1) =: B_t
(top-index monotonicity + symmetry), so the even-sign cells satisfy
w_t >= 2C(d,t) > 2C(d-1,t-1) (overflow above) and the odd-sign cells
w_t <= -2C(d,t+1) < -2C(d-1,t) (overflow below).
(II) boundary t = R-2-d (sign +1 at even R): w = 2C(d,t) - C(d,t+1) and
overflow-above is EQUIVALENT to d(3R-4d-4) > 2(R-2-d)(R-1-d), a quadratic
with asymptotic roots d/R = 1/2, 2/3 — strictly true in the gamma* window.
(III) for t >= R-1-d: A_t <= C(d-1, t-1), C(d,t) >= C(d,t+1), and the two
in-box bounds reduce to (d-t)(d-t-1) >= 0 (upper) and d <= 3t+3 (lower),
both automatic in this region. QED (modulo writing the O(1) boundary cases,
all machine-checked).

This converges EXACTLY with the independent Angle-B1 session's "edge law"
edge(R, D0) = 2 d_0 - (R-2) (their capacity branch, verified on 30+ scales,
including non-dyadic parity branch 2^{v2(R)} via Lucas): two independent
derivations, one identity. T4b was their pre-registered "B2 target lemma";
it is now proved.

## 14. The growth law: D0(R) thresholds to R = 8192 (VERIFIED-exact)

All by the certified fast engine; every threshold pinned by adjacent
die/close pairs; every closure has minus2_rows = (R-2)/2 exactly (T8 debt
paid) and capture at ~0.62 R; every death obeys the (now closed-form) T6b
bound with margin 5..32 rows and marches at exactly 1 row/row after the
dawdle:

| R | D0* (exact) | D0*-1 dies at row | T6b bound | capture row at D0* |
|------|------|------|------|------|
| 128 | 0 | — | — | 81 |
| 256 | 1 | 61 | 52 | 159 |
| 512 | 5 | 121 | 110 | 312 |
| 1024 | 15 | 250 | 230 | 639 |
| 2048 | 38 | 508 | 476 | 1271 |
| 4096 | 89 | 1014 | 980 | 2537 |
| 8192 | 192 | 2045 | 1988 | 5064 |
| 16384 | in [401, 416] | 4055 (D0=400) | 4012 | 10126 (D0=416) |

Pre-registered first-sitting predictions came true: R = 1024 D0 = 0 died at
207 vs bound 202 (margin 5, inside the predicted 5-15); D0 = 8 died at 227
vs bound 218; no death below a T6b bound anywhere (falsifier absent).

**Two conjectures REFUTED by this table:**

1. Angle B1's polylog fit D0*(128·2^n) = C(n+3,4) (exact on 0,1,5,15;
   predicted D0(2048) = 35, D0(4096) = 70). Refuted: 35, 36, 37 all die at
   R = 2048 (rows 494/496/508); 70 and even 88 die at R = 4096 (rows
   962/1014). D0*(2048) = 38, D0*(4096) = 89.
2. First-sitting's open reading "data does not separate polylog from
   linear-with-log corrections" is superseded by:

**The measured law.** Normalized thresholds 1024·D0*/R (R = 256..8192):
**4, 10, 15, 19, 22.25, 24** — vs the pure arithmetic-decay integer pattern
4, 10, 15, 19, 22, 24 (increments 6, 5, 4, 3, 2): the ONLY deviation on six
doublings is +1 in D0 at R = 4096 (89 vs 88). The pattern's R = 8192
prediction (D0 = 192) was verified EXACTLY (191 dies at row 2045, 192
closes). Continuing the arithmetic decay (increments 1, then 0) predicted
**D0*(16384) = 400** and saturation D0*(R) = 25R/1024 for R >= 16384, i.e.
eps_inf = 25/1024 ~ 0.0244. RESULT: **D0 = 400 DIES at R = 16384** (row 4055
vs T6b bound 4012, margin 43 — and the T9' two-phase law predicted the death
row exactly ~700 rows in advance from the mid-march state). So the exact
saturation-at-25/1024 reading is ALSO refuted: **D0*(16384) in [401, 416]**
(400 dies verified; 416 CLOSES verified, capture 10126, debt 8191 exact;
D0 = 408 was mid-march when its run was reallocated — unverified). Increments
of the normalized threshold: 6, 5, 4, 3.25, 1.75, then in (1.06, 2.0].
Threshold ratios D0*(2R)/D0*(R): 5, 3, 2.53, 2.34, 2.16, ~2.1 — declining to
2 from above, the signature of LINEAR growth with a decaying correction.
Best reading: eps_inf >= 401/16384 ~ 0.0245, estimate 0.025-0.028.
CONJECTURE (revised): **the plain rule A needs LINEAR slack:
D0*(R)/R -> eps_inf ~ 25/1024**, i.e. Estimate E (o(R) for plain rule A) is
FALSE. This is a statement about rule A, not about epoch feasibility (hazard
discipline: rule negatives never prove infeasibility; the LP at R = 128
already beat greedy elsewhere). The eerily clean arithmetic decay of the
normalized increments (6,5,4,3,2 exactly, one +0.25 blip) is itself an
unexplained structure — it smells of a second exact edge law one level up.

**Scaling-limit reading (CONJECTURED, now precise).** With D0 = eps·R the
whole initial data + caps + feed become R-homogeneous: junk cells live at
relative positions t/R, loads and caps are exp(R·(entropy densities)), and
the march/stall dichotomy becomes a deterministic property of eps alone in
the R -> infinity limit — predicting D0*(R)/R -> eps_inf with slowly
convergent corrections, exactly as measured. Proving Estimate E was the
wrong target for the plain rule; the right unconditional statements are:

- (E-lin) for every eps > eps_inf there is R_0(eps) with rule A closing all
  dyadic R >= R_0 at D0 = eps·R — would give C* <= 1 + gamma* + eps
  unconditionally for every eps > eps_inf, i.e. **C* <= 1 + gamma* +
  eps_inf** (with the finite scales 128..8192 closed by table above);
- the golden target C* <= 1 + gamma* needs a BULK-modified rule (see 15).

## 15. Phase-transition anatomy and where the golden rule must differ

Exact trajectories at adjacent D0 (R = 2048: D0 = 37 dies, 38 closes):

- DEATH: gap(i) = d_i - front(i) dawdles ~16-24 rows near its initial value
  2d_0-R+2, then closes at EXACTLY 1/row to 0. Junk L1 first SHRINKS
  (~ -1.5 bits/row) to a minimum near half the death row, then regrows
  (+1.6 bits/row) as the marching spike amplifies; the bottom-cell overflow
  at death is astronomically large (2^557..2^4921 across R) — never a
  near-miss.
- CLOSURE at D0*: the gap NEVER goes below its row-0 value (min over the
  whole run is at i = 0 in every threshold closure trace); the front recedes
  from the start while junk L1 decays smoothly to O(1) bits by capture at
  ~0.62 R. The dichotomy is decided in rows ~20-70: whether the
  feed-repumped kernel spread at the frontier is absorbed or not.

Combined with the Angle-B1 steering-invariance lemma (modifying steered
cells k <= K changes later residuals only at absolute positions
>= i + d_i - K, while every die reads position ~ (2 gamma* - 1) R + 2 D0
<< d_0): the die is a BULK phenomenon of the cellwise clamp in the high
cells; no ballot/low-cell schedule can reach it. Any rule that beats
eps_inf must change the clamp/target on the high-index cells — e.g.
anticipatory (tent) shaping of the junk BAND (not the k <= 1 cells), or a
two-row lookahead clamp that solves the 2x2 kernel locally. That is the
sharp, testable design frontier the fast engine can now sweep at
R = 4096-16384 in seconds-to-minutes per run.

## 16. Updated records and status ledger (supersedes sec. 10-11 partially)

- New attainment records (deterministic rule A, exact, fast-engine verified;
  slow-path cross-checks running): R = 1024 closes at D0 = 15, R = 2048 at
  38, R = 4096 at 89, R = 8192 at <= 210. Via THM-3329 assembly the
  n <= 2047 attainment becomes T(n) <= n + 1 + floor(gamma* n) + 15.
  (Witness JSON for R = 1024 D0 = 15 saved; sparse-c format + ballot
  reconstruction documented in file.)
- PROVED (new): T4b exact edge lemma (initial junk block = cells
  [0, R-2-d_0]); closed forms t_lo = R-1-d_0, tau* = (1-gamma*)/gamma*,
  c1 = 2 gamma* - 1; T6b bound = 2d_0 - R + 2 exactly.
- VERIFIED-exact (new): D0* table to R = 8192; T8 debt = (R-2)/2 in every
  closure to R = 8192; phase-transition anatomy; refutation of the two
  growth conjectures above.
- CONJECTURED (revised): D0*(R)/R -> eps_inf in [0.024, 0.04] > 0 (plain
  rule A); (E-lin) closure at every eps > eps_inf; bulk-modified rule with
  o(R) slack exists (the golden route).
- The all-R golden theorem now = (bulk rule design) + (transient bound for
  that rule). The reduction chain, the engine, the exact initial data, the
  death mechanism, and the debt/drain endgame are all in place and proved.

## 17. Lemma T9 (march structure) — proved parts, falsified part, exact law

Script mode `T9` in the fast engine; output
`amm12592_transient_T9_march_check_boxeph.json`.

- **T9a (edge transport; PROVED + verified all runs).** The load at the new
  extreme cell t' = tmax + 1 + delta equals the old extreme junk value
  EXACTLY (only the edge kernel term reaches it).
- **T9b (phi condition; exact inequality, verified at every front of every
  run).** Caps along the march diagonal decay iff d(d-t-1) < (t+1)(t+2),
  asymptotically tau = t/d > 1/phi = 0.618. Since the front starts at
  tau_0 = (1-gamma*)/gamma* = 0.672 > 1/phi and only deepens, cap decay
  holds along every march — the front, once genuinely marching, sheds cap
  protection geometrically. (Golden ratio appears AGAIN, as the stability
  threshold of the march.)
- **T9c (geometric-tail trigger; FALSIFIED as a sufficient condition).** My
  first formalization ("extreme junk > sum of remaining caps implies
  permanent march") fires at row 1 on the R = 256..1024 deaths while the
  front still dawdles. Reason: the OVERFLOW tapers to small values at the
  edge cells (T4b boundary), so the front's own junk is absorbable at its
  landing cell; the march is driven by the INTERIOR bulge pushing from
  behind, not by edge transport. Kept on record per hazard discipline.
- **T9' (two-phase death law; VERIFIED-exact at 6 scales incl. all
  threshold-minus-1 deaths).** Every death trajectory is: a dawdle of
  L rows (L = 5..27, weakly growing with R and D0), then a PURE march —
  the gap d_i - front(i) decreases by exactly 1 EVERY row, zero recoveries,
  and  **death row = L + gap(L)** exactly:
  (R, D0, L, gap(L), death) = (256,0,11,50,61), (512,0,5,102,107),
  (512,4,13,108,121), (1024,0,5,202,207), (1024,14,20,230,250),
  (2048,37,27,481,508). With T4b's exact initial gap 2d_0 - R + 2 this
  pins death rows to [bound, bound + L + few] and reduces the whole
  march/stall question to the dawdle-phase decision.

Cross-session confirmation: the Angle-B1 independent referee (sha-pinned,
from-the-theorems checker) has verified ALL-PASS witnesses for R = 512
D0 = 5 (rule A AND rule B files) and R = 1024 D0 = 15; B1's slow runner
independently found the same D0 thresholds 0/1/5/15 at R <= 1024 and the
same die rows 227..250 at R = 1024 D0 = 8..14 — bit-level agreement with
the fast engine across three implementations and two sessions. The slow
error-dynamics cross-check at R = 1024 D0 = 14 is COMPLETE: death at row 250
with const_bits 599, bit-identical to the fast engine, all 251 rows matching
on every compared stat (d, nclamped, tmin, tmax, junkL1_bits, c0, e_in0).

Design note for the bulk rule (the golden route): the initial junk profile
ALTERNATES in sign cell-to-cell (T4's (-1)^{d-t}), and the kernel (1,2,1)
annihilates alternating profiles — this is why the stall decays the junk at
~1.5 bits/row at D0*. The march-driving bulge is the non-alternating
residue. A bulk clamp that SHAPES junk toward alternation (choosing c below
the edge on selected cells to restore sign alternation) attacks the
mechanism directly and is sweepable with this engine at R = 4096+ in
minutes. That, not ballot steering (refuted by B1's invariance lemma), is
the concrete candidate for killing eps_inf.

## 18. Correction to Section 7 (feed-phase death claim) + refined criterion

Sec. 7's "all deaths occur in the feed phase" was an artifact of the small
scales then available (R <= 512 near-threshold deaths, small-D0 R = 1024).
Exact i_feed vs near-threshold death rows (fast engine):

| R | D0 | i_feed | death row | excess past feed-end |
|------|-----|------|------|-----|
| 2048 | 37 | 491 | 508 | 17 |
| 4096 | 88 | 974 | 1014 | 40 |
| 8192 | 191 | 1940 | 2045 | 105 |
| 16384 | 400 | 3870 | 4055 | 185 |

Near threshold, the march OUTLIVES the feed (the post-feed flow is
autonomous and, by T9' as measured, never de-marches). Refined criterion:
the die/stall decision is whether the interior bulge triggers the front
march AT ALL during the feed phase; once triggered the march runs to death
regardless of feed-end (excess roughly doubling per scale). Deaths at small
D0 (far below threshold) do still land inside the feed phase.
