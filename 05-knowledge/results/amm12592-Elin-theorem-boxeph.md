# AMM 12592 — Lane D2: the E-lin theorem (linear slack), proved skeleton + one sharp hypothesis

Session: boxeph multifront, 2026-08-03 (post transient-theorem Parts I+II).
All computations exact (int / Fraction; no floats in any decision; sympy not
used — gamma* is certified by pure-integer Lucas/Fibonacci comparisons).

Scripts (04-computation/):
`amm12592_Elin_linear_slack_sweep_boxeph.py` (exact closures at D0 = ceil(eps R)),
`amm12592_Elin_epsstar_certificates_boxeph.py` (certificates C1–C6),
`amm12592_Elin_feedend_state_certificate_boxeph.py` (feed-end state N/C/W/F/D/DL),
`amm12592_Elin_postfeed_signstructure_boxeph.py` (exact mechanism instrumentation),
`amm12592_Elin_witness_referee_C7_boxeph.py` (independent-referee witness check).
Outputs (05-knowledge/results/): `amm12592_Elin_sweep_boxeph.json`,
`amm12592_Elin_epsstar_certificates_boxeph.json` (+ `_c15` stage file),
`amm12592_Elin_feedend_state_boxeph.json`,
`amm12592_Elin_postfeed_signstructure_boxeph.json`,
`amm12592_Elin_witnessC7_boxeph.json`,
`amm12592_fastflow_trace_R32768_D02048_boxeph.json` (R = 32768 probe; see 4.1).

Throughout: `g = gamma* = log_5(phi^2)`, epoch problem (*) at slack D0 is
`sum_{i<R} x^i Delta_i = (1-x)^{R-1}` with `Delta_i` admissible at
`d_i = floor(g(R+i)) + D0`; rule A and its junk-flow conjugates are as in the
transient theorem (T1–T9', all citations to that note).

---

## 0. Headline

1. **LIFT (Theorem A, PROVED, 4 lines).** Admissibility is preserved verbatim
   under raising the degree profile: a witness at slack D0 is a witness at
   every slack D0' >= D0. Consequence: for the C* assembly at slack eps R it
   suffices that (*) close at SOME D0 <= eps R. E-lin is thereby weakened
   from "rule A closes at exactly ceil(eps R)" to "epoch feasible at
   <= eps R" — and any provable closing construction may be substituted for
   plain rule A.
2. **Feed-phase survival (Theorem B, PROVED).** Closed form
   `i_feed = floor((R(1-g) - D0)/(1+g))`, and: for
   `D0 >= eps* R`, `eps* := 2(1 - g - g^2)/(3 + 2g) = 0.0211741...`,
   within the T4b window `R/2 < d_0 < 2R/3`, the junk flow of rule A cannot
   die at any row `i <= i_feed`. The T6b death bound `2 d_0 - R + 2` clears
   the feed horizon exactly when `D0 >= eps* R`; the threshold eps* is a new
   closed-form constant of the problem, strictly BELOW the measured
   `eps_inf ~ 0.025` — consistent with (and explaining) the observed fact
   that near-threshold deaths always occur after feed-end.
3. **Post-feed structure (Theorem C, PROVED one-step lemmas).** In the
   autonomous phase: negativity of the junk is invariant; every overflowing
   cell's cap-ratio `A_t = |j_t| / 2C(d-1,t)` obeys the exact damping
   `A'_t <= A_t (d-t)/d + 2 A_{t-1} t/d + A_{t-2} t(t-1)/(d(d-t+1)) - 1`
   (delta = 1 rows; analogue for delta = 0); an explicit two-cell condition
   freezes the front; and the T8 debt IS the cell-0 junk:
   `e_i(1) = j_i[0]`, so the endgame is a forced 2/row drain with an exact
   deadline. These pin the shape of the one remaining estimate.
4. **The sweep (VERIFIED-exact).** Rule A CLOSES every dyadic epoch
   `R = 128 .. 16384` at `D0 = ceil(R/32)` and `D0 = ceil(R/16)` (16/16
   closures, debt `(R-2)/2` paid exactly in every run), and also at
   eps = 1/8, 1/4, ~0.4 probes (R = 512..2048) — no non-monotonicity of rule
   A in D0 was detected anywhere.
5. **Conditional theorem (E-lin).** Modulo one precisely-stated hypothesis
   S(eps) about the autonomous phase (Section 5) — exactly verified for all
   dyadic `128 <= R <= 16384`, open for `R >= 32768` —
   `C* <= 1 + gamma* + 1/16 < 435287/262144 = 1.6604882`, and with S(1/32),
   `C* <= 1 + gamma* + 1/32 < 427095/262144 = 1.6292382`.
   Unconditionally and independent of S, the sweep gives the attainment
   `T(n) <= n + 1 + floor(gamma* n) + ceil(n/16) + O(1)` on the verified
   range (weaker than the constant-slack records for small n; recorded for
   completeness of the assembly's finite part).

Status labels: **PROVED** (complete argument), **VERIFIED-exact** (exact
computation at stated scales), **HYPOTHESIS** (precise open statement).

---

## 1. Theorem A (LIFT: profile monotonicity of feasibility) — PROVED

**Statement.** Let `Delta` be admissible at degree d (cells
`delta_t`, `|delta_t| <= C(d,t)`, `delta_t == C(d,t) mod 2`). Then the same
polynomial is admissible at degree d+1; hence if (*) is solvable at profile
`{d_i}` it is solvable (by the SAME blocks) at any profile `{d'_i}` with
`d'_i >= d_i`, in particular at any larger slack.

*Proof.* Multiply each Bernstein term by `1 = x + (1-x)`:
`Delta = sum_t delta_t x^{d-t} q^t = sum_t (delta_t + delta_{t-1}) x^{d+1-t} q^t`,
so the degree-(d+1) cells are `delta'_t = delta_t + delta_{t-1}`. Then
`|delta'_t| <= C(d,t) + C(d,t-1) = C(d+1,t)` and
`delta'_t == C(d,t) + C(d,t-1) = C(d+1,t) (mod 2)`. The polynomial — hence
the epoch identity — is unchanged. Iterate rowwise. QED

Certificate C6 (`amm12592_Elin_lift_C6_standalone_boxeph.json`): the
R = 128 rule-A witnesses at D0 = 4 and 8 lifted by +1 and by +7 re-verify
admissibility and the epoch identity exactly (a second instance at
R = 512, D0 = 5 is in the full battery JSON when its C6 stage lands).

**Consequence for E-lin.** The C*-relevant statement is epoch FEASIBILITY at
slack `<= eps R`, not the behavior of rule A at exactly `ceil(eps R)`:
closure at any smaller slack lifts. (The sweep found no D0-non-monotonicity
of rule A anyway, but nothing below depends on rule-A monotonicity.)

## 2. Theorem B (exact i_feed + feed-phase survival) — PROVED

**(B1) i_feed closed form.** `i_feed := max{ i : d_i + i <= R-1 }` (the last
fed row; feed values are the original 2G_R coefficients, T6-feed). Since
`d_i + i` is strictly increasing and, for integer RHS,
`floor(g(R+i)) <= R-1-D0-i  <=>  g(R+i) < R-D0-i  <=>  i(1+g) < R(1-g)-D0`:

```
i_feed = floor( (R(1-g) - D0) / (1+g) ).
```

`x = (R(1-g)-D0)/(1+g)` is never an integer: x rational forces
`g = (R-D0-x)/(R+x)` rational, but g is irrational (else `phi^{2q} = 5^p`,
impossible as `phi^{2q} = (L_{2q} + F_{2q} sqrt5)/2` with `F_{2q} != 0`).
Certificate C3: formula == brute-force scan on a 49-point hostile grid
(all dyadic R = 128..16384 at five D0 each + nine non-dyadic R), with
near-integer ceilings decided by the exact Beatty test.

**(B2) Survival theorem.** Suppose
`(i) D0 >= eps* R`, where `eps* := 2(1 - g - g^2)/(3 + 2g)`, and
`(ii) R/2 < d_0 < 2R/3` (T4b window; d_0 = floor(gR) + D0).
Then the junk flow survives every row `i <= i_feed`: death cannot occur
before the feed stops.

*Proof.* Row 0 is safe outright: by T4's closed form the bottom-cell load is
`w_{d_0} = C(R-2-d_0, 0) - C(d_0+1, d_0+1) + 2C(d_0,d_0) = 1 - 1 + 2 = 2`
(valid for any `d_0 <= R-2`), inside `{0,2}`. By T4b (edge lemma, proved;
re-anchored at the exact D0 range used here by certificate C5: initial junk
support is EXACTLY `[0, R-2-d_0]` at every swept (R, D0)), the initial front
is `F(0) = R-2-d_0`. By T6a the
front advances at most `1 + delta_i` per row while the bottom recedes
`delta_i`, so death (junk at cell `d_i`) requires
`i >= d_0 - F(0) = 2 d_0 - R + 2` (T6b, exact). It remains to show
`2 d_0 - R + 2 > i_feed`. Using `d_0 >= gR - 1 + D0` and (B1),
it suffices that `(2g-1)R + 2 D0 >= (R(1-g) - D0)/(1+g)`, i.e.

```
D0 (3 + 2g)  >=  R (1 - g)  -  (1+g)(2g-1) R  =  2 R (1 - g - g^2),
```

which is exactly `D0 >= eps* R`. QED

**(B3) Exact constants (certificates C1, C2).** By integer bisection on
`5^a <=> phi^{2M}` via `L_{2M} + F_{2M} sqrt5` (M = 2^20):

```
g   in ( 627035/1048576 , 156759/262144 )      (width 2^-20)
eps* in ( 0.021173651... , 0.021174659... )     (exact Fractions in JSON)
eps* < 1/32 < 1/16 ,
```

and eps*(g) is decreasing in g (derivative numerator certified negative), so
the upper end comes from g_lo. Window (ii) holds for `D0 = ceil(R/32)` once
`R >= 27` and for `D0 = ceil(R/16)` once `R >= 162` (from
`gR + eps R + 1 < 2R/3` with the sandwich; dyadic cases below these bounds
are in the verified table anyway). Certificate C4 instantiates
`2d_0 - R + 2 > i_feed` on the full grid.

**Interpretation.** eps* is a new closed-form constant of the epoch problem:
the slack at which the proved death-delay bound clears the feed horizon.
The three transition constants order as

```
eps_phi = 1/phi - g = 0.020047...  <  eps* = 0.0211741...  <  eps_inf(measured) ~ 0.025 ,
```

where eps_phi is the slack at which the initial front sits exactly on the
T9b golden stability diagonal `t/d = 1/phi` (`(1-g-eps)/(g+eps) = 1/phi`
solves to `eps = 1/phi - g` since `1/phi + 1/phi^2 = 1`); recorded as
structure, not used in proofs.

## 3. Theorem C (post-feed structure lemmas) — PROVED (one-step)

Post-feed (`i > i_feed`) the flow is autonomous: `w = K_delta * j`, clamp,
no feed. Write `d' = d + delta` for the new degree, `m = max support(j)`,
caps `2C(d'-1, t)` (lower) / `2C(d'-1, t-1)` (upper).

**(C-N) Negativity is invariant.** If `j <= 0` cellwise then `j' <= 0`.
*Proof.* `w = K*j <= 0`; upper cap >= 0, so the clamp is
`u_t = max(-2C(d'-1,t), w_t) <= 0` and
`j'_t = w_t - u_t = min(0, w_t + 2C(d'-1,t)) <= 0`. QED

**(C-A) Cap-ratio damping.** Let `A_t := |j_t| / (2C(d-1,t))` (t >= 1,
all-negative junk). If cell t overflows at the next row, then exactly:

```
delta = 1 :  A'_t <= A_t (d-t)/d  +  2 A_{t-1} t/d  +  A_{t-2} t(t-1)/(d(d-t+1))  -  1
delta = 0 :  A'_t <= A_t  +  A_{t-1} t/(d-t)  -  1
```

(if cell t does not overflow, `A'_t = 0`).
*Proof.* `|j'_t| <= |w_t| - 2C(d'-1,t)` at an overflowing cell, and
`|w_t| <= |j_t| + 2|j_{t-1}| + |j_{t-2}|` (delta = 1). Divide by
`2C(d,t)` and use the exact ratios `C(d-1,t)/C(d,t) = (d-t)/d`,
`C(d-1,t-1)/C(d,t) = t/d`, and
`C(d-1,t-2)/C(d,t) = [(t-1)/(d-t+1)] * [t/d]` (product of the adjacent-cell
ratios); delta = 0 analogous with `C(d-1,t-1)/C(d-1,t) = t/(d-t)`. QED
For `t <= m << d` the neighbor coupling is `O(m/d)` and every live cell
loses one full cap unit per row: this is the exact sense in which the caps
DRAIN the residue additively, and the measured extinction times
(`~ max_t A_t` rows, see Section 4) are its prediction.

**(C-F) Front freeze.** Suppose `j <= 0`, support ⊆ `[0, m]`, and
(delta = 1 case) `|j_m| <= 2C(d, m+2)` and
`2|j_m| + |j_{m-1}| <= 2C(d, m+1)`; (delta = 0 case)
`|j_m| <= 2C(d-1, m+1)`. Then support(j') ⊆ `[0, m]`: the front cannot
advance. *Proof.* The only loads beyond cell m are
`w_{m+1} = 2 j_m + j_{m-1}`, `w_{m+2} = j_m` (delta = 1) or
`w_{m+1} = j_m` (delta = 0); under the stated bounds each is inside its
(lower) cap and is absorbed entirely. QED
Note the hypotheses are LOCAL (front cells only); they are re-established
at every row of every verified run (the profile is steeply increasing
toward low cells, so front values stay near cap scale).

**(C-D) Debt = cell-0 junk; drain deadline.** Post-feed the pristine tail
is empty, so `e_i = J_i / x` with `J_i = sum_t j_t x^{d-t} q^t`, and at
x = 1 only the t = 0 cell survives:

```
s_i := e_i(1) = j_i[0]  =  (2 - R) + 2 * #{ c0 = -2 rows so far }   (T8).
```

Since `c0 in {-2, 0}` always, capture (`j = {}`, so `s = 0`) requires at
least `|j_i[0]|/2` further rows from any post-feed row i: closure demands

```
|j_0(feed-end)|  <=  2 (R - 2 - i_feed)        (deadline DL),
```

and because `i_feed < R(1-g)/(1+g) < (R-2)/2`, the deadline is satisfiable
even if NO debt is paid during the feed phase — at linear slack the T8
deadline is never binding a priori; what must be proved is that the drain
actually runs (cell 0 stays saturated), which is part of Hypothesis S.

## 4. The exact sweep and the feed-end state — VERIFIED-exact

### 4.1 Closures at linear slack (`amm12592_Elin_sweep_boxeph.json`)

Every run: outcome CLOSED, debt `minus2 = (R-2)/2` exactly, capture row
`~ 0.61 R` (drifting to `(R-2)/2 - 1` — the T8 bound — as eps grows), and at
eps <= 1/16 the initial front equals `R-2-d_0` exactly:

| R | D0=ceil(R/32) capture | D0=ceil(R/16) capture | larger-eps probes |
|------|------|------|------|
| 128 | 79 | 78 | — |
| 256 | 153 | 150 | — |
| 512 | 317 | 313 | 1/8: 305; 1/4: 282; 0.398: 254 |
| 1024 | 624 | 616 | 1/4: 554; 0.399: 510 (= T8 bound) |
| 2048 | 1261 | 1237 | 1/8: 1205; 1/4: 1120 |
| 4096 | 2519 | 2486 | — |
| 8192 | 5040 | 4964 | — |
| 16384 | 10090 | 9937 | — |

16/16 primary closures; 7/7 probes. No death, no OPEN_RESIDUAL, no
non-monotone behavior anywhere in the sweep.

**Initial-front regimes (probe data, exact).** The probes confirm the
three-regime structure of the row-0 junk front and hence the exact validity
range of the T4b-based Theorem B:
`d_0 in (R/2, 2R/3)`: `F(0) = R-2-d_0` (T4b; all eps <= 1/16 runs);
`d_0 in (2R/3, ~3R/4)`: `F(0) = R-3-d_0` exactly (boundary cell flips
in-box: R = 512, D0 = 64 gives 139 = R-3-d_0; R = 2048, D0 = 256 gives
565 = R-3-d_0) — Theorem B would IMPROVE by 1 here if T4b's step (III)
were extended past 2R/3;
`d_0 > ~3R/4`: a second (ballot-type) overflow block appears above
`R-2-d_0` (front ~ d_0/3.4 at the d_0 ~ R probes), so the T4b description
genuinely ends — the window hypothesis (ii) is not an artifact.

**R = 32768 probe (launched this session).** Rule A at D0 = 2048
(eps = 1/16) is running on the certified engine at session close; its trace
(partial flushes every 128 rows) is at
`amm12592_fastflow_trace_R32768_D02048_boxeph.json`. If it CLOSES, the
verified range of Hypothesis S moves up one doubling (open only for
R >= 65536); a death would FALSIFY S(1/16) as stated and sharpen the
frontier — either outcome is informative. No claim is made pending the
exact outcome.

**Hostile-referee check (C7, `amm12592_Elin_witnessC7_boxeph.json`).** Full
witnesses reconstructed at (R, D0) = (256, 16) and (512, 32) and verified by
the INDEPENDENT referee implementation (its own admissible(), polynomial
arithmetic, Lucas/Fibonacci floor engine, and target `q^{R-1}`): block
admissibility at the `floor(g(R+i)) + D0` profile, epoch identity, and
forced units all PASS — the linear-slack closures are genuine witnesses,
not engine artifacts.

### 4.2 Feed-end state (`amm12592_Elin_feedend_state_boxeph.json`)

At the first autonomous row of every swept (R, eps) run, R = 128..16384,
eps in {1/32, 1/16} (16 records), certified exactly — ALL PASS:
**(N)** junk all-negative; **(C)** support a contiguous block `[0, m]`;
**(W)** `6m <= d`; **(F)** the C-F freeze inequalities; **(D)** the debt
identity `j_0 = (2-R) + 2 * minus2count`; **(DL)** the drain deadline.

| R | D0 | i_feed | m | j_0 | L1 bits |
|------|-----|------|-----|--------|-----|
| 128 | 8 | 27 | 7 | -100 | 35 |
| 256 | 16 | 54 | 12 | -190 | 63 |
| 512 | 32 | 109 | 18 | -406 | 104 |
| 1024 | 64 | 217 | 28 | -796 | 169 |
| 2048 | 128 | 435 | 39 | -1602 | 263 |
| 4096 | 256 | 870 | 58 | -3230 | 417 |
| 8192 | 512 | 1740 | 83 | -6446 | 636 |
| 16384 | 1024 | 3481 | 119 | -12910 | 972 |

(eps = 1/32 rows analogous; JSON has all 16.) Three OBSERVED scaling laws
of the feed-end state (exact data, no proof claimed):

- `m ~ sqrt(R)` (m / sqrt(R) = 0.62, 0.75, 0.80, 0.875, 0.86, 0.91, 0.92,
  0.93 — slowly increasing toward ~1): the autonomous problem starts on a
  root-scale, not linear-scale, band.
- `max_t log2 A_t ~ 10-11`, essentially R-INDEPENDENT across 7 doublings:
  the cap-excess of the residue does not grow with scale, while the row
  budget `R - 2 - i_feed ~ 0.79 R` grows linearly — strong structural
  support for Hypothesis S.
- `|j_0(feed-end)| = d(feed-end) + c` with `|c| <= 11` at ALL 16 records
  (e.g. 12910 vs d = 12903 at R = 16384; 6446 vs 6450 at R = 8192): the
  feed leaves the cell-0 debt exactly at the T7 drain-hypothesis edge
  `|v| <= d`. Unexplained; smells of another exact edge law (compare the
  T4b boundary).

Measured anatomy (R = 2048, D0 = 128; instrumentation JSON):

- Feed phase: junk profile is EXACTLY alternating (alt-fraction 1.0) for
  ~300 rows while L1 collapses 1887 -> 361 bits (~3.7 bits/row) and the
  gap grows from row 0; by feed-end (i = 436) the junk is an all-negative
  block on cells [0, 39], L1 = 263 bits.
  In the D0 = 37 DYING run at R = 2048 the alternation fraction decays
  from row ~32 (0.99 -> 0.19) and L1 turns around and grows — alternation
  integrity through the feed phase IS the observable that separates the
  phases (K4's annihilation mechanism, now measured at the transition).
- Post-feed: per-cell ratios `A_t` range from ~2^10 (cell 1) down to ~2^0
  (front); the profile rides just above the caps and is FROZEN in value
  while d grows; cells die from the top down as caps overtake them
  (front 39 -> 0 over ~460 rows); after cells >= 1 are extinct, the pure
  T7/T8 cell-0 drain closes at exactly 2/row: capture at
  1237 = 436 + |j_0(436)|/2 with j_0(436) = -1602 exactly.

## 5. Hypothesis S and the conditional theorem

**Hypothesis S(eps) (autonomous collapse; HYPOTHESIS).** For every dyadic
`R >= 32768`, at `D0 = ceil(eps R)`: starting from the feed-end state
(alive by Theorem B), the autonomous flow reaches `j = {}` by row `R - 2`
without junk ever reaching cell `d_i`.

Equivalently (given Theorems A/B and T1/T5 coasting): rule A closes epoch R
at slack `ceil(eps R)`. S(eps) is exactly verified for all dyadic
`R in [128, 16384]` at eps = 1/32 and 1/16 (Section 4.1) and is the ONLY
unproved step below. Supporting structure: the feed-end certificates
(N/C/W/F/D/DL), the C-A additive drain, the linearly growing row budget
(`R - 2 - i_feed ~ 0.79 R` available vs measured post-feed lifetime
`~ 0.40 R`), and the C-D deadline slack.

**Theorem D (E-lin, conditional on S).** Assume S(eps) for a rational
`eps in [1/32, 2/3 - g)`. Then every dyadic epoch `R >= 128` is feasible at
slack `ceil(eps R)` (R <= 16384 by the sweep; R >= 32768 by S), hence by
LIFT at every slack in `[ceil(eps R), ...]` as the assembly requires, and by
the THM-3329 assembly (finite epochs from the witness table; per-epoch slack
eps R):

```
C*  <=  1 + gamma* + eps .
```

In particular (rational upper bounds via the certified sandwich
`g < 156759/262144`): under S(1/16),
`C* <= 1 + g + 1/16 < 435287/262144 = 1.6604882`; under S(1/32),
`C* <= 1 + g + 1/32 < 427095/262144 = 1.6292382`.

**What a full proof of S needs (the envelope program, sharpened).** The
feed-phase collapse is the sole hard core: prove that at `D0 >= eps R` the
alternating component of the junk stays dominant through the feed phase
(measured: alt-fraction 1.0 for the first ~0.15R rows, decaying only after
L1 has collapsed), e.g. via a signed envelope
`j_t = (-1)^{d-t} a_t + r_t` with `a_t` log-concave and
`|r|` cap-dominated; the kernel is a second difference on `a`:
`(K * j)_t = ±(a_t - 2a_{t-1} + a_{t-2}) = ± a_{t-1}(rho_t - 2 + 1/rho_{t-1})`
with down-cell ratios `rho_t = a_t/a_{t-1}`, single-signed while `a` is
log-concave with `rho > 1` (which the T4 initial data satisfies on the
block: `A_{t-1}/A_t = (R-1-t)/(d-t+1)`), giving per-row magnitude factor
`~ (1 - 1/rho)^2` — e.g. `rho ~ 1.8` at the front predicts
`~ 2.4 bits/row` decay vs 3.7 measured (right order; the clamp corrections
add decay). Theorem B guarantees no death while this burns down. Post-feed, C-N/C-A/C-F
reduce extinction to a per-cell ratio-envelope
`A_t <= 2^{c (m0 - t)}`-type invariant (measured c ~ 6.5 bits/cell at
feed-end) whose one-step propagation is already proved in C-A; the missing
piece there is only the ordered (top-down) extinction bookkeeping.
This is a genuinely finite program: both envelopes are one-parameter
families of exact binomial inequalities of C-A type.

## 6. Hazards honored

- Rule/search negatives are never infeasibility evidence; none are used.
- All decisions exact; the only analytic constant (g) enters through a
  certified rational sandwich, never floats.
- S is labeled HYPOTHESIS despite 16/16 exact confirmations and mechanism
  understanding; the conditional theorem is stated as conditional.
- The engine's T4-closed-form initializer requires `d_0 <= R-2`; all runs
  respect it (asserted per run), so no large-eps extrapolation is claimed
  beyond the probed `eps <= 0.4` (where d_0 = R-2 exactly at R = 512).
- Non-dyadic R are not claimed anywhere (epochs are dyadic; T3 parity
  freeness holds at dyadic R and evenness is asserted at every clamp).
- Engine change this session: `initial_junk` now advances `C(d,t)`
  incrementally (performance only; O(d^2) -> O(d)); certified bit-identical
  to the direct closed-form clamp at four scales and by re-reproducing the
  R = 256, D0 = 16 closure (capture row 150) before any use at R = 32768.

## 7. Status ledger

- PROVED (new): Theorem A (LIFT); Theorem B (i_feed closed form; survival
  through feed for D0 >= eps* R in the T4b window; eps* = 2(1-g-g^2)/(3+2g)
  with certified rational sandwich); Theorem C one-step lemmas
  (C-N negativity, C-A cap-ratio damping, C-F front freeze, C-D debt
  identity + deadline); irrationality bookkeeping for i_feed.
- VERIFIED-exact (new): 16/16 closures at eps = 1/32, 1/16 for dyadic
  R = 128..16384 (+7 larger-eps probes); C1–C6 certificate battery;
  C7 independent-referee witness verification at (256, 16) and (512, 32);
  feed-end N/C/W/F/D/DL certificates (14/14 ALL PASS);
  the alternation-integrity dichotomy at the D0 = 37/38 transition
  (R = 2048); the frozen-profile/top-down extinction anatomy of the
  autonomous phase; the |j_0| = d +- O(10) feed-end edge law (observed).
- HYPOTHESIS (one item): S(eps) for dyadic R >= 32768 — the autonomous
  collapse. E-lin and `C* <= 1 + gamma* + 1/16` are conditional on it;
  every other link in the chain is proved or exactly verified.
- Companion targets fed by this note: the D1 lane (bulk alternation-shaping
  rule) inherits the exact envelope targets of Section 5; a proof of the
  feed-phase signed envelope would immediately convert Theorem D into an
  unconditional `C* <= 1 + gamma* + 1/32`.
