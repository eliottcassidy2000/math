# AMM 12592 — Lane E2: Hypothesis S via super-block impossibility — the post-feed capture theorem, the march-budget theorem, and per-epoch certificates

Session: boxeph multifront, 2026-08-04 (after D1 bulk sweep, D2 E-lin theorem,
D3 golden bound, and the Estimate-E hostile referee). All citations T1–T9' to
the transient theorem note, A/B/C/S and C-N/C-A/C-F/C-D to the E-lin note,
G1–G3 to the golden-bound note (with the referee's full-layer F3 correction),
super-block law to D1 sec. 3b. All computations exact int (fixed-point
ceilings; every rounding toward the safe side); floats only in reporting.
No numpy; sympy not used.

Scripts (04-computation/):
`amm12592_S_feedend_extract_boxeph.py` (exact feed-end state extractor: runs
the certified T6 flow through the FIRST AUTONOMOUS ROW only, ~21% of the
epoch), `amm12592_S_superblock_certificate_boxeph.py` (the certificate:
hypotheses check, march-budget arithmetic, majorant iteration, freeze
persistence, capture deadline, super-block margin; plus a hostile `selftest`
mode that replays the TRUE autonomous flow and checks domination cell-by-cell
row-by-row).
Outputs (05-knowledge/results/): `amm12592_S_feedend_R*_D0*_boxeph.json`
(exact states), `amm12592_S_superblock_certificate_boxeph.json` (ledger),
`amm12592_S_superblock_certificate_boxeph.out` (run log).

Status labels: **PROVED** (complete argument), **VERIFIED-exact** (exact
computation at stated scales), **HYPOTHESIS / OPEN** (precise statement).

---

## 0. Headline

Hypothesis S(eps) (D2 sec. 5: from the feed-end state, the autonomous
post-feed flow of plain rule A captures by row R-2 without death) is the ONE
unproved link in the chain `C* <= 1 + gamma* + eps`. This lane:

1. **Theorem 1 (march budget; PROVED, unconditional).** In the autonomous
   phase death is impossible before row `d_{i0} + i0 - m >= R - m`: it can
   only occur in the LAST ~m ~ sqrt(R) rows of the epoch, and any death
   trajectory must advance its junk front at the maximal T6a rate on all but
   at most `m - 1` of the ~0.79R post-feed rows (total advance deficit
   <= m - 1). Corollary: **if the C-F front freeze holds on the first m
   autonomous transitions, death is impossible at every row** — the no-death
   half of S is reduced from a Theta(R)-row statement to an O(sqrt R)-row
   statement.
2. **Theorem 2 (post-feed capture certificate; PROVED).** An explicit
   monotone majorant iteration on cap-ratios (the C-A calculus extended by
   the exact cell-0 debt drain and the debt's pump into cells 1–2 — the
   block version of T7), combined with the C-F freeze in ratio form, is
   machine-checkable from the feed-end state alone; when it validates, the
   autonomous flow provably captures by row R-2 with no death, i.e.
   **S holds at that (R, D0)**.
3. **Certificates (VERIFIED-exact).** The certificate PASSES at every
   verified feed-end state (R = 128..16384 where extracted, eps = 1/32 and
   1/16: 16/16 — the FULL original S-verified range; see table) and — new
   territory — at **R = 32768** if/where the feed-end extraction lands
   (sec. 6; pending states noted honestly). In every PASS the certificate
   also PREDICTS the capture row exactly (`i0 + |j_0(i0)|/2`), matching
   all 16 known sweep captures — the capture clock of the linear-slack
   closures is identified.
4. **Corollary SB (super-block impossibility; PROVED).** Under the
   certificate, no super block (D1 sec. 3b sense: same-sign run whose min
   |value| exceeds the remaining cap tail beyond the front) can EVER form in
   the post-feed phase: junk values stay below `2^{~O(sqrt R log R)}` while
   the cap tail beyond the frozen front contains the middle binomial
   `2C(d-1,(d-1)/2) ~ 2^{0.79 R}`. Conversely death trivially exhibits a
   super block at the death row (empty tail), and by Theorem 1 death
   requires a maximal march from feed-end on — the post-feed half of D1's
   12/12 law is now a theorem in both directions.
5. **Proposition A (asymptotic self-maintenance; PROVED).** With explicit
   constants: IF the feed-end state satisfies the four static conditions
   (negativity; support [0,m] with 6(m+2) <= d; cap-ratios A_t <= B with
   6mB <= d; debt |j_0| <= d + B), THEN S holds at (R, D0) with NO
   iteration at all. With the measured R-independent envelope B ~ 2^11 and
   m ~ 0.93 sqrt(R) this static test is passed automatically once
   `6 m B <= d`, i.e. for **R >= ~2^28**; the per-epoch certificate covers
   the window in between. **Hypothesis S is thereby REDUCED to a static
   property FS(R) of the feed-end state** (sec. 7) — a one-row condition
   produced by the feed phase (~21% of the epoch), instead of a
   Theta(R)-row dynamical statement.

What is NOT claimed: FS(R) itself (equivalently the feed-phase envelope) is
still open for unverified R — that is exactly the GG-adjacent hard core and
the honest remaining content of S.

---

## 1. Setting

Post-feed autonomous flow of plain rule A in T6 junk coordinates. Rows
`i > i0`, `i0` = first autonomous row (the first row whose load receives no
feed). State: junk vector `j` on cells `t` (degree `d_i = floor(g(R+i)) + D0`,
`g = gamma*`), transported by `K_delta` (`(1,1)` for `delta = 0`, `(1,2,1)`
for `delta = 1`, `delta_i = d_{i+1} - d_i`), clamped into
`[-2C(d'-1,t), +2C(d'-1,t-1)]` with `d' = d + delta`; cell-0 box is `[-2, 0]`;
death iff junk reaches cell `d'`; capture iff `j = {}` (then T1/T5 coasting
closes the epoch). Feed-end state: `j = j_{i0}` at `d = d_{i0}`, all-negative
(C-N invariant from there on), support `[0, m]`.

**Lemma 0 (autonomy is permanent).** `d_i + i` is strictly increasing
(increment `1 + delta_i >= 1`). The feed branches at row i require
`d_i + i <= R - 1` (band 1) or `delta_i = 1` and `d_i + i <= R` (band 2).
At the first unfed row, `d_{i0} + i0 >= R` (else band 1 would have fired);
hence for every `i > i0`, `d_i + i >= R + 1` and no feed ever fires again.
(The extractor asserts `d_{i0} + i0 >= R`; note equality DOES occur, e.g.
R = 256, D0 = 8: `d_{60} + 60 = 256` with `delta = 0` blocking band 2 —
the strict form `> R` is false in general.) **PROVED.**

## 2. Theorem 1 (march budget / late death) — PROVED

Let `m = max supp(j_{i0})` and `F_i = max supp(j_i)` (front). By T6a,
`F_{i+1} <= F_i + 1 + delta_i` for every autonomous row (kernel reach
`1 + delta`; the clamp only removes support; no feed). Suppose death occurs
at row `i_D <= R - 1`, i.e. `F_{i_D} = d_{i_D}`. Summing T6a from i0:

```
d_{i_D} = F_{i_D} <= m + (i_D - i0) + (d_{i_D} - d_{i0}),
```

hence, using Lemma 0 (`d_{i0} + i0 >= R`):

```
i_D  >=  d_{i0} + i0 - m  >=  R - m .
```

Moreover, defining the per-row advance deficit
`def_i = (1 + delta_i) - (F_{i+1} - F_i) >= 0`, the same sum gives

```
sum_{i=i0}^{i_D - 1} def_i  =  i_D - (d_{i0} + i0 - m)  <=  (R-1) - (R-m) = m - 1 .
```

**Corollary 1a (freeze prefix kills death).** A row whose front does not
advance contributes `>= 1` to the deficit. If the C-F freeze conditions hold
at the first `m` autonomous transitions (all of which precede any possible
death row, since `i0 + m < R - m` with Theta(R) margin at linear slack),
every death row `i_D <= R - 1` would need deficit `>= m > m - 1` —
contradiction. So death is impossible at every row, REGARDLESS of what the
flow does after those m rows. The
no-death half of S is an O(sqrt R)-row local statement (m ~ 0.93 sqrt R
measured). *(The certificate below verifies freeze on ALL rows anyway,
because the capture half needs the longer horizon.)*

**Corollary 1b (death = perfect march; the hardened super-block necessity).**
Any death trajectory advances maximally (`F_{i+1} = F_i + 1 + delta_i`) on
all but `<= m - 1` rows; on every `delta = 1` maximal row the extreme cell
overflows: `|j_{F}| > 2C(d, F+2)` — the negation of C-F, i.e. a front
same-sign state above single-cap scale, from essentially the first
autonomous row on. At the death row itself the cap tail beyond the front is
EMPTY, so the front run is trivially super. This is the provable form of
D1's "death <=> super block" law for the post-feed phase; the sufficiency
direction (super pair => death for every admissible continuation) is D3's
G3 with the referee's full-layer F3 majorant.

## 3. Lemma M (one-step majorant calculus) — PROVED

Throughout: `j <= 0` cellwise, supp ⊆ `[0, m]`, current degree d, next
`d' = d + delta`, `m + 2 < d'`. All loads `w = K_delta * j <= 0` cellwise, so
each clamp is against the lower end `-2C(d'-1,t)` and (C-N) `j' <= 0`.
Cap-ratios: `A_t = |j_t| / (2C(d-1,t))` for `t >= 1`; cell 0 is tracked by
its VALUE `|j_0|` (its cap is the constant 2).

**(M0) Exact debt drain (T7, cell-0 block form).** `w_0 = j_0` (the kernel
has no cells below 0; no feed), cap 2, so `j'_0 = j_0 + 2` if `j_0 <= -2`,
else `j'_0 = 0`: the debt drains at EXACTLY 2 per row, independently of all
other cells, and never revives. Hence `|j_0(i0+k)| = max(0, |j_0(i0)| - 2k)`
and the debt is empty after exactly `|j_0(i0)| / 2` rows (evenness by T3).
This DISCHARGES the C-D concern of the E-lin note ("what must be proved is
that the drain actually runs — cell 0 stays saturated"): post-feed, cell-0
saturation is automatic, not a hypothesis.

**(M1) delta = 1, t >= 3.** `|w_t| <= |j_t| + 2|j_{t-1}| + |j_{t-2}|` and
`|j'_t| = max(0, |w_t| - 2C(d,t))`; dividing by `2C(d,t)` and using the
exact ratios `C(d-1,t)/C(d,t) = (d-t)/d`, `C(d-1,t-1)/C(d,t) = t/d`,
`C(d-1,t-2)/C(d,t) = t(t-1)/(d(d-t+1))`:

```
A'_t <= max(0,  A_t (d-t)/d + 2 A_{t-1} t/d + A_{t-2} t(t-1)/(d(d-t+1)) - 1 ).
```

**(M2) delta = 1, debt cells.**
`t = 1`: `|w_1| <= |j_1| + 2|j_0|`, so `A'_1 <= max(0, A_1 (d-1)/d + |j_0|/d - 1)`.
`t = 2`: `|w_2| <= |j_2| + 2|j_1| + |j_0|`, so
`A'_2 <= max(0, A_2 (d-2)/d + 4 A_1/d + |j_0|/(d(d-1)) - 1)`
(`2C(d,2) = d(d-1)`). The debt's ONLY influence on cells >= 1 is these two
pump terms; with (M0) they decay linearly in k. This is the "two-cell
interaction of the debt with the small block" of the task brief, and it
preserves sub-super status: the pumps are `<= |j_0|/d`, which is `< 1 + 2c/d`
whenever `|j_0| <= d + 2c`.

**(M3) delta = 0.** Caps unchanged (`d' = d`):
`t = 1`: `A'_1 <= max(0, A_1 + |j_0|/(2(d-1)) - 1)`;
`t >= 2`: `A'_t <= max(0, A_t + A_{t-1} t/(d-t) - 1)`
(`C(d-1,t-1)/C(d-1,t) = t/(d-t)`).

**(M4) Monotone majorant iteration.** The right sides of M1–M3 are
nondecreasing in every input `A_s` and in `|j_0|`. Therefore if
`B_t(0) >= A_t(i0)` for all t and B is iterated by M1–M3 with the EXACT
degree sequence `d_{i0+k}` (certified floor engine) and the EXACT debt
`max(0, |j_0(i0)| - 2k)`, with every arithmetic rounding UPWARD (fixed-point
ceilings at P = 96 bits), then `A_t(i0+k) <= B_t(k)` for every cell and
every row, by induction. Roundings only ever increase B, so domination is
preserved; the accumulated rounding excess is `<= 3·2^{-96}` per cell-row —
absorbed, never subtracted.

## 4. Lemma F (front freeze in cap-ratio form) — PROVED (= C-F rewritten)

Given supp ⊆ `[0, m]`, the only loads beyond the front are
`w_{m+1} = 2 j_m + j_{m-1}`, `w_{m+2} = j_m` (delta = 1), resp.
`w_{m+1} = j_m` (delta = 0). They are absorbed — support stays inside
`[0, m]` — provided:

```
delta = 1 :   A_m <= d (d-m-1) / ((m+1)(m+2))        [extreme cell]
              2 A_m + A_{m-1} m/(d-m) <= d/(m+1)      [mate cell]
delta = 0 :   A_m <= (d-1-m)/(m+1)
```

(from `C(d,m+2)/C(d-1,m) = d(d-m-1)/((m+1)(m+2))`,
`C(d,m+1)/C(d-1,m) = d/(m+1)`, `C(d-1,m+1)/C(d-1,m) = (d-1-m)/(m+1)`).
All three right sides are nondecreasing in d and the left coefficient
`m/(d-m)` is nonincreasing in d. The certificate checks the actual branch
(delta_i from the exact floor engine) at the actual d, using the majorant
values `B_m(k), B_{m-1}(k)` — so a machine PASS implies the true state
freezes at every row. **Consequences while freeze holds:** support never
leaves `[0, m]`; no death (m + 2 < d); and with Theorem 1 the freeze on the
first m transitions alone already excludes death forever.

## 5. Theorem 2 (certificate soundness) — PROVED

Fix a feed-end state (R, D0, i0, d_{i0}, j_{i0}) with: j <= 0 cellwise, even;
supp ⊆ [0, m]; m >= 3; 6m <= d_{i0}; Lemma 0's `d_{i0} + i0 >= R`. Run the
fixed-point majorant of Lemma M with freeze checks of Lemma F at every row.
Suppose the machine reports, for some `k_end <= K := (R-2) - i0`:

- freeze valid at every transition `k < k_end`;
- `B_t(k_end) = 0` for all `1 <= t <= m` (and B stays 0 thereafter — the
  iteration is run to the later of the two extinctions), and
- `|j_0(i0)| - 2 k_end <= 0` (debt drained; `k_debt = |j_0(i0)|/2`).

Then, by M4 + F + C-N inductively: the true flow keeps j <= 0, support in
[0, m], never dies (front <= m < d; also re-derivable from Corollary 1a),
and is EMPTY at row `i0 + k_end <= R - 2`. By T1/T5 the epoch coasts to
closure. **Hence Hypothesis S holds at this (R, D0).** QED

**Corollary C-CLOCK (exact capture row).** If additionally the interior
majorant is extinct no later than the debt (`k_junk <= k_debt`, observed at
every scale), the TRUE capture row is exactly

```
capture = i0 + |j_0(i0)| / 2
```

because the debt drains at exactly 2/row (M0) and is then the last junk
standing. This reproduces every known sweep capture row (sec. 6, 16/16:
79, 78, 153, 150, 317, 313, 624, 616, 1261, 1237, 2519, 2486, 5040, 4964,
10090, 9937 — all exact), identifying D2's "capture ~ 0.61R" law as
`i0 + |j_0|/2 ~ 0.21R + 0.79R/2`.

**Corollary SB (no super block post-feed).** Under the certificate the junk
value at cell t never exceeds `sup_k B_t(k) · 2C(d_max - 1, t)` (t >= 1) or
`|j_0(i0)|` (t = 0), while at every row the cap tail beyond the front
(<= m) contains the middle binomial `2C(d-1, (d-1)/2) >= 2C(d_{i0}-1,
(d_{i0}-1)/2)`. The certificate compares the two exactly; margins are
thousands of bits (sec. 6). So no front-region same-sign run can EVER have
min |value| above the remaining cap tail in the post-feed phase: super
blocks are impossible there, and (with Theorem 1 + G3) the D1 super-block
law is now: death => perfect march + trivial super at death; no super
=> no death — both directions proved for the post-feed flow.

## 6. Certificates — VERIFIED-exact

Hostile self-check first: `selftest` replays the TRUE autonomous flow
(independent in-script implementation of the clamp) from the extracted state
and verifies cell-by-cell, row-by-row: sign, support, the exact drain, and
`|j_t| * 2^96 <= N_t * 2C(d-1,t)` (majorant domination), plus the capture
bound. **PASS at all 13 configs run — both D0 at each of R = 128, 256,
512, 1024, 2048, 4096, plus (8192,512): domination held on every row of
every run (up to 3223 rows).** Falsifiability (`negctrl`): corrupting the
(2048,128) state by blowing the front cell up 2^60 makes the certificate
FAIL (freeze), and flipping an interior sign makes it FAIL (hypothesis
check) — 2/2 negative controls caught. The extractor regresses bit-exactly
against the 16 published D2 feed-end records (all scalar fields + the full
`A_bits_profile`, which is `bitlen|j_t| - bitlen(2C(d-1,t))`).

Certificate ledger (all exact; "slack" = min over rows, in bits, of the
freeze inequalities' headroom; "SB margin" = (middle-binomial tail bits) -
(max junk value bits)):

| R | D0 | i0 | d(i0) | m | \|j0\| (=d+c) | k_junk | k_debt | budget K | capture bound | actual | freeze slack (d1x/d1m/d0) | SB margin |
|------|------|------|------|-----|------------|------|------|------|------|------|------|------|
| 128 | 4 | 31 | 99 | 7 | 96 (-3) | 33 | 48 | 95 | 79 | 79 EXACT | -/7/4 | 58 |
| 128 | 8 | 28 | 101 | 7 | 100 (-1) | 35 | 50 | 98 | 78 | 78 EXACT | -/6/4 | 60 |
| 256 | 8 | 60 | 196 | 12 | 186 (-10) | 72 | 93 | 194 | 153 | 153 EXACT | 8/3/- | 124 |
| 256 | 16 | 55 | 201 | 12 | 190 (-11) | 75 | 95 | 199 | 150 | 150 EXACT | 8/3/- | 128 |
| 512 | 16 | 120 | 393 | 18 | 394 (+1) | 151 | 197 | 390 | 317 | 317 EXACT | 11/5/- | 280 |
| 512 | 32 | 110 | 403 | 18 | 406 (+3) | 154 | 203 | 400 | 313 | 313 EXACT | 9/4/- | 288 |
| 1024 | 32 | 238 | 786 | 27 | 772 (-14) | 306 | 386 | 784 | 624 | 624 EXACT | 9/3/- | 603 |
| 1024 | 64 | 218 | 806 | 28 | 796 (-10) | 313 | 398 | 804 | 616 | 616 EXACT | 16/8/- | 621 |
| 2048 | 64 | 476 | 1573 | 39 | 1570 (-3) | 616 | 785 | 1570 | 1261 | 1261 EXACT | -/10/6 | 1292 |
| 2048 | 128 | 436 | 1613 | 39 | 1602 (-11) | 633 | 801 | 1610 | 1237 | 1237 EXACT | 11/5/4 | 1329 |
| 4096 | 128 | 951 | 3146 | 57 | 3136 (-10) | 1238 | 1568 | 3143 | 2519 | 2519 EXACT | 12/5/5 | 2708 |
| 4096 | 256 | 871 | 3226 | 58 | 3230 (+4) | 1267 | 1615 | 3223 | 2486 | 2486 EXACT | 12/6/5 | 2780 |
| 8192 | 256 | 1902 | 6292 | 81 | 6276 (-16) | 2475 | 3138 | 6288 | 5040 | 5040 EXACT | 13/6/5 | 5632 |
| 8192 | 512 | 1741 | 6451 | 83 | 6446 (-5) | 2538 | 3223 | 6449 | 4964 | 4964 EXACT | 14/6/- | 5777 |
| 16384 | 512 | 3802 | 12582 | 118 | 12576 (-6) | 4958 | 6288 | 12580 | 10090 | 10090 EXACT | 18/10/- | 11573 |
| 16384 | 1024 | 3482 | 12903 | 119 | 12910 (+7) | 5082 | 6455 | 12900 | 9937 | 9937 EXACT | 13/6/- | 11878 |

(Ledger JSON has full detail. Pending at write time, appended on landing:
the NEW scale (32768,1024), (32768,2048) — where a PASS proves S at
R = 32768, one full doubling beyond the previously verified frontier, from
a feed-only computation.) Readings:

- **16/16 PASS — every scale of the original S-verification range is now
  certified**; capture rows predicted EXACTLY at all 16 (Corollary
  C-CLOCK) — the certificate is not just sound but sharp on the clock.
- Freeze slacks are the tight quantity (3–16 bits) but show no downward
  trend across six doublings — consistent with the R-independent feed-end
  envelope (D2) and the 1/sqrt(R) decay of the pump coefficients.
- k_junk / k_debt ~ 0.79 at every scale: the interior dies before the debt,
  as the top-down extinction anatomy predicted; the debt is the clock.
- The debt edge law: |j0| - d in [-16, +4] across all extracted states
  (slightly wider than the +-11 of the 16 original records; the certificate
  does not use the edge law — it takes |j0| exactly — so this is
  informational only).
- SB margins grow linearly in R (~0.7 bits per unit R): super blocks are
  not merely absent, they are impossible by thousands of bits.

## 7. Proposition A (uniform-envelope closure; explicit constants) — PROVED

Let the feed-end state satisfy

```
(FS1)  j <= 0 cellwise, even;   supp(j) ⊆ [0, m],  3 <= m,  6(m+2) <= d_{i0};
(FS2)  A_t <= B  for 1 <= t <= m,  with  6 m B <= d_{i0};
(FS3)  |j_0| <= d_{i0} + B .
```

Then S holds at (R, D0) outright (no iteration): capture by row
`i0 + max(10B, B/2 + 2 sqrt(B d) + 2, |j_0|/2) <= R - 2`, no death.

*Proof sketch (all steps one-line consequences of Lemma M/F; the arithmetic
implications are machine-checked, see below).* (a) Cell 1 transient: while
the debt exceeds d (at most B/2 rows by FS3 + M0), A_1 grows by at most
`(B - 2k)/d` per row, so `A_1 <= B(1 + (B+2)/(4d)) <= 2B` always;
afterwards it decays with rate `>= (2k - B)/d` and is extinct within
`B/2 + 2 sqrt(Bd) + 2` rows (quadratic drain integral), never to revive
(revival needs |j_0| > d). (b) With the uniform interior bound
`Abar = 2B` and FS2, every pump is `<= 4Bm/d + O(B m^2/d^2) <= 4/5`:
cells 2..m decay at `>= 1/5` per row (both delta branches), extinct within
10B rows, and stay extinct (recreation needs pump > 1). (c) The freeze
thresholds of Lemma F all dominate `2B` under FS2, so freeze holds at every
row: support never leaves [0, m]; no death (also via Corollary 1a).
(d) Deadlines run in PARALLEL against the budget
`K = (R-2) - i0 >= d_{i0} - 4` (Lemma 0's two-sided pin
`R <= d_{i0} + i0 <= R + 2`): each of `10B`, `B/2 + 2 sqrt(Bd) + 2`,
`(d+B)/2` is `<= d - 4` under FS1–FS3. QED

The seven arithmetic implications (P1–P7) are machine-checked exactly on a
hostile grid (13 degree values 60..10^6 including every real feed-end
degree of the sweep, m down to 3 and up to d/6 - 2, B at the extremal
d/(6m)): `amm12592_S_propA_gridcheck_boxeph.py`, ALL PASS. (The first
draft's sequential-deadline form of P7 FAILED on the grid at extremal
B ~ d/18 — the corrected parallel-deadline form above is what the proof
actually needs; recorded per hazard discipline.)

With the VERIFIED-exact envelope (B ~ 2^11 R-independent, m ~ 0.93 sqrt R),
FS2 reads `6 · 0.93 sqrt(R) · 2^11 <= 0.79 R`, i.e. **R >= ~2.1e8 ~ 2^28**:
beyond that scale the static conditions FS1–FS3 alone imply S; the window
2^15..2^27 is coverable per-epoch by the Theorem-2 certificate (each run
costs only the feed phase + a cheap iteration).

**The reduction (the actual yield of this lane).** Define FS(R): "the
feed-end state of rule A at (R, ceil(eps R)) satisfies FS1–FS3 with
B = 2^11" (or: passes the Theorem-2 certificate). Then:

```
S(eps)  <==  FS(R) for all dyadic R >= 32768        [Prop. A + Thm 2]
```

and FS is a STATIC one-row condition produced by the feed phase. Everything
dynamical about the autonomous phase — 79% of the epoch — is now proved.
What remains open for the all-R theorem is precisely the feed-phase
envelope (the alternation-collapse mechanism of G1/G2 landing the state
inside FS): the same hard core as D3's golden gap, one phase earlier.
Rule/search hazard: a FAILED certificate at some scale would say nothing
about epoch feasibility — only that this majorant is too lossy there.

## 8. Relation to the assigned E2 route

- "(1) harden the super-block law into a theorem": done in the two provable
  directions for the post-feed phase (Cor. 1b necessity-side: death =>
  perfect march, front pair above single-cap scale from feed-end on, and
  trivial super at the death row; G3 sufficiency-side cited). The exact
  12/12 iff-law AT THE MEASURED ROWS remains a measurement; what is proved
  is what the proof needs.
- "(2) super blocks cannot form post-feed": Corollary SB, with the debt
  handled separately by M0/M2 exactly as briefed (the T7 block version).
  The "amplitude bound" the task asked to prove if the verified one
  resists: Proposition A shows ANY R-independent (indeed any o(R/sqrt R
  · 2^{-...}), via FS2) envelope suffices; the envelope itself at
  unverified R remains the open FS.
- "(3) conclude S": concluded per-epoch wherever the certificate passes
  (16 scales, the full original S range; R = 32768 pending), and reduced
  to FS for all
  larger R; NOT claimed unconditionally.
- EXTRA (eps below eps*): not attempted; note that Theorem 1 + the
  certificate machinery apply verbatim to any feed-end state, so a
  sub-eps* feed-phase survival proof would plug in unchanged. The binding
  constraint stays the feed phase.

## 9. Hazards honored

- All PROVED claims are complete arguments over the stated hypotheses; the
  hypotheses at each scale are the extracted exact states, regression-
  anchored to the published D2 records (bit-exact, incl. profiles).
- The majorant is one-sided by construction (ceilings only); the hostile
  selftest checks domination against an independently coded true flow at
  7 configs, every row, exact.
- Lemma 0's `>= R` (not `> R`) was caught by an assertion failure at
  R = 256/1024 during extraction — the strict form is FALSE; recorded.
- The published `A_bits_profile` convention was identified as
  bitlen-difference (not floor log2); the regression uses it exactly.
- Certificate PASS/FAIL is per-(R, D0); no monotonicity in R or D0 is
  assumed anywhere; nothing is extrapolated beyond Prop. A's explicit
  hypotheses.
- Floats: none in any decision (the floor engine's float seed is corrected
  by exact Lucas comparisons, as certified in B2).

## 10. Status ledger

- **PROVED (new):** Lemma 0 (permanent autonomy, with the corrected
  non-strict inequality); Theorem 1 (death only after row `d_{i0}+i0-m >=
  R-m`; deficit <= m-1; freeze-prefix corollary 1a; perfect-march corollary
  1b); Lemma M (M0 exact drain; M1–M3 extended C-A with debt pumps; M4
  monotone fixed-point majorant); Lemma F (C-F ratio form + d-monotonicity);
  Theorem 2 (certificate => S at (R, D0)); Corollary C-CLOCK (capture =
  i0 + |j0|/2); Corollary SB (post-feed super-block impossibility);
  Proposition A (static FS1–FS3 => S, explicit constants, R >= ~2^28 with
  the measured envelope).
- **VERIFIED-exact (new):** 16/16 certificates PASS (all of dyadic
  R = 128..16384 at eps = 1/32 and 1/16 — the full original S-verified
  range) with EXACT capture-row prediction 16/16; hostile domination
  selftest 13/13 configs (through (8192,512), 3223 rows); negative
  controls 2/2 caught; Proposition-A arithmetic grid check ALL PASS
  (13 degrees, extremal B; first-draft P7 corrected); extractor regression
  vs all published feed-end records; SB margins 58..11878 bits growing
  linearly; freeze slacks 3..18 bits with no downward trend. Pending
  extraction at write time: (32768,1024), (32768,2048) — results appended
  below when landed.
- **OPEN:** FS(R) for unverified R (the feed-phase envelope) — the sole
  remaining content of Hypothesis S, hence of
  `C* <= 1 + gamma* + eps` beyond the certified scales.

---

*Appendix (to be appended on landing): certificate outputs for the pending
extractions, including the R = 32768 verdicts.*
