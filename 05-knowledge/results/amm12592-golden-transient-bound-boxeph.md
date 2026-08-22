# AMM 12592 — Lane D3: the golden transient bound — verdict, the alternation-transport calculus, and the super-pair death certificate

Session: boxeph multifront, 2026-08-03 (after Lanes D1 bulk-rule sweep and D2
E-lin theorem; all citations T1–T9', THM-3329, Theorems A/B/C/S to those
notes). All computations exact int; no floats in any decision; no numpy;
sympy not used.

Scripts (04-computation/):
`amm12592_golden_alternation_calculus_boxeph.py` (Lemmas G1/G2 certificates +
mechanism traces), `amm12592_golden_superblock_certificate_boxeph.py`
(Lemma G3 forward death certificate + threshold grid + hostile controls).
Outputs (05-knowledge/results/):
`amm12592_golden_alternation_calculus_boxeph.{out,json}`,
`amm12592_golden_superblock_certificate_boxeph.{out,json}`.

Status labels: **PROVED** (complete argument, machine-checkable algebra),
**VERIFIED-exact** (exact computation at stated scales), **CONJECTURED /
OPEN** (precise statement, open).

---

## 0. Verdict and headline

**Verdict on the assigned target.** The task was: prove f(R) = o(R) closure
for D1's best bulk rule. D1's best rule (desc1) does NOT have flat or
sublinear D0*: its exact thresholds 0/0/4/14/36/87/<=188 at R = 128..8192
(die/close-pinned) carry the same declining-increment LINEAR signature as
plain rule A, shifted down ~2%. **No o(R) transient bound exists to prove on
the measured behavior of any rule in the D1 class**, and per the task's
fallback this note delivers the sharpest provable structure instead:

1. **Lemma G1 (alternation-transport calculus; PROVED).** The T6 Pascal
   kernel acts on exactly-alternating junk as the (1+delta)-th backward
   difference of the magnitude profile — an exact conjugation. Corollaries,
   all exact: phase-preservation iff the difference profile is single-signed;
   the zero-sum law; the **defect-mass law** `L1(K*j) = 2 * one-signed mass
   of D^{1+delta}a` (the transported L1 *equals* twice the convexity-defect
   mass — the exact budget any golden schedule must manage); the ratio law
   (per-row magnitude multiplier `(rho-1)^{1+delta}`, rho = 2 exactly
   marginal for (1,2,1)). The naive envelope form suggested in the task
   brief ("|j| <= M implies |K*j| <= |D^2 M|") is FALSE — counterexample
   recorded (hazard discipline).
2. **Lemma G2 (the initial data are exactly alternating with an exact
   ratio-2 crossover; PROVED + VERIFIED).** Row-0 junk obeys
   `sign(j_t) = (-1)^{d-t}` on its full support [0, R-2-d_0] (T4b
   corollary), its magnitude profile is additively CONVEX at every interior
   cell (9/9 exact grids to R = 8192), and its envelope-ratio crosses 2
   **exactly at t_2 = 2d_0 - R + 3** = (T6b death-delay bound) + 1: cells
   below t_2 are kernel-contracting, cells above are kernel-expanding.
   The T6b length scale is thus ALSO the boundary between the annihilating
   and amplifying bands of the initial data — a new exact identification.
3. **Lemma G3 (super-pair death certificate; PROVED) — D1's immortality
   conjecture upgraded to a theorem.** If after some row the two leading
   junk cells are same-signed and large enough (exact conditions below),
   then **every admissible continuation** — any in-box even choices at
   every later row, not just the rule that produced the state — dies, at an
   exactly computed row. Machine-certified on 6 near-threshold deaths
   (plain 512/1024/2048/4096/8192, desc1 1024): the certificate validates
   at rows 15/21/20/39/39/57 — within a few rows of D1's observational
   first-super — and **predicts the actual death row exactly, 106..1988
   rows in advance, in all 6 runs**. Hostile controls: on the adjacent
   CLOSING runs the precondition never fires (0 rows in 5/5 runs,
   including 8192 D0=192).
4. **The tau = 2/3 threshold (PROVED + VERIFIED).** The closed-form
   sufficient condition for G3's invariant propagation is the binomial tail
   inequality `C(d, L-1) >= sum_{t>=L} C(d,t)`, whose exact threshold is
   `L*(d) = ceil(2d/3) + {0,1}` at every tested degree — the front-position
   threshold tau = 2/3 joins gamma*, 1/phi, and (1-g)/g in the constants
   ledger of the problem.
5. **The golden gap (the single remaining inequality), with margin
   profile.** Everything reduces to one design inequality (Section 6): keep
   `sup_rows m_i / Theta'_i < 1` (no super-pair) through the feed phase at
   slack o(R), then D2's Hypothesis-S machinery closes. Exact margins: at
   D0* the order parameter stays < 1 (never fires, 5/5 scales); at D0*-1 it
   crosses 1 at rows = 0.029R, 0.021R/0.020R, 0.019R, 0.0095R, 0.0070R — a
   declining fraction of R along the critical line.

---

## 1. Lemma G1 (alternation-transport calculus) — PROVED

Setting: T6 junk flow. Junk vector `j` on cells t (degree d), transported to
the next row (degree d' = d + delta, delta in {0,1}) by
`(K_delta * j)_c = sum_{s=0}^{1+delta} C(1+delta, s) j_{c-s}` before feed
(cells {0,1} only) and clamp.

**(G1-a) Transport identity.** If `j_t = (-1)^{d-t} a_t` (any integers a_t,
finitely supported), then, with `D` the backward difference
`(Da)_c = a_c - a_{c-1}`:

```
(K_delta * j)_c = (-1)^{d-c} (D^{1+delta} a)_c .
```

*Proof.* `(K*j)_c = sum_s C(1+delta,s) (-1)^{d-c+s} a_{c-s}
= (-1)^{d-c} sum_s (-1)^s C(1+delta,s) a_{c-s}`. QED

**(G1-b) Phase law.** In the new phase `(-1)^{d'-c}` the image is
`(-1)^delta (-1)^{d'-c} (D^{1+delta}a)_c`: the image is again EXACTLY
alternating (with a global sign) **iff `D^{1+delta} a` is single-signed** —
for delta = 1: additive convexity (or concavity) of the magnitude profile;
for delta = 0: monotonicity. Every violation cell is an alternation defect
(a same-sign adjacent pair seed) in the transported load.

**(G1-c) Zero-sum and defect-mass law.** For finitely supported a:
`sum_c (D^{1+delta} a)_c = 0` (telescoping). Hence, writing x^+ / x^- for
positive/negative parts,

```
sum_c |(K_delta * j)_c|  =  2 sum_c ((D^{1+delta}a)_c)^+  =  2 sum_c ((D^{1+delta}a)_c)^- .
```

The transported L1 of an exactly-alternating profile EQUALS twice its
one-signed difference mass. In particular a profile whose difference profile
is single-signed apart from boundary cells collapses to boundary scale in
one row — the kernel's annihilation of alternation is total, and the entire
survival budget of the flow is (boundary terms) + (convexity defects). This
is the exact form of K4's "alternation-shaping IS the cancellation scheme".

**(G1-d) Ratio law.** If `a_{t-1}/a_t = rho` (constant, inward growth) then
at interior cells `|(K_delta*j)_c| = (rho-1)^{1+delta} a_c` exactly; rho = 2
is exactly magnitude-preserving for (1,2,1) (`(rho-1)^2 = 1`), rho < 2
contracts, rho > 2 expands. (D1's marginality law, now an identity.)

**(G1-e) Hazard record.** The task-suggested lemma "j alternating with
|j_t| <= M_t implies |(K*j)| <= |D^2 M|" is **FALSE**: a = M on even cells,
0 on odd (M = 100, N = 8) gives image magnitude 200 vs |D^2 M| = 100 at
interior cells. The correct exact statements are (G1-a)–(G1-d); the
saturated case a = M is the special case of (G1-a).

**Certificates** (`..._alternation_calculus_boxeph.json`, key `G1`):
identity 400/400 pseudorandom exact trials (both kernels, values up to
10^9), zero-sum 400/400, defect-mass L1 law 400/400, ratio law 8/8
(rho in {2,3,4,5} x delta in {0,1}, exact integer geometric profiles),
counterexample verified (9 violating cells).

## 2. Lemma G2 (initial data: exact alternation, convexity, and the ratio-2 crossover) — PROVED + VERIFIED-exact

Let d = d_0 = floor(g R) + D0 in the T4b window, j = row-0 junk (clamp
overflow of the T4 closed-form load), a_t = |j_t|.

**(G2-a) Sign law (PROVED = T4b corollary).** On the full support
[0, R-2-d]: `j_t = (-1)^{d-t} a_t`, a_t > 0. (T4b step (I): even-phase
cells overflow above, odd-phase below; boundary cell t = R-2-d has even
phase at even R.) The initial data is EXACTLY alternating.

**(G2-b) Ratio-2 crossover (dominant term PROVED; full profile
VERIFIED-exact).** The dominant magnitude term A_t = C(R-2-t, d-t)
satisfies `A_{t-1}/A_t = (R-1-t)/(d-t+1)`, which equals 2 exactly at

```
t_2 = 2d - R + 3   ( = T6b death-delay bound (2d - R + 2) + 1 ).
```

Verified on the full junk profile at 9 exact grids
(R, D0) = (256,0),(256,1),(512,0),(512,5),(1024,0),(1024,15),(2048,38),
(4096,89),(8192,192): the last sub-2 cell and first super-2 cell bracket
t_2 within 1 cell, with ZERO mixed cells anywhere (the comparison
`a_{t-1} vs 2 a_t` flips exactly once). Cells t < t_2: contracting band
(rho < 2); cells t > t_2 (up to the front R-2-d): expanding band (rho > 2).
Note D0 moves t_2 UP at rate +2 while moving the front DOWN at rate -1:
slack shrinks the expanding band from both sides at combined rate 3.

**(G2-c) Global additive convexity (VERIFIED-exact, 9/9 grids).**
`(D^2 a)_t > 0` at EVERY interior cell 2 <= t <= front (zero concave
cells at all 9 grids). Consequences via G1: one kernel step preserves exact
alternation at every cell (the front boundary term
`a_{L-1} - 2 a_L > 0` too, since the front sits in the super-2 band) —
**the pure transport NEVER creates defects on the initial data**. Defects
are created only by the clamp (cap subtraction at partially-absorbed cells
breaks convexity where junk ~ cap, i.e. in the front taper) and by feed
(cells {0,1}). This localizes the defect-creation operator exactly where
D1 measured nucleation (2–3 cells behind the front).

**(G2-d) Defect birth vs die/close (VERIFIED-exact, plain replays).**
First interior alternation defect (adjacent same-sign junk pair):

| run | outcome | first defect (row, cells, front) |
|-----|---------|----------------------------------|
| 512, D0=4 | DIE 121 | row 10, cells 205-206, front 207 |
| 512, D0=5 | CLOSED 312 | row 70, cell 102 (deep interior) |
| 1024, D0=14 | DIE 250 | row 15, cells 401-402, front 403 |
| 1024, D0=15 | CLOSED 639 | row 139, cell 202 (deep interior) |

Dying runs nucleate a FRONT-REGION defect at rows 10/15 (2 cells behind the
front, before D1's first-super at 13/17); closing runs keep the front band
defect-free for Theta(R) rows (their first defects are deep-interior,
feed/drain-related, and harmless). The dichotomy of K4/D1 is now anchored
at the defect-birth event. (These replays are a THIRD independent
implementation of the plain clamp — die rows 121/250 and capture rows
312/639 agree with the certified fast engine and the D1 sweep engine.)

## 3. Lemma G3 (super-pair death certificate) — PROVED

This upgrades D1's CONJECTURED immortality threshold to a theorem, in the
strongest quantifier order: the death conclusion holds for EVERY admissible
continuation of the flow, not only for the rule that produced the state.

**Setting.** State after the clamp of row i_0: junk j at degree d = d_{i_0},
L = max supp j >= 4, front pair `j_{L-1}, j_L` nonzero, same sign s, values
b_0 = |j_{L-1}|, a_0 = |j_L|. (By definition of L there is no junk above L.)
An admissible continuation is ANY choice, at every later row, of even
in-box values u_t (c-boxes [-2C(d'-1,t), +2C(d'-1,t-1)]), junk = load - u
(T2 generalized-clamp conjugacy, D1 sec. 0). Feed acts on cells {0,1} only
and never reaches the tracked cells.

**Three exact ingredients.**

- **(F1, sign forcing.)** At any cell c, if the load w has |w| beyond the
  box end on its own side, `end_s(c) = 2C(d'-1,c)` for w < 0 /
  `2C(d'-1,c-1)` for w > 0, then every admissible u leaves junk of sign(w)
  with `|junk| >= |w| - end_s(c)`; and at ANY cell, |junk| <= |w| + (box
  width), box width = `2C(d',c)`. In particular junk created at an in-box
  cell is at most the box width.
- **(F2, one-directional transport.)** `(K_delta*j)_c` depends only on
  cells c-s, 0 <= s <= 1+delta: junk strictly above a cell never
  influences its load. Only the <= 2 cells immediately above the front can
  ever touch the advancing pair.
- **(F3, layer majorant.)** Define psi_k on the layer (cells > front) by
  psi_{i_0} = 0 and
  `psi'(c) = sum_s K_s psi(c-s) + 2C(d',c)` (seed = box width).
  Then for every admissible continuation, |junk| <= psi cellwise on the
  layer. *Proof by induction:* a new-layer cell c >= L_{k+1}+1 has
  contributors c-s >= L_k + 1 (old layer), so |w_c| <= sum K_s psi(c-s);
  by F1, |junk| <= |w_c| + width <= psi'(c). QED

**The forced diagonal.** Given the state, define delta_k from the Beatty
word (d_{k+1} = d_k + delta_k), L_{k+1} = L_k + 1 + delta_k, and lower
bounds (a_k, b_k):

```
delta=1:  front cell L+2: load lb  fl = a_k - 2 psi(L+1) - psi(L+2)
          mate  cell L+1: load lb  ml = 2 a_k + b_k - psi(L+1)
delta=0:  front cell L+1: load lb  fl = a_k - psi(L+1)
          mate  cell L  : load lb  ml = a_k + b_k
step OK iff  fl > end_s(front cell)  and  ml > end_s(mate cell);
then a_{k+1} = fl - end_s(front),  b_{k+1} = ml - end_s(mate).
```

*(Load identities: front load = j_L + [<= 2 layer cells]; mate load =
2j_L + j_{L-1} + [<= 1 layer cell] (delta=1) or j_L + j_{L-1} (delta=0);
the interior below the mate never reaches these cells by F2.)*

**Lemma G3.** If every step is OK until the front cell reaches d', then at
that row the bottom-cell junk is nonzero under EVERY admissible
continuation — death (T2/T6: survival requires zero junk at cell d) — at
row `i_0 + (d_{i_0} - L_{i_0})` exactly (the gap d-L decreases by exactly 1
per step).

*Proof.* Induction on k, with hypothesis: junk at cells (L_k - 1, L_k) has
sign s and magnitudes >= (b_k, a_k), and junk at every cell > L_k is
<= psi_k cellwise. (At k = 0 the layer is empty since L = max supp; after
the first step the tracked diagonal may have psi-bounded junk above it —
that is exactly what F3 controls.) Given the hypothesis, the two load
lower bounds hold, with sign s (fl > 0 forces the s-part to dominate the
layer disturbance). By F1 the junk at both diagonal cells is forced, sign
s, with magnitudes >= (a_{k+1}, b_{k+1}) — for every admissible u; and by
F3 the new layer is psi_{k+1}-bounded. At the final step the front cell
equals the bottom cell d', whose box is [0, 2]; fl > end forces nonzero
junk there — death at that row or earlier (junk reaching the bottom cell
by any other route only hastens it). Feed never interferes (cells {0,1};
L_k >= 4 and increasing). QED

**Variant (no-injection class NI; PROVED).** For continuations that take
u = w at every in-box cell (plain rule A is in NI), the layer stays EMPTY
from an empty start (zero loads at layer cells stay in-box, junk 0), so
the same chain with psi == 0 applies. desc1/tr* are NOT in NI (vent); they
are covered by the seeded chain.

**Certificates** (`..._superblock_certificate_boxeph.json`; every number
exact). "Precond" = first row with a same-sign front pair with
`min > Theta' = sum_{t>=L-1} 2C(d,t)` (the cheap gate); "noseed"/"seeded" =
first row the psi==0 / full-psi chain validates to death:

| run | actual death | precond row | noseed valid | seeded valid | predicted death | match | D1 first-super |
|-----|------|------|------|------|------|------|------|
| plain 512 D0=4 | 121 | 15 | 15 | 25 | 121 | EXACT | 13 |
| plain 1024 D0=14 | 250 | 21 | 21 | 42 | 250 | EXACT | 17 |
| desc1 1024 D0=13 | 245 | 20 | 20 | 38 | 245 | EXACT | 16 |
| plain 2048 D0=37 | 508 | 39 | 39 | 70 | 508 | EXACT | 27 |
| plain 4096 D0=88 | 1014 | 39 | 39 | 87 | 1014 | EXACT | 31 |
| plain 8192 D0=191 | 2045 | 57 | 57 | 113 | 2045 | EXACT | — |

Margins: noseed min front margin 3–8 bits (the forced overflow beats the
box end by >= 8x at the tightest row); seeded chains validate with 13–24
bit margins and the layer majorant psi stays 3–8 bits BELOW the pair value
at every step (the compounding adversary loses: mass lingering above the
front is eaten at the front's forced speed — the (2r)^offset path-weight
geometry converges for tau > 2/3). **Death rows predicted exactly, 106,
229, 225, 469, 975, 1988 rows in advance — 6/6.** T9's two-phase law
"death = L + gap(L)" is now a per-run THEOREM (for all admissible
continuations), not a measurement.

**Hostile controls (5/5).** On the adjacent CLOSING runs — plain 512 D0=5,
plain 1024 D0=15, desc1 1024 D0=14, plain 2048 D0=38, plain 8192 D0=192 —
the precondition fires on ZERO rows (sup m/Theta' < 1 throughout). No
contradiction is possible (a validated closing run would refute the
lemma); none occurs. Replay outcomes reproduce the known table (die
121/250/245/508/1014/2045, capture 312/639/628/1271/5064) through D1's
regression-certified choose_row.

## 4. The closed-form corollary and the tau = 2/3 threshold — PROVED + VERIFIED

For readability (the chain certificate above is the authority), a
single-quantity invariant: with `m > Theta' = sum_{t=L-1}^{d} 2C(d,t)`,
layer empty, pair same-signed, the chain conditions hold with
`m' >= m - 2C(d, L+1)`, and the invariant propagates
(`Theta' - 2C(d,L+1) >= Theta'_new`) provided

```
(*0)  delta=0 :  C(d, L-1) >= C(d, L+1)      (automatic for 2L >= d+1)
(*1)  delta=1 :  C(d, L-1) >= sum_{t=L}^{d} C(d, t)    [sufficient form]
```

**The exact threshold of (*1)**: `L*(d) := min{L : C(d,L-1) >=
sum_{t>=L} C(d,t)}` satisfies `L*(d) - ceil(2d/3) in {0, 1}` at
d = 60, 100, 153, 306, 612, 1262, 2538, 5090 (VERIFIED-exact; the
geometric majorization `tail <= C(d,L)/(1-r)`, r = (d-L)/(L+1), gives
tau >= 2/3 asymptotically: `tau/(1-tau) >= tau/(2tau-1) iff tau >= 2/3`).
The front position tau = L/d only increases along the march (L grows
1+delta, d grows delta), so the condition is absorbing. Structural
remarks: the initial front tau_0 = (1-g)/g = 0.67228 > 2/3 at D0 = O(1);
the initial front sits at tau >= 2/3 exactly when
d_0 <= (3R-6)/5, i.e. D0 <= (3/5 - g) R + O(1) with
3/5 - g = 0.0020131...; at larger slack the certified runs enter the
tau >= 2/3 regime a few rows into the march (the chain checks every row
exactly, so no asymptotic step is used anywhere).

## 5. What is proved about the D1 rule class (the negative half, made sharp)

Combining G3 with D1's exact sweeps: at D0*(R) - 1, the plain and desc1
flows reach a certified super-pair state at row ~0.007R-0.03R, after which
NO admissible modification of all later choices — bulk, ballot, vent,
lookahead, anything in the T2-licensed class — can avoid death. All
shaping must act BEFORE the certificate validates. The no-return window
(between D1's observational first-super at rows 13..31 and certificate
validation at 15..39 for the wider class) is 2..21 rows across scales.
This proves the D1 sec.-3b reading: the die row is a delayed consequence;
the battle is lost in the first ~2-4% of the epoch. (Hazard: this is a
statement about flow states, not about epoch feasibility; a different
SCHEDULE can avoid ever entering the state — that is exactly the golden
gap below.)

## 6. The golden gap: the single remaining inequality, with margins

**Conditional chain (all links proved or exactly certified except the two
named hypotheses).** Let f(R) = o(R). Suppose:

- **(GG — the golden gap, OPEN).** For every dyadic R >= R_0 there exist
  admissible choices at slack D0 = f(R) keeping, for every row i <= i_feed,
  every same-sign front pair below the G3 threshold —
  `sup_i m_i / Theta'_i < 1` — while reaching a feed-end state satisfying
  D2's certified conditions (N/C/W/F/D/DL);
- **(S', post-feed collapse, OPEN = D2's Hypothesis S adapted).** From that
  feed-end state the autonomous flow captures by row R-2.

Then the epoch closes at slack f(R) for all dyadic R, and by D2 Theorem A
(LIFT) + THM-3329 assembly: **C* <= 1 + gamma* = log_5(5 phi^2)**.

**Why GG is exactly the right target (proved calculus).** By G1, keeping
the flow alternating is not a heuristic: the kernel annihilates alternation
exactly, with survival budget = 2 x (convexity-defect mass) per row
(G1-c); the initial data is exactly alternating, globally convex, with the
contraction/expansion boundary at the T6b scale t_2 = 2d-R+3 (G2); defects
are born only at the clamp taper and the feed cells (G2-c/d); and every
observed death passes through the G3 state, reached via defect growth
(G2-d: the first front defect precedes first-super by 2-3 rows in every
dying run; the converse "death only via G3" is not claimed or needed).
An admissible schedule beats eps_inf iff it manages the defect-mass budget
through the feed phase — cross-row scheduling of WHERE the clamp absorbs,
as D1's handoff anticipated, now with the exact objective function.

**Margin profile of the order parameter** `sup_i m_i / Theta'_i`
(VERIFIED-exact):

| R | at D0* (closes) | at D0*-1 (dies): first crossing row | row / R |
|------|------|------|------|
| 512 | < 1 (never fires) | 15 | 0.0293 |
| 1024 (plain) | < 1 | 21 | 0.0205 |
| 1024 (desc1) | < 1 | 20 | 0.0195 |
| 2048 | < 1 | 39 | 0.0190 |
| 4096 | < 1 | 39 | 0.0095 |
| 8192 | < 1 | 57 | 0.0070 |

At the crossing the state sits within 0-5 bits of Theta' (m_bits vs
theta_bits: 292/292, 600/599, 596/595, 1214/1214, 2452/2452, 4938/4933): the
threshold crossing is razor-thin at every scale, then margins grow
monotonically — the sup is attained on a knife edge, which is why
one-unit D0 changes flip the outcome. A future golden-schedule attempt
has a concrete acceptance test: keep this exact ratio below 1 for the
first ~0.03R rows; everything after is already covered by proved or
D2-certified structure (up to S').

## 7. Hazards honored

- Rule/search negatives are never epoch-infeasibility evidence; G3
  statements quantify over continuations of a REACHED STATE, not over
  schedules.
- The task-suggested envelope lemma is false as phrased; recorded with
  counterexample (G1-e), replaced by the exact calculus.
- The G2-b crossover is proved for the dominant term and certified within
  1 cell for the full profile; the note does not claim the full-profile
  crossover as a closed form.
- desc1's D0* linearity is a measured signature (7 doublings), not a
  theorem; the o(R)-bound target is reported impossible *on the evidence*,
  not proved impossible.
- All decisions exact; the certificate's every inequality is integer
  arithmetic; the only analytic input is the Beatty word via the certified
  floor engine (Lucas/Fibonacci integer comparisons).
- The seeded G3 chain covers every admissible continuation; the noseed
  chain covers the no-injection class only (plain is in it; desc1 is not —
  its run is covered by the seeded chain, 6/6 validated).

## 8. Status ledger

- **PROVED (new):** G1-a transport identity; G1-b phase law; G1-c zero-sum
  + defect-mass L1 law; G1-d ratio law ((rho-1)^{1+delta}; rho=2
  marginality); G2-a initial sign law (T4b corollary); G2-b dominant-term
  ratio-2 crossover at t_2 = 2d-R+3; G3 super-pair death lemma (all
  admissible continuations, seeded layer majorant; NI class with empty
  layer), with exact death row i_0 + gap; F1 sign forcing; F2
  one-directionality; F3 layer majorant; closed-form Theta' corollary with
  (*0)/(*1) propagation.
- **VERIFIED-exact (new):** G1 batteries (400 trials x 3 laws, ratio 8/8,
  counterexample); G2 sign/support/crossover/convexity at 9 grids to
  R = 8192 (crossover within 1 cell, 0 concave cells); L*(d) = ceil(2d/3)
  + {0,1} at 8 degrees; G3 death certificates 6/6 with EXACT death-row
  prediction (512/1024/1024-desc1/2048/4096/8192),
  margins 3-24 bits, psi below pair value throughout; hostile controls 5/5
  (precondition never fires on closures); defect-birth dichotomy (front
  defects rows 10/15 in dying runs; none for Theta(R) rows in closing
  runs).
- **CONJECTURED / OPEN:** GG (golden-gap schedule at o(R) slack); S'
  (post-feed collapse, = D2's S); desc1 D0* linearity (measured);
  "death only via super-pair" converse (not claimed, not needed).

**The reduction chain now stands:** (*) all-R golden theorem
= GG + S' (both precisely stated, with exact margin data), everything
else — transport, initial data, death mechanism, endgame drain, assembly —
PROVED or exactly certified. Unconditional standing results unchanged:
C* <= 1 + gamma* + 1/32 under S(1/32) (D2); attainment
T(n) <= n + 1 + floor(gamma* n) + 14 for critical n <= 2047 (D1 + THM-3329).
