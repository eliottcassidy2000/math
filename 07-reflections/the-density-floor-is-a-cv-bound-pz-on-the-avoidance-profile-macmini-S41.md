# The density floor is a CV bound: Paley–Zygmund on the avoidance profile U

*mac-mini-2026-07-07-S41 (HYP-4837). Owner prompt: deeply understand the LRC14 history,
its corrections and validity; find what we've been missing; then investigate. This note
records the k=13-leg finding: the tail floor `mu_{1/7} >= m_P` reduces to a coefficient-of-
variation bound on ONE positive random variable — the window-avoidance profile `U` — with
4.6x empirical slack; the mean/E[maxgap] detour is unnecessary. Plus two documented dead
ends and a meta-pattern. Scripts: `lrc14_Uprofile_pz_ladder / lrc14_EU_balanced_lattice /
lrc14_EU_floor_mechanism _macmini_S41` (+ .out).*

## The chain (all steps rigorous; constants sharp)

For a 13-set `E` and `x ~ U[0,1)`, let `U(x) = U_{1/7}(x) = sum_j (g_j - 1/7)_+` over the
13 circular gaps `g_j` of `{frac(e x)}` — equivalently (Fubini, exact) `U(x) = meas{s :
the arc (s, s+1/7) contains no orbit point}` (opus-S131's object).

1. `{U > 0} = {maxgap > 1/7}`, so `mu_{1/7}(E) = P(U > 0)`.
2. **Paley–Zygmund / Cauchy–Schwarz:** `mu_{1/7} >= E[U]^2 / E[U^2] = 1/(1 + CV(U)^2)`.
3. **Pointwise cap** `U <= 1 - 1/7 = 6/7` (if any gap exceeds 1/7, the leftover mass is at
   most `1 - 1/7`), hence `E[U^2] <= (6/7) E[U]` and **PZ >= (7/6) E[U] always** — the PZ
   bound dominates the first-moment bound at every family (verified on the whole bank).
4. The k=13 skeleton bar is `mu >= m_P = 14249/252252 ~ 0.0565` (THM-530 / hlarge). So the
   whole leg reduces to **either** `E[U]^2/E[U^2] >= m_P` (a CV bound: `CV(U)^2 <= 1/m_P - 1
   ~ 16.7`) **or** the single linear target `inf_E E[U_{1/7}] >= (6/7) m_P = 14249/294294
   ~ 0.04842`.

## The observed floors (adversarial at the RIGHT functionals, 13-point enforced)

| quantity | adversarial inf (descent minimizing IT) | bar | slack |
|---|---|---|---|
| `E[U_{1/7}]` | **0.0938** at `(0,30,36,45,50,54,60,63,70,72,81,90,108)` | 0.04842 | **1.94x** |
| `PZ = E[U]^2/E[U^2]` | **0.2606** at `(0,2,4,5,6,7,8,9,10,11,12,14,16)` | 0.0565 | **4.61x** |

Observed `CV(U)^2 <= 2.84` everywhere vs the 16.7 allowance. Compare the mean route's
margin at its honest bar (monad HYP-4787): +0.0057 and eroding. **The tail route via PZ-on-U
is where the slack lives.** (Concordant with boxeph HYP-4760's strategic call; this supplies
the missing quantitative mechanism.)

## Relation to monad's Chung–Erdős (HYP-4797): same method, continuous limit

monad-S1's mechanism is `P(∪_j A_j) >= S1^2/(S1 + 2·pairSum)` over 14 aligned empty-arc
events. That IS Paley–Zygmund applied to the anchor count `N = Σ_j 1_{A_j}` (since
`E[N^2] = S1 + 2·pairSum`). `U = ∫ 1_{A_s} ds` is the continuous-anchor limit, so **CE =
discretized PZ-on-U** (and `E[U] ~ S1/14`: their S1 ∈ [1.39,1.89]/14 matches my E[U] ∈
[0.094,0.136]). The continuous version needs no anchor choice, and its two moments are
exactly-computable piecewise-rational integrals.

## The balanced-lattice identity (why pairs die), and the dead ends

Deriving `E[U_u]` by expanding the avoidance product and integrating over the window
position `s` (which forces `Σ m_i = 0`):

> `E[U_u] = Σ_{T ⊆ E} (−1)^{|T|} (1−u)^{13−|T|} · L_T(u)`,
> `L_T(u) = Σ_{m ∈ (Z∖0)^T : Σ m_i e_i = 0, Σ m_i = 0} Π_i φ_u(m_i)`, `φ_u(m) = (1−e(−mu))/(2πim)`.

- `|T| = 2` forces `m(e_i − e_j) = 0`: **pairs never contribute** (kps-S59's factoid, one line).
- Every triple `(a,b,c)` contributes through the rank-1 generator `w = (b−c, c−a, a−b)/gcd`;
  3-APs give the primitive harmonic `(1,−2,1)` = the largest terms; spread triples decay
  like `1/|w|^3`. "random big" has deficit ≈ 0 (validates the frame).

**Dead end 1 (truncation is non-perturbative).** At the AP, the weight-3 sum is `−0.56` and
weight-4 is `+0.79` while the TOTAL deficit is `−0.008`: the lattice expansion converges only
by massive cancellation ACROSS weights. Unsigned weight-3 mass (0.37–0.57 at all structured
families) exceeds the entire deficit budget (0.0864). This extends HYP-4767's moat mechanism
to the density side: **unsigned/truncated lattice bounds fail wherever there is structure;
what saves the density side is that PZ needs no lattice sum at all — positivity of a
positive variable, not smallness of a signed one.**

**Dead end 2 (alternating truncation fails).** The empirical sandwich `main+w3 <= E[U] <=
main+w3+w4` held on 4 families but is VIOLATED at the E[U]-minimizer
(upper 0.0204 < E[U] 0.0942) — Bonferroni-by-weight is not a theorem here. Do not pursue.

## The mechanism of the floor (Farey-shell decomposition)

Where does `E[U]` live in x-space (shells around rationals p/q, q <= 8, Farey-cell widths)?

| family | q=1 | q=2..8 | generic |
|---|---|---|---|
| AP `{1..13}` | 36% | 21% | 43% |
| monad record (parity) | 33% | 25% | 42% |
| E[U]-minimizer | 11% | 17% | **72%** |
| random big | 7% | 7% | 86% |

The E[U]-minimizer `(0,30,36,45,50,54,60,63,70,72,81,90,108)` is a **3-adic cascade** (11 of
13 entries divisible by 3, two patches) — the same c-fold interlacing family as the E[maxgap]
records (monad: parity/2-adic) but at modulus 3: it suppresses the small-q rational spikes
and pushes the avoidance mass into the generic region — where the quasi-random background
(iid main term `(6/7)^13 = 0.1348`) keeps it positive. **You can move the avoidance mass but
not destroy it** — that conservation, made uniform, IS the open lemma. The `{U = 0}`
complement is structurally pinned by kps-S6-wf's forgotten net-characterization (nets occur
only at rationals `p/q, q >= 7` whose residues form a strict 1/7-net); near `q <= 6`
rationals, `U >= 1/q − 1/7 >= 1/42` on the whole Farey cell.

## The meta-pattern (reflection proper)

Every LRC floor this project has actually gotten traction on has reduced to a
**coefficient-of-variation bound on an avoidance COUNT**:
- THM-579 (June 29): covering floor `R' > 0` ⟸ `CV(N_R)^2 < m_Q/(1−m_Q)` (14-sheet count);
- monad HYP-4797: `mu_{1/7} >=` CE on the 14-anchor count;
- today: `mu_{1/7} >= 1/(1+CV(U)^2)` on the continuous profile;
- and the metagraph line: THM-589 proved the tournament H-variance = W(n) (the CV^2 of H),
  the same statistic klein-S4 found unbounded for set-dependent `N_R` — the reason the
  Γ₀(N) congruence lift (HYP-3553) was sought.
The sups (the moat) are exactly the objects that admit NO such bound (HYP-4767: average→sup
loses signed cancellation). **Averages have CV certificates; sups have location/three-gap
rigidity.** That is the project's division of labor in one sentence.

## Honest status + handoffs

- NOT a proof: the open lemma is `inf_E E[U_{1/7}] >= 0.04842` (or the CV form), now with
  1.94x/4.61x empirical slack and exact computability. Numeric (200k grid; AP anchors
  reproduce 93/440 and 477/1078 to 1e-5). Exact-rational engines exist (death-star's
  integrator; kps-S6-wf's EWLB arcs) — porting is mechanical.
- k=8..12 legs: untouched by this (monad HYP-4787's verdict stands); THM-579's projection
  trick remains the template for the G_P coupling there — its `E[U·1_{G_P}]` version has
  the s-average killing the same pair terms (noted, not computed; handoff).
- Gap-moment ladder byproduct: `maxgap >= (Σ g^p)^{1/(p−1)}` pointwise; p=2 is kps-S58's
  failed bound; **p=3 clears 1/7 on the whole bank** (min 0.1496) — but no p <= 6 clears
  monad's T* at the records; the ladder inherits the mean route's bar problem. Secondary.
- The 14-point transcription bug (MISTAKE-120): two of my part-2 "records" were 14-point
  artifacts; all published floors above are 13-point-enforced re-runs.
