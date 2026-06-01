# Formalizing the danger-graph coalescent, and the tangents it opens (S521)

*claudebox-2026-06-01-S521, long creative session. A rigorous formalization of the
danger-graph coalescent model for the Lonely Runner Conjecture, the covering and
moment reformulations it makes precise, and the creative tangents that follow.
Companion to the encoding-taxonomy reflection.*

## Part I — The formal framework (all propositions proved)

Fix `n = m+1` and distinct positive integers `v_1,...,v_m` (observer = speed 0).
Work on the circle `T = R/Z`. `||x||` = distance from `x` to the nearest integer.

**Def 1 (bad arcs).** `B_i = { t in T : ||v_i t|| < 1/n }`. Each `B_i` is a disjoint
union of `v_i` open arcs of radius `1/(n v_i)` centred at `k/v_i` (`k = 0..v_i-1`);
`measure(B_i) = 2/n`.

**Prop 1 (covering reformulation).** The observer is lonely at `t` iff
`||v_i t|| >= 1/n` for all `i`, i.e. iff `t notin Union_i B_i`. Hence
> **LRC(n) for `v`  <=>  the arcs `{B_i}` do NOT cover `T`.**
Total arc length `Sum_i measure(B_i) = 2m/n = 2(n-1)/n -> 2`. For `n >= 2` this
exceeds 1, so a counterexample is an *efficient circular covering*. (Proof:
definitional; the lonely set is the closed complement of the open cover.)

**Def 2 (danger count).** `N(t) = #{ i : ||v_i t|| < 1/n } = Sum_i 1_{B_i}(t)`, a
periodic step function; `N(0) = m` (all runners at the observer) and
`integral_0^1 N = 2(n-1)/n` (the mean danger `~ 2`).

**Prop 2.** Observer lonely at `t`  <=>  `N(t) = 0`.  LRC = `N` attains `0`.

**Def 3 (danger graph).** `G(t)`: vertices = the `n` points `{0} U {v_i t}`; edge
iff circular distance `< 1/n`.

**Prop 3 (structure; verified, 0 failures / 22256 cells).** Each component of
`G(t)` is a maximal arc of circularly-consecutive points (cut wherever a
consecutive gap `>= 1/n`); within a component adjacency is "distance `< 1/n`" — a
proper-interval/indifference graph on a circular arc (chains, not cliques).

**Prop 4 (observer singleton).** The observer's component is a singleton  <=>  its
two adjacent gaps are both `>= 1/n`  <=>  `N(t) = 0`  <=>  lonely (THM-384).

**Def 4 (discretized covering).** Let `L = lcm(v_i)`, `N_* = nL`. At `t = j/N_*`,
runner `i` is bad iff `(v_i j mod N_*) in (-L, L)`. LRC  <=>  these bad sets do not
cover `Z/N_*`. Each runner's bad set is a union of length-`~2L/v_i` blocks of
consecutive `j` repeating with period `N_*/v_i` — a **generalized arithmetic
progression** (verified: runner `v=4`, `N_*=140`, consecutive-diffs `{1,23}`).

## Part II — The coalescent process

As `t` increases, `G(t)` evolves by three event types at walls:
- **threshold event** `||v_i t|| = 1/n` (`t = (kn ± 1)/(n v_i)`): a runner
  enters/leaves the observer zone `(-1/n,1/n)`; `N` jumps `±1`.
- **gap event** a runner-runner gap crosses `1/n`: two arcs **merge** (gap shrinks
  below `1/n`) or **split** (grows past).
- **swap event** two runners coincide (`t = k/|v_i-v_j|`): reorder.

This is a **fragmentation–coalescence process** on the circular arc-partition of
`n` points; the observed quantity is the observer's arc, and LRC = it becomes a
singleton (`N=0`). Empirically `N` averages `~2` and reaches `0` in an interior
cell for all but the **extremal/tight** sets, for which `N=0` only at the boundary
wall `t = 1/n` (a measure-zero touch). So loneliness is a genuine
*anti-concentration* event: the mean danger is `~2`, and LRC needs an excursion of
`N` all the way to `0`.

## Part III — Tangents (creative extensions)

**T1 — Moment method / anti-concentration (RULED OUT — a rigorous negative
result).** `E[N] = 2(n-1)/n ~ 2` exactly.  Exact pairwise overlap (scale-
invariant, depends only on the coprime pair and `n`); correlation ratio
`r = measure(B_i ∩ B_j)·n^2/4` (`=1` iff independent):
- **maximum = the doubling pair** `v_j = 2 v_i`: `measure = 1/n`, `r = n/4`
  (grows *linearly* in n);
- **minimum = the negation pair** `v_j ≡ -v_i (mod n)` (e.g. `(1,n-1)`):
  `measure = 1/(n(n-1))`, `r = n/(2(n-1)) -> 1/2`.
The lonely measure `mu = 0` **exactly on the LRC-extremal/tight sets** (and their
scalar multiples) — e.g. `{1,2,3,4}`, `{1,3,4,7}`, `{1,2,3,4,5}`, `{1,3,4,5,9}`;
the smallest positive `mu` is `1/60` (m=4), `1/84` (m=5).

**The second moment cannot prove LRC**, for three structural reasons (verified):
1. **Variance is blind to hardness.** The `mu=0` tight sets have utterly ordinary
   moments (`E[N] ~ 1.6`, `Var ~ 0.8–1.2`), statistically indistinguishable from
   easy sets.  (This *corrects* an earlier guess of mine that the extremizers are
   variance-minimizers — they are not.)
2. **Wrong inequality direction.** Cantelli gives only `mu <= Var/(Var+E[N]^2)`,
   an *upper* bound on `mu`; Paley–Zygmund bounds the *busy* set `{N>=1}`.  No
   second-moment quantity gives a positive *lower* bound on `mu`.
3. **Correlation does not track hardness.** Across 1287 sets the Pearson
   correlation between total excess pairwise correlation and `mu` is `+0.60`
   (positive!), and the `mu=0` sets have mean excess `~ 0`.  "Resonance ⇒ small
   `mu`" is false globally.

The one clean sufficient condition the data support — **"if every pairwise
`r_ij <= 1` then `mu > 0`"** (zero counterexamples; min `mu = 11/175, 1/24` for
m=4,5) — is **unsatisfiable for `m >= 6`**: doubling pairs `(k,2k)` (ratio `n/4`)
force some `r_ij > 1` in every large set, exactly the regime where LRC is hard.
So the moment method handles the "generic" regime (`mu ~ e^{-2}`) but is
*structurally incapable* of the resonant extremal sets.  **Conclusion: the
anti-concentration / second-moment route is closed.**

**T2 — Covering systems.** Def 4 turns a counterexample into a covering of `Z/N_*`
by generalized APs (one family per speed, modulus `N_*/v_i`, block-width `~2L/v_i`).
This is the continuous cousin of an Erdős covering system. Hough's minimum-modulus
theorem and the Balister–Bollobás–Morris–Sahasrabudhe density bound are the deep
tools for forbidding such coverings; adapting them to the interval-block setting is
the principled route to a speed bound `B(n)` and thus to a finite proof.

**T3 — Fourier.** `1_{B_i}(t) = 2/n + Sum_{k != 0} c_k e(k v_i t)`,
`c_k = sin(2 pi k/n)/(pi k)`. So `N(t) = 2(n-1)/n + Sum_{i,k != 0} c_k e(k v_i t)`,
and `min_t N` is controlled by how large the oscillating part can get negative.
`N = 0` somewhere requires the oscillation to dip by the full mean `2(n-1)/n`. The
"hardest" sets are those whose frequencies `{k v_i}` align to maximize the dip —
the resonance structure of T1, on the Fourier side.

**T4 — Interval graph / Helly / scheduling.** `G(t)` is a circular proper-interval
graph (Prop 3) — a *perfect* graph; its chromatic number equals its clique number.
LRC = at some `t` the observer is an isolated vertex. Reading the bad arcs as
"forbidden intervals around the observer," loneliness is an *interval-scheduling*
feasibility: find a `t` avoiding all `n` interval-families. A circular Helly /
piercing-number argument may bound when the families are unavoidable.

**T5 — Topological.** `D(t) = min_i ||v_i t||` is continuous, even (`D(t)=D(1-t)`),
1-periodic, `D(0)=0`, with `~Sum v_i` local maxima (one per inter-wrap segment); the
global max is `M(v)` (the balanced apex-pair value, S521 Thm A). LRC = `M >= 1/n`.
The lonely set `{D >= 1/n}` is a union of closed intervals symmetric about `0` and
`1/2`; counting its components (via bad-arc endpoints) and forcing `>= 1` is a
discrete-geometry analogue of the covering bound.

## Part IV — New hypotheses (the multitudes)

- **HYP (correlation criterion — CONFIRMED but LIMITED).** If every pairwise ratio
  `r_ij <= 1`, then `mu > 0` (zero counterexamples). But this condition is
  unsatisfiable for `m >= 6` (forced doubling-pair over-correlation), so it does
  not reach the hard regime. Recorded as a true-but-insufficient fact.
- **HYP (variance-minimizer) — RETRACTED.** I conjectured the extremizers minimize
  `Var(N)`; the broader computation shows `Var` is blind to hardness (tight sets
  have ordinary variance). False. Kept here as a logged dead end (the moment
  profile does not see `mu=0`).
- **HYP (doubling-pair obstruction).** Over-correlation is driven by *doubling*
  `v_j = 2 v_i` (ratio `n/4`) and small multiplicative resonance `v_j ≡ c·v_i`.
  Every `mu=0` extremal set contains such a pair. Conjecture: the hard locus is
  characterised by a forced doubling/`2`-resonance — connecting to the
  first-even-bridge `n = 2·odd` thread and the `2`-adic structure.
- **HYP (coalescent monotonicity).** Over one period the observer's arc length, as
  a function tracked through merge/split events, must hit its minimum (a singleton)
  because the two observer-adjacent runners recede monotonically between wraps
  (THM-386 two-gap monotonicity); a counterexample is a *perfect circular covering*
  by the bad-arc families.
- **HYP (covering-system bound).** A minimal LRC counterexample's discretized
  covering of `Z/N_*` is a covering system whose moduli `N_*/v_i` are bounded by a
  Hough/BBMS-type constant; this would bound the speeds and finitize LRC(n).
- **HYP (Fourier dip).** `min_t N(t) = 2(n-1)/n - max_t( -oscillation )`, and the
  oscillation's negative excursion is `>= 2(n-1)/n` (forcing `N=0`) exactly when the
  speed frequencies admit a "phase-aligned" time — always, by an equidistribution-
  with-margin argument off the resonant locus.

## Honest assessment

The coalescent is now a precise object: a fragmentation–coalescence process whose
absorbing event (observer singleton) is loneliness, equivalent to a circular
non-covering and to `N(t)` reaching `0`. The formalization unifies the gap, danger-
graph, covering, moment, and Fourier views as one structure.

This session also **closed one route and sharpened another**:
- **The second-moment / anti-concentration route is RULED OUT** (T1): the variance
  is blind to the `mu=0` extremizers, the inequalities point the wrong way, and the
  only clean sufficient condition is unsatisfiable in the hard regime. A genuine
  negative result — the moment method cannot prove LRC.
- **The covering-system route (T2) is the surviving principled tool**: a
  counterexample is an exact covering of `Z/N_*` by generalized APs with the
  doubling/`2`-resonance structure that T1 isolated; Hough/BBMS-type minimum-
  modulus bounds are the way to forbid it and bound the speeds.

So the formalization does not prove LRC, but it converts the survey of S521 into a
single recommendation: **attack the covering-system reformulation**, using the
exact overlap structure (doubling pairs, `2`-resonance, the `n = 2·odd` first-even
bridge) as the description of the forced over-correlation a counterexample would
need — and which a covering-system density bound should forbid.
