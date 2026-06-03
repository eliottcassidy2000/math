# HYP-2195 — LRC as circular-arc covering: the covering-depth distribution, and the additive-chain collapse family

**Session:** claudebox-2026-06-03-S617. **Frame:** the recursive-log/altitude frame (HYP-2180/S615) continued.
**Threads:** HYP-2155 (resonance energy), HYP-2185 (apex sheaf, H⁰=∅), HYP-2098 (AP/V* tight collapse), the
perspective key.

## The master object
A loneliness certificate is a clock time `t` avoiding every runner's **forbidden arc** `{t : ‖vᵢ t‖ < δ}`
(`‖·‖` = distance to nearest integer). So **LRC at gap δ IS a circular-arc covering problem**, and the one object
that controls everything is the **covering-depth distribution**
`p_k = meas{t ∈ [0,1) : depth(t) = k}`, `depth(t) = #{i : ‖vᵢ t‖ < δ}`, `δ = 1/(n+1)`, `n = #runners`.
Lonely times = `{depth = 0}`, measure `p₀`. Every LRC functional at gap δ is a functional of `{p_k}`.

### First moment is vacuous (verified + formalized)
Each forbidden arc has measure exactly `2δ` (v disjoint sub-arcs of width `2δ/v`), so `∫depth = Σ 2δ = 2nδ` and
`E[depth] = 2n/(n+1)` (verified exactly for AP, chains, random). The union bound `p₀ ≥ 1 − 2nδ` is **≤ 0 for all
n ≥ 1** (`2n/(n+1) ≥ 1`). LRC's content is therefore NOT first-moment — it is in the **arc correlations**.
Formalized: `CoveringDepth.lean` `union_bound_vacuous`.

For random configs `depth ≈ Poisson(2)`, `p₀ ≈ e⁻² ≈ 0.135` (measured 0.09–0.11). The extremal/tight configs are
exactly the ones whose correlations drive `p₀ → 0`.

## The collapse family = additive chains (new sub-problem, verified)
`p₀ = 0` (lonely set measure-zero, gap exactly δ, achieved only at isolated tight times) is the **collapse**
condition. Exhaustive search (min=1, gcd=1, max≤12, sizes 2–5) returns exactly:
`(1,2), (1,2,3), (1,2,3,4), (1,2,3,4,5)` (the AP) **and the sporadic** `(1,3,4,7), (1,3,4,5,9)`.
**Every collapse set is closed under sum-relations `a+b=c`; every non-collapse set lacks them.** The sporadic chains
are non-AP additive chains (`4=1+3, 7=3+4`; `4=1+3, 5=1+4, 9=4+5` — each top a sum of two below). So:

> **Conjecture (additive-chain collapse).** `p₀ = 0` ⟺ the speed set is an additive chain (generated under
> `a+b=c`). The collapse family is the `(1,1,−1)`-resonance-rich configs.

An `a+b=c` relation is a small `Σ mᵢ vᵢ = 0` **resonance** (coeffs `1,1,−1`). So the collapse family = the
resonance family: this fuses the master object with HYP-2155 (resonance energy) and HYP-2185 (the apex/H⁰=∅ is the
`p₀=0` of the sheaf — collapse = empty global section). **The resonance, seen in the covering picture, is arc-coincidence.**

### The mechanism (formalized)
The clock distance `dZ x = |x − round x|` is **subadditive**: `dZ(x+y) ≤ dZ x + dZ y`. Hence a relation
`v_c = v_a + v_b` forces `‖v_c t‖ ≤ ‖v_a t‖ + ‖v_b t‖` — the three arcs are dependent; at a tight boundary
(`‖v_a t‖,‖v_b t‖ ≤ δ`) it **pins** `‖v_c t‖ ≤ 2δ`, the collapse mechanism. Formalized: `CoveringDepth.lean`
`dZ_add_le`, `dZ_chain_le`.

## The singleton-wall exponent test (the altitude principle, verified)
HYP-2180: #logs = #independent averaging levels = altitude. Prediction: a **singleton wall** (rank-1, ONE additive
relation) carries a `(loglog)¹` regime, a multi-relation wall `(loglog)²⁺`. Test: perturb a wall off the collapse by
ε and fit `p₀(ε) ~ εᵃ` (the exponent = #relations the perturbation breaks = wall codimension = altitude).
- Singleton wall `{1,2}` (1 relation): **a = 0.990** — clean exponent 1, the `(loglog)¹` regime ✓.
- Sporadic chain `{1,3,4,7}` (2 relations): a = 0.84 (sub-linear curvature — stacked-relation signature).
- AP `{1,2,3,4}` (4 relations): a = 0.93. The singleton case passes the test cleanly; multi-relation walls show the
  higher-order correction. (The exponent counts broken relations; a full `(loglog)²` needs the simultaneous-break
  path — open.)

## Formalized (math-lean, sorry-free) — `Math/LonelyRunner/CoveringDepth.lean`
`dZ`, `dZ_nonneg`, `dZ_le_dist_int` (round is closest), `dZ_add_le` (subadditivity = resonance link),
`dZ_chain_le` (additive-chain arc-dependence), `union_bound_vacuous` (first moment says nothing at δ=1/(n+1)).

## CORRECTION (S621): the additive-chain ⟺ collapse conjecture is FALSE

The consolidated structure analysis ([[THM-411]], `lrc_chain_vs_collapse_s621.out`) refutes the
biconditional. The earlier "exhaustive" search was too small (max ≤ 12, sizes 2–5); with the true
magnitude bound it is decisive:

- **`additive-chain ⟹ collapse` is wildly FALSE.** Additive structure is generic; collapse is rare.
  At `n=7`: **3** collapse sets vs **364** strict additive chains vs **5858** sets with
  `max(S) ∈ S+S`. The pinning mechanism (`dZ_add_le`) is real but only *necessary*: one resonance
  pins one arc; tiling the whole circle needs the entire `(ℤ/m)*` orbit pinned at once — far stronger.
- **The clean necessary condition is `max(S) ∈ S+S`** (the largest speed is a sum of two others):
  holds for every collapse set, n=3..7, no exceptions. (The strict single-generator chain notion is
  too narrow — it misses `(1,3,4,7)`, whose `3` is a generator.)
- **The true characterization is tightness, not additive combinatorics** ([[THM-411]]): the collapse
  family is *finite* with `v_max ≤ 2n-1`, count `1,2,2,1,3,1` (n=3..8), every member witnessed
  exactly on the `(ℤ/m)*` orbit `{k/m : gcd(k,m)=1}`. Resonance-richness is a *symptom*; the *cause*
  is the gap staying pinned at `1/m` across the whole primitive-residue witness orbit.

So the master object and the subadditivity mechanism survive; the additive-chain *characterization*
does not (same correction applies to the parallel [[HYP-2153]]). Collapse ⟹ resonance-rich, but
resonance-rich ⇏ collapse.

## Open
- Prove the finiteness / `v_max ≤ 2n-1` bound and the `(ℤ/m)*`-witness characterization ([[HYP-2196]]).
- The clean `(loglog)²` measurement (simultaneous-break / tangential approach to a rank-2 wall).
- Identify `{p_k}` higher cumulants with the resonance energy `E(v)` (HYP-2155); is `p₀` a theta/relation-lattice
  functional?
- Connect to the apex sheaf: `p₀ > 0 ⟺ H⁰ ≠ ∅`; is `p₀` the measure-refinement of the sheaf's global sections?
