# THM-939 — The relation trap on the chain-dense core (death-star-2026-07-16-S34)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCDenseCoreRelationTrap.lean`,
axioms propext/Classical.choice/Quot.sound; verify the build report in the session
log). The formal bridge from kps's THM-934 generic certification law to codex's
`DenseCoreDissociatedB5Supply`. Source: HYP-7171.

## Statement

On any `ChainDenseCore` family (sorted; ratios ≥ 3 at every pair index above the
dense pair j):

**(A1, `no_low_mass_relation_above_pair`)** no vanishing integer relation
`∑ coeff·w = 0` has top support element `t ≥ j+2` with below-mass
`∑_{i<t}|coeff i| ≤ 2`. Kills every doubling-dilate (`w_t = 2w_s`) and unit-pair
(`w_t = w_s + w_r`) shape above the pair: the ladder's factor 3 beats mass 2.

**(A2, `no_unit_relation_high`)** no vanishing unit-coefficient relation
(all `|coeff i| ≤ 1`) has top support element `t ≥ j+4`: the geometric invariant
`2·∑_{j+1≤i<s} w(i) ≤ w(s) − w(j+1)` (three-line ladder induction) plus the
bottom-block bound `≤ 12·w(j+1)` and the 3-step magnification `w(t) ≥ 27·w(j+1)`
give `∑_{i<t} w(i) < w(t)`.

**Adapter** (`denseCoreDissociatedB5Supply_of_trapped`): the B5 supplier may assume
both trap facts on the sorted speeds — `TrappedDenseCoreB5Supply ⟹
DenseCoreDissociatedB5Supply`, hence (`lrc14_from_four_detuned_and_trapped_B5`)
LRC(14) from citation + four-detuned dispatch + B5 on the TRAPPED core.

## Why this is the kps-connection

THM-934's mechanism (THM-930 leverage identity): B5 dies only under SYSTEMATIC
coordinated relations — dilate chains, AP cores — never incidental ones. A1/A2 are
exactly that mechanism made kernel-checkable on the dense core: above its dense
pair the family is PROVABLY relation-free in the killing shapes (relations are
trapped in the bottom j+4 elements; support-2 shapes die from j+2). The remaining
analytic obligation is the singular-series identification (relation masses → B5),
codex-S28's named separate obligation — now to be discharged only UNDER the trap.

Boundary honesty: mass-3 relations at exact ratio 3 (`w_{t} = 3w_{t−1}`) are the
first NOT excluded — precisely kps's "small-ratio incidental blocker" stratum,
which their law prices at ~5×10⁻³ against μ₀ ≈ 0.135.

## Referee

Same script as THM-938: 4000 random dense-core families — proved inequalities PASS,
A1 brute-force PASS (all coefficient patterns of mass ≤ 2), A2 sampled PASS.
