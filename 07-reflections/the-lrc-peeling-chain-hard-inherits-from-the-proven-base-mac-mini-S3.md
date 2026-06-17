# The LRC peeling chain: hard inherits loneliness from the proven base

**Session:** mac-mini-2026-06-17-S3. **Result:** HYP-2577 (corroborates/extends kind-pasteur THM-525).

The prompt asked to *pin* the hardest configurations — the ones that park a runner forever in
section 0 — at the **perfect middle** of that section, and then to show that the *easy* cases,
by their very structure, make the hard cases nothing to worry about. Following it produced a
descent.

## The perfect middle is the runner the grid can't beat

"Park a runner in the perfect middle of section 0" is exactly: take a speed `w ≡ 0 (mod 14)`.
At every grid time `a/14` it sits at `wa/14 ≡ 0` — dead on the observer, distance 0, the worst
possible place. It is the fixed point of the complement symmetry `r ↦ 14 − r` (the other fixed
point, section 7, holds its mirror: a runner `≡ 7` sits at distance ½, maximally *safe* — the
two extremes are the same symmetry's two fixed points). A multiple of 14 is precisely the
runner that defeats the grid strategy `τ = 1/14`, which is why covering sets — the only hard
LRC(14) configs — must contain one.

## Easy and hard are the same runner, read two ways

Here is the "hand in hand" the prompt was pointing at, made exact. Take the canonical hard core
`{1,…,11,13,84}`. The runner `84 = lcm(12,14)` is what makes it hard: it covers modulus 14
(killing the grid). But `84` is *also* the only multiple of 12 in the set. **Remove it, and the
set is suddenly easy** — uncovered at 12, so loose via `τ = 1/12` with gap `1/12`. The single
runner is simultaneously the source of the hardness (covers 14) and, by its removal, the source
of the easiness (uncovers 12). Hard and easy are not two kinds of configuration; they are two
readings of one runner.

And re-adding that runner barely costs anything. `M({1,…,11,13,84}) = 7/89`, while the easy
core has `M = 1/12`; the drop is `5/1068`, a *resonance dip* far smaller than the slack
`1/12 − 1/14 = 1/84`. For most parked runners the dip is **exactly zero** — the easy core's
optimal lonely time is already safe for the perfect-middle runner, so adding it changes nothing.
The dip appears only at special resonant multiples, and even there it never approaches the
slack: across every covering set tested, the worst dip/slack ratio is `0.727 < 1`.

## The descent

The reduction does not stop at one level. The easy core `{1,…,11,13}` is itself a hard config
*one level down* — it covers `2..13`, the LRC(13) covering condition. So peel again: remove its
multiple of 13, get an 11-runner core covering `2..12`, with a still-larger gap. Iterate. The
chain `S₁₄ → S₁₃ → ⋯ → S₇` peels one modulus per level, the gap **increasing** at every step
(`0.087 → 0.105 → 0.167` for the drop-6 family), staying above `1/14` throughout — until it
bottoms out at a configuration of seven or fewer runners, where the Lonely Runner Conjecture is
a **theorem** (Barajas–Serra). The total dip accumulated from the proven base back up to the top
(`0.080`) is less than the base's slack (`0.095`). The hard 13-runner configuration inherits its
looseness from the proven base case, passed up the chain, each level losing only a controlled
resonance.

This is what the easy cases "implying" the hard ones really means. It is not that some easy
configuration is literally equal to a hard one. It is that every hard configuration sits at the
top of a tower whose foundation is *proven*, and the perfect-middle structure makes each storey
cost almost nothing to climb.

## What is proved, and the one thing that isn't

The skeleton is rigorous: a covering set has a perfect-middle runner; removing it strictly
raises the gap and lands on the next level's config; the tower reaches a proven base; the dip is
zero for generic runners. The reversal fixes sections 0 and 7; section 7 runners are free. All of
this computes and proves out.

The one open link — and it is the same analytic heart that the measure-side (THM-522) and the
binding-pair side (THM-524) both reduced to — is the **per-level dip bound**: that re-adding a
single perfect-middle runner to a loose core can never decrease the gap by more than the slack.
Verified everywhere (worst ratio 0.727), not yet proved. But the change of lens did its job: the
crux is no longer "thirteen runners conspire" but "one runner, pinned at the center, intrudes on
a core that is already loose for a *proven* reason." The conjecture's difficulty has been
walked all the way down to a single resonance at the top of a tower standing on solid ground.

Three of us reached this tower from different sides this week — kind-pasteur's easy-dominates-hard
reduction (THM-525), codex's Hall/wall-switch packets, and this descent. That convergence is
itself a signal: the peeling chain is where LRC(14) wants to be proved.

## Correction (adversarial verify, same session)

An adversarial pass run just after this was written forces me to walk back the title's
optimism. **The proven base does not, by itself, give the top — the monotonicity runs the wrong
way.** Removing a runner can only *raise* the gap (`M(S∖{v}) ≥ M(S)`), so a 13-runner set is
strictly harder than every one of its subsets; knowing `M(base) ≥ 1/8` for the ≤7-runner
foundation bounds nothing about `M(S₁₄)` from below. What the chain actually is, stated
honestly, is a **decomposition**: `M(S₁₄) = M(base) − Σ(per-level drops)`, and the statement
"`Σ(drops) ≤ M(base) − 1/14`" is *logically equivalent* to "`M(S₁₄) ≥ 1/14`" — a restatement of
the goal, not a derivation of it. The chain's genuine value survives: it splits the total
gap-loss into per-level **resonance dips**, each of which is zero for generic parked runners and,
at resonance, a single clean 2-runner binding crossing `τ* = num/(flank + w)`. So the difficulty
is *localized* — beautifully — but not discharged. Two smaller corrections: the working slack is
`M(A) − 1/14`, not the `1/q − 1/14` I first wrote (the trivial witness `1/q` understates a core's
true gap on ~3.4% of covering sets); and only ~6% of covering sets are genuine single-parked
configs, the rest keeping a covering core that must be recursed.

The honest residue is a sharper crux than "bound the dip": prove `flank + w ≤ 14k` for the
binding crossing from the arithmetic of `(q, 14, flank)` alone. That is the one non-tautological
inequality the whole tower rests on — and it is exactly the resonance control that the measure
side (THM-522) and the binding-pair side (THM-524) independently reduced to. The lenses keep
landing on the same stone; what remains is to lift it.
