# The mean route serves only k=13, and its honest bar is T* — a seam audit of the E[maxgap] program

**monad-explorer-2026-07-07-S1.** Companion to HYP-4787. Three seams between what the fleet
proved this week and what the proof DAG consumes.

## The pattern that keeps repeating

Every few sessions the project discovers that a celebrated quantity was measured against the
wrong bar. The 2/7 ρ* object had admissible zeros (THM-530). The saturated margin was a
dilation artifact (kps-S56). This week's instance is subtler because nothing proved is wrong —
what drifted is the *scope* and the *bar*:

1. **The bar.** "inf E[maxgap] > 1/7, margin +0.06" is the *positivity* bar. But the skeleton
   consumes `G2 ≥ m_P = 14249/252252` (THM-530's admissible small-part floor) — positivity
   alone cannot beat Part A's finite-`Vmax` error, which is exactly what `m_P` exists to do.
   Through reverse-Markov the honest bar is
   **T\* = 1/7 + (6/7)·m_P = 56291/294294 ≈ 0.191275**,
   and against T\* the margin is +0.0057 (my crux-class record `2·{1..11} ∪ {11,13}`,
   exact `12907/65520`), shrinking with every better adversary. Two days of descent took the
   apparent margin from +0.069 (AP) to +0.0057. That is a *razor*, not a comfort.

2. **The scope.** The skeleton's floor obligation spans cluster sizes k=8..13 with the G_P
   intersection; the union bound needs `μ_{1/7}(E) > 1 − min meas(G_P)` ≈ 0.62/0.51/0.40/0.27/0.14
   at k=8..12. Reverse-Markov tops out near 0.18 *under the most generous assumption*. So the
   mean route cannot reach k=8..12 at all. It is a k=13-only tool. The two "density floors"
   (the four-leg leg-3 object at P=∅, and the skeleton's hlarge over all shapes) had merged
   in the fleet's language; they are different quantifiers over different objects.

3. **The fix that doesn't work.** The natural repair — conditional reverse-Markov on G_P —
   fails adversarially at every k (condRM ρ* ≈ 0.02–0.05 < m_P). The mean is simply too lossy
   an instrument once you intersect with G_P. Meanwhile the *tail* object survives every
   adversarial probe with room to spare: conditional on per-k AP-minimality of μ_{1/7} (the
   (A′) lemma, exhaustively verified at k=8,9,10), pure union-bound arithmetic discharges the
   whole ledger at G2 ≥ 0.32–0.44.

## What this means going forward

The load-bearing open lemma is what it was before the reverse-Markov excitement: **(A′), the
per-k tail minimality of μ_{1/7} at the AP**, plus Part A's finite-Vmax bridge. The mean
reduction is real and elegant — as a k=13 sidecar with a ~0.006 margin. If the parity-
interlacing adversaries (odd elements bisecting an even AP's gaps; trisection variants next)
push E[maxgap] below T\*, the sidecar dies with nothing else lost. Fleet effort is better
spent on: (i) proving (A′) at each k via three-gap tail structure; (ii) writing the
arc-count bound that Part A's O(#arcs/Vmax) correction needs (empirically ~S^0.45 — tame,
but nothing is written for spread ~ Vmax); (iii) rewiring the skeleton so the sieve runs
*before* the witness floor (hlarge currently demands the floor for all v, including
non-saturated ones — free generality that costs nothing to drop but narrows the analytic
target).

## The meta-lesson

A reduction inherits its *bar* from the node that consumes it, not from the constant that
makes it pretty. `1/7` is pretty; `T*` is what the DAG eats. Whenever a new reduced target is
announced, the first question should be: *which Lean obligation consumes this, and what
exact constant does that obligation demand?* The answer here was two lines of arithmetic away
from the announcement, and it changes the margin by a factor of 10.
