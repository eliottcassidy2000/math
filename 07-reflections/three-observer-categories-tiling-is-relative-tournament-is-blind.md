# Three Observer-Categories: the Tiling Model Is Observer-Relative, the Tournament Model Is Observer-Blind

*mac-mini-2026-06-22-S41. The owner asked: observer-relative vs observer-independent, a possible third
category, and the difference between tournament and tiling models. These are one question. The answer
is a FINENESS HIERARCHY under the affine group, and it explains why kps's H-level analogy for the LRC
breaks (S31m) and why coverage ≠ additive energy (S39). Builds on [[the-even-graph-is-the-tournaments-cycle-half]],
[[the-cut-side-is-classical-clebsch-and-the-permutohedron]], kps S31l/m, codex HYP-2891.*

## The affine group acts; quantities split by what they ignore

The runners live on the circle; the symmetries that matter are the affine group on the speeds:
**translation** s → s+c (= MOVE the observer / change the anchor) and **scaling** s → λs (= change
units / dilate time). A quantity falls into one of three categories by which it is blind to. VERIFIED
(meas-safe, additive energy, difference-multiset on {1..13} vs +5 vs ×3):

| category | translation | scaling | the invariant | the model | fineness |
|---|---|---|---|---|---|
| **1. Observer-relative** | SENSITIVE | invariant | meas(safe) / coverage | **tiling** (anchored) | FINEST |
| **2. Metric-difference** (THIRD) | invariant | SENSITIVE | the gap multiset / winding metric | "metric tournament" | middle |
| **3. Observer-blind / affine** | invariant | invariant | additive energy, H=I(Ω,2), the order, cross-ratio | **tournament** | COARSEST |

- **Observer-relative** pins the observer at 0 and measures distances FROM it; free to dilate (so
  scaling-invariant) but moving the observer changes everything (translation-sensitive). This is the
  LRC safe condition. meas(safe): {1..13}→0, ×3→0 (scale-blind), +5→0.138 (sees the moved observer).
- **Observer-blind / affine** ignores both the observer's position AND the scale — only relations
  among runners up to affine survive: the cyclic ORDER, the tournament H, the additive energy.
  A(E): 1469 for {1..13}, +5, AND ×3. This is the tournament model: S_n-symmetric, no anchor.
- **The THIRD, metric-difference**, is between: blind to the observer (translation-invariant, like the
  tournament) but NOT to scale (the actual gap sizes matter). The difference multiset: same under +5,
  different under ×3. This is the winding tournament's METRIC — its order PLUS the gap widths — which
  the bare tournament (order only, scale-blind) throws away.

## Tournament vs tiling = observer-blind vs observer-relative

The two models the project uses are categories 3 and 1:
- **Tournament model** (full arc-flip, S_n-symmetric, no distinguished vertex) is **observer-blind**:
  it keeps only the order / H / Ω — affine-invariant data. Coarsest.
- **Tiling model** (fixed base path = a chosen anchor) is **observer-relative**: the base path is the
  origin, exactly like the LRC observer. The S38 fact "fixing the base path chooses the cut summand"
  is the same act as "fixing the observer" — the tiling model is the anchored, observer-relative one.

So the project already had the observer-relative/observer-blind split built in as tiling vs tournament;
the LRC observer is the tiling anchor.

## Why this resolves the LRC confusion (S39 + kps S31m)

kps S31m (independently): "the LRC coverage is a FINER invariant than the tournament's score/cycle
structure — the H-level analogy breaks," extremal = the anchor-free consec. Exactly predicted here:
**the LRC coverage is category 1 (observer-relative, finest); the tournament H is category 3
(observer-blind, coarsest).** The H-analogy breaks because H lives two levels up — it cannot see the
observer that the LRC is defined against. Likewise S39's "coverage ≠ additive energy": coverage is
category 1, additive energy is category 3 — different invariance classes, so neither determines the other.

## The finishing redirect

The LRC must be attacked in **category 1 (observer-relative / tiling)**, where the coverage actually
lives — not in category 3 (tournament H / additive energy), which is provably too coarse (kps's broken
analogy; my {2..14} counterexample). Concretely the safe condition is:

> **safe at t ⟺ the observer 0 sits in a gap of {frac(s_i t)} with both bounding runners ≥ 1/14 away.**

This couples category 1 (where the observer falls) with category 2 (the gap widths — the metric-
difference / three-gap structure, kps THM-565). The under-used category 2 is the bridge: the bare
tournament order (cat 3) is too coarse, but the order-WITH-gap-widths (cat 2) plus the observer's
placement (cat 1) is exactly the safe condition. So the right object is the **metric winding tournament**
(order + gaps), anchored at the observer — not the combinatorial tournament, and not the affine additive
energy. The S31n correction sharpens the category-1 target: diff-closed tilers are the AP/dilate
family, but the Goddyn-Wong row is a sporadic tight tiler.  The live rigidity statement is therefore
classification of the observer-relative tight locus, with the three-gap metric (cat 2) as the tool.
codex's Clebsch/Bruhat/unital carriers
(HYP-2891) sit on the cut = observer-relative side, consistent with this placement.

Verified core (S41): the affine-invariance table (the three categories). Synthesis: the fineness
hierarchy, the tournament/tiling = blind/relative identification, the resolution of the broken H-analogy.
Related: HYP-2888 (coverage), HYP-2885 (kps additive-energy), kps S31m (tiling rigidity), HYP-2605 (winding).
