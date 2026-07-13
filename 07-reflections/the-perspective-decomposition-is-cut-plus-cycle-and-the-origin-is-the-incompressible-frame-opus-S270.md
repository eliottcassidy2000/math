# The perspective decomposition is cut ⊕ cycle — and the origin is the incompressible frame

**opus-2026-07-13-S270.** Prompt (owner): work LRC(14) creatively; explore past repo threads even if
seemingly unrelated; note that a tournament is made of n−1 perspectives, each of n−1 arcs (the star of
one vertex), and among these (n−1)² arc-slots, T(n−2) are double-counted, T(x)=x(x+1)/2.

## 1. The owner's count is the cut ⊕ cycle split of K₁₄

Take the 14-runner system: origin + 13 movers, vertices of K₁₄. Each **mover** i is a *perspective*:
its n−1 = 13 arcs are its relative speeds {v_j − v_i : j ≠ i} ∪ {−v_i}. Summing over the 13 mover
perspectives: (n−1)² = 169 arc-slots. Of these:

- **T(n−2) = T(12) = 78 = C(13,2) pair-arcs are double-counted** — each difference |v_i − v_j| is seen
  from both of its endpoints;
- **13 spoke-arcs are single-counted** — each v_i is seen once (from mover i; the origin's own
  perspective is excluded from the count).

Now the repo's foundation ("everything is the triangle"): for a tournament on n = 14 vertices, the
**base path is the cut space** (n−1 = 13 arcs) and the **tiles are the cycle space**
(m = C(n−1,2) = 78 tiles of the staircase δ₁₂). The owner's decomposition is *literally* this split
transported to LRC(14): **spokes = cut sector (13, single-counted), pair-differences = cycle sector
(78 = tile count, double-counted).** The perspective sum is the GF(2) Cut ⊕ Cycle identity
169 = 13 + 2·78 — the same identity that makes tilings ↔ even graphs work (the cycle-space bijection).

## 2. The peel family is the perspective family

The peeling identity **L = (6/7)|G'_{~v}| − ε_v holds for every peel v** — thirteen identities, one
per perspective. The fleet had been using its two ends separately without naming the family:

- **Core peel (v = 1):** the runner-1 lemma (opus-S265, kps cont.70). Weak — see §4.
- **Far peel (v = max):** THM-731/732 (klein S287/288). Strong — fine v-grid, small disc_v.

The S270 survey (all 13 peels × 5 covering bodies) shows the spectrum: certificate margin M731 is
roughly increasing in v; the min-L residue {1..11,13,84} is certificate-poorest (**only** the far
peel certifies; all twelve others fail); and the mid-band body 2{1..12}∪{13} — no far element at
all — still certifies at 8/13 peels under the true disc (best v=16), though the crude r²/(3v²)
bound fails everywhere on it. The needed sharpening is quantified: ~19× at v=16 where the true
cancellation provides ~31× (~24× vs ~63× at v=24). kps's Dedekind-collapse (their HYP-6495) has room.

## 3. Frame collapse: the double-counted sector degenerates; the origin cannot

**Micro-theorem (rigorous; LRC(≤13) cited per project policy).** If mover i's frame has ≤ 12
*distinct* relative speeds {|v_j − v_i| : j ≠ i, j ∈ 0..13}, then LRC(≤13) applied to the frame
system makes runner i **1/13-lonely** in the original system (gap ≥ 1/13 > 1/14 to all 13 others,
origin included, at some time).

Census (S270, five covering bodies): **11–12 of the 13 movers are degenerate in every tested
covering body.** The deep-well pack center v=6 sees only 7 distinct speeds — LRC(8) (Rosenfeld)
makes it **1/8-lonely**. The deeper in the near-AP pack, the lonelier, provably.

Why the origin is different: frame collapse means two arcs of the *same star* carry equal values —
possible only through the double-counted sector (v_j − v_i = v_i − v_k: an AP-triple centered at
v_i; or v_j − v_i = v_i: a doubling pair). The **origin's star is all spokes**: coincidence would
need v_j = v_k, impossible. So:

> **Collapse lives only in the double-counted (cycle) sector. The origin is the unique perspective
> with no double-counting available — the incompressible frame — and LRC(14) is precisely the
> loneliness of that frame.**

THM-730 (the AP uniquely maximizes additive triples) now reads: the covering-min body *maximizes
double-counted collapse* — nearly every mover's world degenerates to ≤ 13 runners where the settled
conjecture applies, while the origin's world stays irreducibly 14-dimensional. This sharpens
mac-mini-S82's "structure forgets measure": the structure (additive coincidences, relation lattice,
even-graph shadows) lives on the *pairs*; the measure question lives on the *spokes*. The four
severed facets were all pair-sector objects; the residue is a spoke-sector quantity.

## 4. The frozen fan: why the core peel is weak — exactly, not heuristically

Inside runner 1's danger window D₁ = [0, 1/14] ∪ [13/14, 1], a speed w places a lattice point k/w
in the open window iff w ≥ 15. So for the near-AP pack (w ≤ 13):

**Frozen-fan lemma (proved, 3 lines).** For any W ⊆ {2,…,13}: |D₁ ∩ ⋂_{w∈W} S_w| = (s−1)/(7s),
s = min W. Small runners only trim the window edge, and the trims are nested under the smallest
speed. The consecutive carving curve is exactly flat: T_j = 1/14 for all j ≤ 13.

So kps's guessed mechanism ("each extra small runner carves the arc further") is refuted — but the
refinement survives in corrected two-runner form: the carving is done **entirely by the far
element** (the only runner that *laps* within the window; on runner-1's clock the pack is a frozen
fan, the far element a sweeper), and on the crack body it is tight:
|D₁ ∩ S₂ ∩ S₈₄| = 3/49 = |S_rest ∩ D₁| **exactly** (runners 3..11,13 carve zero even jointly);
margin = L = 563/105105. The (6/7) in 3/49 = (6/7)(1/14) is the same 6/7 as Arg B's ε₁ < 6/7
threshold: Arguments A and B are one mechanism — the far-element sweep — seen as measure vs coverage.

Perspective moral: a peel certifies exactly when its grid is *fine relative to the peeled-away good
set* — the far peel is strong because the sweeper's clock resolves the frozen pack; the core peel is
weak because the pack's clock cannot resolve itself. The user-facing summary of §2–4: **certificates
are perspectives, and the good perspectives are the sweepers'.**

## 5. What transcends

- One lattice, three cosets (klein-S285) said the routes share an object; the perspective count
  says *why the origin's coset is the hard one*: it is the single-counted sector.
- The repo's central duality (cut = hierarchy/scores, cycle = even graphs/tilings) has a metric
  incarnation in LRC: cut = spokes = what must stay lonely; cycle = pairs = what collapses. The
  even-graph/tiling machinery (78 tiles at n=14) is not decoration — it indexes exactly the
  double-counted redundancy that makes covering possible at all.
- Exact rational certificates (disc_v < 6|G'|² in ℚ, Bernoulli form over jump pairs) mean THM-731
  is now *kernel-pure per body* — the remaining analytic work is only ever about *uniformity over
  the class*, never about any single body again.

**Files:** `04-computation/lrc14_multirunner_argA_exact_opus_S270.py`,
`lrc14_peel_survey_perspectives_opus_S270.py`, `lrc14_exact_rational_certificates_opus_S270.py`
(+ .outs in 05-knowledge/results/). **HYP-6505, HYP-6510.**
→ THM-730/731/732, HYP-6248/6485/6495 (×2, collided), MISTAKE-141, everything-is-the-triangle,
even-graphs-as-first-class.
