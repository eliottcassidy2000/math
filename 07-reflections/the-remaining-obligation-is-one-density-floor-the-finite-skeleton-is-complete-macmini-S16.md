# The remaining obligation is one density floor — the finite skeleton is complete

*mac-mini-2026-07-06-S16b (HYP-4432). Owner: work the remaining obligation
creatively, pull from the fleet many times. This session wove together the
fleet's fresh bricks (opus coverer_height, kps gap_candidate/residue-split) with
my witness-denominator lever and multi-scale collapse into a single picture: the
combinatorial skeleton of (G) is complete, and exactly one analytic estimate — the
density floor — remains. Verified: `lrc_fastM_highscale_probe…out`,
`lrc_targeted_gap_search…out`.*

## The skeleton, assembled from four agents' bricks

A gap member would be a covering primitive 12-family with `M ∈ (1/13, 2/25)`. Every
lever built this week constrains it, and they compose:

1. **Divisibility-rich** (kps HYP-4417, `gap_candidate_has_multiple`, GREEN): it
   contains a multiple of *every* `k ∈ {2,…,12}` — the covering-system structure;
   the AP is the minimal such family.
2. **Coverer dichotomy** (opus HYP-4406, `coverer_height`, GREEN): for the
   unique-coverer moduli `r ∈ {7,…,12}`, a pinned coverer is `= r` (unlifted) or
   `≥ 14r ≥ 98`.
3. **Witness denominator** (mac-mini HYP-4432, elementary): `M = c/q ⇒ q ∣ v_i±v_j`
   or `2v_i`, so `q ≤ 2·max(v_i)`.
4. **Multi-scale collapse** (mac-mini HYP-4402): a scale gap ⇒ two parts each of
   size ≤ 11 ⇒ each `M ≥ 1/12 > 2/25` ⇒ safe by decorrelation ⇒ not a gap member.
   **Gap members are single clusters.**

Together: a gap member is a **single-cluster, divisibility-rich, primitive** family
whose witness denominator is **bounded by its height**. So the whole problem is a
finite system *once the height is bounded* — and (1)–(4) say precisely where to
look. This is the "tight finite system" kps flagged; I ran it.

## The searches close every bounded regime

Using my witness-denominator lemma as a **fast exact-`M`** engine (the witness `q`
is among the `O(n²)` pairwise sums/differences, so exact `M` in `O(n²·max)` with no
profile-solver blowup):

- **Bounded-height near-AP residual** (`…targeted_gap_search…out`): all one-, two-,
  and three-element covering-preserving perturbations of `{1,…,12}` (singles to 200,
  doubles to 70, triples to 30) — **59,633 covering primitive families, zero in the
  gap.** `2/25` is *attained* (the doubled-apex `{1..11,24}`) as a hard barrier;
  9,042 families sit in `[2/25, 3/25)`, none below.
- **High-scale near-AP** (`…highscale_probe…out`): the dilated AP `N·{1,…,12}` is an
  **isolated tight point at every scale**. A *generic* integer perturbation
  (`|ε|≤3`) makes it **loose** (`M ≥ 1/8`, all 400/400) — a single integer step
  destroys the exact roots-of-unity resonance; a *structure-preserving* ε gives a
  scale-stable **rung** (`+1` on runner 12 → `1/12` at `N=5,20,100`; `+1,+2` → `1/11`;
  …), always `≥ 2/25`. **Zero in the gap.** (Honest: this is isolation, not exact
  scale-invariance — the values drift slightly with `N` but stay loose/rung.)

So every regime a gap member could live in — multi-scale (rigorous), bounded-height
near-AP (59k clean), high-scale near-AP (isolated) — is accounted for. The
empirical (G) is as tight as computation reaches.

## What actually remains: one density floor

The one thing none of this *proves* is the floor itself — that the safe measure at
`2/25` is bounded strictly below its would-be zero for every non-AP covering family,
uniformly in height. In the language of the week's three faces:

- **Additive/three-gap (mine):** near-tight ⇒ the witness is `{kα}` (small `g`) ⇒
  `M` is a CF rung ⇒ no value inside the Farey cell. Needs the converse-three-gap
  rigidity.
- **Multiplicative (kps):** the tight locus is the `(ℤ/13)*`-orbit of roots of unity,
  an isolated discrepancy minimum. Needs the quantitative isolation.
- **Sum–product (opus):** `safe(S,β)` is a theta-sum over the relation lattice `L(S)`,
  maximized at the AP. Needs the **Riesz-product all-order estimate** — the genuine
  analytic core (opus-S106's renormalization contraction rate).

These are one statement: **the AP uniquely zeroes the safe measure at `2/25`, with a
gap.** Every finite/combinatorial reduction is now in place (four green Lean bricks +
the multi-scale theorem + the fast-`M` searches); the residue is a single
uniform-in-height analytic floor. The height bound and the density floor are the same
open fact, and by my lemma each is equivalent to the other.

## Net

- The (G) skeleton is **complete and largely formal**: divisibility-rich +
  coverer-height + witness-denominator + multi-scale collapse. A gap member is a
  bounded, single-cluster, divisibility-rich finite object — and the bounded
  searches find none (59k near-AP; 400 high-scale; the exhaustive small-`n`).
- The **sole remaining obligation** is the density floor (= the height bound = the
  contraction rate = the Riesz-product estimate), the same fact in four languages.
- Engineering byproduct: **fast exact `M`** via the witness-denominator lemma —
  `O(n²·max)`, no solver blowup — now the fleet's workhorse for gap searches.

## Pointers

- `lrc_fastM_highscale_probe_macmini_S16.py` (the fast-`M` engine + high-scale probe),
  `lrc_targeted_gap_search_macmini_S16.py`-style search (`…out`).
- kps HYP-4417 (`LRCGapCandidate.lean`), opus HYP-4406 (`LRCCovererDichotomy.lean`),
  mac-mini HYP-4432 (witness-denominator), HYP-4402 (multi-scale), HYP-4412
  (three-gap), opus-S106 (renormalization flow / floor).
