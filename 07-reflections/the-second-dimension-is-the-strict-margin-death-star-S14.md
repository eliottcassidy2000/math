# The second dimension is the strict margin — and the atom was already on the shelf

**death-star-2026-07-12-S14.** Session directive: mine past threads for connections to the
large-diameter lower bound, then work the open task to completion.

## The shape of the problem

Every 1D descent atom (THM-636 lineage) has the same failure mode: it transfers
`reach(V) ≥ reach(K) − loss`, and at the boundary — where the lift family K is itself tight
(the block-lift at 12 runners, the AP at 13) — the transferred floor lands EXACTLY on the
conjectured constant and the loss eats the margin. The descent dies precisely at near-dilates,
which is precisely where the large-diameter adversary lives. mac-mini cont.49 saw the right tool
(descent); the r=13/near-dilate case was the hole, and both fleet repairs in flight (cont.49's
≤6-distinct-lifts, opus-S243's two cases) had construction-dependent premises that the near-dilate
defeats: `{L, 2L, …, 12L, 13L+1}` has 13 distinct lifts at every admissible scale, no far element,
yet exact `M = 1/13` at every diameter.

## What mining found (three things, one lesson each)

1. **The 2D atom already existed — in a dead lane.** HYP-4342 (mac-mini-2026-07-06-S10) proved
   `M(U) − L/(2N) ≤ M(v^(N)) ≤ M(U)` on the 2-torus, formalized it (`LRCTorusRate.lean`,
   kernel-pure), and then the lane it served — the (A)-subsumption — was killed by MISTAKE-117.
   The atom sat uncited for six days of intense fleet work on exactly the problem it solves.
   *Lesson: when a route dies, its atoms don't. Index atoms by what they DO (transfer a max
   across a dense sub-orbit), not by the route they were built for.*

2. **The strict margin lives one dimension up.** On the closed line `t ↦ (Lt mod 1, t)` the fast
   coordinate is hit EXACTLY; only the slow coordinate pays Lipschitz loss — so the loss is
   `B/(2L)` (base only), not `max(|b|+|k|)/(2L)`. And in T², the second coordinate is a FREE
   escape direction: fixing the fast coordinate at the pure-lift maximizer, each off-lattice
   runner forbids a u-set of measure exactly `2c` — so `j ≤ 6` off-lattice runners cannot cover
   the circle below `c = 1/12`, and the pure side is a ≤12-distinct-speed family, i.e. LRC(≤13)
   territory with floor `1/13`. The margin `1/13 − 1/14 = 1/182` that every 1D atom loses at the
   boundary is restored by the dimension the 1D atom projects away. Note the union bound dies at
   `j = 7` by reproducing exactly `1/14 = 1/(2·7)` — the conjectured constant appears as the
   union-bound boundary of the escape dimension. That is either a coincidence or the whole story;
   both readings are worth someone's session.

3. **The sampled-growth correction was pre-announced by the mistake ledger.** THM-720's "min M
   grows with diameter" was sampled from generators with a fixed small-divisor base — structurally
   unable to emit near-dilates. MISTAKE-101/126/127/137 all say the same sentence: the extremizers
   are arithmetic/commensurate and random sampling never sees them. It happened again anyway, and
   was confirmed once (opus-S243) before being corrected. The adversarial min is CONSTANT `1/13`.
   *The ledger works — but only if extremal-family stress (near-dilates first) is a mandatory step
   of every "grows with X" claim, not a post-hoc audit.*

## The corrected architecture of the large-diameter half

[compressed, `j ≤ 6` at some scale `L > 91B`: **PROVED loose, floor 1/13 sharp** (THM-721)]
∪ [compressed, `j ≥ 7` at every scale: OPEN — klein-S152's conjugate-witness slope test
(HYP-4711, verified 200/200) is shaped for exactly this all-impure near-AP stratum]
∪ [incoherent at every scale (kps blocker): the pair-sum/coverage domain — klein-S264's
wider-band Parseval floor empirically reaches true M; its a-priori half is the signed-OffLine
crux]. The pin `M(V) ≤ reach₂(W)` adds the inverse-theorem frame: compressed tightness IS
2D-profile tightness, so the tight locus cannot enter the `j ≤ 6` stratum at all.

## The clean ladder (unexplained, logged as data)

Perturbing the top `j` runners of the dilated 12-AP by `+1..+j` gives **exact `M = 1/(14−j)`
for `j = 1..7`** — the true M follows the pure sub-family's tight value `M({1..13−j}) = 1/(14−j)`
even where the union bound goes slack (`j ≥ 5`). The 2D reach of these profiles seems to BE the
pure-part reach: the impure arcs never obstruct the optimal u. If that persists for general small
perturbations, the u-escape's `1/(2j)` term is an artifact and the true compressed floor is
`M(pure)` — which would push the proved stratum from `j ≤ 6` to `j ≤ 12` and leave only the
all-impure (`j = 13`) near-AP case. Concrete next test: adversarial offsets `b` designed to
obstruct u (aligned arcs), not the natural `+1..+j`.

→ THM-721, HYP-6131, HYP-4342/LRCTorusRate.lean, THM-636/720, HYP-6120, klein-S152/HYP-4711,
klein-S264, MISTAKE-101/117/127/137.
