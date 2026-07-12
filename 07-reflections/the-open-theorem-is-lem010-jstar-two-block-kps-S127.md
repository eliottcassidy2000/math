# The single open theorem is LEM-010's j* = O(k) — and my extremal is its two-block hard case

*kind-pasteur-2026-07-11-S127. Owner: "work the single remaining open theorem; search back across random past
threads for inspiration." I searched, and the theorem I have been circling for a dozen sub-sessions turned out
to already have a name, a reduction, and a proved special case — in LEM-010, from months ago.*

---

## Two threads, one theorem

I spent cont.28–42 building a route to LRC(14) through the discrete Bonferroni `B5` / clean rulers / the
bounded-clearing window, and pinned the hard core to: *every divisor-complete family clears at a bounded
non-14 modulus.* opus-S235 made the margin a free corollary of that clearing. So the whole hard core is one
statement: **bounded-clearing for divisor-complete families.**

Searching the past threads, LEM-010 (the deterministic good-period existence, THM-527-A) is the *same*
statement in the covering-reduction framing, and it is much further along than I realized:

- It reduces the covering case to a **good period** `j` (a clearing multiplier) via two elementary lemmas —
  the **wraparound** (`spread < 6·Vmax/7 ⟹ j=1`) and **Dirichlet** (`Vmax > 3^{k-1} ⟹ some j ≤ 3^{k-1}`).
- The sharp form is `j* = O(k)`, and LEM-010(iii) calls it *"the single remaining analytic lemma of the whole
  covering case."*
- The **AP case is PROVED** (mac-mini-S59): for `E = {0,d,…,(k-1)d}`, `j* ≤ ⌈7(k-1)/6⌉` by Dirichlet — one
  `j` makes `‖jd/Vmax‖` tiny, so the phases are a `k`-term AP of step `< 6/(7(k-1))`, spanning `< 6/7`,
  leaving a `> 1/7` gap.
- General `j* = O(k)` has zero counterexamples over 90k+ adversarial clusters.

So my "bounded-clearing for divisor-complete" and LEM-010's "`j* = O(k)`" are the same open theorem, reached
from opposite ends — the discrete B5/margin side and the good-period/covering side.

## And my extremal is exactly its hard case — the ~1% corner, not the bulk

LEM-010 names where the reduction is not yet a proof: *"structured **2-block**/AP clusters, where the soft
discrepancy bound is vacuous."* My M-floor extremal from cont.41 is
`{1,2,3,4} ∪ {10,…,18}` — a **two-block** family (verified: blocks = 2), clearing at `q=24` with multiplier
`p=5 ≈ k`. The worst-clearing family `{1,2,3} ∪ {10,…,19}` is also two-block (`q=25, p=6`). The shift-AP
`{2,…,14}` is one block (an AP) and clears at `p=1` — exactly LEM-010's wraparound lemma.

**But this corner is small — and opus-S237 (same day) pins how small.** The divisor-complete class is
**99–100% spread** (longest-AP ≤ 7): divisor-completeness *forces* spread, because it requires multiples of
`8, 9, 11, 13, 14` — specific, spread-out speeds incompatible with a tight interval. The near-AP / 2-block
families (longest run ≥ 8, like my extremal's `{10..18}` 9-run) are a **~1% corner**. So the honest picture is
a split, and the two are opposite in every way:

- **The 99% spread bulk** (opus-S237): loose (`M ≥ 3/29 > 1/12`), clears at `q ≤ 29` with *zero* exceptions,
  reduced to *one uniform anti-concentration statement*. This is the same object as my cont.36 decoupling —
  "the window-hard covering cores are loose/spread, not near-AP" — which opus-S237 independently confirms.
- **The ~1% near-AP corner** (this reflection, via LEM-010): carries the *M-floor* `1/12` at the two-block
  extremal, and is exactly where LEM-010's soft discrepancy bound goes vacuous. LEM-010's explicit good-period
  construction — `j = 1` wraparound for one block, Dirichlet for the AP — is precisely the tool built for this
  corner.

So the extremal I found and the bulk opus reduced are the two halves of the same residual, and they need
*different* tools: anti-concentration for the loose spread bulk, explicit good-period for the tight near-AP
corner. My cont.40–42 work — the exact M-floor `1/12`, the extremal `{1,2,3,4,10..18}`, the bounded window
`[15,~27]` — is the sharp quantitative face of the corner: the worst near-AP cluster still clears at
`p = O(k)`, giving margin `≥ 1/12`.

## The proof strategy the search hands us — for the corner

The AP proof is a **one-dimensional Dirichlet**: make the single step `‖jd/Vmax‖` small. A two-block cluster
has *two* directions — the within-block step (`1`) and the between-block gap (`b`). The natural extension is a
**two-dimensional simultaneous Dirichlet**: find one `j ≤ C` with `‖j/Vmax‖` *and* `‖jb/Vmax‖` both small.
Then block 1 sits near a short small-step AP at `0`, block 2 near one at `jb/Vmax ≈ 0`, so the two blocks
overlap into a single `< 6/7` span and leave a `> 1/7` gap. Two-dimensional Dirichlet gives `j = O(k²)` —
weaker than the AP's `O(k)`, but still **bounded**, which is all bounded-clearing needs. Iterating to
`m`-block clusters gives `O(k^m)` with `m` bounded (a 13-element family has `≤ 13` blocks, and the extremal
uses `2`). That is a concrete, elementary path to the general lemma via the tool LEM-010 already established
for one block.

## Why the search mattered

For a dozen sub-sessions I treated the divisor-complete core as a fresh frontier and built new machinery for
it — the residue-periodicity, the dilation-invariance, the tiered covering, the exact `1/12`. All of it was
real, and all of it was re-deriving the quantitative skin of a lemma the fleet had already reduced and
half-proved. The lesson is the plain one: before building, search the corpus for the theorem you are actually
proving. LEM-010 had the reduction, the AP proof, and the exact name of the corner — 2-block — sitting there
the whole time; and opus-S237, working the same day from the anti-concentration side, reduced the 99% spread
bulk to one uniform statement. The right next move is not more machinery. It is two named, bounded pieces:
opus's uniform anti-concentration for the spread bulk, and the two-dimensional Dirichlet that turns LEM-010's
one-block good-period proof into a two-block one for the corner. The residual is no longer a frontier; it is a
punch-list of two.

*Files: this session's verification. Connects LEM-010 (good-period, `j*=O(k)`, AP proved, mac-mini-S59) and
opus-S237 (the 99% spread bulk is one uniform anti-concentration statement; divisor-completeness forces
spread) with HYP-6055/6070 (the M-floor `1/12`, the two-block extremal, the bounded window), opus-S235
(band-edge margin), THM-527 (covering reduction). Extends
[[the-window-is-decoupled-from-tightness-kps-S127]] and reconciles with
[[the-residual-is-99pct-spread-the-AP-corner-was-small-opus-S237]].*
