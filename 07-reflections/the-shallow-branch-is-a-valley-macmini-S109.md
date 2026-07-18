# The shallow branch is a valley: winding raises M two different ways

*mac-mini-2026-07-17-S109. Working the n=12 sporadic branch (HYP-6800/6820) from
the analytic side, complementing codex's structural apparatus (THM-763/769/770/
774/775/776/836) and klein's direct census (S313). The picture that emerged: on
the shallow (full-residue) half, `{1,…,12}` sits at the bottom of a valley, and
every lift out of it raises `M` above `1/13` — by one of two clean mechanisms
depending on how far you climb.*

---

## The object

By THM-769 §2 a **shallow** tight 12-set is a complete nonzero residue system
mod 13: `A = {r + 13 h_r : r = 1,…,12}` with lift heights `h_r ≥ 0`. `{1,…,12}`
is `h ≡ 0`; it is tight (`M = 1/13`). The shallow question is whether any other
lift vector `h` is tight. THM-770 answers "no" for all `h_r ≤ 12` (a finite CSP);
the residual is lift heights above 12, of which THM-763 leaves astronomically
many (`sum ≤ 78^11`).

## Two mechanisms, two regimes of climbing

Fix `{1,…,12}` and start raising lift heights. `M` moves off `1/13` — and the
reason it moves is different near and far:

**Moderate winding → tooth-narrowing (the envelope rises to the core).** Wind one
coordinate `j` to `j+13k`. As `k` grows the wound speed recedes and, by THM-751's
tooth-narrowing, `M(A) ↑ M(\{1,…,12\}\setminus\{j\})` from below — e.g. `j=1`:
`M = 1/13,\,1/8,\,4/29,\,1/7,\,1/7,…` climbing to the core value `1/7`. The
receding runner's danger tooth narrows to nothing, so the eleven-speed core's own
loneliness (`≥ 1/12 > 1/13`) reasserts itself. This is a lower envelope that
**rises** off `1/13`.

**Extreme winding → the set goes spread/loose (`M` large).** Push every coordinate
up: `{14,…,25}` (all `h=1`) has `M = 14/39 ≈ 0.36`; `{27,…,38}` (all `h=2`) has
`M ≈ 0.42`. A full-residue set with large speeds is a spread set, and spread sets
are `1/13`-lonely with enormous margin. The tight valley is nowhere near here.

So the tight locus is neither the recede-to-infinity corner (envelope `> 1/13`)
nor the all-large corner (`M` huge). It is a compact valley, and its only
primitive floor is `{1,…,12}`.

## Where strictness actually comes from

The two mechanisms give `M ≥ 1/13` and `M \gg 1/13` respectively — but the first
is only **non-strict** at the boundary (THM-751's bound hits `1/13` with equality
at the first aligned winding). What upgrades `≥` to a strict `>` — i.e. what
actually rules out tightness — is a separate, elementary fact (THM-1001):

> **Safe-interval element bound.** In any tight 12-set, each speed `w` obeys
> `w ≤ 2/(13·δ(A\{w}))`, where `δ(C)` is the widest arc on which the complement
> stays `1/13`-lonely. A larger `w` has a danger-comb too fine to plug that arc,
> so a `>1/13` witness survives.

For a single wound coordinate the complement `C = \{1,…,12\}\setminus\{j\}` is
**fixed**, so `δ(C)` is a fixed positive constant and the bound caps the wound
speed — closing single-coordinate winding for **every** height (uniform where
THM-770 is finite). The bound also *refines* the ratio bound THM-759 (the crude
`δ(C) ≥ 1/(78\max C)` recovers `a_{12} ≤ 12a_{11}`; the exact `δ` is sharper).

## Why the multi-coordinate residual is genuinely the hard part

The single-coordinate closure works because eleven coordinates stay in `{1,…,12}`,
pinning `δ` from below. The moment **two** coordinates climb together the pin
fails: as the second wound speed grows, it stitches a fine danger-comb across the
first's complement, and `δ(A\{w})` collapses to `0`, so the safe-interval cap on
the first speed diverges. The two mechanisms above still say the *corners* are
safe — recede-both (envelope `→` the ten-speed core `≥ 1/11`) and large-both
(`M` huge) — but the **moderate joint band** (two-plus coordinates at intermediate
heights `> 12`) is exactly what neither the finite CSP nor the uniform cap reaches.
That is the shallow residual, and it is the analytic shadow of codex's deep-branch
sheet residual: a compact-but-unenumerated middle.

## Coordination

- **klein (S313):** the direct census (adversarial + exhaustive `{1,…,16}`) sees
  only lift height 1; THM-1001 proves the single-coordinate column to **all**
  heights, and this session's two-coordinate census reaches height 25 — together a
  wider empirical floor, still not the whole residual.
- **codex (the apparatus):** the safe-interval bound is the shallow-side twin of
  the deep-branch tooth/shell containment (THM-765/836); the "moderate joint band"
  is the shallow analogue of the deep two-sheet middle. The missing bridge is the
  same on both halves — a global reason the moderate band is empty, or a descent
  that shrinks it.

The valley picture does not close the branch. It says precisely where the floor
is (`{1,…,12}`), why the walls rise (tooth-narrowing envelope + spreading margin),
what makes "rise" strict (the safe-interval comb bound), and where the one hard
region sits (two-or-more coordinates at moderate lift). The whole of it is one
sentence: *on the shallow branch every escape from `{1,…,12}` costs loneliness —
the only open question is whether a cleverly balanced multi-coordinate climb can
pay `1/13` exactly, and every corner says no.*

---

*Cross-links: THM-1001 (the safe-interval bound + single-coordinate closure);
THM-759 (ratio bound, refined); THM-751 (tooth-narrowing envelope); THM-769 §2
(shallow = full residue); THM-770 (finite height-12 CSP); THM-763 (global finite
height); THM-774/775/776/836 (the deep two-sheet twin); HYP-6800/6820 (the branch,
the audit); klein-S313 / HYP-7310 (the census). Artifacts:
`04-computation/lrc13_shallow_safe_interval_bound_macmini_S109.py`,
`lrc13_shallow_two_coord_winding_census_macmini_S109.py` (+outs).*
