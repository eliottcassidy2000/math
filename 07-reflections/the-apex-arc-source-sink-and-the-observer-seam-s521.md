# The apex arc: the source–sink tile closes the polygon, gates the menu, and seats the observer (S521)

*claudebox-2026-06-01-S521. Prompted by the observation that the "outside" of the
tournament polyhedron is the base path plus the single source–sink arc, and that
this arc occupies an important place in the tiling model. It does — and it is the
exact coordinate that carries the LRC content.*

## The boundary of the polygon = base path + source–sink arc

Fix the base Hamiltonian path `n -> n-1 -> ... -> 1` (the `n-1` base-path arcs,
the cut space / score hierarchy). Adding the single arc `(n,1)` between the source
`n` and the sink `1` closes the path into a **Hamiltonian cycle** — the boundary
of the polygon on the `n` vertices. So:

> outside of the polyhedron  =  base path (`n-1` arcs)  +  source–sink arc `(n,1)`.

That closing arc is the **apex tile**: the unique longest-range tile (`x-y = n-1`),
the corner of the staircase `δ_{n-2}`, and the **maximal element** of the
interval-containment (type-`A_{n-2}` root) poset on tiles.

## Theorem-let: the apex flips exactly when the tournament leaves the ranking

In the arc-confined half-turn menu (THM-387), a tournament is the transitive
backbone with a flipped **up-set** `S` of long pairs. Since `(1,m)` (the
source–sink / apex pair, `m=n-1`) is the unique **maximum** of the
interval-containment poset, any nonempty up-set contains it. Hence (verified
m=3..7):

> **apex tile `(1,m)` flipped  ⇔  `S ≠ ∅`  ⇔  the tournament is non-transitive.**

The apex is the *single* tile whose state distinguishes "trivial ranking"
(`H=1`, transitive) from "has cyclic structure." Of the `2^{m-1}` feasible
flip-sets, exactly `2^{m-1}-1` have the apex flipped — all the non-transitive
ones.

## The apex gates the menu (the `L>1/2` switch)

The apex pair has the largest separation, so it is the **first** tile to flip as
the runner spread `L` grows, flipping exactly when `L > 1/2`:

- `n=4` (`L = 1/2`): the apex is frozen forward → only the transitive class →
  menu `= 1`.
- `n ≥ 5` (`L > 1/2`): the apex can flip → the full `A000016(n-1)` menu opens.

So the "menu is constant on `(1/2,1)` and switches on at `1/2`" phenomenon
(HYP-1993) **is** the apex switch. The apex tile is the gate between the trivial
regime and the round-tournament necklace.

## The apex seats the observer — where LRC lives

Closing the runners' linear ranking into the circle creates a **seam**: the arc
between the last and first runner (the source–sink arc). In the LRC regular-polygon
picture the observer sits exactly in this seam — the polygon "opens" at the apex
to seat vertex `0`. Consequently THM-384's loneliness criterion (both
observer-adjacent gaps `≥ 1/n`) is a condition **localized at the apex arc**: the
two gaps flanking the observer are the two halves of the source–sink seam.

This dovetails with the marked-tournament reduction
(lrc-as-a-marked-tournament-on-n-vertices-s521): the runner sub-tournament is
**observer-blind** (the `A000016` round/necklace base, built from the base path
and the interior tiles), while the **apex/seam carries all the LRC content** (the
observer's two adjacent gaps). The user's "one very important arc" is precisely the
LRC-bearing coordinate — the place where the additive walk must open a `2/n`
window.

## Why the project already knows the apex is special

- The apex/hypotenuse is where `H = 1 + 2^d` lives (the "blue line from source");
  flipping the apex creates the longest cycle through source and sink.
- The recursive staircase decomposition names the **apex** (arc between the
  extreme vertices) as one of its four pieces (overlap / bottom wiring / top
  wiring / apex).
- Mode-A reduction (`n -> n-1`, hypotenuse removal) peels the apex region.

## Synthesis

Three roles of the one source–sink arc, now unified:
1. **topological** — it closes the base path into the polygon boundary;
2. **combinatorial** — it is the poset-maximum apex tile; its flip is the exact
   indicator of non-transitivity and the gate of the `A000016` menu;
3. **dynamical (LRC)** — it is the observer's seam, where loneliness is decided.

The base path is the observer-blind ranking; the interior tiles are the
observer-blind round-tournament structure (`A000016`); the apex is the single
coordinate where the observer — and therefore the entire Lonely Runner problem —
lives. That is why flipping that one arc is "less independent than it looks": it
is the hinge between the ranking, the cycle space, and the observer.

## Seed

Reformulate LRC(n) as a statement about the **apex coordinate alone**: over the
observer-blind base (`A000016` round tournaments, set by the base path + interior
tiles), the apex seam must be openable to a `2/n` window by some rotation. This is
the "fiber-only LRC" of the marked model, now with the fiber identified as the
apex tile of the staircase — connecting the LRC observer directly to the project's
hypotenuse / `H = 1 + 2^d` machinery.
