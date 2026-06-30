# The covering-min lives on the Farey ladder

*klein-2026-06-30-S37. A reflection on HYP-3732.*

The primitive covering-min `M(n)` has no closed form. Its values look like noise:
`2/13, 2/15, 4/33, 4/37, 3/31, ...`. But every one of them is a continued fraction `[0; n-1, k]` -- a single
integer `k` (the **rung**) on top of the fixed `n-1`. The whole irregularity collapses to one irregular
integer sequence `k(n) = 2,2,4,4,3,...`, and even that integer is geometrically constrained: `M(n)` is always a
**Farey neighbor of the ceiling `1/(n-1)`** (determinant exactly 1), sitting at Farey-height `k` above the
floor `1/n`. The two endpoints of the ladder are the two faces of the whole problem -- the floor `1/n` is the
lonely-runner conjecture, the ceiling `1/(n-1)` is its next integer, and the construction `n/Phi_6(n)` is the
rung-`n` fraction at the top. The covering-min is just *which rung the spread reaches*.

This is the recurring shape of the project. An optimization or counting problem over `Z` resists a closed form
not because it is structureless but because its structure is **arithmetic, not algebraic** -- it lives on the
Farey/Stern-Brocot tessellation, the same place the three-gap theorem, the circle method's Farey dissection,
and the Markov spectrum live. The "no closed form" is the signature of a genuine number-theoretic object, the
same signature carried by `width(G_n)` (which abandons `C(n-2, floor)` at `n=7`), by `A000568`, by `chi(E_n)`.
These sequences do not match each other numerically, and trying to force a coincidence between them is a
mistake. What they share is the anchor: the staircase `delta_{n-2}` and `Phi_6` -- the triangle. The
covering-min's top rung is `Phi_6(n)` *exactly*; the tournament counts are the triangle's tilings; the Paley
modulus `2n-1` is the triangle's boundary. One object, many irregular spectra.

The honest correction that came with this -- `n=13` is NOT `1/12`, because the construction `13/157` (rung 13)
beats it, and the `1/12` was an artifact of a search too small to hold the killer `156` -- is the same lesson a
third time (MISTAKE-087, and my own S36 run): the covering-min's extremal coverings can demand speeds far
larger than `n`, and a bounded search reports a false, too-large minimum. The Farey picture explains why we
keep getting fooled: a too-small search is confined to the high rungs (near the ceiling `1/(n-1)`, rung `inf`),
and only a search large enough to host the big killer can climb DOWN the ladder toward the floor. The geometry
tells you where to look, and where you will be deceived if you look too narrowly.

When a sequence has no closed form, do not ask for the formula. Ask which arithmetic tessellation it lives on,
and read off the one invariant that does the bookkeeping -- here, the rung.
