---
source: oracle-2026-06-01-S542
status: synthesis + verified computation (the p-adic tree unifies channels, apex, moat, center; channel rank = omega(n/2))
tags: [lonely-runner, p-adic, bruhat-tits, tree, building, channels, apex, affine-symmetric-group, annular-braid, CRT, rank]
---

# The p-adic tree as the common home: channels are its branches, the apex is the root-branch, the moat is the boundary

**Prompt (user):** continue the themes, integrating the p-adic tree.

The Bruhat-Tits / p-adic descent already lives in the LRC program (codex S399/S400/
S410: forbidden endpoints `e=(nm±1)/(nv)` are rational, and their reduced-denominator
valuations `d_p=v_p(·)` are boundary points of a **product of Bruhat-Tits trees**;
"gate repair = debt export to child vertices"). This session shows that the objects
my recent thread built — the **channels** (S532), the **apex** (S530), the **annular
center = shift** (S541), the **covering character sum** (S526) — are all features of
that *same* p-adic tree. One picture holds them.

## The tree, and what each LRC object is on it

For each prime `p | n`, the speeds and times have a `p`-adic tree (the rooted
`p`-ary tree `Z_p`; the SL_2(Q_p) Bruhat-Tits tree). The LRC system sits on the
**product** of these trees over `p | (n/2)` (channels) and `p | n` (moat).

| LRC object | on the p-adic tree |
|---|---|
| **channels** (S532) = residues mod `n/2` | the **vertices at depth `v_p(n/2)`** of the product tree (CRT); count `= n/2` |
| **apex / source-sink** (S530), speed `n/2` | the **residue-0 branch** (toward the root/gate) — the singleton channel |
| **channel coupling** (S531/S532, why `n≥6` is hard) | the **tree depth + rank**: the channels split only at deeper levels, and several prime-towers couple |
| **endpoint moat** (S410) | the **boundary** at depth `d_p = v_p(n)`; "debt export" = descent to child vertices |
| **annular center = shift** (S541), `v→v+c` | **translation toward the boundary** of the tree; LRC is therefore a **boundary condition** (basepoint-free) |
| **covering character sum** (S526) | **harmonic analysis on the tree** (spherical functions / the residue-field characters) |

Verified (`lrc_padic_tree_channels_s542.py`, cleanest at **n=14**): the 7-adic tree's
level-1 has 7 nodes = the 7 channels; **residue 0 = `[7]` = the apex** (speed `n/2`);
residues 1..6 = the 6 live pairs `{i,i+7}`; and **6/7 of those nodes split at level 2**
(the coupling is the depth). The moat sits at `{2:1, 7:1}`.

## The new invariant: the channel RANK is `ω(n/2)`

Channels `=` residues mod `n/2` `=` the product `p`-adic tree at depth `v_p(n/2)`.
So the **number of independent prime-towers is `ω(n/2)`** (the count of distinct
primes of `n/2`):

```
 n     n/2   factor   rank ω(n/2)   tree
  4     2    2          1           rank-one, depth 1
  6     3    3          1           rank-one, depth 1
  8     4    2^2        1           rank-one, depth 2
 10     5    5          1           rank-one, depth 1
 14     7    7          1           rank-one, depth 1   <- simplest non-trivial
 16     8    2^3        1           rank-one, depth 3
 18     9    3^2        1           rank-one, depth 2
 12     6    2·3        2           rank-two (product of two trees)
 20    10    2·5        2           rank-two
 30    15    3·5        2           rank-two
```

So **`n = 2p^k` are the rank-one cases** (a single `p`-adic tower); among the
genuinely-hard `n` (beyond the proved small range), **`n=14` is the simplest** —
one prime, depth one. The composite-`n/2` cases (`n=12,20,30,…`) are **products of
trees**, where the difficulty is the *coupling between prime-towers* (the rank-`≥2`
form of the S532 channel coupling). This is a sharper "why some `n` are harder":
not the size of `n`, but the **prime-power shape of `n/2`** — its position in the
tree.

## The deep reading: LRC is a recursive moat on a boundary

Put the three frames together:
- the **center = shift** (S541) means LRC only sees the **boundary direction** of
  the tree (basepoint-free) — it is a condition at infinity;
- the **moat** (S410) is the obstruction at each boundary layer, **exported to
  child vertices** when a gate (a protector speed) clears a layer;
- the **channels** (S532) are the first branch level, the **apex** the root-branch,
  and the **coupling** the descent into the tree.

So **LRC at composite `n` is a recursive descent on the product `p`-adic tree**:
clear the moat at the channel level (depth `v_p(n/2)`), and the obstruction either
vanishes or is pushed to deeper vertices (the "debt"); the conjecture is that the
recursion always leaves a clear boundary ray (a fat collar / a lonely time). The
annular braid (S541) provides the group acting: the affine Weyl group `S̃_{n-1}` is
the **apartment** of the building, its alcove walk (S525) a path toward the
boundary, and the center the translation along it.

## Open (→ HYP)

- Is the rank-one structure (`n=2p^k`) enough to make LRC provable for those `n` by
  a single `p`-adic tower argument (the way `n=14`'s `n/2=7` gives one clean 7-adic
  tree)? `n=14` would be the first test: a rank-one, depth-one moat descent.
- Does the rank `ω(n/2)` predict the *order* of the irreducible coupling — i.e., is
  the S526 covering inside-debt a sum over the `ω(n/2)` prime-towers, vanishing
  exactly when the towers decouple?

## Anchor
`04-computation/lrc_padic_tree_channels_s542.py` (+ `.out`): channels = p-adic tree
level (n=14: 7 nodes, apex=residue-0=[7], coupling=depth); rank `ω(n/2)` table.
Builds on S410 (Bruhat-Tits descent), S532 (channels), S530 (apex), S541 (annular
center), S526 (covering), S525 (alcove walk).
