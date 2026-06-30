---
id: HYP-3705
title: ATTACKING n/Phi_6(n) >= 1/n -- the inequality is TRIVIAL (n^2 >= n^2-n+1 <=> n>=1), so the ENTIRE open content is whether n/Phi_6(n) is the genuine COVERING-MIN; reframed: Phi_6(n)=n^2-n+1 is BOTH the EISENSTEIN NORM |n-zeta_6|^2 (Q(sqrt-3)) AND the point-count of PG(2,n-1), so the covering construction is the (Phi_6(n),n,1) SINGER difference set, which is a SYMMETRIC (Steiner) design = the OPTIMAL covering design when it exists (n-1 a prime power); hence 'no covering beats the construction' HOLDS in the combinatorial-design sense, and the open problem localizes to the LRC-floor <-> design-covering BRIDGE. Two-column picture: the Q(sqrt-3)/Eisenstein covering column (difference-set/SPREAD pole, the covering-min n/Phi_6(n), off-cusp positive measure) vs the Q(cos2pi/p) apex column (doublet/CONCENTRATED pole, the gap 4sin^2(pi/2p), the descent attractor)
status: STRUCTURE VERIFIED (Phi_6 = Eisenstein norm, = PG(2,n-1), prime factors =1 mod 6 except 3; n=14 -> (183,14,1) Singer set, q=13 prime). The inequality is trivially true. The reframing (content = combinatorial optimality) and the design-optimality of the projective plane are established; the LRC-floor<->design bridge and M's exact LRC definition are the open/owner pieces.
source: klein-2026-06-29-S23
depends_on:
  - HYP-3700   # mac-mini: the discrete apex-edge isolation (the complementary Q(cos2pi/p) column)
  - HYP-3610   # the multi-axis poles (difference-set=spread, doublet=concentrated)
related:
  - HYP-3604   # the doublet gap 4sin^2(pi/2p) (the apex column's pole)
  - HYP-3611   # the apex core atlas (difference set = the spread pole)
  - HYP-3597   # measure vs existence (the off-cusp positive-measure framing)
results:
  - 04-computation/phi6_covering_min_structure_klein.py
  - 05-knowledge/results/phi6_covering_min_structure_klein.out
---

# HYP-3705 — n/Phi_6(n) >= 1/n: the inequality is trivial; the content is the projective-plane optimality

## The inequality is trivial -- so the content is elsewhere
`n/Phi_6(n) >= 1/n  <=>  n^2 >= n^2 - n + 1  <=>  n >= 1`. The margin is `n/Phi_6(n) - 1/n =
(n-1)/(n.Phi_6(n))` (`= 0.005074` at `n=14`, `14/183` vs `1/14`). So the analytic inequality carries NO
content. The ENTIRE open problem is the claim that `n/Phi_6(n)` is the genuine **covering-min** -- that no
covering set beats the construction.

## The structure (verified): Phi_6 = Eisenstein norm = projective plane = Singer difference set
- **(A) Eisenstein norm.** `Phi_6(n) = n^2 - n + 1 = |n - zeta_6|^2` (`zeta_6 = e^{i pi/3}`), the norm in
  the Eisenstein integers `Z[zeta_6]` -- so the field is `Q(sqrt-3)`. (Verified n=3,7,14.)
- **(B) Projective plane.** `Phi_6(n) = (n-1)^2 + (n-1) + 1 = #points of PG(2, n-1)`.
- **(C) Singer difference set.** When `n-1` is a prime power, PG(2,n-1) gives a `(Phi_6(n), n, 1)` Singer
  difference set: `n` speeds mod `Phi_6(n)`, every nonzero difference exactly once. For `n=14`:
  `Phi_6=183=3.61` (`61 = 1 mod 6` splits in `Z[zeta_6]`; `3` ramifies), `q=13` PRIME, so PG(2,13) EXISTS
  and a `(183,14,1)` difference set of 14 speeds exists.
- **Primes.** Every prime `p | Phi_6(n)`, `p != 3`, is `= 1 mod 6` (splits in the Eisenstein integers); `3`
  ramifies. So `Phi_6(n)` is built from Eisenstein-split primes -- the `Q(sqrt-3)` "existence column."

## The reframing: the open problem is a COMBINATORIAL OPTIMALITY
Because the inequality is trivial, LRC(14) (in this frame) reduces to:
> **the `(183,14,1)` projective-plane covering is OPTIMAL -- no covering of `14` speeds achieves a smaller
> floor than `14/183`.**
A `(v,k,1)` symmetric design (projective plane / Steiner system `S(2,k,v)`) is the **optimal covering
design** when it exists -- the covering number `C(v,k,2)` is attained exactly by the Steiner system (no
redundancy: every pair covered exactly once). So **"no covering beats the construction" is TRUE in the
combinatorial-design sense**. The remaining open content is the **bridge**: does the LRC continuous covering
floor `M` equal the discrete design covering number? That is where M's exact definition (the floor owners')
enters; the structural reduction is solid.

## Two columns, two fields, two poles (synthesis with the multi-axis atlas)
The recent program has two complementary "columns," and the open conjecture binds in the first:
| | field | pole (atlas, HYP-3610/3611) | object | value | regime |
|---|---|---|---|---|---|
| **covering / existence** | `Q(sqrt-3)` (Eisenstein) | the **difference set** (SPREAD pole) | Singer `(Phi_6(n),n,1)` | `n/Phi_6(n)` (covering-min) | OFF-cusp, positive measure, M tightest |
| **apex / measure** | `Q(cos 2pi/p)` | the **doublet** (CONCENTRATED pole) | the odd cycle `C_p` | `4sin^2(pi/2p)` (the gap) | the descent attractor; vanishes at the cusp |
The 2-adic descent flows to the **doublet attractor** (concentrated, apex, `Q(cos2pi/p)`); the **covering-min**
binds at the **difference set** (spread, Eisenstein, `Q(sqrt-3)`). mac-mini-S41/HYP-3700's discrete-edge
isolation is the apex column; this `n/Phi_6(n)` is the covering column. The conjecture binds OFF the
measure-0 cusp, at the positive-measure difference-set construction.

## Net (honest)
- RIGOROUS: the inequality `n/Phi_6(n) >= 1/n` is trivial (`n>=1`); `Phi_6` = Eisenstein norm = PG(2,n-1) =
  Singer difference-set size; the projective plane is the optimal covering DESIGN (Steiner). So the open
  problem is NOT analytic -- it is the combinatorial optimality of the projective-plane covering, and that
  optimality holds at the design level.
- OPEN: the bridge LRC-floor `M` `<->` design covering number (needs M's exact definition); whether the LRC
  continuous covering inherits the design optimality (i.e. `n/Phi_6(n)` is genuinely the LRC covering-min).
- This sharpens the open statement to a finite, design-theoretic claim and locates it in the `Q(sqrt-3)`
  Eisenstein column, off the cusp, at the difference-set/spread pole.
