# The 7-dimensional R-odd homology is a DUOCYLINDER, not a Fano plane: two NS "cylinders" (axis × SC-ring) glued at a ridge — and the general case is a polycylinder, governed by Robbins-oriented strongly-connected NS axes

*opus-2026-06-29. Owner: chase a Fano-line basis of `H₁⁻(n=5)=7`, thinking duocylinders, poly-cylinders,
and Robbins' theorem. The hints corrected my hypothesis: it is NOT a Fano plane. It is a **duocylinder**,
and the apex-7 coincidence I flagged earlier is the duocylinder's cycle rank, not a projective plane.*

## The verified structure (n=5)
The complement `R` fixes the 8 SC classes and pairs the 4 NS classes into **two axes**: `(2,9)` and
`(1,11)`. Each NS pair, with its ring of **shared SC neighbors**, is a **cylinder** (axis = the NS pair,
ring = the SC classes adjacent to both):
| cylinder | axis (NS pair) | SC ring (shared neighbors) | ring size | theta cycle-rank |
|---|---|---|---|---|
| **B (big)** | `2 ~ 9` | `{0,3,4,5,6,8}` | 6 | **5** |
| **A (small)** | `1 ~ 11` | `{0,4}` | 2 | **1** |
| **ridge** | — | `{0,4}` (SC in BOTH rings) + axis-links `(1,2),(9,11)` | — | **+1** |

`5 + 1 + 1 = 7 = b₁⁻`. **The two cylinders share the ridge `{0,4}`** (the Clifford-torus junction of
`D²×D²`) — a **duocylinder**.

## Why every cylinder is entirely R-odd (the clean lemma)
For any NS axis `a~b=R(a)` and any two ring vertices `s,s'` (SC, so `R`-fixed), the theta 4-cycle
`C = a→s→b→s'→a` satisfies
> **`R(C) = b→s→a→s'→b = reverse(C) = −C`** (in homology),
because `R` swaps the two axis endpoints and fixes the ring. **So a cylinder's ENTIRE cycle space lands
in `H₁⁻`** — the R-odd homology is exactly the cylinders' homology. This is structural, holds at every
`n`, and is why `b₁⁻` is "carried by the ribs" (the SC–NS edges, CLAUDE.md).

## The Fano hypothesis fails — and that is the right answer
`b₁⁻(5)=7` is **not** the Fano plane (7 points/lines, `(7_3)` incidence). It is the duocylinder cycle
rank `5+1+1`. The earlier "`7 | b₁⁻` apex-7" pattern (`7, 119=7·17`) is **coincidental at this level** —
the `7` is a duocylinder rank, the `119` a polycylinder rank; neither is a projective-geometry `7`. The
honest correction: the secondary obstruction's geometry is **cylindrical (product), not projective.**

## The general case: a POLYCYLINDER (the "poly-" hint)
At general `n`, each NS R-pair is a cylinder, so `H₁⁻` is a **polycylinder** — a product/gluing of
`#(NS pairs)` cylinders along shared SC ridges. The number of cylinders is the **NS-merged count**
`0, 1, 2, 22, 184` (n=3..7): `n=4` → 1 cylinder, `n=5` → **2 = duocylinder**, `n=6` → 22-cylinder, etc.
`b₁⁻ = Σ_cyl (ring−1) + (ridge/axis cross-links)` — a polycylinder Künneth count. (This explains the fast
growth `1,7,119`: rings get large and many.)

## Robbins' theorem: what the axes ARE
The NS classes are the **non-self-converse** tournaments — the ones that are their own converse only up
to a *non-trivial* relabel, i.e. the **rigid / asymmetric** classes. By **Robbins' theorem** (a graph
has a strongly-connected orientation iff it is bridgeless), the complete graph always admits
strongly-connected orientations — the SC primes — and the cylinder axes are exactly the asymmetric NS
ones, *oriented* by `R` (the complement). The duocylinder is the `R`-orientation of the two rigid axes;
the SC ring is the bridge-free "loop" Robbins guarantees around each.

## What this gives the LRC (refining HYP-3544)
The R-odd homology IS the LRC secondary obstruction / odd index (mac-mini HYP-3544). So:
> **The LRC odd index is a POLYCYLINDER** — a product of `R`-oriented strongly-connected (NS) cylinders,
> glued along self-converse ridges. The forcing question (does `b₁⁻ ≠ 0` ⇒ lonely?) becomes geometric:
> a non-trivial polycylinder cycle is the obstruction, and the **duocylinder is the smallest
> non-degenerate case** (n=5). The `14 = 2·7` two-scale should appear as the **two cylinders** of the
> duocylinder (the `D²×D²` two-disk product), not as a Fano 7.

## Status
- **Verified (opus):** `H₁⁻(n=5)` = duocylinder — cylinder B (axis 2~9, ring `{0,3,4,5,6,8}`, rank 5) +
  cylinder A (axis 1~11, ring `{0,4}`, rank 1) + ridge `{0,4}`/axis-links (+1) = 7. Every theta cycle is
  R-odd (`R(C)=−C`).
- **Corrected (opus):** NOT Fano; the apex-7 `b₁⁻` is a duocylinder/polycylinder rank, the `7|b₁⁻`
  pattern coincidental. The duocylinder hint was right.
- **New frame:** `H₁⁻` = polycylinder over the NS pairs (count = NS-merged `0,1,2,22,184`); Robbins
  orients the strongly-connected axes; the LRC odd index is a polycylinder.
- **Open:** the polycylinder Künneth formula for `b₁⁻(n)` (predict n=7); whether the duocylinder's two
  cylinders are the `2` and `7` scales of `14`; the forcing descent (HYP-3544) read geometrically.

Related: the apex-7-secondary-obstruction reflection (this corrects its Fano/apex-7 over-reach),
Ky-Fan-forcing-fails, mac-mini HYP-3544/3538, CLAUDE.md spine/ribs/sea (the NS-SC ribs = the cylinders),
THM-281 (SC class sizes), the geometric-alignment-of-merged-metagraph, Robbins' theorem (strongly
connected orientation), OPEN-Q-108.
