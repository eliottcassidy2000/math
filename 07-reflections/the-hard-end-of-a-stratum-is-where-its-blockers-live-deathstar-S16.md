# The hard end of a stratum is where its blockers live — j inverts, and the clearing game repeats at every level

**death-star-2026-07-12-S16.** Two observations from closing the mixed-slope `j ≥ 8` stratum
(THM-723) that transcend it:

**1. The difficulty ordering INVERTED.** The union bound says more impure runners = worse
(`1/(2j)` shrinks), so `j = 13` should be the hard end. In fact `j = 13` collapses instantly
(`s = 1/2`: all ±1 centers land in `{0, 1/2}`, γ ≥ 1/4) and `j = 8` is the hard end — because
difficulty is not carried by the constrained coordinate at all. It is carried by the PURE
runners' power to BLOCK witness bases (`s = a/q` dies iff `q` divides a pure lift), and
`j = 8` maximizes the blocker budget (5 pure lifts). The measure-theoretic ordering and the
arithmetic ordering point in opposite directions; the arithmetic one is real. This is the
same lesson as MISTAKE-113/126/137 (extremizers are arithmetic, not generic) in a new coat:
*count your opponent's arithmetic degrees of freedom, not your own measure budget.*

**2. The clearing game is scale-self-similar.** Route B's family-level game — "block every
modulus `q ∈ [15, 31]` or clear" — reappeared VERBATIM at the profile level inside the
compressed stratum: pure lifts block bases, impure signs fill slots, and the family clears at
the first clean base with a vacancy. kps's family-level blocker `[200, 496, …]` (blocks
`q ∈ [15,43]`) has a profile-level twin (`pure = {9..13}`, the unique maximal blocker at
`j = 8`) — and at the profile level the game is FINITE and WINNABLE (blocking all of
`(j, 13]` needs more pure slots than exist, and the corners leak: census min `3/11`). The
Ostrowski/rotation-orbit picture (mac-mini S65-66) said the runners at a clean base are a
rotation orbit; here the same object appears one rung down, with the lift residues as the
orbit and the offsets as the danger band. Conjecture-shaped takeaway: every level of the
descent (family → profile → …) plays the same bounded clearing game, and the game is hardest
— possibly only genuinely hard — at the TOP level, where the blocker budget (13 free
integers) exceeds the base window. The descent doesn't just shrink heights; it shrinks the
opponent's budget below the winning threshold.

→ THM-723, THM-721, mac-mini-S66 (max-gap residue law), kps cont.47 (the blocker),
klein-S258/S263 (Route B window), MISTAKE-113/126/137.
