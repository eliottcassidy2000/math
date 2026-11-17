# A sevenfold square that is strictly increasing and never crosses

Replacing the exterior solves the obstruction from the [five-arc commuting square](../chanpak_commuting_six_20261117/README.md). A different fixed exterior makes both vertical paths strictly increasing at every single-arc edit, and it preserves the horizontal order at every step. There is also a two-insertion route to the same count pairs if vertical edges need not be arc reversals.

## A physical six-arc route

Use outgoing masks (318,8,2,20,38,78,159,63,254) on labels 0 through 8; replace label 8 by any nonempty tournament B. Keep its external contacts uniform. In order, reverse the six exterior edges

    (1,6), (2,6), (2,7), (3,6), (5,7), (3,7).

There is exactly one external entrant into the slot: vertex 0. Therefore any full path can visit B at most twice (one initial run and at most one entered run). The two-visit contribution reads the number P₂−2H of nonempty two-slot covers. The entire read is αP₂+βH; substituting transitive blocks of orders one and two determines both coefficients. All seven contracted quotients are strong, which implies every nonempty substitution is strong. The exact replay obtains:

| stage | α | β | lower example | upper example |
|---:|---:|---:|---:|---:|
| 0 | 68 | 5 | 476249 | 476279 |
| 1 | 143 | 19 | 1005639 | 1005753 |
| 2 | 216 | 25 | 1517213 | 1517363 |
| 3 | 366 | 61 | 2579873 | 2580239 |
| 4 | 417 | 57 | 2933301 | 2933643 |
| 5 | 437 | 47 | 3067811 | 3068093 |
| 6 | 476 | 35 | 3333743 | 3333953 |

The last coefficient pair is seven times the first. Each Δα is positive and each 2Δα+Δβ is positive, so P₂≥2H proves strict vertical increase for every possible source B, not only the example. If an internal switch changes H but leaves P₂ invariant, the horizontal difference is β·ΔH. All seven β's are positive. A literal example uses the strongly connected modules (30,204,200,160,142,87,153,33) and (30,204,72,160,142,87,153,37). They differ only on arc (2,7), have equal P₂=6968, and have H=485 and 491. The two lower endpoints 476249 and 476279 are prime. The six vertical reversals and the horizontal reversal are disjoint, and the square with all its intermediate rungs exists on one labeled vertex set.

## Two ways to bypass the five-edge requirement

1. At the original endpoints, allow a single simultaneous two-arc edit: combine (1,8) with (3,4) in the original five-edge route. The read then skips β=−9, moving from (256,5) to (475,109). Every listed state is strong, both tracks strictly increase, and all remaining β are positive. This deliberately makes no claim about the one-arc intermediate.
2. If one-vertex insertions are allowed, use the [existing induced scaffold](../chanpak_eight_split_20261117/LITERAL_LIFT.md) with exterior rows (38,20,32,7,77,26,47), taking its last label as the variable slot. Add a vertex beating only old 0 and 1; then add a vertex beating only the previous newcomer. The original, intermediate, and final reads are (34,5), (117,13), and (238,35), with the same criterion for strict increase and no crossing. All three quotients are strong. For the new prime example, appending one ordinal singleton inside each variable slot doubles its P₂, so this route realizes the same endpoint counts, with middle counts 1636817 and 1636895. This route changes order but does not reverse any old arc.

## The old square remains obstructed with its module contacts held fixed

The old pair at ordinary counts 6103423 and 6103429 does not admit an increasing no-crossing route obtained by arbitrarily changing exterior edges one at a time while keeping its variable module and all of its external contacts fixed, using the same exterior in both columns and insisting on strong states. This statement is finite and specific to the old labeled endpoints. Unlike the previous 120-order check, it permits detours, repeated edits and every exterior arc.

The exhaustive replay uses each 36-bit orientation on the fixed nine exterior vertices as a state. The single exit of the module remains fixed at label 8; the full substitution is strong iff label 8 reaches every exterior vertex (the first time any path from it enters the module it can only leave back at 8). For each state the all-order two-run formula supplies the exact read (α,β); count monotonicity is checked with the actual common P₂=47680 and h=383,389. Both counts must increase and β must stay positive. A forward/backward search starts at the two labeled exteriors; it terminates with an empty forward frontier, not just a depth cutoff. This rules out the attempted class, not a detour using changing module contacts, a batch move, or a different representation of either endpoint.

## Replay

Run `python replay.py` to check the literal masks, all intermediate coefficients and the alternate insertion route. It also checks direct raw counts at stages 0,3,6 and endpoint primality. Run `g++ -std=c++20 -O3 concrete.cpp -o concrete && ./concrete` for the exhaustive old-route search. No result here places a square at every prescribed prime.
