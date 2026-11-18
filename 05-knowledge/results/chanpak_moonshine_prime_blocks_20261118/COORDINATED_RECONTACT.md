# Two edits can translate an exact +2 window when neither edit can

**All-order:** mixed-difference accounting identity. **Finite-exact:**
all common one- and two-arc edits of the indicated 131/133 pair away
from its two distinguished return arcs. No arbitrary-target selector
or independent composition is proved.

## How the random archive pass suggested the experiment

A deterministic sample (seed 717989, eight files from the sorted
`07-reflections/*.md` list excluding filenames containing tournament,
chan-pak, chanpak, or metagraph) surfaced a theta/code support comparison,
an LRC exact covering-min paper, an LRC star/triangle refinement, and
the odd-square ordinal/Berggren descent. The theta/code note
`theta-code-lattice-gates-and-the-72-support-split.md` distinguishes
scalar cancellation from validity of the support that realizes it.
The LRC reflection
`lrc-nine-relation-slots-as-k4-edges-xor-axes-and-star-triangle-selector-codex-20260815.md`
suggests looking at coordinated edits rather than treating edges in
isolation. The Berggren note
`odd-square-ordinals-and-operation-compatible-natural-addresses-opus-20260823.md`
warns that a scalar address must retain enough information for the next move.
These are sources of the experiment, not claims of transfers between theorems.

## Physical conservation by cancellation, not neutrality

For a labelled tournament T and distinct unordered arcs e,f, let T^e be
the result of reversing e, and define

    J_T(e,f)=H(T^{ef})−H(T^e)−H(T^f)+H(T).

This integer is exactly the signed sum over vertex orderings that use
both edges consecutively (each edge somewhere in the ordering) and
satisfy all remaining arc tests. The sign is +1 if their required
directions originally both agree or both disagree, and −1 otherwise.
Indeed if either edge is unused, the four-corner sum cancels. If both
are used and all other tests pass, the ordering occurs at just one
corner. This argument does not require the edges to be disjoint.

For a pair A,B with gap g=H(B)−H(A), let d_e be the change in g when
e is flipped in both. Then flipping both e and f preserves g **iff**

    d_e+d_f+J_B(e,f)−J_A(e,f)=0.

This is a physical test, not an independence approximation: find two
individually nonneutral edits whose residual is compensated by the mixed term.

## Exact nine-vertex controls

Use the earlier 131/133 pair from EQUAL_GAIN_RECONTACT.md, with outgoing
bitmasks on labels 0..8:

    A=[252,189,504,128,8,24,58,112,251]
    B=[244,253,504,129,8,24,56,112,251]

Exclude their two distinguished return arcs (0,3) and (1,6) from any
common edit. Every one of the remaining 34 single flips and 561
double flips is checked by vertex-subset DP. No single flip translates
the window while preserving the signed +2 gap. Exactly four double
flips do:

| e | f | single-edit gaps | (J_A,J_B) | translated pair |
|---|---|---|---|---|
| (0,5) | (2,6) | −64, 66 | (12,14) | 325/327 |
| (0,5) | (5,8) | −64, 52 | (−46,−30) | 441/443 |
| (1,2) | (3,4) | 22, −12 | (14,8) | 57/59 |
| (3,8) | (6,8) | −142, 122 | (−32,−8) | 391/393 |

The same edits are made to each side, so the resulting witnesses
still differ at exactly the two original return positions. The proof
does not assume their former equal-gain interpretation survives.

There is a sharp scope limit: all five disjoint combinations of these
successful edit pairs lose the signed +2 gap. One cannot count the
two-edit moves as independent digits. The next constructive question
is whether changed sources can make the *higher* mixed differences
cancel while keeping enough such moves to span an interval at a
prescribed absolute height.

Run `python3 coordinated_recontact_replay.py` to reproduce the
complete search, signed differences, and disjoint-composition check.
