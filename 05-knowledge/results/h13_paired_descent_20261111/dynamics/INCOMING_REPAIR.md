# The H13 switch closes to a count-neutral repair at401

**FINITE-EXACT:** two literal nine-vertex rooted sources and their full
four-coordinate states. **PROVED consequence:** a two-arc H-neutral repair
that lowers every future positive construction response, with arbitrary
ordinal contexts. This is not an all-target or repeatable-return theorem.

The incoming
[uniform return](../../uniform_return_metagraph_20261111/reset/UNIFORM_RETURN.md)
supplies a positive model for using a temporary lower source and a specified
birth schedule rather than demanding an immediate same-count return.
Its new all-k branch is `34·2^k+5 -> 102·2^k+17`, retaining both ports and
all but one old vertex, and reaching a residual S5a family. Its separate
next-return problem remains open. The incoming
[birth atlas](../../infinite_birth_families_20261111/README.md) retains
actual common sources for101/103, rather than just the scalar availability
of both targets. Neither incoming claim is re-audited here.

That source-level question led to the following small native probe of the
present H13 twins. Start from Q=(18,8,3,5,14), insert root5 with cuts6
and3, and retain root5 throughout. Append:

1. a global sink6;
2. vertex7 beating only original vertex1 (outmask2);
3. vertex8 beating only the retained root5 (outmask32, the common Z).

The exact response trace is

| Stage | lower-cut source `(H,E,U,d)` | upper-cut source |
|---|---|---|
| original | `(31,7,45,13)` | `(33,7,45,13)` |
| sink | `(31,0,90,13)` | `(33,0,90,13)` |
| specified birth | `(95,7,313,41)` | `(101,7,307,41)` |
| common Z | `(401,41,585,41)` | `(401,41,573,41)` |

The final outgoing rows are

    (498,328,451,485,494,198,384,258,32),
    (466,328,483,485,494,195,384,258,32),

both rooted at5. Only the edges0--5 and2--5 differ. Thus the actual
H13 cut switch becomes an H-neutral two-edge repair after this fixed
three-birth wrapper. It retains every vertex and every other arc.

The middle step moves the difference to `(6,0,−6,0)`. Common Z then
cancels its H component and leaves `(0,0,−12,0)`. This is the exact
mechanism; scalar31/33 alone would not reveal it. In the common comparison
frame, lower minus upper after the wrapper is

    (0,12,12,12).

The upper-cut graph is therefore cheaper after every O/Z/J word containing
Z or J, and no more expensive after arbitrary fixed right-root amalgams.
The unchanged count itself is401. The q order and amalgam formulas are
proved in [PROOF.md](PROOF.md) and [ADDITIVITY.md](../modular_cycles/ADDITIVITY.md).

For any ordinal prefix R, retaining the root in the401 factor multiplies
all four coordinates by H(R). The same two-edge patch therefore preserves
count401H(R), with U decrease12H(R). For an ordinal suffix R, H and d
multiply by H(R), E becomes zero, and U multiplies by Q_R(2); the U decrease
is12Q_R(2). These follow by splitting the two path bins across the ordinal
factors. Both give actual infinite families of count-neutral native repairs.

The literal constants are checked by the standard subset path recurrence,
with one count for every singleton and transitions only along actual arcs.
The root-start table and complementary subset counts give U independently
of the last-step linear formula. The receipt checks the two changed edges
and the exact coordinate cancellation. Run

    python3 -B 05-knowledge/results/h13_paired_descent_20261111/dynamics/incoming_probe.py

The wrapper does not send every H13 right-address pair to a ternary branch.
For example at B1 its middle stage is265/265, rather than209/215. The
construction is attached to the named31/33 source; its arbitrary ordinal
context corollary does not supply an arbitrary count entrance.
