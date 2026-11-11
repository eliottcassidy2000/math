# Exact count collisions lie outside the positive carry cone

**PROVED, scoped.** The fixed-chart lemma below identifies a necessary signed
cancellation in an address collision. **FINITE-EXACT:** two literal named
tournaments at addresses 37 and 39 both have 611 Hamiltonian paths. This
does not rule out a different coordinate frame or a richer signed compiler.

The incoming Session255
[finite support/residue theorem](../../session255_modular_collatz_20261111/cross_frontiers/PROOF.md),
Section 1, gives an exact finite equality compiler for nonnegative operation
matrices and a coordinatewise nonnegative initial difference. The latter
hypothesis is essential. Our [carry matrices](PROOF.md#c3-an-exact-carry-lift-for-every-odd-multiplier)
are nonnegative, but this alone does not make two starting states comparable.

## Fixed-chart collision lemma

Fix a named core K, cut O and left address A. Write the inherited count as

    Phi(A,B) = d*g(B) + c*B + v*(B mod 2) + k,
    d=H(K)>0, c=q*A+s,

where k is constant in B. Suppose `B2>B1`, the addresses have the same
parity, `c>0`, and `Phi(A,B2)=Phi(A,B1)`. Subtraction gives exactly

    d*(g(B2)-g(B1)) = -c*(B2-B1) < 0.                   (1)

Thus the second address increases while its g-coordinate decreases. In the
multiplier-one carry frame

    V1(B) = (g(B), B, [B even], [B odd]),

the difference `V1(B2)-V1(B1)` has both signs. Reversing the order does not
make it nonnegative. More generally, if `B_i=m*n_i+r` for one fixed channel
r of an odd-multiplier lift, the height coordinate n increases and channel
`G_r(n)=g(m*n+r)` decreases, so that frame is also incomparable.

The exact equality additionally requires
`d/gcd(d,c)` to divide `B2-B1`; equal parity adds the factor 2 through the
least common multiple. This arithmetic filter is necessary, not sufficient.
If `c=0`, (1) forces only equality of g; the mixed-sign conclusion does not
hold. Different parities introduce the extra signed term v and are outside
this particular lemma.

## A literal collision retaining core and cut

Take the directed three-cycle `K=(2,4,1)`, named out-cut `O={0}` (mask 1),
and left address `A=1`. The inherited coefficients give

    Phi(1,B) = 3*g(B)+15*B+10+(B mod 2).

The two right addresses have the same six-bit width:

| B | Binary address | g(B) | Phi(1,B) |
|---|---|---|---|
| 37 | 100101 | 15 | 611 |
| 39 | 100111 | 5 | 611 |

Both graphs have eleven named vertices: left vertex 0, core vertices 1–3,
right vertices 4–9, and apex 10. They preserve this naming, the core and the
core cut; their right-flank attachments differ. The receipt retains all
adjacency rows and checks the tournament axioms, the inherited forward path
DP, and an independent predecessor DP. Thus 611 here is an actual count,
not just a formal value of a scalar expression.

Their carry-state difference is `(-10,2,0,0)`. The count equality is the
cancellation `3*(-10)+15*2=0`. The equal counts do not persist under even
one common binary suffix:

| Common appended digit | New addresses | New exact counts |
|---|---|---|
| 0 | 74, 78 | 1279, 1315 |
| 1 | 75, 79 | 1181, 1211 |

These suffix values are exact compiler readouts, not extra graph census
claims. They isolate why equal scalar counts are not an interchangeable
starting source for a future construction.

## What the incoming theorem does and does not transfer

For comparable carry states, nonnegative transitions do give exact support
propagation. A nonnegative output row can be handled by requiring the final
support to avoid all its positive columns, or by adding it as an output
coordinate. This is a legitimate application of the incoming theorem.

The 37/39 pair fails its initial-comparability hypothesis. Erasing magnitudes
or keeping only signed support would discard the cancellation in (1).
Modulo-49/343 phase and other finite congruence controls do not restore that
discarded integer equality by themselves. A usable extension must retain
the needed signed difference or prove comparability in another actual
common frame. No claim is made here about comparability of these graphs in
the different rooted O/Z/J response frame, nor about impossibility of a
signed finite compiler for a restricted family.

## Reproduction

    python3 -B 05-knowledge/results/central_modular_chanpak_20261110/compiler/incoming_height_check.py

The checker constructs exactly the two displayed eleven-vertex graphs and
checks the four suffix readouts. Normal and `-O` receipts must agree. The
incoming proof is pinned to commit `bf945ba947` and its content hash; this
comparison can be reproduced before or after integration of that commit.
