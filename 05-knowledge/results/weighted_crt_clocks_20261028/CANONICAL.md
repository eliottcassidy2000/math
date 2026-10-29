# A complete recursive attachment carrier, with an explicit complement bit

**PROVED** K1-K3 for the recursive ternary family; **FINITE-EXACT** verification
through depth2 by [check_canonical.py](check_canonical.py) and [JSON](check_canonical.json).
This extends the inherited folded joint-port action to arbitrary finite
palettes. It is not a universal canonicalizer for arbitrary strong atoms.

## K1. Joint colourings have a complete recursive normal form

Let f:V(W_k)->P be a vertex colouring by a finite totally ordered palette.
Write f_0,f_1,f_2 for its three child arrays, and alpha for digit negation.
For A_k use the canonical key

    K_A,0(f)=f(0),
    K_A,k(f)=min_cyclic_rotation(K_A,k-1(f_0),
                               K_A,k-1(f_1),K_A,k-1(f_2)).       (K1)

The minimum is over three cyclic rotations, not all six permutations.
For C_k form the paired palette array

    h(v)=(f_1(v),f_2(alpha(v))),
    K_C,0(f)=f(0),
    K_C,k(f)=(K_C,k-1(f_0), K_A,k-1(h)).                        (K2)

Two keys are equal iff the colourings lie in the same indicated group orbit.
Proof is induction using A_k=A_(k-1)^3 semidirect C3 and
C_k=C_(k-1) times diagonal A_(k-1). At each stage equality supplies exactly
the permitted root rotation and child transporters; conversely every such
transporter leaves the key unchanged. The paired palette is what retains
the common action on the mirrored children. Chosen minimizing rotations
can be retained if an explicit transporter, rather than only its orbit,
is needed.

For r named ports, colour vertex v by the subset of port names placed at v.
This permits repeated ports. For multiple binary attachments, colour it by
their ordered bit vector. Thus(K2) classifies their **joint** data at all
depths. This is the recursive counterpart of the inherited LCA-depth and
cyclic-direction classifier; no new group classification is required.

Every palette map applied uniformly to vertices descends to these keys:
apply it to a reconstructed representative and recanonicalize. If different
maps are specified on the immutable C_k single-port strata, they also descend.
An arbitrary vertex-specific operation need not descend.

## K2. Native complement has a one-bit selector

Use a binary colouring f(v) to specify the arc new_vertex->v when f(v)=1.
Reversing all arcs, then using alpha to identify the parent with W_k again,
acts on the attachment by

    S f(v)=1-f(alpha(v)).                                (K3)

Because alpha commutes with C_k, S is well-defined on its attachment orbits.
It is an involution and has no fixed orbit. Indeed the all-zero word o is
fixed by every C_k element and by alpha, while S flips f(o). This proves
that every complement-pair contains exactly one orbit with f(o)=0.

Consequently the exact complement/marking representation can be stored as

    (canonical complement-pair represented with root bit0, original root bit).

Native complement toggles the last bit. The pair key is unchanged. The
additional selector needs exactly one bit, because both orientations occur
and differ on the observed arc to the distinguished fixed vertex. This is
a scoped selector theorem for this parent family, **not** a solution to
the earlier unrelated seven-vertex-selector question.

The state counts under C_k and after quotienting complement are respectively

    2 and1; 8 and4; 192 and96; 887808 and443904

at depths0,1,2,3. The last row follows(W1), not enumeration of2^27 subsets.
For k>=1 the parent is strong. The one-vertex extension is strong iff f
is nonconstant: having both an outgoing and incoming neighbour suffices by
strong connectivity of the parent. Thus strong rooted extensions have
`c_k(2)-2` C_k-orbits, paired into `(c_k(2)-2)/2` complement-orbits. These
retain the parent and new root; rerooting the child can merge them further.

Every base-mark correction used here is realizable by the inherited full
centralizer holonomy. The incoming
[linear loop protocol](../session94_recursive_ports_crt_clocks_20261029/loop_protocol/LINEAR_LOOP_PROTOCOL.md)
improves the quoted upper bound to3(3^k-2k-1)/4 legal re-marking edges.
This is credited concurrent work, with explicit endpoint-construction cost.
The canonicalizer does not equate those edges with native complement.

## K3. What cannot be compressed away

Separate canonical keys cannot generally be joined without a shared frame.
At depth2, single ports at3 or4 have the same type, as do ports6 or7. But
the ordered pairs(3,6) and(3,7) lie in distinct C_2-orbits: the mirrored
blocks share one cyclic shift, so the folded relative offset is invariant.
The complete census is5 one-port types,33 two-port types and261 three-port
types. In particular33 exceeds the25 independent pairs of marginal types.
Canonicalize the product palette **before** forgetting the joint source.

There is also an information lower bound for unrestricted binary patterns.
Let N=3^k and c=(N-2k-1)/4, so |C_k|=3^c. Any complete orbit carrier must
distinguish at least2^N/3^c states, and hence requires at least

    log_2(number of orbits) >= N-c log_2(3)               (K4)

bits in the counting sense (round up for a fixed-length binary encoding).
This grows linearly in N. The explicit canonical keys give a finite recursive
description; they do not promise constant-size storage as depth increases.
Even after identifying complement the lower bound drops by only one bit.
Weights, joint sources, and the allowed operation menu determine which
additional quotients are permissible; orbit counts alone do not settle that.
