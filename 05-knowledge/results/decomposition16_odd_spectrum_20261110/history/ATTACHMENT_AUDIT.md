# Independent audit of the H3 one-birth compiler

**Verdict:** the mixed-contact formula, its finite target inverse and the
marked generating function are sound. The draft's uniform-contact
collapse needed the reentry correction below. The independent verifier
checks the repaired formula for all four triangle-contact cardinalities.
The unrestricted odd tournament spectrum remains OPEN.

## The error caught before publication

The first draft claimed that if the new vertex x has uniform contact
with C3, then C3 is contiguous in every Hamilton path, so replacing it
by one vertex and multiplying the count by3 is valid. This is false.

Let the source be `0 -> C3(1,2,3) -> 4`, where the triangle is
`1 -> 2 -> 3 -> 1`. The new vertex5 beats0,1,2,3 and loses to4.
The adjacency row masks are

    (30,20,24,18,32,15).

The full count is27. The collapsed four-tournament has rows
`(6,4,8,3)` and count5. Moreover

    (1,2,4,5,0,3)

is an explicit Hamilton order which splits the triangle into two runs.
It leaves the triangle through the right flank and reenters through the
left flank. The new vertex need not itself touch the triangle on that
excursion. Thus uniform contact with the module is not a contiguity
certificate in a cyclic exterior quotient.

The first failed implication was `uniform module -> one module run`.
The strongest survivor is divisibility by3, also following directly from
the triangle-rotation automorphism acting freely on ordered Hamilton paths.
The missing coordinate is the exterior return route and hence the
ordered two-path cover of the triangle. This is the same operation-level
hazard that the old fixed-H and unmarked-cover obstructions warn about.

## The repaired all-contact formula

Retain A,B,alpha,beta,g and k from
[the producer proof](../attachment/PROOF.md). For every k=0,1,2,3 the
same odd-cycle calculation gives

    H = 12AB + 3g(A) + 3g(B) - 3
        +3(3-k)A + 3kB + (3-k)alpha + k beta
        +2*[k=1 or k=2].

For mixed contact this is exactly the producer's C(A,B) or C(B,A).
For k=0 put

    U(A,B)=12AB+9A+3g(A)+3g(B)+3alpha-3;

for k=3 exchange A and B. If w_0,w_1 are the collapsed one-apex words,
the exact price of reentry is

    U(A,B) - 3g(w_0) = 6A(B+1),
    U(B,A) - 3g(w_1) = 6B(A+1).

The original multiplication by3 holds precisely when A=0 for k=0 or
B=0 for k=3. This is an equality boundary, not a blanket ban on module
collapse. The displayed hostile has A=B=1 and excess12.

The proof's odd-cycle enumeration has been checked term by term. Every
non-old odd cycle uses x exactly once; the only old nontrivial cycle is
the triangle. If the new cycle avoids the triangle, its residual
collection weight is3; otherwise it is1. Increasing flank paths give
the A,B evaluations at1 and the last-bit evaluations at-1. For uniform
contact the internal out-to-in triangle path term is absent; retaining
the other terms yields U, rather than3g of the collapsed word. No further
old cycle collections are omitted.

The finite inverse of C remains valid because `g(A),g(B)>=1` and the
mixed expression has lower bound `12AB+6A+3B+5`. Uniform outputs are
divisible by3, so this correction does not invalidate the all-order
exclusion of an H3 parent for target11. The marked-class formula also
retains its required hypothesis `3 does not divide target`: deleting
the marked vertex recovers the unique triangle and both ordered flanks,
and k=1 versus2 cannot be identified by a mark-preserving isomorphism.

The uniform inverse also terminates: A=0 reduces to3g(B), and every
one-apex value has an even representative2q by repeatedly deleting a
trailing one bit; `g(2q)>=q+1` bounds q by the target. For A>=1,
`U(A,B)>=12AB+9A+3` gives the producer's finite box. Independent target
inversion through55 then gives exactly7,21,55 missing from the combined
H1/H3 parent lanes. Thus every vertex deletion of an H55 tournament
has count at least5. The seven-vertex witness and all seven deletion
counts are independently checked by literal permutations, including
the deletion with count5. This is a sharp **parent-lane** obstruction,
not a global gap at55.

## Independent finite universe and controls

[attachment_audit.py](attachment_audit.py) imports no producer code. It
constructs dense adjacency matrices directly, reads A and B with explicit
integer sums, computes g by endpoint-gap sums instead of the producer's
binary recurrence, and counts paths by a backwards suffix recursion
instead of the producer's forward subset DP. Literal permutations supply
a third check on every child with at most six vertices.

- All1,032 cuts on every `TT_a -> C3 -> TT_b` with a+b<=4.
- All136 children of order at most6 additionally checked by permutations.
- All258 uniform-contact cuts checked against the exact reentry correction.
- Positive zero-excess cases and the explicit excess12 hostile retained.
- Independent complete target inverse through55, plus literal H55 witness
  and all seven deletion counts.

These are exact checks of a declared finite universe. The all-order
claims depend on the displayed proof, not extrapolation from this census.
Normal and optimized outputs must agree byte for byte; each pins the
producer proof, producer script and independent script hashes.

    python3 -B 05-knowledge/results/decomposition16_odd_spectrum_20261110/history/attachment_audit.py
    python3 -B -O 05-knowledge/results/decomposition16_odd_spectrum_20261110/history/attachment_audit.py
