# Native suspension already loses a harmonic at order five

24 October 2026. **FINITE-EXACT:** complete order-five and order-seven native
line metagraphs, exact old harmonic bank, suspension readouts and dual
certificates. **PROVED from the certificates:** the full old harmonic bank
does not generally transport through literal source suspension, even modulo
two. The target-570 order-eight-to-ten question remains **OPEN**.

This is a separate next-size packet. The earlier complete order-four-to-six
positive result remains true and unchanged. Neither result is extrapolated
to all tournament orders.

## 1. The tiny decisive certificate

Use the canonical merged node indexing fully specified in `normal.json`.
The complete original order-five graph has10 nodes. Its vector

    x = 1 on old nodes {3,4,6,7}, and 0 elsewhere

satisfies L5 x=0 modulo2, including the old root0.

In the complete272-node order-seven metagraph, an exact sum of four FULL
Laplacian rows is

    row10 + row41 + row196 + row244 = e10 + e41 mod2.           (1.1)

Thus every new mod-two harmonic z satisfies z10=z41. This is not a reduced-
root artifact: (1.1) is checked on all272 coordinates.

The two relevant readouts are literal source suspensions:

    old node3  -> S(T3) at new node10;
    old node9  -> S(T9) at new node41.

They can be reconstructed without trusting class labels. The respective
old forward words are55 and0, in the six free-chord-bit convention, with
canonical merged keys16 and258. Suspending their actual matrices gives
canonical new keys256 and8208. Both old matrices, both suspended matrices,
and the exact node maps are in the receipt.

Since x3=1 and x9=0, the requested two evaluations contradict (1.1).
Therefore no new harmonic can retain even these two old source values.
Literal target-memory complications are not needed for this failure.
There are no conflicting record labels on any selected new canonical node.
The obstruction is an actual additional modular electrical relation.

One old source used here has a source and sink. Accordingly this does not
assert a counterexample inside a restricted nonextreme-source family such
as the particular incoming source2393. It disproves unrestricted harmonic
functoriality of suspension, not the specific target-570 construction.

## 2. Exact extent of survival

All64 old forward words are suspended with their actual complements. Their
source and target records give48 distinct readout pairs. The old grounded
mod-two harmonic bank has dimension 3 and 8 elements. Exactly 4 of these 8
grounded characters extend to the new graph subject to all 48 readouts.

These counts use root 0 in each carrier. Allowing constant potentials doubles
the old bank to 16 elements, of which 8 extend: constants always extend as
constants, and grounding removes exactly that common constant direction.

The extendible subspace is exactly

    span over F2 of e_old1 and e_old2.

Both generators have actual new harmonic extensions checked against every
new row and every readout; their four linear combinations are retained.
Together with x above they form a basis of the entire old grounded harmonic bank.
Relation (1.1) kills precisely the remaining direction. Thus the surviving
subspace has dimension2, not merely two successful sample vectors.

The old two-primary critical group is (Z/2)^3. The common-lift theorem in
`../TRANSFER.md` therefore says that full all-dyadic transport would already
be equivalent to mod-two transport. Its failure here is decisive for full all-level transport. This does not
claim that every separately scaled character at a higher modulus must fail: a
higher-modulus lift can exist without supplying the missing mod-two lift, as
the nonsplit example in TRANSFER.md demonstrates.

## 3. The actual prime eleven, with its separate meaning

Exact integer elimination gives

    det A5 = 41096 = 2^3 * 11 * 467.

The old mod-two nullity is3, so all three dyadic elementary factors have
order2. The old kernels at11 and467 each have dimension1. Together these
facts give the complete old Smith factors

    1,1,1,1,1,1,2,2,10274.

In contrast, the new271-by271 reduced Laplacian is invertible at both odd
primes. Exact field elimination gives

    det A7 = 8 mod11,      det A7 = 397 mod467.

Thus its grounded harmonic bank is zero at11 and467, and the old nonzero
characters at either prime cannot transfer through the source readouts.
The checker also retains explicit left-annihilator obstructions for the
full readout systems at both primes.

Here eleven has a genuine, authenticated role: it is an old electrical
critical-group factor that disappears at the new order. It is not the
half-cost of a22-clock move, a selected clock, or a numerical explanation
of the22-bound. No such cross-carrier identification is asserted.

## 4. Completeness and what the certificate does not claim

The full carriers are independently reconstructed from all forward words:
64 words and32 literal complement lines at order5;32,768 words and16,384
literal complement lines at order7. Every unordered free-bit-complement
pair is counted once as one literal line. Actual reflected black partners
are not collapsed into one line; loops retain their records and contribute
zero to the Laplacian. The original base reconstruction is described and
independently audited in `../SUSPENSION.md`.

Canonicalization exhausts permutations inside sorted score classes and
merges with actual transpose/converse. Each degree-preserving candidate
permutation is checked by its complete matrix, not just by its scores.
The receipt includes every word owner, every weighted edge, the complete
Laplacians, old basis vectors, positive lifts, and negative dual witnesses.

The result is not an order-eight-to-ten calculation and does not identify
which of its much larger old bank transfers. It does not alter source2393's
literal three-path suspension, release an old owner pin, invent a kernel
resource, or imply a new tournament wiring. It establishes that a general
suspension-only argument cannot bypass the harmonic-transfer equations.

## 5. Reproduction

Run from this directory:

    python3 -B check.py --out normal.json
    python3 -O -B check.py --out optimized.json
    cmp normal.json optimized.json

The checker imports the hash-pinned finite-carrier builder from
`../native_suspension_probe.py`, but builds both new carriers from scratch.
It uses exact integer arithmetic, finite-field Gaussian elimination and
fraction-free determinant elimination. Its checks use explicit exceptions,
not assertions. The tiny four-row certificate is checked independently of
the solver's returned dual row, and positive extensions are multiplied back
through the complete original system. `MANIFEST.json` freezes this extension
separately from the earlier core and order-four-to-six packet.
