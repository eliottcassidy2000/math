# Independent root review of the recovered arithmetic and norm lanes

Status: **PASS after the Session321 correction**, 2026-11-13. This is a mathematical and source-code review,
not a formal proof-assistant certificate. The author of this review did not
write the reviewed bridge or norm implementations. Status promotion, if any,
is separately bound to exact bytes in `promotions.json`.

The initial review snapshot `2e49f617f827913dc53014d863822d1ab77cc1ff187a6ae01437d81a746c15f1` incorrectly repeated the named-kernel gauge caveat. The stronger incoming identity repairs that clause; `bridges/CORRECTION.md` preserves the producer lineage, and the dated MISTAKES entry records the failed inference. The other reviewed conclusions remain unchanged.

## Reviewed candidates and evidence

| Artifact | Reviewed SHA256 |
|---|---|
| bridges/RESULT.md |29fc9ca464171bd1f5ef30602786cac2bbe222bda2b5b7cab9deb34e7e9699b8|
| bridges/check_bridges.py |89fae3be3a595d5de0b7427485092c24a5099e9cba3894718a08746775f3da95|
| bridges/normal.json and optimized.json |1d98df0f4bbc99b7f4dc6969eb8f6291a22b6ac21a6998c40b4e5ca3b28a82f9|
| census/NORM_CLOSURE.md |ceb1d2279bf5be30a2b51f13efb2e90de11facd417d739bb05c06dcb8382ae6b|
| census/check_norm_closure.py |fac0fb7efe361756dfc90e6bb1c36dc920e6fd0d477a221359c869b48011b05d|
| census/norm.normal.json and norm.optimized.json |bc3e560a9eefc6b1fc380ef8fccd6b81aa02d3a380fbddcd56d1cc045e266ed3|

The root replayed both producers into separate temporary outputs and compared
them byte for byte with the saved receipts: PASS, respectively 19,447 and
12,434 explicit gates. Those replays check reproducibility, not algorithmic
independence. The different kernel, Legendre-symbol and path-cover algorithms
inside the bridge checker were separately inspected.

## B1: declared path observer to a declared first descent — PASS

The two-slot mixed-incidence bit uses opposite local arc directions. Its
integration recovers either global orientation; transpose plus a scalar
identity shift gives exactly the same Phi3 kernel dimension under converse.
The incoming identities additionally make Phi symmetric. It therefore determines
the full labelled Phi and its embedded kernel, including named vectors, without
picking a preferred orientation. M retains the orientation choice. On the
Phi kernel, the doubled Monsky coordinate map is converse-invariant for even n
and changes by (0,x) for odd n, since every kernel vector has even coordinate
sum. The corrected producer replay
checks full matrix and kernel equality, not merely equality of dimensions.

The cited input was independently reopened at
[Watkins et al., section4.1, printed page8](https://dms.umontreal.ca/~andrew/PDF/RankQuadTwist.pdf).
Its odd squarefree parameter matrix is
`[[M+D2,D2],[D2,M+D(-2)]]`, where the row-zero M uses quadratic nonresidues.
For primes3 mod8, D2=I and D(-2)=0, giving precisely the displayed
specialization. Solving the first block equation gives
`y=(M+I)x`; the second becomes `(M^2+M+I)x=0`. Thus the claimed map and
its pure 2-Selmer interpretation have the correct direction and scope.
The Kummer distinction between rank and Sha[2] is retained. No low-corank
preprint assumption is used to infer full BSD.

## B2: smallest arithmetic loss and all-order extension — PASS

The explicit code convention matches the two literal arc lists, binary
matrices, kernels, H=9 and full Q=t^5+6t^3+2t. Primality, residues, and all
Legendre symbols of both supplied marked realizations are tested directly.
The two different integers N define two different elliptic curves.

The ordinal product of numbered optional slots is a bijection, so it proves
Q multiplication for all orders. In the block matrix, M(U)J=0 and
JM(V)=(|V|-1)J hold over F2; adding the parity shift and the linear term
cancels the off-diagonal block of Phi. Thus the claimed direct sum is exact,
not just determinant factorization. Transitive suffixes have zero nullity.
The construction consequently proves the loss at every n>=5; finite
enumeration through4 establishes the stated minimal starting order.
The conclusion concerns (n,Q,scores); it is not an impossibility theorem for
all scalar A000568 recurrences.

## B3: exact Frobenius restriction and full-group obstruction — PASS

Both reciprocal cubic coordinates and their semilinear action are retained.
Phi3(Frobenius) is the trace, whose nonzero kernel is the exponent orbit
{1,2,4}. The F4 module is an operator module, not a subfield of F8. The
displayed alpha^2 to alpha^3 counterexample destroys C7 invariance. Distinct
irreducible cubics split every C7-invariant subspace into the two simple axes,
which rules out linear quotient dimensions2 and4. This proof concerns the
full G21 action on W, and does not rule out nonlinear maps or restriction
to the Frobenius C3 subgroup. No arithmetic realization is obtained merely
from the repeated polynomial Phi3.

The HN recovery retains all three pair channels, the choice of channel
assignment, and geometric/coloring realization. The thirty-code controls
verify pairwise blindness and cubic faithfulness in that declared binary
carrier; they do not prove a new chromatic number for the plane.

## Norm closure — PASS

The cyclic group of order63 has common kernel of the powers21 and9 of
order3, exactly F4*. The multiplication quotient is therefore already
closed. Translation is the new operation. On every scalar orbit outside
F4, independence of1,x over F4 forces distinct translated observations.
On F4* the only surviving pair is phi,phi^2. Multiplication by b outside
F4 sends these two points into already distinct singleton cells. This
proves exactly22->63->64 and the separate stable controls, with no assumed
tournament realization. The field operations, actual order63 of b and
hostile order9 of the polynomial root are explicitly checked.

## Recovery scope

The 33-row catalog has26 inherited index rows and seven marked supplements.
It is an artifact and route recovery, not an audit promoting every historical
mathematical sentence in the role index. The actual census/request evidence
and missing request provenance are separated. Root inspected those scopes
and all seven supplemental routes; no universal numerical interpretation
is inferred from their adjacency in the table.
