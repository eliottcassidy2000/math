---
id: THM-5117
title: "Physical cut switches, opposite-facing response cancellation, and primitive but exponentially separated affine cubes"
status: "PROVED, independently audited. Scoped constructions and exact spacing; no all-target coverage claim."
session: Session245, 2026-11-10
source: 05-knowledge/results/session245_physical_spectrum_20261110/interfaces/PROOF.md
audit: 05-knowledge/results/session245_physical_spectrum_20261110/audit/AUDIT.md
---

# A unit switch outside the previous word-family divisor

The five-vertex deletion tournament with outgoing masks

    [18,8,3,5,14]

and selected cuts {0,1} or {1,2} has seven-coordinate interfaces

    (13,104,32,41,6,5,6),
    (13,104,32,44,8,4,6).

The coordinates are (h,q,A,G,T,S,E), as defined in Session244's exact ordinal
law. Add an unselected sink and then a root beating precisely the selected
cut. The resulting strong seven-vertex graphs have rooted states

    (H, root-ending count, root response, root-deletion count)
       = (71,13,77,13), (69,13,77,13).

This is a pure total-count gap of 2 with deletion count coprime to 3. It
escapes the restricted pure-gap word compiler's common divisor. It is NOT a
full-interface pure-T switch: G and S also differ. Independent enumeration
recounts all interface coordinates and both actual graphs.

## Reversal cancellation

For an interface x, graph converse with the same selected cut sends
(h,q,A,G,T,S,E) to (h,q,G,A,T,E,S). Use Session244's raw JJJJ/JOZJ modules
from physical sources i=1,2; their switch difference is
(0,0,E_i,-4W_i,0,0,-U_i). Here E_i,U_i are source coordinates and W_i is
its rooted ordered-three-path response.

Place module 1, a transitive unselected separator of t>=1 vertices, and the
converse of module 2. If

    2^t E_1 q_2 = 4 h_1 W_2,

simultaneous switching preserves h,q,A,S,E and changes

    G by -4 h_2 W_1 + 2^t q_1 E_2,
    T by 2^t(A_1 E_2 + E_1 A_2 + E_1 E_2).

Protected closure therefore changes only its total count. Explicit strong
examples have counts 589901/589089 at order 19 and 15669451/15654639 at order
22; the latter has deletion 47515, coprime to 3. The order-19 total counts
are literally subset-recounted. The order-22 counts use independently
recounted constituent interfaces and the proved gluing law, not a full
order-22 subset enumeration. This is a positive escape from the previous
same-orientation cancellation obstruction, not arbitrary gap placement.

## An all-length primitive cube with no small spacing

A second deletion core [32,17,3,7,13,30], with cuts {0,4}/{0,5}, has interfaces

    (37,352,96,102,12,12,10),
    (37,352,96,96,32,10,10).

Append one unselected sink to each core. Ordinal-compose one h13 module and
k-1 h37 modules, then take protected closure. The resulting strong graphs
have 7k+1 vertices. All 2^k independent choices are an exact affine scalar
cube. Its coefficient gcd is 2 for EVERY k>=1. Nevertheless, for k>=2 its
2^k distinct counts have minimum spacing exactly

    50 * 37^(k-2).

Proof: the first coefficient is -2*37^(k-1). For later positions i use
c_i=37^(k-i)(6A_(i-1)-20h_(i-1)), with
h_1=13,A_1=64 and h_(j+1)=37h_j,A_(j+1)=704A_j+192h_j.
The boundary A_2=47552 is 7 mod 37; later A_j stay 7 mod 37. Thus the last
coefficient is not divisible by 37, while all coefficients are even.
For spacing, b_j=6A_j-20h_j satisfies b_1=124 and
b_(j+1)=704b_j+14492h_j. After the initial pair -74S,124S,
S=37^(k-2), successive positive coefficients grow by more than 19 and
cannot cancel the earlier sums. The initial minimum 50S is attained.

The exact spacing is the strongest hostile to replacing a physical subset
sum by its gcd. Removing divisibility is real progress but does not create
an interval or a predecessor for an arbitrary odd target.

[Full proof and reproduction](../../05-knowledge/results/session245_physical_spectrum_20261110/interfaces/PROOF.md).
