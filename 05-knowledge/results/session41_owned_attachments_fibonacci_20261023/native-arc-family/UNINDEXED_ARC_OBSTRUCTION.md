# A native order-31 obstruction: an arc occurrence is not an indexed flag

Status: FINITE-EXACT native source/path certificate with elementary proof.
This is about an explicitly exhibited 21-path source-automorphic family, not
all Hamiltonian paths or all tournament classes at order 31.

Let A=D(D(P7)), using the supplied literal doubling convention. Its vertices
are four seven-copies 8i+x (0<=i<=3, 0<=x<=6) and centers 7,15,23. The 21
permutations fixing all centers and acting simultaneously by x->a*x+b,
a in {1,2,4}, b in F7, are actual automorphisms of A. Consider the path

23,0,8,24,16,1,9,25,17,2,10,26,18,3,11,27,19,4,12,28,20,5,13,29,21,6,14,30,22,15,7.

The standalone checker constructs A, verifies this Hamiltonian path P, verifies
every one of the 21 automorphisms, and checks that their images of P are
21 distinct paths with the same actual endpoints (23,7).

## The precise observation distinction

An unindexed arc observation asks whether a literal ordered arc (u,v) occurs
SOMEWHERE in a candidate path. An indexed flag observation asks whether its
vertices at prescribed positions i,i+1 are exactly u,v. A source automorphism
g fixes the latter observation exactly when it fixes u and v, but it may satisfy
the former by carrying a DIFFERENT arc of P to (u,v).

For the displayed P, every one of its 30 arcs occurs in at least two of the
21 candidate paths. Therefore no single unindexed arc observation identifies
P. This is not a marginal synthetic example: all candidates are actual paths
in one fixed source.

There are exactly 13 unordered pairs of arcs whose joint occurrence identifies
P. The first pair in path order is {(23,0),(16,1)}. Thus two unindexed positive
arc observations are necessary and sufficient within this explicit family.
The first arc fixes b=0 because it is the initial center-to-copy transition.
The second then distinguishes the remaining multipliers.

By contrast, the single indexed flag at zero-based positions 4 and 5 is
(16,1). Fixing these two distinct F7 coordinates forces a=1,b=0, so only P
has that indexed flag. One indexed flag and one unindexed arc are genuinely
different query units. The extra information is alignment/position, not an
extra edge color or a parity bit.

## Full inherited-control probe

A separate probe checks the 933 inherited actual-path controls using the
attachment lane's complete source-group solver. Five originally selected arcs
fail the unindexed-occurrence test. Control 697 is the displayed path and has
no occurrence-unique arc at all. Every other nontrivial control has some
occurrence-unique arc. These are finite inventory statements only.

The standalone certificate does not depend on that imported solver or the
inherited compressed receipts. Run check_unindexed_arc_hostile.py normally
and with python -O; both produce identical PASS receipts. It retains actual
vertex names, the actual source and the full ordered path. No statement about
complete native metagraph wiring follows from this single family.
