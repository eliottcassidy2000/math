# Independent coverage continuation audit

**Verdict: ACCEPT the scoped structural results and selected exact controls.**
The complete large quotient censuses remain producer-only computations;
this audit does not silently promote them to independent replays.

The final continuation source `coverage/SYMMETRY_EDIT.md` is pinned at
`f02526c196353835e4740fc765017d79ffe0f4c61a479bb27002f3bbea645c39`.
Final inspection corrected one dimensional typo: the k=6 deletion cards
are order-seventeen parents of order-eighteen children, not order-eighteen
parents. The mathematics and numerical receipts were unchanged.

## 1. Canonical completion of a deleted triangle block

Reviewed `coverage/TRIANGLE_CARD_RIGIDITY.md`, SHA256
`fb8a4f8e7cb81a94d8a8d4db0bb4a0178e9758981956abe276ed022138f224f1`.

For every tournament Q, including nonstrong Q and order one, delete v from
T=Q[C3,...,C3]. The surviving two vertices of its block are the unique
homogeneous pair. Every other pair either lies in an intact triangle and
is distinguished by its third vertex, or meets such a triangle once and
is distinguished by one of the other two blockmates. Orient the unique
pair u->w. The original beaten set is exactly

    S0={u} union (N_B^+(u) minus {w}).

Any completion belonging to any C3-substitution family must use this pair
as the new vertex's two triangle blockmates. Cyclic orientation and the
external homogeneous orientation force precisely S0. Thus uniqueness is
family-wide, not merely uniqueness up to the original isomorphism class.

Consequently card isomorphisms extend canonically, Aut(B) is the stabilizer
Aut(T)_v, and two cards of T are isomorphic exactly when their deleted
vertices are in one automorphism orbit. Nonisomorphic members of the
family have disjoint deletion-card supports. Original blocks are exactly
the homogeneous directed C3 triples, also yielding
Aut(T)=(C3)^k semidirect Aut(Q).

This structural theorem does not prove uniqueness of an H-restoring cut.
It proves that any nonoriginal restoring cut genuinely exits the entire
structured family. The all-order scalar uniqueness question remains open.
The independent literal finite controls are in `audit_card_rigidity.py`:
all 1,099 labelled quotients through five, their 16,215 deletion cards, and
18,828 completion cuts for every deletion through quotient order three.
No producer imports are used. Normal/-O receipts match at SHA256
`cd49b7174bb4f505e0a712775998ca4b5668d256e04dacf242f45cf383d3b53a`.

## 2. Selected order-eighteen row, without orbit compression

`audit_large_blowup.cpp` independently builds the six-vertex quotient of
binary arc code 16, replaces each vertex by a directed triangle, and
deletes vertex zero. It computes forward endpoint subset DP, reverse
endpoint DP, and the quadratic insertion polynomial from endpoint gluing.
This is not the producer's cardinality-profile compression route.

The complete 131,072 literal cuts of the order-seventeen parent give

    original H = 1,049,305,833
    parent H   =    28,273,779
    response gcd = 54
    unique restoring cut = 114689.

There is no response differing from the original H by +/-2. Four selected
order-eighteen children are also counted by direct Hamilton subset DP.
Signed 64-bit arithmetic is safe here: Hamilton paths and ordered
two-path-cover totals are bounded by 18!, and a conservative bound of
4*18! on the insertion-polynomial intermediates is below 2^63. The
endpoint kernels are actual nonnegative path counts; the subsequent
linear/quadratic polynomial is evaluated with signed arithmetic.

    g++ -std=c++17 -O2 -Wall -Wextra -Werror audit_large_blowup.cpp -o audit_large_blowup
    ./audit_large_blowup > large_blowup.normal.json
    g++ -std=c++17 -O1 -g -Wall -Wextra -Werror -fsanitize=undefined audit_large_blowup.cpp -o audit_large_blowup_ubsan
    ./audit_large_blowup_ubsan > large_blowup.ubsan.json

Receipts match with SHA256
`577ccc83c5ef273e2f5c248582da749660fd366bc4358986bad8f4bf302629e2`.
This is one independently reproduced row, not the complete k=6 or k=7
census. The producer's k=7 claim is 353 quotients, 2,471 deletion-orbit
rows, 40,484,864 profiles representing 2,591,031,296 literal cuts, with no
zero escape. Its all-row computation is not independently reproduced here.

The final producer also supplies a separate forward-DP coefficient control
for Paley7 quotient code 85298: 212 literal order-twenty-one children recover
all 211 coefficients of one order-twenty deletion parent. The stated parent
H is 11,160,377,792,361 and original child H is 123,156,134,945,235. This is
a producer cross-route control, not this auditor's independent replay.
The elementary range bound `H(T)<=n!/2^floor(n/2)` follows by retaining only
disjoint consecutive-position arc constraints. Together with the endpoint
exposed-word bound `19H(B)+20!` at parent order twenty, it validates the
documented signed-64-bit range; it does not authorize larger undocumented
orders of the same C++ implementation.

## 3. Formal-capacity firewall

`audit_formal_capacity.py` independently reconstructs endpoint kernels for
the actual order-eleven parent with rows

    [1794,1792,11,19,7,95,159,63,764,1276,508].

Its H is 24,651 and actual restoring cut for target H=241,137 is uniquely
1793. Perturbing the six specified capacities by +/-498 preserves capacity
total 582,822, response gcd six, every cut-cardinality layer sum, the three
intact block rotations, positivity, pair-deletion floors (maximum 4,194),
and the original child. Nevertheless the formal network's restoring cuts
are 480,736,1248,1793. The first new cut rises from actual 239,643 by 1,494.
Six actual children are additionally checked by direct path counting.

This is a deliberately FORMAL capacity tensor. It is not the actual path
kernel of the fixed parent and is not asserted realizable by any other
tournament. It refutes sufficiency of the listed numerical invariants for
the proposed uniqueness deduction, not the tournament uniqueness theorem.
An actual nonstrong quotient control retains its multiple restoring cuts,
so strongness is not silently removed from the scalar conjecture.

Normal/-O receipts have SHA256
`5a7cb8f9b028f387f2e3d38c71b39312fcb42666b71a58445df7687e1c44f7f6`.

## 4. Proof firewall retained

The marked-boundary rotation proof gives the all-order divisor
`3^max(k-3r,0) | H(T')` when only arcs incident to r marked vertices are
altered, and parity doubles the divisor for H(T')-H(T). This does not imply
that zero-response walks remain in the C3-substitution family. Canonical
completion sharpens the missing condition but does not supply it. Neither
finite zero-uniqueness nor formal-capacity positivity closes that gap.
