# Exact first-owner correction at total-mass deficit 2k

**PROVED for k>=4, pending independent audit.** The finite k=2,3 tests are
separately labelled below. This note stops at delta=2k as requested.

Use the model, parameters and notation of `two_unit_cells.md`. Let
L_k(delta) denote its right-hand side (1), evaluated with either the
ordinary or the strong parameter row. It is a formal binomial expression
outside the previously proved range; at delta=2k it is not yet a class
count. Then the exact first-owner correction is

    Z_{k,2k} = L_k(2k) + K_k - 2 H_k,                    (1)

where K_k counts pairs of marked adjacent matching edges and H_k counts
marked matching edges whose endpoint positions differ by two. In the
strong version both are restricted to proper-cut-crossing matchings.

Generally,

    K_k=(k-1)F_{k-1},   H_k=(2k-2)F_{k-1},
    Z_{k,2k}=L_k(2k)-3(k-1)F_{k-1}.                      (2)

Strongly, with H^s_1=0,

    K^s_k=binom(2k-4,2) I_{k-2},                        (3)
    H^s_k=(2k-2)F_{k-1}
          -sum_{j=2}^{k-1} H^s_j sum_{a+b=k-j}F_a F_b,  (4)
    Z^s_{k,2k}=L^s_k(2k)+K^s_k-2H^s_k.                 (5)

Formula (3) is stated for k>=3; K^s_2=0. The first H^s values at
k=2,3,4,5 are 2,8,60,592. These count literal marked edges, not classes
modulo an unrecorded symmetry.

## 1. Four record families, before the owner quotient

At delta=2k, the exact deficit identity still has r=k, because a loss of
one backward edge costs at least 2k+1. Every endpoint is unit. The only
possibilities are

    A: d=0,h=0,U=2k;
    B: d=0,h=1,U=k;
    C: d=1,h=0,U=k-1;
    D: d=0,h=2,U=0.

Thus D, a matching with two empty gaps, is the only additional family.
A matching with an empty gap between two matched consecutive endpoints
is rejected: those vertices form a homogeneous pair. For one or two
empty gaps this is the complete clone obstruction. In D, if the holes
are adjacent, the remaining three consecutive endpoints have at most
one internal matching edge; an adjacent internal edge is exactly the
rejected case, the long edge makes C3, and no internal edge makes a
transitive triple. These exhaust the possible new close pairs. A and C
have all gaps occupied and remain pair free by the original argument.

Set m,E,J,P as in the corridor formula. The exact retained presentation
counts are

    A = E[binom(m+2k-1,2k)
          -m binom(m+k-1,k)+binom(m,2)],
    B = (mE-J)[binom(m+k-2,k)-(m-1)],
    C = P binom(m+k-3,k-1),
    D = E binom(m,2)-(m-1)J+K.                           (6)

For A, inclusion-exclusion now allows two parts at least k, each exactly
k, giving the new binom(m,2) term. For B, the residual deficit k has one
forbidden part of size k; there are m-1 choices. For C, no part can exceed
k-1. For D, if a matching has a adjacent matched pairs, its admissible
hole pairs number binom(m-a,2). The identity

    binom(m-a,2)=binom(m,2)-a(m-1)+binom(a,2)

and summation give the displayed D count.

A comparison with L_k(2k) shows

    A+B+C+D = L_k(2k)+K,                                (7)

because the other new terms cancel:
2E binom(m,2)-mE(m-1)=0 and the two J(m-1) terms cancel.
This cancellation concerns presentation counts only. The owner correction
still has to be paid; it must not be inferred from the scalar cancellation.

## 2. Complete ownership analysis

For k>=4, A has at most two unit gaps since 3(k-1)>2k. Its exact quotient
is the full-matching cell proof in `two_unit_cells.md`. B has at most one
unit gap, since its residual deficit is k<2(k-1); that note's missing/unit
cell analysis applies, even if another large gap has been lowered. Thus
A and B retain one distinct class per admissible record.

C has either no unit gap or exactly one. With none, every gap is large
and its complete P3-plus-matching graph is recovered intrinsically. With
one, there is one two-vertex exterior cell or one three-vertex interior
cell. All other unit vertices lie in singleton cells and are fixed by
the intrinsic large-vertex order. The internal cell cases are:

- No internal backward edge: a transitive cell fixes all its roles.
- An internal edge from an isolated matching component: an unmarked C3;
  its rotations change no record or outside edge.
- An internal edge of the P3: a C3 with exactly one distinguished vertex,
  namely the P3 centre incident with its other, external backward edge.
  That actual external port has two possible endpoint positions in the
  ordered endpoint/unit-gap presentation. These are two different C
  records of the same pointed C3.

The last two C records have exactly the D record constructed in section 5
of `two_unit_cells.md`: all three local vertices are matching endpoints,
the marked vertex is the middle one and has the external matching edge,
and the other two form the long internal matching edge. Both local gaps
are absent. The actual tournament, weights and external arc are identical.

D itself has no internal identifications. With separated interior holes,
the two two-endpoint cells have forward arcs, so their roles are fixed.
With adjacent interior holes, the three-endpoint cell is transitive, or
is the just-described pointed C3. The latter's external port must be its
middle endpoint, giving exactly one D record. Exterior holes create only
an exterior singleton or a transitive two-endpoint cell; their literal
outside incidences fix the endpoint roles.

These are also all C--D collisions. A C record with no unit gap has one
more large-weight vertex than D, so cannot collide. A C record with a
unit gap can match D's cell-size pattern only at an overlapping pair of
holes or at an exterior overlap. In the transitive interior three-cell,
D gives every vertex one distinct external matching edge, while C has
a gap vertex with no external backward incidence; the uniquely ordered
cell retains that difference. At an exterior two-cell, D gives both
vertices distinct external matching partners, whereas C has a gap vertex
with none and possibly a P3 endpoint with two. Again the ordered external
incidences differ. An unmarked C3 from an isolated matching edge has no
external matching port and cannot be D's marked C3. Therefore only the
marked C3 case collides.

A and B have different vertex numbers from C and D, so no other cross-
family collision is possible. The complete nontrivial owner fibre is
therefore precisely

    one D record + two C records -> one weighted class.  (8)

No arbitrary automorphism division is made. The local external port is
retained, and no clone is merged after the counts. Strongness is preserved
by the literal isomorphism and the existing proper-cut criterion, so this
same owner law applies inside the strong family.

## 3. Counting exactly the marked owner fibres

Every D owner fibre is specified by a matching edge (i,i+2), with the two
intervening gaps empty. Its middle endpoint has an external matching
partner automatically. Conversely every such marked edge gives one owner
fibre, as (8). Thus there are exactly H owner fibres, and (7) overcounts
each by two. This proves (1).

For all matchings there are 2k-2 possible positions for (i,i+2), followed
by F_{k-1} matchings on the remaining vertices. Hence H=(2k-2)F_{k-1}.
For K, two marked adjacent edges have disjoint endpoints: choose two
nonadjacent positions among the 2k-1 adjacency slots, in
binom(2k-2,2) ways, then match the remaining 2k-4 points. Thus

    K=binom(2k-2,2)F_{k-2}=(k-1)F_{k-1},

which proves (2).

For strong K, delete the two marked adjacent pairs. The remaining matching
is irreducible, by the earlier adjacent-pair deletion lemma applied twice.
Conversely, distribute two identical marked pair-blocks among the 2k-5
interior gaps of an irreducible matching on 2k-4 points. Both may occupy
the same gap. This gives binom(2k-4,2) insertion choices and proves (3).
At k=2 there are no such irreducible marked-pair choices.

For strong H, split an arbitrary matching at every closed cut. A marked
edge of span two, along with its intermediate vertex, lies entirely in
one irreducible block. Arbitrary matchings precede and follow that block.
The unique decomposition consequently gives the formal-series identity

    sum_{k>=2}(2k-2)F_{k-1} x^k
      = (sum_{k>=0}F_k x^k)^2 (sum_{k>=2}H^s_k x^k).

Taking coefficients proves (4), and hence (5). No analytic convergence
or unexplained sequence fit is used.

## 4. Controls and stopping boundary

`first_owner_check.py` completely generates the four record families at
k=2,3 and delta=2k. It retains every actual vertex, rejects and checks the
clone obstruction, canonicalizes all unit-cell permutations and counts
whole fibres. It tests that every nontrivial fibre is exactly (C,C,D),
with exactly H such fibres, rather than checking just the final sum.
The k=2 totals independently agree with the prior complete compiler.
The code also enumerates matching incidences through k=6 to test K,H and
the two strong formulas. A separate owner_four_check.py exhausts just the
C/D owner families at k=4: 40,605 general and 14,922 strong classes, with
90 and 60 three-record owner fibres respectively. This is not advertised
as a full k=4 mass-layer census. All of these use explicit runtime exceptions.

The k=2,3 all-unit exceptional cases lie outside the at-most-two-unit-gap
proof in section 2; their successful exact census is a FINITE-EXACT
extension, not a silently generalized proof. Normal and optimized replay
receipts identify the exact finite universe.

This resolves the first cross-family owner collision at delta=2k. No
claim, computation or extrapolation for delta>2k is made in this addendum.
