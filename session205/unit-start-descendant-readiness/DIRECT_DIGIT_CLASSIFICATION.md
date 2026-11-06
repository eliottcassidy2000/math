# Direct binary-ready counts in the O/Z/J unit-start grammar

Status: PROVED for the stated edge-seeded grammar, independently audited. This does not settle eventual scalar coverage.

Throughout this note, “reachable” means reachable from the rooted directed edge, whose state is (h,R,e,d)=(1,2,0,1), using only the following operations

- O: (h,2R,0,d).
- Z: (h+R-e,2R-d,d,d).
- J: (h+R-e,h+R,h-e,h).

This is an edge-seeded grammar theorem. It does not independently prove that this grammar exhausts every physical unit-start source. Any such interpretation requires the inherited complete-source-grammar theorem. The singleton source is excluded: its response is (1,1,1,1); including it would add h=1 to the minus-ready count set, but would not change the h>1 statements below.

Call a state minus-ready or plus-ready when R-e=h-1 or h+1, respectively. Either Z or J then gives the desired direct child 2h-1 or 2h+1. O alone preserves h.

## Classification

The set of counts having a minus-ready reachable representative is exactly

    {2^m+1 : m>=1} union {3*2^m+3 : m>=1}.

The set having a plus-ready reachable representative is exactly

    {1}
    union {3*2^t-1 : t>=1}
    union {5*2^t-1 : t>=1}
    union {3*(2^m+1)*2^t-1 : m,t>=1}.

Thus the direct one-letter binary-ready sets are sparse, even though the complete scalar support appears cofinite. This rules out an induction claiming that every sufficiently large represented parent can be replaced by a same-count representative admitting a direct Z/J child of each sign. The corollary below extends the obstruction to arbitrary descendant macros. It does not rule out target-dependent descent with a different scalar parent, or rebuilding outside this grammar.

## Invariants

Every state satisfies

    R>=h>=d+e, d>=e>=0, d odd, h odd,
    C:=R-h+e-d>=0.

These follow by simultaneous induction. Under Z, R'-h'=C; under J, R'-h'=e. The slack h-d-e becomes h-d under O, h+R-e-2d under Z, and R-h under J, all nonnegative using the previous inequalities. C is preserved by Z, reset to zero by J, and changes to C+R-e under O. Both h and d remain odd because R-e is even (also preserved by the recurrences).

Write E=h-2e. Then E>=d-e>=0 and is odd, so E>=1. The transformations relevant below are

    C(O(state)) = d+2C+E,
    E(O(state)) = h,
    E(Z(state)) = e+2x+y, where x=h-d-e and y=R-h,
    E(J(state)) = e+y.

Also h=d holds only for states with h=d=1: every Z/J creates h'>d', and O preserves h,d.

## All E=1 states

They are exactly the following three types:

    A_k = (1,2^(k+1),0,1), k>=0;
    B = (3,3,1,1);
    D_m = (2^(m+1)+1,2^(m+1)+2,2^m,2^m+1), m>=1.

To prove completeness, inspect the last letter.

- After O, E'=h=1 forces type A.
- After Z, e+2x+y=1 forces x=0 and either (e,y)=(0,1) or (1,0). The first gives h=d, hence the initial state (1,2,0,1), and its image is B. The second would give h=d+1, impossible since both h,d are odd.
- After J, e+y=1. If e=0,y=1, an input ending in O would have R>=2h, forcing h=1 and R=2, which is not an O output; thus the input is the initial state, giving B. If e=1,y=0, the input must have d=1. Indeed an e=1 state arising from Z has d=1; arising from J requires the previous h-e=1, which with h>=d+e forces d=1, and then parity and e<=d force previous h=1,e=0, so the new d is again 1.

The reachable states with e=d=1 and R=h are exactly (2^m+1,2^m+1,1,1), m>=1. To see this, a J producing d=1 must have h=1 before it and yields h=1+2^(k+1). From there Z maps h to 2h-1 and preserves this section. Any O makes C positive, which subsequent Z preserves; any J from h>1 makes d>1. The other entry is Z on the initial state. This proves the asserted section classification. Applying J to it produces exactly D_m.

All displayed states are attainable: A_k=O^k(initial); B=Z(initial); the canonical d=1 section is obtained by Z powers, and D_m is its J image.

## Readiness classification

Since R-e-h=d+C-2e, inspect the final operation producing a ready state.

An O output has e=0 and hence R-e-h=d+C>=3, so cannot be ready. The initial state is plus-ready with h=1.

For a Z output, readiness is C-d=+/-1. After the most recent J (or at initialization), C=0. At the first subsequent O, C becomes d+2C+E>=d+1; subsequent Z preserves it, and later O increases it strictly. Therefore C=d-1 is possible only when d=1,C=0. The resulting canonical Z outputs have h=2^m+1.

For C=d+1, there must be exactly one O since the most recent reset of C, applied to a C=0,E=1 state, then one or more Z operations. Among the E=1 states these possible origins are A_0, B, and D_m; A_k for k>0 has C>0. Starting at an origin (h0,R0,e0,d), applying O then Z^t, t>=1, gives

    h_t = h0+2d-2R0 + 2^t*(2R0-d).

The three origins respectively give 3*2^t-1, 5*2^t-1, and 3*(2^m+1)*2^t-1. All are attained.

For a J output, readiness difference is 2e-h=-E. Thus only minus-readiness is possible, precisely from the classified E=1 inputs. Their output counts are respectively 2^(k+1)+1, 5, and 3*2^m+3. This proves both exact sets.

## Computational scope

The inherited exhaustive pruning at bound M=50000 suffices to capture every ready parent h<25000: its required direct child is at most49999. A terminal O omitted by the production pruning cannot create a ready parent in this range, since the required following Z/J child would then also be in range. The finite lists agree with the formulas. This check is a regression control; the proof is the invariant and last-operation classification above.

## Corollary: arbitrary descendant macros do not enlarge the sets

For each sign sigma in {-1,+1} and every odd h>1, there exists a reachable state of count h having an O/Z/J descendant of count 2h+sigma if and only if h lies in the corresponding classified ready set. This is existential over the reachable representatives of h, not a claim about every representative. Fresh same-count representatives are explicitly allowed.

Indeed every non-O operation changes the count to

    h+R-e >= h+h-(h-1)/2 = (3h+1)/2.

Two non-O operations therefore give count at least (9h+5)/4, strictly greater than 2h+1. Any O before the sole non-O operation also prevents success: the state immediately after O has R-e-h=d+C>=3, and further O only increases this difference, so the following non-O count is at least 2h+3. Thus a successful macro consists of one immediate Z or J followed by zero or more O's. No-non-O macros preserve h, which differs from both desired values when h>1. The h=1 minus case is the trivial exception: O preserves 1=2*1-1.

Both ready sets have zero natural density: below X the minus set has O(log X) members, and the plus set has O((log X)^2). Therefore an eventual scalar occupancy proof for this grammar cannot proceed by a universal descendant binary step from all sufficiently large same-count parents, even with arbitrary fixed or target-dependent descendant macros. This is an obstruction to that proof strategy, not to the scalar occupancy statement itself.
