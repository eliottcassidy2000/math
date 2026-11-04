# Why adding a four-vertex atom needs more than one shared port

2026-11-03. **PROVED scoped obstruction and exact two-port hostile.**
This concerns induced grafts, not arbitrary rewrites or native complements.

The recovered mechanism is the one-vertex merge inequality in
[THM-4745, exact gluings of interval chains](../../../01-canon/theorems/THM-4745-exact-gluings-of-interval-chains-and-the-reduction-of-busch-equality-to-a-one-way-lemma.md),
statement 1. If two induced subtournaments A and B cover T and meet in one
vertex, then H(T)>=H(A)H(B). The orientations between the two remaining
vertex sets are arbitrary. A pair of Hamilton paths splits at the shared
vertex; merge the two initial paths and the two final paths while preserving
their orders. Restriction of the chosen merged path recovers the pair, so
this is an injection. No strongness hypothesis is needed for the inequality.

Consequently, grafting the strong four-vertex tournament Q onto a source S
through exactly one shared vertex, while keeping both induced graphs,
forces H(out)>=5H(S). For H(S)>=3 this exceeds both 3H(S)+2 and 3H(S)+4.
Thus a fresh C4 with a single shared port cannot itself be either desired
digit gadget. Extra arc changes can invalidate the induced-source premise;
two or more shared vertices have a different boundary law.

The sharp hostile to extending the inequality to two shared vertices is Q
itself. Use rows [8,1,3,6]. Its induced directed triangles on {0,1,3} and
{0,2,3} share the directed arc 0->3. Each has H=3, while their union has H=5,
not at least 9. Their separate path choices are not independently mergeable
once both shared endpoints must occur consistently. The required sidecar is
their joint endpoint order, not just the two marginal counts.

This explains a useful operational distinction. The inherited C2 operation
grafts a triangle through one old auxiliary and has the exact excess law

    H(out)-3H(source)=2*sigma.

The start defect sigma measures its surplus over the merge lower bound.
The operation moves to a new zero-start-defect auxiliary, so repeating it
there attains the lower bound and gives the zero digit. A four-vertex
replacement must either overlap more of the old boundary or alter its old
relations if it is to avoid the factor-five obstruction.

The new homogeneous-block exclusion is stricter in a different direction:
even a larger fixed exterior around the unchanged TT_t block cannot realize
the missing uniform 24/44 families. Overlapping atoms therefore remain only
a candidate operation, not a solution. A successful construction must also
escape that fixed-block hypothesis, for example by promoting variable
prefix vertices or changing its reader structure with the clock.

`check_overlap.py` independently lists the Hamilton words of Q and both
triangles and checks the graft inequality on every way to join two directed
triangles through one vertex. This is a positive control for the proved
one-port statement and a minimal two-port failure, not an all-graph census.
