# Incoming work and the next useful composition

Reviewed through origin/main `17b4acad2e` before publication. These incoming
computations were read, not rerun or claimed as new here.

The [apex-wedge constructor](../apex_wedge_constructor_20261103/README.md)
constructs every even H edge jump, retains literal sources and ports, and
has an independently proved six-coordinate jump carrier. Its
[PROFILE](../apex_wedge_constructor_20261103/PROFILE.md) lifts the edge
difference to all bin counts and derives a Lucas recurrence with roots
m+1 and m-1. Our six-state tail-word automaton has a different state and
operation: it remembers the last before-tail and first after-tail guards.
The useful common object is the full bin polynomial, not the number six.

There is an immediate exact way to transport a *rooted* edge difference
through the new seam. For fixed Y put u_Y=h_Y+F_Y. Direct subtraction of
CALCULUS C1 gives

    Delta H_(X diamond Y)
       = (u_Y-s_Y) Delta h_X + u_Y Delta F_X
           - (h_Y-s_Y) Delta s_X.

This is an algebraic corollary of the proved seam formula. The incoming
unrooted profile supplies Delta h, but transporting a marked attachment
also needs Delta E_X(2)=Delta F and Delta S_X(1)=Delta s at the same root.
On two sources in the unique-start class, Delta s=0. That is a concrete
possible bridge between edge-jump constructors and the new recursive
classification; source membership and endpoint paths must still be kept.

[Session160's adaptive codes](../session160_ordered_programs_20261103/adaptive_codes/NOTE.md)
encode integers injectively into H values or other statistics. Our rank
instead encodes every rooted class of an exactly classified family. Its
decoder reconstructs the full actual tournament, while its rank need not
be H. The two coding tasks are compatible. A promising next operation is
to carry a source-changing program together with this canonical class
address while it remains in the family, and expand to joint-port
polynomials when an exit leaves it.

The earlier incoming source-selection polynomial carrier and Session159
shift operator are credited in README. Their arbitrary polarity and
triangle-substitution scopes are retained; our sole-new-root formulas
and unique-start grammar do not replace their larger port carriers.
