# Independent square-lane audit

**Status: VERIFIED independent proof audit, 2026-11-12.** Reviewer:
cluster_leech_ports agent. Source: [squares/PROOF.md](../squares/PROOF.md),
sections 3, 7, and 8. No external priority judgment is made.

The multi-arm coefficient formula is correct for every declared arm length
and residual-demand boundary. The internal degree-one equations force
alternation; odd arms contribute (0,0)/(1,1), even arms contribute
(1,0)/(0,1). Solving a=i+j, b=i+e−j gives exactly the two binomial indices.
Length-one arms are handled by an empty even-position product, and the
simple-graph restriction prevents a false length-two cycle from two
parallel length-one arms. The exchange graph is exactly the Cartesian
product of the two Johnson graphs under the defined same-parity exchanges.
Its vertex count, degree, diameter, and connectedness follow from the
explicit subset coordinates. No connectivity assertion for the Hamiltonian
induced subgraph is smuggled into this degree-sector result.

The endpoint-relocation criterion is necessary and sufficient. A finite
simple support with exactly two degree-one vertices and all other degrees
two consists of their path plus disjoint cycles. Removing the marked edge
opens its cycle, or splits the endpoint path. The first case succeeds iff
exactly one cycle was present. In the second case, orienting the old path
from b to t distinguishes the successful a-before-v splice from the
opposite order, which closes a new cycle. Any other old cycle survives and
is fatal. The stated distinctness hypotheses remove degenerate endpoint or
loop cases. The inverse deletion of z and reinsertion of a–v restores all
declared degrees and is unique, so the conditioned count identity is a
bijection rather than merely an injection.

The arithmetic-host family is valid: v+z=(m+1)^2 and a+v=b+z=m^2 by direct
substitution. The inequality 2r<m^2−2m−3 simultaneously guarantees a>b and
v<z; b<v and a<z follow from their differences 2m+1. The explicit a≠v
condition handles the remaining coincidence. The r=0 subfamily excludes
a=v because (m−2)^2=8 has no integer solution. The actual m=6,r=9 values
are (a,b,v,z)=(13,10,23,26). The source-support existence caveat is
necessary and is expressly retained: these are infinitely many graph
hosts, not infinitely many proved Hamiltonian paths.

No mathematical repair was needed in the audited sections. The exact
finite controls are corroboration; the all-size proof does not depend on
extrapolation from their bounds.

## Directed transfer addendum

Section 9 also passes independent proof audit. The allowed old arc a→v is
explicitly retained, which is essential for the inverse. Removing that arc
from a cycle produces v→…→a; inserting a→z→b joins it to b→…→t. Exactly
one old cycle is necessary and sufficient. If the old arc is on the
endpoint path, the new arcs instead close b→…→a→z→b and cannot give a
Hamiltonian path. The new start is v and the end remains t. The inverse
of any target path in the selected (a,z,b) sector cuts its unique order
v→…→a→z→b→…→t into exactly one cycle and one path. No directed segment
reversal is used. The subset convolution counts each source once because
its cycle vertex set is unique.

I independently reconstructed all eight stated six-vertex tournament
extensions and counted H by subset DP, and the conditioned sector by
all 720 vertex permutations for each extension. Both exactly reproduce
the ordered totals 9,3,17,5,17,5,23,5 and conditioned counts all one. This
hostile example correctly shows why the conditional identity cannot be
promoted to control of the total tournament path count.
