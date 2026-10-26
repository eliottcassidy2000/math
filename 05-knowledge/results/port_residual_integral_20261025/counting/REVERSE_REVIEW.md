# Independent review of the reverse valuation potential

2026-10-25. **REVIEW PASS:** the algebraic construction in
`REVERSE_POTENTIAL.md` is correct on the stated structurally reachable,
terminal-reaching ordinary transfer graph. **FINITE-EXACT:** its embedded
generic DAG control passes normally and with `python -O`, with identical
output. Integration into the production tournament graph remains **OPEN**;
this review does not promote that implementation claim.

## Structural graph and state merges

A vertex must contain the stage together with its exact (S,D) key. The
stage strictly decreases even on a zero-cycle skip, so the graph is
acyclic. Restricting to structurally reachable and terminal-reaching
vertices makes every suffix minimum finite. Early completed S=0 states
can be joined to a common terminal by zero-cost unit edges.

The existing ordinary defect quotient preserves every future continuation
with its normalized weight. Thus a prefix reaching a merged state can
be followed by any suffix available from that state. This is essential:
a quotient that merely preserves the final sum, without continuation
compatibility, would not by itself justify the path argument. Here stage,
remaining degree, and retained response coordinates provide precisely
that compatibility. No hidden cycle-count recovery is required.

The graph cannot be built by retaining only states whose prematurely
reduced absolute weights are nonzero. Structural reachability is the
correct input. Once the transformed recurrence and its precision are
fixed, ordinary modular evaluation is legitimate on that graph.

## Pathwise integrality and the scope of the lattice

For every edge v->w, the defining minimum implies
V(v)<=J_e+V(w), giving a nonnegative reduced cost. If a prefix path to v
has cost c, appending a minimizing suffix produces a complete path of
cost c+V(v)>=B. The transformed prefix term is therefore its positive
odd factor times2^[c+V(v)-B], an integer before any merging. Summing at
the same state uses the same V(v), so mergers commute with the change
of coordinates. Initial and terminal normalizations have the claimed
values. The proof does not infer integrality merely from a final integer
answer or from cancellation among prefixes.

This establishes an invariant reachable lattice, not a gauge on arbitrary
independent integer assignments to every Y_v. A useful deterministic
hostile is a four-vertex chain with edge costs0,2,1 and unit factors.
Then B=3, the third vertex has V=1 and reachable Y=4, and its transformed
value is1. Replacing that reachable value arbitrarily by Y=1 would give
1/4. Every reduced edge cost on the actual chain is0. The same example
shows both other boundaries: subtracting B only from the first edge
creates exponent-3, while its absolute terminal value8 is0 modulo2
despite normalized value1 modulo2.

No power-of-two division of an already truncated absolute state is used.
The transformed recurrence must be run from Z_initial=1. The extra reverse
pass computes costs from graph keys and local exponents; it does not
require the unknown tournament count or its full integer coefficients.

## Odd factors and identification of the baseline

For the ordinary recurrence beta=1, and each degree-local transfer factor
is an integer ratio of odd numbers, hence positive odd. Along a complete
path, degree-local odd factorials telescope to

    O(N)/product_d[d^k_d O(k_d)].

The moment valuation gauges telescope as well: the initial exponent is
zero because m=b=0, and the terminal exponent is zero because S=0.
The resulting contribution is O(N) times the original rational partition
summand. Since O(N) and the displayed product factor are odd, total edge
cost equals exactly the summand valuation, including multiplicity-factorial
denominators. There is no residual binary denominator omitted from B.

Consequently the separately audited termwise bound gives B=floor(N/2)
for N>=3. The single N-cycle attains it for odd N; the distinct coprime
pair1,N-1 attains it for even N>=4. Orders0,1,2 all have minimum cost0.
This B is the minimum summand valuation. Even when all contributions are
positive integers after scaling, binary carries among minimum-cost terms
can increase the valuation of their sum; B must not be identified with
v2(a(N)) in general. The order50 example correctly records this distinction.

At termination Z=O(N)*a(N)/2^B. Inverting O(N) modulo2^K uses no additional
binary digits. Thus the construction answers the normalized K-bit
question directly once the structural graph and V are available. It does
not contradict the earlier absolute-residue hostile, which uses a different
coordinate system before truncation.

## Replayed finite control and remaining implementation boundary

The embedded control was executed unchanged normally and under `-O`.
Both runs report204 finite acyclic graphs,2,040 prefix-coordinate checks,
1,224 normalized modular comparisons, and168 failures of naive initial
baseline subtraction. The mandatory consecutive edges ensure every
generated vertex is reachable and terminal-reaching. All gates use
explicit exceptions. The intermediate checks use exact Fractions, while
the transformed recurrence uses integer shifts.

These controls exercise the generic DAG identity, including actual
mergers and nonconstant suffix potentials. They do not construct the
ordinary tournament state graph or demonstrate an implementation speedup.
The document correctly leaves that integration and its cost as open work.
The ordinary positivity of arithmetic factors is also correctly kept
separate from an objectwise positive tournament model and from the signed
clone recurrence.

No proof gap or correction was found.
