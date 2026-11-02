# Independent final-refinement audit

**ACCEPTED:** both separately frozen producer refinements pinned in
`REVIEWED_PRODUCERS.json`, and the exact synthesis in `CONNECTION_REVIEW.json`.
Older independent audit directories are unchanged. This directory contains
all inputs needed for its normal/-O replay and hash validation.

## Connected Erdős extremizers

The translation argument is valid for every n>=1. Move one whole component
A while keeping its complement fixed; the translations making a chosen cross
pair unit form a circle. Only finitely many translations collide vertices.
A permitted translation preserves all original unit pairs (none joined A to
its complement) and adds a new pair, contradicting extremality. This concerns
a final placement, not collision-free motion. For n>=3, triangle plus one-unit
attachments gives u(n)>=n. Connectedness then implies cyclicity, so epsilon
=u(n)-n and beta=u(n)-n+1 exactly. The n2 exception is genuine. Dated interval
tables are arithmetic consequences of the already cited bounds, not newly
proved global extrema.

## Reversal and exact binary content

Integer substitution z=w+1 and its inverse give an invertible integral basis
change. Thus gcd of ordinary coefficients equals gcd of binomial moments;
M0=2^N forces the gcd to be a power of two. This is not a modular division.
For n>=2 and a nonempty reversal-closed bank, every nonempty compatible family
is paired with a distinct reversal: a fixed family would contain a path and
its reverse, prescribing opposite arcs. Union size and family size agree in
the pair. Umax>=1 also includes the empty contribution in the factor
2^(N-Umax+1). A directed non-reversal-closed bank and n1 are genuine hostiles.
The leaf-half doubling explicitly requires j>=1; the empty family is not
doubled. The source proof now states this quantifier.

Independent complete path/clique enumeration proves content2^276 at n25
and2^296 at n26, both by ordinary coefficient gcd and binomial-moment gcd.
The n25 reversal lower factor2^275 is not sharp; four compatible pairs
supply one additional power. At n26 the max union30 occurs among pairs,
not triples; it makes the2^296 bound sharp. Primitive coefficients and the
exact means/probabilities all match. A larger maximum response need not
mean a larger expectation or success probability.

The exact valuation increase is20. Its decomposition is
(325-300)-(30-26)+(0-1)=25-4-1=20, comparing the two spectrum contents.
This is not the incorrect radius subtraction:704-694 remains10.

## Common-block transfer and the indispensable orientation distinction

For a nonempty compatible Hamilton-path family on n>=1 vertices, its arc
union has a connected underlying simple graph. With b common directed blocks
and J ordered joins, m=n-b+|J|, hence beta=|J|-b+1 and the REORIENTABLE
pseudoforest deletion defect is epsilon=max(0,|J|-b). The tree/one-path case
b1,J0 is included. Contract internal block paths without simplifying the
quotient: opposite nontrivial-block joins become distinct parallel edges.
Dropping that multiplicity would lose cycle rank and deletion cost. No
Euclidean unit-distance realization of the square-sum graph is asserted.

If original arrows must be retained, each non-sink can independently retain
one outgoing arc. Therefore the exact fixed-arrow deletion budget is
sum_v max(0,dplus(v)-1)=m-n+z, where z counts sinks. Since the union contains
a Hamilton path, z<=1, with z1 precisely when all paths share a terminal.
A cyclic common-terminal union costs ONE MORE deletion than the reorientable
budget. These two budgets must not be conflated with arithmetic dynamics.

Reversal swaps zero-indegree and zero-outdegree counts. Thus
epsilon_arrow(reverse D)-epsilon_arrow(D)=z_in-z_out, with each count0or1
and corresponding respectively to a common start or a common end. Reversal
changes this budget by at most one while preserving beta and epsilon.

The specifically displayed n25 pair and n26 triple STARTING at18 have
(b,|J|)=(3,4),(5,8), z0, and budgets1,3. Their reversed families END at18,
so z1 and their fixed-arrow budgets are2,4, while reorientable budgets remain
1,3. This distinction was caught by the independent all-start enumeration.
The small compatible paths(1,2,3,4,5),(3,4,1,2,5) have n5,m6,z1, yielding
reorientable budget1 but fixed-arrow budget2. The exact control is retained
in normal.json and under -O.

## Verification scope

`python3 -B verify.py` validates this manifest and producer proof snapshots,
then runs the independent checker normally and optimized from another cwd.
No producer Python code is imported. The connected-extremizer argument is a
proof review, not a new global extremal census. The finite controls do not
replace the general proofs, and no broader conjecture or novelty claim is
certified.
