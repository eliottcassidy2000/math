# Independent audit: dyadic sharpness and least period

2026-10-31. **Mathematical audit PASS** for `../ramsey/DYADIC_PERIOD.md`.
This is a separate automated mathematical review, not a human refereeing or
Lean claim. No unrestricted Ramsey bound is asserted.

## Exact scope

D contains nonempty graphs obtained from a singleton using binary joins,
disjoint unions, and full C5 substitution with five nonempty children.
For positive integer budgets A,B, F(A,B) is the largest possible order with
clique number at most A and independence number at most B. The inherited
G2 source in `ramsey_neighbourhood_lifts_20261030/GRAPH_CAPACITY.md` supplies
its exact slope w(A), upper bound, sharp constructions, and weighted C5
inequality. Session121 supplies eventual translation at L(A).

The conclusions audited here are:

1. F(A,B)=w(A)B exactly when L(A)=2^ceil(log2 A) divides B.
2. The least positive eventual translation period is exactly L(A).
3. The displayed recursive criterion is necessary and sufficient for sharp
   construction trees. It does not count graph isomorphism classes.

## Proof review

**Strict superadditivity.** The producer's parity induction is valid,
including the case x=y=1 and zero half-indices in the odd cases. A separate
verification uses D(x,y)=w(x+y)-w(x)-w(y), with D(x,0)=0:

- D(2a,2b)=(5/2)D(a,b).
- D(2a,2b+1)=(3/2)D(a,b)+D(a,b+1).
- D(2a+1,2b+1)=(3/2)D(a,b+1)+D(a+1,b)
  +(1/2)(w(b+1)-w(b)).

The last term is positive. In the even-odd case a>=1 and the second D-term
has two positive arguments with strictly smaller sum. Non-strict
superadditivity and strict monotonicity are inherited. Thus no induction
case uses the desired strict assertion circularly. Independently,
w(n+1)-w(n)=(3/2)^popcount(n), proved by the binary recurrence, verifies
monotonicity directly.

**Strictness at all integral polytope vertices.** The normalized feasible
independence-port polytope has eleven indicator vertices and one half
vertex. A two-entry indicator selects adjacent original C5 children, so
strict superadditivity makes its objective strictly below w(A). A singleton
indicator is strictly below because every child has clique number at most
A-1; this uses the nonempty-child hypothesis. The zero vertex is strict too.
The half vertex has objective at most w(A). Therefore if the upper bound
is attained, its maximizing face is the singleton half vector. The argument
does not infer unique optimization merely from the existence of one optimal
vertex: it explicitly excludes every other vertex.

**Equality propagation.** The upper bound is the sum of nonnegative child
slacks. Equality forces each child to be sharp, every independence port to
be B/2, and the clique-weight sum to be 2w(A). For A=2, B even already gives
the required lattice. For A>=3, writing k=ceil(log2 A), if all child clique
numbers were at most 2^(k-2), their weight sum would be at most
5w(2^(k-2))=2w(2^(k-1))<2w(A). Hence a child has clique number above that
threshold. Its own clique number is below A, so the outer induction is
legitimate. Its sharpness forces divisibility of B/2 by a power of two at
least 2^(k-1), hence 2^k divides B.

**Unions and budget slack.** Structural induction at fixed A handles union
children that may themselves have clique number A. Their actual independence
numbers add. Equality of the total bound forces equality of every child
bound and exact use of the total independence budget; thus divisibilities
add. Joins cannot be sharp by strict superadditivity. The singleton exception
is separately handled at A=1. Strict monotonicity ensures that a graph sharp
at the displayed budget actually has clique A and independence B.

**Least period.** The inherited finite L(A)-period tail exists. Any eventual
translation must have increment w(A)d by the limiting slope. Choose a
sufficiently large multiple B of L(A); its deficit is zero. Translation
preserves the zero deficit at B+d, and the exact sharpness lattice forces
L(A)|d. Thus no proper period, including one justified only by scalar
denominators, is possible. No bound on the start of the tail follows from
this proof, and none is claimed.

## Independent finite controls

`check.py` imports no producer or inherited checker. It builds w by the
binary digit-sum formula rather than the producer's recursive implementation.
It reconstructs the port-polytope vertices by exact Gaussian elimination of
all 252 five-active-constraint systems. It tests all positive feasible
clique patterns through A=12, indicator strictness, and the deep-child
witness on every sharp pattern in that finite universe.

A separate exact capacity dynamic program checks the sharpness lattice in
A<=6,B<=16. It is cross-checked against unpruned full pattern enumeration in
A<=4,B<=6. Monotonicity permits coordinatewise-maximal pattern reduction in
the larger check. Exact recursive construction statistics verify the inherited
attainment for A<=64 at B=L,2L,3L. These are graph-construction statistics,
not a census of literal graphs or isomorphism classes.

Reproduce:

    python3 -B check.py > result.json
    python3 -B -O check.py > result.optimized.json
    cmp result.json result.optimized.json

The finite controls are evidence about the implementations and proof
boundaries, not substitutes for the all-parameter argument above.
