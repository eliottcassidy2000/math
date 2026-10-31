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

### Completed dyadic controls

Normal and optimized outputs are byte-identical. They certify 511 binary
weight-recursion comparisons, 32,640 strict-superadditivity pairs, all 12
polytope vertices (11 indicators), 51,661 positive feasible clique patterns,
568,271 strict indicator tests, 71 sharp-pattern deep-child tests, 105,264
maximal-DP candidates, 55,131 unpruned-DP candidates, and 192 exact sharp
construction statistics. Integral-slope hostiles include F(3,2)=7<8 and
F(6,4)=35<40.

## Extension audit: complete rows five through eight

**Mathematical audit and independent exact replay PASS** for
`../ramsey/ROWS_5_TO_8.md` and `../ramsey/rows.normal.json`.
For B=8t+r, the claimed rows are respectively

- 62t + (0,5,12,19,29,35,43,51)[r];
- 80t + (0,6,15,23,35,43,55,65)[r];
- 98t + (0,7,17,27,41,51,65,79)[r];
- 125t + (0,8,20,32,50,62,80,98)[r].

They hold for every positive B, hence have threshold1 and least period8.

### Infinite-tail reasoning checked

1. Coordinatewise maximal positive clique patterns suffice because every
   other pattern can be increased finitely and child capacities are
   monotone in the clique budget. No positive child becomes empty.
2. For a nonsharp pattern, the *entire* port-polytope vertex maximum is
   strictly below w(A). The bound B<max_deficit/slope_gap correctly uses a
   strict inequality for a hypothetical violation. The integer cutoff is
   therefore ceil(max_deficit/slope_gap)-1. No positive dual is assumed for
   such a pattern; this avoids an actual negative-dual boundary in row8.
3. For a sharp pattern the independently solved dual is strictly positive.
   The exact identity subtracts its weighted nonnegative integer port
   slacks from w(A)B. A hypothetical violation has total slack at most
   S=ceil(max_deficit/min_dual)-1. Solving the odd-cycle port system yields
   2b_i=B-t_i with |t_i|<=S.
4. The downward reduction chooses B' congruent to B modulo8 at or above
   2+max t_i. The original B meets that bound because each b_i>=1. Thus
   B'<=B, all reconstructed ports remain positive integers, and
   B'<=S+9. Each sharp pattern's child translation law holds at *all positive
   arguments*, with period dividing4. This stronger property, not merely
   eventual periodicity, is required and is satisfied by the exact base
   rows1..4. The downward move therefore preserves a violation exactly.
5. Joins have a strictly positive slope gap. All 64 residue inequalities
   per row prove superadditivity of the proposed H, so structural induction
   covers union roots. The row is positive for all positive arguments and
   bounds the singleton.
6. Each of the eight finite seed constructions per row uses smaller clique
   budgets or the same clique budget at a smaller positive independence
   argument. Thus the construction recursion terminates. Actual child
   budget formulas, not mere objective values, establish feasibility.
   Disjoint unions with the B=8 seed establish every larger argument.

The finite reduction is therefore a genuine proof of the infinite tail,
not an observed-plateau rule. The explicit cutoff maxima are17,18,14,14.

### Independent executable replay

`check_rows.py` reads the producer receipt as untrusted data. It imports only
this auditor's own `check.py`, never producer code. It independently rebuilds
all maximal clique patterns, solves each sharp dual by Gaussian elimination,
checks all gaps and cutoffs, and enumerates the complete maximal
independence-port vectors. The producer instead optimizes with a max-plus
path DP, so the finite objective calculation uses a different algorithm.

All157 maximal clique patterns, all63 finite root-objective rows, all32 seed
construction plans and256 union residue inequalities pass. Normal and
optimized receipts are identical. `mutations.py` rejects six hostile edits:
a missing row, missing clique pattern, false dual, false cutoff, false finite
root value, and false seed order. Mutation results also agree under `-O`.

Reproduce:

    python3 -B check_rows.py ../ramsey/rows.normal.json > rows.result.json
    python3 -B -O check_rows.py ../ramsey/rows.normal.json > rows.result.optimized.json
    cmp rows.result.json rows.result.optimized.json
    python3 -B mutations.py ../ramsey/rows.normal.json > mutations.json
    python3 -B -O mutations.py ../ramsey/rows.normal.json > mutations.optimized.json
    cmp mutations.json mutations.optimized.json

The first four exact rows and slope theorem remain inherited dependencies.
No general Presburger compiler, all-row nonzero-residue classification,
unrestricted Ramsey improvement or graph-isomorphism classification follows.
