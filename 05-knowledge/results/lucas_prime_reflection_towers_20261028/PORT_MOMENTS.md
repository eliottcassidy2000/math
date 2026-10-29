# All ordered-port counts and an exact critical-tree law

2026-10-28. **PROVED** from the inherited wreath actions and elementary
Burnside counting. Same-agent audit; the branching-process interpretation is
an identification of the stated distribution, not a transfer to Collatz.

## 1. The next compatibility target has an all-order recurrence

Retain A_k=Aut(W_k) and C_k=C_(A_k)(t_k) from `FULL_HOLONOMY.md`.
Let a_(k,r) and c_(k,r) count their orbits on **ordered r-tuples of vertices**,
allowing repeated vertices. In particular r=0 has one orbit and both groups
at k=0 act on a singleton, so

    a_(k,0)=c_(k,0)=1,       a_(0,r)=c_(0,r)=1.

For r>=1 the exact recurrences are

    a_(k,r) = (1/3) sum_(i+j+l=r) multinomial(r;i,j,l)
                        a_(k-1,i) a_(k-1,j) a_(k-1,l),

    c_(k,r) = sum_(s=0)^r binomial(r,s) 2^(r-s)
                        c_(k-1,s) a_(k-1,r-s).

The first sum distributes the named tuple positions among the three top
blocks. Cyclic rotation acts freely on every nonempty assignment of positions
to blocks, and the child groups act independently. Divide by3.

For the second, choose s positions in block0; the rest choose block1 or2.
C_(k-1) acts on the former. The latter two blocks share one A_(k-1) action,
conjugated by t_(k-1) on block2. Applying that conjugation to the block2
coordinates makes it a single diagonal action on all r-s remaining positions.
Their combined orbit count is a_(k-1,r-s), **not a product of two marginal
counts**. This is the precise shared-port compatibility term.

The fixed-r computation uses O(k r^3) integer arithmetic operations if all
counts through r and depth k are retained. Integer bit costs are additional.
The recurrence counts orbits; by itself it does not select canonical tuple
representatives or compute weighted responses on those orbits.

## 2. A polynomial hierarchy instead of enumerating every group element

For fixed r>=1, a_(k,r) is a polynomial in k of degree r-1 with leading
coefficient r!, and c_(k,r) is a polynomial of degree r with leading
coefficient2^r r!. In particular,

    a_(k,1)=1,             a_(k,2)=2k+1,
    a_(k,3)=6k^2+2k+1,
    c_(k,1)=2k+1,          c_(k,2)=8k^2+1.

Proof by induction on r: in the a recurrence the three terms with all r
positions in one block reproduce a_(k-1,r). Every other term uses smaller
moments. The terms with exactly two nonempty blocks have degree r-2 and
combined leading coefficient(r-1)r!; the three-nonempty-block terms have
smaller degree. Discrete summation gives degree r-1 and leading r!.
The case r=1 is constant1.

In the c recurrence, remove the s=r term. Every remaining term has degree
r-1 by the preceding result and induction on c. For each s=0,...,r-1 its
leading coefficient is2^r r!, including s=0. Summing r such terms and
discretely integrating gives2^r r! k^r. This also proves polynomiality,
not merely an asymptotic bound.

Thus the exponentially large centralizer has a polynomial number of orbit
types for each fixed number of ordered ports. The degree grows with the
number of jointly retained ports. There is no uniform constant-size quotient.

For two ports, independent single-port classes give(2k+1)^2 labels, whereas
the actual joint count is8k^2+1. The missing refinement has

    (8k^2+1)-(2k+1)^2 = 4k(k-1)

additional cells. At orders9 and27 the joint counts are33 and73, compared
with marginal products25 and49. This is an exact count of the failure of
independent port compression. For **distinct** ordered ports subtract the
2k+1 diagonal orbits, obtaining8k^2-2k.

## 3. Why the counts are moments of a critical tree

For a uniformly chosen element of A_k, let X_k be its number of fixed
vertices. A nonzero macro rotation has no fixed vertex. With probability1/3
the macro rotation is zero and the three uniform child automorphisms are
independent. Therefore, in distribution,

    X_0=1;
    X_k=0 with probability2/3,
        X_(k-1)^(1)+X_(k-1)^(2)+X_(k-1)^(3) with probability1/3.

This is exactly the population size at generation k of a Galton-Watson tree
with0 offspring with probability2/3 and3 offspring with probability1/3.
Its mean offspring count is1 and variance2. The probability-generating
function is iterated by f(s)=(2+s^3)/3. Elementary conditioning gives

    E X_k=1,         Var(X_k)=2k.

For a uniform element of C_k, let Y_k be its fixed-vertex count. Its linked
child actions have the same fixed count, not independent counts:

    Y_k = Y_(k-1) + 2 X_(k-1),     Y_0=1,

with the two displayed summands independent. Iterating uses independent
copies of X_j at the different depths. Thus

    E Y_k=2k+1,      Var(Y_k)=4k(k-1),
    E Y_k^2=8k^2+1.

Burnside says a_(k,r)=E X_k^r and c_(k,r)=E Y_k^r, because an ordered
r-tuple is fixed exactly when each entry is fixed. This proves the moment
interpretation with a specified probability law: uniform group elements.
The variance4k(k-1) is **exactly** the extra joint-port count above the
product of the single-port counts. No decorrelation assumption is made.

The growth of higher moments explains a genuine rare-large-fixed-set effect
in these uniform symmetry samples. It does not establish arithmetic heavy
tails or cancellation rates for actual Collatz sources. Those have different
sampling laws and clocks. The precise common research move is to derive
the joint source law before using its moments.

## 4. The small one-bit and one-trit controls fit the same action

At order9 the residual C3 acts on paired blocks by(j,l)->(j+a,l-a),
leaving j+l mod3. On two distinct ports in the same block it leaves their
nonzero difference, one of two values. On all three distinct ordered ports
in one block it has exactly two orbits, the two cyclic orientation classes.

This recovers the **type of joint invariant** in the incoming Session92
[six-order witness](../session92_joint_order_carriers_20261029/runs/RECOVERED_RUN_WITNESS.md):
three vertices individually lie in one C3 orbit, but their six orders split
into two classes and their native-complement Hamilton counts differ. The
underlying action on those three vertices is literally the same rotation
action. The surrounding thirteen-vertex tournament and its Hamilton counts
are not identified with W_2.

The port recurrence now closes the counting question left open by the
single-port calculation. Remaining targets concern canonical representative
construction, weighted native operations on these tuple classes, and the
largest operation menu under which their polynomial fixed-r size survives.
