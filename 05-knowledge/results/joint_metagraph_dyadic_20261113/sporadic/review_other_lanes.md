# Independent review of the joint clock and native metagraph checker

**Verdict: PASS within the reviewed scope.** This is an independent
mathematical/code review by the sporadic lane; it does not promote this
lane's own RESERVED candidates. No producer module was imported. The reviewed
bytes are pinned in [review_other_lanes.json](review_other_lanes.json).

Reviewed targets:

- [pairings/PROOF.md](../pairings/PROOF.md), definitions and claims P1/P2.
- [audit/check_small.py](../audit/check_small.py), including its later rational
  column-rank, hidden-operator, and scalar-memory additions.

## 1. P1: the uniform clock is valid with its stated cutoff and onset

The nontrivial SCC inventory is finite for fixed H because a strong component
has at least3 Hamilton paths and at most H vertices: cutting a spanning
directed cycle gives distinct Hamilton paths. For factors at least3 their
sum is at most their product. Thus the all-order core sum has no unbounded
component-order variable hidden in its use.

For the converse word U,K,U^op, stripping singleton components in U leaves
t nontrivial components and t+1 padding slots. Formula (3)'s binomial is
exactly the stars-and-bars count for those slots; its parity restriction and
center types prevent double counting. The estimate2m+d<=h is valid also for
empty/singleton centers. Consequently n>=B is indeed a common safe onset.

For R>=1, ceil(log2(R+1))=1+floor(log2 R), including powers of two.
With L=k+floor(log2 R), every1<=j<=R lies below2^L, and
binom(2^L,j) has valuation L-v2(j)>=k. The binomial factor in the proof is
odd because2^L-1 has all its binary digits equal to1. Vandermonde therefore
gives the ordinary period without an unjustified polynomial-coefficient
division. The even step also preserves the converse-center parity; the
half-step works since t<=floor(R/2) and, when t>0,
floor(log2 floor(R/2))=floor(log2 R)-1. For R=1, t=0 and no logarithm of
zero is used. B=1,2 are explicitly separated. No failure was found.

The independent checker probes R=1..48, k=1..5, all allowable r,t, and
41 consecutive admissible upper arguments. These are controls for the
all-parameter proof, not replacements for it.

## 2. P2: reduced pole order gives the least eventual period

After numerator cancellation over F2, a positive denominator exponent e
means the reduced numerator is nonzero at1. Eventual period p is equivalent
to (1-x^p)A(x) being polynomial, since its sufficiently high coefficients
are exactly a_n-a_(n-p). For p=2^j q with q odd,1-x^p has multiplicity
exactly2^j at1. It cancels the sole pole iff2^j>=e. Therefore the least
period among all positive integers, not just powers of two, is the claimed
least sufficient power of two. The e=0 polynomial case has eventual-zero
tail and period1. Converse preserves H, so its nonfixed pairs prove ordinary
and fixed-locus parity agreement at every order.

The H9 numerator cancellation and the even H15 numerator are appropriate
hostiles to degree-only or support-only reasoning. The finite vector's
least period is the maximum of its component dyadic periods. No such
minimality is transferred to the sum of that vector, which could cancel.

Our independent test considers all2040 binary numerator/denominator cases
with denominator(1+x)^d, d=1..8, and numerator degree<=d+1. It computes
coefficients directly using binomial series, factors the numerator by
independent GF2 division, and tests all candidate periods1..16 on a safe
tail. Polynomial quotients, complete cancellation, zero series, and poles
on both sides of power-of-two boundaries are included.

## 3. The native checker enumerates the correct objects

The adjacent transpositions generate S_n. Its breadth-first traversal starts
from every previously unseen labelled orientation, so the orbit partition
is complete and its representatives are minimal masks. The automorphism
order n!/orbit_size is the correct orbit-stabilizer quotient. Every Hamilton
path is tested literally among all n! permutations, without a recurrence
import from another lane.

An easily misread line is `u=owner[t^flip]`, where flip contains the path
edges. The native complement operation flips the *nonpath* edges. These
two resulting tournaments are converses, so they give exactly the same
merged target used here. Reversing the path then gives the inverse incidence.
Thus the use of path-only flipping is valid specifically because the target
is the converse-merged metagraph, not the ordinary unmerged graph.

Dividing each path incidence by Aut(T) correctly passes to path orbits:
an automorphism fixing an ordered spanning path fixes every vertex. Summing
over both ordinary classes of a nonself-converse pair explains its factor2.
The blue test asks whether reversal of that very path is an anti-automorphism,
so it retains the actual frame rather than merely testing whether T happens
to be self-converse. The symmetry, integrality, doubled-loop, and total
tiling-mass gates are independent consistency controls.

The folded-score partition takes the lesser of the score sequence and its
converse, exactly as needed after converse merging. Every degree is positive.
Writing D=diag(d_i), P=D^-1 A, and block masses M_a, the code computes
conditional expectation with weights d_i and verifies

    M_a [(E P^2 L)_(ab) - (E P L)^2_(ab)]
      = sum_j d_j (p_j(a)-Q_(block(j),a))
                    (p_j(b)-Q_(block(j),b)).

Reversibility follows from A=A^T. The diagonal is a weighted sum of squares,
so all its diagonal entries vanish iff every within-block response is
constant. This justifies both the positive return defect and the exact
lumpability equivalence tested by the program.

## 4. The added linear-algebra checks are correctly typed

`column_rank` performs rational pivot elimination with normalized pivot
entries. Insertion of an earlier pivot does not invalidate older basis
columns because their earlier entries are zero. The accumulated columns
are L,PL,...; once their rank stops growing, their span is P-invariant, so
the early stopping criterion is sound. At orders5,6 the reported full ranks
10,34 certify reachable dimension. Turning that into a *minimal moment
realization* uses the reversible weighted Gram/Hankel argument in the
metagraph proof, rather than reachability alone; the proof supplies it and
excludes a separate direct-feedthrough term.

At order5, each hidden basis vector has one unit nonanchor coordinate and
a compensating last coordinate -d_i/d_anchor. Hence it is weighted-mean-zero,
and the selected nonanchor coordinates recover its exact expansion. Applying
the hidden projection after P gives the correct3-by-3 matrix C. The trace,
principal-minor sum, and determinant give

    det(I-zC)=1-(6590/186615)z-(16714/186615)z^2
                      +(2126/186615)z^3.

The transitive block is a singleton. The initial hidden vector is
(I-LE)P e_T; each iteration applies C, and the observation is the T coordinate
of P times that vector. Thus the stored moments are the stated scalar memory
kernel, with constant1/18. Cayley-Hamilton proves the displayed recurrence
at every index; testing seven further identities is a supplemental control.

The root receipt inspected before the explicit Hankel gate had137379 active
checks. Its first five scalar
memory moments additionally have3-by-3 Hankel determinant
**1/2059408867230**, computed from the receipt using independent rational
arithmetic. Thus that particular scalar response also needs dimension3;
the cubic is not merely an unreduced possible denominator. The final reviewed
source also adds that determinant as an explicit gate and output field; its
Leibniz formula uses precisely the moments indexed i+p(i), so it is correct.
The receipt is
root-produced data; this review is not a third independent full metagraph
enumeration.

## 5. Receipts and remaining limitations

`python3 check_other_lanes.py` and `python3 -O check_other_lanes.py` produced
byte-identical [normal](review_other_lanes.json) and
[optimized](review_other_lanes.optimized.json) receipts with382887 active
checks. This reviewer also independently verifies the Lorentzian Gram
matrix [[0,-1],[-1,0]] and the explicit projection used in the recovered
23/24/25/26 bridge; that check is separate from P1/P2.

No defect was found in the reviewed statements. The review does not audit
pairings P3/P4's external strong-inventory cutoff, the metagraph lane's
all-order non-lumpability construction, or a complete tournament-generation
recurrence. It certifies the stated finite checker logic and P1/P2's
all-order modular claims, with their existing precision, cutoff, and onset.
