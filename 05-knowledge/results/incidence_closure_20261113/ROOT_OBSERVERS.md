# Three different meanings of a closed observer

**Status: PROVED IN NOTES, O1–O3, after independent review; see [AUDIT.md](AUDIT.md).**
Finite controls are in `audit/independent.py`; there is no all-order tournament
separation claim. The inheritance and hostile examples are in `BOARD.md` and
`closure/PROOF.md`. This note distinguishes three minimality objectives that
must not be interchanged.

## O1. The algebra retained by an exact finite-state observer

Let X be a finite set, let M_c be specified real matrices on functions on X,
and let V_0 be exactly the block-constant algebra of an initial partition Pi_0. The smallest
unital pointwise algebra containing V_0 and invariant under every M_c exists:
intersect all such algebras, including R^X. It is exactly the block-constant
algebra of the coarsest stable refinement in `closure/PROOF.md`, Section 2.

To prove the only nontrivial identification, take any unital subalgebra A of
R^X and declare x equivalent to y iff every f in A has f(x)=f(y). For each
equivalence class C and each other class D, choose f_D separating C and D.
The product of the factors

    (f_D - f_D(D))/(f_D(C) - f_D(D))

is 1 on C and 0 on every other class. Each f_D is constant on each class
by definition. The resulting class indicator belongs to A, so A is precisely
the block-constant algebra. Matrix invariance is equivalent to constancy of
the sums of each row over each target class. This proves O1.

Thus partition closure is a closure operator: monotone, extensive in the
information order, and idempotent. In the poset of observation algebras, its
universal property is the least invariant algebra above the initial one.
One can describe this as a reflector into the invariant algebras. This is a
statement about a poset category; it does not turn an arbitrary graph quotient
into a Grothendieck site.

The invariant **linear** span of all operator words applied to V_0 is another
minimum. It need not be a pointwise algebra. For example, on three points the
span of 1 and (0,1,2) is invariant under the identity, but does not contain
(0,1,4), the pointwise square. A minimal linear realization can therefore
retain different information from a minimal state partition. Neither minimum
is the number of allowed experimental columns: that is a separate set-cover
problem after the admissible measurements have been specified.

## O2. A second-order diagnostic is not a two-pass repair theorem

Let P be reversible with positive stationary weights D, let L lift partition
functions, put M=L^T D L, E=M^(-1)L^T D, and R=I-LE. Then R is the
D-orthogonal projection onto the hidden coordinates. With B=RPL,

    M(EP^2L - (EPL)^2) = B^T D B.

Indeed the difference is EPRPL; multiplying by M and using DP=P^T D and
R^T D=DR gives the displayed Gram matrix. It vanishes iff B=0, precisely
strong lumpability. This is the inherited exact memory theorem of the prior
packet, with a useful continuous-time form. If G=P-I and F(t)=E exp(tG)L,

    M(F''(0) - F'(0)^2) = B^T D B.

The same calculation applies because RGL=RPL. No factor 1/2 belongs here:
this is the derivative, not its Taylor-series coefficient. Consequently a
nonclosed observer is detected at derivative order two, even if repairing
that observer takes many refinement steps.

For an explicit reversible hostile, let P be simple random walk on the path
0--1--...--m, m>=2, and observe {0} versus {1,...,m}. After k strict
refinements the classes are {0},...,{k},{k+1,...,m}; at k=m-1 the
partition is discrete. The induction follows because the first point of
the remaining tail alone has positive probability into the preceding
singleton. All other tail points send total probability 1 to the tail.
Thus there are m-1 strict rounds, with no universal two-pass bound.

Nevertheless its omitted-state defect is already nonzero at two steps.
The actual return probability from 0 is 1/2; the degree-weighted quotient
predicts 1/(2m-1). This example is checked exactly for m=2,...,14, while
the preceding induction proves it for every m>=2. Compare the genuine
two-plus-pass sheafification theorem in `bridges/PROOF.md`: its bound follows
from site axioms and matching families, not from reversibility or a derivative.

## O3. What recurrence is actually closed

For the stable partition, M_c L=L B_c. Hence every word of operators descends:

    M_c1 ... M_ck L = L B_c1 ... B_ck.

For one operator, the formal generating matrix is (I-zB)^(-1). This closes
walk enumeration at the fixed order, including specified colors and endpoints.
It does not construct the matrices for order n+1, produce new tournament
isomorphism classes, or establish that every allowable Hamiltonian count is
realized. The empirical discreteness at orders 5,6,7 rules out a smaller strongly
lumpable partition quotient refining the specified initial observer. It is
not a lower bound on every linear realization of selected word statistics
or every conceivable A000568 recurrence.

The exact scalar earlier-value recurrence already recorded in
`../joint_metagraph_dyadic_20261113/README.md`, final section, remains available with
its independently computed cycle-inventory coefficients. A compact formula
for those coefficients, or a growth rule retaining joint incidence and
realization, is still an additional target.
