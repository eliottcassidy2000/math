# Multi-return seams, physical profile-neutral trades, and a rooted dipole

11 November 2026. **PROVED locally:** the reverse-matching seam formula,
its complete literal-cycle decoder, its square-zero determinant expression,
and the explicit six-vertex selector exchange. **FINITE-EXACT:** the stated
small constructor tests. **OPEN:** arbitrary target realization and Chan–Pak.
No new canon ID or literature-priority claim is made.

## 1. The constructive gain

Session263's one return edge creates a clique of new cycle addresses. Several
return edges create both compatible new cycles and new cycles using multiple
returns. Those latter cycles need not be treated as uncontrolled leakage.
Here they have an exact cyclic-word decoder and a finite matrix formula.

The first useful output is already concrete. Two different ways to reverse
two edges between two triangles give modules A and B with

    Z_A(t,u) = 1 + 6tu^3 + 4tu^5 + 2t^2u^6,
    Z_B(t,u) = 1 + 6tu^3 + 6tu^5 + t^2u^6.

Their complete unmarked ordered-cover profiles agree:

    Q_r(A) = Q_r(B) = r^6 + 12r^4 + 16r^2.

Yet their exact selector graphs are respectively K10 minus a two-edge
matching and K12 minus one edge. At the same marked vertex they also agree
on the full deletion profile and on both alternating endpoint responses,
while their two ordinary endpoint responses move oppositely by two. This is
a **physical rooted dipole**, not just two tables with the same scalar H.
The [shared-root fusion theorem](../fusion/PROOF.md) uses it as an actual
sorting input: exchanging BA for AB strictly raises H in every physical
left/right fusion context, giving ab+1 different counts from a copies of A
and b copies of B. This cross-lane consequence is proved there.

The formula below is a compiler for these objects and for arbitrary numbers
of matched returns. It retains literal paths, vertex exclusions, cycle
addresses, and all packing orders. It does not pretend that ordinary numeric
matrix multiplication remembers disjointness.

## 2. Inheritance and the deliberately unusual transfers

Closest proved mechanism: Session263 RJ1/RJ2 in
`../../session263_spin_flip_20261111/selectors/RETURN_JOIN_SELECTORS.md` parses a cross cycle at its unique return
edge and retains induced-deletion packing data. Canonical hostile: replacing
literal directed cycles by vertex supports merges distinct addresses. The
corrected near miss in this lane is multiplying arbitrary endpoint-path
counts without deleting the other reversed forward seams. The least-used
sidecar is a commuting square-zero variable for each actual vertex.

Anchor: construct new physical selector gates. Niche: two-return profile
trades invisible to every unmarked cover count. Wildcard: convert a hard-core
partition function into a finite logarithmic determinant over a ring that
remembers every occupied vertex.

Live concept board: return-word parsing; occupied vertices; odd length;
packing components; profile-neutral trades; directed rooted response.
The old sources sampled, their exact paths, and what was transferred are
recorded in `SOURCES.md`. In brief:

- The independence-polynomial viewpoint supplies exclusion by intersection,
  but literal cycles remain addresses, including cycles on the same support.
- GMC first-return work suggests cutting at return events. Unlike its signed
  scalar moments, our objects have a positive literal parse before any
  aggregation; we do not infer atom isolation from an equation of numbers.
- Paley series/parallel kernels suggest choosing a coefficient algebra
  closed under composition. Here that algebra retains occupied vertices,
  rather than quadratic-character terminal differences.
- Spine automata suggest indexing a seam by its outstanding ports. The new
  carrier is a matrix of actual path legs with square-zero vertex memory,
  not a claim to improve the already known finite-cutwidth transfer theorem.
- Electrical owner-code work highlights why a small observable is not an
  allocation. Our determinant retains vertex-use capacity algebraically;
  it is not an electrical flow theorem or a polynomial-time capacity solver.

## 3. Literal reverse-matching construction

Let A and B be disjoint tournaments. Choose distinct vertices

    alpha_1,...,alpha_k in A,  beta_1,...,beta_k in B.

Orient all cross edges A -> B, then reverse precisely

    beta_i -> alpha_i,  1 <= i <= k.

Call the result R_k(A,B). Both shores and all their internal edges survive
literally. The distinctness condition says that the reversed cross edges are
a matching. It is part of the theorem, not an implicit genericity assumption.

Every cross cycle uses some nonempty set of these returns. Cut immediately
after each return. Each resulting **leg** has the unique form

    p followed by q,

where p is a nonempty simple path in A starting at some alpha_i, q is a
nonempty simple path in B ending at some beta_j, and the forward seam
end(p) -> start(q) actually exists. Thus

    (end(p), start(q)) is not (alpha_l,beta_l) for ANY l.       (F)

The leg starts at port i, ends at port j, and its next edge is the return
beta_j -> alpha_j. A complete cross cycle is a cyclic sequence of such legs.
Its legs must be pairwise vertex-disjoint, including all marked vertices
traversed incidentally inside a leg. A return port can occur at most once.

### Literal decoder theorem

Cross cycles of R_k correspond bijectively to cyclic leg words satisfying
(F), with pairwise disjoint vertex sets. Their size is the sum of the leg
sizes. Odd cross cycles are exactly the words with odd total size.

**Proof.** The only B-to-A edges are the listed returns. Between two
consecutive returns a cycle stays first in A, crosses A-to-B exactly once,
and then stays in B. This proves the unique decomposition. Conversely the
legs and return edges form a directed cycle; disjointness makes it simple.
Cyclic rotation of the leg list is precisely the choice of where to cut.
No singletons or directed two-cycles arise: the only two-vertex diagonal
leg would use the forbidden seam alpha_i -> beta_i. QED.

This is the actual new address dictionary. No scalar identity is being used
to postulate a cycle that the tournament does not contain.

## 4. The square-zero path matrix

Work in the finite commutative algebra

    R_V = Q[x_v : v in V] / (x_v^2 : v in V).

Write x_S for the product over S. Multiplication has the exact rule

    x_S x_T = x_(S union T) if S and T are disjoint, and 0 otherwise.

This is a commuting square-zero algebra, NOT an exterior/Grassmann algebra.
There are no permutation signs from reordering vertices.

Let P_i^A be the sum of x_V(p) over all simple A-paths starting at alpha_i,
and E_j^B the analogous beta_j-ending sum. Let P_il^A count paths starting
at alpha_i and ending at alpha_l, and E_lj^B count paths starting at beta_l
and ending at beta_j. A path with equal endpoints is the singleton path.
Use the same actual vertex variables in every entry. Define

    M_ij = P_i^A E_j^B - sum_(l=1)^k P_il^A E_lj^B.             (M)

Despite this subtraction, M_ij is coefficientwise nonnegative: its monomials
count exactly the legal legs. Each pair (p,q) has one forward seam, so the
forbidden classes in the subtracted sum are disjoint. In particular it is
not enough to remove only the diagonal singleton terms.

To keep literal addresses, attach an independent auxiliary variable
Y_(i,p,q,j) to every legal leg. Formula (M) then means the same leg list with
those labels attached. Distinct paths on the same vertex support remain
separate summands. The computations aggregate their coefficients after
independently checking the full literal-cycle decoder.

Let sigma replace every x_v by -x_v, leaving auxiliary variables unchanged,
and define Odd(f)=(f-sigma(f))/2. Then Odd selects total vertex degree,
NOT the number of returns.

### Multi-return cycle formula

    C_cross = Odd( sum_(ell=1)^k trace(M^ell)/ell ).             (C)

Its coefficient at a vertex support counts the literal odd cross cycles
with that support. With the leg labels retained, each literal cross cycle
has its own decodable monomial.

**Proof.** A term of trace(M^ell) is a cyclic leg word with a chosen starting
port. A repeated vertex makes its product zero. Each surviving word is a
literal simple cross cycle by the decoder theorem and uses ell distinct
return ports. It is counted exactly ell times, once per starting port.
Division by ell therefore leaves a nonnegative integer sum. Odd filters its
actual length. More than k legs would repeat a return port. Indeed the entire
matrix satisfies M^(k+1)=0 because any product of k+1 entries repeats a row
port alpha_i. QED.

One may instead enumerate each cycle starting at its smallest return-port
index. This gives the same formula without rational arithmetic and is the
literal enumeration used in the test script.

## 5. Positive packing formula and determinant expression

Let Z_X(t;x) be the multivariate packing polynomial of X: a packing of j odd
cycles on vertex set S contributes t^j x_S. With auxiliary cycle addresses,
retain their products too. Write Z_empty=1. Because intersecting cycles
multiply to zero,

    Z_X(t;x) = exp(t C_X(x)).                                  (E)

The exponential is a finite polynomial; its j! cancels the order of the j
chosen disjoint cycles. It is not a signed or asymptotic expansion.

All old cycles lie in A or B, and (C) is the complete new list. Therefore

    Z_Rk(t;x) = Z_A(t;x) Z_B(t;x) exp(t C_cross).                (P)

This replaces the exponentially indexed deletion list by one reusable
algebra-valued interface. It does NOT eliminate exponential storage: R_V
has dimension 2^|V|, and constructing all path entries can be expensive.
The advantage is exact multi-port composition and a literal address parser.

Since M is nilpotent, the formal trace-log determinant identity is finite:

    log det(I-M) = -sum_(ell=1)^k trace(M^ell)/ell.

It also follows directly from the determinant's permutation-cycle
expansion; the logarithm retains connected cyclic index words. Applying
sigma and subtracting yields

    exp(2 C_cross) = det(I-sigma(M)) / det(I-M).                (D)

The denominator has constant term 1, so it has a unique finite inverse in
R_V. At general t the formal expression is the t/2 power of this ratio,
interpreted through (P). At t=2r for an integer r>=1 it is its ordinary
r-th power.

### What the ratio does NOT permit

Setting every x_v=1 is NOT a homomorphism from R_V: it would turn x_v^2=0
into 1=0 and erase the disjointness condition. To get a scalar, sum the
coefficients linearly after all square-zero multiplication, or equivalently
extract a top coefficient after adjoining unused-vertex factors.

Let n=|A|+|B|. Using the inherited odd-cycle formula,

    Q_r(R_k) = [x_V] (product_v (1+r x_v))
        Z_A(2r;x) Z_B(2r;x)
        (det(I-sigma(M))/det(I-M))^r.                         (Q)

A selected cycle contributes 2r, an unused vertex contributes r. Thus a
packing F contributes 2^|F| r^(n-|V(F)|+|F|), exactly the ordered r-path-cover
formula with empty slots allowed. For r=1 this is H(R_k).

Labels are preserved until deliberately forgotten. A selected leg-label
monomial decodes to its index graph; its directed connected components are
exactly the chosen cross-cycle addresses. The old cycle labels decode the
old packing. This recovers the precise independent-support predicate and
its two signs per selected cycle, as in Session263.

## 6. The two-return formula: no determinants needed

For k=2 put

    a = Odd(M11),   b = Odd(M22),
    c = Odd(M12 M21).

Every term of a uses alpha_1,beta_1, every term of b uses alpha_2,beta_2,
and every term of c uses all four marks. Consequently a^2=b^2=c^2=ac=bc=0.
The exact positive formula is

    Z_R2(t;x) = Z_A(t;x) Z_B(t;x)
                 [1 + t(a+b+c) + t^2 a b].                  (2R)

The first channel chooses one one-return cycle; c chooses a single cycle
using both returns; ab chooses two compatible one-return cycles. Square-zero
multiplication also removes old cycles that intersect the selected legs.

For parity-separated matrices M_ij=E_ij+O_ij,

    c = E12 O21 + O12 E21.

Replacing c by O12 O21 would count even cycles as odd. This is a genuine
parity sidecar, independently of the number of crossed seams.

The k=1 case reduces to Z_A Z_B(1+t Odd(M11)), exactly RJ2 after coefficient
expansion. The new two-return channel and compatible-pair channel are the
first additional physical ingredients.

## 7. Explicit six-vertex profile-neutral selector exchange

Take both shores to be triangles with outgoing masks (2,4,1), label the
first shore 0,1,2 and the second 3,4,5.

- Aligned ports (alpha,beta)=(0,0),(1,1) give

      A = (50,44,57,17,34,8).

- Swapped ports (0,1),(1,0) give

      B = (42,52,57,18,33,8).

These masks fully specify literal tournaments, not isomorphism classes.
Their complete cycle lists and compatible address pairs are in the receipt.
The aligned list has ten odd cycles, with exactly compatible pairs
{1,10} and {6,7}; the swapped list has twelve odd cycles, with exactly the
compatible pair {1,12}. In each list the first and last cycles are the two
original shore triangles. Hence their exact conflict graphs are

    Gamma(A) = K10 minus two disjoint edges,
    Gamma(B) = K12 minus one edge.

The cycle-size counts give the polynomials stated in §1. Their difference is

    Z_A(t,u)-Z_B(t,u) = t u^5 (t u-2).                         (T)

The full cover substitution is Q_r(X)=r^6 Z_X(2r,1/r). It annihilates the
factor tu-2, proving equality for EVERY r, not just H,Q2,Q3. In particular
(H,Q2,Q3)=(29,320,1845) for both.

### Exact signed-state exchange, including the missing color

There is an explicit bijection on colored signed packing states for every
r. Match six triangle addresses and four five-cycle addresses, retaining
the two original triangles. Match the original compatible doubleton.
The aligned module has one remaining compatible pair of triangles, with
two sign bits and two cycle-color choices. The swapped module has two
remaining five-cycle singleton addresses, with one sign bit, one cycle
color and one unused-vertex color. Send

    first sign bit -> choice of the extra five-cycle address,
    second sign bit -> sign of that cycle,
    first triangle color -> its cycle color,
    second triangle color -> its unique unused vertex color.

Use any fixed bijections on the unused vertex sets for the other matching
states. This preserves all r-colored signed-state counts. It does not
preserve every original cycle address or individual shell vertex; its
precise loss is the exchange of one compatible doubleton for two singleton
addresses. It is not asserted to be a canonical bijection of Hamiltonian
paths, which would require choosing the inherited OCF correspondence.

Ordinally place either module after an unchanged protected core. The core,
its mark, and its zero-shell state survive literally, and the complete
unmarked Q profile is unchanged by this exchange. The selector geometry
and marked continuation data can change. This is a constructive trade,
not a claim that full Q determines an interface.

## 8. The rooted dipole and physical readout

Mark vertex 0 in both A and B. For a marked tournament (T,z) define

    h = H(T-z),  Q = Q_2(T-z),  k = K_(T-z)(-1),
    U = sum_(p starts z) H(T-V(p)),
    V = sum_(p ends z)   H(T-V(p)),
    d = sum_(p starts z) (-1)^(|p|+1) H(T-V(p)),
    e = sum_(p ends z)   (-1)^(|p|+1) H(T-V(p)).

Here K_X(s)=sum_(S subset V(X)) H(X[S])H(X-S)s^|S|, with H(empty)=1.
Direct enumeration of the explicitly listed tiny graphs gives

    shared (h,Q,k,H,d,e) = (9,84,0,29,4,5),
    aligned (U,V) = (52,41),
    swapped (U,V) = (50,43).

Moreover their full root-deleted profiles agree:

    Q_r(A-0)=Q_r(B-0)=r^5+6r^3+2r.

Thus the trade changes only the two directed ordinary rooted responses in
this eight-coordinate carrier, by (+2,-2). It preserves their sum as well
as the entire original and deleted unmarked cover polynomials.

A literal one-vertex readout distinguishes the source. Take T -> b and
reverse b -> 0. RJ2 gives

    H(readout(T)) = H(T) + 2 sum_(even-size p starts 0) H(T-V(p))
                 = H(T)+U-d.

Therefore the two actual seven-vertex counts are 77 and 75. No abstract
linear functional or unattainable test vector is involved.

For reference, deleting vertices 1 or 3 gives r^5+4r^3 in both; deleting
vertex 4 gives r^5+6r^3+2r in both. At 2 or 5 the aligned profile is
r^5+8r^3+2r and the swapped profile is r^5+8r^3+4r. Equality is asserted at
the chosen root, not at every vertex.

## 9. Hostile controls and precise finite scope

Reproduce from repository root:

    python3 05-knowledge/results/session264_constructive_20261111/returns/multi_return_probe.py

The script constructs actual outgoing masks, independently enumerates
literal simple odd cycles, and compares every cross cycle to its leg-word
parse. It then checks all packing coefficients by occupied vertex set and
number of cycles, and separately counts Hamiltonian paths by subset DP.

The complete small universe consists of every labelled pair of shores of
orders 1 through 3, prefix-selected marks, and every permutation pairing
of those marks: 121 one-return, 200 two-return and 384 three-return cases.
There are also 24 seeded constructors of orders 7 and 8. For those and the
explicit examples it checks the determinant ratio and r=1,2,3 cover counts
against independent subset convolution. This is not an exhaustive census
of all marked tournaments or all choices of marks at larger orders.

One seven-vertex three-return control has 18 one-return, 14 two-return and
4 three-return odd cycles. Hence the multi-return terms are exercised;
testing only two triangles would leave the three-return odd channel empty.
Its rows are (114,52,105,83,36,8,50), with H=109.

A concrete seam hostile takes two transitive triples, with ports
(0,2),(2,0). Multiplying uncorrected endpoint-path lists predicts ten cross
odd cycles; only eight exist. The spurious supports are masks 47 and 61.
The first false implication is treating a forward cross edge as present
after it has itself been reversed. Condition (F) and formula (M) remove it.
This probe motivated the final exact seam correction; it is not a rejection
of the inherited one-return formula.

The receipt records source-code SHA256, literal examples, all cycle lists
for the rooted trade, matrices, and finite counts. The proof is the cycle
bijection; computation supplies concrete initial controllers and hostile
checks, not an inference from a finite sample to all k.

## 10. Scope and next physical target

**Preserved:** literal shores; every permitted cycle with its return-word
address; disjointness and parity; exact signed packing states; all Q_r;
protected external core under ordinal attachment.

**Not preserved by scalarization:** vertex capacity, cycle labels, rooted
responses, or future arbitrary seam behavior. The square-zero carrier is
exact but can be exponentially large. The reverse-matching restriction is
explicit. A general nonmatching return graph would need a different port
parse and is not covered silently.

**New constructive targets:** build a library of profile-neutral rooted
trades with several independent response directions; compose them with the
root-sharing carrier; and find a finite physical controller bank whose
positive response cone permits prescribed inverse steps. The present
(29,29) trade supplies one genuine direction and an immediate readout, not
that whole bank. Arbitrary odd targets and Chan–Pak remain open.
