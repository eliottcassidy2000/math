# Uniform ternary schedules: explicit square-sum prefixes and cycle multiplicity

**PROVED from complete independently audited certificate catalogs.** This
constructs a specified family of prefixes and rays. It does not solve arbitrary
prescribed-path embedding, classify all square-sum rays, or claim novelty for
ordinary interval Hamiltonicity or general cycle-count asymptotics.

## Three typed prefix operators

| Symbol | Radix Q | Offset c | Fixed prefix H | Added vertices |
|---|---:|---:|---|---:|
| A |25|12|(1,3)|2|
| B |25|8|(1,3,6,19)|4|
| C |49|4|(1,3,6,19)|4|

All positions are one-based. Define

    F_a(p1,...,pk)=H_a || (Q_a*p_j+(-1)^j*c_a)_(j=1..k).

For a finite outer-to-inner word w=a1...al over{A,B,C}, put

    P_w=F_a1∘...∘F_al(1,3),
    R_w=∏ Q_ai,
    S_w=(85*R_w−1)/2.

For the empty word use P_empty=(1,3), R_empty=1 and S_empty=42.

### Uniform all-size extension theorem

For every nonempty word w and **every** integer N>=S_w, the square-sum graph
on [1,N] has a Hamiltonian cycle whose oriented traversal begins with the
original unscaled P_w. If w has a,b,c occurrences of A,B,C, respectively, then

    |P_w|=2+2a+4b+4c,
    S_w=(85*25^(a+b)*49^c−1)/2.

These are individual existential embeddings, not simultaneous embeddings of
all coded prefixes in one cycle. In particular, pure B has4l+2 vertices and
threshold(85*25^l−1)/2; pure C has4l+2 vertices and threshold(85*49^l−1)/2.
They are **different rays**, so this is not a lowered threshold for the old
C-ray itself. Pure A has2l+2 vertices at the same25-radix threshold.

## Proof: a common input type, not a common historical construction

A weak B pair has Hamiltonian paths on[1,n] and[1,n+1], both starting1,3,
ending at8 for an even-length member and10 for an odd-length member, with
position-parity agreement for every common label. All three paired catalogs
consume and produce exactly this type. A and C impose stronger sufficient
within-block shifts; B allows independent split-head colors. B's head colors
and tail shifts each agree between outputs, which is exactly the compatibility
needed by the next operator. No stronger provenance condition is required.

Use the certified weak bases42..2081. For B(n,w), n>=42, one may use any fixed
audited all-size weak-B constructor when w is empty (for example A alone).
If w=a v, set Q=Q_a and h=(Q−1)/2. For n<42Q+h, return the weak base B_n;
this lies in42..1061 for A/B and42..2081 for C. Otherwise put

    m=floor((n−h)/Q), d=n−Qm in[h,3h],

recursively construct B(m,v), then apply the paired certificate for a at
state(m mod2,d). This terminates and returns a weak B pair.

Since S_(av)=Q*S_v+h, induction on word length gives the common prefix P_w
whenever n>=S_w. Apply the terminal consumer for a to B(m,v). It has the
same fixed-prefix/intact-block start and closes to1: at8 for A/B, at8 or15
for C. This proves the cycle theorem.

The threshold formula telescopes because composing x↦Qx+(Q−1)/2 gives
x↦Rx+(R−1)/2; start from x=42. The base cutoffs1061 and2081 are essential:
using arbitrary larger weak bases instead of recursive outputs would lose
the required deeper prefix.

## Uniform vertex formulas and every square equation

For an infinite schedule a1,a2,... define

    R_t=∏_(i=1..t) Q_ai,
    C_t=Σ_(i=1..t) c_ai*∏_(k<i) Q_ak.

At depth t, the block contributed by operator a_(t+1) has vertices

    (R_t*h_j+(-1)^j*C_t)_(j=1..|H_(a_(t+1))|).

For a finite word of length l, append the transformed final pair
(R_l−C_l,3R_l+C_l). All block lengths are even, so these parity signs agree
with the global ordering. Within an A block the root is2√R_t; within a B/C
block the roots are2√R_t,3√R_t,5√R_t. The outgoing root is respectively
4√R_t,6√R_t,8√R_t for A,B,C, because

    3+25−12=16, 19+25−8=36, 19+49−4=64.

Thus every adjacent sum is explicitly a square.

## Simplicity and exact decoding of the coded family

Each F_a preserves simplicity: a collision of lifted labels would force
Q_a(x−y) to be0 or±2c_a, while0<2c_a<Q_a. Its image is disjoint from H_a,
as follows from the residues±12 modulo25 for A, ±8 modulo25 for B, and±4
modulo49 for C. None of the corresponding fixed-prefix labels has one of
those residues. Positivity is immediate from Q_a−c_a>0.

Every F_a(1,3) starts1,3, so prefix preservation shows P_w is a prefix of
P_(wa) for each next symbol a. Operators can be decoded without ambiguity:
A has third vertex13; B/C have third vertex6; then B has fifth vertex17 and
C has fifth vertex45. Remove the decoded H_a and invert its affine lift to
continue decoding. Therefore different words of the same operator length,
or more generally any antichain in the word-prefix order, give vertex-prefixes
none of which is a prefix of another. Different finite words are distinct.

Each infinite schedule yields a simple infinite square-sum ray, and decoding
makes distinct infinite schedules give distinct rays. This family has
cardinality continuum. Every finite prefix of every one of its rays extends
to a Hamiltonian cycle on every sufficiently large interval. These sparse rays
are **not claimed Hamiltonian on the positive integers**.

## Distinct finite cycles, not simultaneous embedding

Every constructed cycle contains edge1–3. An undirected Hamiltonian cycle
containing this edge has a unique traversal beginning1,3; reverse orientation
cannot also begin1,3. Thus incompatible prescribed prefixes require distinct
undirected cycles.

Consequently, for every l>=1 and N>=(85*25^l−1)/2 there are at least2^l
distinct Hamiltonian cycles, by the words{A,B}^l. For N>=(85*49^l−1)/2,
all3^l words of length l yield distinct cycles. More sharply, words of length l
with exactly c occurrences of C give

    binomial(l,c)*2^(l−c)

distinct cycles for all N>=(85*25^(l−c)*49^c−1)/2. These are explicit
certificate-derived lower bounds, not claims of optimality or literature priority.

### A weighted prefix-code bound

Let α be the unique positive root of

    2*25^(−α)+49^(−α)=1.

For N>=2082 put X=(2N+1)/85. Starting with the empty word of product1,
expand a node of product R into its three children exactly when49R<=X.
Every leaf is nonempty and has X/49<R<=X. The frontier is finite and is an
antichain. Each expansion preserves the weight sum ΣR^(−α), so that sum is1
on the frontier. Hence the number L of leaves satisfies

    L >= (X/49)^α.

Each leaf gives a different Hamiltonian cycle by the theorem and decoding.
This is a constructive lower bound with an exact frontier algorithm; the
irrational exponent is only a compact summary of its weighted count.
The proof does not infer any major open conjecture from the cycle multiplicity.

## Certificate and replay boundary

A:50 paired+50 terminal states, independently audited with51,027,430 checks.
B:50 paired+50 terminal states, independently audited with51,028,631 checks,
including separate-head parity checks and an explicit hostile head-color test.
C: the already frozen Q49 core plus threshold42 base extension. No old A/final
transition tables are needed. The old compressed base bundle and eleven small
bases42..52 are the only input data, beyond these finite transition catalogs.

The portable constructor and independent auditors replay actual cycles and
all universal states. Mixed controls cover every ordered pair of operators
at its exact threshold and threshold+1, all six ABC permutations, and AAA,
BBB,CCC at depth3. Finite controls supplement the universal composition proof;
they are not its justification. The piece-overlap compatibility theorem explains
precisely why the relaxed B output remains a valid input to both other operators.

## Exact optimal frontier inside this certificate-coded family

The previous full-ternary frontier gives a clean exponent proof. One can count
a larger, exact maximum antichain in the finite schedule tree
{w:R_w<=X}. Expand a node whenever its A/B children fit; add its C child also
when it fits. Every internal node has two or three children. Its leaves are
an antichain, and any other antichain injects into these leaves by assigning
a descendant leaf to each node (the assigned subtrees are disjoint).
Thus this leaf set is a maximum-cardinality antichain **within this coded
family and budget**, not an optimum count of all square-sum cycles.

For1<=X<25 the count f(X)=1; otherwise

    f(X)=2*f(X/25)                         if25<=X<49,
    f(X)=2*f(X/25)+f(X/49)                 ifX>=49.

For N>=1062, the root has children and every selected leaf is nonempty, so
f((2N+1)/85) is a valid distinct-cycle lower bound. This improves or equals
the complete-ternary frontier bound above. Exact examples from the recurrence:

|N|Certified distinct cycles at least|
|---:|---:|
|1062|2|
|2082|3|
|26562|5|
|102042|9|
|16601562|35|
|1000000000|145|

`schedule_frontier.py` computes this count using only integer comparisons and
memoization on products R=25^a*49^b. There are O((log N)^2) such states; this
is a bound on arithmetic-state count, not an assertion that all cycles can be
materialized in polylogarithmic time. Listing the entire frontier or producing
its cycles incurs the corresponding output cost. The count does not require
enumerating Hamiltonian cycles or searching the graph.

### Exact growth exponent of the coded antichain count

For any finite word antichain, the full ternary-tree weights R_w^(−α) have
sum at most1: extend all its nodes to a common operator depth and use the
weight-conservation identity. Since every admitted product is at mostX,
each node's weight is at leastX^(−α). Thus f(X)<=X^α. Combining this with
the complete-frontier lower bound gives, for X>=49,

    (X/49)^α <= f(X) <= X^α.

Therefore the optimal antichain count **inside this coded family** is
Θ(X^α), with α the unique root2*25^(−α)+49^(−α)=1
(approximately0.3205317162). This upper bound is not an upper bound on the
number of Hamiltonian cycles of the square-sum graph; it only bounds this
particular prefix-certificate counting method.
