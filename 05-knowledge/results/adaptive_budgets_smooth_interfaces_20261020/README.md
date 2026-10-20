# Adaptive budgets, a sharp prime-batch scale, and a three-wing boundary

**Status: PROVED scoped propositions with same-agent proof review and
FINITE-EXACT controls; no independent-agent audit or priority claim.**
The smooth-number input is CITED. Universal Collatz entry and an improved
Hadwiger–Nelson bound remain OPEN. These are research propositions, not
new numbered canon entries.

## Inheritance and portfolio

Anchor: source-preserving Collatz budgets and the exact A000568 response
interface. Niche: a six-colour boundary relation beyond triangle strips.
Wildcard: commutative operation fibres, Farey/braid clocks and half-planes.

Closest mechanisms are Propositions 1, 7 and 8 of
[prime batches and slack six](../prime_batch_triangle_20261020/README.md):
the smooth-partition upper bound, exact source cylinders, and monotone
contiguous-block envelopes. Hostiles are source 63's five-block first exit,
source 27's failure of every one-block prefix, and the full four-colour
relation of a two-link triangle strip. The corrected near miss is the
shifted affine anchor: `T_1^s(1)=2(3/2)^s-1`, not the old expression with
coefficient 3. The least-used useful sidecar is the divisor-incidence
coordinate of [THM-516, Burnside core kernel](../../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md),
from June: it upgrades a bank of prime tests to a complete smooth bank.

The live board has six objects: source cylinders; ordered block debt;
smooth divisor profiles; full colour-extension tensors; swap-fixed operation
fibres; and integer lifts of quotient histories. See [CONNECTIONS.md](CONNECTIONS.md)
for incoming work, map/loss contracts and next tests.

## 1. The prime-batch width has a matching leading logarithmic scale

Let W(N) denote peak occupied `(degree, future gcd response)` states in the
ordinary positive-weight tournament recurrence through degree N, processing
odd cycle lengths by decreasing largest prime factor, with length 1 last.
Within a prime batch the ordering is arbitrary. The response at a cut is

    R_j = sum_(a in past cycle multiset) gcd(a,j)

for each unprocessed odd j affordable at the committed degree. This is the
same W and the same recurrence as the preceding packet, not a count of
isomorphism classes or a general algorithm complexity invariant.

**Proposition 1 (PROVED, using the cited smooth-number estimate).** Put
`B_N=sqrt(log N log log N)`, with natural logarithms. Then

    log log W(N) = (1/2) log N - (1/2+o(1)) B_N.          (1)

Equivalently,

    W(N) = exp(sqrt(N) exp(-(1/2+o(1)) B_N)).             (2)

The same leading double-logarithmic scale holds for the maximum rank of
the wordwise weighted past/future interface, allowing a cut and fixed past
and future degrees whose sum is at most N. Thus changing only a linear
encoding of this particular interface cannot remove its width scale.
No lower bound for arbitrary methods of evaluating A000568 is asserted.
No lower bound after signed core cancellation is asserted.

### A complete smooth probe bank

Take a real y>=3, an odd prime q>y, and an integer z>=1. Let

    S = {c<=z: c odd and every prime factor of c is <=y},
    D = |S|,       M = sum_(c in S) c.

Assume `(q+1)M<=N`. At the cut just after the q-batch, choose a history
indexed by any subset A of S, using cycle lengths `{qc:c in A}`. Every
such length has largest prime factor q; every c in S is still a future
length. All these histories and even the entire future bank are affordable.
For j in S their responses are

    R_j(A) = sum_(c in A) gcd(qc,j) = sum_(c in A) gcd(c,j). (3)

Let G be the matrix `G_(c,j)=gcd(c,j)` on S. The set S is divisor-closed.
Writing `Z_(c,d)=1_(d|c)` on S gives

    G = Z diag(phi(d):d in S) Z^T.                       (4)

This uses `sum_(d|g) phi(d)=g`. In increasing numerical order Z is lower
triangular with diagonal 1. Thus G is positive definite and invertible,
and (3) distinguishes all `2^D` histories. Ordinary recurrence weights
are positive, so these occupied states cannot disappear by cancellation.
Different degrees only help distinguish states; every response entry in
the bank remains present even at the largest committed degree qM.

The old prime-only bank is a strict restriction of this construction.
Composite probes carry independent divisor coordinates, including prime
powers. Replacing c by its radical would invalidate (4) as a recovery map.

### A rank lower bound at common degrees

Pigeonhole subsets by `sum A`, which takes at most M+1 values. There is a
family F of at least `2^D/(M+1)` subsets with a common sum s. Use past
lengths qA and future lengths B for A,B in F. Past degree is qs, future
degree s, and their sum is at most N. After positive row and column
factors from internal edge-orbits and cycle-index denominators are removed,
the weighted continuation matrix is

    K(A,B) = 2^(sum_(c in A,d in B) gcd(c,d)).             (5)

This matrix is strictly positive definite on any distinct subset vectors.
Indeed set `x_A=sqrt(log 2) G^(1/2) 1_A`; then K(A,B)=exp(x_A dot x_B).
For any real coefficients t_A its quadratic form is the nonnegative sum

    sum_(multi-index alpha) (sum_A t_A x_A^alpha)^2/alpha!.

If it vanishes, every polynomial moment vanishes. For each one of the
finitely many distinct x_A there is a polynomial equal to 1 there and
zero at all the others: multiply the normalized linear forms separating
that point from each other point. Hence every t_A is zero. This proves
strict positivity without a numerical eigenvalue argument.

Therefore the interface rank is at least `2^D/(M+1)`. Its upper bound by
the occupied response-state count follows because all future transition
weights factor through that state, up to a past-only scalar. These are
wordwise interfaces; summing all future words into one scalar is a different
operation and has no such rank lower bound.

### Choosing the bank: the previously missing lower bound

The only new external input is the classical estimate

    log Psi(z,y) = log z - (1+o(1)) u log u,
    u=log z/log y,

as y,u tend to infinity with `u<=y^(1-epsilon)` for a fixed epsilon>0.
This is Corollary 1.3 of Hildebrand–Tenenbaum,
[Integers without large prime factors](https://tenenb.perso.math.cnrs.fr/PPP/Survey_friables.pdf),
printed page 417 (PDF page 7); it is unconditional in this range.
The result's own uniformity condition is checked below. The more restrictive
Psi/rho asymptotic in their Theorem 1.1 would also suffice here.

Discarding even integers does not change the required logarithmic estimate:
the odd-part map has fibres of size at most `1+floor(log_2 z)`, so

    Psi(z,y)/(1+floor(log_2 z)) <= D <= Psi(z,y).           (6)

For any fixed eta>0, put v=log N and B=sqrt(v log v), and choose

    y=exp(B/2),   y<q<2y prime,
    z=floor(exp(v/2-eta B)).                             (7)

Bertrand's prime bound suffices to choose q; a least-prime asymptotic is
unnecessary. Here `u=(1+o(1))sqrt(v/log v)`, which is smaller than
`y^(1/2)` for all sufficiently large N, so the cited estimate applies.
It gives `u log u=(1/2+o(1))B` and, by (6),

    log D = v/2-(eta+1/2+o(1))B.                         (8)

Since M<=zD,

    log((q+1)M) <= v-(2eta+o(1))B < v.                  (9)

This verifies the full degree budget, rather than merely the size of an
abstract matrix. Equations (3)--(9) give

    log log W(N) >= v/2-(eta+1/2+o(1))B.

The rank loses only `log(M+1)=O(v)` from its logarithm; this is negligible
compared with D, so it has the same lower bound. Since eta>0 is arbitrary,
combine with Proposition 1 of the previous packet, which proves the upper
bound `log log W<=v/2-(1/2-o(1))B` uniformly inside every batch. This proves
(1) and the rank assertion. QED.

This closes the leading-scale gap for the **specified** prime-batch
interface. Sharper second-order asymptotics and genuinely different
interfaces remain meaningful targets; another prime-only lower bound
would now be weaker than what is already available.

## 2. Collatz block budgets: source arithmetic and clock requirements

Use the odd map `T_a(n)=(3n+1)/2^a`, with a the exact valuation, and keep
the original positive odd source n>1. For a word w of length m and clock
K=sum w, the inherited cylinder and endpoint are

    n = (2^K-B(w)) (3^m)^(-1) mod 2^(K+1),
    T_w(n) = (3^m n+B(w))/2^K.

For a contiguous partition P, U_P denotes composition of decreasing
rearrangements **only in the upper-bound calculation**. The physical word
is never sorted. The inherited theorem proves monotonicity under splitting
P and preservation of n,m,K. The unaccelerated clock convention is m+K,
counting each odd `3n+1` operation and every subsequent halving separately.

### Integer dynamic programming

**Proposition 2 (PROVED, implementation refinement of the inherited DP).**
Let `K_j=sum_(i<=j) a_i`, and let B_(i,j) be the numerator of the decreasing
rearrangement of positions i+1,...,j. If d_(b,j) is the minimum numerator
at the common denominator `2^K_j` over b contiguous blocks, then

    d_(0,0)=n,
    d_(b,j)=min_(b-1<=i<j)
             [3^(j-i) d_(b-1,i) + 2^K_i B_(i,j)].        (10)

Only defined predecessor states occur. A candidate pays exactly when
`d_(b,m)<n 2^K_m`. All operations and comparisons are integers.

Proof: each segment has slope `3^(j-i)/2^(K_j-K_i)`; substitute the
preceding numerator at denominator `2^K_i`. All slopes are positive, so
retaining the least previous value is valid. This is precisely the earlier
rational DP, with denominators aligned rather than repeatedly reduced.
For repeated exponent a with multiplicity r its segment numerator is
`(3^r-2^(ar))/(3-2^a)`, which is an integer geometric sum. QED.

**FINITE-EXACT extension.** Every odd source 3<=n<=200001 was followed
to its first strict descent, with a declared cap of 1000 odd steps. No cap
was hit. Exact DP and direct rational replay of the selected cuts give:

| Minimum blocks on that first-exit prefix | Sources |
|---|---:|
|1|99,562|
|2|432|
|3|4|
|4|1|
|5|1|

The four- and five-block sources remain 27 and 63. This does not assert a
universal five-block first-exit bound, an asymptotic frequency, or universal
existence of first exits. The maximum observed first-exit length is 85,
at 35655. The receipt retains record words, cuts, exact upper endpoints
and clocks, not just the histogram.

### Two blocks suffice after certified terminal arrival, with a paid clock

**Proposition 3 (PROVED; a certificate-compression theorem).** Suppose an
authenticated finite word w takes n>1 to 1. Put `U=U_{one block}(n)`.
After appending r actual exponent-two steps at 1, the two-block envelope is

    V_r = 1 + (3/4)^r (U-1).                             (11)

Thus it pays strict descent once `(3/4)^r(U-1)<n-1`. If it already pays,
r=0 is allowed and the empty second block is omitted. Otherwise the least
r is found by exact integer comparison, with no floating logarithm needed.
One may always take `r<=2m+4`. The retained integer clock is K+2r and
the unaccelerated clock is m+K+3r.

Proof: `T_2(x)=1+(3/4)(x-1)`, and the actual orbit at 1 has exponent 2
on every repetition. Every T_a is increasing and bounded above by T_1 on
positive arguments; hence `U<=(3/2)^m(n+1)-1`. For n>=3,
`(n+1)/(n-1)<=2`, while `(3/4)^(2m)(3/2)^m=(27/32)^m<=1` and
`2(3/4)^4<1`. This proves the explicit bound. The authenticated appended
word retains the original cylinder, now at modulus `2^(K+2r+1)`; it does
not reset the source to 1 in the certificate. QED.

For 27, its 41-step route to 1 followed by 25 exponent-two repetitions
gives two blocks. For 63 the corresponding numbers are 39 and 22.
These are the least r for the particular split at first arrival at 1,
not the globally shortest two-block certificates. Source 27's inherited
all-future one-block obstruction shows two is sharp for unrestricted
certificate length among convergent sources.

**Exact reformulation, not a proof of convergence.** Universal positive
Collatz convergence is equivalent to every odd n>1 admitting some paid
at-most-two-block prefix. Forward implication is Proposition 3; reverse
implication follows by strong induction from actual descent below n.
Consequently a small block budget without a clock/arrival requirement
does not resolve the open problem. Productive strengthening must prove
entry or bound the authenticated clock independently of already knowing
terminal arrival. First-exit budgets remain a different, stronger target.

## 3. A geometric three-input six-colour relation

**Proposition 4 (PROVED, explicit gadget; no novelty claimed for the wheel).**
Let `r_j=(cos(j*pi/3),sin(j*pi/3))`, j mod 6, let h=0, and put
`t_i=r_(2i)+r_(2i+1)` for i=0,1,2. The ten points consist of a regular
unit hexagon, its hub, and three outward equilateral tips. Their induced
unit graph has exactly:

* six spokes h--r_j;
* the six rim edges r_j--r_(j+1);
* six tip edges t_i--r_(2i), t_i--r_(2i+1).

The nine ports are the three ordered wings
`(r_0,r_1,t_0), (r_2,r_3,t_1), (r_4,r_5,t_2)`.
For every proper assignment of six colours to the induced port graph,
the exact number of extensions is

    6 - |{colour(r_0),...,colour(r_5)}|.                  (12)

Thus the relation excludes exactly rainbow rims. Every projection onto
two complete wings is the full proper relation of their induced graph.
The full three-wing relation is not recoverable from those pairwise
projections or from the three unordered wing palettes.

Proof: all coordinates lie in the triangular lattice. Checking pairwise
squared distances gives the listed 18 edges and no coincidences. The sole
uncoloured vertex h can use exactly the colours absent from its six
neighbours, giving (12).

For the projection, retain the first two wings. Their only cross edge
is r_1--r_2. The four assigned rim vertices use at most four colours.
Choose a colour e outside that set and require the hub to use e. There
are five other colours; choose the new r_4 to avoid e and r_3, then r_5
to avoid e, r_4 and r_0. Choose t_2 distinct from r_4 and r_5. These
choices always exist, proving extension. Rotation covers all three pairs.

For palette loss, the following triples have identical unordered palettes:

    ((0,1,2),(2,3,0),(0,5,4))   extends;
    ((0,1,2),(2,3,0),(4,5,0))   does not extend.

Both obey every port edge. The first repeats rim colour 0, leaving colour
4 for the hub; the second puts all six colours on the rim. QED.

### Complete boundary control

Fix the first wing to colours (0,1,2), which is valid by global palette
symmetry. Enumerating the other two ordered proper triangles gives 8,336
proper port assignments. Exactly 384 are blocked, and 7,952 extend. The
extension-count histogram at 0,1,2,3,4 is respectively
`384,3456,3840,640,16`.

Independent counting: the rim has `5^6+5=15630` proper six-colourings.
Each tip has four choices, so there are `15630*4^3/120=8336` normalized
port assignments. Rainbow rims give `6!*4^3/120=384` blocked assignments.
The two-wing projection has all 100 ordered proper choices of the second
triangle compatible with its one cross edge to the fixed first wing.

### Why this is a boundary gadget, not six-chromatic forcing

The entire gadget is three-colourable. In coordinates `a+b exp(i*pi/3)`
the rule `a+2b mod 3` properly colours every unit edge of its triangular
lattice. In particular this construction does not contradict
[THM-4614](../../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md):
both its field and its fixed eight-port configuration differ here; our hub
has six port neighbours, not at most two. None of that theorem's 570 patterns
has been eliminated by this new-port calculation.

Moreover, gluing these lattice gadgets along a shared unit edge, by plane
isometries and with all cross distances retained, keeps them in one common
triangular lattice. Indeed a unit lattice segment generates the lattice
with its 60-degree rotation; any isometry identifying two such segments
identifies their lattice cosets. Such edge-gluing therefore cannot create
a graph requiring four or more colours, even if new cross edges appear.
This is an all-depth obstruction for that construction class.

The useful reroute is to retain the nine named ports and the full ternary
relation (12), then seek a compatible context outside this common lattice
that forces a forbidden rim assignment. Pairwise triangle transitions
cannot certify that obligation. A one-port rotation may leave the lattice,
but only authenticated cross contacts can make its new relative angle
relevant; arbitrary abstract gluing supplies no coordinates.

## 4. Half-planes and recombination

**Proposition 5 (standard consequence, proved here for the requested lens).**
The unrestricted unit-distance chromatic number of every nonempty open or
closed Euclidean half-plane equals that of the whole plane.

Every finite point configuration can be translated into the open half-plane,
preserving all distances. Hence a k-colouring of that half-plane gives a
k-colouring of every finite plane unit-distance graph. The
[de Bruijn–Erdos compactness theorem, Theorem 1](https://www.renyi.hu/~p_erdos/1951-01.pdf)
then gives a plane k-colouring. The reverse direction is restriction.
This is for unrestricted colourings; no measurable or periodic assertion
is imported. A bounded-width strip is a different object. QED.

The translation moves the whole marked configuration. It does not construct
a new context at prescribed absolute port positions, so this global equality
does not discharge a fixed-port forcing problem.

The Banach–Tarski prompt suggests testing piecewise congruent
recombination, but such a map preserves distances only within each piece.
Cross-piece unit edges are not invariant: two points at unit distance,
placed in separate pieces, lose their edge after translating only one
piece. Thus any colouring-transfer claim needs a complete cross-piece
edge ledger. The relevant resource is compatible boundary information,
not a paradoxical cardinality argument. No equidecomposition theorem is
used in Propositions 1--5.

## Reproduction and remaining targets

Run from repository root:

```text
python3 05-knowledge/results/adaptive_budgets_smooth_interfaces_20261020/verify.py --limit 200001
python3 -O 05-knowledge/results/adaptive_budgets_smooth_interfaces_20261020/verify.py --limit 200001
python3 05-knowledge/results/adaptive_budgets_smooth_interfaces_20261020/separate_controls.py
```

The producer uses only standard-library exact arithmetic. The separate
checker imports no producer code; see [AUDIT.md](AUDIT.md) for its universe,
proof review, limits and receipts. The finite checks supplement, rather
than establish, the all-order arguments.

Next obligations are: an arrival/clock estimate for Collatz not paid from
known convergence; second-order prime-batch asymptotics or a different
interface; and coordinates for a non-lattice context forcing a forbidden
three-wing assignment while preserving every marked port.
