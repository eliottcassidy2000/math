# Four binary blocks and a caterpillar route to balanced coverage

Session 199. **PROVED:** a binary-block compiler and a sparse physical graph
interpretation of bounded root changes; the entire unit-deletion source class omits a
positive fraction of odd targets. **FINITE-EXACT:** at most three root changes
produce every odd balanced target from 27,783 through 1,000,001. This is a
strictly bounded-root-change result, not an unrestricted source search.
**OPEN:** the explicit four-block coverage lemma in Section 5. Neither its
large finite test nor the bit interpretation proves the all-order lemma.

## 1. Inheritance and exact domain

A source is a nonempty rooted tournament with exactly one directed Hamilton
path starting at its root. Its response is (D,E,U,d), with the path, endpoint,
reserved-cover and deletion conventions of the frozen positive-source proof.
Every nontrivial source is generated from t=(1,0,2,1) by

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

Chronological words act from left to right. Only J changes the marked root.
The number of J letters is a complexity of a supplied word, not an invariant
of an unmarked tournament. The least such number among words for a response
is well-defined, but need not equal the minimum for any specified graph.

The inherited canonical balanced construction is

    s . O J,       H=D+2U=:G(s).

It has one Hamilton cycle and a balanced root, and every balanced source is
obtained this way from a physical source, with the singleton giving H=3.
The companion finite repair table handles the 42 full balanced gaps through
50,001 other than 7 and 21. No new nonexistence proof for 7 or 21 is claimed.

The new restriction is to source words with at most three J letters. It keeps
arbitrarily long O/Z blocks: it is not a finite word menu. A canonical hostile
is zero root changes, which already has an all-order density obstruction below.
The missing sidecar in a scalar approximation is the ordered binary block,
not merely its length, digit count, modular image, or slope.

## 2. Exact binary-block compiler

Let epsilon=(epsilon_0,...,epsilon_(k-1)) be a binary word. Its chronological
source word uses O for zero and Z for one. For k>=1 define

    B=2^k,
    A=sum_i epsilon_i 2^i,
    R=sum_i epsilon_i 2^(k-1-i),
    C=sum_(0<=j<i<k) epsilon_j epsilon_i 2^(i-j-1)
        +sum_(1<=i<k) epsilon_(i-1) epsilon_i.

**Theorem 1.** On every physical input source, and indeed as an identity on
four formal coordinates, the block has transfer

    D'=D-epsilon_0 E+A U-Cd,
    E'=epsilon_(k-1)d,
    U'=BU-Rd,
    d'=d.                                                   (1)

The empty block is the identity, not the k=0 specialization of E' above.
Every such transfer has a literal O/Z graph recipe.

**Proof.** Immediately before position i,

    U_i=2^i U-sum_(j<i)epsilon_j 2^(i-j-1)d.

The old E is used only at i=0. At later positions E_i=epsilon_(i-1)d.
Summing the Z increments epsilon_i(U_i-E_i) proves the D identity; the
other three follow directly from the source laws. QED.

The coefficient C is even: for adjacent positions its two contributions
sum to 2, and every nonadjacent pair has coefficient 2^(i-j-1), also even.
Thus the unit-deletion seed gives the explicit integer expressions

    D=1+2A-C, U=2B-R,
    G=4B+1+2(A-R)-C.                              (2)

For a fixed length k, R runs through 0,...,B-1, so U runs bijectively through
B+1,...,2B. The empty block gives U=2. Consequently every integer U>=2 has
one O/Z word from t. Its D is the function g(U) determined by

    g(1)=1,
    g(2n)=g(n),
    g(2n-1)=g(n)+n-(n mod2), n>=2.              (3)

Equation (2) is a weighted binary-pair interaction formula for that exact
physical sequence. It preserves orientation order: the two terms A and R
are opposite bit readings, while C retains the ordered pair interactions.
For example OZ and ZO have the same length and digit count but produce
(5,1,7,1) and (3,0,6,1), hence G=19 and G=15.

## 3. The actual graph is a chronological caterpillar

Place the vertices in their birth order, starting with the transitive edge.
Relative to the forward transitive tournament, O reverses no edge. A Z or J
birth reverses exactly the edge from the newborn to the current root.
A Z newborn is never subsequently the current root; a J newborn is.

**Theorem 2.** The graph of reversed edges of any source word is a forest
with at most one nontrivial component. That component is a caterpillar:
the successive current roots form a path, and each Z newborn is a leaf on
its then-current root. All O newborns are isolated in this reversal graph.
A word with r J letters has a distinguished spine of r+1 vertices.

This follows directly from the birth edges. The word supplies an ordering
constraint beyond the abstract caterpillar: all leaves at a spine vertex
are born before the next spine vertex. Arbitrary orderings of a caterpillar
are not asserted to be sources. The statement is about reversal edges in
this specified ordering, not about an induced forest inside the tournament.

Appending OJ produces a balanced C1 graph and extends that spine by one.
Thus the tested r<=3 class has balanced witnesses with a reversal forest
of this ordered type and a spine of at most five vertices. This is a narrow
physical construction class despite unbounded total graph order.

## 4. A rigorous hostile: unit deletion has positive-density gaps

Let B_0={g(u)+2u:u>=2}, adjoining the singleton value 3 if desired. Every
zero-J source has precisely this readout, by (1)-(3).

In fact the zero-J response class is exactly the physical class d=1.
O/Z preserve d; J sets d to the preceding D, and D>=d, so d cannot decrease.
If a word ending with d=1 ever uses J, its last such birth must have a
predecessor with D=1. The source grammar makes that predecessor a transitive
source with response (1,0,2^k,1), k>=1. Its J output is
(2^k+1,1,2^k+1,1), which is exactly t . Z^k. Replace the prefix through J
by this zero-J word, then retain its O/Z suffix. The replacement is also a
rooted graph isomorphism: send the new J root to the initial root of the
star word, and the old transitive chain to the remaining vertices in order.
Thus no deletion-one graph or response is omitted by the zero-J image.
The singleton boundary only adds G=3.

**Theorem 3.** Put

    c=sum_(k>=0) 1/(2^(k+2)+1) = 0.431166447... .

Then

    |B_0 intersect [1,N]| <= c N+O(log N).

In particular at least (1/2-c)N-O(log N) odd integers through N are absent.
This applies to every physical deletion-one source, not only to a selected
O/Z construction. The absent proportion among odd integers is at least
1-2c = 0.13766710... asymptotically. These constants are bounds, not the
actual density of B_0; collisions can only make the image smaller.

**Proof.** Write u=2^k m with m odd. For m>=3 the unit source with cover m
is strong, so the inherited physical cone gives g(m)>=m/2+1. Therefore

    g(u)+2u >= (2^(k+1)+1/2)m+1.

At most N/(2^(k+2)+1)+O(1) odd m can qualify for each k. Only O(log N)
values of k occur. The excluded m=1 cases are the powers of two and
contribute O(log N). Sum and bound the finite coefficient sum by c. QED.

A coarse rational version already suffices: the k=0 term is 1/5 and all
remaining terms total less than 1/4, so c<9/20. Thus at least N/20-O(log N)
odd integers are missed without evaluating the series numerically.

This is a direct boundary on the preceding session's exact signed-response
coverage: even though suitable aU-bD is surjective over this same d=1 class,
the positive readout D+2U has a positive-density exceptional set. Signed
surjectivity therefore cannot be converted into this positive coverage while
keeping deletion one. The obstruction disappears only after allowing larger
deletion counts; it does not refute other positive readouts or source pairs.

No corresponding obstruction for one, two or three J letters is proved.
The exact census below shows two J still missing large examples, but that
finite observation does not establish infinitely many two-J exceptions.

## 5. The explicit new arithmetic target

Let B_epsilon be the four-by-four matrix in (1), including the identity
for the empty word, let J be the displayed source matrix, and let
l=(1,0,2,0), t=(1,0,2,1)^T. Then the restricted scalar language is exactly

    union_(0<=r<=3) {
        l B_epsilon_r J B_epsilon_(r-1) ... J B_epsilon_0 t
        : each epsilon_i is an arbitrary finite binary word }.          (4)

The at most four binary words are chosen afresh for the target. They are
not fixed words, fixed slopes or a fixed extension of a supplied parent.
Every choice in (4) is physically realizable, by Sections 1-3; conversely
splitting a word at its J letters proves completeness for this restriction.

> **Four-block coverage lemma (OPEN).** Every odd n>27,781 belongs to (4).

This is a stronger, more explicit sufficient lemma than the unrestricted
balanced-cofinite statement. It asks for at most three root changes and
four arbitrary binary blocks, with the full exact bit-reversal and quadratic
coefficients retained. It makes no claim that a formal positive tuple or a
modular solution has a physical preimage.

**Conditional consequence.** If the four-block lemma holds, every odd n>=3
other than 7 and 21 is a Hamilton-path count of a strong tournament with
at most three Hamilton cycles. For n>27,781 use (4) followed by OJ, obtaining
C1. At smaller n, the inherited exact balanced bank and its physical repair
table apply. The finite census here independently supplies all odd n from
27,783 through 1,000,001 inside the proposed restricted class.

The remaining all-order obligation is precisely to solve the integer equation
(4) for every sufficiently large odd n. Calling the matrices a spin system,
a transducer, or a self-affine system does not discharge that obligation.

## 6. Complete finite experiment and controls

The finite universe is every source response obtained from t by any O/Z/J
word, with at most three J letters and G<=1,000,001. The inherited positive
coordinates prove that G strictly increases under every letter (indeed it
increases by at least a factor 3/2). Thus the target cutoff loses no accepted
word, regardless of graph order or length. Zero-one shortest-path relaxation
records the smallest J cost of each response and then of each scalar.
Discarding a higher-cost copy of an identical response is safe because every
future transfer depends only on that response, and costs are nonnegative.

The complete census gives 19,852,699 distinct nontrivial source responses.
With the singleton value 3 adjoined at cost zero, the numbers of represented
odd scalar values by minimum J cost are

    cost 0: 213,873; cost 1: 266,697;
    cost 2: 19,137; cost 3: 220.

There are exactly 73 gaps among the 500,000 odd integers from 3 through
1,000,001. Their largest is 27,781; their full list is in GAPS_THREE.json.
The same gap list already occurs at the cutoffs 50,001 and 500,001.
At most two J letters miss 293 values through 1,000,001, largest 879,001.
This is evidence for (4), not an induction or an all-order conclusion.

The checker independently compares the raw-response search with a positive-
coordinate layered closure at cutoff 50,001; checks (1) on a full coordinate
basis for every binary word through length ten; checks (2)-(3) through
u=20,000; checks the all-k unit-deletion prefix rebase through k=40,
including literal rooted isomorphisms through k=12; and tests literal source graphs, their caterpillar reversals,
balanced readouts and unique Hamilton cycles on all source words through
length four. Explicit exceptions remain enabled under python -O.

Run python3 verify.py and python3 -O verify.py serially. Each full run uses
several GiB of memory; do not run them simultaneously. The normal and
optimized receipts must agree. The earlier exploratory scripts are retained
separately and are not the verification entry point.

## 7. Historical bridges and their precise scope

This lane was prompted by three older threads rather than a new scalar census
alone. The source paths and limitations are recorded in CONNECTIONS.md:

* THM-853's closed continued-fraction return monoid distinguishes a short
  present summary from a true continuation state. Here the entire source
  response is retained when minimizing J cost; scalar G alone is not used
  as a recursive state.
* The Collatz parity/continued-fraction boundary note separates exact binary
  source guards from coarse growth exponents. Here every binary block has
  an explicit finite source recipe and the affine endpoint term in (1).
* THM-3435's dyadic-fibre analysis retains affine phase beyond fibre counts.
  Here opposite bit readings and their ordered pair correction play the
  analogous role, but no LRC mask or Collatz integer is identified with a
  tournament. The transported object is the finite-word transfer discipline.

The genuinely new physical object is the bounded-spine chronological
caterpillar, and the new tested arithmetic obligation is (4). The global
{7,21} coverage theorem remains open.
