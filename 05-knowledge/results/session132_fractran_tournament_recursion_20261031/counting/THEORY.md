# One arithmetic counter, three tournament sequences

**Proved:** the all-order identities and finite recurrences below.
**Finite-exact:** `check_counts.py` computes orders 0 through 25 by distinct
integer algorithms, checks literal signed permutation actions through order
8, and enumerates every tournament through order 5. It independently
enumerates ordered pairs of tournaments, with and without vertex colors,
through order 4. `normal.json` and `optimized.json` are identical.

These are counting recurrences, not a claim that the naked metagraph node
count determines its wiring. The FRACTRAN program itself is a separate
session artifact. This file supplies a compilation target and its proof.

## 1. Definitions and inherited scope

Let A_n count n-vertex tournaments up to vertex isomorphism. Let S_n count
the classes isomorphic to their converse (all arcs reversed). Let M_n count
classes under vertex isomorphism **and** converse. Then

    A_n = A000568(n),                    offset 0,
    S_n = A002785(n),                    published offset 1,
    M_n = A059735(n) = (A_n+S_n)/2,       published offset 1.

We extend S_0=M_0=1 for the empty object. This is an explicit convention,
not an alteration of the published offsets. Complementary pairs in the
last sequence include singleton fixed orbits, not only chiral pairs.

The identifications and offsets were checked on [A000568](https://oeis.org/A000568),
[A002785](https://oeis.org/A002785), and [A059735](https://oeis.org/A059735).
The first page warns that the historical Harary--Palmer values at orders
24 and 25 were incorrect. The external controls here use the corrected
[A000568 table](https://oeis.org/A000568/b000568.txt).

The ordinary odd-cycle formula is classical; see Moon, *Topics on
Tournaments*, Section 29, especially the fixed-point lemma and Theorem 40.
Self-converse enumeration and anti-involution structure are likewise
classical; see Eplett, *Self-Converse Tournaments*, Section 1. Full citations,
retrieval limits, local control values, and source pins are in
`SOURCE_PINS.json` and `sources/`.

The repository already contains a more general defect-graded cycle grammar
in Session98 and parent-only re-rooting/projector enumeration in Session50.
Those files are byte-pinned here. This packet does not claim either mechanism
as a new discovery. Its purpose is an integer arithmetic specialization,
one parameterized engine for the three sequences, and independent checks.

## 2. Odd-cycle Burnside formula

For an odd-part partition lambda=(d_1,...,d_r) of n, put

    z_lambda = product_d d^(m_d) m_d!,
    E(lambda) = sum_i (d_i-1)/2 + sum_(i<j) gcd(d_i,d_j).

Equivalently,

    E = sum_d m_d(d-1)/2
        + sum_d d binom(m_d,2)
        + sum_(d<e) gcd(d,e) m_d m_e.

Then

    A_n = sum_(lambda odd-part partition of n) 2^E / z_lambda.       (1)

**Proof.** An even vertex cycle has an antipodal pair that a power of the
permutation swaps. An invariant orientation of that pair is impossible.
For an odd d-cycle, the unordered pairs internal to the cycle form
(d-1)/2 independent pair orbits. Pairs between distinct d- and e-cycles
form gcd(d,e) independent orbits. None reverses orientation on return.
Each pair orbit admits two choices, so the fixed-point count is 2^E.
There are n!/z_lambda permutations of this cycle type. Burnside averaging
gives (1), including the empty partition at n=0. QED.

For compilation there is no need for rational accumulators:

    N_n = sum_lambda (n!/z_lambda) 2^E,
    A_n = exact_divide(N_n,n!).                                    (2)

Every summand in N_n is an integer because n!/z_lambda is a conjugacy-class
size. The final division is exact by the orbit proof, not by rounding.

## 3. Twisted Burnside and a half-order engine

For a permutation pi, let tau reverse every arc. On an unordered-pair
orbit of length L, let s be 1 when pi^L swaps its endpoints and 0 otherwise.
The equation tau*pi(T)=T is consistent precisely when

    s = L mod 2.                                                  (3)

Each consistent pair orbit contributes two choices. Equation (3) implies
that the only permitted vertex-cycle types have at most one singleton and
all other cycle lengths congruent to 2 modulo 4:

* An odd cycle of length greater than 1 has an internal odd-length pair
  orbit with no endpoint swap, contradicting (3).
* Two fixed vertices give the same contradiction.
* For an even d-cycle, its antipodal pair orbit has length d/2 and swaps
  endpoints on return. Thus d/2 must be odd.
* These conditions suffice: all non-antipodal internal orbits have even
  length, as do cross-orbits between distinct allowed even cycles and
  between the singleton and an even cycle.

Writing n=2m+epsilon, epsilon in {0,1}, the permitted type is
(2d_1,...,2d_r,1^epsilon), where the d_i form an odd-part partition of m.
Its denominator is 2^r z_lambda and its number of free pair orbits is

    E_tw = sum_i d_i + 2 sum_(i<j) gcd(d_i,d_j) + epsilon*r
         = 2E(lambda)+(1+epsilon)r.

Define the rational polynomial

    F_n(q,x) = sum_(lambda odd-part partition of n)
                    q^E(lambda) x^r / z_lambda.                    (4)

Then the entire count family is

    A_n       = F_n(2,1),
    S_(2m)    = F_m(4,1),
    S_(2m+1)  = F_m(4,2),
    M_n       = [F_n(2,1)+F_floor(n/2)(4,2^(n mod 2))]/2.           (5)

Why does the twisted average count self-converse classes, rather than
converse-merged classes? Within a self-converse labelled orbit, each
tournament has exactly |Aut(T)| anti-isomorphisms and there are
n!/|Aut(T)| labelled representatives. Its contribution to the sum over
pi is therefore n!. Non-self-converse orbits contribute zero. Dividing by
n! counts exactly S_n. The further C2 orbit quotient gives M_n=(A_n+S_n)/2.

For q=2^k with k>=1 and an integer x>=1, F_n(q,x) has a direct interpretation:
it counts k ordered tournament structures on the **same** vertex set,
together with an x-coloring of the vertices, up to simultaneous relabeling.
The fixed-point count is (2^E)^k x^r. Thus S_(2m) equals the number of
simultaneous isomorphism types of ordered tournament pairs on m vertices;
S_(2m+1) adds two distinguishable vertex colors. Arbitrary numerical q
need not give an integer, so the checker does not assert that stronger claim.

The session's separate folding artifact develops the actual half-order
bijection using an anti-involution and an internally winning transversal.
The counting identity (5) here does not depend on that construction.

## 4. Explicit integer recursion with retained interaction state

The partition sum is finite but the following recursive form makes the
missing state especially clear. For every possible future odd cycle length
j, retain a nonnegative integer h_j. If previous cycles have lengths
d_1,...,d_t, its meaning is

    h_j = sum_(i=1)^t gcd(d_i,j).

Define B(0,h;q,x)=1 and, for k>0,

    B(k,h;q,x) = sum_(odd d<=k)
       [(k-1)!/(k-d)!] x q^((d-1)/2+h_d)
       B(k-d, (h_j+gcd(d,j))_(odd j<=k-d); q,x).                   (6)

Then

    F_n(q,x) = B(n,0;q,x)/n!.                                    (7)

**Proof.** In a permutation of the remaining k labelled vertices, the
cycle containing the least label has odd length d. Choosing its other
vertices and their cyclic order gives binom(k-1,d-1)(d-1)! choices, the
factorial quotient in (6). Its internal pair orbits contribute
q^((d-1)/2), its interactions with previous cycles contribute q^h_d,
and its color contributes x. Each future j-cycle acquires gcd(d,j)
new interaction orbits, exactly the displayed update. Recursing enumerates
every permitted permutation once and counts all pair choices. This proves
(6)-(7). QED.

Only nonnegative-integer operations appear before the final exact division.
No precomputed A_n, S_n, or M_n is used as a program input. A direct finite
register implementation may instead use (2): enumerate multiplicities
m_1,m_3,... with sum d*m_d=n, compute factorials and gcds by bounded loops,
accumulate the integer numerator, and exact-divide. Multiplication is
repeated addition, fixed-base powers are repeated multiplication, gcd is
Euclidean subtraction, and exact division is repeated subtraction with a
zero-remainder check. These observations specify a register-machine target;
the concrete ordered FRACTRAN fractions and their machine invariant must
still be supplied and verified separately.

The h-state is not optional: interactions depend on gcd of cycle lengths,
not only remaining vertex number. Equation (6) is not a scalar recurrence
in the bare sequence values. In the checked implementation, n=25 uses 142
odd partitions or 445 memoized cycle-insertion states. This is a finite
observation, not an asymptotic complexity or FRACTRAN runtime claim.

## 5. Structural recursion by strong components

Let C_n count strongly connected tournament types and D_n count the
self-converse strongly connected types. Set C_0=D_0=0. Let A(z),S(z),C(z),D(z)
denote their ordinary generating series. Then

    A(z) = 1/(1-C(z)),
    S(z) = [1+D(z)] A(z^2).                                      (8)

Consequently, for n>=1,

    C_n = A_n - sum_(k=1)^(n-1) C_k A_(n-k),
    D_n = S_n - [n even] A_(n/2)
          - sum_(1<=j<n, n-j even) D_j A_((n-j)/2).                (9)

**Proof.** The strongly connected components of a tournament are uniquely
totally ordered by dominance: opposite directions between two components
would force them into one strongly connected component, and a directed
cycle in the quotient would do the same. Thus a tournament is a unique
ordered word of strong tournament types. This gives the first identity.
Converse reverses this word and converses every letter. A fixed word has
an arbitrary left half, its reversed-converse right half, and either no
middle letter or one self-converse strong middle letter. This gives the
second identity. Extracting coefficients gives (9). QED.

These are genuine structural decompositions, but (9) does not independently
close a scalar recursion for A_n without strong counts. The companion
strong merged sequence is (C_n+D_n)/2, also checked through order 25.

## 6. Vertex growth: keep stabilizers, marks, and the converse sign

For each ordinary parent representative T on n vertices, let G_T=Aut(T),
with its actual permutation action. Adjoining a new vertex v whose
out-neighborhood in T is X gives a rooted extension. Precisely

    disjoint_union_[T] (2^V(T)/G_T)

is the set of (n+1)-vertex tournament types with a distinguished vertex.
It has size sum_[U] o(U), where o(U) is the number of vertex orbits of
Aut(U). It is not A_(n+1), and dividing by n+1 is generally wrong.
For example at n=2 the four rooted extensions give 4/3, not A_3=2.

One exact unpointing identity is

    A_(n+1) = sum_(rooted extension classes [(U,v)]) 1/o(U).        (10)

All root classes of a given U occur once, so its o(U) summands add to 1.
Equation (10) requires child orbit information. The inherited Session50
method is stronger algorithmically: parent-deletion transport constructs a
matrix R with R^2=(n+1)R and A_(n+1)=tr(R)/(n+1), without child-order
canonicalization. The source is pinned, not claimed anew here.

For converse-merged growth, use the signed stabilizer

    G_T^pm = {(g,epsilon): tau^epsilon g(T)=T}.

Its action on attachment subsets is

    (g,0): X -> gX,
    (g,1): X -> g(V(T)\X).                                      (11)

Choose one representative of each converse-merged parent class and quotient
its subsets by (11). These states are exactly rooted converse-merged child
classes. Formula (10) then holds with signed vertex-orbit counts and M in
place of A. For a chiral parent there are no sign-1 stabilizers; for a
self-converse parent they form an Aut(T)-torsor. Replacing all of this by
the unsigned group loses the sheet identification already at n=1: the two
attachments are converses and must merge.

The checker verifies both extension identities, each child's root
multiplicity, signed complement action, and the anti-isomorphism torsor for
every tournament through order 5. These states retain object-level data;
none of the count formulas recovers line colors, multiplicities, loops, or
metagraph adjacency by itself. M_n counts the repository's converse-merged
tournament nodes, but not its tilings or lines.

## 7. Reproduction and exact coverage

    python check_counts.py > normal.json
    python -O check_counts.py > optimized.json
    cmp normal.json optimized.json

The checks deliberately use explicit exceptions, not removable `assert`s.

* A_n: partition formula and integer cycle insertion for 0<=n<=25;
  external corrected OEIS controls for every one of these 26 values.
* S_n: half-order formula, integer recursion, and direct doubled-part
  twisted formula for 0<=n<=25; external displayed controls for 1<=n<=23.
* M_n: converse quotient through n=25; external displayed controls for
  1<=n<=19. The inaccessible S/M b-files are not falsely counted as checked.
* Literal signed permutation actions on unordered pairs through n=8,
  including every one of the 40320 permutations at n=8, independently
  validate both admissibility rules and fixed-point counts.
* Complete labelled orbit census through n=5, including all 1024 labelled
  five-vertex tournaments, validates A,S,M,C,D and automorphism histograms.
* Direct simultaneous pair-tournament census through m=4 yields
  F_m(4,1)=1,1,2,12,176 and F_m(4,2)=1,2,8,88,2752.
* Unsigned and signed pointed extension checks run from n=0 through n=4.

For orientation, at n=0..10 the three main rows are

    A: 1,1,1,2,4,12,56,456,6880,191536,9733056
    S: 1,1,1,2,2,8,12,88,176,2752,8784
    M: 1,1,1,2,3,10,34,272,3528,97144,4870920.

All 26 rows, the strong companions, and work counts are in `normal.json`.
