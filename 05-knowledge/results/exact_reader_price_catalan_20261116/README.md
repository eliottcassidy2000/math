# Exact reader prices, the 8/9 carry, and prime selection

Status: all-order reader and arithmetic statements below are proved; the
particular masks and local spectra are finite exact, with **no size cutoff
on the extensions** in the stated absolute interval. No global spectrum
or Goldbach/twin-prime implication is claimed.

## A. Pay the contact before trying to land

Fix a nonempty strong tournament T with count h=H(T). Write the edges
from old vertices to the new apex as X. For old vertices a,b, let w_ab
sum H(T-V(P)) over all positive even-vertex directed paths P from b
to a. Here H(empty)=1. The exact increment is a nonnegative directed cut:

    H(T+x)=h+2 K_T(X),   K_T(X)=sum(w_ab : a in X, b outside X).

Indeed, in the inherited odd-cycle packing identity, the packs avoiding
x give h. In the rest, the unique cycle using x is x,b,...,a,x, with
weight 2; all remaining disjoint cycles give the induced complement count.
This is the positive directed-cut presentation of the inherited quadratic
attachment law (THM-2527 and the one-head bank work). It is not counting
Hamilton paths by unweighted adjacency powers. Testing whether a
particular mixed cut releases p is arithmetic; deciding whether its output
can fit an absolute budget uses its actual nonnegative capacity.

**Subtripling grammar.** The entire spectrum of tournaments containing T
as an induced subtournament and having H<3h is obtained as follows.
Start at T, repeatedly add a vertex with a nonempty proper cut, stop
before crossing 3h, and allow neutral ordinal padding. In particular a
finite exhausted search of this grammar proves an all-orders nonhit
for a prescribed N<3h relative to this T.

Proof: let C be the strong component containing T. All other SCCs must
be singletons: any nontrivial SCC has H>=3 and SCC counts multiply.
Build up the current strong subset S inside C starting with T. If it is
still proper and no vertex mixes with S, then S is a homogeneous module
of C. Its proper contraction is a nonsingleton strong tournament:
project any directed reachability path to the quotient and erase repeated
vertices. It has H>=3. Expanding each Hamilton path of the quotient by
a Hamilton path of S gives at least 3H(S)>=3h, a contradiction.
A mixed vertex preserves strong connectivity and increases H strictly:
some old arc goes from outside its cut into the cut, furnishing a
two-vertex positive even path. The general induced monotonicity needed
for the cap follows by assigning unchanged weights to all odd-cycle
packs on the old vertices. The converse uses every displayed mixed
extension and padding is optional. This proves the grammar. In the
implementation subsets are labelled, not deduplicated by isomorphism.

One immediate practical lower bound: an extension that escapes a prime
dividing h costs at least min(3h, h+2 min_{proper X}K_T(X)).
If the cost-minimizing mixed cut escapes p and lies below 3h, that is
the exact first release among extensions of *any order* of this host.
The inner minimum of the ordinary capacities is a directed global min-cut.
The minimum for a *specified residue* need not be the same; after
filtering by a capacity or residue, exact hitting still requires a
physical cut. General congruency-constrained submodular minimization
for cardinality does not directly solve congruence of the *cut value*.

Examples, with explicit masks in the verifier:

| Source | Divisor | Cheapest release keeping the core induced, any order | First homogeneous release in exterior-vertex budget |
|---|---:|---:|---|
| rigid, permanently opaque H=49 | 7 | 83 | never |
| second rigid H=49 | 7 | 85 | 873 at 2 |
| rigid H=63 | 3 or 7 | 103 | never for 3; 1935 at 2 for 7 |
| rigid H=343 with delayed homogeneous release | 7 | 579 | at 4; matching reader 112983 |

The H=49 locked host reaches 83 and 85 but **not 87 or 89** below its
tripling threshold; 89 is prime, so release does not imply even a
tiny prime-target window. These exclusions are about this host, not
any tournament realizing 89. A second wave below 147 is needed: for
the first host count 143 has two-attachment witnesses.

## B. Repeated transparent triangles: exact vertex budget, large count price

Let Q_T(t) be the path-colour polynomial (THM-3166): distribute vertices
to t named colours and multiply the induced path counts in the colour
cells (empty count 1). An ordered cover with j nonempty runs contributes
to Q_T(t) with coefficient binom(t,j). Let U(T)=C3[T,1,1]. Then

    Q_U(T)(t) = A Q_T(t),   A Q(t) = t^2 Q(t)+t(Q(t+1)-Q(t-1)).

When the new vertices have distinct colours the contribution is
t(t-1)Q(t), as they simply extend their respective paths. With the
same colour the old cell count is unchanged in the one-run portion;
the added split portion is a pair of independently coloured cells,
except both empty. Summing gives t[Q(t)+Q(t+1)-Q(t-1)].

Thus H(U^m(T))=(A^m Q_T)(1), involving Q(1),...,Q(m+1) with top
coefficient m! by the unique all-up walk and the zero row at t=0.
Writing D_j=(j-1)! times the number of ordered j-run covers, the walks
are lower unitriangular in D through depth m+1. The lower-coordinate
integrality follows from freely permuting internal nonempty exterior
gaps. Any homogeneous reader with at most 2m+1 exterior vertices
only reads these coordinates, because each internal gap uses at least
two vertices. Therefore actual nested triangles match the best
possible homogeneous *vertex* budget for detecting any residue
disagreement or releasing any prime, even though they need not
minimize the counts. (The 343 host's second wrap reaches 118639.)

For odd p, the signs Q(-t)=(-1)^n Q(t), and Q(0)=0 at n>0, follow
from the odd-cycle pack polynomial. The first (p-1)/2 evaluations
determine the function mod p, and the leading m! is still a unit.
Every releasable p has a homogeneous release by p-3 added vertices;
for 7, four is sharp by the 343 example. For the first four primes:
2 is excluded by count parity; 3 has one channel and cannot release
once in the core; 5 has at most two channels and costs at most two
vertices; 7 has three channels and costs at most four vertices.
These are individual tests, not simultaneous-release guarantees.

Modulo 7 the coordinates Q(1), Q(2), Q(4) are indexed by the doubling
orbit 1,2,4, representatives of the sign pairs in F_7^*. This is
the exact same exponent cycle as Frobenius on F_8^*, hence on the
seven nonzero points used in Fano/octonion models, not an octonionic
or Bott model of the reader itself. For odd n, on coordinates
Q(1),Q(2),Q(3), A has rows (1,1,0),(-2,4,2),(0,-3,6) modulo 7.
Its square is (1,3,1)^T(6,5,2), and A acts on that line by 4.
After two wraps the count is either forever 7-divisible or forever
7-free. Even n changes the terminal diagonal to 5 and the matrix
is invertible mod 7. There is a genuine parity bifurcation, not a
sin/cos analogy promoted to period eight. The odd-order H75 example
with vector (5,2,1) goes to H1155 and permanently loses the ability
to escape 7 homogeneously. It also illustrates why a reader can
destroy a useful phase. In fact mod 3, A on Q(1) is 1+(-1)^n, so every
wrap of any odd-order core puts it into the permanent mod-3 channel.

## C. A literal 8 to 9 carry and its cutoff

On an ordered spine of five old vertices 0,...,4, take X={4}.
The weights landing at 4 are 4,2,1,1. Grouping starts {0,2} and
{1,3} gives 8=5+3. Allow vertex 2 to beat the apex as well. The old
landing port loses one path and now contributes 7; the new landing
port contributes two, giving 9=7+2. Thus this gives twin counts 17/19
on the same host with one flipped arc; the decomposition changes
from groups on one port to two individual ports. It is not a
statement about uniqueness of that grouping.

The same move has an exact cutoff on a transitive spine: for initially
X={n-1}, switching i to the cut side adds 2^(i-1) (i>=1) and
takes away 1 if n-1-i=1, otherwise 2^(n-3-i).
Reason: for b<a the even paths from b to a select an even subset
of the intervening vertices, so w_ab=1 for adjacent positions or
2^(a-b-2) otherwise; everything else is zero. For i=0 the deposit
is zero. A net +1 requires 2-1, so i=2 and n=4 or 5 only.
The other witness is 9/11 at n=4. A repeated shallow local bit
is not an unbounded prime seam; it needs a changing reader or a
coordinated edit which preserves a small signed discrepancy.

## D. Arithmetic mirrors and a conditional two-seam route

See ARITHMETIC.md for the older summand/multiplicand, Catalan and
Zsigmondy context, a wheel conjugacy, and the exact Goldbach
conditional. After ordinal products reduce the target atoms to
p,7p and 49,63,343, all large atoms are +1 or -1 mod 6. Two seams
p+q+1 and p+q-1 with compatible ready primes of opposite types
would therefore suffice, assuming compatible Goldbach for large
multiples of six and the finite startup cases. A fixed finite bank
of left primes cannot do this even if all its seams work: CRT and
Dirichlet force infinitely many prime targets and infinitely many
7p targets on which every offered right prime is composite.
Thus the inverse step must choose an adaptive left prime, or admit
a nonprime right source. The Goldbach/twin wheel equivalence carries
local obstructions, not primality, exact port readiness, or a
computability reduction. In fact there is an elementary pointwise
decision algorithm without a sharp birth bound. Strip singleton SCCs
(their count factors are 1). Every remaining vertex lies in a directed
triangle: otherwise all of its in-neighbours beat all of its
out-neighbours, impossible within a strong component. Hence the
remaining vertex count is at most 3 times the number of triangles,
which is at most 3(H-1)/2 by the odd-cycle identity. Enumerate to this
finite bound (handle H=1 separately). No computable reduction of the
halting problem to the yes/no question at a *single* supplied count
can exist. That observation says nothing about independence or the
truth of the universal statement.

The 8/9 perfect-power collision is Mihailescu's theorem; see
[Mihailescu's account](https://www.math.kth.se/4ecm/program/mihailescu.pdf).
The tournament 7/21/63 relation here uses cyclotomic arithmetic, not
an inference from Catalan to Chan–Pak.
