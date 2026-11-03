# A smaller closed attachment carrier: pool internal vertices, retain parity

**PROVED all-order compression and closure; PROVED exact monomial-identification
quotient; FINITE-EXACT minimum 571 at order six with independent rank audit.
Minimum linear dimensions at orders seven and above remain OPEN.**

Source snapshot: `7b76044ca13e6684d1746923285b83dbf83828d2`.
The inherited theorem is
`05-knowledge/results/primegame_port_spectrum_20261103/FOREST_RESPONSE_RECURRENCE.md`.
Its all-order upper bound remains correct. Its exact minimum dimensions at
orders 3, 4, 5 remain 7, 28, 126. What changes is the next-order target:
**The exact minimum is 571 at order 6, not 661.** A closed 571-coordinate
carrier supplies the upper bound; 571 independently audited exterior
coefficient rows supply the matching lower bound. The exact dimensions at
orders 3, 4, 5, 6 are therefore **7, 28, 126, 571**.

## 1. Universe, operations, and inherited boundary

Fix named old vertices `[n]`. A source is an arbitrary rational signed measure
on all `2^binom(n,2)` labelled tournaments, or a probability law when positivity
is stated. For `i<j`, write `y_ij=a_ij-1/2`, and `y_ji=-y_ij`. For any enlarged
label set, the Hamilton polynomial is

    H_N(y) = sum_(vertex words w) product_i (1/2+y_(w_i,w_(i+1))).

The value on a Boolean orientation is the number of directed Hamiltonian
vertex words. An allowed attachment fixes every arc with at least one new
endpoint, independently of the old orientation, and leaves every old arc
and vertex name unchanged. There may be any finite number of new vertices.
Every such attachment is realizable by a sequence of fixed one-vertex cuts.
Updates, readouts, and dimensions here are rational linear ones. Source-
dependent attachments, legality conditioning, source copying, and products
of two responses on the same source are not single-response operations.

The previous coordinates `Q_pi` specify each path component's vertex set and
two endpoints, while forgetting the internal order with the correct signs.
Their number `b_n` has EGF `exp(x+x^2 exp(x)/2)`. We quotient further by
forgetting the allocation of internal vertices among components, retaining
only each component's length parity. This is not forgetting names: the set
of all internal vertex names is retained, just not their component allocation.

## 2. The signed parity signature

Let F be a labelled undirected linear forest on all old vertices, including
isolates. Define its signature `sigma(F)` by:

1. its set S of isolated vertices;
2. the perfect matching M of the degree-one vertices, pairing the two ends
   of each nontrivial path component;
3. for each pair in M, the parity epsilon of the number of internal vertices
   on that component.

The remaining old vertices form the pooled set U of degree-two vertices.
Thus S, the endpoints, and U are determined as named sets. Equivalently one
may retain parity of the edge count, which is `epsilon+1 modulo 2`.

Traverse each nontrivial path from its smaller named endpoint to its larger
one. Let `s(F)` be the product of +1 for a traversal along increasing labels
and -1 for a traversal along decreasing labels. Define

    R_sigma(y) = sum_(F with signature sigma) s(F) product_(e in F) y_e,

where each undirected variable `y_e` uses increasing-label order. Distinct
R coordinates have disjoint monomial supports, and every forest occurs in
one. Hence the R coordinates are linearly independent as polynomials and
as functions on the full Boolean cube. This is independence of the proposed
coordinates, not yet independence of their actual response columns.

Each R is the sum of the old Q coordinates with that signature. All forests
in one signature have the same degree `n-|S|-|M|`.

## 3. Exact full-Hamilton coefficient law

Let L be an edge subset of the complete graph on N vertices. Write `c(L)`
for the number of components including isolates and `r(L)` for the number
of nontrivial components. The coefficient of its squarefree monomial in H_N
is zero unless L is a linear forest with **every nontrivial path having an
even number of edges**. If this condition holds, the coefficient is

    s(L) * 2^r(L) * c(L)! / 2^(c(L)-1).                 (1)

For N=0 one may separately set H_0=1; formula (1) is used for N>=1.

Proof. Any contributing Hamilton word must contain every edge of L, so a
cycle or degree at least three is impossible. For a linear forest each
component must occur as a contiguous path block. Its two orientations
have opposite signs if its edge count is odd, causing cancellation. If its
edge count is even, both orientations have the same sign. There are `c!`
orders of the blocks, `2^r` choices of their orientations, and `c-1`
unselected edges between blocks, each contributing 1/2. This gives (1).
No assumption about a particular tournament has entered the calculation.

For fixed old vertices, let G consist only of edges with at least one new
endpoint. Expanding an attachment response in the exterior centered edge
variables gives old coefficient columns

    K_G(F) = [y_(F union G)] H_(n+k).                     (2)

Every actual attachment is an evaluation of these exterior variables at
+/-1/2. Conversely, all coefficients in the exterior variables are rational
linear combinations of these Boolean evaluations, by unique multilinear
interpolation. Thus the span of all attachment responses is exactly the
span of the rows (2), over all finite exterior sets and exterior edge sets.
More explicitly, if m is the total number of exterior-edge variables, the
coefficient polynomial p_G(y) is

    p_G(y) = 2^(|G|-m) sum_(s in {+1,-1}^m)
             (product_(e in G) s_e) H_(n+k)(y,s/2).           (2a)

This is ordinary Boolean multilinear interpolation, with no limit or
positivity assumption. Every summand is an actual fixed-cut attachment.
Some coefficient-extraction contexts are signed combinations of actual
counts, not individual graph constructors.

## 4. Compression theorem

**Theorem.** For F and F' with the same parity signature, and every exterior
edge set G,

    K_G(F)/s(F) = K_G(F')/s(F').                            (3)

Consequently every fixed attachment response is in the R span, at every
order and every attachment size.

Proof. An exterior edge touching an old degree-two vertex makes degree at
least three on both sides, hence both coefficients vanish. Otherwise the
old components act as two-ended blocks and the old isolates as named
singletons. Their endpoints and exterior incidence are identical. Thus
F union G is a linear forest if and only if F' union G is one: degrees and
cycles are determined by the same contracted incidence structure. Its new
path endpoint matching is identical. Each traversed old block contributes
the same length parity, so all enlarged path parities also agree. The
numbers of vertices, edges, components, and nontrivial components agree.
Finally, the sign of traversing an old block is its canonical sign times
`(-1)^(edge count)` if traversed backwards; this factor is determined by its
retained parity. Thus dividing out s(F) makes the enlarged canonical sign
identical on the two sides. Equation (1) now proves (3).

This proof also explains the information loss: an exterior context never
enters the interior of an old path selected in a monomial. Even the actual
numbers of internal vertices per component can change by two, provided all
old names remain in the pooled set and every component parity is unchanged.

## 5. Closed-update theorem

**Theorem.** The vector `(E[R_sigma])` is closed under every allowed fixed
attachment, fixed relabelling, and ordinary converse.

Proof. Expand an enlarged coordinate R_tau and substitute all fixed exterior
y values. The coefficient of an old forest monomial F is a sum over exterior
edge sets G for which F union G is a forest of signature tau, weighted by
its canonical sign and the product of those fixed y values. For F,F' with
the same old signature, the preceding contraction proof shows that the
same G are admissible, that the enlarged signatures agree, and that their
canonical signs divided by s(F) or s(F') agree. Thus each substituted R_tau
has constant signed coefficients on the old signature classes and belongs
to the old R span. Taking expectations gives a rational linear update.
Relabelling permutes named data and changes only canonical sign gauges.
Converse sends every y to -y; all terms in one R have the same degree.

The coordinate vector need not be nonnegative and is not a probability
vector on signatures. The kernel of this linear map is stable under the
specified updates, but need not be stable under conditioning or taking
shared-source response products.

## 6. Exact maximal monomial aggregation, not yet a rank theorem

The signature is the **coarsest possible pairwise identification of signed
forest monomials** for the entire exterior-coefficient language (2). In
particular two different signatures never have proportional nonzero
coefficient columns: an exterior G makes exactly one coefficient zero.

Proof of separation. All attachments used below use fresh, distinct vertices.

* If degree profiles differ, interchange F and F' if necessary and choose
  a vertex whose degree is lower in F than in F'. Attach enough fresh
  leaves to saturate its degree to two in F. This raises its degree above two in F'. Complete every nontrivial
  path of F union G to even edge count by adding a leaf at an available
  endpoint when needed. If the selected vertex had degree one, choose
  the repair at the opposite endpoint; if it had degree zero, its two
  leaves already form an even path. Hence F has nonzero coefficient and
  F' has zero coefficient.
* If degree profiles agree but endpoint matchings differ, take endpoints
  u,v paired in F but not F'. Add the two-edge exterior path u-z-v.
  This closes a cycle in F, but joins two distinct paths in F'. Complete
  the latter and all other paths to even edge count by leaf additions.
  Thus F has zero coefficient and F' nonzero.
* If degree profiles and matchings agree but parity decorations differ,
  extend each F path to even edge count with zero or one endpoint leaf.
  At least one corresponding F' path remains odd. Formula (1) separates.

These are valid coefficient contexts, each obtainable by a signed linear
combination of actual fixed-cut responses. The argument proves a precise
classification of equal/proportional columns. It does **not** exclude a
linear relation among three or more distinct columns. Consequently it
does not prove that the number of signatures is the minimum linear
response dimension. Such a claim requires a full-rank response certificate
or an additional general independence theorem.

## 7. Counting the new carrier

Let `c_n` count signatures. Choose `s=|S|`, `r=|M|`, and `t=|U|`, where
`s+2r+t=n`. For r>=1, a parity vector on the r endpoint pairs is feasible
exactly when its number j of odd entries satisfies `j<=t` and `j=t mod 2`.
Indeed allocate one internal vertex to each odd pair, then allocate the
remaining even number to any one pair. All labels in U can then be assigned
and ordered, so every such signature is realized by a forest.

Therefore

    c_n = 1 + sum_(s+2r+t=n, r>=1)
          n!/(s! t! 2^r r!)
          * sum_(0<=j<=min(r,t), j=t mod 2) binom(r,j).       (4)

The initial 1 is the all-isolated forest. There is no signature with
r=0 and t>0; imposing this boundary is essential.

For orders 0 through 12, old carrier dimensions b_n and new dimensions c_n:

| n | b_n | c_n |
|---:|---:|---:|
| 0 | 1 | 1 |
| 1 | 1 | 1 |
| 2 | 2 | 2 |
| 3 | 7 | 7 |
| 4 | 28 | 28 |
| 5 | 126 | 126 |
| 6 | 661 | 571 |
| 7 | 3928 | 2668 |
| 8 | 25628 | 13028 |
| 9 | 181828 | 65908 |
| 10 | 1396261 | 347626 |
| 11 | 11517606 | 1902231 |
| 12 | 101321320 | 10777900 |

At n=6 the old 661 endpoint-block coordinates split into 481 singleton
classes and 90 classes of size two, giving 571 new classes. Choose the two
pooled internal labels in 15 ways and pair the remaining four endpoints
in 3 ways. For each of the 45 choices, two old-coordinate collisions arise:
(1) one internal vertex on each path, with their assignments swapped;
(2) both internals on one path or both on the other. Both preserve each
path's parity. These are the first possible collisions, so the inherited
sharp dimensions 7, 28, 126 at orders 3, 4, 5 are unchanged.

For normalized source laws the constant coordinate is known, leaving at
most `c_n-1` varying coordinates. At order six exactly 570 are necessary,
not 660, by the independently audited rank certificate below. The rational
linear source-law ambiguity has dimension exactly `32768-571=32197`
(also 32197 within the normalized affine hyperplane).
These are linear dimensions, not state cardinalities or bit costs.

## 7a. Independently certified exact minimum at order six

The finite certificate `independent/rank_contexts.txt` contains 571 exterior
edge sets, each with its total vertex count and explicit edges. The old
vertices are 0,...,5; there are between 0 and 12 new vertices. Every edge
has at least one new endpoint. Their file SHA-256 is
`19f3e917d64e15d6cc79d15bf11d5f1fa3d3cac78215096226715cb21636987f`.

`independent/verify_exterior_rank.py` independently enumerates all 32,768
old edge subsets, reconstructs the 571 signatures and canonical
representatives, evaluates (1) for every representative plus each saved
exterior edge set, and divides its coefficient by the old canonical sign.
This is the matrix of the 571 coefficient polynomials in the R basis.

The resulting matrix has rank 571 over both prime fields:

- modulo 101: determinant 40;
- modulo 103: determinant 84.

All rational coefficients have denominators that are powers of two, so a
nonzero determinant modulo either odd prime proves the rational matrix is
invertible. Equation (2a) places each row in the span of actual fixed
attachments. Hence the response span has dimension at least 571. The
closed-carrier theorem bounds it above by 571, proving equality.

The independent verifier gives byte-identical normal and optimized
receipts, checks the coefficient law by direct Hamilton-word expansion at
orders 1 through 5, and drops to rank 570 after either deleting one row or
replacing one row by a duplicate. Its parser checks vertex ranges,
absence of loops, distinct undirected edges, and the absence of old-old edges.
The discovery checker separately checks every full-Hamilton coefficient through
order 6, without importing auditor code.

This is a finite exact rank proof at n=6. It is not an exhaustive listing of
the actual attachment evaluations in (2a), whose number can be large;
multilinear interpolation proves that span inclusion exactly. It uses all
fixed source-oblivious attachments, not a claimed small bank of words in
only S and the original B ports. It makes no all-order sharpness claim.

## 8. Explicit strictly positive new blindness at order six

Use labels 0,...,5, lexicographic pair-bit order, and forests

    F1 = {04,14,25,35},     F2 = {05,15,24,34}.

Their edge masks are 10376 and 5392. Their two components have endpoint
pairs {0,1},{2,3}; internal labels 4,5 are swapped. Both canonical signs
are +1. They belong to different old Q coordinates but the same new R.
Let `chi_F(a)=(-1)^(sum_(e in F) a_e)`. Define probability laws

    mu_+(a) = (1+(chi_F1(a)-chi_F2(a))/4)/32768,
    mu_-(a) = (1-(chi_F1(a)-chi_F2(a))/4)/32768.             (5)

Both are normalized and strictly positive, with weights 1/65536,
1/32768, or 3/65536. By (3), all fixed source-oblivious attachment counts
have the same mean under these two laws. But the old coordinate for F1
has expectations +1/64 and -1/64 respectively, since its degree is four
and `y_F1=chi_F1/16`. Thus some of the old 661 coordinates really contain
information no allowed mean response can reveal.

This blindness is neither cyclic-support nor branching-support blindness:
both changed Fourier supports are genuine disjoint unions of paths.
Shared-source second moments can distinguish the laws. The independent
actual-graph audit finds the two-cut context with masks 17 then 69 (a set
bit means the new vertex points to that old vertex) has common mean 627/2
but second moments 15043377/128 and 15038127/128, a gap of 2625/64.
This numerical second-moment witness is FINITE-EXACT; the all-context mean
identity follows from the general theorem, not from sampling contexts.

## 9. Connection contract and stopping boundary

Source: labelled-tournament laws under fixed attachment mean queries.
Target: the signed pooled-parity forest carrier. Map: sum old Q coordinates
with equal isolated set, endpoint matching, and component parity. Preserved:
all fixed response means and all permitted linear updates. Destroyed:
internal-label allocation and even changes of individual component lengths.
Required sidecar to recover old Q data: their allocation-resolving moments.
Cheapest decisive test: the six-vertex forests and positive laws (5).

The parity signature is a decorated matching with a pooled interior set.
This supplies a precise common combinatorial boundary object for matching-
based gluing. It does not identify tournament attachment means with tangle
isotopy, planar embedding, pseudoforest enumeration, fence area, product
packing, or multicolour Ramsey feasibility. Those problems need their own
operations, guards, and retained structure.

Closed here: a stronger all-order upper bound, closed updates, a sharp
pairwise monomial-aggregation classification, and the first strict loss.
Closed by the independent finite certificate: exact minimum response
dimension 571 at order six. Open here: exact minimum response dimensions
at order seven and beyond, and whether further linear relations arise
between distinct signatures at higher orders. The
well-known distinction between pairwise distinguishability and linear
independence is part of the stopping boundary, not a technicality.
