# Elliptic contacts and the first mixed eleventh root interfaces

2026-10-16. **PROVED, scoped and self-audited:** the infinite contact-tree
theorem, the low-width Cartesian theorem, the two-port linkage theorem,
and the common-parameter specialization theorem below.
The real quadratic-seed restriction G6 and its odd cubic extension repair
G7 are also proved in their stated scopes.
**FINITE-EXACT:** the declared geometric and elliptic controls in
the companion receipt. **INHERITED:** the coordinate decoder and the
stronger Collatz rank obstructions. Six-colour forcing, the uniform paired
reset bound, universal ROOT entry, and planar Jacobian remain **OPEN**.
No independent referee or global chromatic improvement is claimed.

## Inheritance and the live objects

The closest proved mechanism is the positive-rank contact curve in
[recursive interfaces, R2b](recursive_interfaces_20261016b.md). Its hostile
is the 49-point induced tree: infinitely many eligible contact parameters
need not supply compatible cross-edges. The corrected near miss is forgetting
that those elliptic contacts still lie in a smaller field. The stronger
[five-colour ceiling and odd-degree barrier](rhombus_overlap_elliptic_20261011.md)
and [one-layer port theorem](eleven_port_layers_20261011.md) remain in force.
The least-used coordinate here is the width in an eleventh-root basis over
the smaller field, together with local denominators of rotation ratios.

Anchor: contacts that can change a marked boundary relation. Niche: the
positive fruit cubic and its coordinate decoder. Wildcard: Paley periods
and the first field-coordinate interaction. The concept board consists of
elliptic contact parameters, quadratic phase factors, field width, joint
port relations, and literal Collatz clocks.

Write

\[
 L=\mathbb Q(\zeta_{55}),\quad K=\mathbb Q(\zeta_5,\sqrt{-11}),
 \quad j=\sqrt{-11},\quad k=\sqrt{-55},\quad \xi=\zeta_{11}.
\]

Thus `[L:K]=5`. Golden arithmetic supplies the real quadratic factor:
`jk=-11 sqrt(5)`. The old pentagon and rhombus ports lie in K and are
retained literally, in their original order, in the finite linkage test.
The infinite tree theorem concerns its explicitly named vertices, not a
claim that the old eight ports are vertices of that tree.

At startup and the first research checkpoint, fetching `origin/main`
found no incoming commits beyond `7a2b8a743a`. The earlier incoming
[dense census preregistration](collatz_census_2p42_preregistration_20261016.md)
still supplies proposed tests, not a new uniform bound. The inheritance
pass recovered the September fruit decoder and the stronger two-place
rank obstruction before proposing another height correction.

## 1 The two quadratic phases cannot mix freely

**G1 (PROVED).** Let `A,C` be nonzero elements of `Q(j)` and let `B,D` be
norm-one elements of `Q(k)`. If `|AB-CD|²` is rational, then either
`A/C` is rational or `B/D` is rational. In the latter case `B/D=+1` or `-1`.

Indeed write `A bar(C)=h0+h1 j` and `B bar(D)=t0+t1 k`, with rational
coefficients. Their product and its conjugate sum to

\[
 2h_0t_0+2h_1t_1jk.
\]

The two individual squared norms are rational. Since `jk` is an irrational
real multiple of `sqrt(5)`, rationality of the distance forces `h1 t1=0`.
Dividing by the nonzero rational norms gives the stated alternatives.
This is a necessary condition, not a sufficient unit-edge test.

In particular, if `|A|²/|C|²` is not a rational square, only the second
alternative remains. The extra real coordinate supplied by golden
arithmetic is an exact obstruction to mixing the two phases; discarding
it would invent cross-edges.

## 2 The earlier finite tree extends to an infinite tree

Set

\[
 a=(1+j)/2,\qquad g=(5+j)/6,
\]

and take the two points of the contact quartic already proved to exist:

\[
 (r_1,y_1)=(1/5,8/5),\qquad
 (r_2,y_2)=(91/135,10552/3645).
\]

For `t=1,2` and `sigma=+1,-1`, define

\[
 s_t=(1+55r_t^2)/4,\qquad
 b_{t,\sigma}=\frac{s_t+2+\sigma(y_t/4)j}{2\bar a},
 \qquad d_t=\frac{55r_t^2-1+2r_tk}{1+55r_t^2}.
\]

The inherited identities give `|a|²=3`, `|b_(t,sigma)|²=s_t`,
`|a-b_(t,sigma)|=1`, and `|b_(t,sigma)-d_t b_(t,sigma)|=1`.
Here `s1=4/5` and `s2=23684/3645`.

**G2 (PROVED).** On the vertex set

\[
 \{g^i a:i\in\mathbb Z\}
 \cup
 \{g^i d_t^h b_{t,\sigma}:i,h\in\mathbb Z,
          \ t=1,2,\ \sigma=\pm1\},                 \tag{1}
\]

all vertices are distinct. Its **complete induced unit-distance graph**
has exactly these edges:

* consecutive spine vertices `g^i a, g^(i+1) a`;
* `g^i a` joined to `g^i b_(t,sigma)`;
* consecutive vertices along each fixed `(i,t,sigma)` branch, with `h`
  changing by one.

It is an infinite tree: four doubly infinite paths, each attached once,
at every vertex of a doubly infinite spine. In particular it is bipartite.
An explicit colouring is `i mod 2` on the spine and `i+h+1 mod 2` on a
branch. No enlargement within this entire family can make it four-chromatic.

**Proof.** The only roots of unity in each of `Q(j), Q(k)` are `+1,-1`.
The displayed nonreal norm-one rotations therefore have infinite order.
On a fixed circle a point has precisely two possible unit neighbours;
the displayed rotations give them. This proves the claimed spine and
within-branch edges.

The ratio `3/s_t` is not a rational square for either t. G1 forces any
spine-to-branch edge to have `d_t^h=+1` or `-1`, hence `h=0`. For `h=0`
the equation `|a-z b_(t,sigma)|=1`, `|z|=1`, has the two rotation roots
`1` and

\[
 \lambda_{t,\sigma}
 =\frac{s_t+2-\sigma(y_t/4)j}{s_t+2+\sigma(y_t/4)j}.
\]

Their nontrivial traces, independent of sigma, are respectively

\[
 \lambda_{1,\sigma}+\bar\lambda_{1,\sigma}=19/15,
 \qquad
 \lambda_{2,\sigma}+\bar\lambda_{2,\sigma}
       =110354899/64746135.                         \tag{2}
\]

Both have a denominator divisible by 5 in lowest terms. Every power
`+/-g^m` is integral at all primes over 5, so has trace integral at 5.
Consequently neither lambda is `+/-g^m` for any integer m. Only the root
1 is attained, giving the indicated attachment with the same i.

For two branch vertices on the same t-circle, a unit edge makes their
ratio equal `d_t` or `d_t^-1`. Their ratio also factors as an element of
`Q(j)` times a power of `d_t`. The intersection of `Q(j)` and `Q(k)` is Q,
so the first factor must be `+1` or `-1`. When the signs sigma agree,
infinite order forces equal i. When they differ, that first factor is a
power of g times a lambda from (2), so it cannot be either sign. This
proves that distinct branches on the same radius do not meet or cross.

The two different radii have separation greater than one. An entirely
rational check is `s2>s1+1` and

\[
 (s_1+s_2-1)^2-4s_1s_2=250681849/13286025>0.
\]

Thus there are no unit edges between those circles at any angles.
Distinctness follows from the same phase argument and the three distinct
radii. Connectivity and the tree description now follow from the complete
edge classification. QED.

The finite receipt tests both signs, `-2<=i,h<=2`, rather than only the
positive signs of the earlier patch. It checks every pair, with exact field
norms on every survivor of two exact modular filters.

## 3 Three eleventh root coordinates still give a Cartesian graph

**G3 (PROVED).** Suppose K is a conjugation-stable subfield of C, xi has
norm one and degree d over K, and `2m<d`. If

\[
 v=\sum_{h=0}^m a_h\xi^h,\quad a_h\in K,\quad v\bar v=1,
\]

then exactly one coefficient is nonzero and that coefficient has norm one.

Multiply `v bar(v)-1` by `xi^m`. It is a polynomial of degree at most
`2m<d`, so its vanishing at xi means it is identically zero. If two or more
coefficients were nonzero, the coefficient at the largest positive Laurent
exponent would be the product of the last coefficient with the conjugate
of the first, and could not vanish. The constant coefficient then gives
the remaining norm. QED.

For `[L:K]=5`, the map

\[
 K^3\longrightarrow L,\qquad (x_0,x_1,x_2)
       \longmapsto x_0+x_1\xi+x_2\xi^2              \tag{3}
\]

is injective and identifies the complete unit graph on its image with
the Cartesian product of three copies of the unit graph on K. An edge
changes exactly one coordinate by a unit. This is an all-coordinate
classification, not a generic-angle statement.

The inherited `chi(K)<=5` therefore gives a five-colouring of this whole
image: sum any three K-colours modulo five. More precisely, for finite
coordinate sets `A0,A1,A2`, with `0 in A1,A2`, the full product has exactly
the same q-colour boundary relation on `B subset A0 x {0} x {0}` as A0,
provided A1 and A2 are q-colourable. Extend an arbitrary valid colouring
of A0 by adding the other two colourings modulo q, normalized to zero
at their marked zeros. Restriction supplies the reverse inclusion.
For an arbitrary subset of that product this gives a sufficient lift,
not an equality with its own projection relation.

Thus a mixed edge in this fixed basis needs exponent width at least three.
Width three is a necessary threshold, not an existence theorem. The
norm-one positive control `(2-xi)/(2-xi^-1)` has width four: its coefficient
of `xi^4` in the quintic basis is `12/(34+9j)`, which is nonzero.

This threshold is a real connection with the tournament calculation.
The quadratic-residue orbit `{1,3,4,5,9}` of the Paley tournament produces
the minimal polynomial

\[
 z^5+(1-j)z^4/2-z^3+z^2-(1+j)z/2-1.                \tag{4}
\]

Its degree five is the degree used in G3. Golden arithmetic adds the
other quadratic factor, while the remaining degree-five extension
governs when coordinate interactions can start. No arbitrary orientation
is imposed on the geometric contact graph. The tournament supplies a
Galois representation, not a chromatic forcing certificate.

## 4 A two-port linkage enters the larger field but has only two seam edges

A direct rotation about a K-centre has no noncentral cross-edge to a K-piece
when its unit multiplier has degree five: the cross-edge equation is
quadratic over K. To escape that restriction, move the second piece as
well as rotating it, while keeping the first piece and its ports fixed.

Let `a in K`, `a!=0`, and `A=|a|²!=1`. For a norm-one xi of degree at least
five over K, set

\[
 u=\frac{a-\xi}{\bar a-\xi^{-1}},\qquad
 F_\xi(q)=\xi+u\bar q.                             \tag{5}
\]

The denominator is nonzero because xi is not in K. We have `|u|=1`;
F is a reflected isometry. The two exact unit attachments are

\[
 |F_\xi(0)-0|=1,\qquad |F_\xi(a)-a|=1.
\]

**G4 (PROVED).** For **all** `p,q in K`,

\[
 |p-F_\xi(q)|=1
 \quad\Longleftrightarrow\quad
 (p,q)=(0,0)\text{ or }(a,a).                     \tag{6}
\]

Also `F_xi(K)` is disjoint from K. Thus one such linkage never creates
an unaccounted seam, regardless of how many elliptic points are placed
in either K-piece.

**Proof.** Multiply `p-F_xi(q)` by `bar a-xi^-1`. Its numerator is

\[
 R\xi+S+T\xi^{-1},\quad
 R=\bar q-\bar a,\quad S=p\bar a+1-a\bar q,
 \quad T=-p.
\]

The unit equation after multiplication by `xi²` has degree at most four.
Since xi has degree at least five, all its coefficients must vanish.
The outer coefficient is `R bar(T)`, so `p=0` or `q=a`.

For `p=0`, put `q=ah`. The coefficient of xi gives
`bar h+A h-A|h|²=0`. As `A!=1`, subtracting its conjugate shows h is real.
If h is nonzero it gives `h=(A+1)/A`. The constant coefficient gives
`h=4/(A+1)`, so `(A-1)²=0`, a contradiction. Therefore q=0.

For `q=a`, put `p=ah`. The xi coefficient gives
`A|h|²+(1-A)bar h=1`, hence h is real and
`(h-1)(Ah+1)=0`. The alternative `h=-1/A` violates the constant
coefficient unless `(A-1)²(A+1)=0`. Thus h=1 and p=a.
Both advertised pairs do give unit distances by (5).

A collision would make the numerator itself zero. Its three coefficients
force `p=0,q=a,1-A=0`, again impossible. QED.

The exclusions `a=0` and `|a|=1` are necessary: the transformation then
degenerates into a rotated reflected K-piece and can have many seam edges.
Neither case is silently covered by (6).

**Boundary consequence.** A finite q-colourable new piece attached by only
these two edges cannot restrict an already coloured old piece for `q>=3`.
Among all q! permutations of a new-piece colouring, each seam edge forbids
exactly `(q-1)!`; two such sets cannot exhaust q!. This remains true for
arbitrary old boundary assignments that already extend to the old piece.
It does not assume that an arbitrary assignment to eight ports extends to
all of K.

The new test therefore composes five linkages, for
`xi,xi^3,xi^4,xi^5,xi^9`, retaining **every** edge between the new pieces.
These are the conjugates fixing K. The K-seed includes the eight original
ports, both signs of both elliptic parameters with `-1<=i,h<=1`, and the
unit path

\[
 0,g,2g,3g,3g-1,3g-2=a.
\]

The saved receipt contains the coordinates, the source and copy indices,
every edge, the ten source-to-copy seams, and a colouring checked on every
edge. This assembly has **306 vertices and 322 edges**, with no collisions
between pieces and no cross-edges between the five new pieces. Its
51-vertex seed has 52 edges; the six copies contribute 312 edges and the
ten seams supply the rest. The source-to-copy seams are additionally checked
by direct norms on every pair, independently of the modular filters.

## 5 A shared linkage parameter can be specialized back into K

**G5 (PROVED).** Fix finitely many anchors `a in K` with `a!=0,|a|!=1`,
and finitely many source points in K. Form any finite labelled assembly
from those source points and the points `F_(a,xi)(q)` of (5), with one
**shared** norm-one parameter xi of degree at least five over K. Assume
K contains a nonzero purely imaginary element j.
There is a norm-one parameter w **in K** giving the identical labelled
unit graph, including precisely the same coordinate identifications,
while fixing every original K-port pointwise.

**Proof.** For an old-to-new pair the edge numerator was computed in G4;
after multiplying by `xi²`, its norm equation has degree at most four.
For a new-to-new pair with possibly different anchors a,b, subtract the
two expressions; the initial xi cancels. After multiplication by
`(bar a-xi^-1)(bar b-xi^-1)`, the numerator is

\[
 (a-\xi)(\bar b-\xi^{-1})\bar q
 -(b-\xi)(\bar a-\xi^{-1})\bar r.
\]

This has Laurent exponents only `-1,0,1`. Its squared norm and the
squared norm of the denominator have exponents only `-2,...,2`.
The unit equation again has polynomial degree at most four. Thus every
edge at xi is a polynomial identity and remains an edge at every
admissible norm-one w. Every nonedge excludes only finitely many w.
Collision equations have degree at most two after clearing a Laurent
power; the same reasoning preserves identifications and avoids new ones.

There are finitely many pairs and denominators. The infinitely many
distinct values `(1+t j)/(1-t j)`, for rational t, lie in K and have norm
one. Choose one outside the finite exceptional set. All old coordinates
were fixed throughout. QED.

Consequently, arbitrarily many elliptic anchors with a single shared
parameter do not yield a graph or boundary relation unavailable over K.
This is stronger than observing that each individual linkage is harmless.
It does **not** apply to the several conjugate parameters `xi^h` together:
substitution increases the degree, so vanishing need not be a polynomial
identity. Nor does it say that K has no nontrivial six-colour boundary
relations. Its five-colour ceiling concerns chromatic number, while a
prescribed six-colour port assignment is a different predicate.

The finite control specializes the two-anchor shared-parameter assembly
and compares its entire labelled graph, not just its edge count. The
foreign assembly has **152 vertices and 159 edges**. Choosing t=1 is a
hostile: it creates a collision and additional edges, giving 151 vertices
and 176 edges. Choosing t=2 gives exactly the original labelled graph,
with all eight original ports unmoved. Thus the finite-exception clause
performs real work. The separate two-anchor/five-conjugate test escapes
G5's shared-parameter hypothesis. With anchors `a,b_(1,+)`, it has
**556 vertices and 587 edges**, and an exact three-colouring. Five pairs
of presentations coincide at the expected points `xi^h`; no new edge
joins vertices that do not already belong to a common piece, apart from
the source seams. Thus changing the anchor creates actual shared vertices,
but does not yet yield a forcing context. By
[THM-4614, the fixed-port extension theorem](../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md),
both tested assemblies preserve all 570 proper six-colour port patterns.

### Real quadratic seeds miss the width three target

**G6 (PROVED).** Put `K+=K intersect R`, of degree four over Q. For
`a,b in K+`, consider the norm-one element

\[
 v=\frac{a+b\xi+\xi^2}{a+b\xi^{-1}+\xi^{-2}}.
\]

If its coefficient of `xi^4` in the quintic basis vanishes, then it is
one of the monomials `xi²,-xi²,xi³,-xi³`. Thus allowing arbitrary golden
and other real K-coefficients in this quadratic seed does not produce
the desired nonmonomial of width three.

**Proof.** Multiply by the denominator and use (4). The multiplication
matrix is invertible: a nonzero polynomial of degree at most two cannot
vanish at xi. The Cramer numerator of the coefficient of `xi^4` is
`R+jI`, where

\[
 R=-\tfrac12(a-1)(a+b+1)
     (2a^2b-a^2-2ab^2+3ab+2a-2b+2),
\]
\[
 I=\tfrac12a(a-1)(a+b)(a-b+1).
\]

Both R and I are real, so both must vanish. The branch `a=1` gives
`v=xi²`. The branch `a=0` requires `b=+/-1` and gives `+/-xi³`.
The branch `b=a+1` leaves `a=0,1,-1`; the additional point `a=-1,b=0`
gives `-xi²`. Finally `b=-a`, apart from the already handled `a=1`,
requires `2a³+2a²-2a-1=0`. This cubic is irreducible over Q by the
rational-root test and cannot have a root in the degree-four field K+.
These exhaust the factors of I. QED.

The script verifies the entire symbolic numerator identity, rather than
inferring it from a coefficient scan. Complex coefficients in K, higher
degree seeds, or another support pattern are outside G6. This identifies
a specific missing phase to retain in the next width-three search.

### An odd cubic extension supplies a mixed direction and preserves the ports

**G7 (PROVED; alternative carrier).** Let theta be the real root in `(4/5,1)`
of

\[
 2\theta^3+2\theta^2-2\theta-1=0,
\]

and put `K'=K(theta), L'=L(theta)`. This root exists and is unique in that
interval by the endpoint signs and the positive derivative there. The
cubic is irreducible, so `[K':K]=[L':L]=3` and `[L':K']=5`.
The explicit unit

\[
 v_\theta=
 \frac{\theta-\theta\xi+\xi^2}
      {\theta-\theta\xi^{-1}+\xi^{-2}}             \tag{8}
\]

has **exactly exponent width three** in its quintic basis over K'.
It is a genuine mixed direction, rather than a monomial in xi.

**Proof.** The numerator and denominator are conjugate and nonzero, so
the norm is one. Substitute `a=theta,b=-theta` in the Cramer identities:
I vanishes and R is `(theta-1)(2theta³+2theta²-2theta-1)=0`. Thus the
coefficient of `xi^4` vanishes. Suppose the result were `c xi^h` for one
`h in {0,1,2,3}`. For h=0 or 1, cross multiplication gives a nonzero
polynomial of degree four below the minimal degree five. For h=2,
comparing coefficients gives `c theta=1` and `c theta=theta`, forcing
theta=1. For h=3, the degree-five difference would be `c theta` times
the minimal polynomial (4), but its constant coefficient is zero whereas
that multiple has constant coefficient `-c theta`. All are impossible.
G3 excludes a nonmonomial of width at most two, so the width is exactly
three. The degree statements follow because the irreducible cubic cannot
have a root in either degree-eight K or degree-forty L. QED.

An [independent exact calculation](../../04-computation/experiments/elliptic_cubic_direction_20261016c.py)
constructs (8) by a 30-by-30 rational multiplication matrix in
`Q(theta,j,xi)` and checks its norm by conjugation. Its
[receipt](elliptic_cubic_direction_20261016c.json) gives all five
coefficients: the coefficient of `xi^4` is zero, those of 1 and `xi³`
are nonzero, and their common rational denominator is 1541. This
independent path verifies the positive witness rather than merely
reusing the vanished Cramer numerator.

This changes the ambient field explicitly; it is not a construction
inside the requested L. However it preserves a surprisingly strong part
of the old interface. Any new common unit neighbour of two old K-ports
would have degree at most two over K, while `[L':K]=15` is odd.
It must therefore already lie in K. Consequently the **same nine exterior
double-port points and no triple-port point** remain the complete census.
Likewise L' has no `sqrt(-3)`: a quadratic extension of L cannot embed in
its degree-three extension. Its unit graph is still triangle-free.

This gives a precise optional reroute: enlarge by the particular cubic
that obstructed the real seed, gain a non-Cartesian direction, and keep
the original port coordinates and contact census. The source boundary
predicate has not been altered. Four-chromaticity or six-colour forcing
for this enlarged carrier still requires an actual graph certificate;
the explicit new unit by itself supplies neither. The fixed-L route
continues to require a different seed or genuinely complex coefficients.

## 6 The huge positive triple and the actual Collatz arrow

The user recalls **one positive integer triple**, with coordinate lengths
79, 80, and 81. The primary
[Bremner--Macleod paper, Remark 2.2](https://ami.uni-eszterhazy.hu/uploads/papers/finalpdf/AMI_43_from29to41.pdf)
identifies it with `9G` on

\[
 E:y^2=x^3+109x^2+224x,\qquad G=(-4,28).
\]

The [2025 corrigendum](https://publikacio.uni-eszterhazy.hu/8862/1/AMI_62_from26to27.pdf)
retains those digit lengths. Both primary papers were reopened for this
session. The repo already has the stronger
[coordinate-controlled decoder](creation_elliptic_20260925.md), including
the torsion phase, the intrinsic rational-half test, and a terminating
binary creation certificate. The exact triple is regenerated and its
rational sum verified again in the receipt.

That certificate is `9G -> 4G -> 2G -> G -> O`. Actual plus-Collatz sends
`9G -> 14G`. Both operations can be performed directly on curve coordinates,
but they certify different predicates. The positive fruit locus is also
lost under the creation decoder. The large decimal height is not a long
creation word; the short word is not an authenticated Collatz descent.

For a native finite Collatz controller, replace each integer endpoint n
by `nG`. Since G has infinite order, this is an injective relabelling.
Use the inherited square-class test to choose the actual arrow

\[
 n\text{ even}:P\mapsto\operatorname{selectedHalf}(P),\qquad
 n\text{ odd}:P\mapsto\operatorname{selectedHalf}(3P+G),
\]

with ROOT `G` stopped. Keep the original source, its coefficient
`w(n)=3[n=2 mod 3]-1`, and the **integer** clock as additional coordinates.
Then every signed terminal row and every Gram expression in
[refinement debt](collatz_refinement_debt_20261016a.md) is unchanged by this
relabelling. Elliptic addition of different endpoints is not a native
merge: `nG+mG=(n+m)G` supplies no common-future receipt.

A useful richer representation of an authenticated source row is

\[
 \mathcal V_I(z)=\sum_{n\in I}w(n)z^{d(n)}[e(n)G],   \tag{7}
\]

in the free module on the retained elliptic endpoints, with polynomial
clock variable z. A native extension of terminal e by exactly tau steps
acts by `[eG] -> z^tau [T^tau(e)G]`, uniformly on all its occurrences.
Composition adds the actual clocks; evaluation at `z=-1` recovers the
existing signed endpoint row exactly. No source coefficients are replaced.

The two valid singleton rows for sources 4 and 1 are `-z²[G]` and `-[G]`.
They agree at `z=-1` but differ before specialization. Their difference
`1-z²` also vanishes at `z=1`; even both parity specializations lose the
integer clock. This is not an objection to the parity-energy observable,
which intentionally uses that specialization. The full polynomial is the
extra certificate data needed when chronology is part of the claim.

The tempting next step of correcting elliptic height by a bounded state
has already been treated more strongly by
[two-place Lyapunov ranks](procgen_rank_20260926_two_place_lyapunov.md),
Theorem A/Corollary A1 and the stated adaptive-counter obstructions.
The elementary Mersenne control explains its relevance: for `k>=2`, the
first k shortcut steps from `2^k-1` are all odd and end at `3^k-1`, producing
unbounded expansion ratio. Since canonical height on `nG` is quadratic
in n, no globally comparable positive height can be nonincreasing at
every native step. Recovering this older theorem prevents re-proving a
weaker finite-state version. Uniform paired-reset weighted cancellation
is still the separate sufficient target; (7) does not prove it.

## 7 Which next constructions escape these proved restrictions

| Source and target | Map | Preserved | Lost unless retained | Decisive next test |
|---|---|---|---|---|
| Contact quartic to unit graph | `(r,y)` to the marked point b and rotation d | Field and unit incidences | Which circles can share a new edge | Full pair distances, then the boundary relation |
| Paley periods to eleventh-root coordinates | Quadratic-residue orbit to (4) | Exact minimal polynomial and degree | Vertex colours and marked ports | A unit of width at least three plus a genuine coupled graph |
| Product seed to field geometry | (3) | Every unit edge at width two | No interaction is generated | Leave the Cartesian region before increasing its size |
| K-piece to reflected translated copy | (5) | Geometry, two named attachments | Relations between several copies | Vary anchors and conjugate parameters together; a shared parameter specializes back into K |
| Integer controller to elliptic controller | `n -> nG`, with (7) | Source coefficients, native endpoints, full clocks | Descent if the wrong halving rule is used | The original weighted overlap bound on the actual transported rows |

The next useful geometric target is therefore **interference between
multiple linkages**, or a genuinely non-Cartesian direction with a coupled
port relation. Increasing the two-family elliptic tree or a width-two
Cartesian grid cannot provide it. A width-three norm-one element inside L
remains a concrete arithmetic helper question. G7 constructs one after
an explicit odd cubic extension while preserving the old port census;
compatibility and four-chromaticity still need separate certificates.
The existing width-four direction shows that L itself contains mixed
directions.

The next useful Collatz target is to use the full marked module to verify
a changing-terminal controller while estimating its actual opposite-sign
overlap. Coordinate compression is now exact; the all-scale quantitative
estimate remains unpaid. The geometric obstruction and the Collatz
obstruction both concern compatibility after composition, but no numerical
bound or convergence claim is transferred between them.

META-PATTERNS used: recover stronger prior results; retain the product and
shared labels before scalar totals; audit a representation under its next
operation. No new universal meta-pattern is claimed from these examples.

## Reproduction

Run from the repository root:

```bash
python3 04-computation/experiments/elliptic_crossedges_20261016c.py
python3 -O 04-computation/experiments/elliptic_crossedges_20261016c.py
python3 04-computation/experiments/elliptic_cubic_direction_20261016c.py
python3 -O 04-computation/experiments/elliptic_cubic_direction_20261016c.py
```

The [script](../../04-computation/experiments/elliptic_crossedges_20261016c.py)
uses exact rational arithmetic, the existing cyclotomic engine, and the
existing elliptic decoder. The [receipt](elliptic_crossedges_20261016c.json)
states the finite universes and supplies complete geometry and colour
witnesses. The finite tests do not replace the all-index proofs of G1--G7.
Ordinary and optimized outputs agree byte for byte for both programs.
The main receipt records 257 always-active checks, including complete
geometry and explicit colouring witnesses. Repository documentation and
relative-link checks also pass.
