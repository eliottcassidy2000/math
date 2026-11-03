# Nine-cycle carry and complete degree-three magic relabeling classification

**Status:** PROVED (self-contained arguments below), with FINITE-EXACT controls.
The unrestricted problem of a 3-by-3 magic square of nine distinct integer
squares remains **OPEN**. No originality or global-conjecture closure is claimed.

## 0. Complete classification of degree-at-most-two numerical assignments

**Theorem.** For real `A,B,C`, let `f(k)=Ak^2+Bk+C`, `k=-4,...,4`. Its nine
indexed values can be arranged in a 3-by-3 magic square if and only if `A=0`.
Repetitions and negative entries are allowed. If distinct entries are required,
the complete condition is `A=0` and `B!=0`.

**Proof.** Every real 3-by-3 magic square is centrally symmetric as a multiset
about its center/mean `mu`: its four opposite pairs sum to `2mu`. Consequently
all odd centered moments vanish, in particular `sum((x-mu)^3)=0`.

For the values of `f`, the mean is `mu=C+20A/3`. Using

    sum(k^2)=60, sum(k^4)=708, sum(k^6)=9780,

and the vanishing odd power sums, exact expansion gives

    sum((f(k)-mu)^3) = (44/3) A (65A^2+63B^2).

The mixed `A^2 B` and `B^3` terms vanish. The `A^3` coefficient is
`9780-20*708+3*(400/9)*60-9*(8000/27)=2860/3`; the `A B^2` coefficient is
`3*(708-(20/3)*60)=924`. The expression is nonzero whenever `A!=0`, so such
an arrangement is impossible. If `A=0`, arrange the indices by the centered
Lo Shu grid

     3 -4  1
    -2  0  2
    -1  4 -3.

Every line sums to zero; applying `k->Bk+C` gives eight line totals `3C`.
These values are distinct exactly when `B!=0`. QED.

**Square-root corollary.** With `f(k)=(a+dk)^2`, the cubic centered moment is

    (44/3) d^4 (65d^2+252a^2) > 0   when d!=0.

Thus no such squared-root arithmetic progression is magic. Section1 supplies
an independent center/maximum proof not using the moment expansion.

### 0.1. General centrally symmetric packets and the sharp small boundary

For `k=-m,...,m`, put `n=2m+1` and `q=m(m+1)`. For the same quadratic `f`,
its mean is `C+Aq/3`, and the elementary power-sum identities

    average(k^2)=q/3,
    average(k^4)=q(3q-1)/15,
    average(k^6)=q(3q^2-3q+1)/21

give

    sum((f(k)-mean)^3)
      = n A q(4q-3)/945 * ((4q-15)A^2+63B^2).

Thus a genuine real quadratic sampled on any odd number of equally spaced
arguments at least five cannot be centrally symmetric as a multiset. Nonzero
argument spacing is understood and can be absorbed into the coefficients.
This does not exclude arbitrary higher-order magic squares, which need not
be centrally symmetric. For three arguments the boundary is sharp: with
`A!=0`, central symmetry holds exactly when `B=+A/3` or `B=-A/3`; for example
`3k^2+k` gives the three values `2,0,4`. The real-coefficient hypothesis is
essential to the positivity argument.

### 0.2. Complete odd-cubic classification: a positive degree-three boundary

**Theorem (finite-exact linear classification).** Let `a,b,C` be real, with
`a!=0`. The nine values

    f(k)=a k^3+b k+C,  k=-4,...,4,

can be arranged in a 3-by-3 magic square if and only if `b=-15a`. All nine
values are distinct in that case. This classifies odd cubics plus a constant,
not arbitrary cubics with a quadratic term, and not magic squares of squares.

**Reduction and finite proof.** Applying the same nonzero affine transformation
to all entries preserves magic. Subtract `C` and divide by `a`, writing
`t=b/a` and `p_i=i^3+t i` for `i=1,2,3,4`. The mean is zero. Every centered
3-by-3 magic square has multiset

    {0, +/-u, +/-v, +/-(u+v), +/-(u-v)},

as seen by solving the eight line equations (or by the displayed general
magic-square formula

    [u,    -u-v, v]
    [-u+v,  0,   u-v]
    [-v,    u+v, -u]).

Therefore a solution exists
exactly when a signed permutation of the four indexed `p_i` equals the
ordered quartet `(u,v,u+v,u-v)`. Repeated values cause no exception: within
each absolute-value class, match the indexed opposite pairs arbitrarily;
if a pair is zero, its two copies can likewise be paired.

For every permutation `(i,j,h,l)` of `(1,2,3,4)` and sign tuple
`(s_i,s_j,s_h,s_l)` in `{+1,-1}^4`, impose

    (s_i i+s_j j-s_h h)t + s_i i^3+s_j j^3-s_h h^3 = 0,
    (s_i i-s_j j-s_l l)t + s_i i^3-s_j j^3-s_l l^3 = 0.

These are two explicit linear equations, not a bounded search over `t`.
The checker enumerates all `24*16=384` branches using exact rational arithmetic.
No branch has both equations identically zero; precisely eight branches
survive, all with `t=-15`. Hence the classification covers every real `t`,
not merely integer or bounded parameters. The output records all survivors.
The finite reduction is independently audited; the normal and optimized
checker outputs byte-match.

Conversely the root's positive witness is

     18 -22   4
    -14   0  14
     -4  22 -18,

whose multiset is exactly `{k^3-15k:-4<=k<=4}` and whose eight totals are zero.
Scaling by `a` and adding `C` gives the entire family. This also proves the
degree-two obstruction must not be extrapolated to all nonlinear assignments.
Making these nine numbers positive by adding a constant is easy; making all
nine distinct integer squares is a separate, unsolved condition in this packet.

### 0.3. Full real cubic classification

**Theorem (finite-exact linear classification, independently audited).** For
real coefficients, the nine values

    f(k)=A k^3+B k^2+Ck+D,  k=-4,...,4,

can be arranged into a 3-by-3 magic square **if and only if**

    B=0 and (A=0 or C=-15A).

Thus, apart from linear/constant assignments, the single nonlinear family is
`A(k^3-15k)+D`. Every member with `A!=0` has distinct values. `D` is arbitrary.

**Completeness reduction.** The mean is `D+20B/3`. Pick the input index `j`
whose value occupies the center. The remaining eight indices must be paired
into four opposite pairs. There are exactly nine center choices and
`7*5*3*1=105` pairings per center, with no distinctness assumption. The center
equation and each pair `(u,v)` give homogeneous linear equations in `(A,B,C)`:

    3j^3 A + (3j^2-20)B + 3j C = 0,
    (u^3+v^3-2j^3)A + (u^2+v^2-2j^2)B + (u+v-2j)C = 0.

The constant `D` cancels. Exhausting these **945 exact five-by-three matrices**
classifies every real coefficient vector compatible with opposite pairs;
this is not a bounded coefficient search. Rational Gaussian elimination and
an independent rank-by-minors path give the rank census

    rank3: 914; rank2: 30; rank1: 1.

The sole two-dimensional kernel is the plane `B=0`. Every one-dimensional
kernel with `B=0` lies in that plane. The only kernels outside it are the four
coefficient rays

    (A,B,C)=(1,+/-3,-18), (2,+/-3,-39).

Their positive centered deviation quartets, up to common nonzero scaling and
sign, are respectively

    (20,20,34,36), (20,54,56,70).

Neither can be a magic quartet. For a real quartet
`{|u|,|v|,|u+v|,|u-v|}`, the largest member is `|u|+|v|`, the sum of two
members, and hence is at least the sum of the two smallest members. Here
`36<20+20` and `70<20+54`. Thus all nonzero-`B` kernels fail the required
magic arrangement. The remaining plane `B=0` is exactly Section0.2's
odd-cubic-plus-linear family (with linear assignments handled separately).
The converse is supplied by its displayed magic grids. QED.

**Independent paths.** `full_cubic_audit.py` enumerates all 945 systems and
checks every rank by both Gaussian elimination and minors. The independent
`audit/full_cubic_nonodd.py` normalizes `B=1`, eliminates `C` using the center
equation, and solves candidate pair-edge equations followed by exact perfect
matching. It finds precisely the same four nonodd rays and rejects all four.
A separate root audit also independently replayed the 945-system route
without importing producer code. Normal and optimized outputs agree.

Do not confuse the center's **value** with its input index: for general
odd cubics, a nonzero index can also evaluate to the mean. The 945-system
enumeration covers all index choices, including repeated zero values.

## 1. All signed arithmetic-progression roots are excluded

**Theorem.** Let `a,d` be real numbers with `d != 0`. The indexed multiset

    {(a+kd)^2 : -4 <= k <= 4}

cannot be arranged into a 3-by-3 magic square. The conclusion does not require
the nine squared values to be distinct, integral, rational, or strictly positive.

Here magic means that all three row sums, all three column sums, and both main
diagonal sums are one common number `T`.

**Proof.** In any 3-by-3 magic square with center `c`, adding the central row,
central column, and the two diagonals counts every noncentral entry once and
the center four times. Since the total of all entries is `3T`, this gives
`4T=3T+3c`, hence `T=3c`. Thus the center equals the mean of all nine entries,
and every opposite pair sums to `2c`.

Changing the sign of `d` or of `a` only permutes the squared multiset, so assume
`d>0` and `a>=0`. Its mean is

    mu = a^2 + (20/3)d^2,

because `sum(k)=0` and `sum(k^2)=60` on `k=-4,...,4`.

* If `a>=4d`, all roots are nonnegative and increasing. The two entries `a^2`
  and `(a+d)^2` are adjacent in their sorted squared list, whereas

      a^2 < mu < (a+d)^2,

  since `2ad+d^2 >= 9d^2 > (20/3)d^2`. The required center is not an entry.
* If `0<=a<4d`, the largest squared entry is `M=(a+4d)^2`, and

      M-2mu = a(8d-a) + (8/3)d^2 > 0.

  All entries are nonnegative. Since `mu>0`, this maximum is not the center;
  its required opposite entry would be `2mu-M<0`, impossible.

These alternatives cover all parameters. QED.

**Boundary.** At `d=0`, the constant square-entry grid is magic, so the excluded
parameter is essential. For signed roots, repeated squares do not provide an
escape. The proof is all-real; no finite numerical search is being extrapolated.

**Independent finite reduction.** If the center has root index `j`, then `j=0`
is impossible, and for `j in {-4,-3,-2,-1,1,2,3,4}` the center equation forces

    a/d = (20-3j^2)/(6j).

The checker exhausts these eight exact rational candidates. Each has maximum
greater than twice the mean; moreover no pair of distinct noncentral root
indices has squared sum twice the mean. This is a second, finite algebraic
route, not the source of the two-case real proof.

## 2. What the cycle actually is

The displayed word

    4,8,3,7,2,6,10,5,9,4

has **nine distinct states and nine transitions**; the final `4` is the return.
On labels `2,...,10`, its transition is

    T(l) = 2 + ((l-2+4) mod 9).

It is the translation by `4` in `Z/9Z`, and has order nine because `gcd(4,9)=1`.
In ordinary integer labels it has five jumps of `+4` and four jumps of `-5`.
It is not an eight-state object.

Write `l-2=u+3v` with `u,v in {0,1,2}`, then regard the digits as elements of
`F3`. Its faithful digit formula is

    T(u,v) = (u+1, v+1+1_{u=2})
           = (u+1, v+1+u+2u^2),

with arithmetic in `F3`. After three transitions the low digit returns and
exactly one carry has occurred:

    T^3(u,v) = (u,v+1).

After nine transitions both digits return. Deleting the carry gives
`(u,v)->(u+1,v+1)`, which has order three and loses the actual cycle.

**No affine repair by relabeling.** There is no bijection of these nine states
with `F3^2` that makes this permutation affine. Indeed `AGL_2(F3)` embeds into
`GL_3(F3)` using homogeneous coordinates. A matrix with order nine would have
minimal polynomial dividing `X^9-1=(X-1)^9`, so it is unipotent. A unipotent
3-by-3 matrix is `I+N` with `N^3=0`; in characteristic three,
`(I+N)^3=I+N^3=I`, a contradiction. Permutation order is preserved by relabeling.
The displayed degree-two formula is therefore minimal among polynomial-degree
realizations in these two digit coordinates, and no choice of coordinates
makes degree one sufficient.

**Cube boundary.** The cycle graph is `C9`. Every hypercube graph is bipartite,
so there is no edge-preserving graph map of this odd cycle into a hypercube
(in particular no injective adjacency-preserving embedding). This does not
forbid drawing the states near a cube or embedding the nine points as a set;
such drawings do not preserve its transition graph.

The carry mechanism is a direct `p=3`, two-digit specialization of the
repository's previously proved odometer mechanism (THM-2788, exact slug in
SOURCE_PINS.json); it is not a new general extension theory.

## 3. The precise bridge to magic squares, and its limit

The cycle's labels are `l=6+k`. Section0 completely classifies **all real
polynomial numerical assignments of degree at most three** to the centered
labels: constants and linear assignments admit magic arrangements; every
genuine quadratic fails; the only genuine cubic family is
`A(k^3-15k)+D`. This is stronger than excluding affine root assignments
followed by squaring. It does not classify degree-four or higher assignments.

**Corollary.** Assign a real root `r(l)=alpha*l+beta`, `alpha!=0`, to each of
the nine distinct labels in this cycle. No permutation of the nine squared
roots `r(l)^2` is a 3-by-3 magic square.

The roots are an arithmetic progression; apply Section1 with
`a=6alpha+beta` and `d=alpha`. This proves an obstruction independent of the
cycle's starting point, direction, chosen square-grid placement, translation,
or nonzero scaling of the roots. It does **not** cover a nonlinear root
assignment, arbitrary nine integer squares, or decorated-biclique rational
points.

For the literal roots `2,...,10`, the squared sum is `384`, so the required
center is `128/3`, not even integral. The unsquared labels do admit the shifted
Lo Shu arrangement

    9  2  7
    4  6  8
    5 10  3

with all eight totals `18`. Squaring its entries destroys that property.

**Connection contract.** Source: the cyclically ordered label set with its
actual wrap carry. Target: an unordered squared-root multiset and a proposed
magic-square placement. Map: apply one common affine root assignment and
square. Preserved: nine indexed labels, root AP structure, sum and sum of
squares independently of placement. Lost: cyclic order and carry history.
Required sidecar: the proposed center and four opposite pairs. Decisive tests:
center equals mean, the third centered moment vanishes, and (for square
entries) opposite pairs are nonnegative complements. Retaining the whole
matching, rather than only its moments, completes the degree-three numerical
classification. This is not a classification of arbitrary integer squares.

## 4. Variant firewall and inherited status

The recovered Session57 work already distinguishes:

1. **Nine distinct positive integer squares, eight equal totals:** OPEN.
2. **Integer-square entries, equal total also an integer square:** impossible
   for a nonzero such grid, since `T=3e^2` has odd 3-adic valuation.
3. **Nine distinct positive square entries, eight square totals allowed to
   differ:** explicitly realized by the inherited root grid

       2655 2976  320
       2624 2520 1665
       1440  895 3624.

   Its six row/column totals are `4001^2`; its diagonals are `5151^2` and
   `2920^2`. It is **not magic**. The checker independently verifies the
   arithmetic with explicit exceptions (not disabled by `python -O`).
4. **Bare square-sum graph/biclique:** positive constructions are known.
   Session57 proves the actual MSOS bridge needs reflected parallelogram
   decoration and nine distinct square edge labels. A bare `K4,4`, an
   original-entry Hamiltonian path, or a cycle alone is not enough.

The prior geometry file identifies the decorated locus with the usual
six-quadric surface in `P8`, distinct from the six-dimensional bare biclique
variety. This new packet does not reprove all of that geometry or claim its
open rational-point problem is solved.

## 5. Current literature check and scope

Checked 2026-11-03: Auel and Singer, *The algebraic geometry of 3-by-3 magic
squares of squares*, [arXiv:2609.09351v1](https://arxiv.org/html/2609.09351v1),
submitted 2026-09-08. The Introduction and Section1 explicitly retain the
distinct rational/integer-square existence question as unknown. Their
surface and automorphism results do not supply a prize solution. This is a
primary-source status check, not an independent audit of that entire paper.
The bounded source search did not find a later resolution; no universal
literature-exhaustion or novelty claim follows from that search.

## 6. Reproduction

    python3 verify.py > results.json
    python3 -O verify.py > results.optimized.json
    cmp results.json results.optimized.json
    python3 full_cubic_audit.py > full_cubic_audit.json
    python3 -O full_cubic_audit.py > full_cubic_audit.optimized.json
    cmp full_cubic_audit.json full_cubic_audit.optimized.json
    python3 audit/full_cubic_nonodd.py > audit/full_cubic_nonodd.json
    python3 -O audit/full_cubic_nonodd.py > audit/full_cubic_nonodd.optimized.json
    cmp audit/full_cubic_nonodd.json audit/full_cubic_nonodd.optimized.json
    python3 verify_manifest.py
    python3 verify_hostile_controls.py

`verify.py` checks the cubic-moment polynomial coefficient by coefficient,
all 384 signed linear branches of the odd-cubic classification,
all 432 affine permutations of `F3^2`, exact cycle/digit identities, all eight
center-equation candidates, and semantic controls.
The affine element orders are `1:1, 2:45, 3:80, 4:54, 6:144, 8:108`; none is9.
No source checkout, network, package installation, or large finite search is
required. Source pins and the portable manifest freeze the published bytes.
The hostile checker changes only temporary copies: one source byte, one cycle
state, and one cubic-witness entry. All twelve normal/optimized rejection runs
must fail for the intended reason while the original bundle stays unchanged.
