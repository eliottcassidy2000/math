# Odd-square coordinates, boundary families, and an inradius permutation

**PROVED**, by the elementary arguments below; **FINITE-EXACT** for the
declared companion computation. This is a recovery and extension of existing
repo mechanisms, with no literature-priority claim. In particular, the
odd-square chart and its Berggren descent are inherited results, not new
discoveries. No Collatz convergence claim is made.

## 1. Inheritance and the objects being compared

The closest proved mechanism is
[THM-3756, odd-square ordinal Berggren affine descent](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md).
It already establishes the lossless chart, totient fibres, three parent cones,
and the nonprimitive content forest. Related mechanisms are
[THM-3334, Berggren parabolic spine and Gaussian collision torsor](../../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md)
and [THM-4057, Stern--Brocot depth pullback](../../01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md).
The [prior semicircle note](collatz_mod6_20260917_pythagorean_semicircle.md)
also gives the chart and Gaussian squaring, while the
[prior geometry audit](arithmetic_braids_20260917_geometry.md) corrects the
integer-versus-rational parameter distinction and proves the sole integer
parameter duplication at 2 and 3.

The canonical hostile example is the collision of the two triples
`(5,12,13)` and `(15,8,17)` at the same upper odd square `25`. The corrected
near miss is that integer parameters in `(k²−1,2k,k²+1)` do not give all
primitive shapes; rational parameters do. The less-used sidecar here is
the fixed-square fibre as a unit group modulo sign, rather than as a level
of an expanding tree. This is distinct from the state-dependent Pell return
in [THM-3819, Chebyshev--Pell Berggren return conductor staircase](../../01-canon/theorems/THM-3819-chebyshev-pell-berggren-return-conductor-staircase.md).

The live concept board is: primitive triangles; odd-square roots; Berggren
descent; integer/rational shape parameters; unit doubling and cyclotomic
traces; Collatz edges. The anchor is lossless indexing, the niche is a
permutation within one fibre, and the wildcard is its exact Chebyshev
realization. Sections 4 and 7 state the preserved and destroyed information.

Throughout, `(a,b,c)` means positive integral sides with `a` odd, `b` even,
`a²+b²=c²`, and gcd one, unless explicitly stated otherwise.

## 2. What the odd square does and does not index

**OS-T1 (inherited chart).** Every such triple has a unique representation

\[
 a=uv,\qquad b=(u^2-v^2)/2,\qquad c=(u^2+v^2)/2,
 \qquad u>v>0\text{ odd},\quad\gcd(u,v)=1.                 \tag{1}
\]

The inverse is `u=√(c+b)`, `v=√(c−b)`. Indeed, the two positive odd factors
`c+b` and `c−b` are coprime and have product `a²`, so each is a square.
Conversely (1) gives integral sides of the required parity. An odd prime
dividing all sides would divide both roots; a factor 2 is excluded by `a`
being odd. Thus coprime roots are sufficient as well as necessary.

The Euclidean parameters are `m=(u+v)/2`, `n=(u−v)/2`. They are coprime and
of opposite parity. Notice the typed distinction

\[
c+b=u^2,\quad c-b=v^2,\qquad c+a=2m^2,\quad c-a=2n^2.
\]

It is the **even** leg that produces odd squares. The other two sums are not
positive integer squares. Keeping both roots retains the complete ordered
triple; keeping only `u` loses a finite fibre.

**OS-T2 (inherited fibre).** For every odd `u≥3`, that fibre has `φ(u)/2`
members. The units `1,…,u−1` pair as `v,u−v`; exactly one member of each pair
is odd. The fibre is therefore the unit quotient
`(Z/uZ)×/{±1}`, represented by its unique odd integer in `(0,u)`.
For `u=1` it is empty. The first ambiguity is `u=5`, with `v=1,3` giving
the two triples above. An index by odd squares is a grouping, not a unique
natural-number address, unless the inner coordinate is included.

Without coprimality, put `g=gcd(u,v)`. The coordinate gcd in (1) is exactly
`g²`, since division by `g²` leaves a primitive triple. Thus the same chart
contains odd-square multiples, such as `(u,v)=(9,3)` giving
`(27,36,45)=9(3,4,5)`. It does not contain every nonprimitive scale:
`2(3,4,5)` has `c+b=18`.

## 3. A tree of all triples, already proved in the repo

**OS-T3 (inherited Berggren root chart).** The three children are

\[
 L(u,v)=(u+2v,v),\quad M(u,v)=(2u+v,u),\quad R(u,v)=(2u-v,u).
                                                               \tag{2}
\]

Their matrices have determinant `1`, `−1`, and `1` respectively, so they
preserve the root gcd; parity and the positive chamber are immediate.
The unique parent, except at a root, is

| Child chamber | Parent | Branch |
|---|---|---|
| `u>3v` | `(u−2v,v)` | L |
| `2v<u<3v` | `(v,u−2v)` | M |
| `v<u<2v` | `(v,2v−u)` | R |

The boundary `u=2v` is impossible for odd roots. The boundary `u=3v` is
terminal: primitively it is `(3,1)`, and generally it is `(3g,g)`.
Each parent strictly decreases `u`. Hence descent terminates, gives a
unique branch word, and proves that (2) generates every primitive triple
exactly once. No statistical or Collatz descent assumption is used.

The depth at a fixed `u` is at most `(u−3)/2`, with equality precisely on
`v=1` or `u−v=2`; these are the pure L and pure R rays. This equality
classification, and conversion to the odd ordinals `(u+1)/2,(v+1)/2`, are
proved in THM-3756 cited above. The integer family below is exactly these
two extreme rays.

## 4. A permutation that replaces the even leg by four inradii

**OS-T4 (fixed-fibre operation).** For fixed odd `u≥3`, define

\[
 D_u(v)=|u-2v|.                                                \tag{3}
\]

It preserves `0<v<u`, oddness and `gcd(u,v)`, and permutes the admissible
fibre. Given an odd output `w`, the two formal preimages are
`(u−w)/2` and `(u+w)/2`. They sum to the odd integer `u`, so exactly one is
odd; that is the unique admissible inverse. Both lie in `(0,u)` and have
the same gcd with `u` as `w`. The real folded map `t↦|1−2t|` has two
branches, but the retained parity coordinate chooses one inverse on this
finite arithmetic set.

In sides the same operation is

\[
\mathcal D(a,b,c)=
\bigl(|b+c-2a|,\ 2(a+b-c),\ 3c-b-2a\bigr).                  \tag{4}
\]

This follows by substituting (1) and replacing `v` by `|u−2v|`.
All sides remain positive, the triple remains primitive, and `b+c=u²`
stays fixed. If `r_in=(a+b−c)/2` is the inradius, then the new even leg is
exactly `4r_in`. This is an integral operation inside one fixed-square
family, although the hypotenuse usually changes.

There is also an exact quadratic-form certificate: before taking the
absolute value of the first coordinate, (4) has matrix

\[
 Q=\begin{pmatrix}-2&1&1\\2&2&-2\\-2&-1&3\end{pmatrix},
 \qquad Q^{\mathsf T}\operatorname{diag}(1,1,-1)Q
       =4\operatorname{diag}(1,1,-1).
\]

**OS-T5 (complete cycle law).** Under the quotient identification in
OS-T2, (3) is exactly multiplication by 2. Consequently all its cycles
have the same length

\[
d_u=\min\{k\ge1:2^k\equiv 1\text{ or }-1\pmod u\},
\qquad \#\text{cycles}=\frac{\varphi(u)}{2d_u}.                \tag{5}
\]

To see the uniformity, a unit class `[v]` returns after `k` steps exactly
when `[2^kv]=[v]`; cancellation of `v` makes this independent of the point.
Dropping coprimality invalidates that uniform period statement: `(9,3)`
is a fixed point, whereas the primitive `u=9` fibre is a three-cycle.

**OS-T6 (an actual Chebyshev conjugacy on a finite set).** Let `ζ_u` be a
primitive `u`th root of unity. The map

\[
 (u,v)\longmapsto \zeta_u^v+\zeta_u^{-v}
                 =2\cos(2\pi v/u)                           \tag{6}
\]

is a bijection from the primitive fixed-`u` fibre to these real cyclotomic
conjugates, and intertwines `D_u` with `x↦x²−2`. Equality of traces is
equivalent to equality of root-of-unity exponents modulo sign; the
intertwining identity is `(z+z⁻¹)²−2=z²+z⁻²`.

This retains the entire finite dynamical system. It does **not** identify
the cyclotomic angle `2πv/u` with the triangle angle
`ψ=2 arctan(v/u)`, opposite the odd leg. Gaussian squaring doubles the
oriented complex argument and squares the hypotenuse; converting back
to positive ordered legs folds the angle. For this acute odd-leg angle,
the resulting value is `|π/2−2ψ|`. Operation (4) instead preserves
`b+c` and is periodic. Berggren generation increases `u`. These are three
different operations on triangles, not three descriptions of one map.

Examples, always in odd-leg/even-leg order:

| u | v cycle | Corresponding triangle cycle |
|---|---|---|
| 3 | 1 | `(3,4,5)` |
| 5 | 1,3 | `(5,12,13)`, `(15,8,17)` |
| 7 | 1,5,3 | `(7,24,25)`, `(35,12,37)`, `(21,20,29)` |
| 9 | 1,7,5 | `(9,40,41)`, `(63,16,65)`, `(45,28,53)` |

The exact period-three conductors are `u=7,9`: divide `2³−1` or `2³+1`
and exclude shorter returns. For period six the same test on `63,65`
gives exactly `u=13,21,63,65`, with respectively `1,1,3,4` cycles, hence
nine six-cycles on 54 triangles. These conclusions hold for all odd `u`,
not merely the finite test cap, because any such conductor divides
`2^k−1` or `2^k+1`.

**OS-T7 (power-of-two neighbours).** For every `m≥1`,
`d_(2^m+1)=m`. The return at `m` is negative, and for `j<m`, both positive
integers `2^j−1` and `2^j+1` are smaller than `2^m+1`. Likewise
`d_(2^m−1)=m` for `m≥3`; the size comparison holds again. The excluded
case `m=2` has `u=3,d_u=1`. Thus `63` and `65` furnish respectively three
and four six-cycles for an explicit modular reason. This is a genuine
role for `63=2^6−1`; it does not explain its separate exceptional status
in a primitive-prime-divisor theorem.

If `u=2^(2^r)+1` is a Fermat prime, (5) gives `d_u=2^r` and
`2^(2^r−r−1)` cycles. The examples `u=3,5,17,257,65537` therefore have
`1,1,2,16,2048` cycles. Primality is needed only to replace `φ(u)` by
`u−1`; the signed order formula itself is valid for composite neighbours
as well. This makes the first split at 17 precise without invoking an
unrelated periodicity theory.

## 5. Integer parameters are two rays; rational parameters cover all shapes

**OS-T8 (exact parameter classification).** Reduce
`(k²−1,2k,k²+1)` to its primitive triple and place its odd leg first.
For integer `k≥2`, the gcd is one when `k` is even and two when `k` is odd.
Its roots are respectively

\[
 (u,v)=(k+1,k-1)\quad(k\text{ even}),\qquad
 (u,v)=(k,1)\quad(k\text{ odd}).                              \tag{7}
\]

The two boundary rays intersect only at `(3,1)`, corresponding to
`k=2,3`. The first primitive triple omitted by these integer parameters
is `(21,20,29)`, with `(u,v)=(7,3)`. Thus a family whose raw hypotenuse is
one more than an integer square is not a full primitive-triple chart.

For a reduced rational `k=p/q>1`, clear denominators to obtain
`(p²−q²,2pq,p²+q²)`. Its gcd is one when `p,q` have opposite parity, and
two when both are odd. After reduction and leg ordering, every primitive
shape occurs at exactly two rational parameters:

\[
 k=\frac{u+v}{u-v},\qquad k^*=\frac uv,
 \qquad k^*=\frac{k+1}{k-1}.                                 \tag{8}
\]

They are the two choices of which leg to place first in the usual rational
circle chart. The involution has no rational fixed point greater than one
(its positive fixed point is `1+√2`). If the raw first leg is required to
be odd after primitive reduction, the opposite-parity Euclidean chart
chooses just the first sheet. For the omitted integer example the two
rational parameters are `5/2` and `7/3`.

**OS-T9 (exact sparse count).** Let `N_int(X)` count distinct primitive
triples of this integer family with `c≤X`, for integer `X`. It is zero
for `X<5`, and otherwise

\[
 N_{\rm int}(X)=
 \left\lfloor\frac{\lfloor\sqrt{X-1}\rfloor}{2}\right\rfloor
 +\left\lfloor\frac{\lfloor\sqrt{2X-1}\rfloor-1}{2}\right\rfloor-1
 =\frac{1+\sqrt2}{2}\sqrt X+O(1).                            \tag{9}
\]

The first term counts the even `k≥2`, the second the odd `k≥3`, and the
subtraction removes the repeated `(3,4,5)`. The companion
[shape-law note](odd_square_20260921_shapes.md) proves the full count
`N_PPT(X)=X/(2π)+O(√X log X)`. Hence the relative share of this integer
family is asymptotic to `π(1+√2)/√X`, and tends to zero.

## 6. Closing the two rays under the permutation still misses triples

**OS-T10.** For each fixed `u`, the two boundary nodes `v=1,u−2` lie in
the same `D_u` cycle, since `D_u(1)=u−2`. Therefore closing the entire
integer family under (4) fills one cycle per fibre, of size `d_u`.
It fills a fibre if and only if `d_u=φ(u)/2`.

The first failure is `u=17`, where the cycles are

\[
 (1,15,13,9),\qquad (3,11,5,7).
\]

The triple `(51,140,149)` at `(17,3)` is the first missing one by
hypotenuse after that closure. For `u<17` all fibres are covered (the
finite list is `3,5,7,9,11,13,15`); for `u≥19`, even the least hypotenuse
is `(u²+1)/2≥181`. Within `u=17`, the first missing `v` is 3.

The full quotient group of units acts transitively on a fibre; the
single multiplier 2 need not generate it. Whole-fibre distribution
results in the shape-law note consequently do not imply an individual
`D_u` orbit has the same shape distribution. This is the precise lost
coordinate in a claim that one doubling orbit generates every triangle.

## 7. The Collatz edge test fails at a small explicit example

For the accelerated signed family
`U_b(n)=(3n+b)/2^v2(3n+b)`, the positive `b=−1` map has `5↔7`.
Taking the unordered odd pair `(u,v)=(7,5)` gives `(35,12,37)`, but (4)
sends it to `(21,20,29)`, whose root pair is `(7,3)`.
Neither direction between 7 and 3 is an edge for `b=−1` or `b=+1`:

\[
 U_{-1}(7)=5,\quad U_{-1}(3)=1,\qquad
 U_{+1}(7)=11,\quad U_{+1}(3)=5.
\]

Thus (4) preserves primitive triangles and the odd square, but destroys
the Collatz edge predicate. An edge label and a direction are necessary
sidecars if the source is a Collatz graph. The period-three or period-six
fibres above cannot be identified with Collatz components merely from
their sizes. The successful bridge is (6), where the map, inverse,
preserved dynamics and information loss have all been specified.

## 8. Exact verification and scope

Run from the repository root:

```text
python 04-computation/experiments/odd_square_20260921_triangles.py
python -O 04-computation/experiments/odd_square_20260921_triangles.py
```

The [standard-library script](../../04-computation/experiments/odd_square_20260921_triangles.py)
writes a [JSON certificate](../../04-computation/experiments/odd_square_20260921_triangles.json)
with an LF-normalized source hash, explicit universes, every fixed-fibre
cycle through `u=301`, and active exception-based checks that survive `-O`.
It verifies 9,242 primitive root pairs, all nonprimitive odd pairs in that
range, the independent Berggren traversal, and the rational two-sheet
reconstructions. An independent Euclidean-parameter census agrees with
the root census on 1,593 triples through `c=10000`; a direct side-length
census agrees on 158 triples through `c=1000`. It checks integer
parameters `2≤k≤2000`, exact count cutoffs through `10^6`, power neighbours
`2^m±1` through `m=16`, and the positive and hostile controls above.

The Chebyshev check uses exact Laurent polynomials modulo `z^u−1`, not
floating-point agreement of cosines. The uniform cycle law, global chart,
descent, counting formula and conductor statements are proved above;
finite checks do not stand in for those proofs. No assertion about
Collatz convergence, density along Collatz orbits, rational-quadratic
global conjugacy, or single-cycle equidistribution follows from this note.
