---
id: THM-3396
title: "Four-bit pairwise-independent Fourier cone"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  An unbiased pairwise-independent
  law on four Rademacher bits is uniquely determined by its four cubic Walsh
  moments a_i and quartic moment d.  Its sixteen atom masses are
  2^-4(1+q(x)(d+a dot x)); nonnegativity is equivalent to one sharp absolute
  inequality on each parity class.  The integer form classifies every
  OA(N,4,2,2) cell-count packet, and coordinatewise multiplication of
  independent packets multiplies the five surviving moments.  Thus quadratic
  Gram data are completely blind on this cone, while the cubic-quartic packet
  is a complete minimal sidecar.  The moment polytope is the polar of the odd
  5-demicube: it has 26 vertices, namely the ten signed coordinate vectors and
  the sixteen even sign vectors divided by three, and f-vector
  (26,120,160,80,16).  Its sharp absolute functional norm is
  max(l-infinity,l-one/3), and its vertices are exactly ten uniform eight-run
  and sixteen twelve-run OA packets.  Repeated independent coordinatewise
  multiplication either converges exponentially to the uniform law, with an
  exact five-mode chi-square formula, or starts at one of the ten signed
  half-cubes and is fixed or period two.  Three exact Hadamard-puzzle packets of
  sizes 48, 120, and 896 are verified; the 896 packet omits exactly two atoms,
  lies in the relative interior of their triangular-bipyramid ridge, and has
  two different exact H8/H12 vertex decompositions.  This is a
  four-column/local realizability theorem, not a Hadamard completion theorem
  or a Grothendieck bound.
source: root-2608-sign-puzzle-2026-08-14
audit: independent Fourier/OA/convolution and raw-signword reconstruction; exact rational polar face lattice; separate support-function, edge-incidence, and 896-face replay; normal/optimized agreement
related:
  - THM-3394-twelve-formerly-missing-hadamard-orders-through-2000
  - THM-3392-bipartite-sign-lift-and-synchronization-loss
script: 04-computation/four_bit_pairwise_independent_fourier_cone_thm3396.py
output: 05-knowledge/results/four_bit_pairwise_independent_fourier_cone_thm3396.out
script_sha256: 43bc4e9bdc99bf0723ace6d124128a00fd8f89bd43e8738cf7b4bb8ce2139bcc
output_sha256: 47bc7bda869b6b5ce799fee770e8ca56dc1fbbcb1c9a33145678d411f3edfa0b
semantic_sha256: 3b23fb3d84c1d191c975384b3a789779311e19f5136e346743786679d97acc8e
hash_basis: LF-normalized bytes
---

# THM-3396 -- the exact higher-chaos sidecar after four-bit whitening

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. The cone

Let `X=(X_1,X_2,X_3,X_4)` take values in `{+-1}^4`.  Assume that its
coordinates are unbiased and pairwise independent, equivalently

```text
E X_i=0,                 E X_i X_j=0  (i != j).               (1)
```

For `x in {+-1}^4` put `q(x)=x_1x_2x_3x_4`, and define

```text
a_i = E[q(X)X_i] = E product_(j != i) X_j,
d   = E q(X).                                                   (2)
```

Then the law is uniquely determined by `(a_1,a_2,a_3,a_4,d)`, with

```text
P(X=x) = (1/16)[1+q(x)(d+a dot x)].                            (3)
```

Conversely, a real packet `(a,d)` is the packet of such a probability law
if and only if

```text
max_(q(x)=+1) |a dot x| <= 1+d,
max_(q(x)=-1) |a dot x| <= 1-d.                                (4)
```

Thus (4) is an exact five-dimensional polytope, not merely a necessary
moment relaxation.  It also makes the equality boundary literal: a facet is
reached precisely when at least one of the sixteen masses in (3) vanishes.

## 2. Proof by Fourier inversion

For `S subset {1,2,3,4}`, write `chi_S(x)=product_(i in S)x_i`.  The sixteen
characters form an orthogonal basis for functions on the four-cube.  If
`mu_S=E chi_S(X)`, Fourier inversion gives

```text
P(X=x)=2^-4 sum_S mu_S chi_S(x).                               (5)
```

The empty coefficient is one, while (1) kills every coefficient of degree
one or two.  The degree-three character missing coordinate `i` is exactly
`q(x)x_i`, and the only degree-four character is `q(x)`.  Substitution in
(5) proves (3) and uniqueness.

For fixed parity `q(x)=+1`, the two atoms `x,-x` have numerators
`1+d+a dot x` and `1+d-a dot x`.  Their simultaneous nonnegativity is the
first inequality in (4).  On the odd parity class the corresponding
numerators are `1-d-a dot x` and `1-d+a dot x`, giving the second inequality.
This proves necessity.  Conversely, (4) makes every value in (3)
nonnegative; character orthogonality makes their sum one and recovers exactly
the coefficients prescribed above.  This proves sufficiency.

The uniform law is an interior point, so the affine cone has full dimension
five.  Consequently the four cubics and the quartic are not only sufficient:
after the degree-at-most-two data have been fixed, no one of the five can be
deleted from a globally complete linear sidecar.

## 3. Exact polar geometry

Put `y=(a_1,a_2,a_3,a_4,d)` and associate to an atom `x` the sign vector

```text
n(x)=-q(x)(x_1,x_2,x_3,x_4,1).                              (4a)
```

Its five coordinates have product `-1`, and `x -> n(x)` is a bijection from
the four-cube to the sixteen odd sign vectors in dimension five.  Pointwise
nonnegativity in (3) is therefore `n dot y<=1` for every odd sign vector `n`.
If

```text
D=conv{n in {+-1}^5: product_i n_i=-1},
```

then the moment polytope is exactly

```text
P=D^circ.                                                      (4b)
```

Here is a self-contained facet proof, rather than an appeal to a named
demicube.  One has

```text
D={z in [-1,1]^5: s dot z<=3 for every even sign vector s}.    (4c)
```

Indeed, an odd sign vector differs from an even one in at least one
coordinate, so every vertex of `D` satisfies the displayed system.  For the
reverse inclusion, maximize an arbitrary linear functional `c dot z`.  Put
`b_i=|c_i|`.  If some `b_i=0`, choose the signs of `c` on its zero coordinates
so that their product is odd; the cube bounds give `c dot z<=sum_i b_i`, with
equality at that odd sign vector.  The same argument works if the sign vector
of `c` is already odd.  Otherwise its sign vector `s` is even and all `b_i`
are positive.  If `b_m=min_i b_i`, combine `s dot z<=3` with the five cube
bounds `s_i z_i<=1` to obtain

```text
c dot z
 = b_m(s dot z)+sum_i (b_i-b_m)s_i z_i
 <= 3b_m+sum_i(b_i-b_m)=sum_i b_i-2b_m.                       (4d)
```

Flipping coordinate `m` of `s` gives an odd sign vector attaining the last
quantity.  Thus the two compact convex sets in (4c) have the same support
function in every direction and are equal.

Every coordinate equation `z_i=+-1` cuts out a four-dimensional crosspolytope:
the remaining four coordinates are the eight sign vectors of one parity,
four independent antipodal pairs.  These give ten facets of `D`.  For each
even sign vector `s`, equality in `s dot z<=3` holds at exactly the five odd
sign vectors obtained by flipping one coordinate of `s`; they are affinely
independent and give a four-simplex facet.  Thus these are all `26` facets.
Polarity turns them into the complete vertex list

```text
vert(P)={+-e_i: 1<=i<=5} union {(1/3)s: s in {+-1}^5,
                                            product_i s_i=+1}. (4e)
```

For completeness, two odd sign vectors have Hamming distance two or four.
A distance-two pair is an edge in a common simplex facet; a distance-four
pair is an antipodal nonedge in the coordinate crosspolytope determined by
their common coordinate.  Thus `D` has `16 choose(5,2)/2=80` edges.  Its
sixteen simplex facets have five tetrahedral ridges each, and its ten
crosspolytope facets have sixteen each; double counting gives `120` ridges.
Euler's relation gives `160` two-faces.  Reversing dimensions under polarity
therefore yields

```text
(f_0,f_1,f_2,f_3,f_4)(P)=(26,120,160,80,16).                  (4f)
```

This geometry identifies the full boundary stratification hidden by the two
compressed maxima in (4); it adds no realizability claim beyond (3).

## 4. Sharp functional norm and the H8/H12 atoms

The vertex classification gives a useful exact norm. For every real
five-vector c,

~~~text
sup_(y in P) |c dot y| = max{||c||_infinity, ||c||_1/3}.       (4g)
~~~

The signed coordinate vertices realize the first term. For the second,
choose signs matching c. Since the ambient dimension five is odd, exactly
one of that sign vector and its negative has even sign product; its
one-third vertex realizes ||c||_1/3 in absolute value. No interpolation
argument is needed because a linear functional attains its extrema at a
vertex.

The vertices also have literal orthogonal-array meanings. At each of the
ten vertices +-e_i, formula (3) is uniform on one eight-atom half-cube. At
each even-sign vertex s/3, the sixteen atom masses have multiplicities

~~~text
{0^5, (1/12)^10, (1/6)^1}.                                  (4h)
~~~

Thus it is the empirical law of a twelve-run strength-two array, with ten
atoms occurring once and one atom twice. Consequently every law in P is a
convex mixture of ten H8-type and sixteen H12-type four-column packets.
This is a local mixture statement; it does not supply a compatible square
Hadamard completion.

For the 896-row point y=(-1/4,0,0,-1/4,1/2), put

~~~text
s_- = (-1,-1,-1,-1,+1),       s_+ = (-1,+1,+1,-1,+1).
~~~

The triangular-bipyramid face exposes two different exact decompositions:

~~~text
y = (1/2)e_5+(1/4)(-e_1)+(1/4)(-e_4)
  = (1/4)e_5+(3/8)(s_-/3)+(3/8)(s_+/3).                      (4i)
~~~

At the row-count level these are respectively block sizes

~~~text
448+224+224 = 224+336+336 = 896.
~~~

Fourier inversion verifies that both sums give exactly the same sixteen
cell counts. Hence even the complete local moment law forgets which
extremal OA compiler produced it.

## 5. Exact integer and orthogonal-array form

Let a multiset of `N` four-bit rows have cell counts `n_x`, and define the
unnormalized moments

```text
A_i=sum_x n_x q(x)x_i,              D=sum_x n_x q(x).          (6)
```

The rows form an orthogonal array `OA(N,4,2,2)` if and only if their four
columns are balanced and pairwise orthogonal.  Applying (3) to the empirical
law gives the exact cell ledger

```text
n_x = [N+q(x)(D+A dot x)]/16.                                  (7)
```

Hence a proposed integer packet `(N,A,D)` is realized by such an array if
and only if all sixteen values in (7) are nonnegative integers.  There is no
search or hidden matching condition: those values are the unique cell
counts.  If `R` is the `N x 4` sign matrix of rows, the entire degree-two
shadow is simply

```text
[1 R]^T[1 R]=N I_5.                                            (8)
```

All arrays in the theorem therefore have the same normalized Gram data.
The five coordinates `(A,D)` are exactly what that Gram shadow discards.
In particular, every polynomial statistic of degree at most two has the same
expectation on every law in the cone.

Equation (8) is only four-column data.  It neither completes `R` to a square
Hadamard matrix nor asserts that an arbitrary local packet occurs in one.

## 6. An operation law

Let `X,Y` be independent laws satisfying (1), and multiply them
coordinatewise:

```text
Z_i=X_iY_i.                                                     (9)
```

For every Walsh character,

```text
E chi_S(Z)=(E chi_S(X))(E chi_S(Y)).                           (10)
```

Thus `Z` again satisfies (1), and its surviving coordinates obey

```text
a_i(Z)=a_i(X)a_i(Y),               d(Z)=d(X)d(Y).              (11)
```

For integer row multisets, the Cartesian product of runs gives the exact
compiler

```text
(N,A,D) star (M,B,E)=(NM,(A_iB_i)_i,DE).                       (12)
```

This operation preserves the whitened quadratic shadow while multiplying
the higher-chaos sidecar.  It is the lawful transfer; arbitrary multiplication
or projection of cell-count tables need not preserve strength two.

## 7. Complete convolution dynamics

Identify a law with its moment packet
`mu=(a_1,a_2,a_3,a_4,d)`, and let the superscript star `r` denote the
`r`-fold law obtained by independent coordinatewise multiplication. Equation
`(11)` gives

~~~text
mu^(star r)=(a_1^r,a_2^r,a_3^r,a_4^r,d^r).                    (12a)
~~~

Let u be the uniform law on the four-cube. Parseval applied to the density
relative to u gives the exact identity

~~~text
chi^2(mu^(star r) || u)
 = sum_i a_i^(2r)+d^(2r).                                    (12b)
~~~

Consequently

~~~text
TV(mu^(star r),u)
 <= (1/2)sqrt(sum_i a_i^(2r)+d^(2r))
 <= (sqrt(5)/2) rho^r,                                      (12c)

rho=max{|a_1|,...,|a_4|,|d|}.
~~~

The first bound is sharp on every signed one-coordinate segment
`mu=t e_i`, `|t|<=1` (where `rho=|t|`).
The vertex description makes the convergence dichotomy complete. Every
coordinate of P has absolute value at most one. If one has absolute value
one, any convex decomposition into (4e) can use only vertices attaining
that same coordinate equality; the unique such vertex is +e_i or -e_i.
Therefore

~~~text
rho<1:       mu^(star r) converges exponentially to u;
mu=+e_i:    every power is +e_i;
mu=-e_i:    odd powers are -e_i and even powers are +e_i.     (12d)
~~~

These ten exceptional packets are exactly the signed eight-run half-cubes.
All sixteen H12 vertices have rho=1/3 and hence mix. For an integer packet,
the r-fold Cartesian compiler has N^r runs and unnormalized high moments
(A_1^r,A_2^r,A_3^r,A_4^r,D^r); no limiting or probabilistic construction is
needed.

## 8. Exact packets exposed by the sign-word certificate

Three equal-four-row sidecars inside the reconstruction underlying THM-3394
have the following exact packets.  The two margins are the right side minus
the corresponding maximum in (4), after division by `N`.

| `N` | `A` | `D` | even margin | odd margin | multiset of the 16 cell counts |
|---:|---|---:|---:|---:|---|
| 48 | `(-16,-8,0,0)` | -8 | `1/3` | `2/3` | `{1^2,2^4,3^4,4^4,5^2}` |
| 120 | `(0,0,-8,8)` | 24 | `16/15` | `2/3` | `{5^2,6^4,7^2,8^2,9^4,10^2}` |
| 896 | `(-224,0,0,-224)` | 448 | `1` | `0` | `{0^2,28^4,56^4,84^4,112^2}` |

For the `N=896` packet, the absent atoms are

```text
(-1,-1,+1,-1),             (-1,+1,-1,-1).                    (13)
```

The corresponding active odd normals are

```text
n_1=(-1,-1,+1,-1,+1),      n_2=(-1,+1,-1,-1,+1).
```

Their common face in `P` has the five vertices

```text
e_5,  -e_1,  -e_4,
(1/3)(-1,-1,-1,-1,+1),  (1/3)(-1,+1,+1,-1,+1).               (13a)
```

They span dimension three.  Dually, `n_1,n_2` form an edge of `D`; the six
ways to choose one of their two differing coordinates and one of their three
common coordinates give the six triangular two-faces through that edge.
Consequently (13a) is a triangular bipyramid.  Direct substitution shows that
the packet has no active atom inequality besides `n_1,n_2`, so it lies in the
relative interior of this ridge, not merely somewhere on one odd-parity
facet.

Its degree-one and degree-two Walsh data still vanish identically.  This is
a sharp finite demonstration that whitening is not uniformity: sparse higher
interactions can survive after every linear and quadratic test is zero.

The familiar parity laws give the extremal controls `(a,d)=(0,+1)` and
`(0,-1)`, supported on only the even or odd half-cube.  At the other end, the
uniform law has `(a,d)=0` and full support.

## 9. Hostiles and transfer boundary

Both parity inequalities in (4) are load-bearing.  For

```text
a=(1,0,0,0),                d=1/2,                             (14)
```

the even-parity bound passes (`1<=3/2`) while the odd-parity bound fails
(`1>1/2`), and (3) has a negative atom.  Replacing `d` by `-1/2` exchanges
the roles.

The quadratic hypotheses are also load-bearing.  The uniform law on the
eight atoms satisfying `X_1=X_2` has balanced one-coordinate marginals and
has all five moments `(a,d)` equal to zero, but `E X_1X_2=1`.  Formula (3)
would incorrectly reconstruct the uniform sixteen-atom law if the degree-two
condition were omitted.

This theorem supplies an exact finite low-chaos filter: quadratic/Gram
observables are destroyed, and `(a,d)` is the complete repair.  It is useful
beside cubic-Hermite arguments and sign synchronization, but it does not by
itself transfer Gaussian positivity, improve the Grothendieck constant,
complete a Hadamard matrix, or solve any LRC/FC/JC frontier.

## 10. Exact companion

Run

```bash
python3 04-computation/four_bit_pairwise_independent_fourier_cone_thm3396.py
python3 -O 04-computation/four_bit_pairwise_independent_fourier_cone_thm3396.py
```

The standard-library companion:

1. checks exact Walsh orthogonality on all sixteen characters;
2. verifies the facet iff on `3,125` rational half-grid packets;
3. exhausts all `65,535` nonempty binary atom supports, finding exactly ten
   pairwise-independent half-cubes and the full cube;
4. reconstructs all `26` polar vertices from rank-five rational intersections
   and all faces from the `2^16` facet subsets, obtaining (4f);
5. verifies the sharp norm (4g) on a complete five-dimensional integer grid
   and reconstructs the ten H8 and sixteen H12 vertex packets;
6. verifies (12) on all `11^2=121` pairs of the binary laws;
7. checks the Parseval mixing identity and sharp total-variation bound for
   six convolution powers of all 61 feasible half-grid laws;
8. reconstructs the three displayed integer packets and their margins;
9. identifies the exact five-vertex, nine-edge, six-facet face containing the
   `N=896` packet; and
10. verifies both row-level decompositions in (4i) and checks the one-parity
   and missing-quadratic hostiles.

Normal and optimized runs are byte-identical.  The computation is an exact
referee for the proof and frozen packets; the theorem itself follows from
Fourier inversion and the pointwise inequalities above.
