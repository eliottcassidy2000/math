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
  is a complete minimal sidecar.  Three exact Hadamard-puzzle packets of sizes
  48, 120, and 896 are verified; the 896 packet lies on a sharp odd-parity
  facet and omits exactly two atoms.  This is a four-column/local realizability
  theorem, not a Hadamard completion theorem or a Grothendieck bound.
source: root-2608-sign-puzzle-2026-08-14
audit: independent Fourier/OA/convolution derivation, raw-signword packet extraction, facet/hostile reconstruction, and normal/optimized replay
related:
  - THM-3394-twelve-formerly-missing-hadamard-orders-through-2000
  - THM-3392-bipartite-sign-lift-and-synchronization-loss
script: 04-computation/four_bit_pairwise_independent_fourier_cone_thm3396.py
output: 05-knowledge/results/four_bit_pairwise_independent_fourier_cone_thm3396.out
script_sha256: a1472ad8ac7506d0513711c21f4be3b6dc3281b82f825f320329336e9a18d7c3
output_sha256: aa2b12178e658ae2aca4d9ba7f5b688e654b0027490febe60618bc1c9dac9d2b
semantic_sha256: 3a13dc0acd6929f3fdba1ec646c045a639dee4c34f729be19d51a214dcc3dee5
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

## 3. Exact integer and orthogonal-array form

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

## 4. An operation law

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

## 5. Exact packets exposed by the sign-word certificate

Three equal-four-row sidecars inside the reconstruction underlying THM-3394
have the following exact packets.  The two margins are the right side minus
the corresponding maximum in (4), after division by `N`.

| `N` | `A` | `D` | even margin | odd margin | multiset of the 16 cell counts |
|---:|---|---:|---:|---:|---|
| 48 | `(-16,-8,0,0)` | -8 | `1/3` | `2/3` | `{1^2,2^4,3^4,4^4,5^2}` |
| 120 | `(0,0,-8,8)` | 24 | `16/15` | `2/3` | `{5^2,6^4,7^2,8^2,9^4,10^2}` |
| 896 | `(-224,0,0,-224)` | 448 | `1` | `0` | `{0^2,28^4,56^4,84^4,112^2}` |

The `N=896` packet is therefore an exact facet witness.  The absent atoms are

```text
(-1,-1,+1,-1),             (-1,+1,-1,-1).                    (13)
```

Its degree-one and degree-two Walsh data still vanish identically.  This is
a sharp finite demonstration that whitening is not uniformity: sparse higher
interactions can survive after every linear and quadratic test is zero.

The familiar parity laws give the extremal controls `(a,d)=(0,+1)` and
`(0,-1)`, supported on only the even or odd half-cube.  At the other end, the
uniform law has `(a,d)=0` and full support.

## 6. Hostiles and transfer boundary

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

## 7. Exact companion

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
4. verifies (12) on all `11^2=121` pairs of those laws;
5. reconstructs the three displayed integer packets and their margins; and
6. checks the one-parity and missing-quadratic hostiles.

Normal and optimized runs are byte-identical.  The computation is an exact
referee for the proof and frozen packets; the theorem itself follows from
Fourier inversion and the pointwise inequalities above.
