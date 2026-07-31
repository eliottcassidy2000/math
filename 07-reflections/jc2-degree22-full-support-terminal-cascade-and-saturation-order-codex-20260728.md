# The full-support obstruction is a two-order terminal cascade

**Scope:** independent exact degree-twenty-two `BCDEW` projection referee
downstream of THM-2411, THM-2671, THM-2683, and the now-promoted THM-2692.
The canonical THM-2692 companion separately proves the stronger identical,
squarefree degree-seven carrier and unit-reconstruction statements; this
referee supplies a direct characteristic-zero unit-ideal proof of terminal
emptiness over both factor fields.  Nothing here closes split/even descent,
integral raising, another degree branch, `JC(2)`, or `DC(2)`.

## Inheritance and the missing coordinate

The support-three calculation closes a chart when the first coefficient past
the eliminant's `t`-degree has full monomial rank.  THM-2683 repairs the same
idea on support four: when the ambient terminal kernel survives, intersect it
with the actual monomial image through its toric equations.

On full support, neither view alone is final.  The order-eleven equations have
nine monomials and rank four.  Their intersection with the four-dimensional
coefficient torus is not visibly empty.  The coordinate missing from the
one-order picture is the **next Hensel jet**.  Order eleven compresses the
torus to a finite projected candidate set; order twelve is incompatible with
that projected set and removes it.

This suggests a reusable ladder:

```text
terminal linear rows
    -> intrinsic monomial/toric reparametrization
    -> pivot-compatible finite eliminant
    -> next-order compatibility projections
    -> exact unit ideal.
```

The full-support case is therefore not a failed rank argument.  It is the
first genuine two-order member of the same terminal-obstruction family.

## The intrinsic four-variable object

The complete order-eleven support is

```text
(E, EW, DE, C, CW, CD, CD^2, C^2E, C^3).
```

On the active torus put

```text
A=E/C,                       S=C^2.
```

After dividing by nonzero `C`, this becomes

```text
(A, AW, AD, 1, W, D, D^2, AS, S).
```

Each of four independent scalar rows consequently has the form

```text
w_i(A) W + s_i(A) S + b_i(A,D)=0,
```

where `w_i,s_i` are affine in `A` and `b_i` has support
`(1,A,D,AD,D^2)`.  This is the structural compression: a nine-coordinate
monomial vector becomes a four-row affine pencil in the two fibre variables
`W,S`, over the two-coordinate base `(A,D)`.

For a row pair `(i,j)`, the Cramer pivot

```text
Delta_ij=w_i s_j-w_j s_i
```

is quadratic in `A`.  On `Delta_ij != 0`, Cramer solves `W,S`.  The other two
rows become two compatibility polynomials `R_ij,R'_ij` of `D`-degree two and
`A`-degree at most three.  Their exact `D`-resultant has degree eleven in
`A`; saturation by `Delta_ij` leaves degree seven on every one of the six
pivots, in both the root and unordered-pair fields.

Thus order eleven is best understood as a finite spectral atlas of at most
seven physical `A`-candidates per pivot, not as an ambient kernel.

## Why order twelve is exactly shaped to kill the atlas

The order-twelve support is

```text
1, W, W^2, E^2, D, DW, D^2, D^3,
CE, CDE, C^2, C^2W, C^2D, C^4.
```

After `E=AC`, `S=C^2`, and the Cramer formulas

```text
W=N_W/Delta,                S=N_S/Delta,
```

every monomial has denominator exponent at most two.  Multiplication by
`Delta^2` therefore gives an honest polynomial numerator of `D`-degree four
with exactly twenty-eight support terms.  Each of the five order-twelve rows
is projected against both compatibility equations.  The resulting ten
univariate polynomials have exact degree twenty one.

For every pivot, the exact gcd of the degree-seven saturated terminal
resultant and all ten order-twelve projections is one.  A hypothetical
factorization would supply a common `(A,D,W,S)` point, hence a common zero of
all eleven univariate polynomials, contradicting that unit ideal.

The six quadratic pivots themselves have exact gcd one.  Hence they have no
common `A`-zero, so every order-eleven solution lies in at least one Cramer
chart.  This closes the atlas without assuming a preferred row pair.

## The reduction-order trap

The most informative hostile control was a failed residue.  In the root
field at `F_109`, generator `43`, pivot `(0,2)`, the exact saturated resultant
has degree seven, but reduction creates an additional common factor with
`Delta`.  Saturating **after** reduction produces degree six.  That smaller
polynomial has a unit-ideal certificate in the special fibre, but it is not
the reduction of the characteristic-zero saturated polynomial and cannot
prove the desired statement.

The repaired order of operations is:

```text
exact resultant -> exact saturation -> coefficientwise reduction,
```

never

```text
reduction -> special-fibre saturation.
```

At the good simple residues `F_113` with root generator `24` and `F_103`
with pair generator `61`, direct reductions of the exact saturated
polynomials do give full coefficient-convolution rank.  All denominators are
checked nonzero, independently recomputed residue polynomials agree up to
units, and named square minors are nonzero.  These are lifted determinant
witnesses, not an inference that an empty special fibre alone forces an empty
generic fibre.

## What this sharpens

The useful object is not merely the Hensel matrix, the toric variety, or a
resultant.  It is a **filtered terminal scheme**: its order-eleven truncation
has a finite projected support after toric localization, while the
order-twelve equations have no common point with that support.  This view
explains both why one order was insufficient and why only one additional
order was needed here.  It does not assert a Jacobian-transversality or
reduced-intersection calculation; those stronger carrier properties belong
to the canonical THM-2692 companion.

For future branches, the cheapest decisive test is now clear:

1. quotient the terminal monomials by the active scaling torus;
2. solve the smallest affine fibre using all pivot charts;
3. saturate in characteristic zero by every chart denominator;
4. project the next nonzero jet against every base compatibility equation;
5. test the resulting exact univariate ideal, retaining a named lifted minor
   as an independent certificate.

The equality boundary is also sharp: a common zero of the saturated terminal
polynomial and next-order projections is exactly where this method stops and
a third jet or a different sidecar would be required.  In the present
`BCDEW` chart no such zero survives in either fixed-factor field.
