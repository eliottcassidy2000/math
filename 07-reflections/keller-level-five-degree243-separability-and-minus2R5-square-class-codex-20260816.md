# The 243-sheet gate is now an eliminant theorem, not a tower dimension

**Status: PROMOTED as THM-3525 after two exact representations.**  At the
lawful target `(1,1,1)` over `F_251`, the fixed Keller fifth x-eliminant has
degree `243` and derivative gcd one.  The generic fifth eliminant therefore
has 243 distinct roots, and its discriminant square class is `[-2R_5]`.
This does not prove that `R_5` is a fifth image prime.

## Inheritance pass

- **Closest proved mechanism:** THM-3498 used a degree-81 good fibre to make
  the fourth THM-2582 odd-block step lawful and obtained `[Delta_4]=[2G]`.
- **Canonical hostile:** THM-3523's 243-dimensional evaluation algebra is a
  complete norm universe but says nothing by itself about the degree or
  squarefreeness of a polynomial in a fresh fifth x-coordinate.
- **Corrected near miss:** a full tower dimension is not a primitive-element
  or separability certificate; the new variable must survive with degree 243
  and nonzero discriminant.
- **Least-used sidecar:** multiplicative Fourier inversion on `F_p^*` gives
  an interpolation method disjoint from consecutive-point Newton recovery.

## What changed

The new object is

```text
P_5(X)=Norm_(A_4/K)(L_4X^3+T_4X-2z_4),
dim_K A_4=81.
```

At `p=251`, both exact routes recover a 244-coefficient polynomial with hash

```text
912f32ec0b9b375d9db2ba71d7fdf224456c86862871ae0f8f92bac5038f00ab.
```

The first route uses nested coefficient tuples, transitive cubic norms, and
Newton interpolation at `0,...,243`.  The second uses FLINT regular matrices,
all 250 nonzero field values, and character inversion on `F_251^*`.  It sees
zero coefficients in degrees `244,...,249`, a nonzero coefficient in degree
243, and derivative gcd one.  The fibre factors with degree multiplicities

```text
1^2,2^6,3^4,4^4,6^4,9,12^3,24,36^3,
```

all exponent one.  This is a squarefree reducible fibre, exactly the right
kind of witness: irreducibility would be stronger but is unnecessary.

## Why the consequence is generic

All four inverse graphs are checked and every `L`, chart denominator, and
cubic derivative used by the tower is a unit.  The reduction is therefore
inside the finite-etale chart.  Nonzero reduction of the leading coefficient
and discriminant makes both generic characteristic-zero functions nonzero.
After splitting `A_4`, the norm is a product of 81 cubic blocks; a squarefree
degree-243 product forces each block to remain cubic and separable and every
pair to be coprime.

This distinction is the main conceptual gain:

```text
dimension 243 evaluation algebra
  != degree-243 eliminant;

degree-243 squarefree norm polynomial
  => 81 cubic blocks are generically disjoint and separable.
```

## Square-class closure

The THM-2582 recursion can now be applied one more time:

```text
[Delta_5]
 =[N(Delta_4)][Delta_1]
 =[N(2G)][-L]
 =[-8LN(G)]
 =[-2R_5].
```

The alternating sign and the constant class `[2]` continue.  The exponent
`270` of the old `L` denominator is even, so the old component again vanishes
from the square class.  What remains is the newest cleared norm, but no
irreducibility of that norm is smuggled in.

## New frontier

The degree-243 gate is closed.  The fifth-image program now has a cleaner
separation:

1. find a tractable slice or other witness showing `R_5` is a first power of
   its image prime and coprime to `L,H,J,G`;
2. prove absolute irreducibility or otherwise identify its image divisor;
3. only then promote a fifth nonproperness component;
4. independently, a sixth discriminant step would require a new degree-729
   separability witness.

The current theorem proves none of those four later statements.  It removes
the generic block-separability obstruction that previously sat in front of
them.
