# Braid-arrangement shadows of the missing-region program

> **SCOPE CORRECTION — MISTAKE-223 (2026-07-21).** S208 found exact
> arrangement and determinant shadows, but promoted them to transfers that do
> not follow. Braid localization does not factor the NC2 hyper-Bessel boundary;
> the cake--bagel and finite-shadow deficits have not been realized by one Euler
> valuation; and THM-2023 had already proved the general `Phi_(p,q)`
> Laguerre--Polya statement. Read the separated facts below, not the original
> “four pillars verified” claim.

*boxeph-2026-07-21-S208. Original aim: look for geometry beneath the
finite-rank missing-region law and test whether it supplies algebraic leverage.*

## 1. Exact braid-arrangement fact

The polynomial

```text
V(a_1,...,a_n) = product_(i<j) (a_j-a_i)
```

is the defining polynomial of the braid arrangement `a_i=a_j`. Its real
complement has `n!` chambers, one for each total order of the coordinates and
hence one for each labelled transitive tournament. Flats correspond to set
partitions recording which coordinates coincide. The complex complement is
connected, so the `n!` chamber statement must not be attached to it.

THM-2033 identifies a Vandermonde in one special factorial moment matrix. In
that scope, distinct nodes put the determinant off the arrangement. It does
not follow that an arbitrary scalar NC2 moment is nonzero exactly on the braid
complement; ACTIVE-GUARDRAILS item 18 and MISTAKE-215 forbid that upgrade.

## 2. Exact coalescence formula, with the right exponent

Let coincidence blocks be `B_1,...,B_m`, and write the coordinates in block
`B_i` as `c_i+epsilon*delta_j`, with distinct block centres. Put

```text
N_X = sum_i C(|B_i|,2).
```

Then, up to the sign determined by coordinate ordering,

```text
V(c+epsilon*delta)
 = epsilon^N_X
   * product_i V(delta restricted to B_i)
   * product_(i<j) (c_j-c_i)^(|B_i||B_j|)
   + O(epsilon^(N_X+1)).
```

This is an exact leading expansion. `N_X` is the Vandermonde's vanishing order,
not the flat codimension: for block sizes `3+2`, they are `4` and `3`,
respectively. In arrangement localization the within-block braid arrangements
are the local factors; the cross-block product is a nonzero unit near a generic
point of the flat.

This formula is valuable wherever an actual determinant has already been
reduced to `V`. It does not produce a factorization of a different analytic
object merely because that object occurs near a “wall.”

## 3. Why the hyper-Bessel transfer failed

THM-2017's

```text
Phi_(p,q)(x) = sum_(k>=0) x^k / ((pk)! (qk)!)
```

arises as a univariate channel limit. S208 supplied no map from braid blocks or
the preceding polynomial asymptotic to this series, so the claimed product of
“one hyper-Bessel per block” was a type jump. Sampling the necessary Laguerre
inequality also cannot prove Laguerre--Polya membership.

The desired conclusion is nevertheless known for an independent reason:
THM-2023 uses Gauss multiplication to express `Phi_(p,q)` as a positive-
parameter generalized hypergeometric function, and the Baricz--Singh theorem
puts all its zeros on the negative real axis. Arrangement product closure is
neither needed nor established.

## 4. Exact companion determinant, without a topology claim

For the recurrence

```text
a_d = a_(d-1) + a_(d-g-1),     a_0=1,
```

take the `(g+1)`-state companion matrix whose first row is `(1,0,...,0,1)`
and whose subdiagonal entries are one. Its characteristic polynomial is
`lambda^(g+1)-lambda^g-1`, hence

```text
det(I-x*M_g) = 1-x-x^(g+1).
```

Thus the full-rank gap sequence has generating function
`1/det(I-x*M_g)`, a legitimate Bowen--Lanford-style determinant/zeta syntax.
The companion digraph has a loop and is not a tournament; common syntax does
not identify the underlying objects.

The finite-rank shadow's first deficit and the bagel identity remain separate.
For positive `n`, the stronger direct relation is

```text
bagel(n) = cake(n+1)-2,
bagel(n)-cake(n) = (T_n+1)-2 = T_n-1.
```

[OEIS A003600](https://oeis.org/A003600) describes the solid-torus correction
as an extra cake slit that leaves two pieces unseparated. A solid torus is
`S^1 x D^2`, not `T^3`. No chain complex, intersection poset, Mobius function,
or valuation has been supplied that sends this two-piece correction and the
shadow's first missing cell to the same class.

## 5. The repaired research target

The topology idea remains testable, but the missing data are now explicit:

1. Construct a cell or arrangement model for the bagel cut count and derive
   `cake(n+1)-2` by a named valuation.
2. Construct the finite-shadow cell complex and derive its first deficit by
   the same kind of valuation.
3. Give a source-to-target map, show what cells and incidence data it preserves,
   and account for every lost boundary term.
4. If no such map exists, isolate the first incompatible operation; that
   obstruction is more informative than another matching scalar.

Potential tournament vertices were chambers, flats, coincidence blocks,
companion states, and proof obligations. Only real braid chambers carry the
canonical total-order tournament. Orienting the other objects would discard
the incidence or analytic data at issue, so no new tournament quotient is
asserted.

## Reproduction scope

Run:

```bash
python3 04-computation/arrangement_topology_leverage_boxeph_S208.py
```

The repaired script checks the braid count, one exact coalescence family, the
explicit companion determinant, and small recurrence coefficients. It prints
the failed transfers as non-results. It computes no Euler characteristic of a
bagel-cut complex and proves no hyper-Bessel factorization.

Related: HYP-8825, MISTAKE-222, MISTAKE-223, THM-2023, THM-2033, and the
corrected S207 binomial-atlas reflection.
