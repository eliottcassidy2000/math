# The scalar barrier survives restricting to strong witnesses

The [padding proof](SCALAR_DESCENT_BARRIER.md) treats all tournament
witnesses. Restricting to strong ones removes the free ordinal tail,
so that proof alone does not apply. A different proof uses the
[homogeneous switch from Session 319](../session319_arc_usage_and_signed_cycles_20261113/RESERVOIR_POLYNOMIAL.md).

For every m>=2 it gives strong A_m,B_m of the same order and the
same count 9(2^m+1), with ordered nonempty cover responses
```
F_1(A_m)=F_1(B_m),   F_j(A_m)<F_j(B_m)  (2<=j<=m+2),
```
and equality at higher j. Fix a homogeneous context with any
positive number of copies of its input and all other factors fixed.
The itinerary expansion has only nonnegative coefficients, so
every term is weakly larger for B_m. If any marked vertex lies
on a quotient directed cycle, the [insertion
argument](SCALAR_DESCENT_BARRIER.md) gives an itinerary using that
mark twice and every other vertex once. That term is strictly
larger, since all its other factors are positive and F_2 differs.
This works for any number of copies: all higher coordinates
are weakly ordered, so no bound on the visit count is needed.

Thus the full monomial classification holds even when the
input is promised strong and both compared inputs have the
same order. This rules out a putative uniform strong scalar
lift of 2h+1 even after eliminating all ordinal padding.

For a context chosen solely by h modulo a fixed M, some residue
class contains infinitely many counts 9(2^m+1) from this family.
Every class hit by a member of this family with m>=2 must use a
singleton-SCC context if it is to be a uniform scalar rule on
that class. The corollary is only for these residues; it does
not constrain a branch using additional port information.
