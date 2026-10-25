# The twist average, its formal seven, and a positive-torus obstruction

2026-10-25. **PROVED** by finite-field character orthogonality and labelled
line incidence. This applies the already proved negative-torus-coefficient
test from Session59's `arithmetic-configurations/RESULT.md`; that test is
inherited, not rediscovered here.

Let V be the complement of the seven labelled hourglass lines from
`KUMMER_AND_F1.md`, and let q be any odd prime power. Every point of V has
E!=0. For each of the 64 vectors delta in (F_q*/F_q*^2)^6, let X_delta^o
be the unramified cover given by

    u_i^2=delta_i L_i/L_0,  i=1,...,6.

Then the exact twist average is

    (1/64) sum_delta #X_delta^o(F_q) = #V(F_q).                (A1)

At a fixed source P, exactly one twist vector makes all six radicands
squares; its fibre has 64 points. The other 63 fibres are empty. Summing
first over twists proves A1. This does not assert that any one fixed twist
has the average number of points, or produce a rational point over Q.

For either complementary norm-seven pair, its six twist coordinates are
bijective linear coordinates on this same twist set. The product of the
two quotient fibre sizes equals the full fibre size **at each common P**.
That product cannot be moved outside the sum over P. Separate marginal
point counts would instead count pairs of possibly different base points.

## The seven-line polynomial, including its bad characteristic

The seven primitive line coefficient rows are

    (1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(1,1,1),(1,-1,-1).

Their 35 triple determinants lie in {0,+/-1,+/-2,3}. In characteristic
different from 2 and 3 the only concurrent triples are 012,034,056.
There are therefore twelve ordinary double points and three triple points,
all rational over the prime field. No four lines meet. Inclusion-exclusion
on the union, or subtracting the excess multiplicity at each intersection,
gives

    #V(F_q)=q^2-6q+12,      char(F_q)>3.                     (A2)

In characteristic three there are exactly two more concurrent triples,
136 and 245. There are six ordinary double points and five triple points,
still no fourfold point, giving the repaired formula

    #V(F_q)=q^2-6q+10,      char(F_q)=3.                     (A3)

These are all-extension formulas, not merely fits to prime-field data:
all intersections are already defined over the prime field, and the
displayed determinant list identifies all exceptional characteristics.

For A2 the formal value at q=1 is seven. But with t=q-1 the polynomial is

    t^2-4t+7;                                                (A4)

A3 similarly gives t^2-4t+5. A finite decomposition of V over a finite
field into split tori would give sum_j a_j(q^r-1)^j with a_j nonnegative
integers for every extension degree r. Infinitely many extension sizes
force equality of polynomials, contradicting the coefficient -4.
Thus neither complement admits such a split-torus decomposition.

This is an obstruction to that precise positive torus model, not to every
F1 framework. Nor does it contradict the F1 **ambient torus isomorphism**
in K3. The faithful object is the pullback

    X_delta^o = V x_(T_b, b_delta) T_u,   T_u -> T_b: u -> u^2.

The root-power map and repaired coordinate transformations are multiplicative;
the branch map b_delta contains the indispensable additive equations E+x,
E-z, and so on. Omitting the branch map changes the object. The formal seven
counts a virtual combination of torus pieces, not seven actual rank-zero
pieces. The same positivity warning is visible in the incoming tournament
example: a high root multiplicity need not leave a quotient polynomial with
nonnegative coefficients. These are two explicit obstructions to interpreting
a scalar specialization as a decomposition; they are not isomorphic objects.
