# Hostile: binary Pfaffian coefficients are odd, not necessarily unit

Status: REFUTED proposed termwise lemma; strongest survivor and open repair
recorded below. This later hostile packet is outside frozen manifests.
No central mistakes ledger is modified here.

## Proposed claim and first failed implication

A proposed refinement selected principal Pfaffians containing exactly one
entire uncertain pair and otherwise only binary vertices. It asserted that
such a Pfaffian has the form +/-x+an even integer, and hence its square is
at least x^2 for |x|<=1. The first assertion is false: the coefficient is
the Pfaffian of a binary complementary tournament. Its parity is odd, but
its magnitude need not be one.

## Exact witness

Let the only uncertain skew edge be K01=x. On vertices2,3,4,5, set

    (K23,K24,K25,K34,K35,K45)=(1,1,1,1,-1,1).

Its Pfaffian is3. Set the eight cross entries

    (K02,K03,K04,K05,K12,K13,K14,K15)
       =(-1,-1,-1,-1,1,-1,1,1).

Complete by skew symmetry and zero diagonal. The full six-vertex Pfaffian
is exactly3x-2. At x=2/3 it vanishes, although x^2=4/9>0. Converting
K to the weighted tournament A=(J-I+K)/2 gives a valid matrix whose only
fractional unordered edge has probability5/6.

The exact determinant coordinates are

    u=32/9, v=-5/3, u-|v|=17/9.

The suggested TOTAL matching-support lower bound in this c=1 case would
be1-q=31/36. That total inequality is NOT refuted by this matrix: the
actual total margin17/9 is larger. Other principal Pfaffians carry excess
energy, so the missing sidecar is precisely how that excess can compensate
the failed individual square bound.

## Strongest survivor

If the complementary binary Pfaffian is known to equal+/-1 (for example,
for a transitive complementary tournament), the coefficient is a unit and
the intercept is even. Then |+/-x+2z|>=|x| for |x|<=1 and integer z, so the
termwise square bound is valid in that scoped case.

For arbitrary binary complements, the unit-coefficient argument is
invalid. A compensating total-margin refinement remains OPEN. The accepted
multiplicative Pfaffian margin proof uses Gram/Cauchy and the entire
sum-of-squares expression, not this unit assumption, and is unaffected.

## Reproduction

    python3 -B pfaffian_coefficient_hostile.py > pfaffian_coefficient_hostile.json
    python3 -B -O pfaffian_coefficient_hostile.py > /tmp/hostile.optimized.json
    cmp pfaffian_coefficient_hostile.json /tmp/hostile.optimized.json

The script uses independent Fraction Gaussian elimination and a recursive
Pfaffian expansion. It checks coefficient3, intercept-2, the exact zero,
and the full determinant margins. No numerical approximation is used.
