# Fifteen pure switches and a thirty-tick Paley count clock

This combines the prior 5–6–7 curvature/rotation-carrier work, THM-4762's
Fibonacci cat map, and the signed Pythagorean-conic and projective-switching
notes. The concrete new control below is a finite five-vertex orbit; the
determinant law applies at all orders. No Chan–Pak selector is asserted.

## Three different meanings of 15 and 30

1. The inherited (2,3,q) triangle has angular excess
   pi(1/2+1/3+1/q−1). At q=5,6,7 this is pi/30, 0, −pi/42:
   spherical, Euclidean, and hyperbolic respectively. At 5 there are
   120 reflected chambers on the sphere and 60 in the orientation-preserving
   group A5. The 15 is also the number of involutions in A5 (choose the fixed
   label and one of the three double transpositions on the other four).
2. Put M=((1,1),(1,0)), the classical Fibonacci matrix (Wall,
   https://doi.org/10.1080/00029890.1960.11989541). Its matrix order is the
   Pisano period. Modulo 15 it is 40; modulo 30 it is 120, by CRT from
   the orders 3,8,20 modulo 2,3,5. If we instead stop at a unit scalar
   matrix, the orders are 20 and 60, respectively, and the scalar is
   11 in both rings. The corresponding local scalar orders are 3,4,5.
   Since det M=−1, no ordinary Pisano period15 exists for any modulus
   above2. Ordinary and scalar-quotiented periods are different.
3. Modulo 31, M has both linear and projective order30, while M² has
   order15 in both. Modulo 61, M has linear order60, scalar order15,
   and M^15=11I, with 11²=−1. In particular M² has linear order30
   but scalar order15 there. These order facts can be verified by the
   Fibonacci power formula and testing the divisors of the displayed
   periods, or by the exact replay.

The recurring numbers measure different things. The common structural
question is what additional data the return must preserve: orientation
of a spherical chamber, a signed vector, or just a projective line.

## Conic and hyperbolic transition

The form Q(x,y)=x²−xy−y² satisfies Q(Mv)=−Q(v), so M² preserves it.
Over the reals its signature is (1,1); M² has eigenvalues phi² and
phi^(−2), the hyperbolic boost from THM-4762. Independently, applying
the prior Pythagorean quadratic lift V(u,v)=(u²−v²,2uv,u²+v²)
carries a projective spinor to the conic X²+Y²=Z²; its Lorentz pairing
squares the alternating determinant. Thus the conic alone forgets the
sign needed for the tournament. The real triangle curvature above
is not a curvature assigned to a finite-field graph.

## What the Fibonacci clock actually does to tournaments

Over an odd field with q≡3 mod4, use the projective chart with representative
a(t)=(1,t), and a(infinity)=(0,1). Choose distinct points x_i and orient
i→j iff chi(det(a(x_i),a(x_j)))=+1. This is Paley with an added sink
when all projective points are used. Define the normalization factors
lambda_i by M^k a(x_i)=lambda_i a(M^k x_i). Then

    chi(det(a(M^k x_i),a(M^k x_j)))
      = chi((-1)^k) chi(lambda_i) chi(lambda_j)
        chi(det(a(x_i),a(x_j))).

Thus even steps are pure cut switches; odd steps are switches followed
by converse. Every step gives a real tournament in this chart, including
when a point hits infinity. If one keeps the unnormalized signed
vectors instead, the cut factors disappear and the count is constant
(converse preserves H). No coordinate shift in this orbit changes the
switching class up to converse. This is the inherited projective
determinant law specialized to the Fibonacci clock.

Over F31, take (0,1,2,3,4), with projective update t↦1/(1+t) and
the usual values at the pole and infinity. In the normalized chart its
Hamiltonian counts over one period are exactly

    13,15,5,3,13,9,9,15,5,5,5,15,9,3,5,
    5,5,5,5,9,5,5,15,3,5,5,9,13,3,5.

The count clock itself has minimal period30 (compare each proper
divisor). Its even ticks, the 15-step M² orbit, already see 3,5,13.
At 61 the determinant character is symmetric because −1 is a square;
the superficially similar order15 clock does not orient a Paley
tournament in that field.

For Chan–Pak this is a bounded catalogue of physically reachable
switches, not an absolute selector. The switching-sheet obstruction
still applies, and one cyclic subgroup need not even exhaust its sheet.
A productive next move would change triple products, then use a
short clock of this kind to optimize which cut reads the new source.

## Prior repository work and replay

- THM-4762: Fibonacci matrix, cat map, and Pisano periods.
- geometry_567_flag_carriers_20261004.md and
  ogg_triangular_chowla_20260927.md: excess and the role of rotation data.
- CONIC_SIGN_LOSS_AND_TARGET_LIFT.md and
  PROJECTIVE_SWITCHING_SELECTOR.md: signed conic lifts and cut switches.

Run python3 pisano_switch_replay.py. The script checks matrix and scalar
orders and enumerates the literal 120 permutations at each Paley tick.
