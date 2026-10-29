# Session99: Collatz clocks and tournament singularities

**Status:** exact bridges, independently audited, with explicit hostile controls.
Source pin: `d9d7e33349653b6b65dbad477841c9ed518c081a`.
This session improves the mathematical connection; it does not prove Collatz
convergence. The native integer source and its arithmetic guards are retained.

## 1. More than the same numeral: a pole amplitude reads the halving clock

The normalized tournament GF from Session98 has pole p_k=4^k of order3k+1.
Its leading scaled amplitude is

\[
 A_k=\frac{(3k)!}{3^k k!}\,2^{3k^2+k}.
\]

The exact ladder is

\[
 A_{k+1}=2^{6k+4}(3k+1)(3k+2)A_k.
\]

For **odd k**, the cofactor3k+2 is odd. Therefore the actual accelerated
Collatz division clock is recovered by

\[
 v_2(3k+1)=v_2(A_{k+1})-v_2(A_k)-(6k+4).
\]

Together with the pole order, this determines the next Syracuse source.
Equivalently, the explicitly defined pole-order readout
4^k ->4^oddpart(ord_(4^k)B) is conjugate to odd Collatz. This is a precise
arithmetic encoding, not an analytic-flow or descent theorem.
[Proof, binary-digit formula and parity guards](readout/POLE_ORDER_READOUT.md).

The residual identity is
v_2(3k+1)=2+H(k+1)-H(k), H(k)=s_2(k)-s_2(3k).
It does not telescope along Collatz trajectories, whose next index is generally
not k+1. This identifies exactly which alignment a potential argument lacks.

## 2. Why both constructions have 3k+1 geometric slots

A maximal tournament core3^k has3k labeled support positions and3k+1 gaps
for fixed vertices. An ordered full ternary tree with k branch vertices has
3k edges and3k+1 vertices. Fixing a core label order and one tree shape gives
an explicit bijection between gap occupancies and vertex residual clocks.

After removing prefactors and matching variables, both give the negative
binomial generating factor(1-t)^-(3k+1). At the critical parity controller's
residual law t=1/4, this identifies the exact clock distribution, not just
its mean. The native clock variable still has poles at +/-2, not4^k.
There is no corresponding equivariant identification of all labeled core
symmetries with ordered tree symmetries.
[Geometric-slot theorem and obstruction](structural/GEOMETRIC_SLOT_BRIDGE.md).

## 3. Genuine analytic growth is a different operation

For the top-core tower T(z,u), with D=z*d/dz,

\[
 \partial_u T(z,u)=\frac{z^3}{12}(D+1)(D+2)(D+3)T(z/4,u).
\]

This transports pole locations by a factor4 and orders by+3, while the
core index moves k->k+1. It is different from the arithmetic readout above.
The full core insertion law retains gcd-response coordinates; cores(3,7)
and(5,5) have the same size and defect but different responses to insertion3.
[Exact insertion operators and all-ell towers](poles/POLE_TRANSPORT.md).

The same construction for ell=5 has pole orders5k+1, and its analogous
odd-part readout contains13->33->83->13. Thus meromorphic continuation,
positive coefficients and geometric-slot laws cannot alone establish native
convergence. This is a hostile example built by the same mechanism.

## 4. Keeping the additive carry changes the pole order

For the cheaper controller, a valuation divisible by3 creates three children;
other valuations are leaves. A stopped tree with k branches has J=3k+1
valuation slots. Its affine map is

    n -> (3^J n+C)/2^D.

The raw carry moment GF of order r, for fixed k>=1, has positive-real poles

- z=1, exact order3k+1;
- z=2^-e, e=1,...,r, exact order3k.

The missing slot has a precise cause: the terminal valuation enters the
final division clock D but never enters the numerator carry C. The proof
tracks marked chronological prefixes and excludes cancellation of the
highest pole. Two words with the same branch count and total clock have
carries179 and211, and different original-source cylinders.
[Carry moment theorem and source controls](transfer/CARRY_CLOCK_POLES.md).

The critical parity rule in section2 and this cheaper mod3 rule are DIFFERENT
controllers. Their branching probabilities, residual clocks and pole
coordinates must not be interchanged.

## 5. A new moment boundary for the cheaper rule

Under normalized odd Haar source measure, write U=3^J/2^D and V=C/2^D.

- All slope moments exist: U<=(3/2)(27/32)^k<=3/2.
- At fixed k, every V moment is finite, although every positive integer raw C moment
  is infinite when k>=1.
- After mixing over all tree sizes, **EV, EV^2 and EV^3 are finite, but EV^r
  is infinite for every integer r>=4**.

The exact criterion is the child amplification factor
(1+m_r+m_r^2)/7<1, with m_r=EU^r the least nonnegative root of a cubic.
The proof uses completed-tree truncation, rather than assuming finiteness.
The mean7/12 is inherited; new algebraic second moments give approximately
EU^2=1.30746709249, E(UV)=0.53890572060, EV^2=1.19925614077.
[Moment boundary and exact equations](transfer/INTERCEPT_MOMENT_BOUNDARY.md).

These are source-measure statements, not permission to resample a fixed
positive integer's future or to replace its exact carry by a mean.

## Reproduction

Run `python3 verify_packet.py` in this directory. It runs each mathematical
checker normally and under `python3 -O`, compares stdout and stored receipts,
and leaves checks active in both modes. The two symbolic carry/moment
checkers require **SymPy1.14.0** and mpmath1.3.0 (see requirements.txt); the remaining checkers use the standard
library. Independent proof review and source-cylinder tests are in
[audit/AUDIT.md](audit/AUDIT.md).

Finite computation supports but does not prove the all-order identities,
analytic continuation or moment-divergence statements. The notes give their
proofs separately. No claimed progress depends on a numerical coincidence
between tournament counts and the values of a Collatz orbit.
