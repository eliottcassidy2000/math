# A quarter-price source bound and an eighth-power ROOT error

2026-10-11. **PROVED, scoped and self-audited:** Q1--Q5 below.
**FINITE-EXACT:** the companion replay. **OPEN:** an independently established
positive lower readout at every source, and universal Collatz. This changes
the measured weight; it is not an eighth-power estimate for the old W-readout.

## Inheritance and the retained source

The [adaptive mixture](collatz_adaptive_mixture_flow_20261005.md) already
constructs every fixed-price flow f_r, its critical incoming inequality,
summability, and exact positive support. The [refuel bill](collatz_refuel_bill_code_20261005.md)
already singles out r=1/4: its base-edge size comparison charges one bit to
a valuation-one base edge. These are inherited, not new weight constructions.
The [height envelope](collatz_height_envelope_20261011.md) gives the sharp
logarithmic cap for the beta mixture W. The [deadline theorem](collatz_floor_transport_deadlines_20261005.md)
already contains the stronger source-height bound behind that envelope.
The [measurement audit](collatz_measurement_independence_20261005.md) forbids
calling a certified off-source bank an independent marked-source measurement.

The present additions are an exact *ordinary-source* envelope for the
quarter-price flow, its affine normalization and defect product, an explicit
q^-8 alias bill, and a quantitative passage back to W. The closest hostile
is still a missing marked atom, not a random long orbit. The least-used
sidecar is the actual next odd integer in the multiplicative carry correction.

Portfolio: ROOT measurement / forced graph lists / algebraic source changes.
Board: **marked integer, fixed-price flow, normalized defect, full-label
coupling, palette compatibility, divisor transport**. The supplied Kakeya
paper suggests retaining weights under restriction; its hypotheses and the
source-specific obstruction are recorded in the [intake](kakeya_marked_context_intake_20261011.md).

Let U(n)=oddpart(3n+1), S(n)=4n+1. A rooted n>1 has strict first-hit
word a_0,...,a_(tau-1), states n_0=n,...,n_tau=1, and counters

    K=sum floor((a_i-1)/2),   L=tau-1,
    f(n)=f_(1/4)(n)=4^-K (3/4)^L,
    W(n)=2 K! tau!/(K+tau+1)!.

Both scores are zero on unrooted positive odd integers and one at1.
This definition assumes no universal convergence.

## Q1. A sharp rational envelope and the exact remaining defect

For every positive odd n,

    0 <= f(n) <= 4/(3n+1).                               (1)

Equality holds precisely at n=1 and the one-step ROOT ray
n=(4^h-1)/3, h>=2. In particular the coefficient and exponent in this
source envelope cannot be improved uniformly.

Put Z(n)=(3n+1)f(n)/4. For a rooted n>1,

    Z(n)=2^-o product_(i=1)^(tau-1) [3n_i/(3n_i+1)],      (2)

where o counts the odd valuations among a_0,...,a_(tau-2).
The terminal valuation is even, so it could equivalently be counted too.
The empty product equals1.

**Proof.** A terminal source has a=2h, K=h-1 and f=4^(1-h),
giving equality in(1). At a nonterminal edge n->y>1 write
j=floor((a-1)/2). Its exact fixed-price multiplier is

    f(n)=(3/4)4^-j f(y).

Using 3n+1=2^a y gives

    Z(n)=epsilon(a) [3y/(3y+1)] Z(y),
    epsilon(a)=1 if a even, and 1/2 if a odd.             (3)

Iteration gives(2). Every nonterminal factor is strictly below1;
unrooted scores are zero. This proves the bound and equality statement.

This is the corrected potential: Z is not a new freely assignable score.
It retains the same ROOT support. A lower bound Z(n)>=eta would bound
both the number of odd-valuation charges and the accumulated reciprocal
carry cost in(2). No bound on that total cost is asserted here.

There is also an unconditional interpretation of the infinite product.
Every nonrooted positive orbit must contain infinitely many valuations1:
if eventually all valuations were >=2, each later integer >1 would strictly
decrease. Therefore it contains infinitely many odd valuations, and the
product of the factors in(3) tends to zero. A nontrivial positive cycle is
covered by the same argument. This proves that terminating at the one-step
ray or taking the infinite product gives exactly Z; it does not prove that
every source takes the terminating alternative.

The normalization also removes the entire sibling coordinate:

    Z(S^j n)=Z(n).                                      (3a)

Indeed 3S^j n+1=4^j(3n+1), cancelling the exact multiplier4^-j of f.
This is an exact quotient preserving ROOT support, not a forgotten
source label: the original integer is still retained for its measurement.
Equations(3), the value1 on the one-step ray, and boundedness uniquely
determine Z. On an infinite orbit their products tend to zero, so the
difference of two bounded solutions tends to zero after iteration.
Existence and uniqueness of this bounded solution do not imply it is
strictly positive at every source; its zeros encode the unresolved case.

## Q2. The new source-specific measurement threshold

For n=2m+1 and j>=1 keep the inherited native modulus

    q=(3m+2)(4^j-1)/3,
    nu(a)=(a+1)^-6-(a+2)^-6,
    M_q^f=sum_(k>=0) nu(m+kq) f(n+2kq).

The first alias is S^j(n), and its multiplier is exactly4^-j.
Set

    A_q^f=nu(m)+4^-j nu(m+q),
    E_q^f=9/(448 q^8).                                  (4)

For any **independently justified** L_q^f<=M_q^f,

    f(n)>=max(0,(L_q^f-E_q^f)/A_q^f).                    (5)

**Proof.** Equation(1) and the integral expression for nu give, for k>=2,

    nu(m+kq) <= 6/(kq)^7,
    f(n+2kq) <= 4/[3n+6kq+1] <= 2/(3kq).

Thus the residual is at most

    4 q^-8 sum_(k>=2) k^-8
      <=4 q^-8 [2^-8 + integral_2^infinity x^-8 dx]
      =9/(448 q^8).

Separating k=0,1 now gives(5). A sharper valid bill retains the sum
sum_(k>=2) 4 nu(m+kq)/(3n+6kq+1), with a certified integral tail.
The numerical constant in(4) is convenient, not claimed optimal.

This replaces the q^-7/(log q)^2 capacity only after changing from M_q^W
to M_q^f. An old lower bound for M_q^W is not automatically a lower bound
for M_q^f. As before, q_(2j)=(4^j+1)q_j gives nested native moduli.
No L_q^f crossing(4) is supplied for an arbitrary unknown source.

## Q3. A genuine lower-bound transfer back to the beta mixture

For 0<r<1 and p>1 define

    C(r,p)=(1-r^p)^2-(1-r)^(2p)>0.

For every positive odd source, with the inherited f_r,

    W(n)>=C(r,p) f_r(n)^p.                              (6)

In particular,

    W(n)>=(9/16) f(n)^2.                               (7)

**Proof.** On the nonempty interval
[r^p, 1-(1-r)^p], both t>=r^p and 1-t>=(1-r)^p.
For a rooted source, t^K(1-t)^L>=f_r(n)^p there. Integrate against
2(1-t)dt; its mass on this interval is exactly C(r,p). Unrooted sources
give zero on both sides and n=1 is immediate. At r=1/4,p=2 the interval
is[1/16,7/16] and its mass is9/16.

Consequently the independently positive readout in(5) would imply

    W(n)>=(9/16) [max(0,(L_q^f-E_q^f)/A_q^f)]^2.        (8)

This lower transport is distinct from the sharp inherited inverse-section
inequality SF1 with coefficient160/441. The two can be composed, retaining
the actual source and the section's native guard.

## Q4. Why the Kakeya information move needs a new arithmetic input

For probabilities P,Q and an event A with P(A)=1,

    D_KL(P||Q)>=log(1/Q(A)).                             (9)

Indeed Q conditioned on A has density Q/Q(A); splitting the logarithm
gives D(P||Q)=D(P||Q(.|A))+log(1/Q(A)). This is the elementary support
test used in the supplied paper's Section8.4. It suggests the following
**OPEN helper**, rather than assuming independence of two Collatz views:

*Construct, without a completed ROOT word at n, a same-source coupling P
and a comparison law Q such that unrootedness implies Q(A)<=delta^c,
P(A)=1, and an independent argument proves D(P||Q)<c log(1/delta).*

Neither the coupling with the required support nor the small information
bound has been constructed here. The paper gets them from geometric
nonconcentration and averaging over offsets, which our fixed integer has
not been shown to supply.

The deletion hostile is exact. If Q has a marked atom of mass p<1, deleting
that atom and renormalizing gives P with D(P||Q)=-log(1-p). Retaining a
constant fraction of the total mass counts as an extensive restriction in
the paper's exponent sense; it can erase the marked atom entirely.
Even a loss tending to zero can erase a marked atom whose mass tends to
zero. Thus a useful Collatz import must preserve the marked row under every
restriction, or explicitly pay its lost mass. Typical-index conclusions
alone do not establish an independently positive ROOT score.

## Q5. An exact hostile to a uniform lower bound with exponent below two

The sharp n^-1 upper envelope should not suggest a comparable uniform
lower envelope. For every integer h>=1 put n_h=2^(h+1)-1. Its first h
odd steps have valuation1 and states

    n_(h,i)=3^i 2^(h+1-i)-1,       0<=i<=h.

The last of these is 2*3^h-1>1. Applying the fixed-price recurrence h
times, followed by Q1, gives unconditionally

    f(n_h) <= (3/4)^h * 4/(6*3^h-2)
            =8/[(n_h+1)^2 (3-3^-h)].                   (10)

If the source were unrooted both it and the endpoint have score zero,
so this argument assumes no unproved termination. In particular, for
every p<2, n_h^p f(n_h) tends to zero. There is **no** constant c>0 with
f(n)>=c n^-p for every positive odd n when p<2. This does not exclude a
source-dependent floor, an exponent at least2, or another choice of weight.
The elementary all-valuation-one orbit formula is not claimed new; its
role here is to calibrate this particular lower-bound target before
attempting an entropy or renewal estimate.

## Reproduction and next decisive probes

[Checker](../../04-computation/experiments/marked_context_transfers_20261011.py)
and [output](marked_context_transfers_20261011.json). All odd sources through
20001 are replayed with a10000 odd-step cap; hitting the cap raises rather
than fabricating a zero score. The first1000 source range has six sibling
checks each. Twenty-four native source/modulus pairs compare forty literal
residual terms with the proven bound for the entire infinite tail.
The all-valuation-one formulas and bound(10) have separate finite controls
for h=1,...,18; their all-h scope follows from the displayed induction.

1. Seek a lower readout for **M_q^f**, not the old mixture readout.
2. Test same-source coupling candidates against marked-atom deletion before
   attempting entropy regularization.
3. Use the defect product to compare actual excursion blocks; any claimed
   budget must control odd-valuation charges, not just total density.

These are narrower obligations than an unspecified positivity proposal,
but none is asserted solved for every positive integer.
