# Pell scales, odd-zeta capacity differences, and the information they retain

2026-10-11. **PROVED, scoped and self-audited:** the norm-one Pell family,
its field-changing limit, and the capacity expansion/difference below.
**CITED:** the classical Apéry mechanism. **FINITE-EXACT:** the declared
controls. **OPEN:** new odd-zeta irrationality and independent universal
ROOT positivity. A positive capacity is not a positive actual ROOT readout.

## Inheritance before analogy

[THM-3335, square-triangular Pell/Markov selector](../../01-canon/theorems/THM-3335-square-triangular-pell-markov-pythagorean-selector.md)
already gives the full orbit and its Pythagorean, Markov and tournament
interfaces. Its recurrence is inherited here, not rediscovered. The
[odd-zeta audit](odd_zeta_parallels_petersen_20261001.md) already rejects
the cited three-zeta preprint's denominator gain, and the
[credit/zeta bridge](collatz_credit_zeta_bridge_20261005.md) already isolates
the missing arithmetic lattice for the ROOT readout. Neither gap is removed
by a matching recurrence root.

The nearest mechanism is an exact weighted tail bill. The hostile is a
positive smoothing kernel with a missing marked atom. The corrected near
miss is promoting rational coefficients of an integrand to integrality of
its infinite expectation. The useful sidecars are the actual measure,
denominators, and the field containing the limiting point.

## P1. A literal square-triangular family of unit directions

Start with (x_0,q_0)=(1,0) and use the inherited recurrence

    (x_(k+1),q_(k+1))=(3x_k+8q_k,x_k+3q_k).

Then x_k^2-8q_k^2=1 and T_((x_k-1)/2)=q_k^2. With w=sqrt(-55), define

    u_k=(x_k+q_k w)/(x_k-q_k w)
       =[1-47q_k^2+2x_k q_k w]/[1+63q_k^2].              (1)

Every u_k lies in Q(sqrt(-55)) and has norm one. This follows by multiplying
the conjugate factors; the exact identity is

    (1-47q^2)^2+55(2xq)^2=(1+63q^2)^2
    whenever x^2-8q^2=1.

Here **63=8+55** is produced by combining the two norm forms. It is not a
claim that the Mersenne exception at 63 causes this identity. The first row
(x,q)=(3,1) gives u_1=(-23+3w)/32=bar(eta)^2, with eta=(-3+w)/8.

Let lambda=3+2sqrt2. The Pell orbit has x_k/q_k -> 2sqrt2 and hence

    u_k -> u_infinity=(-47+4sqrt(-110))/63.              (2)

The error decays on the scale lambda^(-2k), since
x_k/q_k-2sqrt2=1/[q_k^2(x_k/q_k+2sqrt2)]. All finite directions in (1)
belong to Q(zeta55), but the limit in (2) does not: its quadratic field
Q(sqrt(-110)) is not one of that cyclotomic field's three quadratic
subfields Q(sqrt5), Q(sqrt(-11)), Q(sqrt(-55)).

This is a precise warning about a fractal or limiting construction: metric
convergence preserves the norm equation, but need not preserve the arithmetic
carrier. It also does not preserve injectivity or a finite graph's nonedges.
One must retain exact finite coordinates or explicitly change the field.
The [degree-eight ceiling](rhombus_overlap_elliptic_20261011.md) applies to
each finite direction in (1), so this family by itself supplies no six-colour
forcing context.

## P2. The Pell connection to Apéry has a precise scope

The square-triangular recurrence has characteristic roots 3+/-2sqrt2.
Their squares are 17+/-12sqrt2, the roots of X^2-34X+1. These are the
limiting characteristic roots of Apéry's recurrence for zeta(3).
[Brown, sections 1.2 and 10.2.2](https://www.ihes.fr/~brown/IrratModuliMotivesv8.pdf)
records the recurrence and the successful balance between decay
(sqrt2-1)^4 and denominator growth. It also separates smallness,
nonvanishing, denominator control, and unwanted-zeta coefficient vanishing.

The identity of these quadratic units explains a shared growth scale. It
does not identify the solutions of the constant-coefficient Pell recurrence
with solutions of Apéry's polynomial-coefficient recurrence. In particular,
it gives no rational approximation to zeta(7). The original recurrence,
initial data, represented linear form, and arithmetic denominator bound are
the missing sidecars in such a proposed transfer.

## P3. A positive two-scale difference isolates an odd-zeta coefficient

For a fixed integer m>=0 and positive integer q, take the sharp universal
envelope of the quarter-price aliases after the known first alias:

    B_m(q)=sum_(k>=2) g(m+kq),
    g(x)=2[(x+1)^-6-(x+2)^-6]/(3x+2).                   (3)

It bounds the actual unknown tail because f(2x+1)<=2/(3x+2). This B is a
capacity computed from the envelope; it is not the actual weighted measure.

There is an exact positive-integral representation

    g(x)=28 integral_(1<=t<=2) integral_(0<=v<=1)
               v^6 [x+2/3+v(t-2/3)]^-8 dv dt.          (4)

Indeed the difference of sixth powers is 6 integral_1^2(x+t)^-7 dt;
the elementary parameter identity for A^-7 B^-1 gives (4). If desired,
that identity follows by substitution in the derivative of
v^7/[B+(A-B)v]^7, whose derivative is
7B v^6/[B+(A-B)v]^8.

Taylor's theorem applied to (1+z)^-8 for z>=0 gives alternating finite
upper and lower bounds. With

    C_m=32m+134/3,
    D_m=144m^2+402m+2620/9,

one obtains, for every positive q,

    4(zeta(8)-1)/q^8 - C_m(zeta(9)-1)/q^9
       <= B_m(q)
       <= 4(zeta(8)-1)/q^8 - C_m(zeta(9)-1)/q^9
          + D_m(zeta(10)-1)/q^10.                      (5)

The lower bound can be negative at small q; the claim is the inequality,
not positivity of that truncation. The coefficients are moments of the
bounded shift a=m+2/3+v(t-2/3), weighted by 28v^6. Higher alternating
enclosures follow by the same formula, without an infinite-series exchange.

Now define

    Delta_m(q)=256 B_m(2q)-B_m(q).

It is **strictly positive for every q**: in (4), multiplying the 2q term
by 2^8 replaces (kq+a)^-8 by (kq+a/2)^-8, which is larger because a>0.
The leading even-zeta coefficient cancels, leaving the exact asymptotic

    q^9 Delta_m(q) -> (16m+67/3)(zeta(9)-1).             (6)

The remainder bounds in (5) justify this limit directly. Thus an odd zeta
value arises from a **difference between two scales of the new eighth-power
capacity**, whereas the older unweighted sixth-power prior produced zeta(7)
at seventh order. This explains the change of zeta index instead of treating
all occurrences of odd zeta as interchangeable.

More generally, for any fixed c>1 the difference
c^8 B_m(cq)-B_m(q) has the same sign and a next-order coefficient multiplied
by 1-1/c. Integer c keeps the original integer moduli. Pell-generated pairs
of integer scales can use the exact ratio of their eighth powers, but their
ratio tending to lambda adds no denominator or nonvanishing theorem.

## P4. The two gates this construction does not cross

First, B_m(q) is an infinite sum. Rational summands and exact rational
moment coefficients do not make B_m(q), Delta_m(q), or a denominator-cleared
version rational. Therefore (6) is not an irrationality proof. The useful
extension question is whether an independently constructed linear-form
carrier can reproduce its cancellation while paying **all** coefficient
denominators and eliminating unwanted constants. The inherited failed
odd-zeta audit makes this a mandatory arithmetic calculation, not a heuristic
bonus from rapid decay.

Second, the positivity of Delta_m is specific to the envelope kernel (4).
It does not imply positivity of 256 M_(2q)-M_q for the actual ROOT measure.
A measure with its only mass at the odd alias m+q has M_q=1 and M_(2q)=0,
so that signed difference equals -1. This is a hostile for transferring a
kernel property to an arbitrary law; it is not claimed to be the Collatz
law. A Collatz-specific coupling or lower bound must be proved separately.

The fixed-source progress is instead the [charge-floor theorem](collatz_shadow_charge_20261011.md):
bounded odd charge pays the carry loss and then the marked readout. The
capacity expansion reduces an error bill; it cannot pay the missing source.

## Connection ledger and reproduction

| Source -> target | Map / preserved quantity | Lost information and decisive test |
|---|---|---|
| Square-triangular Pell point -> unit direction | Formula (1), norm one and exact coordinates | Integer size becomes a denominator bill; verify the quadratic norm |
| Pell sequence -> its metric limit | Formula (2), unit length | The field changes; test quadratic subfield membership |
| Pell unit -> Apéry limiting roots | Square the two characteristic roots | Full recurrence and linear form are lost; retain polynomial coefficients and denominators |
| Quarter-price envelope -> odd-zeta scale difference | Equations (4)--(6), sign and cancellation | Actual ROOT weights are replaced by a capacity; the odd-alias point mass defeats that transfer |
| Changing ROOT controllers -> requested source | Shared bits plus bounded total odd charge | Unbounded charge defeats coverage; use the finite precision receipt in C4 |

The [checker](../../04-computation/experiments/changing_controllers_rhombus_20261011.py)
verifies twelve Pell rows, exact norm identities, the first three rational
moment coefficients, and pointwise Taylor inequalities at the declared
finite shifts and scales. The infinite expansion bounds and limit are
proved above. The [saved output](changing_controllers_rhombus_20261011.json)
also retains the unrelated geometry and ROOT controls without blending their
statuses. Universal Collatz, six-colour forcing, and new odd-zeta
irrationality remain OPEN.
