# Trinomial channel evolution: strict local interlacing, an actual-mass hostile, and the remaining sign test

Status is separated by claim. The exact Euler identities, their strict
interlacing consequence, and the real-gauge lemma are **PROVED / INDEPENDENTLY
AUDITED** separately by root and observer. The displayed noninterlacing witness is an
**EXACT REFUTATION** of a global transfer. The higher-multiplicity two-rung
sign and adjacent-level coprimality boxes are **FINITE-EXACT**, not proofs
of unbounded assertions. No theorem identifier or navigation is changed.

## 1. Inheritance and the faithful evolving object

The closest proved mechanism is
[THM-4436, complete factorial-row simple negative roots and trinomial phase wall](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md).
It places each individual return polynomial in the simple-negative-root
locus but deliberately does not compare different masses. The complete
two-first-channel sign statement is already
[THM-4432, trinomial two-channel two-rung noncancellation with carries](../../01-canon/theorems/THM-4432-trinomial-two-channel-two-rung-noncancellation-with-carries.md).
The earlier sidecar
[THM-2760, exact-prefix even Faber flux gcd and smooth-boundary exclusion](../../01-canon/theorems/THM-2760-exact-prefix-even-faber-flux-gcd-and-smooth-boundary-exclusion.md)
uses a derivative response to turn simple roots into genuine coprimality.
Here an elementary factorial shift supplies exactly such a response.

The live board is: finite simple-negative roots; coefficientwise Euler
response; complete carry coordinates; first versus later support fibres;
real coefficient gauge; constant-term square versus positive norm.
The map from a factorial parameter step to `a+b*t*d/dt` preserves the
whole polynomial and its zero set. A mass step generally combines several
such moves and a change of left monomial prefactor. Forgetting those
coordinates destroys the claimed interlacing, as the explicit witness
below demonstrates.

Targeted searches for Euler/contiguous/factorial interlacing,
nearest-coordinate motion and the exact hostile root word found the
mechanisms above, not an existing all-mass interlacing theorem. The local
identities are elementary consequences of factorial ratios; no claim of
external priority is made.

## 2. Four exact local response identities

Take integers

```text
0<A<B, h>=0, x>=0, 0<=r<B, 0<=z<A,
delta=B-A, m=x+B*h+r+z,
F_(x,h,r,z)(t)=sum_(j=0)^h
 t^j/[(x+delta*j)!(B*h+r-B*j)!(z+A*j)!].             (1)
```

This is the THM-4436 polynomial divided by the positive constant `m!`.
Write `theta=t*d/dt`. Direct coefficient comparison gives

```text
F_(x,h,r,z)=(x+1+delta*theta) F_(x+1,h,r,z),                         (2)

F_(x,h,r,z)=(z+1+A*theta) F_(x,h,r,z+1),       if z+1<A,             (3)

F_(x,h,r,z)=(B*h+r+1-B*theta) F_(x,h,r+1,z),   if r+1<B,             (4)

F_(x,h,B-1,z)=B*(h+1-theta) F_(x,h+1,0,z).                         (5)
```

All functions displayed remain inside the exact canonical factorial
domain of THM-4436. In `(5)` the leading coefficient on the right is
annihilated, not discarded: the new row really has degree `h+1`.
If the `m!`-normalized moment polynomials are used instead, the left side
of each identity is multiplied by `m+1`; the checker retains that factor.

### Strict root order

Assume `h>=1`, let `P` be the old polynomial in each identity and `Q`
the new one, and denote their roots in increasing order by `p_i,q_i`.
Then

```text
(2),(3):  q_1<p_1<q_2<p_2<...<q_h<p_h<0;
(4):      p_1<q_1<p_2<q_2<...<p_h<q_h<0;
(5):      q_1<p_1<q_2<...<p_h<q_(h+1)<0.            (6)
```

Thus each listed neighbouring pair is coprime. Increasing the left
increasing-count offset `x`, or the residue `z` before its wall, moves
roots strictly left. Increasing the decreasing-count residue `r` before
its wall moves them strictly right. Its carry creates one extra root and
gives the third strict ordering.

Proof: THM-4436 gives simple negative roots for `Q`. At each root,
`P(q_i)=b*q_i*Q'(q_i)` where the slope `b` is respectively
`delta,A,-B,-B`; it never vanishes. In the equal-degree cases the leading
and constant coefficients of `P` are positive. The signs at consecutive
`q_i` alternate. For positive `b`, the final root lies in `(q_h,0)` and
the other `h-1` roots lie between adjacent `q_i`; for negative `b`, the
remaining root lies left of `q_1`. In `(5)`, `P` has degree `h` and `Q`
has degree `h+1`, so the `h` sign changes between consecutive `q_i`
already account for every root of `P`. This proves all strict inequalities.
For `h=0`, the first three statements have no roots to compare; the carry
case creates one simple negative root.

These are precise local statements, not a license to compose arbitrary
root orders. In particular `(2)` and `(4)` push in opposite directions.
No assertion is made here about crossing the `z=A-1` wall by `(3)`.

## 3. Consecutive admissible masses need not interlace

For the primitive support `(-3,1,5)`, the return modulus is `g=4`, and
every positive multiple of four is admissible. The monomial invariant is
`tau=alpha*gamma/beta^2`. At masses eight and twelve the complete scalar
polynomials in that **same** invariant are

```text
P(t)=56t^3+420t^2+280t+28,
Q(t)=3960t^4+18480t^3+16632t^2+3960t+220.             (7)
```

Their gcd over the rationals is a nonzero constant. Both have simple
negative roots, yet the exact increasing root order is

```text
P Q Q P Q P Q.                                      (8)
```

A small independent certificate consists of seven disjoint rational
intervals, in that order:

| Polynomial | Interval | Values at the left and right endpoints |
|---|---|---|
| P | `(-7,-6)` | `-560,1372` |
| Q | `(-4,-3)` | `81532,-40172` |
| Q | `(-4/5,-3/4)` | `-17908/125,1991/32` |
| P | `(-2/3,-1/2)` | `308/27,-14` |
| Q | `(-1/3,-1/4)` | `1012/9,-121/32` |
| P | `(-1/8,-1/9)` | `-35/64,1456/729` |
| Q | `(-1/12,-1/13)` | `-1441/288,157828/28561` |

Each row forces a real root by a sign change. The three `P` intervals
and four `Q` intervals exhaust the respective degrees, so this certificate
alone proves completeness, simplicity and the order `(8)`. It does not
need a numerical root finder or a generic interlacing theorem.

The failed implication is therefore

```text
individual simple-negative roots + consecutive admissible masses
  => interlacing                         FALSE.
```

This is not a counterexample to coprimality or to the sharp moment bound.
It is also not claimed minimal among every possible support or ordering.
The source-to-target loss is the actual sequence of residue carries and
offset increments. A separate small example, `(-1,1,2)` at masses twelve
and thirteen, even has scalar degrees `2` and `1`: normalized row degree
need not increase with actual mass.

## 4. Exact higher-multiplicity first/second sign signal

The following remains **FINITE-EXACT evidence**, except for the inherited
`h=1` theorem THM-4432. Let

```text
g=x+B*h+r+z,
a=A*(B*h+r)+B*z,
b=g*A-a>0, c=g*B-a,
gcd(A,B)=gcd(a,g)=1.                                (9)
```

Then the first nonempty support fibre is at mass `g`: every support
return mass is divisible by `g`, and the explicitly nonempty row `(1)`
occurs there. Its complete channels are

```text
(x+delta*j,B*h+r-B*j,z+A*j), 0<=j<=h.
```

Let `P` be its multinomially normalized polynomial. Set

```text
epsilon_y=floor(2r/B), epsilon_z=floor(2z/A),
h_2=2h+epsilon_y+epsilon_z,
r_2=2r-B*epsilon_y, z_2=2z-A*epsilon_z,
x_2=2x-delta*epsilon_z.                             (10)
```

The complete second polynomial `Q` is the same factorial construction
with these new coordinates. Its left monomial prefactor differs from the
square of the first one. If `M=alpha^x*beta^(B*h+r)*gamma^z`, then exactly

```text
CT(f^g)=M*P(tau),
CT(f^(2g))/M^2=tau^(-epsilon_z)*Q(tau).              (11)
```

The left carry sign in `(11)` is essential. At a negative root of `P`,
the observed second normalized sign is always **strictly negative**, not
the uncorrected sign of `Q`.

The full tested universe is `2<=B<=6`, `1<=A<B`, `gcd(A,B)=1`,
`h=2,3,4`, every canonical `r,z`, `x=1,2,3`, filtered exactly by `(9)`.
There are 638 rows: 213 of degree two, 210 of degree three, 215 of degree
four. Carry counts `(epsilon_y,epsilon_z)` are

```text
(0,0):224, (0,1):118, (1,0):194, (1,1):102.
```

For every row the checker independently reconstructs both full weighted
fibres from the original charge equations. It isolates every first root
over the rationals, divides `Q` by `P`, and refines an interval evaluation
of the exact remainder until its sign is certified. Thus no rounded root
or approximate near-zero value supplies the sign claim.

The tempting unbounded statement is

```text
P(t)=0, t<0 => t^(-epsilon_z) Q(t)<0               for every (9).   (12)
```

It is **OPEN here** for arbitrary `h>=2`. If proved, `(12)` would close
the two-rung obstruction for every collided trinomial first fibre and
give the sharp ceiling `2g<=a+c=gB`. That consequence is why evidence
must not be promoted merely from a large finite box.

### Large-offset hostile-direction controls

A separate test takes `A=1,B=2,r=z=0`,
`h in {3,4,5,6,8,10,12}`, and
`x in {5,7,11,101,503,997}`, filtered by `gcd(2h,x)=1`.
The resulting 40 primitive supports are
`(-2h,x,2x+2h)` with first mass `g=x+2h`.
All first roots again give strictly negative second values, exactly.
To avoid enormous irrelevant common factors, the checker uses

```text
P(t)/P(0)=sum_j (2h)! t^j/[j!(2h-2j)!(x+1)_j]
```

and its doubled row with `(h,x)` replaced by `(2h,2x)`.
This is a stress test along an asymptotic direction, not an asymptotic
proof. In particular the factor `2x` must remain when comparing limiting
polynomials: an unscaled same-argument Hermite comparison would describe
the wrong second row.

## 5. A real coefficient gauge exists, but the square is not a positive norm

For any trinomial coefficient torus point with `tau<0`, its orbit under
common coefficient scaling and Laurent-variable scaling has a real
representative with exactly one negative coefficient.

To see this without assuming real root choices, first normalize
`alpha=beta=1`: choose `sigma^(a+b)=alpha/beta` and common scale
`lambda=sigma^a/alpha`. The new `gamma` satisfies `gamma_new^A=tau`.
As the `g*A` possible choices of `sigma` vary, `gamma_new` is multiplied
by all `A`th roots of unity, because `gcd(A,B)=1`. Thus when `A` is odd
we may choose `gamma_new<0` real. If `A` is even then `B` is odd;
normalize `alpha=gamma=1` instead. The same argument gives
`beta_new^B=1/tau` and allows `beta_new<0` real. This proves the claim
without requiring primitive original charges.

Both `tau` and the normalized quotient in `(11)` are unchanged by these
operations. In a real representative, `M^2>0`, so the proposed negative
sign `(12)` is exactly the sign of the real constant term `CT(f^(2g))`
at a zero of `CT(f^g)`.

But for a real Laurent polynomial `q`,

```text
CT(q(u)^2)=sum_j q_j*q_(-j)
```

is an indefinite bilinear form. The positive norm is instead
`CT(q(u)*q(u^(-1)))=sum_j q_j^2`. Mean-zero examples already separate
the two signs:

```text
CT((u+u^(-1))^2)=2,
CT((u-u^(-1))^2)=-2.
```

These are not trinomial first-fibre counterexamples; they identify the
first invalid norm inference. A successful proof of `(12)` must keep the
specific multinomial fibre constraints or an actual involution that turns
them into a signed norm. The negative phase alone does not supply that.

### Individual opposite-frequency products also fail to have a fixed sign

This failure already occurs inside an actual collided first fibre. Take

```text
f=u^(-4)+u+t*u^6,       q=f^5,
P(t)=CT(q)=10t^2+30t+5.
```

Let `q_e` denote the coefficient of `u^e` in `q`. Modulo `P`, the exact
products at opposite frequencies are

| e | Remainder of `q_e*q_(-e)` modulo P | Signs at the two increasing roots of P |
|---|---|---|
| 5 | `3220t+560` | negative, negative |
| 10 | `0` | zero, zero |
| 15 | `700t+125` | negative, positive |
| 20 | `-35t-25/4` | positive, negative |

These are all nonzero-frequency pairs contributing to `CT(q^2)`.
Their convolution identity is checked literally against the complete
mass-ten moment. Thus the proposed implication from the first moment's
zero to termwise nonpositivity is **REFUTED**, even in this small
three-channel first row. The exact stronger survivor is grouping:
the last two remainders are negatives of one another with factor twenty.
Keeping only the signs loses that decisive coefficient ratio.

## 6. Adjacent-level gcd controls and reproduction

A separate universe takes every primitive support
`1<=a<=8`, `1<=b<=5`, `b<c<=10`, and every mass `1<=m<=40`.
After retaining consecutive **nonempty** support returns with both scalar
degrees positive, it contains 247 supports and 2,759 pairs. Every pair
has constant gcd over the rationals. The degree-jump counts are

```text
-1:211, 0:1909, 1:589, 2:34, 3:11, 4:5.
```

This is a finite control bank, not a universal coprimality theorem. It
includes later levels and carry onsets, unlike the first/second bank.

The local identity/interlacing controls use every `2<=B<=4`, `1<=A<B`,
`h=1,2,3`, all canonical residues and `x=0,2`, without a gcd filter.
They audit 210 x-steps, 90 z-steps, 150 interior r-steps and 60 r-carries.
Root labels are endpoint-safe: an isolating interval can share a rational
endpoint with a different root, so inclusive endpoint counts are corrected
before assigning the labels. A separate rational sign-bracket certificate
checks the actual-mass hostile independently of that labelling routine.

```bash
python3 04-computation/nc2_channel_contiguous_overnight_hexagon_sep05.py
python3 -O 04-computation/nc2_channel_contiguous_overnight_hexagon_sep05.py
```

The script imports only standard-library arithmetic and SymPy, not another
repository producer. All failure gates remain active under optimization.
All 23,045 gates pass; normal and optimized runs are byte-identical to
the stored output. The recorded exact-polynomial backend is SymPy 1.9.
Raw-LF SHA-256:

```text
ea9a138929984480539495ddae7856dcb0d5236bb160ebf425272778dab7a6d6  04-computation/nc2_channel_contiguous_overnight_hexagon_sep05.py
b25aaf8b8f5137a1788b119a92aa3190f7bd7ea17bec72be4a9d3000bbeedcd3  05-knowledge/results/nc2_channel_contiguous_overnight_hexagon_sep05.out
```

Root and observer separately audited all four Euler identities, the
normalization, carry degree change, strict root orders, coefficient-gauge
root choices and quotient invariance. Root also checked the rational
bracket degree count and the carry-normalized two-rung formula. These
proof audits do not promote the finite cross-mass sign or gcd statements.
No new literature claim is needed for this continuation: the only
analytic inheritance is the proved THM-4436 simple-root statement.
