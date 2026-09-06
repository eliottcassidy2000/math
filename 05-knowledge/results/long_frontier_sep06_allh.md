# Exact beta steps and an obstruction to all-channel sign propagation

**Status: PROVED elementary parameter-step identities; REFUTED qualitative
propagation rule; FINITE-EXACT genuine controls; INDEPENDENTLY AUDITED.**
The [independent proof, source and replay audit](long_frontier_sep06_allh_audit.md)
passes all stated identities, both model obstructions, and the bounded
genuine controls. No mathematical correction was required.
The actual all-h doubled-row sign theorem remains **OPEN**. The obstruction
below is explicitly a coefficient model with incorrect interior amplitudes,
not an actual Laurent counterexample. It keeps more than separate
real-rootedness: the genuine first equations, the genuine lower carry and
leading coefficient, and the exact parameter-step operator all survive.

## 1. Inheritance and the precise question

Use the complete A=2,B=3 progression

    h>=1, x>=1 integers, g=x+3h+1,
    f(u)=alpha*u^(-6h-3)+beta*u^(2g-6h-3)+gamma*u^(3g-6h-3).

The coefficients alpha,beta,gamma are nonzero. The designated mass g is
the first support return if gcd(g,6h+3)=1. All coefficient identities
below hold without that separate gcd filter. Put tau=alpha*gamma^2/beta^3
and X=alpha^x beta^(3h) gamma. The exact monic first row and normalized
carried second row are

    P_(h,x)(t)=sum_(j=0)^h
      [(2h+1)! (x+h)_(h-j)/((3h-3j)!(1+2j)!)] t^j,

    Q_(h,x)(t)=sum_(e=-1)^(2h)
      [(2x+2h)_(2h-e)/((6h-3e)!(2+2e)!)] t^e.            (1)

They retain the literal laws

    CT(f^g)=X binom(g,2h+1) P_(h,x)(tau),
    CT(f^(2g))=X^2 (2g)_(4h+2) Q_(h,x)(tau).             (2)

The full count fibres are (x+j,3h-3j,1+2j), j=0,...,h, and
(2x+e,6h-3e,2+2e), e=-1,...,2h. The lower carry is never discarded.

The closest mechanisms are the
[all-h characteristic interface](overnight7_20260906_laurent_midpoint_transport.md),
the [completed-alpha negative response](overnight8_20260906_alpha_completion.md),
and the now-audited [endpoint-39 certificate](long_frontier_sep06_endpoint39.md).
**THM-4436**,
[complete factorial-row simple negative roots](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md),
proves that P_(h,x) has h distinct negative roots. It does not compare
two different designated masses. The coefficient-monotone path comparison
also does not compare values on the negative phase ray.

The canonical hostile is the actual derivative-cone separator in
[the endpoint-33 note](second_20260906_laurent.md); its entire permitted
cone misses a genuine response whose sign is nevertheless negative.
The corrected near miss is the same note's positive phase multiplier
that fails at x=10000. The newer
[sparse-amplitude repair](continuing2_20260906_laurent_sparse_amplitude.md)
requires selected square-root branches and algebraic coefficients at one
actual row. It supplies no uniform propagation theorem.

The retained concept board is: complete carried coefficients; original
first zeros; beta integration in x; separate negative-root geometry;
the actual interior amplitudes; and hit/skip response decomposition.
Targeted searches found no existing statement of the beta-step propagation
rule or the exact obstruction below. No external priority claim is made.

## 2. An exact parameter step for the entire progression

Let Theta=t*d/dt. Direct coefficient ratios in (1) give, for every h,x
in scope,

    (x+1+Theta) P_(h,x+1)=(x+h+1) P_(h,x),

    (2x+1+Theta)(2x+2+Theta) Q_(h,x+1)
      =(2x+2h+1)(2x+2h+2) Q_(h,x).                     (3)

For the first identity the coefficient ratio is
(x+h+1)/(x+j+1). For the second it is
(2x+2h+1)(2x+2h+2)/[(2x+e+1)(2x+e+2)]. These statements include
the leading coefficients and e=-1 lower carry. Every denominator is
strictly positive for x>=1 and the full displayed index range.

Equivalently, coefficientwise integration gives

    P_(h,x+1)(t)=(x+h+1) integral_0^1 s^x P_(h,x)(st) ds,

    Q_(h,x+1)(t)=(2x+2h+1)(2x+2h+2)
      *integral_0^1 s^(2x)(1-s) Q_(h,x)(st) ds.          (4)

The second integral is legitimate even for the carry: its smallest
power of s is 2x-1>=1. Thus the equations preserve the entire Laurent
row, rather than replacing it by its nonnegative-power part.

These identities suggest propagating the sign Q<0 at P's roots from x
to x+1. However the two integration kernels differ. A root of the new P
is a zero of one weighted integral, not a zero of every integrand; the
old pointwise sign at the old P-roots supplies no direct sign to the
new Q integral.

## 3. The qualitative propagation rule fails at h=1

Here the genuine first rows are P_(1,1)=t+2 and P_(1,2)=t+3.
For clarity multiply the genuine Q by the fixed positive number 720.
The complete genuine polynomial rows t*720Q at x=1 and x=2 are

    t^3+20t^2+6t+1/21,
    t^3+30t^2+15t+5/21.                                (5)

Consider instead the model Laurent row M_1 determined by

    t M_1(t)=(t+2/43)(t+1/2)(t+43/21)
      =t^3+(4685/1806)t^2+(2063/1806)t+1/21.             (6)

Apply exactly the second step in (3), with h=x=1, to define M_2.
Writing coefficients of tM in ascending order, the multiplier at index
k is 30/[(k+2)(k+3)]. Thus

    t M_2(t)=t^3+(4685/1204)t^2+(10315/3612)t+5/21.       (7)

The first-row step is also exact and gives t+2 -> t+3. Both model
cubic polynomials have all coefficients strictly positive and three
distinct strictly negative roots. This is immediate from the three
distinct factors in (6). For (7), the intervals

    (-3,-2), (-1,-1/2), (-1/2,0)

each have a strict sign change; degree three exhausts all roots. The
exact endpoint values are included in the source/output certificate.
Independently its discriminant is

    1152093393330275/56737436781312>0.

Nevertheless, at the genuine original first roots,

    M_1(-2)=-3/43<0,
    M_2(-3)=557/5418>0.                                (8)

The lower-carry and leading coefficients in (6)--(7) agree exactly with
(5). Only the two interior coefficients differ. Therefore the following
package is insufficient to propagate negativity: the genuine first
equations and their step, correct carry and leading anchors, strict
negative-root geometry of both doubled polynomial rows, positive
coefficients, the exact doubled parameter step, and the preceding
same-root negative sign.

This is the smallest nonconstant channel count h=1 and the smallest
allowed x-step. At h=0 there is no first cancellation root to propagate.
The model is not substituted into (2) and is not declared an actual
constant-term response. In particular (8) does not refute the proved
actual two-first-channel theorem.

## 4. The obstruction is not caused by a missing first-return gcd

At h=1, the first example's adjacent designated masses are 5 and 6,
and only the first is coprime to 9. A second exact model keeps both
adjacent masses primitive. Take x=3, so g=7 and g+1=8, and set

    t M_3(t)=(t+1/6)(t+100/101)(t+101/25)
      =t^3+(78731/15150)t^2+(73301/15150)t+2/3.

The exact x-step now multiplies coefficient k of tM by
90/[(k+6)(k+7)], yielding

    t M_4(t)=t^3+(78731/12120)t^2+(219903/28280)t+10/7.

The genuine first rows are t+4 and t+5, and the displayed lower carry
and leading coefficients are again genuine after the same factor 720.
Both model cubic rows have three distinct negative roots: the first
is factored, and the second has positive coefficients and discriminant

    859047951135198760589/2467080636572160000>0.

The responses reverse sign:

    M_3(-4)=-874/7575<0,
    M_4(-5)=221/21210>0.

Thus requiring the two designated first masses to be primitive does not
repair the qualitative propagation rule. The missing coordinate remains
the genuine interior factorial amplitudes.

## 5. What survives and what the genuine controls say

Equations (3)--(4) are valid on the actual coefficient family and can
still be useful. The obstruction excludes an invariant based only on
the retained qualitative predicates and carry/top anchors. It does not
exclude a sign invariant restricted to the exact factorial orbit, an
invariant retaining additional interior coefficient ratios, or an
all-h proof by a different representation.

The inherited alpha-completed response is

    W_e=[binom(2g,2+2e)/(2g)_(4h+2)] sum_j B_j B_(e-j),
    B_j=binom(x+3h-2j,x+j), -x<=j<=h,

which uses the same positive normalization as Q. At every first root W<0
is already proved. The residual Q-W is the actual beta-skip response;
its universal nonpositivity would suffice to close the target.

The cheap genuine bank here is exactly h in {7,8}, x in {1,2,10,100}.
For each of these eight rows, full literal charge enumeration reconstructs
P and Q, including the e=-1 carry. The full negative support of B is
retained before convolution. The characteristic polynomials of quotient
multiplication by Q and by Q-W have every coefficient strictly positive.
Since P has only real roots by THM-4436, both responses are strictly
negative at every first root in these eight rows. This gives 60 actual
root phases with both signs certified. Six rows have a primitive first
mass; the other two are coefficient controls only.

The full rational characteristic coefficients are printed in the eight
`GENUINE_CONTROL` output lines. These are fixed-row sign certificates,
not polynomial certificates for an unbounded endpoint parameter and not
an all-h theorem. No endpoint-45 or endpoint-51 table was compiled.

## 6. Reproduction, evidence and stopping reason

[Source](../../04-computation/long_frontier_sep06_allh.py) and
[output](long_frontier_sep06_allh.out) retain the exact models, their
correct and incorrect coefficient coordinates, all sign brackets, both
primitive-mass controls, the symbolic step checks, and the eight genuine
complete-fibre tests. The all-h identities follow from the displayed
coefficient ratios; the symbolic h=1,2,7,8 checks are independent controls,
not an inference from finitely many h.

```bash
python3 -B 04-computation/long_frontier_sep06_allh.py
python3 -B -O 04-computation/long_frontier_sep06_allh.py
```

Both runs pass **142 always-active exact gates**, with byte-identical
frozen output. Raw file hashes:

```text
source SHA256 7cc83c7518c27e70db33936809e119a315dc16e1ceb9573ea40bce939d42b1ad
output SHA256 bb012dd9cf4b00b904d1cf74746e3455dc753db16fb308bd94d67f4e68bdf07e
semantic SHA256 a9d05da8f843b6365851e248308082bfbd5d57a7da0379b91c5bddd3726a891f
```

Connection contract: the source is the exact pair of carried coefficient
rows; the map is the pair of beta integrations in (4); the target is
same-root negativity after the parameter step. The map preserves the
Laurent coefficient law on the actual orbit, but the tested abstraction
forgets the interior amplitudes. The two models prove that root geometry,
anchors and the first equation do not recover that information. The next
useful sidecar must retain more of the interior factorial ratios or provide
a direct coupled inequality. This is the bounded stopping obstruction;
no genuine all-h sign failure has been found.
