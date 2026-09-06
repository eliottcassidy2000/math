# The smallest model phase closes and the uniform tail falls to 75000

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [independent referee](continuing3_20260906_laurent_finite_phase_audit.md)
accepts both uniform theorems and the exact hostile scopes without repair.
Its separate engine passes 663 normal/optimized exact gates and also proves
the full box margin by separate convexity and eight univariate Sturm chains. In the two-anchor model, the entire smallest first-phase branch has
negative full carried response. This first-branch theorem needs no interlacer.
With both interlacers, the inherited fourth-coefficient floor improves the
uniform negative tail from 118163898523 to **75000**. Full model negativity
and general actual Laurent noncancellation remain **OPEN**.

The remaining geometry cannot be replaced by coefficient inequalities and
first-row real-rootedness: an exact hostile below satisfies both anchors,
all four displayed Newton inequalities, and four positive simple first
phases, but its full response is positive at one phase. Its beta polynomial
has only one real root, so this is not a counterexample to the model.

## 1. Inheritance, model, and retained information

The closest proved mechanisms are the
[qualitative quadratic-anchor theorem](open_frontier_sep06_quadratic_anchor.md)
and its audited
[effective fourth-coefficient floor](continuing2_20260906_effective_anchor.md).
The latter proves e4>1/100 and an explicit uniform tail starting at
118163898523. The canonical hostile is the
[linearly anchored escaping-phase model](open_frontier_sep06_laurent.md).
The corrected near miss is treating negative numerical optimization as a
uniform sign proof. The underused sidecar here is the exact positive-monomial
inventory after eliminating e5 by the original first equation.

The board is the two coefficient anchors, the original first zero, both
beta interlacers, the eliminated quadratic response, whole-box Bernstein
coefficients, and the beta-root predicate lost by coefficient relaxation.
The actual fixed-row
[four-window amplitude repair](continuing2_20260906_laurent_sparse_amplitude.md)
provides a useful contrasting mechanism: it has common positive algebraic
coefficients at one actual row, while the present bounds are uniform over
the specified model. No theorem transfers merely from both using root values.

Let a_1,...,a_5 be nonnegative, with

    sum a_i=13,  sum a_i^2=59,  equivalently e1=13, e2=55.

Write x=e3, y=e4, z=e5 and keep

    B(v)=v^5-13v^4+55v^3-xv^2+yv-z,
    C(v)=v^4-12v^3+45v^2-(2/3)xv+(3/7)y,
    D(v)=v^4-11v^3+36v^2-(5/12)xv+(1/7)y.

The full model K additionally requires C and D to weakly interlace B.
Zero and repeated beta roots are allowed. In the same original variable t,
define coefficient arrays and Laurent rows

    G_B=(1,13,55,x,y,z),
    G_C=(1,12,45,(2/3)x,(3/7)y),
    G_D=(1,11,36,(5/12)x,(1/7)y),
    beta=t^-1 G_B, C_raw=t^-1 G_C, D_raw=t^-1 G_D,
    O=sum_j binom(14,2j+1)t^j, E=sum_j binom(14,2j)t^j,
    P=O star beta,
    Q=(O^2+t^-1 E^2) star (beta^2+2t C_raw D_raw).          (1)

Here star is coefficientwise multiplication. All original shifts, binomial
coefficients, and the relative crossing factor 2t remain fixed. In particular
Q has the genuine exponent -1 term 28/t. For s>0,

    P(-s)=182-20020s+2002x s^2-3432y s^3+2002z s^4,
    T(s)=sQ(-s).                                          (2)

The source-to-target map in the first theorem replaces the beta shape by a
containing coefficient box; it preserves (1) exactly but loses real-rootedness
and interlacing. Its stronger whole-box sign certificate makes that loss
harmless only on the stated interval. The tail instead eliminates z using
the exact equality P(-s)=0, retaining this equality and the inherited floor
y>1/100. The last hostile identifies why neither relaxation closes all phases.
No new literature theorem is invoked.

## 2. Elementary bounds from the two anchors

Every nonnegative five-root shape with the anchors satisfies

    0<x<123<130,  0<=y<152,  0<=z<72,  4y<=13x.           (3)

For completeness these bounds can be obtained without importing an extremal
classification. If M is the largest root, Cauchy on the other four gives
5M^2-26M-67<=0. Since M>=13/5 and the left side is increasing thereafter,
its value 9/20 at M=71/10 gives M<71/10. The elementary third-power identity
is x=(sum a_i^3-52)/3, so

    x < ((71/10)*59-52)/3 =1223/10 <123.

There must be at least three positive roots, because otherwise
13^2<=2*59, which is false; hence x>0.

An explicit sum of squares supplies y<=55^2/20=605/4<152. Put

    S=sum_(i<j) (a_i-a_j)^2 e2(a with i,j omitted),
    W=sum_(i<j<k<l) [(a_i a_j-a_k a_l)^2
                    +(a_i a_k-a_j a_l)^2
                    +(a_i a_l-a_j a_k)^2].

Direct expansion gives e2^2-20e4=S+W/3>=0. AM-GM on the ten pair
products gives z^2<=(55/10)^5<72^2. Finally

    e1 e3-4e4=sum_i a_i^2 e2(a with i omitted)>=0.

The source independently checks both polynomial identities in five formal
variables; the inequalities use only nonnegative roots.

## 3. A uniform closed interval contains and settles the smallest phase

**Theorem A.** For every nonnegative five-root shape with the two anchors,
P(-s) has exactly one positive root below 1/90. That root is simple and lies
in (1/110,1/90). At every s in the larger closed interval [1/120,1/80],

    s Q(-s)<-400.                                        (4)

In particular the complete response is strictly negative at the smallest
first phase, without either C or D interlacing B.

For 0<s<=1/110 the constant/linear part of P(-s) is nonnegative. Using
4y<=13x, its next two terms satisfy

    2002x s^2-3432y s^3 >= x s^2(2002-11154s)>0.

The last term is nonnegative, so there is no earlier positive root. At
s=1/90, the box bounds in (3), enlarged to x<=130, give

    P(-1/90)<=182-20020/90+2002*130/90^2+2002*72/90^4<0.

Also the derivative of P(-s), viewed as a polynomial in s, is bounded above
on [0,1/90] by

    -20020+4004*130/90+8008*72/90^3<0.

This proves the existence, uniqueness, and simplicity assertion.

The interval sign (4) is proved on the entire larger box

    (x,y,z,s) in [0,130] x [0,152] x [0,72] x [1/120,1/80].

Substitute x=130X, y=152Y, z=72Z and s=1/120+U/240 in T(s). The resulting
polynomial has coordinate degrees (2,2,2,9). Its **270** tensor-product
Bernstein coefficients are all strictly less than -400. The exact largest is

    -22645374245632441/52254720000000 < -400.               (5)

The Bernstein basis functions are nonnegative on [0,1]^4 and sum to one;
therefore (5) proves (4) on the whole continuous box. This is an identity
certificate, not a sample of parameter values. Explicitly, if c_k are the
power coefficients and d=(2,2,2,9), the certified Bernstein coefficient at
index i is

    b_i=sum_(k<=i) c_k product_j binom(i_j,k_j)/binom(d_j,k_j).

The frozen certificate contains every b_i. The producer also reconstructs
the entire power polynomial by expanding these Bernstein basis functions,
checking the reverse direction of the basis change exactly.

## 4. Eliminating the original zero makes the tail effective at 75000

At an original positive phase, P(-s)=0 gives, with u=1/s,

    z=(12/7)y u-xu^2+10u^3-u^4/11.                       (6)

After substituting (6), write R=Q(-s)/s^7=sum_(j=0)^8 R_j u^j.
The exact full carried coefficients are

| j | R_j |
|---|---|
| 0 | -(26075790/7)y^2 |
| 1 | (153780300/7)xy-(179344800/7)y^2 |
| 2 | -16900975x^2+(647843760/7)xy-1329865290y |
| 3 | -53986980x^2+1122025905x-4282905900y |
| 4 | 3467704710x+(5521932000/11)y-10070260200 |
| 5 | -(3690469830/11)x-9902880y-30313505040 |
| 6 | 6175260x+3654364350 |
| 7 | -1022439600/11 |
| 8 | 565082 |

**Theorem B.** On K, at every positive original phase s>=75000,

    Q(-s)/(s^7 e4^2)<-120000<0.                          (7)

This includes multiple first roots. It inherits y=e4>1/100 from the proved
effective-anchor theorem, and uses x<123 from (3). No interlacer inequality
is silently substituted for the original first equation.

To prove (7), keep the negative constant term and discard every other
negative monomial in the displayed table. Divide by y^2, then use
x/y<12300, x/y^2<1230000, 1/y<100, and 1/y^2<10000. Thus R/y^2 is bounded
above by the following one-variable envelope, with every nonconstant
coefficient strictly positive:

    F(u)= -26075790/7
       +12300[(153780300/7)u+(647843760/7)u^2]
       +1230000[1122025905u^3+3467704710u^4+6175260u^6]
       +(552193200000/11)u^4
       +10000[3654364350u^6+565082u^8].                  (8)

Direct rational evaluation gives

    F(1/75000)
      =-470440303496224131940323407036971853244343
        /3854347229003906250000000000000000000
      <-120000.

Since F increases with u>0, (7) follows for all s>=75000. The source
checks that every positive monomial of the full eliminated row, and no
unrecorded term, contributes to (8). This is a positive-term domination
argument after exact cancellation, not numerical extrapolation of a tail.

## 5. Exact hostiles and the remaining missing predicate

Several separate controls keep the quantifiers honest.

The rational beta roots (1,3,9,22,30)/5 satisfy both anchors and both
strict interlacings, giving a nonvacuous point of K. Its full response at
s=75000 is positive, while P(-75000) is nonzero. Thus the original-root
condition is essential in (7); the theorem does not assert negativity at
every large phase independent of P.

The coefficient-box point (x,y,z,s)=(130,0,0,1/60) has

    T(s)=34734566083/93312000>0.

It shows why the whole-box interval proof cannot simply be extended to
1/60. It is not a beta-root shape and is not a model counterexample.

A sharper hostile survives the original first equation and first-row
real-rootedness. Take

    (x,y,z)=(104,50,37435088/3898125),  s=15/2.             (9)

The first polynomial factors exactly as

    P(-s)=(26/50625)(2s-15)
           *(18717544s^3-26680920s^2+2595600s-23625).

Its other three roots have sign-changing brackets
(1/99,1/98), (3/32,5/53), and (70/53,37/28). These, together with 15/2,
exhaust degree four and give four simple positive first phases. At the
displayed root, however,

    Q(-15/2)=78541969368658673/18480>0.                   (10)

All four usual Newton coefficient inequalities hold strictly; their
explicit margins are

    13^2-(5/2)*55=63/2,
    55^2-2*13*x=321,
    x^2-2*55*y=5316,
    y^2-(5/2)*x*z=2437924/779625.

But its beta polynomial B has exactly one positive real root by exact
Sturm counting. It has no negative real root, since B(-v) is strictly
negative for v>=0. Thus it has four nonreal roots and is outside K.
Neither the positive first-phase roots nor these four coefficient
inequalities restore the lost beta-root predicate. The missing coordinate
is the complete beta-root/interlacer geometry, not another untracked carry.

The first-root elimination also leaves a quadratic in x,y whose Hessian
determinant is negative for every u>0; unconstrained concavity cannot replace
that geometry. The exact hostile (9) already gives a decisive failure of
the stronger proposed coefficient-only surrogate, so no broad census is
needed to establish this stopping boundary.

Combining Theorems A and B, **any still-unresolved original phase in K lies
in (1/80,75000)**. There are at most three such phases per shape, counted
with multiplicity, since the unique simple smallest phase has been settled.
The sign on this remaining interval remains OPEN. A small initial numerical
optimization found no positive response in K, but none of its numerical
maxima or apparent shared boundary points is used as proof.

## 6. Exact scope, reproduction, and artifacts

The [standalone source](../../04-computation/continuing3_20260906_laurent_finite_phase.py)
uses only rational arithmetic and the Python standard library. Its explicit
finite universe consists of the formal five-variable anchor identities;
the full Laurent convolution and all original carried coefficients; every
coefficient of the formal elimination; all 270 Bernstein coefficients; the
complete positive-term envelope; one strict rational model; and the two
precise relaxation hostiles above. Analytic inequalities and the Bernstein
identity supply the uniform quantifiers, rather than these control points.

The [JSON certificate](../../04-computation/continuing3_20260906_laurent_finite_phase_certificate.json)
stores the full T polynomial, eliminated row, 270 coefficients, and exact
tail envelope. Both runs below pass **341 always-active gates** with
byte-identical LF [output](continuing3_20260906_laurent_finite_phase.out).

```text
python -B 04-computation/continuing3_20260906_laurent_finite_phase.py --certificate 04-computation/continuing3_20260906_laurent_finite_phase_certificate.json
python -B -O 04-computation/continuing3_20260906_laurent_finite_phase.py
```

```text
source      8c49a16864d2b17c3fc888df0313e241b6fd5cdfb47c55a06920d709902269d5
output      85de04cb1e32eae863d12c321e0defe35fddd68ddc19182164f2400c27ad3133
certificate afe514688067ab0fdc38f33ca650fbecdb0ddb40e08466b6395c7c5c29659385
```

The proof, certificate and independent audit are filed with their frozen
source/output identities. No new theorem ID is allocated.
