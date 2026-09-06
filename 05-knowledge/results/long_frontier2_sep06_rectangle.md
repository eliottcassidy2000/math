# An open coefficient region with negative full Laurent response

**Status: PROVED + FINITE-EXACT;
[independent analytic and full coefficient audit passed](long_frontier2_sep06_rectangle_audit.md).**
The proof and complete coefficient certificate below extend the previous
one-dimensional four-anchor fibre to a three-dimensional coefficient prism.
No interlacer hypothesis is used. The general two-anchor model remains OPEN.

## Inheritance, connection and scope

The closest proved mechanism is the all-root fibre certificate in
[long_frontier_sep06_residue_tail.md](long_frontier_sep06_residue_tail.md),
with its [exact admissible domains](long_frontier_sep06_fibre_domain.md).
Its hostile is f=6: the full response becomes positive at an original phase,
but the beta polynomial has lost its nonnegative-root predicate. The corrected
near miss is applying the incoming tail2500 without its two interlacers.
The least-used sidecar is the full bivariate dependence of the eliminated
degree-eight response, together with a first-critical-point bound on f.

The board is: coefficient region, original phases, complete carried response,
rooted beta domain, positive polynomial charts, and canceled-root boundaries.
The source-to-target map is elimination of f at the original zero. It preserves
the response there and forgets which coefficients are actual factorial rows.
Phase intervals restore root coverage; a beta-root barrier restores the entire
admissible f range. Cheap hostile probes retain f6 and a positive off-root
response at the genuine centre. No map from every model coefficient to an
actual Laurent row is asserted.

The incoming [degree-eight moment decoder](continuing4_20260906_moments_packet.md)
is exact for the two-interlacer domain. It is not needed for the broader
B-only theorem here. Both results can be used when subdividing the remaining
coefficient region, but their hypotheses must be kept separate.

## 1. The coefficient-prism theorem

Write a=e3, b=e4, f=e5, and define

    B(v)=v^5-13v^4+55v^3-av^2+bv-f,
    C(v)=v^4-12v^3+45v^2-(2a/3)v+3b/7,
    D(v)=v^4-11v^3+36v^2-(5a/12)v+b/7.

Keep the original shifted coefficient carriers

    beta(t)=t^-1(1+13t+55t^2+at^3+bt^4+ft^5),
    c(t)=t^-1(1+12t+45t^2+(2a/3)t^3+(3b/7)t^4),
    d(t)=t^-1(1+11t+36t^2+(5a/12)t^3+(b/7)t^4),
    O(t)=sum_j binom(14,2j+1)t^j,
    E(t)=sum_j binom(14,2j)t^j,
    P=O star beta,
    Q=(O^2+t^-1 E^2) star (beta^2+2t c d),

where star is coefficientwise multiplication. The negative-support carry
of Q is q_(-1)=28. No coefficient or Laurent shift is dropped.

**Theorem 1.** Suppose

    167/2 <= a <= 169/2,
    69/2 <= b <= 71/2,
    0 <= f <= 5.                                           (1)

Every original zero of P is strictly negative and simple. If f>0 there are
exactly four, written -s with one s in each interval

    I1=(99/10000,1/100), I2=(1/9,13/100),
    I3=(1,8/5),           I4=(10,infinity).                  (2)

If f=0, there are exactly three, in I1,I2,I3. At every one,

    Q(-s)<0.                                               (3)

This holds throughout the entire coefficient prism, without requiring B
to be real-rooted and without requiring C or D to interlace it.

**Corollary 2.** For every nonnegative-root B with a,b in the rectangle(1),
the conclusion holds at every original zero, with no further restriction
on f and no interlacer assumption. In fact its admissibility forces f<5.

The rectangle is not merely a formal box disjoint from genuine B geometry:
the whole smaller prism with the same a,b ranges and 1/2<=f<=3/2 has five
distinct positive B roots.

## 2. All original roots are covered uniformly

Direct extraction gives P(-s)=2002 p(s), where

    p(s)=f s^4-(12b/7)s^3+a s^2-10s+1/11.                  (4)

At a fixed s this is affine in a,b,f. Thus its extrema on(1) occur at the
eight corners. The certificate checks every corner at the seven endpoints
in(2), with the following strict uniform signs:

| Interval | Left sign of p | Right sign of p |
|---|---:|---:|
| I1 | + | - |
| I2 | - | + |
| I3 | + | - |
| I4 | - | + at infinity when f>0 |

There are56 finite endpoint checks. For f>0, the positive leading
coefficient supplies the last sign; four disjoint sign changes exhaust
the degree-four polynomial. Hence each root is simple and there are no
additional real or complex zeros. When f=0, the degree is three because
b>0, and the first three intervals exhaust it. No limiting assertion
about a finite fourth zero at f=0 is made.

## 3. The full response after elimination

At an original zero, equation(4) gives exactly

    f=12b/(7s)-a/s^2+10/s^3-1/(11s^4).                     (5)

Substitution in the complete Q yields

    s Q(-s)=-(14/11) H(a,b,s),                             (6)

where H=sum_(j=0)^8 h_j s^j has the following complete coefficients:

    h0 = -443993,
    h1 = 73031400,
    h2 = -4851990a-2871286275,
    h3 = (1845234915a+54465840b)/7+23817753960,
    h4 = -(19072375905a+2760966000b)/7+7912347300,
    h5 = 296928390a^2/7-1763183565a/2+3365140350b,
    h6 = 26558675a^2/2-3563140680ab/49+7314259095b/7,
    h7 = (-845791650ab+986396400b^2)/49,
    h8 = 143416845b^2/49.                                  (7)

This is a19-monomial identity in a,b,s. In particular, at a=84,b=35
it recovers the earlier one-variable polynomial coefficient for coefficient.
The source reconstructs Q directly from literal convolutions before making
this substitution; equation(7) is not inserted as the computation's input.

Let a- and a+ be the rectangle endpoints, and likewise b-,b+. For each
finite phase interval [l,r] in(2), expand the polynomial

    (1+u)^2(1+v)^2(1+w)^8
      H((a-+a+u)/(1+u), (b-+b+v)/(1+v), (l+rw)/(1+w)).    (8)

For the tail use

    (1+u)^2(1+v)^2 H((a-+a+u)/(1+u),
                       (b-+b+v)/(1+v),10+w).             (9)

Every coefficient of u^i v^j w^k, for i,j=0,1,2 and k=0,...,8,
is strictly positive in each of the four charts. Thus there are324
positive exact coefficients. The complete rational lists, rather than
samples of polynomial values, are pinned in the
[certificate](long_frontier2_sep06_rectangle_certificate.json).
Their minima are:

| Chart | Minimum coefficient |
|---|---|
| I1 | 13592097147660088381611435635279441/80000000000000000000000000000000 |
| I2 | 176722546078297280862289/392000000000000000 |
| I3 | 16827346271419/392 |
| I4 | 682807599045/196 |

For finite intervals, this is a strictly positive Bernstein expansion
after clearing positive denominators. Positivity holds on the closed
rectangle, including its boundary faces: each boundary specializes to
a positive face of the homogeneous coefficient array. The tail polynomial
has a positive constant coefficient and all its remaining coefficients
positive for w>=0. Therefore H>0 on all four phase regions. Equation(6)
and s>0 prove Theorem1.

## 4. Recovering the full nonnegative-root beta domain

Put F_(a,b)(v)=B(v)+f. On 0<=v<=2/5 it is bounded above by

    F_upper(v)=v^5-13v^4+55v^3-(167/2)v^2+(71/2)v.

Also F'_(a,b)(2/5)<=F'_upper(2/5)=-81/10<0, while
F'_(a,b)(0)=b>0. The transform of5-F_upper on[0,2/5] is

    (1+w)^5 [5-F_upper((2w/5)/(1+w))]
      =5+(54/5)w+(164/25)w^2+(34/25)w^3
        +(983/625)w^4+(3008/3125)w^5.                    (10)

Every coefficient is positive, so F_(a,b)(v)<5 throughout that interval.

If B has five nonnegative roots, write them beta1<=...<=beta5.
Its derivative has four real roots interlacing them, with multiplicities.
The least derivative root gamma is positive and lies below2/5 by the
two derivative signs. It lies between beta1 and beta2, where B(gamma)>=0.
The same weak inequality holds if these roots coincide. Hence

    f <= F_(a,b)(gamma) <5.

Nonnegative beta roots give f>=0. This proves the full-domain corollary
without borrowing a tail that needs C/D interlacing.

For open nonvacuity, take the whole a,b rectangle and1/2<=f<=3/2.
At v=0,1/10,1,3,5,7 the B signs alternate -, +, -, +, -, +
uniformly. The eight corners give respectively the value ranges

    [-3/2,-1/2], [115871/100000,226871/100000],
    [-17/2,-11/2], [33/2,59/2], [-133/2,-71/2],
    [1117/2,1231/2].                                     (11)

Affineness again extends these signs to the whole smaller prism. Five
disjoint sign changes exhaust B's degree and give five distinct positive
roots. This also shows that neither repeated roots nor a zero product
are secretly required by the result.

## 5. Failure boundaries, verification and next question

At a=84,b=35,f=6 there is an original phase in
(16693/2000,41733/5000). The endpoint p signs are opposite and every
coefficient of the transformed H on this interval is negative, so the
full response is strictly positive there. This is an exact hostile to
discarding the f bound; it is outside the nonnegative-root B domain.
At the genuine centre a=84,b=35,f=1, the value s=4 is not an original
zero and 4Q(-4)=350398552675052>0. Thus the original-zero condition is
essential even inside the valid coefficient region.

The [source](../../04-computation/long_frontier2_sep06_rectangle.py) has
457 always-active gates, including complete coefficient arrays, every
affine corner, the B-root barrier and the two hostiles. Normal and
optimized output and regenerated certificates agree byte for byte:

    python3 04-computation/long_frontier2_sep06_rectangle.py --certificate /tmp/rectangle.json
    python3 -O 04-computation/long_frontier2_sep06_rectangle.py --certificate /tmp/rectangle-opt.json

Raw LF SHA256:

    source 4fb5fe880fbc030b3461889d021f7fa24b6376b23cdc1251afd99b0327591ea5
    output 572d5f8192b0a567efd95b6db25d23f8dcd7277477dacfd27d75a98134668bad
    certificate cd1af1f580b51fc1b556a7cdf44f28203e464a3ba41e3d03e55363278816d412

The universal statements follow from the polynomial identities and positive
coefficient lists, not from numerical sampling. This is a model region;
only the inherited centre f1 is identified here with an actual factorial
row. There is no claim that the prism is maximal. The next useful boundary
question is which failure first obstructs an enlarged region: a phase
chart, the product bound, or the actual response sign. The exact degree-eight
model decoder can distinguish lost geometry from a genuine sign change.
