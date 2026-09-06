# Independent referee of the four-anchor full-response sign theorem

**Status: INDEPENDENT ANALYTIC AND EXACT POLYNOMIAL AUDIT — PASS.**
This audits [the residue-tail theorem](long_frontier_sep06_residue_tail.md),
its [producer](../../04-computation/long_frontier_sep06_residue_tail.py),
and [frozen output](long_frontier_sep06_residue_tail.out). No correction is
needed. The separate [fibre-domain classification](long_frontier_sep06_fibre_domain.md)
is a different proof obligation; its stronger domain assertions are not
inferred from the moment determinant audited here.

## 1. Exact accepted domain and analytic argument

The sign conclusion holds for all real f in [0,5], with the fixed induced
arrays GB=(1,13,55,84,35,f), GC=(1,12,45,56,15), GD=(1,11,36,35,5).
No interlacer hypothesis is needed for the sign conclusion on this box.
If the corresponding B has five nonnegative roots and D weakly interlaces
B, the independently checked necessary determinant bound places f in
[0,14sqrt(2)-16], inside that sign box. The implication is sufficient for
entry into the box, not an iff characterization of the model.

Weak interlacing and cancellation give a positive probability residue
measure for D/B. Repeated roots create no additional pole order after
the common factors are cancelled; cancelled nodes may have zero weight.
The nonnegative-root assumption permits the shifted Cauchy inequality.
Direct series division recovers

    nu0=1, nu1=2, nu2=7,
    nu3=7A/12-19,
    nu4=115A/12-632-(6/7)b.

The shifted two-moment minor gives A>=522/7. The three-moment determinant
gives 2592b<=-343A^2+67788A-3157056. At A=84,b=35, a separate exact
division through nu8 and a direct symbolic determinant give

    det(nu_(i+j))_(i,j=0)^4
      =(f^2-15f-225)(f^2+32f-136).

AM-GM on the five nonnegative summands of e4 gives f^4<=7^5<12^4.
For 0<=f<12 the first factor is negative, so positivity of the Gram
determinant forces f<=14sqrt(2)-16<19/5. This audits the stated necessary
D-domain bound, including zero and repeated-root boundaries.

## 2. Independent reconstruction from the ordinary arrays

The reconstruction did not import the producer or its polynomial package.
It used ordinary dictionaries indexed by Laurent exponent, multiplied
by direct convolution, with

    O_j=binom(14,2j+1), E_j=binom(14,2j),
    beta_j=(GB)_(j+1), Craw_j=(GC)_(j+1), Draw_j=(GD)_(j+1).

Thus beta,Craw,Draw retain their exponent -1. Form O^2+z^-1 E^2 and
beta^2+2z Craw Draw by convolution, and multiply matching Laurent
coefficients. The resulting Q retains exponent -1 with coefficient 28.
The first row is obtained independently by matching O and beta:

    P(-s)=2002[f s^4-60s^3+84s^2-10s+1/11].

Only after both full rows were formed was

    f=60/s-84/s^2+10/s^3-1/(11s^4)

substituted. Exact cancellation gave sQ(-s)=-(14/11)h(s), with all nine
coefficients, in descending order, equal to

    (3585421125, -26087589000, -83518139925,
     343030019640, -234760993560, 46232902140,
     -3278853435, 73031400, -443993).

This verifies the complete lower-carry response at the original equation.
It does not replace the original root by a contiguous-row root or discard
a negative exponent. Since s>0, multiplication by s preserves the response
sign. The f=1 interpretation as the first and doubled rows of (-27,1,15)
is also checked by the producer's fifteen literal factorial channels;
other f values remain coefficient-model deformations.

## 3. Exhaustion, positive transforms, and boundaries

An independent exact check evaluated p_f at both f=0 and f=5 at each of
the seven endpoints

    99/10000, 1/100, 1/9, 13/100, 1, 3/2, 10.

Both sign lists are (+,-,-,+,+,-,-). These fourteen endpoint signs apply
to the full coefficient interval because p_f is affine in f with
nonnegative coefficient s^4. For f>0 there is a root in each of the three
finite sign intervals and a fourth after 10; degree four forces four
distinct simple positive roots and excludes any additional complex roots.
For f=0 the degree drops to three and the three finite intervals exhaust
it, again with simple positive roots. No complex or repeated original
roots are silently omitted.

Every coefficient was independently checked positive in each of

    (1+t)^8 h((a+b*t)/(1+t)),
    (a,b)=(99/10000,1/100), (1/9,13/100), (1,3/2),

and in h(10+t). Each finite transform covers its closed phase interval:
t>=0 covers a<=s<b, and the positive leading coefficient equals h(b).
The tail transform covers s>=10. Therefore h is positive at every
original root in the sign box, proving the claimed strict Q negativity.
As f tends down to zero, the three finite roots remain in their fixed
intervals; the fourth tends to infinity because the total root sum is
60/f. The same tail transform covers this entire escape.

The outside-box f=6 hostile was checked independently as well. Its
polynomial has signs (-,+) at
(16693/2000,41733/5000). All nine coefficients of the finite transform
of h over this interval are negative. Thus this bracket contains an
original root with h<0 and Q>0. The producer's other disjoint brackets
exhaust the remaining roots and identify this as the largest one. Its
D moment determinant is -25668, so it does not contradict the stated
interlacer domain. The coefficient bound is load-bearing.

## 4. Producer inspection and independent replay

The entire standalone source was read. Its sparse rational algebra retains
all four formal variables until their declared specializations, computes
the fifth moment determinant by the permutation formula, and forms full
raw first and doubled rows before elimination. Its second path uses
ordinary binomial carriers and coefficient extraction rather than the
Hadamard construction. The signs of all fifteen named control roots are
certified by rational brackets and interval Horner evaluation. Checks are
explicit runtime gates, so optimized Python does not disable them.

The producer was independently run in normal and optimized modes and
compared byte for byte with its frozen output:

```
python3 -B 04-computation/long_frontier_sep06_residue_tail.py > /tmp/residue-tail-certificate-audit.normal.out
python3 -B -O 04-computation/long_frontier_sep06_residue_tail.py > /tmp/residue-tail-certificate-audit.optimized.out
cmp /tmp/residue-tail-certificate-audit.normal.out /tmp/residue-tail-certificate-audit.optimized.out
cmp /tmp/residue-tail-certificate-audit.normal.out 05-knowledge/results/long_frontier_sep06_residue_tail.out
```

All three outputs agree, with 174 explicit gates. Freshly computed hashes:

* Source: `a89ede31236a8caeaace4e89e7aa8c47efffc97a417944bf3badf6d54cc46acd`.
* Output: `859a7b8a9d6610c45656ebff8b436bfe96a520df43a56325c9db048ed46a08ac`.
* Semantic: `e615355a57495d52795eeeac6fb7dc877eec7adae1166525a987be16f3d39828`.

The independent reconstruction used local symbolic polynomial algebra only
as an exact identity checker; the primary certificate remains standard-library
rational arithmetic. Its formula inputs and all tested identities are
specified in Sections 1--3. No numerical root approximation, shape census,
or separate domain-classification claim is a dependency of this acceptance.
