# Independent audit of the smallest-phase and 75000-tail model theorem

**Verdict: PASS, analytically and by an independent exact certificate path.**
No mathematical repair is requested. The accepted statements are precisely
Theorems A and B of `continuing3_20260906_laurent_finite_phase.md`, together
with its exact relaxation hostiles. Full model negativity and general actual
Laurent noncancellation remain OPEN. This audit does not promote a finite
control bank into an unproved whole-phase statement.

The producer's pre-promotion report SHA-256 is
`b9f335982e8b62e3beac190f365634ae51580f4ed851161902e7b102f722f8ba`.
Its source, output, and JSON are byte-pinned by the independent verifier.
The full report and source were reviewed after independently deriving the
eight-corner sign mechanism below. No producer is imported or executed by
the referee source.

## 1. Model and original normalization

The source object is a nonnegative five-root shape with sum 13 and square
sum 59. These are the exact anchors e1=13 and e2=55. In the full class K,
the two specified monic contiguous quartics C and D weakly interlace the
degree-five polynomial B; zero and repeated beta roots are allowed.

The referee reconstructs the entire original coefficient product from

```text
alpha = t^-1 (t O^2+E^2),
beta packet = t^-2 (G_B^2+2t G_C G_D).
```

It multiplies the literal binomial polynomials O and E before comparing
every alpha coefficient with binom(28,2j+2). The resulting full Laurent row
has exactly q_-1=28; hence T(s)=sQ(-s) has constant term -28. The original
crossing factor 2t and both lower shifts are retained. Independent symbolic
extraction gives exactly

```text
P(-s)=182-20020s+2002x s^2-3432y s^3+2002z s^4.
```

The actual inherited endpoint-27 coefficient control (x,y,z)=(84,35,1)
gives 2002 times its monic first polynomial, with all signs unchanged.
Every coefficient of the JSON's T polynomial agrees with this literal
construction. The model scope is not enlarged to arbitrary real-rooted
carriers, arbitrary multiplier normalizations, or all Laurent supports.

## 2. First-branch proof and a different whole-box sign certificate

The elementary root bounds in producer Section 2 are valid, including the
two five-variable polynomial identities. The maximum-root Cauchy bound
gives x<123; the fourth-coefficient sum of squares gives y<152; AM-GM on
the ten pair products gives z<72. The anchors force at least three positive
roots, so x>0. The identity 4y<=13x proves P(-s)>0 through s=1/110.
The upper endpoint and derivative bounds prove a unique simple root in
(1/110,1/90), and no earlier positive root. Neither interlacer is needed
for these statements.

Here is an independent replacement for the tensor sign proof. On the full
phase interval [1/120,1/80], half the second derivatives of T with respect
to x,y,z are respectively

```text
198835 s^5(66-85s),
(2005830/7)s^7(140-13s),
13123110 s^9.
```

All are positive. At fixed s, T is therefore separately convex in the three
coefficient coordinates. Applying the elementary chord bound successively
in x, y, and z bounds T by the maximum of its eight corner values on
[0,130] x [0,152] x [0,72]. No joint convexity is asserted or needed.

The referee reconstructs each corner polynomial T+400 over the rationals.
Its exact Sturm variations at 1/120 and 1/80 are, in lexicographic corner
order,

```text
(2,2), (4,4), (4,4), (5,5), (4,4), (5,5), (6,6), (7,7).
```

Thus all eight have zero roots in the interval; every one is negative at
the endpoints and at 1/100. This proves **T(s)<-400 on the entire continuous
four-dimensional box**, with the full advertised margin, by a certificate
independent of the producer's tensor Bernstein signs.

Separately, all 270 frozen Bernstein entries are checked to be below -400.
Their entire power-basis polynomial is reconstructed by the inverse
finite-difference identity

```text
c_alpha = product_j binom(d_j,alpha_j)
          * sum_(beta<=alpha) (-1)^(|alpha-beta|)
            product_j binom(alpha_j,beta_j) b_beta.
```

Every coefficient agrees with the independently reconstructed T after the
stated affine box substitution. The claimed exact maximum coefficient is
also verified. Consequently both the original certificate and the different
eight-corner proof pass.

## 3. Original-root elimination and the full uniform tail

The inherited supplier is
`05-knowledge/results/continuing2_20260906_effective_anchor.md`, with its
independent audit: **e4>1/100 on K**. That statement requires both specified
interlacers. It is not inferred here from the coefficient box or from the
first-row polynomial's roots.

Substitution of the exact original equation yields
z=(12/7)yu-xu^2+10u^3-u^4/11, u=1/s. Independent symbolic division and
substitution recover every coefficient of Q(-s)/s^7 in the report and JSON.
After division by y^2, the retained constant is -26075790/7. Every positive
monomial, and only those positive monomials, contributes to the displayed
envelope using x<123 and y>1/100. All other discarded monomials are
nonpositive for x,y,u>=0. The coefficient bounds are therefore in the
correct direction, with no missing high or low carry.

The exact envelope at u=1/75000 is

```text
-470440303496224131940323407036971853244343
/3854347229003906250000000000000000000 < -120000.
```

Its nonconstant coefficients are strictly positive. Monotonicity in u
therefore proves Q(-s)/(s^7 e4^2)<-120000 for every original positive root
s>=75000, including multiple roots. Both the original-root condition and
the interlacer-supplied coefficient floor remain attached to this theorem.

## 4. Hostiles, missing predicate, and remaining scope

The exact shape (1,3,9,22,30)/5 has both anchors and both strict interlacers.
It gives positive Q(-75000) away from the first-root equation. This correctly
refutes any attempted deletion of that equation from the tail assertion.
The coefficient-box point (130,0,0,1/60) has the published positive T value;
it is explicitly outside the beta-root shape class.

For the stronger relaxation hostile
(x,y,z)=(104,50,37435088/3898125), the referee verifies the complete first
factorization, its four distinct positive roots, all three other rational
brackets, and the exact response

```text
Q(-15/2)=78541969368658673/18480>0.
```

All four displayed Newton inequalities hold strictly. Independent Sturm
counts show that B nevertheless has exactly one real root, which is
positive. It has four nonreal roots and is not in K. Thus first-row real
rootedness and these coefficient inequalities do not restore the lost
beta-root predicate. This is a correctly typed obstruction to a stronger
surrogate, not a model counterexample.

The ancillary indefinite-Hessian claim also checks exactly: its determinant
is

```text
-(1031232600/49)u^2(144097056u^2+72692884u+10966105)<0,
```

for u>0. Unconstrained joint concavity therefore provides no shortcut.
Theorems A and B leave only phases in (1/80,75000) unsettled, at most three
per shape counted with multiplicity. No assertion that such phases exist,
or that a sampled search resolves them, is required.

## 5. Frozen reproduction

The independent source and outputs share this report's stem. The verifier
uses Python 3 and SymPy, imports no repository implementation, and performs
**663 always-active exact gates**. Normal and optimized runs are byte-identical
LF. It reads the JSON and producer source beside itself, and reads the
producer output beside itself or from `../05-knowledge/results` when filed.

```text
python3 -B 04-computation/continuing3_20260906_laurent_finite_phase_audit.py
python3 -B -O 04-computation/continuing3_20260906_laurent_finite_phase_audit.py
```

Referee raw-LF SHA-256 hashes:

- source: `f56a5cfa459e24be65b2865379c8e18fd908e3e2b991a4afd002c11b34528167`;
- each output: `30518ee02c155c3e0bd2dc573ffa8330c5ac7ab34dbb49549f67c085e6c16eee`.

Frozen producer pins checked by the source:

- source: `8c49a16864d2b17c3fc888df0313e241b6fd5cdfb47c55a06920d709902269d5`;
- output: `85de04cb1e32eae863d12c321e0defe35fddd68ddc19182164f2400c27ad3133`;
- JSON: `afe514688067ab0fdc38f33ca650fbecdb0ddb40e08466b6395c7c5c29659385`.

No producer file, shared mathematical document, theorem namespace, or Git
state was changed during this independent audit.
