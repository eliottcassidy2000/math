# Endpoint 39: seven genuine first channels and complete doubled separation

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The exact parameter certificate has 258 strictly positive
coefficients. This note extends the established endpoint-33 method to
endpoint 39, with an unbounded opposite endpoint. It asserts no new
all-channel method and no general trinomial two-rung theorem.

## 1. Statement and exact coefficient scope

Let g>=20 be an integer with gcd(g,39)=1, and let alpha,beta,gamma be
arbitrary nonzero complex numbers. For

    f(u)=alpha*u^(-39)+beta*u^(2g-39)+gamma*u^(3g-39),

the first nonzero constant-term moment occurs at exactly g or 2g.
Both alternatives occur for every eligible g. In particular,

    first nonzero mass <=2g<3g=39+(3g-39).

There are seven channels in the first support-return fibre and fourteen
in the doubled fibre. All six distinct canceling phases of the first
moment have a nonzero doubled moment. The sign established below is a
sign of the displayed normalized coefficient response; a complex raw
moment is not assigned an order sign.

The condition gcd(g,39)=1 belongs to the first-mass statement, not to
the polynomial sign certificate. At g=21 the same coefficient formulas
apply to support (-39,3,24), but its first support return is seven.
Thus deleting the gcd condition changes the claimed first mass.

## 2. Inheritance, comparison and the retained concepts

The closest proved mechanism is the
[endpoint-33 actual quotient theorem](second_20260906_laurent.md), with its
[independent full polynomial reconstruction](second_20260906_laurent_audit.md).
The polynomiality and degree argument below specializes the existing
[all-h characteristic interface, Section 5](overnight7_20260906_laurent_midpoint_transport.md).
The negative and simple first roots come from **THM-4436**,
[complete factorial-row simple negative roots](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md),
with its exact parameter substitution stated below.

The canonical hostile is the actual endpoint-15 derivative-cone separator
recorded in the endpoint-33 note: the whole permitted derivative cone
can miss the genuine carried response even when that response is negative.
The corrected near miss is the endpoint-33 square-phase multiplier whose
coefficients lose positivity at x=10000. The newer
[sparse positive-amplitude repair](continuing2_20260906_laurent_sparse_amplitude.md)
works at one endpoint-27 coefficient point, with algebraic coefficients on
selected positive square-root branches. It does not provide a uniform
endpoint-39 representation or a rational identity on all conjugate phases.
The least-used retained sidecar here is the exact lower-carry coefficient
before taking the first-root quotient.

The concept board is: complete integer count fibres; original coefficient
phase; lower carry; multiplication in the first-root quotient; shifted
characteristic positivity; and independent literal-fibre reconstruction.
The incoming amplitude result changes the possible representation of a
fixed response but does not replace its coefficient law. The cheap hostile
probes therefore used the actual h=6 response at x=1,2,3,10,100. All six
characteristic coefficients were positive at each point, motivating the
full symbolic computation rather than an extrapolated parameter claim.

Targeted searches at inherited HEAD f0521b87281f for endpoint39,
seven-first-channel closure and the corresponding supports found only
the explicitly open next step in the endpoint-33 report. The general
degree interface and endpoint-33 mechanism are credited as inheritance;
the new proof object is the endpoint-39 positive certificate. No external
priority claim is made.

## 3. The complete actual first and doubled rows

Put

    x=g-19>=1, tau=alpha*gamma^2/beta^3,
    X=alpha^x beta^18 gamma, K=(2g)_26>0,

where (a)_m is the falling factorial. At mass g the charge equation is
2*n_beta+3*n_gamma=39. It forces n_gamma odd, and its complete
nonnegative solutions, including the alpha count, are

    (n_alpha,n_beta,n_gamma)=(x+j,18-3j,1+2j), j=0,...,6.

At mass 2g, the complete solutions are

    (2x+e,36-3e,2+2e), e=-1,...,12.                       (1)

The lower carry e=-1 has counts (2x-1,39,0) and is valid for every x>=1.
No selected semigroup subfamily or midpoint model replaces these fibres.

Define

    p_x(t)=sum_(j=0)^6
      [13! (x+6)_(6-j)/((18-3j)!(1+2j)!)] t^j,

    q_x(t)=sum_(e=-1)^12
      [(2x+12)_(12-e)/((36-3e)!(2+2e)!)] t^e.             (2)

The first polynomial is monic. Direct multinomial expansion gives the
exact identities, in the same original coefficient phase,

    CT(f^g)=X binom(g,13) p_x(tau),
    CT(f^(2g))=X^2 K q_x(tau).                           (3)

THM-4436 applies with (A,B,h,r,z,x)=(2,3,6,0,1,x): its first row
is exactly a positive multiple of p_x. Therefore p_x has six distinct
strictly negative real roots at every integral x>=1. Also

    p_0(x)=13! (x+6)_6/18! >0,                           (4)

so none of these roots is zero and the inverse phase in (2) is legitimate.

For any support-return mass m, the charge equation modulo g gives
g | 39m. Under gcd(g,39)=1, it follows that g | m. The first fibre
above is nonempty, so g is exactly the first support-return mass. The
remaining obligation is the actual same-root inequality

    q_x(rho)<0 whenever p_x(rho)=0.                      (5)

Off those six phases the first moment is nonzero. Each root is attainable
by taking beta=gamma=1 and alpha=rho. The all-positive coefficient choice
gives the alternative first mass g.

## 4. Polynomial carry cancellation and the degree-bound certificate

Initially work over Q(x)[t]/(p_x), where t is invertible. The lower carry
has the exact polynomial cancellation

    [q_x]_-1 / p_0(x)
      = [128*18!/(39!*13!)] x product_(j=0)^5(2x+2j+1).   (6)

Indeed (2x+12)_13 has seven even factors 2x,2x+2,...,2x+12
and six odd factors 2x+1,...,2x+11. The six factors x+1,...,x+6
cancel (4), leaving exactly (6). This proves polynomiality after
multiplication by the carry coefficient; it does not assert that t is
invertible over the unspecialized ring Q[x][t]/(p_x).

Use ordinary monic division for the nonnegative powers, and use

    t^-1=-sum_(j=1)^6 p_j(x)t^(j-1)/p_0(x)

for the inverse term. The full remainder is

    R_x(t)=sum_(j=0)^5 R_j(x)t^j in Q[x,t],
    q_x(t)=R_x(t) mod p_x(t), deg_x R_j<=12-j.            (7)

For the degree bound, give both x and t weight one. Each relation
coefficient p_j has degree 6-j. Every nonnegative term of q_x has
total weight 12. Monic reduction preserves that bound. For the inverse
term, (6) has degree seven and p_(j+1) has degree 5-j, giving the same
bound in (7).

Let M_x be multiplication by R_x on the quotient basis 1,t,...,t^5.
Its entry at row i, column j has degree at most 12+j-i. Therefore

    det(zI-M_x)=z^6+c_1(x)z^5+...+c_6(x),
    deg_x c_k<=12k.                                     (8)

To see this without sampled-degree inference, express c_k as the signed
sum of k-by-k principal minors. Within any determinant term, the row
and column index sums cancel. The bound is thus 12k for every term.

## 5. The positive certificate proves the actual sign

The [standalone source](../../04-computation/long_frontier_sep06_endpoint39.py)
constructs (2), (6), (7) and (8) over exact rational polynomial rings.
For each k=1,...,6 it obtains

    c_k(x)=d_k sum_(j=0)^(12k) a_(k,j)(x-1)^j,
    d_k>0, and every a_(k,j) is a strictly positive integer. (9)

The full six d_k and coefficient arrays are frozen in the six
`CHAR_CERTIFICATE` lines of the
[matching output](long_frontier_sep06_endpoint39.out). There are exactly
13+25+37+49+61+73=258 strictly positive entries. These complete data,
not a short parameter table, are the finite certificate in the proof.

For every real x>=1, (9) gives c_k(x)>0. Hence the characteristic
polynomial in (8) is strictly positive at every real z>=0.
For a real root rho of p_x, evaluating the quotient multiplication
identity gives

    det(q_x(rho)I-M_x)=0.

Equivalently, apply the verified Cayley-Hamilton identity to the class
of 1 and then evaluate at rho. It follows that the real number
q_x(rho) cannot be nonnegative, proving (5). For integral x>=1,
THM-4436 places all six first roots on that real axis. Equation (3)
therefore gives nonzero doubled moment at every first cancellation.
Together with the support-return gcd argument this proves Section 1.

The certificate proves the sign at every real root for all real x>=1.
The claim that there are six simple negative first roots is used only
for integral x and follows from the stated inherited theorem. No
unproved extension of its parameter domain is required.

## 6. Independent controls, scope and reproduction

The symbolic characteristic coefficients are computed from matrix powers
and Newton trace identities. The program separately checks the complete
symbolic Cayley-Hamilton matrix identity. A second path enumerates every
nonnegative count triple at masses g and 2g for

    x in {1,2,...,12,100}, g=x+19.

It reconstructs both full rational coefficient rows from their literal
multinomial weights, verifies the inverse phase against the original
first equation, and verifies that dropping the carry changes the quotient
response. It then computes all six characteristic coefficients by sums
of principal minors, using rational elimination, and compares them with
the symbolic polynomials. These 13 samples include eight primitive
first-return cases; nonprimitive samples test only the coefficient laws.
The separately enumerated g=21 support has first return seven and records
the exact failure of the omitted-gcd claim.

These controls use a different determinant path and literal fibres but
are not an independent agent audit. The certificate is the universal
proof input; the finite samples are implementation checks. No numerical
root solver, free coefficient model, deleted carry, retuned first phase,
or extrapolation in h enters the argument.

Reproduce from the repository root:

```bash
python3 -B 04-computation/long_frontier_sep06_endpoint39.py
python3 -B -O 04-computation/long_frontier_sep06_endpoint39.py
```

The frozen source passes **421 always-active exact gates**. Normal and
optimized outputs match the saved output byte-for-byte. Raw LF hashes:

```text
source SHA256 a92bab693c728ab98257f20842d53eafd33faef0cd8023c35633b143d07e1c9f
output SHA256 4c686cea98965eb7b4ea785f87a73592dc959ec6d19f47ebc6ad7e7dd0c5a5bb
semantic SHA256 021c9851d4ffa39bc698ae718a02faca35025bc3886773e0b6e84d22935290fb
```

The [independent analytic and full certificate audit](long_frontier_sep06_endpoint39_audit.md)
passes. It reconstructs all six polynomial identities from73 distinct
literal-fibre parameters, using the proved degree bound72 and a separate
Berkowitz characteristic algorithm. Its744 always-active gates pass
normally and with optimization; no mathematical correction was needed.

The source-to-target map takes the complete actual Laurent coefficient
response to multiplication in the quotient by its actual first equation.
It preserves every first-root value and the same-root noncancellation
predicate, while discarding values off that locus. The retained literal
rows, coefficient monomial X, lower carry and original tau supply the
interpretation of that quotient. This is an actual coefficient theorem.
The all-h shifted-characteristic positivity question remains open; this
endpoint certificate supplies no uniform proof as h varies.
