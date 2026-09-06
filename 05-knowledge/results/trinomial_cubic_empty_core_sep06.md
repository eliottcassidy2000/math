# A three-minor certificate for an unbounded cubic first-return family

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
For every integer `g>=11` with `gcd(g,21)=1`, a Laurent polynomial with
nonzero complex coefficients on

```text
(-21, 2g-21, 3g-21)
```

has first nonzero constant-term moment exactly `g` or `2g`. Both orders
are attained on each support. At every one of its three first-return
cancellation parameters, the doubled moment divided by the square of the
first channel's coefficient monomial is strictly negative.

The certificate consists of all three signed leading principal minors
of a weighted trace matrix. Each factors into positive terms after
`g=s+11`. The complete second row has eight channels, with its lower
carry retained. This is one fixed family, not a result for every support
of smaller endpoint twenty-one or a general higher-channel theorem.

## 1. Inheritance and actual coverage of the incoming duplication theorem

The closest proved mechanism is the
[quadratic first-return family](trinomial_width15_empty_core_returns_sep06.md),
where a negative trace and positive norm determine the two signs.
The [Hermite trace formulation](trinomial_trace_sign_empty_core_certificate_sep06.md)
identifies the correct higher-degree replacement: every signed leading
principal minor. Its old control `(-21,1,12)` already has the three
positive minors recovered below; that specialization is inherited.
The corrected near miss is using only trace and determinant when there
are three signs. The least-used sidecar is the middle principal minor.

The incoming
[THM-4440, signed duplication SOS and real-rooted Laurent return](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md),
promoted in `f96ddbf594`, was read with its
[full proof](nc2_signed_duplication_overnight_hexagon_sep05.md).
It closes the arithmetic-progression trinomial case and any applicable
real-rooted compressed core. Here the ordinary compressed core is
`alpha+beta*s²+gamma*s³`. Reversing it and scaling gives

```text
1+s+tau*s³,               tau=alpha*gamma²/beta³.
```

For real `tau<0`, its discriminant is `-tau*(4+27tau)`; it is
real-rooted exactly on `[-4/27,0)`. The scalar `tau` is unchanged under
nonzero scalar and variable gauges. A gauge with real-rooted real core
would normalize back to this same polynomial, so changing gauge does
not enlarge that sector.

This sector covers at most **one** of our three cancellation parameters
for every admissible support; the proof is in Section 5. Thus at least
two roots per support are outside THM-4440's real-core consumer. The
present theorem is an analytic extension on actual uncaptured
coefficient parameters, not a restatement of its AP theorem or an
inference from a nonreal-rooted toy.

Exact family strings, `(21,2,3)`, the named control, and two-rung routes
were searched on the inherited clean `HEAD=b637305db0`. The prior
occurrences recovered the trace control and the explicitly open next
target in the quadratic-family note. The
[incoming finite root/sign bank](overnight3_20260906_moments_root_geometry.md)
already retains the raw lower-carry factor and includes the first two
admissible parameters `g=11,13` here. Those signs and the carry correction
are concurrent inherited evidence. The
[endpoint-eight theorem](overnight_20260906_moments_width8.md) does not
cover these supports: `min(21,3g-21)>=12`. The two-channel theorem
[THM-4432](../../01-canon/theorems/THM-4432-trinomial-two-channel-two-rung-noncancellation-with-carries.md)
does not cover four first channels. No prior unbounded closure of this
family was recovered, and no literature-priority claim is made.

The subsequently read
[Hadamard-transport result](nc2_hadamard_transport_overnight_hexagon_sep05.md)
proves a negative **virtual** doubled row at every first root, including
these cubic rows. It explicitly distinguishes that row from the actual
`tau^-1 Q_g`: coefficientwise domination does not imply domination at
a negative root. Its universal signed-defect transport remains **OPEN**;
its unbounded transport proof is the earlier `(15,2,3)` quadratic
family. Its finite head includes the present `g=11,13` controls, adding
stronger finite defect-sign evidence there. It does not subsume the
unbounded cubic actual-sign or first-detection conclusion proved here.
Conversely, the present direct trace certificate does not prove that
the actual row is below the virtual one.

The live board is: complete affine fibres; carry gauge; ordinary-core
real-root sector; weighted trace inertia; translated symbolic positivity.
The proof couples the first, second, fourth and fifth objects, while the
third identifies the exact overlap with incoming work.

## 2. Full rows and the first-return statement

Let

```text
f(u)=alpha*u^-21+beta*u^(2g-21)+gamma*u^(3g-21),
g>=11 integral, gcd(g,21)=1, alpha*beta*gamma!=0.
```

The support is primitive because its gcd is `gcd(g,21)`. A balanced
multiplicity vector `(x,y,z)` of mass `m` satisfies

```text
g(2y+3z)=21m.
```

Hence `g|m`. Solving `2y+3z=21` and `2y+3z=42` gives the complete
first and doubled rows:

```text
(g-10+j, 9-3j, 1+2j),      j=0,...,3,
(2g-21+j, 21-3j, 2j),     j=0,...,7.
```

All negative-charge counts are positive for `g>=11`. These four and
eight channels are complete. There is no earlier positive support
return and no support return strictly between `g` and `2g`.

Set `u0=g-7` and define

```text
tau=alpha*gamma²/beta³,       X=alpha^(g-10)*beta⁹*gamma,
L(g)=(g)_7/9!,               K(g)=(2g)_14/128>0,

P_g(tau)=72tau³+504u0*tau²+84u0(u0-1)tau+u0(u0-1)(u0-2),
Q_g(tau)=sum_(j=0)^7 (2g)_(21-j) tau^j/[(21-3j)! (2j)!].
```

Here `(v)_k` denotes a falling factorial. The exact moment identities are

```text
CT(f^g)=X*L(g)*P_g(tau),
CT(f^(2g))=X²*tau^-1*Q_g(tau).                       (1)
```

The canonical tuple is `(A,B,h,r,z,x)=(2,3,3,0,1,g-10)`.
The second left anchor is twice the first minus the affine step
`(1,-3,2)`, giving lower carry one and upper carry zero. Its factor
`tau^-1` in `(1)` cannot be suppressed from a sign claim.

## 3. Complete anchored residue and the trace matrix

Modulo `P_g`, the nonzero constant coefficient makes `tau` invertible:

```text
tau^-1= -[72tau²+504u0*tau+84u0(u0-1)]/[u0(u0-1)(u0-2)].
```

Let `R_g` be the degree-at-most-two residue of
`tau^-1 Q_g/K(g)`. Full polynomial division gives

```text
R_g(tau)=r0(g)+r1(g)tau+r2(g)tau²,

r0=-(g-8)(g-7)A4(g)/28510570408320000,
r1=-(g-7)B4(g)/593970216840000,
r2=-C4(g)/4157791517880000,

A4=981739413g⁴-29190405089g³+323857306851g²
     -1589986674295g+2916051235200,
B4=1715049346g⁴-49273046168g³+530071519257g²
     -2530936512265g+4525872382400,
C4=70405697381g⁴-1950797880338g³+20268905917687g²
     -93593992728910g+162061152680400.                (2)
```

Use the quotient algebra over the real numbers, with basis `1,tau,tau²`,
and define

```text
H[i,j]=Tr(R_g*tau^(i+j)),       0<=i,j<=2,
D_k=(-1)^k det H[0:k,0:k],      k=1,2,3.             (3)
```

The trace is the ordinary trace of multiplication in the quotient. It
can be computed with a companion matrix or from the Newton power sums
of the monic cubic. There is no root approximation in these formulas.

The four polynomial factors in the next section give, simultaneously,
the discriminant identity

```text
disc(P_g)=15552(g-8)(g-7)² F3(g)>0.                 (4)
```

Since all coefficients of `P_g` are positive for `g>=11`, `(4)` proves
it has three distinct real roots, all strictly negative. Thus its
Vandermonde evaluation matrix `V[l,j]=lambda_l^j` is invertible and

```text
H=V^T diag(R_g(lambda_1),R_g(lambda_2),R_g(lambda_3)) V.
```

The three inequalities `D_1,D_2,D_3>0` are exactly the leading-minor
criterion for negative definiteness. They therefore force all three
root values negative. This elementary identity is the classical
weighted trace mechanism recovered in the trace companion; no new
abstract positivity theorem is asserted.

## 4. The complete positive polynomial certificate

Define `F_d` by the following descending coefficient arrays of
`F_d(s+11)`, beginning with the coefficient of `s^d` and ending with
the constant coefficient:

```text
F5:
21956511799903, 447248962924222, 3643877376899891,
14842834323259268, 30227960317081116, 24622377792469920

F11:
140892841323838612672637, 5965004072995999656704026,
114714746974061256056387574, 1322756664440858664210176156,
10161118159501638625772346617, 54598896731247155548982886362,
209397347882224306724476576476, 573182745817688552468969130664,
1097399763351738257252121426176, 1399541151141843962715007684992,
1070005363538327691749959345920, 371517273613263339950411923200

F3:
74863, 860377, 3285600, 4167636

F15:
141387288101294140031271809, 7376544414686784800082667655,
179216442364749171385042322678, 2689550519731785864607802436098,
27881104094183316787745471091341, 211468678134040753645990344857219,
1212233732617430261960470627128172, 5347770140720395948826740811577460,
18303547297713040708451344127658688, 48600649295038820551605884973566416,
99288474821888816842823372182321600, 153249236436278889420201068735488192,
172971501480743281449391277285037312, 134765740119228340806711738891048960,
64801772706710330583381523716710400, 14495117169229267635561205094400000
```

These are all **38** coefficients, and every one is positive. The
complete three identities are

```text
D_1=(g-7)F5(g)/d1,
D_2=(g-8)(g-7)² F11(g)/d2,
D_3=(g-8)²(g-7)^4 F3(g)F15(g)/d3,

d1=28510570408320000,
d2=119489335876142491574207692800000000,
d3=41207553558341744647619275991171138077065216000000000000.
                                                               (5)
```

The discriminant factor is equivalently
`F3(g)=74863g³-1610102g²+11532575g-27511000`.
The standalone source checks `(2)`, `(4)` and every identity in `(5)`
over the rational polynomial field. All denominators in `(5)` are
positive. For every real `g>=11`, every remaining factor is positive
by its displayed shift. Consequently `H` is negative definite and

```text
P_g(tau)=0 => tau^-1 Q_g(tau)=K(g)R_g(tau)<0.        (6)
```

There is no equality or common-root case in the domain. Combining `(6)`
with `(1)` and support-mass divisibility proves

```text
m_*(f)=g     when P_g(tau)!=0,
m_*(f)=2g    when P_g(tau)=0.                        (7)
```

Positive coefficients attain the first case. Each of the three scalar
roots is attained by `alpha=tau`, `beta=gamma=1`, so each supplies the
second case. Thus the exact worst detection order per support is `2g`,
strictly below its total exponent width `3g`.

The symbolic sign statement also holds in the real falling-factorial
extension `g>=11`; the first-return interpretation requires integral
`g` and `gcd(g,21)=1`. The excluded parameter `g=12` gives
`(-21,3,15)`, whose first support return is four, not twelve. Its
formal mass-twelve identity is valid; calling it the first row is not.

## 5. Exact real-core boundary and why this adds coverage

On `tau in [-4/27,0]`,

```text
P_g''(tau)=432tau+1008u0>0,
27P_g'(-4/27)=2268(u0-4)²+11844(u0-4)+11216>0.
```

Hence `P_g` is strictly increasing on the entire real-core sector.
At most one of its three roots belongs there. At least two first
cancellation parameters of **every** support are therefore outside
the hypothesis of THM-4440. All three nevertheless satisfy `(6)`.

For precision, substitution at the left boundary gives

```text
P_g(-4/27)=u0³-(139/9)u0²+(2066/81)u0-512/2187.
```

The [independent sector audit](trinomial_cubic_sector_empty_core_sep06.md)
proves this expression negative for
`4<=u0<=13` and positive for `u0>=14`. Thus among admissible integer
parameters, exactly one first root lies in the real-core sector at
`g=11,13,16,17,19,20`; none does for admissible `g>=22`. These sector
counts are a scope diagnostic, not an input to the proof of `(6)`.
The primary exact controls replay all six overlap parameters and two
parameters with no overlap.

## 6. Reproduction, controls and audit

The [standalone source](../../04-computation/trinomial_cubic_empty_core_sep06.py)
and [saved output](trinomial_cubic_empty_core_sep06.out) import no
repository producer or theorem. Symbolic polynomial division retains
the inverse lower carry and reconstructs `(2)` from the full `Q_g`.
Newton sums construct the weighted trace matrix and verify all three
factorizations and all positive coefficients.

The eight named admissible controls are `g=11,13,16,17,19,20,22,29`.
For each, an original charge-equation scan and independent literal
Laurent multiplication agree on both full weighted rows. Every earlier
mass is checked empty. Specialized polynomial inversion and companion
multiplication independently recover the three Hermite minors and
their symbolic values. Exact Sturm counts check all three first roots
and the real-core sector overlap. The nonprimitive hostile above is
also replayed. These are named controls, not a growing support census.

At `g=11`, the inherited signed minors are exactly

```text
752350432547692,
21236578848830718128306804/3,
942083106929885721679784497035400/27.
```

Here the earlier trace note already knew the signs. The contribution
is the uniform parameter certificate `(5)` and its consequence `(7)`.

```sh
python3 -B 04-computation/trinomial_cubic_empty_core_sep06.py
python3 -B -O 04-computation/trinomial_cubic_empty_core_sep06.py
```

All **120 explicit gates** pass; normal and optimized outputs are
byte-identical. The semantic digest is
`c99c8bf9cd37e1e25bcddfaa563b7069942eca2f1bf2b55d213c427d7263b65d`.

Raw SHA-256:

```text
primary source 7ebdffa7cd529cc1d641dc319954c23a594cff5b4e98ae538e21ed18849b77f5
primary output cfdd743006c16e3e33fab99a9b656bb1c7849118c237ca60fba75b234c187ce5
```

The [independent identity and proof audit](trinomial_cubic_audit_empty_core_sep06.md)
uses a separate
[standard-library Fraction source](../../04-computation/trinomial_cubic_audit_empty_core_sep06.py)
and [output](trinomial_cubic_audit_empty_core_sep06.out). It reads only the
literal positive coefficient arrays and denominators from the hash-pinned
primary source, without importing or executing it. Independent companion
recurrences, multiplication traces and direct determinant expansions
reconstruct the residue and all three minors.

The audit first proves cancellation of the inverse term's parameter
denominator. Assigning weight one to `g` and `tau` then bounds the three
residue coefficient degrees by `6,5,4`, and the trace entry degree by
`6+i+j`. The leading `k` by `k` determinant has degree at most
`k²+5k`, giving bounds `6,14,24`. Twenty-five exact evaluations at
distinct `g=11,...,35` certify all polynomial identities, including the
degree-six discriminant. This is a polynomial-identity proof with an
explicit degree bound, not extrapolation from tested parameters.

All **205 independent gates** pass normally and with optimization, with
identical bytes. The complete analytic, type and final-text audit is
**PASS**, including the inverse carry, `K^k` normalization, primitive
domain, root attainment and overlap with THM-4440. The separate sector
source has **55 gates**, with normal/optimized outputs also identical.

```sh
python3 -B 04-computation/trinomial_cubic_audit_empty_core_sep06.py
python3 -B -O 04-computation/trinomial_cubic_audit_empty_core_sep06.py
python3 -B 04-computation/trinomial_cubic_sector_empty_core_sep06.py
```

Independent raw hashes:

```text
audit source 06dcd740dc0f3f133451ede0de2f258c22f34e7ed28393c77fa14baceeb2fa85
audit output 319af5bd209616715c379528d99587a7d537567a3cfee29848787326db348614
sector source 0a4565ea2b852e8342502ac40b0aa75e81d4586bce2be14c0b73945ff973e177
sector output 45d4bad25282b5d8a6985a51e95a96f28ddf474940f24ebe3dd1798903774b82
```

## 7. Preserved information and the stopping boundary

The source consists of complete first/doubled weighted fibres and their
lower carry. The map removes nonzero monomials, divides by the positive
factor `K`, and reduces in the cubic quotient. It preserves the actual
rootwise doubled sign and every common-zero question. The weighted
trace matrix retains all three signs; its determinant alone loses two
sign coordinates. This explains why the third-degree step needed the
middle principal minor rather than only the trace/norm pair used before.

The old method is translated polynomial positivity; its new input here
is the three-dimensional trace form. The incoming duplication theorem
supplies a complementary geometric sector, whose exact overlap was
checked on actual cancellation roots. This is the source-to-target
connection, with neither a discarded carry nor a claimed equivalence
between unrelated mathematical objects.

General higher-channel negative definiteness and general trinomial
two-rung separation remain **OPEN**. The next useful obstruction is
whether the complete shifted-positive-minor pattern persists in a
larger fixed family or fails while noncancellation survives; no larger
family is claimed by these three minors.
