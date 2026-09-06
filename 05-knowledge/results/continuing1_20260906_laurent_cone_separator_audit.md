# Independent referee: the actual all-integer Laurent phase cone obstruction

**Status: independent analytic audit PASS, with an independent exact verifier.**
The corrected main report is accepted. Its conclusion is an obstruction to
a specified sign-certificate cone on the genuine support `(-21,1,12)`.
It neither refutes Laurent noncancellation nor proves a new general
two-return theorem. The optional half-integer extension is outside this
main audit's scope.

## Source normalization and inherited scope

I first read the proof and certificate, reconstructed the objects independently,
and ran the referee successfully before reading the producer source. The
nearest proved inputs are the coupled-window reduction in
`05-knowledge/results/creative_20260906_laurent_bridge.md` (relative to
**THM-4440**, `01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md`),
the complete midpoint identities in
`05-knowledge/results/overnight7_20260906_laurent_midpoint_transport.md`,
and the alpha completion retained in the main report. The correction lineage
is essential: arbitrary separate carrier replacements do not preserve the
selected zero, and coefficientwise positive additions need not have the
desired sign at a negative phase.

For `h=1,2,3`, I enumerate all nonnegative count triples at every positive
mass through `2g`, with `g=3h+2`, and test the literal charge equation for
`(-6h-3,1,3h+3)`. All three charges are `1 mod g`; hence a positive return
mass is divisible by `g`. The triple `(1,3h,1)` witnesses mass `g`.
The complete first fibers are `(1+j,3h-3j,1+2j)`, `0<=j<=h`, and doubled
fibers are `(2+j,6h-3j,2+2j)`, `-1<=j<=2h`. Every multinomial coefficient
matches the report, including the lower carried term.

For general nonzero source coefficients `(a,b,c)`, these counts give the
common first factor `a b^(3h)c` and phase `t=a c^2/b^3`; the doubled factor
is its square. In the report's specialization `a=b=1`, the identities are
exactly `CT(f^g)=c P_raw(c^2)` and
`CT(f^(2g))=c^2 Q_raw(c^2)`. There is no independently adjustable relative
coefficient in the mixed midpoint payment.

My independent path compiler counts north/east paths by a rectangular dynamic
program, retaining the complete beta support beginning at index `-1`.
It constructs `H,H_C,H_D` and verifies the polynomial identities

```
[u^k]H=-s P_raw(-s),
s Q_raw(-s)=s^(-1)[u^(2k)]H^2-2[u^(2k-2)]H_C H_D.
```

In particular the factor `s` inside the unnormalized mixed payment is correct.
The cubic first polynomial is `-330(s^3-28s^2+14s-1/3)`, at the original
multinomial scale. The actual carried target is the certificate's full degree
seven polynomial; its constant coefficient is `-22`, not a deleted endpoint.

## Coupled windows and the actual quotient

The referee obtains each basic window by literal symbolic differentiation
of `H` and ordinary polynomial squaring. It obtains every two-sided window
by reversal, differentiation, reversal, and differentiation with the actual
ambient degrees. A separately inverted paired-grid matrix reconstructs its
nonnegative coordinates in the basic kernels. All **147** windows for the
three controls pass, including all **84** windows and seven basic rays for
the cubic support. Equality is checked in the genuine first-row quotient;
the discarded central term is a multiple of its equation.

The report originally described the alternating displayed beta polynomial
as negative-rooted. Root identified this wording error; the corrected text
now distinguishes the unsigned negative-rooted path polynomial from its
alternating reversal, which has positive roots. I read and accept the repair.
An independent exact Sturm count confirms every displayed beta root is
positive in all three controls. Therefore its `u^2/s` pullback is real-rooted
for `s>0`, with nonzero endpoints, as required for the inherited strict
window theorem. No producer coefficient or certificate changed.

All quadratic constant-cone separation values, its positive affine phase
repair, the cubic remainders, and the signed quadratic repair match exactly.
The nonnegative Laurent cone means finite sums with nonnegative real
coefficients and arbitrary **integer** powers of `s`. It does not mean all
functions positive at the three roots. This distinction is mathematically
substantive, since the latter class contains the reported signed repair.

## All-integer proof and independent root certificates

The three disjoint positive rational intervals have opposite signs of the
monic cubic at their endpoints. Their three distinct roots exhaust degree
three, so each is simple and there are no further roots. The middle interval
lies strictly inside `(1/3,1)`.

Rather than repeat the producer's interval Horner code, the referee changes
variables from each interval to `[0,1]` and computes exact Bernstein
coefficients. Their strict signs certify `p'` with signs `(+,-,+)`, the
Lagrange numerator `N<0`, all seven quotient ray values `<0`, the actual
target `<0`, and the signed multiplier `>0` throughout the respective
intervals. Reduction modulo `p` preserves their values at its roots.

Applying the coefficient functional to the Lagrange numerator gives exactly
`N(z)=y0(z^2-28z+14)+y1(z-28)+y2`. Thus its spectral weights are
`lambda_i=N(r_i)/p'(r_i)` with signs `(-,+,-)`. Invertibility of `s` in the
quotient follows from its nonzero constant term and is also checked by exact
polynomial inversion. Every integer power, including a negative one, therefore
has the stated spectral interpretation.

For each ray, the normalized sequence has the exact form

```
f_d=e_d/beta^d=A(alpha/beta)^d-B+C(gamma/beta)^d,
A,B,C>0.
```

Its consecutive second difference is strictly positive for all integer `d`.
The reconstructed four levels give `e0,e1>=0`, `e(-1)>3e0`, and `e2>e1`.
Since `beta>1/3`, `f(-1)=beta e(-1)>e0=f0`; this remains strict if `e0=0`.
Since `beta<1` and `e1>=0`, `f2>f1`. Increasing first differences then show
that every index outside `{0,1}` has value larger than one of these two
nonnegative central values. Hence **every** shifted basic ray has
nonnegative functional value. This is an analytic all-integer implication;
the referee's independent matrix-power bank `-9<=d<=9` only checks arithmetic.

The target has functional value `-8569669792006914054574`. Since `ell(J0)=0`,
the same value separates `T-J0`. A linear functional is continuous in this
three-dimensional real quotient, so the same separation excludes the
closure of the specified cone. No finite power enumeration is being
substituted for this argument, and no arbitrary-root-function closure is
being asserted.

The exact table also implies that among the basic integer-shifted rays,
the only zero evaluations are `J0` and `s J6`; all others are strict.
This observation is not needed for the main result.

## Boundary controls and reproduction

The independent verifier checks the true lowering identity, all three mixed
pieces and their complete sum, and the nonzero selected coefficient after
deleting the second piece. It reproduces the hostile
`H=(1+u)(1-2u)^2`: the selected coefficient of `H` vanishes, while that of
the squared divided carrier is `24>0`. A mixed piece crossing the separating
wall does not by itself become an admissible sign certificate.

The verifier imports SymPy but no producer code, repository theorem engine,
numerical root solver, or optimization solver. Its complete finite universe
is three source supports, **27** literal return fibers, **147** coupled
windows, and the frozen cubic/quadratic certificate. Every gate remains active
under Python optimization. Source and JSON are resolved beside the verifier;
the pinned producer transcript additionally supports its eventual
`05-knowledge/results` location.

```
python -B 04-computation/continuing1_20260906_laurent_cone_separator_audit.py
python -B -O 04-computation/continuing1_20260906_laurent_cone_separator_audit.py
```

Frozen producer inputs:

```
source 5bd98bc1fe072fcfc88b34af46b2310e1a4e184dcefdff6c5a508a820d82d00b
output 207b86469f1b5e8b5bb18a5cd1cd2a931c4c0c596a5fe66ba354923fe3069fee
JSON   4b1ee5770b484e4164e692fbf2934f4099800b0b85d379ce01d8afef71040cc0
```

The two final runs pass **646 always-active gates** and have byte-identical
LF output. Frozen referee hashes:

```
source 57a90b04ead8c7ae2c068719590e815721f8ef4359f328bf85f1e70af80a0e9b
output 964631c09264dabe9c85371d2e24b2de22eb410b0ce28a390b9ae438d74018cd
```

No repository or Git edits were made in this independent audit.
