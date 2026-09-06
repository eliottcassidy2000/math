# Higher residue moments close the genuine four-anchor coefficient fibre

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
In the q=5,h=4 model, fix `(e1,e2,e3,e4)=(13,55,84,35)` and put `f=e5`.
**Every B with five nonnegative roots on this fibre has strictly negative
full carried response at every positive original first-row root. Neither
interlacer assumption is needed.** The audited
[domain theorem](long_frontier_sep06_fibre_domain.md) places this entire
class below f<22/5. The algebraic sign theorem below holds throughout
the larger coefficient interval `0<=f<=5`. There are four simple positive roots for f>0 and
three for f=0. The D-measure argument that led to the result separately
gives the exact necessary bound `f<=14sqrt(2)-16` when D interlaces B.

The remaining two-anchor model, with e3,e4 free, is **OPEN**.

## 1. Inheritance and the selected operation

The [first-phase theorem](long_frontier_sep06_anchor.md) uses low C/B residue
moments to control the smallest phase uniformly. The incoming
[effective-anchor result](continuing2_20260906_effective_anchor.md) controls
an explicit far tail in the larger two-anchor model. Neither result is an
all-root sign theorem for that model. The [actual sparse-amplitude certificate](continuing2_20260906_laurent_sparse_amplitude.md)
already proves the sign at the single genuine coefficient point f=1 for
support `(-27,1,15)`; that fixed-point sign is inherited, not new here.

The closest mechanism is the positive residue measure of an interlacer.
The canonical hostile is the linearly anchored escaping family in
[open_frontier_sep06_laurent.md](open_frontier_sep06_laurent.md). The corrected
near miss is treating an indefinite quadratic form after coefficient
elimination as a sign certificate. The underused sidecar is the higher
Hankel matrix of D/B. The live board consists of the C measure, the D measure,
original-root elimination, the genuine coefficient fibre, and complete
negative-support transport. Higher moments turn the retained fourth anchor
into an effective bound on the fifth; elimination then removes the fifth
coefficient from the response at the original zero.

This map retains the original phase and all carried coefficients. Passing to
one coefficient fibre loses freedom in e3 and e4; it supplies no uniform
answer when those two coefficients vary. No shape census supplies a proof
quantifier, and no literature-priority claim is made.

## 2. The D residue measure supplies a sharp algebraic necessary bound

For general e3=A,e4=b,e5=f, write

```text
B(v)=v^5-13v^4+55v^3-A v^2+b v-f,
D(v)=v^4-11v^3+36v^2-(5/12)A v+(1/7)b.
```

If D weakly interlaces a B with nonnegative roots, cancellation gives a
positive probability measure on those roots:
`D(v)/B(v)=sum_i w_i/(v-a_i)`, with w_i>=0 and sum w_i=1. This includes
repeated and cancelled nodes. The first moments are

```text
nu0=1, nu1=2, nu2=7,
nu3=7A/12-19,
nu4=115A/12-632-(6/7)b.
```

Consequently two general necessary constraints are

```text
A>=522/7,
2592b<=-343A^2+67788A-3157056.                         (1)
```

They follow respectively from `nu1 nu3>=nu2^2` and positivity of the
three-by-three moment determinant. At A=75 the second right-hand bound on b
is `-259/288`, directly excluding the C-only boundary `(0,0,3,5,5)`.
These necessary inequalities are not asserted sufficient.

Now fix A=84,b=35. Formal division at infinity determines nu0,...,nu8
exactly. Its five-by-five Gram determinant has the factorization

```text
det(nu_(i+j))_(i,j=0)^4
  =(f^2-15f-225)(f^2+32f-136).                        (2)
```

This determinant is nonnegative. The five summands of e4 are nonnegative
and have product f^4, so AM-GM gives `f^4<=(35/5)^5=16807<12^4`.
Thus 0<=f<12, making the first factor in (2) strictly negative. The second
factor must be nonpositive, proving

```text
0<=f<=14sqrt(2)-16<19/5.                              (3)
```

The [independently audited domain companion](long_frontier_sep06_fibre_domain_audit.md)
proves that (3)'s algebraic boundary is attained and is the exact endpoint
of the simultaneous-interlacer domain. It also classifies the larger C-only
and B-only domains. That sufficiency requires its ordered-root argument;
it does not follow from determinant positivity alone.

## 3. Retain the original row and every lower carry

Use the ordinary coefficient arrays

```text
GB=(1,13,55,84,35,f), GC=(1,12,45,56,15), GD=(1,11,36,35,5),
beta=z^-1 GB, Craw=z^-1 GC, Draw=z^-1 GD,
O=sum_j binom(14,2j+1)z^j, E=sum_j binom(14,2j)z^j,
P=O star beta,
Q=(O^2+z^-1 E^2) star(beta^2+2z Craw Draw).
```

The complete Q exponent range is -1,...,8, and its raw exponent -1
coefficient is 28. At f=0 the upper coefficients can vanish. No shift or
phase renormalization is introduced. The original equation is

```text
P(-s)=2002p_f(s),
p_f(s)=f s^4-60s^3+84s^2-10s+1/11.                  (4)
```

At any positive root of (4), solve for f in that same equation. Forming Q
before substitution gives the identity

```text
sQ(-s)=-(14/11)h(s),
h(s)=3585421125s^8-26087589000s^7-83518139925s^6
     +343030019640s^5-234760993560s^4+46232902140s^3
     -3278853435s^2+73031400s-443993.                 (5)
```

It is a complete response identity modulo the original equation, not a
replacement-root computation. This is the useful gain of the additional
anchors: the original zero removes the sole remaining coefficient.

## 4. Four phase intervals exhaust every original root

For every 0<=f<=5, the endpoint signs of p_f are as follows. Dependence
on f is affine with nonnegative coefficient, so the two endpoint values
f=0 and f=5 certify the whole interval.

| Phase interval | Sign at left | Sign at right |
|---|---|---|
| (99/10000,1/100) | positive | negative |
| (1/9,13/100) | negative | positive |
| (1,3/2) | positive | negative |

Also p_f(10)<0. For f>0, the leading coefficient is positive, so a fourth
root lies in (10,infinity). Four disjoint sign changes exhaust the degree
four and force all four roots to be positive and simple. At f=0 the degree
is exactly three; the first three intervals likewise exhaust it. In
particular there are no additional complex or repeated roots to audit.

The exact coefficient certificate proves h>0 on all four closed intervals:
for each finite interval [a,b], every coefficient of

```text
(1+t)^8 h((a+b t)/(1+t))
```

is strictly positive. This covers a<=s<b with t>=0, and the positive leading
coefficient covers s=b. Every coefficient of h(10+t) is also strictly
positive, covering s>=10. Hence (5) proves Q(-s)<0 at every original root.
The proof is uniform all the way to the f=0 boundary, including the root
escaping to infinity as f decreases to zero: the other three roots remain
bounded, while their total sum is 60/f. It does not use the previous
large effective-tail threshold.

## 5. Boundary, control, and scope

The sign theorem in Sections 3-4 holds on the entire coefficient interval
0<=f<=5, even without an interlacer assumption. Section 2 places every
D-interlacing member of the four-anchor model inside that interval.
The domain companion goes further: nonnegative B roots alone imply
0<=f<=kappa<22/5, with the exact double-root boundary retained. Thus the
sign theorem covers the whole nonnegative-root four-anchor class, including
members beyond both interlacer domains.
All relevant higher-row coefficients remain induced by the same GB,GC,GD;
there is no independent scale parameter.

The actual point f=1 recovers the full first and doubled factorial rows for
`(-27,1,15)`. Other f values are model deformations, and are not asserted to
be actual rows of that Laurent support. The f=0 degree-drop control is
essential to the quantifiers.

Dropping the coefficient bound is false even with the four other anchors:
at f=6 the original quartic still has four positive roots, but its largest
root, in `(16693/2000,41733/5000)`, has positive full response. The determinant (2) is then `-25668<0`,
so this is an exact hostile to the unrestricted coefficient extension,
not to the nonnegative-root interlacer model. Minimality is not claimed.

The next unresolved operation is to bound the eliminated response while A
and b vary under both residue measures. Positivity of a handful of moment
minors alone has not been promoted to a complete characterization.
This bounded cycle stops at the complete four-anchor fibre. Its next cheap
test is a rational rectangle around `(A,b)=(84,35)`: recompute the original
phase brackets and the four coefficient transforms uniformly on that
rectangle. The unbounded phase interval must retain its full positive
shifted-polynomial certificate. No such rectangle has yet been certified.

## 6. Exact universe, independent controls, and reproduction

The [standalone source](../../04-computation/long_frontier_sep06_residue_tail.py)
and [frozen output](long_frontier_sep06_residue_tail.out) pass **174 explicit
exact gates**, with byte-identical normal and optimized runs. The universe
contains the general first five D moments and their two displayed minors,
the specialized first nine moments and fifth determinant, the complete
first and doubled rows, the same-zero elimination, all four positive
coefficient transforms, and named controls f=0,1,5,6. Standard-library
`Fraction` arithmetic supplies every decision; there is no numerical
root solver, disabled assertion, producer import, or inferred census bound.

All fifteen factorial channels of the actual f=1 first and doubled rows
are reconstructed directly. For every named control, a second path uses
ordinary carriers `(1+u)^14 s^degree B(u^2/s)` and the corresponding C,D
carriers. It extracts coefficients 9,18,16 and verifies the original first
row and full response at rational test phases. Rational sign brackets
and interval Horner bounds certify all fifteen control roots, including
the degree-three f=0 boundary and the positive-response f=6 hostile.
The formal proof covers the continuous interval; these controls check
normalization and the boundaries rather than providing the quantifier.

The full positive transform coefficient lists and exact control brackets
are in the output. Reproduce from the worktree root:

```bash
python3 04-computation/long_frontier_sep06_residue_tail.py > /tmp/residue-normal.out
python3 -O 04-computation/long_frontier_sep06_residue_tail.py > /tmp/residue-optimized.out
cmp /tmp/residue-normal.out /tmp/residue-optimized.out
cmp /tmp/residue-normal.out 05-knowledge/results/long_frontier_sep06_residue_tail.out
```

SHA-256 manifest:

```text
source  a89ede31236a8caeaace4e89e7aa8c47efffc97a417944bf3badf6d54cc46acd
output  859a7b8a9d6610c45656ebff8b436bfe96a520df43a56325c9db048ed46a08ac
semantic e615355a57495d52795eeeac6fb7dc877eec7adae1166525a987be16f3d39828
```

The concurrent [domain companion](long_frontier_sep06_fibre_domain.md)
retains the ordered-root information lost by the moment determinant and
identifies exact admissible coefficient intervals. Its domain statements
and this note's full-response sign statement are separate proof obligations.
The [independent response audit](long_frontier_sep06_residue_tail_audit.md)
reconstructed the raw Laurent
response with ordinary dictionary convolutions, substituted f from the
original equation, and recovered all nine coefficients of h. It separately
checked all four positive transforms, all fourteen phase-endpoint signs,
the D/B series through nu8 and the fifth determinant. Its analytic review
of the measure argument, AM-GM bound, degree exhaustion, simplicity, and
lower-carry sign transfer is **PASS**. This acceptance concerns the present
sign theorem and necessary D-domain bound; the companion's stronger domain
classification has its own [independent audit](long_frontier_sep06_fibre_domain_audit.md).
The response audit additionally verifies the f=6 hostile through a direct
negative-coefficient transform of h on its entire coarse root bracket,
independently of this producer's bisection and interval-Horner path.
