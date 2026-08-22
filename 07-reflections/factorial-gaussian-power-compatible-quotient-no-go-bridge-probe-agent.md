---
source: bridge-probe-agent-2026-08-03
status: >
  RESEARCH NOTE / PROVED ELEMENTARY LEMMAS + FINITE-EXACT HOSTILE; NOT CANON.
  The factorial-to-Gaussian diagonal embedding is exact and multiplicative.
  The reverse attempt by angular Reynolds projection is disproved at the
  second power on THM-3290's polynomial.  Cross-frontier comparisons are
  typed no-factorization statements, not LRC(14), FC, GMC, or JC results.
tags:
  - factorial-conjecture
  - gaussian-moments
  - torus
  - reynolds-operator
  - quotient
  - lonely-runner
  - planar-jacobian
  - strict-transform
  - operation-response
related:
  - 01-canon/theorems/THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go.md
  - 01-canon/theorems/THM-3301-symmetry-vanishing-is-mathieu-compatible.md
  - 01-canon/theorems/THM-3018-factorial-conjecture-as-a-simplex-moment-problem.md
  - 01-canon/theorems/THM-3290-archimedes-flatness-and-the-gmc3-gvc3-counterexample-family.md
  - 01-canon/theorems/THM-3267-norm-phase-factorization-ladder-and-projective-determinant-blindness.md
  - 01-canon/theorems/THM-3285-same-ancestry-allocation-switching-horn.md
  - 01-canon/theorems/THM-3287-weighted-backbone-dominance-witness-section-and-selector-cut.md
  - 01-canon/theorems/THM-3288-maximizing-witness-lifted-walk-rational-series.md
  - 01-canon/theorems/THM-3279-affine-transverse-clutch-critical-no-go.md
  - 01-canon/theorems/THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient.md
script: 04-computation/factorial_gaussian_reynolds_power_no_go_bridge_probe_agent.py
output: 05-knowledge/results/factorial_gaussian_reynolds_power_no_go_bridge_probe_agent.out
script_sha256_lf: e3aaeecf8f3c5a120cb1aabcdfd25d0a9e911e378604e0cdf9200fbbae317bde
output_sha256_lf: 9eaab54495d512bd367855e8bd27d41355ae7c631402bd4602b16338ecd404cd
---

# Power-compatible quotients: an exact factorial/Gaussian bridge and a cross-frontier no-go

## 0. Outcome

There is a clean **one-way** bridge.  The factorial functional is exactly the
standard-real-Gaussian moment functional on the coordinatewise rotationally
invariant subalgebra.  This map is an injective algebra map, so it preserves
all powers before taking moments.  While this probe was running, that bridge
and the structural Archimedes boundary were independently promoted as
[THM-3300](../01-canon/theorems/THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go.md)
on `origin/main` at `6b0a31712`.  Section 2 below is therefore an independent
rederivation, not a novelty claim.

The tempting reverse move is false.  Haar/Reynolds averaging preserves the
expectation of each already-formed power, but does not commute with forming
powers.  On the explicit polynomial `P` of
[THM-3290](../01-canon/theorems/THM-3290-archimedes-flatness-and-the-gmc3-gvc3-counterexample-family.md),
this fails at `m=2`, with a nonzero symbolic defect and a strictly positive
Gaussian hostile moment.  Thus THM-3290's Archimedes null sequence cannot be
pushed into the factorial/radial subalgebra by angular averaging.

The same exact obstruction type clarifies three current frontiers: LRC phase
and allocation scalarizations, static witness-walk series, and planar-JC
strict-transform/resultant quotients.  In every case the quotient preserves
one observable but does not intertwine the next operation.  This is a
no-factorization theorem, not a transfer of any conjecture.

## 1. Inheritance pass and concept board

Closest proved mechanisms:

- [THM-3018](../01-canon/theorems/THM-3018-factorial-conjecture-as-a-simplex-moment-problem.md)
  gives the homogeneous factorial/simplex identity, with the full-vs-
  homogeneous correction from MISTAKE-350.
- [THM-3300](../01-canon/theorems/THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go.md)
  now proves the torus-invariant Gaussian bridge and structural Archimedes
  no-go; its independent immutable audit is still pending.
- [THM-3290](../01-canon/theorems/THM-3290-archimedes-flatness-and-the-gmc3-gvc3-counterexample-family.md)
  proves that its GMC(3) cancellation is created only after angular constant
  term and Archimedes-coordinate integration.
- [THM-3267](../01-canon/theorems/THM-3267-norm-phase-factorization-ladder-and-projective-determinant-blindness.md)
  computes exactly how projectivization loses the LRC `C_12` phase.
- [THM-3279](../01-canon/theorems/THM-3279-affine-transverse-clutch-critical-no-go.md)
  makes saturation and a unit leading coefficient load-bearing when a
  critical resultant is lifted back to an affine critical point.

Canonical hostile/corrected near misses:

- MISTAKE-215 forbids global de-factorialization from a prime-local face.
- MISTAKE-226/235 forbid identifying LRC and Gaussian sums merely because
  both are sums over kernels.
- MISTAKE-350 forbids reading the homogeneous simplex reduction as full FC.
- [THM-3287](../01-canon/theorems/THM-3287-weighted-backbone-dominance-witness-section-and-selector-cut.md)
  gives an exact static-section/physical-transition separation.

Least-used relevant sidecar: an **operation intertwiner**.  A quotient must
not only preserve the current scalar; its fibres must be stable under the
next operation.  In the present lanes that sidecar appears respectively as a
character grading, an endpoint/current update, or a primitive gradient lift
with cofactor units.

The live concept board was:

| object | representation | predicate | operation | lost coordinate | decisive test |
|---|---|---|---|---|---|
| factorial functional | radial complex-Gaussian subalgebra | all power moments | powers | none inside the image | monomial moment identity |
| THM-3290 witness | `U(1)` character grading | Gaussian null sequence | power then average | reciprocal-charge pairing | compare `R(P^2)` with `R(P)^2` |
| LRC phase/owner | full increment and allocation pieces | phase/current continuation | physical translation | scale, endpoint, owner | test quotient-fibre covariance |
| static response graph | decorated adjacency matrix | relation-walk count | adjacency powering | chronology and one-pole state | check for a physical intertwiner |
| planar-JC infinity row | primitive gradient/cofactor frame | unit Jacobian mate | strict transform/resultant | content and cofactor unit | DVR projective hostile, then saturation |

## 2. The exact factorial-to-Gaussian torus embedding

Let `(U_i,V_i)`, `1<=i<=n`, be independent standard real Gaussian pairs and
put

```text
X_i=(U_i^2+V_i^2)/2.
```

Each `X_i` is exponential of mean one, so

```text
E[X_1^a1 ... X_n^an]=a_1! ... a_n!.                 (1)
```

Consequently the substitution

```text
J:C[x_1,...,x_n] -> C[U_1,V_1,...,U_n,V_n],
J(f)=f((U_1^2+V_1^2)/2,...,(U_n^2+V_n^2)/2)          (2)
```

is an injective algebra map with image the coordinatewise `SO(2)^n`-
invariant polynomial subalgebra, and

```text
L(f^m)=E[J(f)^m]                                      (3)
```

for every polynomial `f` and every `m>=0`.  Equation `(3)` is not merely a
linear-functional coincidence: `J(f^m)=J(f)^m` is what makes the bridge
legal.  Thus the FC nullcone question is exactly the question whether the
Gaussian nullcone meets this invariant subalgebra nontrivially.  This does
not prove FC, SFC, GMC, or a converse implication.

For homogeneous `f` of degree `d`, write `X=r u`, where
`r=sum X_i` and `u` lies in the standard simplex.  Then `(3)` becomes

```text
L(f^m)=(dm+n-1)! int_Delta f(u)^m dA(u),              (4)
```

which is precisely the proved homogeneous scope of THM-3018.  The torus
quotient of the complex Gaussian sphere is the simplex; it does not add a
new homogeneous theorem.

## 3. Operation-faithful quotient lemma

The elementary lemma useful across the portfolio is the following.

**Lemma (descent criterion).**  Let `pi:X->Y` and `T:X->X`.  A map
`Tbar:pi(X)->pi(X)` satisfying

```text
pi T = Tbar pi                                           (5)
```

exists if and only if every `pi`-fibre is stable in the weak sense

```text
pi(x)=pi(x')  =>  pi(Tx)=pi(Tx').                       (6)
```

The proof is immediate: `(6)` makes `Tbar(pi(x)):=pi(Tx)` well-defined, and
`(5)` implies `(6)`.  For a linear quotient, `(6)` is exactly
`T(ker pi) subset ker pi`.

For an algebra quotient used to transport a power-moment problem, the
operation is `a |-> a^m`.  Therefore the quotient class must be a congruence
for powers.  A linear projection need not be one.

THM-363's LRC scalar gauge is the positive control for this rule.  Adding
`m*i` to a residue vector is intertwined exactly by the open-cell translation
`alpha -> alpha+s*m/n`; quotienting by scalar ramps therefore preserves the
full blocker predicate, not merely its count.  This is why the independently
audited `v_1=0` SAT classification in HYP-1823 is a legal quotient, whereas
the Reynolds move below is not.

For a compact abelian group, decompose a commutative algebra into characters,
`a=sum_chi a_chi`, and let `R(a)=a_1` be Reynolds projection.  Then

```text
R(a^2)-R(a)^2=sum_(chi!=1) a_chi a_(chi^-1).            (7)
```

The right side is precisely the reciprocal-charge return information erased
at the first projection.  Higher powers similarly sum all zero-product
character tuples.  Thus `R` preserves invariant expectations but is not an
algebra homomorphism.

## 4. Exact THM-3290 hostile: the failure occurs at the second power

Use THM-3290's variables

```text
rho=t^2+xy,
A=rho+x^2,
C=y rho^2-2x t^2 rho-x^3t^2,
P=A C^2.                                                (8)
```

Give `x,y,t` the `U(1)` weights `1,-1,0`, let `R` keep equal `x,y`
exponents, and write `q=xy`.  The exact weight support is
`{-2,0,2,4,6,8}`, independently matching THM-3301.  Exact sparse expansion
then gives

```text
R(P)=q(q-4t^2)(q+t^2)^4,                               (9)

R(P^2)=q^2(q+t^2)^8(q^2-20qt^2+24t^4),                (10)

R(P^2)-R(P)^2
 =-4q^2t^2(q+t^2)^8(3q-2t^2) !=0.                     (11)
```

For `x=u+iv`, `y=u-iv`, standard real independent `u,v,t` give
`q=u^2+v^2`.  Using `E[q^a]=2^a a!` and
`E[t^(2b)]=(2b-1)!!`, the exact moments are

```text
E[R(P)]       =0,
E[R(P^2)]     =0,
E[R(P)^2]     =3,212,537,328,000,
E[(11)]       =-3,212,537,328,000.                    (12)
```

The first two equalities agree with invariance of Gaussian expectation and
THM-3290's null sequence.  The third is also forced qualitatively because
`R(P)` is a nonzero real polynomial, so its square has positive Gaussian
integral; `(12)` supplies the exact hostile value.

The concurrent
[THM-3301](../01-canon/theorems/THM-3301-symmetry-vanishing-is-mathieu-compatible.md)
proves that an infinite-order character eigenvector cannot refute the
Mathieu conclusion and records that THM-3290's `P` mixes weights.  Equation
`(11)` is the complementary local mechanism: reciprocal nonzero weights
return at the second power and make Reynolds projection nonmultiplicative.

The first failed implication is therefore

```text
E[P^m]=0 for every m
  => E[R(P)^m]=0 for every m.                           FALSE             (13)
```

What is true is only

```text
E[P^m]=E[R(P^m)] for each already-formed power m.       (14)
```

The strongest survivor is a sequence of invariant polynomials
`R(P),R(P^2),...`, not the powers of one invariant polynomial.  The repair is
to retain the complete character grading (or the entire projected-power
sequence); neither is a factorial-conjecture input.  This proves a literal
Archimedes-to-factorial **no-go**, while leaving other genuinely radial
simplex mechanisms open.

## 5. Typed cross-frontier comparison

| lane | source -> target and map | preserved predicate | destroyed information | required sidecar | hostile |
|---|---|---|---|---|---|
| factorial -> Gaussian | `f -> J(f)` in `(2)` | multiplication and every moment | nothing inside the image | none | monomial identity `(1)` |
| arbitrary Gaussian -> radial | `P -> R(P)` | expectation of one fixed polynomial | reciprocal-charge returns, hence powers | full character grading | equations `(9)--(12)` |
| LRC phase | `q in F_169^* -> [q] in P^1(F_13)` | exactly phase parity | the other five values of the same parity and affine covariance | full `q` plus chosen norm gauge; ultimately endpoint origin/current | THM-3267: six phases twice in every direction |
| LRC allocation | `(M,R) -> B=M disjoint-union R` | the undecomposed whole cylinder and scalar mass | the exact `R-M-R` switching word | address chamber and literal contributor labels, then an endpoint/current action | THM-3285 equations `(7)` and `(11a)` |
| response sequence | physical relation data -> decorated adjacency | static relation-walk count and rationality | basepoint, lifetime, one-pole chronology | an actual update intertwiner | THM-3287: every dominance arrow has `L1=2`, physical arrows have `L1=1`; THM-3288 only powers the static adjacency |
| planar-JC infinity | gradient pair -> projective direction or saturated resultant | generic line, or common-root image under the stated unit hypotheses | content, affine-vs-infinity root, cofactor unit, integrability | primitive gradient lift, saturation and a branchwise determinant-one cofactor | DVR hostile below; THM-3279 uses leading coefficient `4` after saturation |

The JC hostile is elementary but load-bearing.  Over a DVR `R` with
uniformizer `pi` and fraction field `K`, the rows

```text
v_0=(1,0),             v_1=(pi,0)                      (15)
```

define the same point of `P^1(K)`.  But `v_0` is unimodular and extends to
an `SL_2(R)` frame, whereas for every `w=(c,d) in R^2`,

```text
det(v_1,w)=pi d
```

is not a unit.  Hence existence of a Keller cofactor does not factor through
the generic projective gradient direction.  The missing coordinate is the
content ideal `(a,b)` (equivalently a primitive lift), followed by the actual
cofactor/integrability data.

This does not contradict
[THM-2230](../01-canon/theorems/THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient.md).
THM-2230 says that **after a unit mate exists**, its global response fibre is
exactly one target-shear orbit.  It does not make a mate exist from a
projective gradient line.  Likewise THM-3279's constant leading coefficient
`4` and owner saturation are precisely what make its resultant roots lift
to affine critical points.

## 6. The actual connection and the no-factorization theorem

The LRC and planar-JC objects share one honest abstraction: both have a
projective direction together with a scale torsor that controls the target.

```text
LRC: full finite-field increment q -> projective direction [q],
     scale sidecar = norm phase / absolute affine lift;

JC:  primitive gradient row v -> generic direction [v],
     scale sidecar = content valuation / determinant-one cofactor lift.     (16)
```

Projectivization preserves incidence but neither target predicate factors
through it: THM-3267 supplies the exact LRC fibre hostile, and `(15)` supplies
the exact JC fibre hostile.  Therefore any proposed common bridge that uses
only projective directions, determinant shadows, or radial torus invariants
cannot preserve both full LRC phase transport and Keller cofactor existence.
The coordinatewise factorial torus quotient is even coarser: it intentionally
erases angular scale/phase and cannot recreate either sidecar.

This is useful positively.  A cross-area construction should work one level
up, on **primitive lifted lines** rather than on their scalar/projective
shadows.  The natural next comparison is between:

- the full `q`/chosen-norm lift over THM-3267's phase backbone, with endpoint
  origin and current attached before any average; and
- the primitive strict-transform gradient line with its content divisor and
  a normalized cofactor section attached before taking a resultant.

The common question is not whether the scalar invariants match.  It is
whether the lifted line admits a global section whose transition maps satisfy
the operation square `(5)`.

## 7. Integration recommendation and honest frontier

1. [THM-3300](../01-canon/theorems/THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go.md)
   is now proved on `origin/main`, with independent immutable audit pending.
   The exact hostile `(9)--(12)` is an independent sharpening suitable for
   that audit: it makes the failure of a Reynolds-quotient equivalence visible
   already at power two, without changing THM-3300's statement.
2. For factorial progress, search directly on the simplex/radial image of
   THM-3018.  Averaging known angular GMC(3)/SIC witnesses cannot manufacture
   an FC counterexample because it loses power compatibility.
3. For LRC, do not feed THM-3288's rational static-walk sequence into a moment
   or torus engine.  First attach THM-3285's proposed endpoint-origin fibre and
   test a literal physical intertwiner on the `R-M-R` horn.
4. For planar JC, record the content ideal and a branchwise cofactor unit at
   every strict-transform step.  A scalar resultant or projective gradient
   direction alone cannot decide Keller entry.  The provisional THM-3289
   remains outside the proof graph and is not used here.

Direct progress on LRC(14) or JC(2): none.  Indirect progress: one exact
factorial/Gaussian embedding, one exact Archimedes transfer obstruction, and
one reusable operation-descent gate that explains the sharp boundaries in
THM-3267, THM-3285, THM-3287/3288, and THM-3279.

## 8. Reproduction

Run

```text
python3 04-computation/factorial_gaussian_reynolds_power_no_go_bridge_probe_agent.py
python3 -O 04-computation/factorial_gaussian_reynolds_power_no_go_bridge_probe_agent.py
```

and compare with

```text
05-knowledge/results/factorial_gaussian_reynolds_power_no_go_bridge_probe_agent.out.
```

Both modes are byte-identical.  The companion uses exact sparse polynomial
arithmetic, integer Gaussian moments, and no random sampling or floating
point.
