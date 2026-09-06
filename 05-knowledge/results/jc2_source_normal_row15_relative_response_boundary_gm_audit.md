# Row-15 relative response over the THM-4426 boundary `G_m`

Status: **FINITE-EXACT / CANONICAL AUDIT ARTIFACT / THEOREM CANDIDATE; NOT THEOREM CANON.**

This report continues only the row-14 partial source and boundary component of
THM-4426. The exact source, frozen output, and this audit note are canonical
evidence, but no theorem identifier is promoted here. It proves a row-15
bracket/projected-depth lift after adjoining one
prefix-preserving exact-valuation-15 source response.  It does **not** complete
weights 15--22, the full source, all rows, termination, chart entry, `B_2`,
`JC(2)`, or `DC(2)`.

## Inheritance pass

- Closest proved mechanism: THM-4410/4415's frozen source-normal bracket and
  Hasse-depth operators, as replayed by the hash-pinned two-memory THM-4426
  companion (`c647bdf...fd52a`).
- Inherited solution component: THM-4426's boundary rational curve
  `Phi=eta=0`, `alpha11=1`, `c51=1087/135`, `Q(z,h)=0`, parametrized by
  `s in K*` over every characteristic-zero field `K`.
- Canonical hostile: `h=0,z=1` separates row-14 bracket from row-14 depth.
  It is not on the simultaneous boundary component and is not reused as a
  row-15 witness.
- Corrected near miss: proportionality-deduplicated source debts are adequate
  for a fixed source but unsafe after new response variables are added; a
  previously proportional family can split.  The row-15 response theorem uses
  all six raw left-nullspace compatibilities.  The deduplicated diagnostic is
  retained in the output and explicitly labelled `not_used`.
- Least-used sidecar: a source coordinate is tracked by `(weight, valuation,
  y-exponent, monomial)`.  Weight or triangular address alone loses the
  prefix-preservation and parity data used below.

## Exact row-15 bracket quotient

After the row-14 depth terminal is retained with its ten free tangent
coordinates, the row-15 bracket coefficient matrix is

```text
16 x 10, rank 10, pivots 0,...,9,
constant pivot minor 56056/390625.
```

Its raw cokernel has dimension six.  For the inherited partial source, its
nonzero source-debt span has representatives `P(z,h)` and `4Q(z,h)`.  The
second vanishes on the boundary curve.  The first is the exact quartic

```text
P(z,h)=
16615075021278830747989951720275819375429868011474609375*h^4
+2135739285093369291251960507427346972581243522949218750000*h^3*z
+110974322471541214146345210168307641503076194394432000000000000*h^3
+102724539571362544875978866741258037594474141236572265625000*h^2*z^2
+10493460117151504104257704219164846041866854743571456000000000000*h^2*z
+221844886843582624860396319525502752912380842795426811113932800000000*h^2
+2191089657459379569867864447311674997794346731962890625000000*h*z^3
+329866595380521166783896906031259893361399772421046016000000000000*h*z^2
+13572917244058089659339812857031051819140831252204016691332761600000000*h*z
+122033725213346235975852735862827724830846648816199012642548647429406720000*h
+17487435751100036783279731625048934728281943033649902343750000*z^4
+3448119974418148472951204235809436813919596354358228992000000000000*z^3
+207461284584216691729064584727505171220617425035861815587993548800000000*z^2
+3695964762068179185414288790438836972120872467162810490870369796143185920000*z
+20091716380508389394799928049046274881795912873098199870629385763077023961448448.
```

Its frozen SymPy-expression hash is
`52cd95ba3a379a4adf1ce2c5a7f9d69bc300916b6254f3d3e6cad33d18935c8b`.

Adjoin the complete exact-valuation-15 packet

```text
sum_(b=0)^7 r_b p^(15-2b)y^b.
```

It preserves every row at most 14.  Its map to the six-dimensional bracket
cokernel has rank one, with observable

```text
Lambda = 145*r0 + 30*r2 + 20*r4 + 24*r6;
r1,r3,r5,r7 are invisible.
```

Writing

```text
C=10852621164972710686787843667734315747451565056000000000000000,
```

the complete response graph is

```text
r0 = -P/C - (6/29)r2 - (4/29)r4 - (24/145)r6,             (1)
```

with `r1,...,r7` free.  Substitution of (1), followed by the exact rational
`G_m` parametrization, kills all 16 bracket coefficients identically over
`Q(s,r1,...,r7)`.

## Sharp least-weight response

For `a+2b=15`, `a,b>=0`, one has `weight=2a+3b=30-b` and `b<=7`.
Consequently there are no prefix-preserving exact-row-15 monomials of weights
15--22.  The unique weight-23 monomial is `p*y^7`; it is the invisible `r7`
direction.  The first visible monomial is therefore

```text
p^3*y^6, valuation 15, weight 24.
```

A one-coordinate section is obtained by setting all responses except `r6`
to zero and taking

```text
r6 = -145*P(z,h)/(24*C).                                  (2)
```

This has constant rational denominator before the `G_m` parametrization.
Equation (2) kills all 16 row-15 bracket coefficients.

The conclusion about weights 15--23 is deliberately relative: weights
15--22 are absent from the prefix-preserving exact-valuation-15 packet, and
weight 23 is invisible.  This does **not** rule out arbitrary lower-weight
deformations, which enter before row 15 and would require re-solving the
frozen prefix and the row-14 component.

## Projected depth

After the bracket response, the row-15 depth selector is

```text
91 x 16, rank 6, pivots 10,...,15,
constant pivot minor 1/64, terminal fibre A^10.
```

Its sole primitive source condition is exactly

```text
-4Q(z,h)=
-39181163108999707367578125*h^2
-2518217618047835370778125000*h*z
-130848131112847741467219640320000*h
-40196592278165226030060937500*z^2
-3962921571539300855745614315520000*z
-43085855533637881520123131891571941376.
```

It is independent of all seven relative responses and vanishes identically on
the inherited `G_m`.  All 91 raw depth coefficients vanish symbolically over
`Q(s,r1,...,r7)`.  Thus every base point has an `A^7` relative-source fibre
and an `A^10` row-15 terminal fibre; the weight-24 section (2) gives a
distinguished lift.

## Hostile and field quantifiers

Without an exact-valuation-15 response, pullback of `P` to the rational
parameter has zeros governed by the primitive quartic

```text
N(s)=
1969580657460775305686670067359225259066326135216064453125*s^4
-2339565479646249347106578234510498427810281609606595342336000000000*s^3
-990094381730294856177750934143522466948785694709196769644219241226112000000*s^2
-3335638504705139514152441104886050951243451910196917001116955688815072182100033536000*s
-402145106762885452045358763958282697115463254376331932971211888786341995230765606457537724416.
```

More precisely, `P(s)` is a nonzero rational multiple of
`N(s)/s^2` (the replay records denominator
`16938406087473901921821*s^2`).  Exact factorization over `QQ` leaves `N`
irreducible.  Hence **every** `s in Q*` on the inherited family is hostile to
row-15 bracket continuation without a new response.  This rational-hostile
statement is only over `Q`; an extension field may contain a root of `N`.

The positive lift has the broader quantifier: for every characteristic-zero
field `K`, every `s in K*`, and every seven relative parameters in `K`, (1)
and the constant bracket/depth pivots give the exact finite row-15 lift.

## Connection ledger

| source | target | map | preserved | lost / sidecar | decisive test |
|---|---|---|---|---|---|
| THM-4426 boundary `G_m` | row-15 bracket quotient | frozen `G_15-predicted_G_15` coefficient map | all rows through 14 and `Q=0` | quotient forgets ten tangent coordinates; retain raw six-cokernel basis | `16x10` exact rank and all 16 coefficients |
| valuation-15 source packet | bracket cokernel | left-nullspace response matrix | source valuation and monomial labels | rank-one image forgets seven directions; retain `(weight,b)` | raw `6x8` response rank one |
| bracket-neutral `A^7` | row-15 projected depth | append row-15 tangent, then Hasse-depth nullspace | every relative parameter | terminal quotient forgets six pivot tangents | all 91 residuals in `Q(s)[r]` |
| row-14 depth conic | row-15 depth debt | inherited specialization | exact equation `Q=0` | scalar normalization only | depth debt equals `-4Q` coefficientwise |
| rational parameter line | no-response hostile | pull back `P` | rational points `s!=0` | extension-field roots not excluded | irreducible quartic `N(s)` over `QQ` |

## Independent canonicalization audit

The canonicalization pass independently reran and line-audited the exact
certificate; it is not a second implementation of the row operators.

- The imported THM-4426 companion has the required raw-LF SHA256
  `c647bdf2c894d60413b041197cff446884191bf1469372e6b4a6a55ae96fd52a`.
- Fresh normal and optimized runs both exit zero and reproduce the frozen
  `18754`-byte output exactly. Verification gates remain live under `-O`;
  the script's own `check` function, rather than Python `assert`, controls the
  verdict.
- The row-15 bracket ranks are exact over the rational-function state. Their
  tangent pivot minor is the nonzero constant `56056/390625`. The response
  graph is constructed from all six raw left-nullspace compatibilities, not
  from the unsafe proportionality-deduplicated pair.
- Direct substitution verifies all `16` bracket coefficients. The projected
  depth selector has constant pivot minor `1/64`; its only source debt is
  exactly `-4Q`, and direct substitution verifies all `91` raw residuals.
- Pullback of `P` has denominator a nonzero rational constant times `s^2` and
  primitive numerator `N(s)`. Exact factorization over `QQ` returns one
  irreducible quartic factor, so the no-response hostile excludes every
  rational `s!=0` but deliberately makes no such assertion over extensions.
- Solving `a+2b=15` in nonnegative integers gives precisely `b=0,...,7` and
  weight `2a+3b=30-b`. Thus weights `15--22` are absent from this
  prefix-preserving diagonal, weight `23` is the verified invisible `b=7`
  direction, and the first visible direction is the verified `b=6`, weight
  `24` monomial. This says nothing about lower-weight deformations that alter
  and re-solve the inherited prefix.
- The positive field quantifier uses only rational identities, constant
  nonzero pivots, and the localization `s!=0`; it therefore base-changes to
  every characteristic-zero field. It does not upgrade the partial source to
  a Keller pair or establish termination, entry, `JC(2)`, or `DC(2)`.

## Replay and hashes

Run from the repository root:

```powershell
python 04-computation/jc2_source_normal_row15_relative_response_boundary_gm_referee.py
python -O 04-computation/jc2_source_normal_row15_relative_response_boundary_gm_referee.py
$env:PYTHONHASHSEED='161803'; python 04-computation/jc2_source_normal_row15_relative_response_boundary_gm_referee.py
```

The replay is hash-pinned to the canonical two-memory companion and includes
both the raw left-nullspace compatibility calculation and direct
coefficientwise substitution.  It finishes with `checks=259;result=PASS`,
including `16/16` row-15 bracket and `91/91` row-15 depth coefficients.

```text
jc2_source_normal_row15_relative_response_boundary_gm_referee.py
  SHA256 be992bd9f3550c7a23334ed6d12f770d6422426412f17d3db2863ac55ce10374
jc2_source_normal_row15_relative_response_boundary_gm_referee.out
  SHA256 a86ba2e50407989f7bc5784880c3267bbdb5d43a9710aadae25d1b119957b1cc
```

Fresh normal and `-O` runs are byte-identical to the frozen output. A fixed
`PYTHONHASHSEED=161803` replay is also prescribed as a hostile determinism
control. The canonical script and output are LF-only; the output has 84 LF /
0 CR.
