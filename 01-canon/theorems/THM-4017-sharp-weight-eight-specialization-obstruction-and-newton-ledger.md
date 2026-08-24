---
id: THM-4017
title: "Sharp weight-eight specialization obstruction and Newton ledger"
status: >
  PROVED CORRECTION + VERIFIED-EXACT + INDEPENDENTLY AUDITED; FINITE-EXACT
  TRUNCATED NEWTON LEDGER; CONDITIONAL GEOMETRIC CONSEQUENCE. At the exact
  THM-4007 sharp 5x5 cancellation, the forced coefficient
  [p^4](R/gamma)=-512/(9*A5^4) is nonzero. Under THM-4008's proposed
  weight-six scaling it has rho-exponent -2, so that family is not integral:
  the direct realization of THM-4008's six-node j=0 face, and hence the
  THM-4016 sharp-face exclusion, is REFUTED. THM-4016's arithmetic
  non-torsion theorem survives for its formal truncation point. For the
  explicitly truncated support through p^4+y^2, the exact lower subdivision
  has a j=1728 elliptic side component, toric rank 7, and abelian dimension
  1. This yields a no-isogeny obstruction only if the complete residual has
  that stable model. The first uncontrolled p*y^2 term destroys the side
  facet, so no conclusion about the full sharp residual, reduced 2:3 cell,
  or JC(2) follows.
source: root + stable-specialization gate probe + no-import audit, 2026-08-24
audit: >
  PASS. The primary exact verifier separates raw and normalized sharp
  coefficients, detects the rho^-2 pole, enumerates the two lower facets,
  checks both Newton-genus ledgers and the j=1728 quartic, and probes p*y^2.
  A no-import verifier reconstructs F_Q coefficientwise and passes 33 exact
  gates, including the shared coefficient, facet initial forms, hostile
  replacement subdivision, and normalization firewall. Normal and optimized
  runs of both scripts match their frozen LF outputs after platform-newline
  normalization.
depends_on:
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
related:
  - THM-4008-pure-p-residual-totally-degenerate-generic-fibre-no-go
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
  - THM-4016-sharp-five-by-five-elliptic-attachment-nontorsion
script: 04-computation/jc2_sharp_weight8_specialization_obstruction_thm4017.py
output: 05-knowledge/results/jc2_sharp_weight8_specialization_obstruction_thm4017.out
script_sha256: a39600e0258c37d33cee45b1e6047ce0d488aab4482578d38596584f7c4c6fe1
output_sha256: c64cdc3359aa1c744958ae007599653b08dc46faed88b54bea2e2633af71addc
independent_audit_script: 04-computation/jc2_sharp_weight8_specialization_obstruction_thm4017_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_sharp_weight8_specialization_obstruction_thm4017_independent_audit.out
independent_audit_script_sha256: e57ed6cf5a2798915556e722b350ca3dacf63d8438a3e8f29fd73b0c712b6b1d
independent_audit_output_sha256: 57da016e20da76f6a5f483837ceed9795f4875bde53a19af65c0530603b16ad4
hash_basis: raw LF bytes
---

# THM-4017 -- the forced weight-eight term changes the sharp model

**PROVED CORRECTION + VERIFIED-EXACT + INDEPENDENTLY AUDITED; FINITE-EXACT
TRUNCATED NEWTON LEDGER; CONDITIONAL GEOMETRIC CONSEQUENCE.**

Work over an algebraically closed field of characteristic zero in the
normalized reduced `(2,3)` cell. Put

```text
G=gamma*u+H(p,y),             H=lambda*p+R(p,y),
y=sp,                         A5=a^5,
Rtilde=R/gamma,               gamma=-a^3/2.              (1)
```

This theorem corrects one application without weakening THM-4008's proved
pure-`p` theorem or THM-4016's proved arithmetic statement.

## 1. The fixed sharp point has a forced weight-eight term

THM-4007's third row and the sharp fixed-support cancellation used in
THM-4016 give

```text
[p^3]Rtilde=2752/(135 A5^3),
c40=[p^4]Rtilde=-512/(9 A5^4),
c02=[y^2]Rtilde=-8128/(135 A5^3).                       (2)
```

In raw coefficients of `R`, write

```text
epsilon=[p^3]R=-1376/(135a^12),
kappa=[y^2]R=4064/(135a^12),
delta=[p^4]R=256/(9a^17).                               (3)
```

In particular, `delta!=0`. The normalized and raw pairs must not be mixed:

```text
delta=gamma*c40,                 kappa=gamma*c02.        (4)
```

THM-4008's formal weight-six scaling is

```text
q=rho^-6,              s=rho^-1 S,
p=rho^-2 P,            y=rho^-3 SP.                     (5)
```

If a residual monomial is `p^i y^j`, its contribution to
`rho^6 H(rho^-2P,rho^-3SP)` has exponent

```text
6-(2i+3j).                                                (6)
```

Thus `p^3` and `y^2` have exponent zero, but the forced `p^4` term has
exponent `-2`. The proposed family is not integral. Clearing its denominator
makes the first special equation

```text
-delta*(S^2-P)P^4=0,                                    (7)
```

not THM-4008's two-component weight-six equation. Therefore:

```text
the exact sharp 5x5 residual does not realize the
THM-4008 weight-six six-node j=0 model.                  (8)
```

This refutes the direct THM-4008-to-THM-4016 sharp-face application at the
first algebraic model gate. It is stronger than saying that stable
specialization is still unproved.

## 2. What survives from THM-4016

If one formally discards `delta*p^4` and every higher term, the two
weight-six coefficients in `(3)` give the six points

```text
X^3=-epsilon/(epsilon+kappa)=43/84,
Y^2= kappa/(epsilon+kappa)=127/84                       (9)
```

on `Y^2=X^3+1`. THM-4016 proves, by incompatible exact reduction orders at
places above `11` and `17`, that all six choices are non-torsion. That
arithmetic statement is unchanged.

What fails is the identification of `(9)` with actual attachment points in
the stable reduction of the fixed sharp residual. Equation `(8)` removes
that identification, so `(9)` excludes no actual sharp face.

This formal point is also distinct from THM-4012's actual
maximum-weight-six branch. That branch has `[p^4]Rtilde=0`, ratios
`(43/224,267/224)`, and a separately proved face-stability calculation.
Neither arithmetic certificate can be substituted for the other.

## 3. The exact conditional stable-model lemma

The geometric implication underlying the old invoice is valid with its
owners stated explicitly. Let a finite nonconstant generic morphism to a
smooth proper elliptic scheme acquire, after finite base change and graph
resolution, a source special fibre with:

1. a smooth genus-one component `D`;
2. every other component rational;
3. `D` the unique positive-genus component; and
4. all displayed branches on `D` meeting one connected rational subcurve.

Every rational component maps constantly by Riemann--Hurwitz, and constants
agree across their connected clutch. For a relatively ample target line
bundle `L`, fibrewise intersection gives

```text
deg(phi_generic^*L)
 =sum_i m_i deg((phi_i)^*(L|E_0))>0.                    (10)
```

All rational summands vanish. Hence `D` owns positive degree and maps
nonconstantly; after translation its restriction is an isogeny. The common
contracted clutch makes all displayed branch images equal.

The two load-bearing sidecars are therefore:

```text
degree owner: D carries positive degree;
clutch owner: the branches meet one connected contracted subcurve.        (11)
```

This is a stable-map argument to a fixed smooth target, not automatically a
classical admissible cover. It proves the conditional implication, but it
does not identify the stable model of the complete sharp residual. Extra
positive-genus components can steal degree, and disconnected rational
clutches need not have one image.

## 4. Finite-exact truncated Newton subdivision

There is a useful lower-face signal if one makes the explicit truncation

```text
H=lambda*p+alpha*p^2+epsilon*p^3
  +kappa*s^2p^2+delta*p^4,             delta*kappa!=0.  (12)
```

Put `Q=q^-1` and

```text
F_Q=(s^2-p)(1-QH)+gamma*Q*s^2.                           (13)
```

Give the `QH` coefficients height one and the base monomials height zero.
Exact lower-hull enumeration gives precisely two two-dimensional facets:

```text
A: z=(x+2y-2)/8,
   in_A(F_Q)=(s^2-p)(1-delta*p^4);

B: z=(x+y-2)/4,
   in_B(F_Q)=s^2(1-kappa*s^2p^2-delta*p^4).             (14)
```

After removing the monomial factor and putting `T=sp`, the nontrivial side
component is

```text
kappa*T^2=1-delta*p^4.                                  (15)
```

It is a smooth genus-one quartic. Its binary-quartic invariants satisfy
`I!=0,J=0`, hence `j=1728`. Its rational endomorphism algebra is `Q(i)`;
the target `j=0` curve has rational endomorphism algebra
`Q(sqrt(-3))`. A nonzero Hom would be an isogeny and identify these fields,
so

```text
Hom(E_1728,E_0)=0.                                      (16)
```

The genus ledger is exact. Facet A has one parabola and four vertical
rational components, eight transverse nodes, and graph rank `4`. The four
roots of the shared edge add three more graph cycles when the one elliptic
facet-B component is glued. Thus

```text
toric rank=7,              abelian dimension=1,
total genus=8.                                                (17)
```

Equivalently, the full truncated Newton polygon has eight interior lattice
points, split as `4+1+(4-1)`.

If `(12)` were the complete relevant support and the nondegenerate toric
model supplied the actual stable reduction, `(16)` would obstruct a target
elliptic quotient. That consequence is **CONDITIONAL**. Equations
`(14)--(17)` themselves are **FINITE-EXACT** facts about the named
truncation.

## 5. The first uncontrolled term destroys the signal

The same-weight monomial `p*y^2= s^2p^3` adds the lifted endpoint
`(4,3,1)`. Its gap below facet B is exactly

```text
1-(4+3-2)/4=-1/4.                                       (18)
```

It therefore destroys the `j=1728` side facet. The independently enumerated
replacement subdivision has a rational side cell and moves the positive
genus into the enlarged primary cell. This is a minimal hostile showing that
complete residual support is load-bearing.

THM-4012 studies the actual highest `(2,3)`-weighted face of the entire
polynomial. At total weight eight that face is `c*p^4+d*p*y^2` and has a
Bolza genus-two normalization when both coefficients are nonzero. That
highest-face observer and the lower `Q`-Newton truncation in `(12)--(18)` are
compatible but different constructions.

## 6. Scope and reproduction

This theorem proves:

1. the direct exact-sharp weight-six realization is **REFUTED**;
2. THM-4016's formal-point non-torsion arithmetic remains **PROVED**;
3. the stable-model implication is **PROVED** under owners `(11)`; and
4. the truncated Newton data `(14)--(18)` are **FINITE-EXACT**.

It does not prove the complete sharp residual, its stable reduction, an
actual elliptic quotient, emptiness of the reduced `(2,3)` cell, or `JC(2)`.

Reproduce both independent exact paths with

```bash
python3 -B 04-computation/jc2_sharp_weight8_specialization_obstruction_thm4017.py
python3 -B -O 04-computation/jc2_sharp_weight8_specialization_obstruction_thm4017.py
python3 -B 04-computation/jc2_sharp_weight8_specialization_obstruction_thm4017_independent_audit.py
python3 -B -O 04-computation/jc2_sharp_weight8_specialization_obstruction_thm4017_independent_audit.py
```

All four streams match their frozen LF outputs after platform-newline
normalization. **QED.**
