---
id: THM-3517
title: "Weighted odd family: all-coordinate primitivity, the m=3 quintic atlas, and sign-blind Jelonek components"
status: >
  PROVED + VERIFIED-EXACT; M3 ATLAS AND ALL-GRADE Z GATE INDEPENDENTLY
  AUDITED.  Every
  actual source coordinate x,y,z is primitive in every member of the explicit
  THM-3448 cyclic weighted family, hence all three share the inverse
  discriminant square class.  At m=3 they are explicit quintics with class
  [D5]=[L5] after target pullback, while the exact Jelonek set is
  V(C) union V(L5).  For every m>=3 in the odd subfamily, the sign class
  misses the genuine C-component because its local inertia cycle is odd and
  therefore even as a permutation.  The all-grade z proof has disjoint
  symbolic-remainder and rational two-root certificates, with exact finite
  hostiles through ell=30 and forward replays through n=256.  No
  identification with THM-1605's unstored E_m is made.
author: codex-2026-08-16
source: codex exact extension of the Gallagher/THM-3438 weighted lift
audit: >
  deterministic exact Sympy resultants, irreducibility factorizations,
  discriminants, square roots, branch/index gcds, source-map replay,
  determinant, Newton residual, m=2 hostile, and byte-identical normal/-O
  transcripts; independent python-flint reconstruction of the map, all three
  resultants, the 191-term z row, 268-term index core, pullback, and an F_31
  separator; symbolic all-n remainder proof and direct ell=1..30 hostiles for
  z-primitivity; independent fractions-only two-root proof and n=3..256
  exact forward replay; structural all-m boundary deduction from audited
  THM-3448
depends_on:
  - THM-3438-weighted-lift-keller-degree-spectrum
  - THM-3448-weighted-keller-cyclic-jelonek-inertia-family
  - THM-3494-weighted-lift-primitive-coordinate-discriminant-atlas
related:
  - THM-1605-infinite-family-extent-vs-mechanism
  - THM-2546-integral-coordinate-dichotomy-and-parity-lens-scope
  - THM-3508-level-two-sporadic-keller-three-coordinate-primitive-discriminant-square-class
  - THM-3519-level-three-sporadic-keller-three-coordinate-primitivity-and-common-discriminant-class
script: 04-computation/jc_weighted_odd_family_m3_coordinate_jelonek_probe_20260816.py
output: 05-knowledge/results/jc_weighted_odd_family_m3_coordinate_jelonek_probe_20260816.out
independent_script: 04-computation/jc_weighted_odd_family_m3_three_coordinate_independent_audit_20260816.py
independent_output: 05-knowledge/results/jc_weighted_odd_family_m3_three_coordinate_independent_audit_20260816.out
script_sha256: ff03f0ac71922f5881bf112ae2c8cbb7f5cd49198500c19d62b5eaf5b286002e
output_sha256: 6d64acd12b3b546e8820a683f38ada905f0a52af96ab3661319eb111c45699e0
independent_script_sha256: e0eaed0b1f034c799cd95a00ede7fdd60ff4c9258330b39845db45ea159e97e2
independent_output_sha256: b458c072d8220403162bd972ab32dc3c4aaacd217cfb2ebf7a0525fc0e14f96d
semantic_sha256: ed6845e743f8554327653521f243817264b08d1ca864c8513c0b2af7ce17ac81
independent_semantic_sha256: bd208dad9732439dfa14a794ef54dbfe57d66360f5f5a39b26353ffb82b6bba3
all_m_script: 04-computation/jc_weighted_cyclic_z_primitivity_all_m_probe_20260816.py
all_m_output: 05-knowledge/results/jc_weighted_cyclic_z_primitivity_all_m_probe_20260816.out
all_m_script_sha256: 0c74f0546444ef28265800e89bbcd0c1161eadfdd75599b6707428c22853424b
all_m_output_sha256: 01e1db90014ecd1f8f5fda0ecc250e23429db94bb04033bc1ab91500b4a38e4a
all_m_semantic_sha256: 6f06d5042a944f817b23e5bc43e3d70f74098ed4640e20eb95a638bec575cbf0
all_m_independent_script: 04-computation/jc_weighted_cyclic_z_primitivity_two_root_independent_audit_20260816.py
all_m_independent_output: 05-knowledge/results/jc_weighted_cyclic_z_primitivity_two_root_independent_audit_20260816.out
all_m_independent_script_sha256: b95756882ba89e51a542858807c939fd80bdd7a0d69fe3e6a7187e3318ef9cd7
all_m_independent_output_sha256: 53c5bb974b37d3f5263247f24e70b46973ce34a286c1f5de1e9cc579b428ae05
all_m_independent_semantic_sha256: 2c8d9f191ae5a8dddd6ef15feeb20d6564491bd0be6ea34c8e662428f74b729b
hash_basis: LF-normalized bytes
---

# THM-3517 -- the odd weighted family beyond the cubic

**PROVED + VERIFIED-EXACT; M3 ATLAS AND ALL-GRADE Z GATE INDEPENDENTLY
AUDITED.**  Sympy and a clean-room python-flint implementation agree on the
full all-coordinate `m=3` ledger, including the 191-term `z` eliminant and
its 268-term index core.  Disjoint symbolic-remainder and rational two-root
arguments prove `z`-primitivity in every cyclic grade.

## 1. Provenance boundary and the explicit object

The legacy Keller file carrying the colliding identifier `THM-1605` contains
no literal definition of its historically reported family `E_m`, no source
artifact from which that formula can be recovered, and no exact `m>=3`
replay.  Its degree list therefore cannot lawfully define an input to this
test.  This theorem instead uses the explicit cyclic weighted subfamily of
THM-3448, itself a specialization of THM-3438's published weighted lift.

For `m>=2`, put

```text
ell=2m-3,
p_ell(w)=(ell+1)w^ell-(ell+2)w^(ell+1),
R_ell(w)=w^(ell+1)(1-w),
q_ell(w)=ell w^(ell+1)-(ell+1)w^(ell+2),
a_m=-(2ell+1)/(2ell)=-(4m-5)/(4m-6).                 (1)
```

For source coordinates `(x,y,z)`, set

```text
u=1+xy,                  gamma=1+a_m xy+x^2z,
w=u gamma,

beta=1+(ell+1)u^ell gamma^(ell-1)
        -(ell+2)u^(ell+1)gamma^ell,
alpha=u+ell u^(ell+1)gamma^(ell-1)
        -(ell+1)u^(ell+2)gamma^ell,                  (2)

E_m^cyc=(alpha/x^2,beta/x,x gamma).                  (3)
```

The apparent quotients cancel.  THM-3438 and THM-3448 prove

```text
det J(E_m^cyc)=1,
[Q(x,y,z):Q(E_m^cyc)]=ell+2=2m-1,
Mon_geom(E_m^cyc)=S_(2m-1).                          (4)
```

The ordinary coordinate degrees and `z`-degrees follow directly from (2):

```text
deg(E_m^cyc)=(10m-13,10m-14,4),
deg_z(E_m^cyc)=(2m-3,2m-3,1).                        (5)
```

Thus this explicit odd family has first-coordinate degrees
`7,17,27,37,...`, not the historical `7,13,26,43,64,...` string in the old
THM-1605 prose.  The common `m=2` numerical degree is not an identification
of the two unstored constructions.  We write `E_m^cyc` precisely to retain
that distinction.

## 2. The first lawful beyond-cubic member

At `m=3`, one has `ell=3` and

```text
p=4w^3-5w^4,       R=w^4-w^5,
q=3w^4-4w^5,       a=-7/6.                            (6)
```

Equations (2)--(3) become

```text
u=1+xy,             gamma=1-(7/6)xy+x^2z,

E_3^cyc=
 ((u+3u^4 gamma^2-4u^5 gamma^3)/x^2,
  (1+4u^3 gamma^2-5u^4 gamma^3)/x,
  x gamma).                                              (7)
```

Exact expansion proves polynomial cancellation, coordinate degrees
`(17,16,4)`, `z`-degrees `(3,3,1)`, and determinant one.  For target
coordinates `(A,B,C)`, put

```text
P=BC,                    Q=AC^2.                       (8)
```

The source invariant `w=u gamma` obeys the irreducible inverse quintic

```text
T(w)=w^5-w^4+Pw-Q.                                    (9)
```

Its exact irreducible reduced discriminant is

```text
D5(P,Q)=256P^5-27P^4-36P^3Q-50P^2Q^2
        -2500PQ^3+3125Q^4+256Q^3.                    (10)
```

The degree is five and the global monodromy is `S5`; in particular this
inverse is not solvable by radicals.  The `C3` below is only a local boundary
inertia subgroup, not the global monodromy.

## 3. One fully displayed coordinate eliminant

In the generic target field `K=Q(P,Q,C)`, reconstruction gives

```text
gamma=P-p(w),                x=C/gamma.                (11)
```

Eliminating `w` from `T(w)` and `X gamma-C` gives the irreducible quintic

```text
E_x(X)= -C^5 +(1-15P)C^4 X
        +(-80P^2+8P+8Q)C^3 X^2
        +(-160P^3+18P^2+28PQ+50Q^2)C^2 X^3
        +D5(P,Q)X^5.                                  (12)
```

There is no `X^4` term.  Its discriminant factors exactly as

```text
Disc_X(E_x)=C^20 D5 I_x^2,                            (13)
```

where the branch-coprime index core is

```text
I_x=625P^6-18750P^5Q-64P^5+234375P^4Q^2+1840P^4Q
    -1562500P^3Q^3-59200P^3Q^2
    +5859375P^2Q^4+564000P^2Q^3+3840P^2Q^2
    -11718750PQ^5-1920000PQ^4-79360PQ^3
    +9765625Q^6+1950000Q^5+153344Q^4+4096Q^3.        (14)
```

Exact factorization gives `gcd(D5,I_x)=1`.  Hence the full square class is

```text
[Disc_X(E_x)]=[D5] in K*/K*2.                         (15)
```

## 4. All three actual coordinates

The other reconstruction relations are

```text
Cy=w-gamma,
C^2z=gamma[gamma(gamma-1+a)-aw].                      (16)
```

Taking the corresponding exact resultants with (9) gives:

| view | eliminant terms | eliminant total degree | index-core terms/degree | forced `C` power in index |
|---|---:|---:|---:|---:|
| `x` | 17 | 10 | 17 / 6 | 10 |
| `y` | 29 | 10 | 26 / 6 | 10 |
| `z` | 191 | 15 | 268 / 26 | 20 |

All three resultants are irreducible of coordinate degree five over `Q`, so
`x,y,z` are each primitive elements of the same degree-five extension.  Their
discriminants satisfy, for nonzero `J_i in Q(P,Q,C)`,

```text
Disc(E_x)=D5 J_x^2,
Disc(E_y)=D5 J_y^2,
Disc(E_z)=D5 J_z^2,                                  (17)
```

and every computed index core is coprime to `D5`.  This is the exact
degree-five continuation of the fixed map's three-coordinate common-class
phenomenon.  The mechanism is not special to cubics: once the coordinates
are primitive, their power bases are bases of one trace algebra, so basis
change squares relate the discriminants.

The independent FLINT audit reconstructs the map before touching the inverse
atlas, obtains the same `17`, `29`, and `191`-term coordinate resultants, and
reproduces all three eliminant, discriminant, and index-core hashes.  In
particular the `z` row is irreducible of coordinate degree five, its index
root has forced factor `C^20`, the remaining `268`-term core is coprime to
`D5`, and the exact discriminant ratio is a square.

As an orthogonal separator, over `F_31` the target

```text
(A,B,C)=(26,23,1),                 (P,Q)=(23,26)       (18a)
```

has inverse roots

```text
w=(8,9,12,16,18),                                     (18b)
```

and the corresponding source-coordinate values have Vandermonde
determinants

```text
(V_x,V_y,V_z)=(28,14,1) in F_31*.                     (18c)
```

Thus all three actual coordinates separate the five sheets in one lawful
good-reduction fibre.  This independently certifies the primitive-element
conclusion without trusting the large symbolic `z` factorization.

The primitivity sidecar is indispensable.  Replacing an actual coordinate by
the flat base-field view `X=P` gives exactly

```text
Res_w(T,X-P)=(X-P)^5,
Disc_X((X-P)^5)=0.                                    (18)
```

Thus neither five visible sheets nor a symmetric comparison certifies a
coordinate discriminant without a primitive-element gate.

## 5. Pullback, nonproperness, and the missed component

Substitution (8) gives

```text
D5(BC,AC^2)=C^4 L5(A,B,C),                            (19)

L5=3125A^4C^4-2500A^3BC^3+256A^3C^2
   -50A^2B^2C^2-36AB^3C+256B^5C-27B^4.              (20)
```

The polynomial `L5` is irreducible and reduced, and `L5(A,B,0)=-27B^4`, so
`V(C)` and `V(L5)` are distinct.  THM-3448 proves the exact effectivity
statement

```text
S_(E_3^cyc)=V(C) union V(L5).                         (21)
```

The two generic component ledgers are:

| component | escaping sheets | local inertia | raw `w` discriminant order | `w` index | sign sees it? |
|---|---:|---|---:|---:|---|
| `C=0` | 3 | `(123)1^2` | 4 | 1 | no |
| `L5=0` | 2 | transposition | 1 | 0 | yes |

Because `C^4` is a square,

```text
[Disc(E_x)]=[Disc(E_y)]=[Disc(E_z)]=[L5]
    in Q(A,B,C)*/Q(A,B,C)*2.                          (22)
```

Equation (22) detects the transposition component but erases the genuine
`C=0` nonproperness component.  This is not a computational accident: a
3-cycle is even, while a transposition is odd.  The sign resolvent records
the parity of local monodromy, not the full effective Jelonek divisor.

The cubic boundary is a sharp hostile.  At `m=2`, the same pullback has a
factor `C^2`, but THM-3448 proves

```text
S_(E_2^cyc)=V(L3),                                    (23)
```

with no separate `C=0` component.  Hence an even discriminant factor may be
either a genuine sign-blind nonproper component (`m=3`) or only index/chart
data (`m=2`).  Effectivity cannot be recovered from square class plus
multiplicity parity.

## 6. All-coordinate primitivity in every cyclic grade

The third coordinate admits an all-grade proof without computing another
large resultant.  Write `n=ell+2`, let

```text
T=w^n-w^(n-1)+Pw-Q,
gamma=P-(ell+1)w^ell+(ell+2)w^(ell+1),
H=gamma[gamma(gamma-1+a)-aw],                         (24)
```

so that `C^2z=H(w)`.  Let `lambda(g)` be the coefficient of `w^(n-1)` in the
remainder of `g` modulo `T`.  If

```text
c_k=lambda(w^k mod T),
```

then multiplication of `T=0` by `w^(k-n)` gives the exact recurrence

```text
c_k=c_(k-1)-P c_(k-n+1)+Q c_(k-n).                   (25)
```

For `n>=4`, the initial block has `c_k=1` for
`n-1<=k<=2n-3`.  Continuing (25) through the four powers occurring at the
top of `gamma^3` gives

```text
lambda(gamma^2)=1-n(n-2)P,
lambda(w gamma)=1,

lambda(gamma^3)
 =1+(n^3-3n^2+3n)P^2+(6-4n)P+4(n-1)Q.              (26)
```

Since `a=-(2ell+1)/(2ell)=-(2n-3)/(2n-4)` and
`H=gamma^3+(a-1)gamma^2-a w gamma`, equations (25)--(26) yield

```text
lambda(H)=((ell+1)^3+1)P^2
          +(4ell^2+ell-2)P/2+4(ell+1)Q.              (27)
```

This is a nonzero polynomial.  Hence the remainder of `H` has exact degree
`n-1`.  If `z` belonged to the target field `K=Q(P,Q,C)`, then

```text
rem_T(H)-C^2z                                         (28)
```

would be a nonzero polynomial of degree `n-1` over `K` annihilating `w`,
contradicting the irreducible degree-`n` polynomial `T`.  Thus `z notin K`.
At the exceptional small index `ell=1`, direct reduction gives the nonzero
leading coefficient `3P(6P+1)/2`, proving the same conclusion.

There is also a coefficient-free independent certificate.  Specialize

```text
P=2^(-(n-1)),                 Q=0.
```

Then `0` and `1/2` are roots of `T_n`.  If `m=n-3` and
`delta=1/(2n-4)`, direct substitution gives

```text
gamma(0)=P,                       gamma(1/2)=-mP,
H(0)=P^2(P-2-delta),
H(1/2)=-mP[mP(mP+2+delta)+(1+delta)/2].
```

For `n=3`, `H(0)<0=H(1/2)`.  For `n>=4`, `m>=1`, `P<=1/8`, and
`delta<=1/4`, so

```text
|H(1/2)|>P/2>=(16/9)|H(0)|.
```

Thus `H mod T_n` cannot be constant in any grade.  The `n=3` half-root has
`gamma=0` and is used only as a polynomial-identity witness.  For every
`n>=4`, both roots lie in the reconstruction chart and give two finite source
points above `(A,B,C)=(0,2^(1-n),1)` with distinct actual `z` coordinates.
This route does not use the recurrence (25) or coefficient formula (27).

THM-3438 gives global `S_n` monodromy, whose point stabilizer `S_(n-1)` is
maximal.  Therefore the extension has no proper intermediate field.  Combining
`z notin K` with THM-3494's `x/y` argument proves, for every `ell>=1`,

```text
K(x)=K(y)=K(z)=K(w).                                  (29)
```

All three generic coordinate eliminants consequently have degree `n`, and
their trace-form discriminants share the full inverse square class.  This is
an all-family three-coordinate theorem; only the explicit 191-term `z`
eliminant is special to the displayed `m=3` atlas.

Now restrict to the odd reindexing `ell=2m-3`.  THM-3448 gives, for every
`m>=3`, a genuine `C=0` component with

```text
2m-3 escaping sheets in one C_(2m-3) orbit,
v_C(Disc_w)=2m-2,
sign(C_(2m-3))=(-1)^(2m-4)=+1.                        (30)
```

Thus every one of the three primitive coordinate sign classes misses `V(C)`
for every `m>=3` in `E_m^cyc`.  The other generic component has transposition
inertia and odd discriminant order, so it remains visible.  If its reduced
equation is normalized as `L_m` and the discriminant unit is `u_m`, then the
common pulled-back field class is `[u_m L_m]`; the constant class must not be
dropped.  At `m=3`, the exact normalization above has `u_3=1`.

## 7. Persistence and failure ledger

| fixed cubic feature | `m=3` verdict | exact reason |
|---|---|---|
| three actual primitive coordinate views | persists for every cyclic grade | remainder (27), `S_n` maximality, and THM-3494 |
| one common discriminant square class | persists for every cyclic grade | trace-form basis-change squares |
| literal factor `-4` and class `[-L]` | fails | the quintic invariant unit is `+1`; class is `[D5]`, then `[L5]` |
| cubic eliminants | fails | all three eliminants have degree five |
| sign class equals effective Jelonek support | fails | even `C3` inertia erases the real `V(C)` component |
| one Jelonek component | fails | exact support is `V(C) union V(L5)` |
| `S3` global monodromy | becomes `S5` | THM-3438's two-root incidence proof |
| radical inverse | fails generically | `S5` is not solvable |
| odd generic fibre grade `2m-1` | persists in `E_m^cyc` | inverse degree `ell+2` |
| identification with historical THM-1605 `E_m` | unavailable | its literal formula/source is absent and its recorded degrees disagree |

Nothing here classifies arbitrary Keller maps, maps within a fixed grade, or
the unstored historical outside family.  Nothing here settles `JC(2)`, and no
LRC, tournament-current, or composition-tower consequence is asserted.

## 8. Exact companions and reproducibility

The `m=3` companion expands (7), checks both apparent cancellations and the full
Jacobian, verifies the inverse identity, factors every resultant over `Q`,
computes all four discriminants, extracts exact rational square roots of the
three ratios, checks every branch/index gcd, proves the target pullback and
Newton residual, and runs the flat-view and `m=2` hostiles.  Truth gates use
`require`, not `assert`; there is no randomness or elapsed-time field.

The independent python-flint companion imports neither Sympy nor the candidate
script.  It reproduces every map/eliminant/discriminant/index hash, the target
pullback, and adds a split `F_31` fibre where the three coordinate Vandermonde
determinants are `(28,14,1)`.

The all-grade companion derives (25)--(27) symbolically and independently
computes the complete `H mod T` remainder for every `1<=ell<=30`.  Its sharp
hostile sets `P=Q=0`, where every tested remainder vanishes even though the
generic leading polynomial (27) is nonzero.  This prevents a degenerate
single-target check from being mistaken for a generic proof.

The all-grade independent companion imports neither Sympy nor the recurrence
companion.  It proves nonconstancy by the displayed two-root specialization,
then uses exact rational arithmetic to replay both roots, reconstruction, and
the full weighted map for every `3<=n<=256`.  At `n=3` it enforces the
`gamma=0` chart boundary; for `n>=4` both finite points return the same target
and different `z` values.

Run

```bash
python -B 04-computation/jc_weighted_odd_family_m3_coordinate_jelonek_probe_20260816.py
python -B -O 04-computation/jc_weighted_odd_family_m3_coordinate_jelonek_probe_20260816.py
python -B 04-computation/jc_weighted_odd_family_m3_three_coordinate_independent_audit_20260816.py
python -B -O 04-computation/jc_weighted_odd_family_m3_three_coordinate_independent_audit_20260816.py
python -B 04-computation/jc_weighted_cyclic_z_primitivity_all_m_probe_20260816.py
python -B -O 04-computation/jc_weighted_cyclic_z_primitivity_all_m_probe_20260816.py
python -B 04-computation/jc_weighted_cyclic_z_primitivity_two_root_independent_audit_20260816.py
python -B -O 04-computation/jc_weighted_cyclic_z_primitivity_two_root_independent_audit_20260816.py
```

Each normal/optimized transcript pair is byte-identical to its stored
LF-normalized output.  The independent FLINT semantic digest is
`bd208dad9732439dfa14a794ef54dbfe57d66360f5f5a39b26353ffb82b6bba3`.
The independent all-grade semantic digest is
`2c8d9f191ae5a8dddd6ef15feeb20d6564491bd0be6ea34c8e662428f74b729b`.
