---
id: THM-4304
title: "Cubic-corner repeated first-face rationality"
status: >
  PROVED RELATIVE TO THM-4301 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  At the exact-M=12 corner D=Lambda=0 with UZ!=0, every genuine repeated-q
  irreducible first-face factor is linear in q, occurs with multiplicity
  exactly two, and has rational reduced normalization. The literal repeated
  locus consists of three quadratic-section regimes and two balanced cubic
  coefficient values U_+/-; a (1,3,16) cubic realizes the latter. Thus every
  reduced first-face carrier is Keller-constant. THM-4307 later corrects the
  local algebra to A1, closes the balanced completed-local refinements, and
  leaves the three Regime-A refinement towers open.
source: root / alternating LRC14-JC2 session, 2026-09-01
depends_on:
  - THM-4301-cubic-corner-first-face-keller-extinction
related:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4299-d-zero-square-face-elliptic-splitting-and-off-corner-extinction
  - THM-4307-cubic-corner-balanced-double-section-first-refinement
primary_script: 04-computation/jc23_cubic_corner_repeated_first_face_rationality_thm4304.py
primary_output: 05-knowledge/results/jc23_cubic_corner_repeated_first_face_rationality_thm4304.out
primary_script_sha256: 029673d009905c0e11d46b9cc02090628ac201fb2c5ab550aadc0123fb13e4f2
primary_output_sha256: 2eb91282121ae842241160834cc198faaae7080822c8dba55e0845a9636c81e7
independent_audit_script: 04-computation/jc23_cubic_corner_repeated_first_face_rationality_independent_audit_thm4304.py
independent_audit_output: 05-knowledge/results/jc23_cubic_corner_repeated_first_face_rationality_independent_audit_thm4304.out
independent_audit_script_sha256: 6059034088a8f671002ac74d0d1fcd9733e899b04147be9a015b7267932cd59b
independent_audit_output_sha256: ddb1a041510aeec7a37898d43aeff039f6d4a4ad637d57a515eb74e80987338a
hash_basis: raw LF bytes
audit: >
  PASS. The primary SymPy path reconstructs the literal
  q-support, all three valuation regimes, the balanced cubic discriminant,
  both repeated values and the exact hostile. A dependency-free quadratic-
  field path independently checks the two coefficient values, double roots,
  nontriple condition, weight ledger, and separable control.
---

# THM-4304 -- Cubic-corner repeated first-face rationality

**PROVED RELATIVE TO THM-4301 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
EVERY REDUCED FIRST-FACE CARRIER IS CONSTANT. THM-4307 LATER CLOSES THE
BALANCED COMPLETED-LOCAL TAILS; THE THREE REGIME-A TAILS, SEAM ENTRY,
`JC(2)`, AND `DC(2)` REMAIN OPEN.**

## 1. Statement

Use THM-4301's exact corner

```text
D=Lambda=0,               UZ!=0,               (U,W,Z)=(U,-2U,U),
z=1/S,                    q=P/S^2-1,            t=sigma*z.           (1)
```

Let `v` be a divisorial valuation centered at `sigma=z=q=0`, put

```text
s=v(sigma)>0,       beta=v(z)>0,       tau=s+beta,       gamma=v(q)>0, (2)
```

and let `f` be its saturated first weighted face. A **genuine repeated-q
factor** means an irreducible factor `K` which depends on `q` and satisfies
`K | gcd(f,f_q)` after torus-monomial saturation.

> **Theorem.** Every genuine repeated-q factor is linear in `q`, occurs in
> `f` with multiplicity exactly two, and its reduced projective curve has
> rational normalization. Hence its reduced first-face map to the good
> elliptic target is constant.
>
> The complete literal repeated locus is the union of the three quadratic
> regimes in Section 4 and the two balanced cubic values in Section 6.
> There are no other repeated first faces.

This theorem alone does **not** say that resolving or refining a nonreduced
section adds only rational curves. Higher rows can create later discriminant
covers. THM-4307 later closes the two balanced completed-local towers and
exhibits a positive-genus Regime-A refinement, so the qualification is
essential.

The inheritance pass is:

- closest proved mechanism: THM-4301 kills every separable first face;
- canonical hostile: THM-4301's irreducible genus-four separable cubic;
- corrected near miss: rationality of the reduced repeated carrier does not
  classify components born above its nilpotent thickening;
- least-used sidecar: `gcd(f,f_q)` locates the section, while the translated
  critical-value/discriminant coordinate retains the transverse deformation.
  THM-4307 corrects the provisional two-coordinate description below.

## 2. Degree-three lemma and rational graph

THM-4301 gives coefficient-independent monomials `Uq^3` and `-t^12/2`.
If `N` is the face weight, then

```text
N<=3gamma,                         N<=12tau.                (3)
```

Every monomial of `q`-degree at least four has weight strictly larger than
`3gamma`, hence strictly larger than `N`. Therefore

```text
deg_q(f)<=3.                                                (4)
```

Let `K` be a genuine repeated-q irreducible factor. In characteristic zero,
if `K` occurred once then reduction of `f_q` modulo `K` would be
`K_q*(f/K)`, which is nonzero because `K` depends on `q`. Thus `K^2 | f`.
Equation `(4)` gives

```text
deg_q(K)=1.                                                (5)
```

After clearing a torus monomial, write the homogeneous factor as

```text
K=a(sigma,z)q+b(sigma,z).                                 (6)
```

Genuineness and saturation give `a,b!=0`; irreducibility gives
`gcd(a,b)=1`. Projection from `K=0` to the weighted line
`P(s,beta)` is the graph of the weighted rational function `q=-b/a`.
A vertical fibre would require `a=b=0`; for two binary weighted-homogeneous
forms this would give a common factor, including at either coordinate
boundary. Hence the graph is birational to the weighted line, whose
normalization is `P^1`. This proves rationality of the reduced carrier.

Algebraically `(4)--(5)` allow a double linear factor, or a triple linear
factor when `deg_q(f)=3`. Sections 4--6 prove that the literal source permits
only the double case.

## 3. Literal coefficient rows

Use THM-4301's notation

```text
c1=alpha_11+beta_11,          c2=upsilon_5+xi_10,
c3=eta+zeta_3,                c4=Delta+Theta,
c5=Phi,                       c=7168/135-(7/6)Delta.       (7)
```

The first three coefficients above the `q`-linear row are

```text
A1=5alpha_11+4beta_11,
A2=5upsilon_5+4xi_10,
A3=4eta+3zeta_3.                                      (8)
```

They multiply `q^2t`, `q^2t^2`, and `q^2t^3`, respectively. The fixed
`q`-linear terminal row is `(8/3)qt^8`.

If neither `Uq^3` nor `-t^12/2` belongs to the face, then after removal of
the common `q` monomial its `q`-degree is at most one, so it has no genuine
repeated-q factor. The remaining three anchor patterns are exhaustive.

## 4. Regime A: constant absent, quadratic section

Suppose `Uq^3` is on the face but `-t^12/2` is not. A genuine repeated factor
requires a `q^2` term (without it the saturated face has no repeated nonzero
root); if its first order is `k`, equality of weights gives

```text
gamma=k tau,                         k in {1,2,3}.          (9)
```

After removal of the common torus factor `q`, the face is

```text
Uq^2+A_k t^k q+B_(2k)t^(2k).                              (10)
```

It is repeated exactly when

```text
A_k^2=4U B_(2k),                                          (11)
```

and is then `U(q+A_k t^k/(2U))^2`. If `12beta=2k tau`, the
term `-qz^12` also lies on the face and contributes a distinct monomial
`4Uz^12` to the discriminant. Thus repetition requires

```text
12beta>2k tau.                                            (12)
```

The three literal rows are:

| `k` | forced cancellations | nonzero transverse row | repeated equation | strict range |
|---:|:---|:---|:---|:---|
| 1 | `c1=0` | `A1=alpha_11!=0` | `A1^2=4U c2` | `s<5beta` |
| 2 | `alpha_11=beta_11=0`, `c2=c3=0` | `A2=upsilon_5!=0` | `A2^2=4U c4` | `s<2beta` |
| 3 | first two coefficient pairs zero, `c3=c4=c5=0` | `A3=eta!=0` | `A3^2=4U c` | `s<beta` |

Every factor is the rational section displayed after `(11)` and has
multiplicity exactly two because the saturated face is quadratic.

## 5. Regime B: constant present, cubic absent

Here `N=12tau<3gamma`. If the first `q^2` row has order `k`, then

```text
gamma=(12-k)tau/2,                    k in {1,2,3}.         (13)
```

No literal `q`-linear `t`-row has the required order `(12+k)/2`; the missing
integer case would be order seven. Hence the only possible `q`-linear term is
`-qz^12`. The saturated quadratic discriminant is

```text
z^24+2A_k t^(12+k),                                      (14)
```

a sum of two distinct nonzero monomials. It cannot vanish identically.
Therefore Regime B contains no genuine repeated factor.

## 6. Regime C: balanced cubic

If both coefficient-independent terms occur, then

```text
gamma=4tau,                         beta>=2s.              (15)
```

Eliminating every lower-weight `q` and `q^2` row forces

```text
alpha_11=beta_11=upsilon_5=xi_10=eta=zeta_3=Phi=0,
Theta=-Delta,                        Delta=2048/45.        (16)
```

If `beta=2s`, the `-qz^12` term is present. In THM-4301's chart
`x=sigma^2/z`, `y=q/z^6`, the cubic discriminant has constant term `4U!=0`.
Thus it is not identically repeated.

If `beta>2s`, the `z^12q` term is above the face. Put `y=q/t^4`; division by
`t^12` gives the scalar cubic

```text
p_U(y)=Uy^3+(2048/45)y^2+(8/3)y-1/2.                     (17)
```

Its exact discriminant is

```text
Disc_y(p_U)=-(820125U^2+141926400U-24696061952)/121500.   (18)
```

It vanishes exactly at

```text
U_epsilon=-315392/3645
          +epsilon*(217088 sqrt(265))/18225,   epsilon=+1,-1. (19)
```

At `U_epsilon`,

```text
gcd(p_U,p_U')=(256y+15+epsilon*3sqrt(265))/256.           (20)
```

The gcd is linear, so the factor is double and the remaining root is simple.
A triple root is impossible: comparison with `U(y-r)^3` forces
`r=9/16`, `U=2048/729`, and then the `y^2` coefficient would be
`-128/27`, not `2048/45`.

## 7. Exact hostile and positive control

Take

```text
(s,beta,gamma)=(1,3,16),
U=U_+,                 W=-2U,                 Z=U,
```

and impose `(16)`. The five relevant weights are

```text
v(Uq^3)=v((2048/45)q^2t^4)=v((8/3)qt^8)=v(t^12/2)=48,
v(qz^12)=52.                                               (21)
```

Thus the face is `t^12p_(U_+)(q/t^4)` and its repeated section is

```text
256q+(15+3sqrt(265))t^4=0.                                (22)
```

It is a literal, nonempty repeated first face with rational reduced
normalization. The same coefficients with `U=1` give the separable control

```text
Disc_y(p_1)=24553315427/121500!=0.                         (23)
```

## 8. Refinement firewall

The theorem classifies the **reduced first-face carrier** only. Translating an
actual repeated section to `Q=0` gives a double root with a distinct simple
factor. After that factor is made a unit, its intrinsic local algebra is

```text
T^1=C[[Q]]/(partial_Q Q^2)=C[[Q]]/(Q)=span{1}.             (24)
```

The earlier two-dimensional quotient `C[[Q]]/(Q^2)=span{1,Q}` is the A2
algebra of a triple root `Q^3`; Section 6 excludes that multiplicity. Raw
constant and linear Taylor rows remain computational data, but critical
recentering combines them into one discriminant coordinate. THM-4307 carries
out that repair: the balanced coordinate has a smooth formal graph and only
rational completed-local refinements, whereas an exact Regime-A section can
produce a genus-one cover. Thus the present theorem proves constancy of every
reduced literal first-face component but does not by itself close later
refinements.

## 9. Reproduction and scope

Run

```bash
python3 -B 04-computation/jc23_cubic_corner_repeated_first_face_rationality_thm4304.py
python3 -B -O 04-computation/jc23_cubic_corner_repeated_first_face_rationality_thm4304.py
python3 -B 04-computation/jc23_cubic_corner_repeated_first_face_rationality_independent_audit_thm4304.py
python3 -B -O 04-computation/jc23_cubic_corner_repeated_first_face_rationality_independent_audit_thm4304.py
```

The universe is the literal exact-`M=12` corner `(1)` over `C`. The result
does not cross `U=0` or `Z=0`, classify later refinements, complete the
component inventory, prove seam entry, or prove `JC(2)` or `DC(2)`. See
THM-4307 for the balanced completed-local classification and the surviving
Regime-A boundary.

**QED.**
