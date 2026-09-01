---
id: THM-4297
title: "General Lambda-zero central and tail Keller extinction"
status: >
  PROVED RELATIVE TO THM-4103/4230/4272/4292/4294 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. On the full exact-M=12 coefficient wall
  Lambda=U+W+Z=0 with U*Z*D!=0, put A=2U+W. Then D=A^2, the raw source has
  one A23 contact, and the good-target differential has exact vertical order
  nine on the genus-seven central component. Relative to the W=0 local model
  with U_eff=A/2, the only new top-face term in the prepared equation is
  -(W/2)r^4q^3. On every repeated scale it enters only at t^6, after the
  complete four-step critical ladder and the terminal b^12 splitter. Thus
  every exceptional tail is Keller-constant as well. Proper-flat degree
  conservation excludes the entire wall. Together with THM-4290, the full
  exact-M=12 region U*Z*D!=0 is closed. The walls U=0, Z=0, D=0, seam entry,
  and JC(2) remain open.
source: root/jc2-planar-20260831
depends_on:
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4272-lambda-zero-a23-contact-and-e0-infinity-jet-obstruction
  - THM-4292-lambda-zero-repeated-face-keller-extinction
  - THM-4294-lambda-zero-central-keller-extinction
related:
  - THM-4289-a23-blowdown-observer-kahler-dualizing-quotient
  - THM-4290-exact-weight-twelve-deck-equivariant-visible-quotient-exclusion
  - THM-4293-lambda-zero-normalized-response-degree-stratification
primary_script: 04-computation/jc23_general_lambda_zero_keller_extinction_thm4297.py
primary_output: 05-knowledge/results/jc23_general_lambda_zero_keller_extinction_thm4297.out
primary_script_sha256: ae0d14f31313faff594f57619af97c1ddd09668ca91268cac27359d0a8df5202
primary_output_sha256: d58b12e1b75259c983092ab9c8e0819ecae8243a1256ff03aa50ffb88ab7ddd6
independent_audit_script: 04-computation/jc23_general_lambda_zero_keller_extinction_independent_audit_thm4297.py
independent_audit_output: 05-knowledge/results/jc23_general_lambda_zero_keller_extinction_independent_audit_thm4297.out
independent_audit_script_sha256: 33637ff886f7de1d966422921591a5f8a0f54ee0452411293f54b9b1504eac5f
independent_audit_output_sha256: 37170c6624fa25c642030bf895d81d7054c7cd0d47068d761f2c96c16d518081
hash_basis: raw LF bytes
audit: >
  PASS. A SymPy path reconstructs the complete general top face, contact,
  uniform central generic-point resultant, repeated local equation,
  moving-critical ladder, nonzero-root discriminant, deepest W-nonzero
  hostile, and all valuation inequalities. A standard-library
  sparse-polynomial path independently rebuilds the pre-t^6 equation from
  binomial coefficients and checks the exact transport and a larger order
  grid. LF-normalized normal and optimized streams byte-match frozen outputs.
---

# THM-4297 -- General Lambda-zero central and tail Keller extinction

**PROVED RELATIVE TO THM-4103/4230/4272/4292/4294 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THE FULL `Lambda=0`, `U*Z*D!=0` WALL IS EXCLUDED.
TOGETHER WITH THM-4290, THE COMPLETE EXACT-`M=12` REGION `U*Z*D!=0` IS
CLOSED. SEAM ENTRY AND `JC(2)` REMAIN OPEN.**

## 1. Statement

Work over `C` in the exact-weight-twelve reduced `(2,3)` seam. Put

```text
D=W^2-4UZ,                         Lambda=U+W+Z.               (1)
```

Assume

```text
Lambda=0,                          U*Z*D!=0.                    (2)
```

No condition is imposed on `W`; in particular the theorem includes the
previously open locus `W!=0`.

> **Theorem.** A nonautomorphic planar Keller pair cannot give an
> exact-weight-twelve response satisfying `(2)`.

The proof makes every component of a proper regular special fibre map
constantly to the good elliptic target. It does not use a wall specialization
of the interior degrees `34/42`, and therefore bypasses MISTAKE-531.

The inheritance pass is:

- closest proved mechanism: THM-4294's exact order-nine calculation on the
  central component;
- canonical hostile: THM-4292's deck-invariant prepared quadratic with a
  genuine genus-two, degree-two tail of differential order zero;
- corrected near miss: the interior response packet cannot be transported
  across a double top-edge root;
- least-used decisive sidecar: on `Lambda=0`, the entire `W`-dependence
  beyond the transverse scalar is one cubic-in-`q` perturbation.

The live concept board was

```text
{D-square gate, A23 contact, cubic q-perturbation,
 critical-value ladder, good differential, component degree}.       (3)
```

The new connection has source the exact top face, target the repeated Newton
model, map the substitution `q=t^6y`, preserved predicate positivity of the
Keller-form order, destroyed information terms of order at least `t^6`, and
sidecar the strict timing inequality for the `b^12q` splitter. The cheapest
hostile test is to force every earlier critical coefficient to zero at a
point with `W!=0`.

## 2. The general wall has the same single `A_23` contact

Set

```text
A=2U+W.                                                       (4)
```

Since `Z=-U-W` on `(2)`, one has the exact identity

```text
D=W^2-4UZ=(2U+W)^2=A^2.                                      (5)
```

Thus the gate `D!=0` is exactly `A!=0`.

THM-4230's multiplicity-one special fibre is

```text
G_0=R C,
R=S^2-P,
C=1-UP^6-WS^2P^5-ZS^4P^4.                                   (6)
```

At top infinity put

```text
b=1/S,                         r=P/S^2.                       (7)
```

The regular closure of `C` and its intersection with `R:r=1` are

```text
Ctilde=b^12-Ur^6-Wr^5-Zr^4,
Ctilde|R=b^12-Lambda=b^12.                                   (8)
```

Moreover

```text
partial_r Ctilde|(b,r)=(0,1)=-(6U+5W+4Z)=-A!=0.              (9)
```

THM-4272 proves that `C` is a smooth projective genus-seven curve throughout
`U*Z*D!=0` and that `C intersect R` is finite flat of rank twelve. Equations
`(8)--(9)` show that on the full wall the intersection is one fat point
`12Q`. The two smooth branches have intersection order twelve, so their union
is analytically

```text
q(q-unit*b^12)=0,                                            (10)
```

an `A_23` singularity. Away from `Q`, the multiplicity-one special fibre is
smooth. Every noncentral component created by normalization or regularization
is therefore either a rational toric/point-blowup component or is centered
over this one contact.

## 3. Exact order nine on the central component

Let `sigma^12=Q_target`. The coefficient-independent scaling calculation of
THM-4103/4294 gives

```text
phi^*eta_0=sigma^9 dS/G_P,                                  (11)
```

where `eta_0=dX/(2Y)` is the nowhere-zero invariant differential on the good
elliptic scheme and `G` is the exact integral source equation.

At the generic point `xi_C` of `C`, differentiation of `(6)` gives

```text
G_P=-C+R C_P=R C_P                    in k(C).               (12)
```

Both factors are nonzero. Indeed,

```text
C|_(P=S^2)=1                                             on Lambda=0, (13)
```

so `R` is not identically zero on `C`. Uniformly for every allowed
specialization, setting `S=0` gives

```text
C=1-UP^6,
C_P=-6UP^5,
Res_P(C,C_P)=46656 U^6!=0.                                 (14)
```

Hence `C_P` cannot vanish identically on the specialized central component.
The primary audit also computes `gcd(C,C_P)=1` over the generic coefficient
field as a hostile control.

Thus `G_P` is a unit in the DVR at `xi_C`, and the implicit-function chart
makes `dS` a relative differential generator. Equation `(11)` yields

```text
ord_C(phi^*eta_0)=9.                                      (15)
```

After finite ramified base change this is `9e>0`. Properness of the good
target extends the map at the generic point. Its special restriction to `C`
pulls back the nonzero target form to zero. A nonconstant morphism of smooth
characteristic-zero curves has nonzero differential, so

```text
phi_C:C->E_0 is constant.                                  (16)
```

This is the good-model differential, not a dualizing residue on the raw
`A_23` curve; THM-4289's Kahler/dualizing firewall is respected.

## 4. The exact `W`-perturbation

Use THM-4292's local variables

```text
t=sigma*b,                       q=r-1,
F=q(Hhat-b^12)-t^12/2.                                      (17)
```

The top part of `Hhat` satisfies

```text
Ur^6+Wr^5+Zr^4
 =(A/2)(r^6-r^4)-(W/2)r^4(r-1)^2.                         (18)
```

The first term is THM-4292's `W=0` top face with `U_eff=A/2`. Because
`(17)` already contains the outer factor `q`, the entire difference from
that model is

```text
Delta_W F=-(W/2)r^4q^3.                                    (19)
```

This locates the new information at one exact Newton order; it is not a
continuity argument.

## 5. Every exceptional tail remains Keller-constant

### 5.1. Nonrepeated faces

Weierstrass preparation gives

```text
F=q ell+q^2V(q,t)-t^12/2,               V(0,0)=A!=0.       (20)
```

Before a repeated discriminant root, THM-4291/4292's divisorial calculation
uses only the nonzero quadratic coefficient in `(20)`. With

```text
s=v(sigma)>0,                  beta=v(b)>0,               (21)
```

every generically separable face has Keller-form order at least

```text
3s+5beta>0.                                                  (22)
```

The possible lower-row/`b^12` cancellation for `beta<s` ends in a quadratic
with discriminant `x^2+2A`, whose roots are simple because `A!=0`. On the
balanced face `beta=s`, after lower rows of order below six vanish, the
discriminant is

```text
(c-X^6)^2+2A.                                               (23)
```

A nonzero multiple root would force both `c-X^6=0` and `A=0`. The only
possible multiple root is `X=0`, occurring precisely when

```text
c^2+2A=0.                                                   (24)
```

### 5.2. Repeated faces and the timing firewall

On `(24)`, set `q=t^6y`. Then `(19)`, divided by `t^12`, is

```text
-(W/2)t^6(1+t^6y)^4y^3.                                   (25)
```

It first appears at order `t^6`. Through order five, the exact critical
equation is therefore identical to THM-4292:

```text
-(1/2)(cy-1)^2
 +alpha_11 t y^2
 +t^2(upsilon_5 y^2+(8/3)y)
 +eta t^3y^2
 +t^4(Delta y^2-3y)
 -epsilon y.                                              (26)
```

Here `epsilon=b^12/t^6`. The four moving-critical coefficients are

```text
C_1=alpha_11/c^2,
C_2=upsilon_5/c^2+8/(3c),
C_3=eta/c^2,
C_4=(Delta+32/9-3c)/c^2.                                 (27)
```

There is no `t^5` coefficient. If `C_j` is first nonzero, its splitter has
gap `j(s+beta)` for `1<=j<=4`. If all four vanish, the `b^12q` splitter has
gap `6(beta-s)`. It always precedes every correction beginning at `t^6`:

```text
6(beta-s)<6(s+beta).                                      (28)
```

Thus the `W`-perturbation can change later normalized equations but cannot
reach the first splitter. THM-4292's order table transports unchanged:

| first splitter | Keller-form lower order |
|:---|:---|
| `b^12q` | `6s+2beta` |
| `C_1t` | `(5s+9beta)/2` |
| `C_2t^2` | `2s+4beta` |
| `C_3t^3` | `(3s+7beta)/2` |
| `C_4t^4` | `s+3beta` |

Every entry is positive. Normalization at a simple splitter cancels the
apparent derivative zero as in THM-4292, ramification scales all orders
together, and a root at `X=0` passes to the next ratio case. This exhausts
every divisorial valuation centered over `Q`. Every exceptional component
above the contact, including every positive-genus or repeated Newton tail,
therefore has constant Keller map to `E_0`.

## 6. A genuine `W!=0` deepest hostile

The transport is nonvacuous. Put

```text
c=5152/405,                    Delta=4672/135,
alpha_11=eta=Phi=0,            upsilon_5=-8c/3,             (29)
```

impose the paired cancellations, and choose

```text
U=1,
W=-13599602/164025,
Z= 13435577/164025.                                       (30)
```

Then `Lambda=0`, `W*U*Z*D!=0`, condition `(24)` holds, and all four
coefficients in `(27)` vanish. At `(s,beta)=(1,6)`, the lower gaps `14` and
`28` cancel, the first splitter is `b^12q` at gap `30`, and the Keller-form
order is

```text
6s+2beta=18>0.                                             (31)
```

This is the canonical hostile to treating `W!=0` as a generic perturbation
which automatically removes repetition. It validates the timing argument.

## 7. Complete component extinction and degree contradiction

After finite base change, normalize and resolve the exact source graph over
the good elliptic scheme. Every special-fibre component is one of:

1. the genus-seven central component `C`, constant by `(16)`;
2. the rational component `R`;
3. a fixed toric or ordinary point-blowup component, hence rational; or
4. a divisorial component centered above `Q`, constant by Section 5.

Thus every component map is constant. Let

```text
L=O_E(O),                         M=phi^*L,                 (32)
```

where `L` has relative degree one. If the special fibre is

```text
X_0=sum_i m_i X_i,                                           (33)
```

proper-flat intersection gives

```text
deg(M_generic)=sum_i m_i deg(M|X_i)=0.                       (34)
```

This uses the actual pullback bundle and retains every fibre multiplicity;
no vertical twist or dual-graph genus can absorb degree. A Keller response in
the inherited seam is nonconstant and has positive degree, contradicting
`(34)`. This proves the theorem.

Combining the wall result with THM-4290 on `Lambda!=0` gives

```text
no exact-M=12 reduced-seam candidate with U*Z*D!=0.           (35)
```

## 8. Failure boundary and generated next tasks

The result does not cross `U=0`, `Z=0`, or `D=0`. At those walls the smooth
genus-seven family, the coefficient `A`, or both can degenerate, so neither
the component inventory nor the timing proof may be imported. It also does
not prove exact-`M=12` seam entry, another reduced cell, `JC(2)`, or `DC(2)`.

The exact `t^6` boundary in `(25)` generates two typed follow-ups:

1. audit whether the later repeated-discriminant order set
   `{13,14,15,16,18}` changes when `W!=0`; its missing `17` is not
   transported by this theorem; and
2. at `D=0`, retain `A=0` rather than dividing by the quadratic coefficient
   and compute the first nonzero cubic face from `(19)`.

These sidecars are not dependencies of `(35)`.

## 9. Reproduction

Run

```bash
python3 -B 04-computation/jc23_general_lambda_zero_keller_extinction_thm4297.py
python3 -B -O 04-computation/jc23_general_lambda_zero_keller_extinction_thm4297.py
python3 -B 04-computation/jc23_general_lambda_zero_keller_extinction_independent_audit_thm4297.py
python3 -B -O 04-computation/jc23_general_lambda_zero_keller_extinction_independent_audit_thm4297.py
```

The first path uses SymPy. The second is a dependency-free sparse-polynomial
reconstruction over exact rational arithmetic. **QED.**
