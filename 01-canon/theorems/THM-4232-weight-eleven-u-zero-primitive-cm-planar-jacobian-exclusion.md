---
id: THM-4232
title: "Weight-eleven U-zero primitive-CM planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3992/3997/4012/4045/4218/4222
  + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In the inherited b=d=0
  reduced (2,3) seam, the complete exact-M=11 U=0 wall with
  A*B*Z*(A+B) nonzero contains no nonautomorphic planar Keller pair.
  Together with THM-4222 this closes the entire exact-M=11 coefficient
  chart A*B*Z*(A+B) nonzero for arbitrary U. The four remaining named
  coefficient walls, other cells, seam entry, JC(2), and DC(2) remain OPEN.
source: codex-planar-jacobian-frontier-synthesis-20260826
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
  - THM-4045-live-two-three-max-seven-hidden-elliptic-tail-no-go
  - THM-4218-exact-weight-ten-hidden-elliptic-tail-degree-three-planar-jacobian-exclusion
  - THM-4222-dense-weight-eleven-primitive-cm-planar-jacobian-exclusion
related:
  - THM-4217-complete-mixed-off-antidiagonal-delta-zero-planar-jacobian-exclusion
  - THM-4220-weight-ten-zeta-zero-genus-two-planar-jacobian-exclusion
  - THM-4227-two-outsider-scale-separated-depth-eight-haar-wedge
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
external: >
  Tim Dokchitser, "Models of curves over DVRs," arXiv:1807.00025v2,
  Definitions 3.7, 3.9, and 3.12 and Theorem 3.14, through the audited use
  in THM-4045, supplies the general face/edge model and rational toric chains.
  J. S. Milne, "Complex Multiplication" course notes (2020), Sections
  1.9--1.10 and Proposition 3.13, supplies the standard theorem that primitive
  CM-pairs classify simple CM abelian varieties up to isogeny;
  https://www.jmilne.org/math/CourseNotes/CM.pdf. The exact support, faces,
  CM characters, regular-model ledger, and planar-Jacobian consequence are
  proved here.
script: 04-computation/jc23_weight11_u0_primitive_cm_exclusion_thm4232.py
output: 05-knowledge/results/jc23_weight11_u0_primitive_cm_exclusion_thm4232.out
script_sha256: 23ca63d6004b8e9b3be2b455a19c336c9598dfff3c7cf1e843615b6fc806ad37
output_sha256: 1fbdbefae40e75e0223fffbba0b589e8135d3448bae10d17cad6ec2987b4943d
independent_audit_script: 04-computation/jc23_weight11_u0_primitive_cm_exclusion_thm4232_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_weight11_u0_primitive_cm_exclusion_thm4232_independent_audit.out
independent_audit_script_sha256: 96f13d674dd9de18f107566ab80dc7a3512c5de7ba6579ce5d1e17db17f5c1ad
independent_audit_output_sha256: 7e22887a28255c3d4b6cb4ae083ee3e7f489900697201647ffbbd01f66ce7f9a
hash_basis: raw LF bytes
audit: >
  PASS/ACCEPT. The primary executes 12,657 exact gates and the standard-library
  clean-room audit executes 12,695. Normal, optimized, fixed-hash-seed, and
  combined runs are deterministic and match the frozen LF transcripts. Both
  include the fixed residual point (2,0,1), the 8,192/4,096 hostile censuses,
  all face/edge/toric/chart/CM/genus ledgers, and the degree-zero inventory.
---

# THM-4232 -- weight-eleven `U=0` primitive-CM planar Jacobian exclusion

**PROVED RELATIVE TO THM-3992/3997/4012/4045/4218/4222 + VERIFIED-EXACT
+ INDEPENDENTLY AUDITED; JC(2) REMAINS OPEN.**

## 1. Statement and scope

Work over `C` in the inherited `b=d=0` reduced `(2,3)` seam and put

```text
A=[p^4y]H, B=[py^3]H, U=[p^5]H, Z=[y^3]H.
```

> **Theorem.** The complete exact-weight-eleven wall

```text
U=0,                     A*B*Z*(A+B) != 0
```

> contains no nonautomorphic planar Keller pair.

The quantifier includes every allowed lower coefficient. In particular, it
includes both `Delta != 0` and `Delta = 0`; no condition on
`Phi,Theta,eta,Xi,K` is added beyond the inherited identity
`K=2848/45-(7/6)Delta`.

> **Corollary (complete exact-M=11 coefficient chart).** THM-4222 proves the
> same conclusion when `U!=0`. Hence their union closes
>
> ```text
> A*B*Z*(A+B) != 0,                 U arbitrary.
> ```

Thus `U=0` is no longer a remaining wall in this chart. The only named
coefficient walls left are `A=0`, `B=0`, `Z=0`, and `A+B=0`. Other cells,
seam entry, `JC(2)`, and `DC(2)` remain open.

## 2. Inheritance and overlap audit

- Closest proved mechanism: THM-4222,
  `01-canon/theorems/THM-4222-dense-weight-eleven-primitive-cm-planar-jacobian-exclusion.md`.
  Its main genus-five primitive-CM component is unchanged here.
- Canonical hostile: THM-4218's hidden elliptic tail.  It does not recur: both
  replacement faces below are rational.
- Corrected near miss: MISTAKE-487 requires the highest entire weighted support;
  the certificate enumerates the full exact-M=11 universe and all lower rows.
  MISTAKE-519 requires scope archaeology before promoting a recomputed wall.
- Exact search of current canon found no theorem already closing `U=0`.
  Before this theorem, THM-4222 section 7 explicitly listed `U=0` as open.
- Harmless inherited certificate omission: THM-4222's `expanded_support` helper
  omitted the fixed residual point `(2,0,1)` coming from `-Q*s^2/2`.  It lies
  strictly above every face and does not affect THM-4222 or this theorem.
  The new certificates include and gap-check it.
- Incoming THM-4227's scale-separated Haar wedge supplies a useful
  cross-frontier *method comparison*: split a parameter wall by its first
  surviving owner and keep the object fixed through the scale change. It
  supplies no algebraic map, genus invariant, or proof dependency here.

## 3. Complete source and coefficient collision ledger

With `epsilon=-1376/135`,

```text
H=-3p+(8/3)p^2+epsilon p^3+K y^2+Phi p^2y
  +Delta p^4+Theta p y^2+eta p^3y+Z y^3
  +Xi p^2y^2+A p^4y+B p y^3,
K=2848/45-(7/6)Delta.
```

The thirteen possible inherited rows through weight eleven are

```text
p,p^2,p^3,y^2,p^2y,p^4,py^2,p^3y,y^3,
p^5,p^2y^2,p^4y,py^3.
```

For `F_Q=(s^2-p)(1-QH)-Qs^2/2`, a row `c p^i y^j` contributes endpoints

```text
(j+2,i+j,1) with coefficient -c,
(j,i+j+1,1) with coefficient +c.
```

Together with `(2,0,0),(0,1,0),(2,0,1)`, the universal valued support has
24 points.  The active supports have 23 points for `Delta != 0` and 22 for
`Delta=0`.  The complete coincident-point coefficient ledger on `U=0` is

```text
(2,3,1): K-epsilon = 1984/27-(7/6)Delta,
(2,4,1): Theta-Delta,
(2,5,1): Xi,
(3,4,1): Z-eta,
(3,5,1): B-A.
```

The first four may disappear harmlessly; `(3,5,1)` is only an edge midpoint,
so `A=B` remains harmless.  The independent gate needed at the main nodes is
`A+B != 0`, not `A-B != 0`.

## 4. Exhaustive lower hull

The two cases exhaust the wall because the first surviving pure-`p` owner is

```text
Delta != 0:  Delta p^4;
Delta = 0:   epsilon p^3, with epsilon=-1376/135 != 0.
```

The common main and tail planes are

```text
M11: nu=(r+2l-2)/11,
T3:  nu=(r-2)/3.
```

The replacement plane is

```text
V4: nu=(-r+l-1)/4                 when Delta != 0,
V3: nu=(-2r+l-1)/3                when Delta = 0.
```

For a retained row `(i,j)`, the two exact gaps are

```text
M11: ((11-2i-3j)/11, (11-2i-3j)/11),
T3:  ((3-j)/3, (5-j)/3),
V4:  ((7-i)/4, (4-i)/4),
V3:  ((8-i+j)/3, (3-i+j)/3).
```

All are nonnegative in their respective active universes.  The fixed residual
point `(2,0,1)` has gaps `(1,1,7/4)` in the `V4` case and `(1,1,8/3)` in the
`V3` case.  Exact plane enumeration then survives every simultaneous lower-row
and aggregate-collision deletion:

```text
Delta != 0: 2^8 * 2^5 = 8,192 patterns;
Delta = 0:  2^7 * 2^5 = 4,096 patterns;
total:                    12,288 patterns.
```

## 5. Faces, polygons, edge schemes, and packets

The face equations, up to torus monomials, are

```text
g_M=(S^2-P)(1-A S P^5-B S^3P^4)=R*C,
g_T/S^2=1-(SP)^3(Z+B P),
g_V/P=-1+c_k P^k+A S P^5,

(k,c_k)=(4,Delta) or (3,-1376/135).
```

`T` and `V` are rational because their displayed cores are linear in `P`
after `T0=SP`, and linear in `S`, respectively.  The exact polygon ledgers are

| case/face | polygon | `(2Area,boundary,interior)` |
|---|---|---:|
| common `M` | `(0,1),(2,0),(5,4),(1,6)` | `(33,5,15)` |
| common `T` | `(2,0),(5,3),(5,4)` | `(3,5,0)` |
| `V4` | `(0,1),(1,6),(0,5)` | `(4,6,0)` |
| `V3` | `(0,1),(1,6),(0,4)` | `(3,5,0)` |
| global `Delta!=0` | `(0,1),(2,0),(5,3),(5,4),(1,6),(0,5)` | `(40,12,15)` |
| global `Delta=0` | `(0,1),(2,0),(5,3),(5,4),(1,6),(0,4)` | `(39,11,15)` |

Both outer packets are

```text
(10,10,5,4,2,2,2,1), full degree 36, defect 28.
```

The unchanged cubic carrier has three index-two branches, so the finite/full
response ledger remains `(30,36)`.  The vertical affine endpoint polynomial
has degree four or three, but lies on the omitted `s=0` affine divisor and is
used only for generic smoothness, not as a new response invoice.

All six outer plus two internal edge schemes are

```text
X-1, 1-ZX^3, -Z-BX, (X-1)(AX+B),
A+c_kX, c_k-X^k, AX-1, 1-BX.
```

Their discriminants are

```text
1, -27Z^2, 1, (A+B)^2, 1,
-256 c_4^3 or -27 c_3^2, 1, 1.
```

Thus exactly `A*B*Z*(A+B)*c_k != 0` makes every edge reduced and avoids
corners; in the `Delta=0` case `c_3` is fixed nonzero.

## 6. Exact common regular model

Take `Q=sigma^132`.  The face heights / primitive normals and exact charts are

| face | height | primitive normal | chart | scale on `F_Q` |
|---|---|---|---|---:|
| `M` | `12r+24l-24` | `(12,24,-1)` | `s=sigma^-12 S, p=sigma^-24 P` | `sigma^24` |
| `T` | `44r-88` | `(44,0,-1)` | `s=sigma^-44 S, p=P` | `sigma^88` |
| `V4` | `-33r+33l-33` | `(-33,33,-1)` | `s=sigma^33 S, p=sigma^-33 P` | `sigma^33` |
| `V3` | `-88r+44l-44` | `(-88,44,-1)` | `s=sigma^88 S, p=sigma^-44 P` | `sigma^44` |

Every face has multiplicity one.  The outer lifted-edge gcds, equal to the
planar lengths, are `(1,3,1,2,1,4)` and `(1,3,1,2,1,3)`; both internal gcds
are one.  The outer slopes are

```text
Delta != 0: (12,-44,-44,-12,-33,-33),
Delta = 0:  (12,-44,-44,-12,-44,-88).
```

The exact determinant-one internal chains are

```text
Delta != 0: AE -24,...,-33 (8 intermediates),
Delta = 0:  AE -24,...,-44 (19 intermediates),
both:        BD -36,...,-44 (7 intermediates).
```

Every inserted toric component is rational and multiplicity one.  In the main
chart the branches `R,C` meet at

```text
P=S^2, 1-(A+B)S^11=0,
```

with determinant `-11(A+B)S^10`.  There are eleven distinct transverse nodes,
and locally

```text
UV=sigma^132/2.
```

Hence each is an `A_131` smoothing, inserting 131 rational multiplicity-one
curves.  The edge schemes cover all compactified boundary points; the torus
derivative determinants are `-11AB`, `A P^5`, and `-B T0^3` for `C,V,T`.

After contracting rational paths, the core has vertices `R,C,T,V`, eleven
`R--C` paths, one `C--T` path, and one `C--V` path.  Thus

```text
b1=13-4+1=10, special genus=5+10=15,
```

matching both global Pick ledgers.

## 7. Primitive CM and degree obstruction

The only positive-genus normalization is the unchanged component

```text
C: 1-A S P^5-B S^3P^4=0.
```

With `x=A S P^5`,

```text
P^11=(B/A^3)x^3/(1-x).
```

It is the same genus-five cyclic degree-eleven curve as in THM-4222.  Its
branch residues are `(3,10,9)`, its Newton interior points are

```text
(1,2),(1,3),(1,4),(2,3),(2,4),
```

and both Newton differentials and Chevalley--Weil give CM type
`{4,5,8,9,10}`.  Its unit stabilizer is `{1}` and its full `H^1` contains
every nontrivial character once.  The cited Milne primitive-CM theorem makes
`J(C)` simple; dimension five then gives `Hom(J(C),E0)=0` for every elliptic
curve `E0`.  All other components are rational.

Since `q=sigma^-132`, scale the inherited target by

```text
A_target=sigma^-44 X, C_target=sigma^-66 Y.
```

Its special fibre is the smooth elliptic curve `Y^2=X^3+1`.  Every source
component, including any additional blow-up exceptional curve, maps
constantly to it.  Proper-flat degree conservation therefore gives special
degree zero, contradicting the positive degree of a hypothetical finite
generic Keller morphism.

## 8. Complete exact-M=11 synthesis and status firewall

The theorem above and THM-4222 partition the coefficient chart by the
exhaustive alternatives `U=0` and `U!=0`. Therefore the union corollary in
Section 1 is exact: under `A*B*Z*(A+B)!=0`, no further `U` stratum remains.

- **PROVED/CITED inherited:** the seam reductions and target equation; the
  Dokchitser face/edge model input as audited in THM-4045; the primitive-CM
  implication from Milne as used in THM-4222; proper-flat degree conservation.
- **VERIFIED-EXACT maintained:** every monomial/support/collision, plane, polygon,
  Pick, packet, face/edge polynomial, discriminant, primitive normal, lifted
  edge gcd, slope/chain, side chart, main local equation, CM-character, target
  scaling, genus, and degree-inventory calculation stated above.
- **PROVED RELATIVE:** the `U=0, A*B*Z*(A+B)!=0` exclusion and, with
  THM-4222, the arbitrary-`U` union corollary.
- **OPEN:** `A=0`, `B=0`, `Z=0`, `A+B=0`, other M11 cells, seam entry, JC(2),
  DC(2), and the unresolved M12 hidden-Hom locus/walls.
- **Later orthogonal M12 advance:** THM-4230 promotes the Prym analysis and
  proves that every hypothetical point in its exact-M12 gate lies on a
  countable proper, nonempty hidden-`E_0` Hom locus. It repairs the naive
  `Hom=0` transfer by saturating the visible integral Hom lattice and keeping
  the connected Prym. It is not a dependency of the present M11 theorem and
  does not close M12 or JC(2).

## 9. Replay and audit boundary

```bash
python3 -B 04-computation/jc23_weight11_u0_primitive_cm_exclusion_thm4232.py
python3 -B -O 04-computation/jc23_weight11_u0_primitive_cm_exclusion_thm4232.py
PYTHONHASHSEED=4232 python3 -B 04-computation/jc23_weight11_u0_primitive_cm_exclusion_thm4232.py

python3 -B 04-computation/jc23_weight11_u0_primitive_cm_exclusion_thm4232_independent_audit.py
python3 -B -O 04-computation/jc23_weight11_u0_primitive_cm_exclusion_thm4232_independent_audit.py
PYTHONHASHSEED=4232 python3 -B 04-computation/jc23_weight11_u0_primitive_cm_exclusion_thm4232_independent_audit.py
```

The primary executes 12,657 exact gates.  Normal, `-O`, fixed-hash-seed, and
combined primary runs are identical.  The standard-library clean-room audit
executes 12,695 gates; its same four runs are also identical.  On Windows the
live streams use CRLF while the frozen artifacts use repository-style LF;
normalizing CRLF to LF gives byte equality.

The raw-LF SHA-256 ledger is

```text
primary.py                  23ca63d6004b8e9b3be2b455a19c336c9598dfff3c7cf1e843615b6fc806ad37
primary.out                 1fbdbefae40e75e0223fffbba0b589e8135d3448bae10d17cad6ec2987b4943d
independent_audit.py        96f13d674dd9de18f107566ab80dc7a3512c5de7ba6579ce5d1e17db17f5c1ad
independent_audit.out       7e22887a28255c3d4b6cb4ae083ee3e7f489900697201647ffbbd01f66ce7f9a
```

The computations establish the in-repo exact claims. Milne's primitive-CM
simplicity theorem, Dokchitser's regular-model theorem through THM-4045, and
proper-flat degree conservation remain the explicitly cited proof inputs.
**QED.**
