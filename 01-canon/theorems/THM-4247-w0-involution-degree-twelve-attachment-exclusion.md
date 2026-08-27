---
id: THM-4247
title: "W=0 involution-projection refinement and degree-twelve attachment exclusion"
status: >
  PROVED RELATIVE TO THM-4230/4241 + VERIFIED-EXACT + HOSTILE-AUDITED. In
  the full W=0 Hom lattice, the visible/hidden involution projections refine
  the 36,288 degree-34 and 16,992 degree-42 vectors into exact bivariate theta
  rows. Every hidden projection has positive degree divisible by twelve. The
  degree-twelve hidden shell has 24 maps in two symmetry orbits, and neither
  orbit can collapse the twelve admissible attachments. Together with the
  inherited fibre-degree and pure-hidden exclusions, this removes 2,112 raw
  degree-34 vectors and 864 raw degree-42 vectors, leaving 34,176 and 16,128.
  The marked-ratio sets remain finite but unenumerated; W=0, the M=12 seam,
  JC(2), and DC(2) remain OPEN.
source: root/sun-odd-cycles/2026-08-26
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4241-w0-hidden-cm-saturation-and-visible-hidden-index-four-gluing
related:
  - THM-4218-exact-weight-ten-hidden-elliptic-tail-degree-three-planar-jacobian-exclusion
  - MISTAKE-521
  - MISTAKE-522
scripts:
  - 04-computation/jc23_w0_involution_projection_histogram_thm4247.py
  - 04-computation/jc23_w0_hidden_degree12_attachment_audit_thm4247.py
  - 04-computation/jc23_w0_hidden_degree12_orbit0_symbolic_audit_thm4247.py
outputs:
  - 05-knowledge/results/jc23_w0_involution_projection_histogram_thm4247.out
  - 05-knowledge/results/jc23_w0_hidden_degree12_attachment_audit_thm4247.out
  - 05-knowledge/results/jc23_w0_hidden_degree12_orbit0_symbolic_audit_thm4247.out
script_sha256:
  - fcfb1334edbc4f1ea897a2f08caca27256631317d29ab60503d14c3b53b60c30
  - 9dc1f8614db388f463acb93951285fb4b8245cf20c976d1d2129cb883ccf9c28
  - 19497a9d81875b64e3edf6f8ca2bfc492acdbfd1aed4582a18fbf65b9a50f41c
output_sha256:
  - da4bb1e48c0f6cedb4e1d26a4ab63a5074cb6101976d7b97544b622ef6150251
  - 219eaf90c319dcd5f0f9ab98cbb7e8b9b1db75b96b4df5c90e0d9b593015fa7c
  - bff39073b8e387a80dac850d231dee68c7f0e04531e4cf3153dde733fab841db
hash_basis: raw LF bytes
---

# THM-4247 -- W=0 involution refinement and degree-twelve exclusion

**PROVED RELATIVE TO THM-4230/4241 + VERIFIED-EXACT + HOSTILE-AUDITED.**

## 1. Statement

Use the notation of THM-4230/4241:

```text
C0: x^6+y^4=1,              E0:Y^2=X^3+1,
M=Hom(J(C0),E0),            iota:(x,y)->(x,-y),
V=M^(iota=+1),              L=M^(iota=-1).
```

For `m in M`, put

```text
v=m+m composed_with iota,
ell=m-m composed_with iota.                           (1)
```

Then `v in V`, `ell in L`, and the Rosati degree form satisfies

```text
q(v)+q(ell)=4q(m).                                    (2)
```

For every raw full-Hom vector in the two response degrees, the exact hidden
projection histograms are

```text
q(m)=34:
q(ell):count =
12:1536, 24:2304, 36:5952, 60:5760, 72:8064,
84:5376, 108:4992, 120:1728, 132:576;                 (3)

q(m)=42:
q(ell):count =
12:672, 24:288, 36:2304, 60:288, 72:3456,
84:3072, 108:3744, 120:1728, 132:576,
156:672, 168:192.                                     (4)
```

Every displayed hidden degree is a positive multiple of `12`.

The hidden degree-twelve shell consists of exactly `24` maps in two free
orbits of size `12` under postcomposition by `mu_6` and precomposition by
`T`. No map in this shell sends all twelve admissible `W=0` gate-interior
attachments to the origin of `E0`.

If a degree-34 or degree-42 full map collapses the twelve attachments, its
hidden projection must send all of them to the origin. Consequently the
following disjoint rows are impossible:

```text
degree 34: q(ell)=12, count 1536;
           (q(v),q(ell))=(4,132), count 576;
degree 42: q(ell)=12, count 672;
           (q(v),q(ell))=(0,168), count 192.           (5)
```

The last two exclusions use, respectively, the degree-four fibre bound and
THM-4230's pure-hidden degree-42 calculation. The new degree-twelve theorem
supplies the first and third rows. The remaining exact raw-vector ledgers are

```text
degree 34:
{24:2304,36:5952,60:5760,72:8064,84:5376,
 108:4992,120:1728},                         total 34176;

degree 42:
{24:288,36:2304,60:288,72:3456,84:3072,
 108:3744,120:1728,132:576,156:672},         total 16128. (6)
```

These are lattice-vector counts before quotienting by source or target
automorphisms. They do not enumerate the finite marked-ratio sets
`S_34,S_42` from THM-4230.

## 2. Integral projections and the exact theta refinement

THM-4241 gives an `O=Z[omega]` basis `[u,f,g,h]` of `M` with

```text
2h=v0+P,                    P=omega^2 f+g,
q(u)=q(v0)=4,
H_L=[[6,-4-2omega],[-4-2omega^2,6]].                  (7)
```

Write

```text
m=a u+b f+c g+d h,              a,b,c,d in O.
```

Since `iota` fixes `u,v0` and negates `f,g`,

```text
ell=(2b+omega^2 d)f+(2c+d)g,
v=2a u+d v0.                                           (8)
```

The eigenspaces are Rosati-orthogonal, so the parallelogram identity proves
`(2)`. Formula `(8)` gives

```text
q(v)=16N(a)+4N(d).                                    (9)
```

Conversely, the hidden coefficients in `(8)` lie in the exact residue class
`(omega^2d,d) mod 2L`, and this class recovers `b,c` uniquely. Enumerating the
bounded positive-definite shells in `(7)`, grouped by that residue class and
then by `(9)`, is therefore a bijective enumeration of `M`, not a superset.
The companion reproduces THM-4241's totals `36,288` and `16,992`; none of
its coordinate boxes touches a boundary.

The visible degree is divisible by four, while the hidden Gram makes every
hidden degree divisible by six. Equation `(2)` makes `q(ell)` divisible by
four as well, hence by twelve. If `ell=0`, then `m` is visible and its degree
is divisible by four, impossible at degrees 34 and 42. This proves positivity
and divisibility independently of the finite census.

## 3. Collapse and the involution sidecar

Choose the Abel-Jacobi base point `Q_*` among the `iota`-fixed points `y=0`.
The twelve attachments satisfy

```text
iota Q_j=Q_(j+6).                                     (10)
```

If their common image under `m` is `P`, then

```text
ell(Q_j)=m(Q_j)-m(iota Q_j)=O,
v(Q_j)=m(Q_j)+m(iota Q_j)=2P.                         (11)
```

This is only a necessary condition. The converse can lose a two-torsion
difference, exactly as the index-four gluing sidecar in THM-4241 predicts.

## 4. Uniform denominator bound in hidden degree twelve

THM-4230 proves the coefficient-form mechanism; the degree-twelve
specialization is recorded here explicitly. For every nonzero `r in L`, the
pair `(X_r/x,Y_r/y)` lies in `t^-1 K[t^2]`, and the elliptic group law
preserves oddness. Both points over `t=0` map to `O` under `f,g`, hence under
`r`. Therefore the reduced denominator `d_r(t)` of `X_r/x` is odd and has
odd degree.

A finite denominator root of multiplicity `e` contributes at least `6e` to
the `X`-pole divisor, except that `t=0` contributes `6e-2`; the `t=+-i`
fibres cost more. Thus

```text
6 deg(d_r)-2 <= 2 deg(r).                             (12)
```

At `deg(r)=12`, this gives `deg(d_r)<=4`; oddness sharpens it to

```text
deg(d_r)<=3.                                          (13)
```

## 5. The complete degree-twelve shell

Exact enumeration in the hidden Gram `(7)` gives `24` degree-twelve vectors.
Scalar `mu_6` and

```text
T(a,b)=(-omega b,a)                                   (14)
```

split them into two disjoint free orbits of size `12`. Representatives may
be taken as

```text
r0=(1-omega)f+g,                 r1=f+omega g.         (15)
```

Post-units fix the origin and do not alter poles. Precomposition by `T` sends
`t` to `-1/t` and preserves the wall pair `{+1,-1}`, so `(15)` is a complete
set of cases.

For an equivalent representative of the `r1` orbit, exact coefficient-form
elliptic addition over the characteristic-zero relations

```text
z^4-z^2+1=0,
p^2-(1+2z-z^3)p+1=0
```

gives

```text
d_r1(t)=C1 t(t^2-1),                                  (16)
```

where the absolute resultant norm of `C1` is `65,536`. Exact numerator
resultants at `t=0,+1,-1` are nonzero, so none of these factors cancels.

For `r0`, direct addition first forms `(1-omega)f` and then adds `g`. After
the same algebraic reductions the raw denominator is

```text
C0 t(t^2-1) R8(t^2),                                  (17)
```

with nonzero leading-coefficient resultant `18,339,659,776`. The raw
numerator resultant norms at `t=0,+1,-1` are respectively

```text
4,477,456,
40,290,721,869,103,654,477,234,176,
40,290,721,869,103,654,477,234,176.                  (18)
```

Thus `t(t^2-1)` survives reduction. Bound `(13)` forces every remaining
factor in `(17)` to cancel, proving

```text
d_r0(t)=C0' t(t^2-1).                                 (19)
```

The orbit-zero certificate explicitly verifies that setting the common
target scale to one is legal: elliptic addition is homogeneous with weights
`(2,3)` in `(X,Y)`.

As hostile controls, the four exact good reductions

```text
(q,z,p,s)=(313,29,135,21),(349,24,246,28),
          (373,69,297,33),(397,157,161,27)             (20)
```

give reduced denominator degree three and monic reciprocal gcd `t^2-1` for
both orbit representatives.

THM-4230's attachment criterion says that a hidden map sends all twelve
attachments to `O` only if the denominator and its `T`-reciprocal have a
common root. Equations `(16),(19)` show that the only common finite roots are
`t=+-1`. But the node coordinate is

```text
Z/U=((t^2-1)/(2t))^2,                                 (21)
```

so both roots lie on `Z=0`, excluded by the gate. This proves the complete
degree-twelve attachment exclusion.

## 6. Exact row deletion and scope

Equation `(11)` deletes both `q(ell)=12` rows. In degree 34, the remaining
row `(q(v),q(ell))=(4,132)` cannot collapse: a nonconstant degree-four map
has fibre divisor of total degree four and cannot contain twelve distinct
gate-interior attachments. In degree 42, `(0,168)` is exactly the 192-vector
pure-hidden shell excluded by THM-4230. The rows are disjoint, and subtraction
from `(3),(4)` proves `(5),(6)`.

The degeneration dual graphs considered along this route do not supply a
shortcut: their exact repeated path lengths are even and their other cycle
blocks are two-edge circuits, so the resolved graphs are odd-cycle-free.
The lawful finite-order information is instead the attachment action:
`tau^8` has four three-cycles and `iota=tau^6` supplies the incompatible
three- and two-torsion conditions used upstream. No tournament is asserted.

The theorem does not quotient the remaining vectors, enumerate or empty
`S_34,S_42`, handle the other walls, close `W=0`, prove seam entry, or prove
JC(2) or DC(2).

## 7. Reproduction

```text
python -B 04-computation/jc23_w0_involution_projection_histogram_thm4247.py
python -B -O 04-computation/jc23_w0_involution_projection_histogram_thm4247.py
python -B 04-computation/jc23_w0_hidden_degree12_attachment_audit_thm4247.py
python -B -O 04-computation/jc23_w0_hidden_degree12_attachment_audit_thm4247.py
python -B 04-computation/jc23_w0_hidden_degree12_orbit0_symbolic_audit_thm4247.py
python -B -O 04-computation/jc23_w0_hidden_degree12_orbit0_symbolic_audit_thm4247.py
```

All normal/optimized output pairs are byte-identical. Frontmatter hashes use
raw LF bytes. **QED.**
