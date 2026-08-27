---
id: THM-4241
title: "W=0 hidden CM saturation and visible-hidden index-four gluing"
status: >
  PROVED RELATIVE TO THM-4230 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  On C0:x^6+y^4=1, the explicit hidden Eisenstein lattice is saturated, but
  its orthogonal sum with the saturated visible lattice has underlying
  Z-index four in the full E0 Hom lattice, with quotient O/2O witnessed by
  an explicit mixed degree-four map. The full lattice has 36,288 degree-34
  and 16,992 degree-42 vectors. Attachment collapse remains open only on
  THM-4230's finite, unenumerated marked-ratio sets S_34,S_42; W=0, M=12,
  seam entry, JC(2), and DC(2) are not closed.
source: root/cross-frontier-bridge/2026-08-26
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
related:
  - THM-4218-exact-weight-ten-hidden-elliptic-tail-degree-three-planar-jacobian-exclusion
  - THM-4222-dense-weight-eleven-primitive-cm-planar-jacobian-exclusion
  - THM-4226-dense-weight-thirteen-primitive-cm-bolza-planar-jacobian-exclusion
  - THM-4232-weight-eleven-u-zero-primitive-cm-planar-jacobian-exclusion
mistake_firewall:
  - MISTAKE-521
  - MISTAKE-522
scripts:
  - 04-computation/jc23_w0_hidden_cm_saturation_gluing_thm4241.py
  - 04-computation/jc23_w0_prime_ideals_independent_audit_thm4241.gp
  - 04-computation/jc23_w0_full_theta_independent_audit_thm4241.gp
outputs:
  - 05-knowledge/results/jc23_w0_hidden_cm_saturation_gluing_thm4241.out
  - 05-knowledge/results/jc23_w0_prime_ideals_independent_audit_thm4241.out
  - 05-knowledge/results/jc23_w0_full_theta_independent_audit_thm4241.out
script_sha256:
  - b2ead203e325b19ee0813b30d1f41c42263a1c86496506bde3a1341352459b7f
  - 7a124feb0fb99e086c0fd2a20c4c293b1de740c540c0454afa908f401e7adb08
  - 18b9366623af1c9ca402204c7f4f20f736c4c549d8c55070b7c984d42e8bfeea
output_sha256:
  - 823b6f59a3d74d52d3a4ebe0a89ab5e85972b0de639762fcc6878d506861b91b
  - 649e2b28d91935e5bf96ca9239ef9a91e872cbb90acaa76013a815f2fb3c5b62
  - 080df5a924aca4ce0d181ceaede5d0ae8bd93d73a362785aa6ed828999f6e1ad
hash_basis: raw LF bytes
audit: >
  Hostile audit conditional PASS after two repaired presentation defects:
  the scalar trace factor is 3+sqrt(3) in the stated zeta_12 embedding, and
  current THM-4230 already confines degree-34/42 attachment equality to
  finite ratio sets. Independent CM-different/overideal, pole-classifier,
  gluing, widened theta, and GP qfrep probes pass; normal/optimized outputs
  byte-match. No computation substitutes for the named geometric inputs.
---

# THM-4241 -- W=0 hidden saturation and full gluing

**PROVED RELATIVE TO THM-4230 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem deliberately separates four layers that MISTAKE-521 warns must
not be conflated:

1. the rational isotypic decomposition;
2. scalar degree-integrality of a CM ideal;
3. integrality of the actual `O`-valued Rosati pairing; and
4. visible-hidden gluing in the full integral Hom lattice.

The fixed-support issue of MISTAKE-522 is upstream in THM-4230 and is neither
reopened nor claimed to be rechecked here.

## 1. Statement

Let

```text
C0: x^6+y^4=1,                 E0: Y^2=X^3+1,
O=Z[omega],                    omega^2+omega+1=0,
iota:(x,y)->(x,-y).
```

Let `M=Hom(J(C0),E0)`, where degree is the degree of the corresponding map
`C0 -> E0`, equivalently the diagonal of the Rosati form.  Let

```text
V=M^(iota=+1),                 M^-=M^(iota=-1).
```

THM-4230 proves that `V` is saturated and has an orthogonal `O`-basis of
degree four, and constructs hidden maps `f,g=Tf` with

```text
L_exp=O f direct-sum O g,
H_L = [ 6,             -4-2omega   ]
      [ -4-2omega^2,    6          ],                (1)
det(H_L)=24,
```

where the Hermitian form is linear in its first slot.

> **Theorem.** The hidden lattice is exactly saturated:
>
> ```text
> M^-=L_exp.                                             (2)
> ```
>
> Nevertheless `D=V direct-sum L_exp` is not the full Hom lattice.  Rather
>
> ```text
> M/D ~= O/2O = F_4,                                               (3)
> |M/D|=[M:D]_Z=4.                                     (4)
> ```
>
> Thus the hidden saturation index is `1`, while the separate full
> visible-hidden gluing index is `4`.

There is an explicit degree-four mixed map producing the nonzero class in
`(3)`.  After normalizing bases, if `u,v` are the visible orthogonal
degree-four basis and

```text
P=omega^2 f+g,                h=(v+P)/2,              (5)
```

then `[u,f,g,h]` is an `O`-basis of `M`.  Its Hermitian Gram is

```text
[ 4,  0,             0,              0          ]
[ 0,  6,            -4-2omega,        2omega-2  ]
[ 0, -2+2omega,      6,               2-2omega  ]
[ 0, -4-2omega,      4+2omega,        4          ],   (6)
```

with determinant `96`.

## 2. Rational rank and the CM ideal problem

This step uses the rational decomposition only to establish rank, not
integral exhaustion.  At `W=0`, THM-4230 gives

```text
A_12 ~ E0^2 x E_1728^2
```

and the other hidden elliptic factor is `Hom`-orthogonal to `E0`.  Therefore
`M^- tensor Q` has dimension two over `Q(omega)`, hence dimension four over
`Q`.  The nonzero determinant in `(1)` shows that `f,g` are independent, so
they span this rational space.

Put

```text
K=Q(zeta_12),       omega=zeta_12^4,       T=zeta_12^5,
T^2=-omega.
```

Then `O_K=Z[zeta_12]=O direct-sum O T` and `L_exp=O_K f`.  Precomposition by
the order-twelve automorphism preserves the actual hidden Hom lattice.
Consequently `M^-` is a fractional `O_K`-ideal containing `O_K`, multiplied
by `f`.  This is the exact source and target of the ideal calculation.

In the identification `x f <-> x in K`, the scalar degree form is

```text
q(x)=Tr_{Q(sqrt(3))/Q}((3+sqrt(3)) x conjugate(x)),   (7)
sqrt(3)=zeta_12+zeta_12^(-1).
```

In the `Z`-basis `(1,omega,T,omega T)`, its coefficient matrix is

```text
Q = [ 6, -3, -3,  0 ]
    [-3,  6,  3, -3 ]
    [-3,  3,  6, -3 ]
    [ 0, -3, -3,  6 ],                                (8)
```

so `q(n)=n^t Q n`, `det(Q)=324`, and

```text
Smith(Q)=(3,3,6,6).                                   (9)
```

The associated bilinear matrix is `2Q`; `(9)` is used only for its exact
prime support.  The cyclotomic different is `p_2^2 p_3`, of norm `144`.
If `I` is any degree-integral overideal, then `I` lies in the bilinear
dual of `O_K`, so only the primes `p_2,p_3` can occur in `I/O_K`.

### Prime convention and the four scalar-integral candidates

Fix the following *principal prime ideals* (absolute norms are used):

```text
pi_2=1-zeta_12^3=1-i,          p_2=(pi_2),
N_K/Q(p_2)=4,                  pi_2^2=-2 zeta_12^3;

pi_3=zeta_12+zeta_12^3,       p_3=(pi_3),
N_K/Q(p_3)=9,                  pi_3^2=3 omega.         (10)
```

Thus `(2)=p_2^2` and `(3)=p_3^2`; each is the unique prime above the named
rational prime and has ramification index two and residue degree two.  The
words “index four” and “index nine” below mean cardinality of the quotient,
equivalently underlying `Z`-index:

```text
[p_2^-1:O_K]_Z=4,              [p_3^-1:O_K]_Z=9.      (11)
```

Complex conjugation gives `pi_2 conjugate(pi_2)=2` and
`pi_3 conjugate(pi_3)=3`.  Hence a generator of
`p_2^-e2 p_3^-e3` has degree

```text
6/(2^e2 3^e3).                                      (12)
```

Together with `(9)`, this proves that the complete list of scalar
degree-integral overideals is

```text
O_K,             p_2^-1,             p_3^-1,
(p_2 p_3)^-1,                                      (13)
```

with minimum displayed generator degrees `6,3,2,1` and underlying indices
`1,4,9,36`.  This is only the scalar degree test.

## 3. What excludes the degree-two `p_3^-1` option

For actual maps `a,b:J(C0)->E0`, principal polarizations define

```text
<a,b>_R = a composed_with b^dagger in End(E0)=O.       (14)
```

The diagonal is `deg(a)`.  Thus every pair of actual Hom vectors has an
`O`-valued Hermitian pairing, a stronger requirement than integrality of all
diagonal degrees.

On the natural `O`-bases `{pi_j^-1,T pi_j^-1}`, exact change of basis in
`(1)` gives

```text
H_(p_2^-1) = [ 3,          -2-omega ]
              [ -1+omega,   3        ],               (15)

H_(p_3^-1) = [ 2,          (-4-2omega)/3 ]
              [(-2+2omega)/3, 2       ].              (16)
```

Matrix `(15)` is `O`-valued, so polarization alone does not remove the
degree-three candidate.  Matrix `(16)` is not `O`-valued.  It contradicts
`(14)` and therefore rules out `p_3^-1`.  The product candidate in `(13)`
contains `p_3^-1`, so it is ruled out at the same time.  After the Rosati
test, the only possible *proper hidden* enlargement is the index-four
`p_2^-1`; the next section rules that out too.

There is also an independent geometric hostile control.  A degree-two map
from a genus-seven curve to an elliptic curve supplies a bielliptic
involution.  Two distinct such involutions would give, by
Castelnuovo--Severi,

```text
g(C0) <= 2*1+2*1+(2-1)(2-1)=5,
```

so a bielliptic involution would be unique and hence central in
`Aut(C0)`.  The group `G=C_6 x C_4` acts by independent scalings of `x,y`.
Its quotient `x^6:C0->P^1` has branch-inertia orders `6,4,12` at
`0,1,infinity`.  An automorphism centralizing `G` descends and must fix those
three differently labelled branch points, hence belongs to `G`.  The three
involutions of `G` are

```text
(-x,y),            (x,-y),            (-x,-y),
```

and their quotient genera are respectively `3,2,4`.  None is bielliptic.
Thus `C0` has no degree-two elliptic quotient, independently confirming the
exclusion forced by `(16)`.

## 4. Exact exclusion of the index-four degree-three candidate

Suppose that `p_2^-1 f` were contained in the actual hidden Hom lattice, and
choose a degree-three vector `h`.  Because `p_2^2=(2)`, one has `2h in
L_exp`.  Put

```text
sigma:(x,y)->(omega x,y).
```

On the hidden rational `K`-line, precomposition gives

```text
h composed_with sigma = [omega]h,
h composed_with iota  = [-1]h                              (17)
```

at the homomorphism level.  The following translation audit makes `(17)`
exact for curve maps.  First translate a representative so that `2h` equals
its `L_exp` representative exactly.  The affine error in the sigma relation
then lies in `E0[2]`.  Multiplication by `1-omega` is invertible on
`E0[2]=F_4`, so a unique translation by a two-torsion point kills that error
without changing `2h`.  The remaining affine error in the iota relation is
both two-torsion (after doubling) and fixed by `omega` (because sigma and
iota commute); it is therefore zero.  Hence `(17)` holds as written.

The fixed field of `sigma,iota` is generated by

```text
t=(1+y^2)/x^3,
x^3=2t/(t^2+1),                 y^2=(t^2-1)/(t^2+1).   (18)
```

Equivariance forces

```text
X_h=x A(t),                     Y_h=y B(t),            (19)
```

with `A,B in C(t)`.  The elliptic equation is exactly

```text
(t^2-1)B(t)^2-2t A(t)^3=t^2+1.                        (20)
```

The branch and `x`-order ledger for `t:C0->P^1` is

| base point | points above | ramification | `ord(x)` | contribution from an `A`-pole of order `m` |
|---|---:|---:|---:|---:|
| `0` or `infinity` | 2 | 3 | `+1` | `6m-2` |
| `+1` or `-1` | 3 | 2 | `0` | `6m` |
| `+i` or `-i` | 1 | 6 | `-2` | `6m+2` |
| generic | 6 | 1 | `0` | `6m` |

If `A` is regular and nonzero at `+i` or `-i`, the intrinsic pole of `x`
contributes two there; a zero of `A` at that base point cancels it completely.
Since degree three means `deg((X_h)_infinity)=6`, any generic or `+/-1`
pole already spends six and would still require two distinct zeros to cancel
the two intrinsic poles, impossible for a rational function with total pole
order one.  A pole at `+/-i` already costs at least eight.  The only surviving
patterns have one simple pole at `0` or at `infinity` and one zero
`a in {i,-i}`:

```text
A=c_0(t-a)/t                 or                 A=c_0(t-a),
a^2=-1.                                                (21)
```

For `B` to be rational, the right side of `(20)` divided by `t^2-1` must
have even valuations.  At `t=1,-1` this forces

```text
A(1)^3=-1,                       A(-1)^3=1.            (22)
```

For the first expression in `(21)`, eliminating the nonzero scalar from
`(22)` requires

```text
(1-a)^3=(-1-a)^3,
```

but the difference is `-4` modulo `a^2+1`.  For the second expression it
requires

```text
(1-a)^3+(-1-a)^3=0,
```

but the left side is `-4a` modulo `a^2+1`.  Both are contradictions.  There
is no degree-three hidden map.  This excludes `p_2^-1` and proves `(2)`.

## 5. Every full class lies in the half-sum, and at most one glue line exists

Set `L=L_exp` now that hidden saturation is proved, and `D=V direct-sum L`.
For any `m in M`, define

```text
v=m+iota m in V,                  ell=m-iota m in L.   (23)
```

Then

```text
2m=v+ell,                  D subset M subset (1/2)D.  (24)
```

This proves the often implicit “every class lies in `(V+L)/2`” assertion;
no rational direct-sum claim is being substituted for it.  Since `M` is an
`O`-module, `(24)` makes `M/D` a vector space over

```text
O/2O=F_4.
```

In particular, if its dimension is `d`, then its module order, its quotient
cardinality, and its underlying `Z`-index are all

```text
|M/D|=[M:D]_Z=4^d.                                   (25)
```

It is incorrect to call a one-dimensional quotient here “index two.”

Projection to the negative part gives an injection

```text
M/D -> L/2L,                  [m] |-> ell mod 2L.      (26)
```

Indeed, if `ell=2ell_0`, then `m-ell_0` is fixed by `iota`, hence lies in
`V`, so `m in D`.

The two eigenspaces are Rosati-orthogonal.  Therefore

```text
q(m)=(q(v)+q(ell))/4.                                 (27)
```

Visible saturation gives `q(V) subset 4Z`; integrality of `q(m)` forces
`q(ell)=0 mod 4`.  Because every entry of `(1)` lies in `2O`, changing
`ell` by an element of `2L` changes its degree by a multiple of four, so this
condition is well-defined on `L/2L`.  An exact four-by-four residue check in
`L/2L=F_4^2` gives

```text
{(a,b):q(af+bg)=0 mod 4}
  ={(0,0),(1,omega),(omega,1+omega),(1+omega,1)}
  ={(omega^2 b,b):b in F_4}.                          (28)
```

It is a unique `F_4`-line, so `(26)--(28)` prove `d<=1` and
`[M:D]_Z<=4`.

## 6. Explicit degree-four map and exact nontrivial gluing

Choose complex constants satisfying

```text
r^2-12r+24=0,                lambda^4(r-9)=1,
c=r lambda,                  alpha^3=c,              epsilon^2=-1,
Q_lambda(y)=y^2-(c/2)y+3lambda^2.                     (29)
```

Define

```text
H_lambda:C0->E0,
X=alpha (y-lambda)/x^2,
Y=epsilon Q_lambda(y)/x^3.                            (30)
```

The load-bearing identity is

```text
c(y-lambda)^3+1-y^4=-Q_lambda(y)^2.                  (31)
```

Direct expansion reduces the `y^2` coefficient to
`lambda^2(r^2-12r+24)/4` and the constant coefficient to
`1-lambda^4(r-9)`; all other coefficients cancel.  Thus `(31)` proves
`Y^2=X^3+1` exactly.

The equation `lambda^4=1` would force `r=10`, which does not satisfy the
quadratic in `(29)`.  Hence at each of the four points `x=0`, the numerator
of `X` is nonzero and `X` has exact pole order two.  These are all its poles:
at either point at infinity, `ord(x)=-2` and `ord(y)=-3`, so `X` has a zero,
not a pole.  Therefore

```text
deg((X)_infinity)=4*2=8,          deg(H_lambda)=4.     (32)
```

Both iota eigenspaces occur.  With `eta=H_lambda^*(dX/Y)`, differentiation
on `C0` gives

```text
eta=(alpha/(3epsilon)) x^-5 F(y)dy,
F(y)=(3+y^4-4lambda y^3)/Q_lambda(y),                 (33)
F(0)=lambda^-2,                 F'(0)=c/(6lambda^4).  (34)
```

Both values in `(34)` are nonzero.  Since
`iota^*eta=-constant*x^-5 F(-y)dy`, `F(0)` gives a nonzero anti-invariant
part and `F'(0)` a nonzero invariant part.  Thus the two projections of
`H_lambda` to `V` and `L` are both nonzero.

If `H_lambda` belonged to `D=V direct-sum L`, orthogonality and the minimum
degrees `4` and `6` would give degree at least `10`, contradicting `(32)`.
Therefore `M/D` is nonzero.  Combined with `d<=1`, this proves `(3)--(4)`.

Let

```text
v_H=H_lambda+H_lambda composed_with iota,
ell_H=H_lambda-H_lambda composed_with iota.
```

The parallelogram law gives `q(v_H)+q(ell_H)=16`.  Both are nonzero,
visible degrees are multiples of four, and hidden degrees are multiples of
six, so necessarily

```text
q(v_H)=4,                         q(ell_H)=12.          (35)
```

Scaling by an `O`-unit and changing representatives by `D` normalizes the
unique residue line `(28)` to `(5)`.  This proves the basis and Gram in
`(5)--(6)`; its determinant is `96` (the direct-sum determinant was
`16*24=384`, reduced by the norm/index factor four).

## 7. Exact degree-34/42 consequence and finite-ratio boundary

Exact enumeration of `(6)` gives the theta coefficients through degree 42:

```text
q : count
0:1, 4:60, 6:72, 8:324, 10:864, 12:276, 14:2592,
16:2868, 18:1152, 20:5832, 22:9504, 24:2556,
26:15552, 28:15456, 30:6480, 32:22356,
34:36288, 36:7836, 38:49248, 40:44280, 42:16992.     (36)
```

In the normalized basis `(5)`, explicit witnesses are

```text
q(h-f+g)=34,
q(-f-omega g+(1-omega)h)=42.                         (37)
```

These counts are raw Hom-lattice vectors, before quotienting by target or
source symmetries.  They have two precise consequences:

1. THM-4230's absence of degree 34 and its finite-exact noncollapse of all
   192 degree-42 vectors **inside `L_exp`** remain correct.
2. Those sublattice results no longer supply an arithmetic obstruction in
   the full lattice: `M` has 36,288 degree-34 and 16,992 degree-42 vectors.

Use THM-4230, section 6.2, for the current attachment boundary.  For
`d in {34,42}`, let `S_d` be the set of marked ratios `U/Z` for which
some admissible attachment orbit on `C0` and some `m in M` of degree `d`
have all twelve attachment images equal.  THM-4230 proves

```text
S_34 and S_42 are finite.                              (38)
```

It does not enumerate them or prove them empty.  The explicit map
`H_lambda` cannot itself collapse twelve distinct attachments, because a
fibre of a degree-four map has total degree four; this says nothing about the
vectors in `(36)--(37)`.

Thus the exact remaining `W=0` task is no longer an unrestricted attachment
locus.  It is the enumeration or emptiness proof for the finite sets
`S_34,S_42`, now using the explicit full lattice `(5)--(6)` and its unique
visible-hidden coset.  A decisive computation should reduce first by scalar
`mu_6`, the order-twelve source action, and every proven
attachment-preserving symmetry, then perform exact attachment-resultant
tests.  No such emptiness claim is made here.

## 8. Preserved and lost data; nonclaims

The maps and operations above preserve exactly the following data:

- the `iota` eigenspace splitting over `Q`;
- the actual integral Rosati pairing, not merely its diagonal degree;
- the labelled CM action needed to make the hidden lattice an `O_K`-ideal;
- the full `O/2O` residue class controlling visible-hidden gluing; and
- exact map degree via the `X`-pole divisor.

They deliberately do not preserve or decide:

- individual attachment images under an arbitrary full-lattice vector;
- emptiness or nonemptiness of THM-4230's finite ratio sets `S_34,S_42`;
- the complete hidden-Hom locus away from `W=0`;
- any wall outside the THM-4230 gate; or
- `M=12`, seam entry, `JC(2)`, or `DC(2)`.

No claim of full `Hom` is made from rational decomposition alone: `(2)` uses
the CM overideal classification, Rosati integrality, and the independent
degree-three geometric obstruction; `(3)` then performs a separate integral
gluing audit.

## 9. Reproduction

From the repository root, run

```bash
python3 -B 04-computation/jc23_w0_hidden_cm_saturation_gluing_thm4241.py \
  > /tmp/jc23-w0-thm4241.out
python3 -B -O 04-computation/jc23_w0_hidden_cm_saturation_gluing_thm4241.py \
  > /tmp/jc23-w0-thm4241-opt.out
cmp /tmp/jc23-w0-thm4241.out /tmp/jc23-w0-thm4241-opt.out
cmp /tmp/jc23-w0-thm4241.out \
  05-knowledge/results/jc23_w0_hidden_cm_saturation_gluing_thm4241.out

gp -fq 04-computation/jc23_w0_prime_ideals_independent_audit_thm4241.gp \
  > /tmp/jc23-w0-primes-thm4241.out
cmp /tmp/jc23-w0-primes-thm4241.out \
  05-knowledge/results/jc23_w0_prime_ideals_independent_audit_thm4241.out

gp -fq 04-computation/jc23_w0_full_theta_independent_audit_thm4241.gp \
  > /tmp/jc23-w0-theta-thm4241.out
cmp /tmp/jc23-w0-theta-thm4241.out \
  05-knowledge/results/jc23_w0_full_theta_independent_audit_thm4241.out
```

The Python certificate uses exact `Fraction`/SymPy arithmetic, reconstructs
the corrected `3+sqrt(3)` trace Gram, and exhausts a bounded theta census
whose retained vectors do not touch the coordinate boundary.  The first
PARI/GP path independently reports discriminant `144`, class number `1`,
the unique primes above `2,3`, their norms `4,9`, the prime-square unit
relations, and the corrected trace sign.  The second converts the full
Hermitian lattice to an eight-dimensional integral quadratic matrix of
determinant `2916` and independently reproduces `(36)`; PARI's `qfrep`
counts each pair `{x,-x}` once, so its nonzero counts are doubled.

The computation certifies the displayed algebra and finite enumerations.  It
does not replace the standard facts that Rosati composition is `O`-valued,
that two distinct bielliptic involutions obey Castelnuovo--Severi, or the
rational CM decomposition inherited from THM-4230.
