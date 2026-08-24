---
id: THM-3932
title: "Infinity-component conic torus sextics: fold-one collapse, fold-two exclusion, and a fold-three family"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. If
  an irreducible torus sextic has affine normalization A1 and its projective
  quadratic coefficient is a reduced conic containing the line at infinity,
  then its Cardano fold degree is 1, 2, or 3. Fold one is necessarily a
  triangular polynomial-coordinate pullback of the cusp, fold two is empty,
  and fold three has a complete trace/norm parameter grammar and is nonempty.
  The explicit sextic 4X^3Z^3-(Y^3-X^2Z)^2 has normalization A1 with one
  infinity address, a nonconstant coefficient Jacobian, and a genuine
  order-three divisor class on its normal quadratic resolvent. Its associated
  normal cubic order is nevertheless globally monogenic, so it is a sharp
  branch target rather than a planar Keller completion.
source: jc_degree6_one_place/infinity_component_conic / post-THM-3928 remaining singular-conic lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23). The audit independently
  reconstructed the Cardano fibre-product irreducibility argument, all three
  fold rows, and the fold-three trace/norm coefficient equations, including
  the total-degree bounds after the removable triangular shift. It also
  checked the projective normalization, the contact-two finite branches and
  (2,9) infinity cusp, both normality arguments, and the distinction between
  a genuine order-three subgroup and an uncomputed full resolvent class
  group. At the singular local ring, the reflexive prime (x,q+w) needs two
  generators, so div(q+w)=3D+ has exact rather than merely formal order
  three. In normal and optimized mode the assertion-free 55-gate companion
  LF-normalizes exactly to the frozen raw-LF output; raw and semantic hashes
  and documentation checks pass.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
related:
  - THM-3879-rational-torus-sextic-c3-packet-one-place-tradeoff
  - THM-3925-fivefold-conic-contact-torus-sextic-one-place-fold
  - THM-3928-split-affine-conic-one-place-fold-degree-barrier
  - THM-3929-root-regular-one-place-linear-color-cubic-is-monogenic
script: 04-computation/jc2_infinity_component_linear_conic_torus_sextic_thm3932.py
output: 05-knowledge/results/jc2_infinity_component_linear_conic_torus_sextic_thm3932.out
script_sha256: 39f362c41cdbf6fd481af85dc91b96ed5f1151810949895ef248b32913a39374
output_sha256: 5fd6493bfc51a36626f2f45be5ba0bfcd209fd415941988ee1b8812d23e27322
semantic_sha256: f21747d66cce559781355bb75789500399f7d8eb8281254e33668ffb0697a849
hash_basis: raw LF bytes
---

# THM-3932 -- the infinity component leaves a genuine cubic fold

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero. Suppose
the projective quadratic coefficient of a torus sextic is a reduced conic
having the chosen infinity line as one component. After a projective change
and a nonzero rescaling, write

```text
Q2=3XZ,             p=Q2(X,Y,1)=3x,
Q3 in k[X,Y,Z]_3,   q=Q3(x,y,1),
Delta=4Q2^3-27Q3^2.                                      (1)
```

Assume that `Delta` is irreducible of degree six and that its affine
normalization is `A1_s`; equivalently, its projective normalization is `P1`
with exactly one point over `Z=0`. Then the Cardano parameter has degree

```text
d in {1,2,3}.                                             (2)
```

The three rows have an exact outcome:

1. `d=1` forces `(p,q)` to be a triangular polynomial automorphism, hence
   the branch is a polynomial-coordinate cusp fold;
2. `d=2` is impossible for a sextic;
3. `d=3` has the complete trace/norm grammar in Section 4 and is nonempty.

In particular, the last singular-conic lane left by THM-3928 is not empty.
A one-parameter family is

```text
Q2=3XZ,
Q3=Y^3+rXYZ-X^2Z,
H=s^3+rs-2,
(x,y)=(H^2,sH).                                           (3)
```

For every `r`, `(3)` is birational onto an irreducible one-place sextic and
has fold degree three. The member `r=0` has equation

```text
F=Delta/27=4X^3Z^3-(Y^3-X^2Z)^2.                         (4)
```

It is not a polynomial-coordinate fold: the affine coefficient map has

```text
Jac(3x,y^3-x^2)=9y^2.                                    (5)
```

Its quadratic resolvent carries a genuine nonzero class of exact order
three. However, the natural depressed cubic has leading coefficient one and
is already its normal integral closure. Thus `(4)` advances the branch and
Kummer geometry but does not evade the monogenic Keller obstruction.

## 1. One place makes the Cardano parameter polynomial

Let `nu:A1_s -> V(Delta)` be the normalization. In its function field put

```text
h=3(q o nu)/(2(p o nu)).                                  (6)
```

The torus equation gives

```text
x o nu=h^2,                    q o nu=2h^3.               (7)
```

The first identity makes `h` integral over `k[s]`. Since `k[s]` is normal,

```text
h in k[s].                                                 (8)
```

It is nonconstant: otherwise both coefficient coordinates would be constant
on the sextic. Write `d=deg_s h`. Because `h` is recovered rationally from
`(x,y)` by `(6)`,

```text
k(s)=k(h,y),                    [k(s):k(h)]=d.             (9)
```

Now form the Cardano fibre product

```text
R(T,Y)=q(T^2,Y)-2T^3.                                     (10)
```

It is irreducible. Indeed

```text
R(T,Y)R(-T,Y)=q(T^2,Y)^2-4T^6=-F(T^2,Y),                 (11)
```

where `F=Delta/27`. If `R=AB` nontrivially, then
`A(T,Y)A(-T,Y)` and `B(T,Y)B(-T,Y)` are nonconstant polynomials in
`T^2,Y`; `(11)` would factor the irreducible polynomial `F`. Therefore `R`
is the minimal polynomial of `y` over `k(h)`. If `m=deg_Y q`, then

```text
d=deg_Y R=m<=3,                                           (12)
```

which proves `(2)`. This is the basic reversal from THM-3928: a product of
two affine line coordinates forced a large fold, while the single affine
line left by an infinity component caps the fold at three.

## 2. Fold one is exactly the coordinate fold

For `d=1`, equation `(12)` writes

```text
q(x,y)=a(x)y+b(x),             deg a<=2, deg b<=3.         (13)
```

Since `h` is an affine coordinate on the normalization, `(7)` and `(13)`
give

```text
y=(2h^3-b(h^2))/a(h^2).                                  (14)
```

If `deg a>=1`, polynomiality of `y` makes the denominator divide the
numerator and gives

```text
deg_s y <= 6-2 deg a <=4.                                 (15)
```

But a basepoint-free polynomial parametrization with one infinity place has
degree

```text
deg Gamma=max(deg_s x,deg_s y).                           (16)
```

Here `deg_s x=2`, so `(15)` contradicts `deg Gamma=6`. Hence

```text
a in k*.                                                   (17)
```

The map

```text
(x,y) |-> (3x,a y+b(x))                                   (18)
```

is triangular with a polynomial inverse. Thus fold one is precisely the
polynomial-coordinate cusp pullback. Its depressed cubic is globally
monogenic, exactly as in THM-3925.

## 3. Fold two is empty by trace parity

For `d=2`, write

```text
q(x,y)=a_2(x)y^2+a_1(x)y+a_0(x),           deg a_i<=3-i.  (19)
```

The polynomial `y(s)` is integral over `k[h]`. Its monic minimal polynomial
therefore belongs to `k[h][Y]`. Since the irreducible polynomial `(10)` is
an associate of that monic polynomial, its leading coefficient
`a_2(h^2)` must be a scalar; otherwise it would be a nonconstant content
factor of `R`. Thus `a_2 in k*`.

After translating and scaling `s`, write

```text
h=alpha s^2+beta,                         alpha!=0.        (20)
```

The trace of `y(s)` from `k(s)` to `k(h)` is its value at `s` plus its value
at `-s`. On the one hand, `(19)` says

```text
Tr(y)=-a_1(h^2)/a_2,                                      (21)
```

a polynomial having only even powers of `h`. On the other hand, `(16)` and
`deg_s x=4` force `deg_s y=6`. Its nonzero `s^6` term contributes a nonzero
multiple of `(h-beta)^3` to the trace. Hence `Tr(y)` has degree three in
`h`, contradicting `(21)`. No fold-two sextic exists.

## 4. The complete fold-three grammar

For `d=3`, the coefficient of `y^3` in `q` is already a nonzero scalar.
Scale `y` to make it one, and make an affine change of the normalization
coordinate so that

```text
h=s^3+r s+beta.                                           (22)
```

Every polynomial `y(s)` of degree at most six has a unique expression

```text
y=E(h)+v,
E=e_0+e_2h^2,
v=(B_0+B_1h)s+C_0(s^2+2r/3).                             (23)
```

Here the absence of an `h` term in `E` is exactly the requirement that the
trace coefficient come from `k[h^2]`. Starting first with
`C=C_0+C_1h`, the `h^5` coefficient of `Norm(v)` is `C_1^3`; the required
constant coefficient of `(10)` has no `h^5` term, so `C_1=0`. A direct
trace/norm calculation then says that `(10)` has all coefficients in
`k[h^2]`, with odd norm part exactly `2h^3`, if and only if

```text
2B_0B_1r-3B_0C_0+3B_1C_0 beta=0,                         (24)

B_1^2(3B_0-B_1 beta)=2,                                  (25)

B_0^3-3 beta B_0^2B_1+(4/3)B_0B_1C_0r^2-B_0C_0^2r
 +beta B_1C_0^2r-2 beta C_0^3=0.                         (26)
```

These equations are necessary and sufficient. To see sufficiency, take the
monic minimal polynomial

```text
M(Y)=Norm_{k(s)/k(h)}(Y-y).                               (27)
```

Equations `(24)-(26)` make `M(Y)+2h^3` a polynomial in `h^2,Y` with the
total-degree bounds of a cubic `q(x,Y)`. Put `x=h^2` to define that `q`.
Equation `(25)` makes `B_1` nonzero, so `v` is not in `k(h)` and `(27)` is
irreducible of degree three. Conversely, coefficient comparison in any
fold-three branch gives exactly `(23)-(26)`. The removable term `E` is a
triangular change `y |-> y-e_0-e_2x`.

The family `(3)` is the especially simple solution

```text
beta=-2,        B_0=C_0=e_0=e_2=0,        B_1=1.          (28)
```

Indeed

```text
h=s^3+rs-2,             x=h^2,             y=sh,
q=y^3+rxy-x^2=2h^3.                                     (29)
```

The inverse on a dense open is explicit:

```text
h=q/(2x),                         s=y/h=2xy/q.             (30)
```

Thus `(29)` is birational. Its polynomial parametrization has degrees
`(6,4)`, so `(16)` gives image degree six. Since its degree-six image lies
in the degree-six polynomial `4x^3-q^2`, that polynomial is the irreducible
equation of the image. This proves all claims about the one-parameter family.

## 5. The explicit projective normalization and address packet

Set `r=0` and homogenize `(29)`. With

```text
H=T^3-2S^3,
nu([S:T])=[H^2:S^2TH:S^6],                               (31)
```

all three coordinates have degree six and substitution gives `(4)`. They
have no common zero: if `S!=0`, the last coordinate is nonzero; if `S=0`,
then `T!=0` and the first coordinate is `T^6`. Thus `(31)` is a
basepoint-free projective normalization. Its only zero of the last
coordinate is `S=0`, which maps to

```text
P_infinity=[1:0:0].                                       (32)
```

Consistently, `F(X,Y,0)=-Y^6`. The branch therefore has one projective
infinity point and exactly one normalization address over it.

The complete singularity ledger has two supports. Affinely, the singular
ideal has reduced support only at `(0,0)`. Its normalization fibre is

```text
s^3=2,                                                     (33)
```

so it contains three distinct points. Near a root `alpha`, `h` is a local
parameter and

```text
x=h^2,                         y=alpha h+O(h^2).           (34)
```

The three smooth branches have distinct quadratic coefficients and each
pair has intersection multiplicity two. Hence

```text
delta_(0,0)=3*2=6.                                        (35)
```

At infinity use the chart `X=1`, with `u=Y/X`, `v=Z/X`, and local parameter
`z=S/T`. Then

```text
u=z^2/(1-2z^3),
v=z^6/(1-2z^3)^2,
v-u^3=-2z^9/(1-2z^3)^3.                                  (36)
```

This is a unibranch `(2,9)` cusp, so

```text
delta_infinity=(2-1)(9-1)/2=4.                            (37)
```

The arithmetic genus of a sextic is ten, and `(35)+(37)=10`; the explicit
`P1` normalization and the singularity packet agree exactly. The affine
normalization is `P1\{infinity}=A1`.

## 6. The natural cubic is normal and monogenic

For `(4)`, the affine coefficient map and its Jacobian are

```text
(p,q)=(3x,y^3-x^2),                    Jac(p,q)=9y^2.      (38)
```

Therefore this fold cannot become the triangular coordinate fold through
pre- and post-composition by polynomial automorphisms: those operations
change a Jacobian only by a nonzero scalar.

The associated binary cubic is

```text
Phi(U,V)=U^3-3xUV^2-(y^3-x^2)V^3.                         (39)
```

It represents `1` at `(1,0)`, and its cubic algebra is

```text
A=k[x,y,u]/(u^3-3xu-(y^3-x^2)).                           (40)
```

The discriminant of `(39)` is `27F`. The generic cubic is irreducible: a
rational root would be a polynomial root by Gauss' lemma; degree comparison
makes it affine-linear, and direct coefficient comparison gives none. The
singular locus of `(40)` is the single closed point `(x,y,u)=(0,0,0)`.
As a hypersurface it is `S2`, and an isolated singularity has codimension
two, so `(40)` is normal.

Consequently `(40)` is already the integral cubic completion in its
fraction field and is globally monogenic. THM-3801 therefore excludes it as
the finite cubic completion of a constant-unit planar Keller map. The
nonconstant Jacobian in `(38)` and monogenicity in `(40)` are distinct
obstructions; neither invalidates the branch geometry proved above.

## 7. The quadratic resolvent has a genuine three-class

Put

```text
B=k[x,y,w]/(w^2-((y^3-x^2)^2-4x^3)),
q=y^3-x^2,
a=q+w,                         b=q-w.                     (41)
```

The resolvent polynomial is irreducible. Its singular ideal has support only
at `(x,y,w)=(0,0,0)`, so the same `S2+R1` argument proves that `B` is normal.
Over `x=0` it has the two height-one primes

```text
D^+=(x,a),                         D^-=(x,b).              (42)
```

Since

```text
ab=4x^3,                                                   (43)
```

and the opposite factor is a unit at each generic point,

```text
div(x)=D^++D^-,                 div(a)=3D^+.               (44)
```

The second relation does not merely produce a formal candidate. The class
`[D^+]` is nonzero. At the singular maximal ideal, eliminate `w` in favour
of `a`; the local equation is

```text
a^2-2qa+4x^3=0.                                           (45)
```

It has no linear term. Hence the prime ideal `(x,a)` has two independent
generators modulo `m(x,a)` and is not locally principal. Moreover
`B/(x,a)=k[y]` is Cohen--Macaulay, so the depth lemma makes `(x,a)` reflexive.
It therefore represents a nonzero Weil class. Equation `(44)` proves

```text
ord([D^+])=3,                         [D^-]=-[D^+].         (46)
```

This is an intrinsic `C3` direction predicted when the affine part of
`Q2` has one component. It is exactly the Cardano class: adjoining a cube
root of `(q+w)/2`, and putting the conjugate cube root so that their product
is `x`, gives the splitting field of `(40)`. Thus this genuine class recovers
the already normal monogenic cubic; it does not manufacture a second,
nonmonogenic completion.

No claim that `(46)` is the whole class group is made. In fact, with local
weights

```text
wt(y)=1,                         wt(x)=2,       wt(w)=3,
```

the initial equation at the resolvent singularity is

```text
w^2+4x^3-y^6,                                            (47)
```

the simple-elliptic weighted sidecar. Therefore the two principal-divisor
rows in `(44)` are sufficient to prove a genuine order-three subgroup but
not to present the full global class group. A nonmonogenic twist would need
additional global class data and compatible involution descent.

## 8. Reproduction and scope

Run

```bash
python3 04-computation/jc2_infinity_component_linear_conic_torus_sextic_thm3932.py
python3 -O 04-computation/jc2_infinity_component_linear_conic_torus_sextic_thm3932.py
```

After platform newlines are normalized to LF, both streams must byte-match
the frozen raw-LF output named in the metadata. The companion verifies the
trace/norm grammar, explicit family, irreducibility,
projective identity and basepoint freedom, rational inverse, full
singularity/delta packet, coefficient Jacobian, cubic and resolvent
normality, Smith relation, non-Cartier local test, and weighted elliptic
sidecar.

This theorem closes the fold-degree classification for the singular conic
with an infinity component. It does not classify every torus sextic, compute
the full class group of `(41)`, produce a nonmonogenic cubic twist, or prove
`JC(2)`.
