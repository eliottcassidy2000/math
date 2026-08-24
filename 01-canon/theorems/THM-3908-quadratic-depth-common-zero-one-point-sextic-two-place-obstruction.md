---
id: THM-3908
title: "Quadratic-depth binary cubics cannot realize a nonmonogenic one-place sextic"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
  independent hostile audit.  Let a binary cubic over k[A,C] have
  coefficient degree at most two, irreducible generic fibre, irreducible
  degree-six discriminant, and globally nonmonogenic Delone--Faddeev order.
  If the degree-six discriminant form is one linear sixth power, then the
  unique projective infinity point has at least two normalization places.
  For a common coefficient zero, every irreducible row is exactly in the
  triple-root THM-3906 grammar and has two places.  Constants do not provide
  an escape: a triple-root leading row is reducible or at least two-place; a
  fixed double root is reducible or monogenic; and the sole moving-double-
  root normal form becomes reducible when its sextic top is a pure power.
  Coefficient depth three, lower-degree unit-ideal discriminants, Keller
  realization, and JC(2) remain open.
source: jc_degree6_one_place / post-THM-3906 leading-stratum classification, 2026-08-23
audit: >
  SELF-AUDITED proof candidate.  The proof removes every coefficient base
  divisor before lifting the leading row to the normalization of the binary-
  cubic discriminant.  It distinguishes the intrinsically marked repeated
  root from the simple root, includes finite triple-root specializations,
  and checks all possible pullback bidegrees.  The exact companion freezes
  the discriminant tangent, twisted-cubic and double-root degree ledgers,
  the common-zero two-edge packet, the constant-row gamma/beta/alpha Newton
  chain, and the moving-root A^2-times-quartic identity in 42 active gates.
  Normal and optimized replays byte-match the frozen transcript.  An
  independent hostile audit remains required before promotion.
related:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3889-maximally-confluent-quadratic-binary-cubic-two-place-obstruction
  - THM-3890-universal-quintic-common-zero-resolvent-class-group-dichotomy
  - THM-3906-degree-six-common-zero-normal-cubic-two-place-boundary
  - THM-3907-unit-ideal-nonmonogenic-cubic-six-place-boundary
script: 04-computation/jc2_quadratic_depth_one_point_sextic_obstruction_thm3908.py
output: 05-knowledge/results/jc2_quadratic_depth_one_point_sextic_obstruction_thm3908.out
script_sha256: 8fd4fefd95bf729934f8e1fdd6b8203d3f17a3025ac310f165cb8d796aeb9231
output_sha256: 175a7a753fac08ed4ab13f1e7aec8fdc8a936d3724aeec60acc9badd44ad6f28
semantic_sha256: dbf51a2ab95c0c45304577625fe3d17609522eb74d0b4df26a9e1cf8dc5a18a7
hash_basis: raw LF bytes
---

# THM-3908 -- quadratic coefficient depth cannot pay the one-place sextic invoice

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
independent hostile audit.**

Work over an algebraically closed field `k` of characteristic zero.  Put
`R=k[A,C]` and let

```text
Phi=aU^3+bU^2V+cUV^2+dV^3,                 a,b,c,d in R,  (1)
```

where every coefficient has total degree at most two.  Write

```text
Phi=Phi_2+Phi_1+Phi_0                                  (2)
```

for the homogeneous target-degree rows, and let `Delta=Disc(Phi)`.  Assume:

1. `Phi` is irreducible over `k(A,C)`;
2. the Delone--Faddeev cubic order belonging to `Phi` is globally
   nonmonogenic, equivalently its binary index form represents no element of
   `k*`;
3. `Delta` is irreducible of total degree six; and
4. its degree-six part is

   ```text
   Delta_6=lambda L^6,                    lambda in k*,  (3)
   ```

   for a nonzero target-linear form `L`.

Then the projective curve `V(Delta)` has one infinity support point, but its
normalization has at least two points above it.  In particular no normal
finite-flat, globally nonmonogenic `S3` cubic order of coefficient depth at
most two realizes a one-place sextic discriminant.

There is a sharper common-zero statement.  If `Phi_0=0`, then after target
and binary-variable linear changes every irreducible row satisfying `(3)` is

```text
a=C^2+alpha A+alpha_1 C,
b=beta A+beta_1 C,
c=gamma A+gamma_1 C,
d=C,                                                       (4)
```

with

```text
gamma^2(beta^2-4alpha gamma)!=0.                           (5)
```

Its unique infinity point has exactly two smooth normalization branches,
of `z`-orders two and four in the chart `x=C/A,z=Z/A`.

The proof also identifies the exact constant-row boundary.  A triple-root
leading row remains reducible or at least two-place.  On the smooth
double-root stratum, a fixed repeated root makes the full row reducible or
the order monogenic.  The only leading row whose repeated root moves is

```text
Phi_2=(AU+CV)^2V,                                         (6)
```

and tangent cancellation forces

```text
Phi_1=(AU+CV)Q_2(U,V),             Phi_0=Q_3(U,V).         (7)
```

Condition `(3)` then forces `V|Q_2,Q_3`, so `V|Phi`, contrary to generic
irreducibility.  Thus constants close rather than open the quadratic-depth
escape.

## 1. The two leading strata and the base-divisor ledger

Since `deg Delta=6`, the target-degree eight and seven pieces of the
discriminant vanish:

```text
Disc(Phi_2)=0,
d Disc_(Phi_2)(Phi_1)=0.                                  (8)
```

The first identity places the generic `Phi_2` on the discriminant
hypersurface of projective binary cubics.  Its singular locus is the twisted
cubic of triple roots.  Therefore the generic leading row is either
double-plus-simple or triple-root.  It cannot vanish identically because
then `Delta` would have degree at most four.

The word "generic" here hides no base-point assumption.  Let `H` be the gcd
of the four homogeneous quadratic coefficients of `Phi_2`, with
`g=deg H`.  Dividing by `H` leaves four degree

```text
e=2-g                                                        (9)
```

forms with no common zero on `P1_[A:C]`, hence a genuine morphism to
projective coefficient space.  Every isolated zero of the original row is
accounted for by `H`; none is silently discarded.

### 1.1. Triple-root leading row

If the residual map lands in the twisted cubic, it factors through

```text
nu_3:P1 -> P3,                     ell |-> ell^3.           (10)
```

Because `nu_3^*O(1)=O(3)`, a nonconstant root map of degree `r` would give

```text
e=3r.                                                        (11)
```

But `e` belongs to `{0,1,2}`.  Hence `r=0`: the binary linear factor is
fixed, including across the zeros removed into `H`.  Thus

```text
Phi_2=H_2(A,C) ell(U,V)^3.                                 (12)
```

After a binary change take `ell=U`.

### 1.2. Double-plus-simple leading row

The normalization of the binary-cubic discriminant surface is

```text
P1_ell x P1_m -> D,                (ell,m) |-> ell^2 m,    (13)
```

and the hyperplane bundle pulls back to `O(2,1)`.  The double root `ell` is
intrinsic on the smooth locus: ramification type `2+1` distinguishes it from
the simple root `m`, so there is no root-swapping ambiguity.  The two generic
root maps are rational maps from the nonsingular complete source curve to
projective lines; properness extends them across the finitely many points
where the row becomes triple-root.  Substitution then extends the lift to
`(13)`.  These nonclosed leading strata therefore add no hidden monodromy
case.

If the two root maps have degrees `r,s`, then

```text
2r+s=e<=2.                                                  (14)
```

The complete list is

```text
(e,r,s)=(0,0,0),(1,0,1),(2,0,2),(2,1,0).                  (15)
```

Thus either the repeated root is fixed, or `(e,r,s)=(2,1,0)` and it moves
linearly while the simple root is fixed.  This exhausts common factors,
base points, ramified root maps, and finite collisions `ell=m`.

## 2. The smooth tangent equation

At a generic double-plus-simple row choose binary coordinates over
`k(A/C)` so that

```text
Phi_2=U^2(rU+sV),                         s!=0,             (16)
Phi_1=a_1U^3+b_1U^2V+c_1UV^2+d_1V^3.
```

The coefficient of `epsilon` in

```text
Disc(Phi_2+epsilon Phi_1)
```

is exactly

```text
-4s^3 d_1.                                                   (17)
```

Hence the degree-seven cancellation in `(8)` is equivalent to

```text
U divides Phi_1.                                            (18)
```

Intrinsically, `Phi_1` vanishes at the repeated root.  This is the tangent
hyperplane to the smooth discriminant, not merely a dimension count.

If the repeated root is fixed, `(18)` and `ell|Phi_2` leave two cases.  If
`ell|Phi_0`, then `ell|Phi` and the generic cubic is reducible.  Otherwise
evaluate the binary index form at the fixed point `ell=0`.  The value is
the nonzero scalar `Phi_0|_(ell=0)`, so the index form represents a unit and
the cubic order is monogenic.  This proves the fixed-root
reducible-or-monogenic dichotomy.

In particular, when `Phi_0=0`, every smooth-leading row is reducible.  This
is why the common-zero branch must enter the triple-root stratum.

## 3. The triple-root UFD reduction

In coordinates `(12)`, write the coefficient of `V^3` in `Phi_1` as the
linear form `d_1`.  Direct expansion gives

```text
Delta_6=-27 H_2^2 d_1^2.                                   (19)
```

Together with `(3)`, unique factorization in `k[A,C]` forces

```text
H_2=rho L^2,                     d_1=sigma L,              (20)
```

with `rho,sigma` nonzero.  Target and diagonal binary changes normalize
`L=C`, `rho=sigma=1`.  The common-zero row is then exactly `(4)`.

At `C=0`, its discriminant is

```text
Delta(A,0)=gamma^2(beta^2-4alpha gamma)A^4.                (21)
```

If the scalar coefficient vanishes, then `C|Delta`, contradicting the
assumed irreducibility.  Hence `(5)` holds.

Let `Delta^h(A,C,Z)` be the degree-six homogenization and put

```text
h(x,z)=Delta^h(1,x,z).
```

Its lower boundary has vertices

```text
(0,2), (2,1), (6,0),                                      (22)
```

with coefficients

```text
[z^2]h=gamma^2(beta^2-4alpha gamma),
[x^2z]h=-4gamma^3,
[x^6]h=-27.                                                (23)
```

Every support exponent satisfies both

```text
i+2j>=4,                         i+4j>=6.                  (24)
```

The two primitive compact edges therefore produce one branch each, of
orders two and four.  Their different orders make the normalization places
distinct.  This proves the sharper common-zero statement.

## 4. Why constants do not merge the triple-root places

Retain the normalization `(20)` but allow arbitrary constants:

```text
a=C^2+alpha A+alpha_1 C+a_0,
b=beta A+beta_1 C+b_0,
c=gamma A+gamma_1 C+c_0,
d=C+d_0.                                                    (25)
```

The top form remains `-27C^6`.  In the same local equation `h`, the
following conditional coefficient chain is exact:

```text
[x^6]h=-27,
[x^2z]h=-4gamma^3,
gamma=0             ==> [xz^2]h=-4beta^3,
gamma=beta=0        ==> [x^4z]h=-54alpha.                  (26)
```

Irreducibility forces an `x`-independent term in `h`; otherwise `x|h`, so
`C|Delta`.  Let `(0,j_0)` be the lowest such exponent.

- If `gamma!=0`, then `j_0>=2`, and `(2,1)` lies strictly below the chord
  from `(0,j_0)` to `(6,0)` because `1<2j_0/3`.
- If `gamma=0,beta!=0`, then `j_0>=3`, and `(1,2)` lies strictly below that
  chord because `2<5j_0/6`.
- If `gamma=beta=0,alpha!=0`, then `j_0>=4`, and `(4,1)` lies strictly below
  that chord because `1<j_0/3`.

In every case the lower Newton boundary has at least two compact edges.
Over the declared algebraically closed characteristic-zero field, each edge
has a nonzero Newton--Puiseux root; distinct edge slopes give distinct
normalization branches.

The only terminal seam is

```text
alpha=beta=gamma=0.                                        (27)
```

Then every coefficient in `(25)` is independent of `A`, so `Delta` is a
degree-six polynomial in `C` alone.  It is reducible over `k`.  This proves
that constants never turn a triple-root quadratic row into an irreducible
one-place sextic.

## 5. The moving repeated root is a forced global factor

It remains to close the sole moving case in `(15)`.  After target and binary
changes, the degree-one repeated-root map is an isomorphism and the fixed
simple factor is `V`, giving `(6)`.  Gauss's lemma applied to the primitive
linear form `AU+CV` upgrades divisibility over `k(A,C)[U,V]` in `(18)` to
divisibility in `k[A,C,U,V]`.  Target homogeneity then forces `(7)` with
constant forms

```text
Q_2=pU^2+qUV+rV^2,
Q_3=sU^3+tU^2V+uUV^2+vV^3.                                (28)
```

An exact discriminant expansion gives

```text
Delta_6=A^2[
 (r^2-4v)A^4+(-2qr+4u)A^3C
 +(2pr+q^2-4t)A^2C^2+(-2pq+4s)AC^3+p^2C^4].              (29)
```

If this nonzero form is a sixth power of one linear form, its visible factor
`A^2` forces that form to be proportional to `A`.  Successively comparing
the `C^4,C^3,C^2,C` coefficients gives

```text
p=0,                 s=0,
t=q^2/4,              u=qr/2.                              (30)
```

Consequently

```text
Q_2=V(qU+rV),
Q_3=V((q^2/4)U^2+(qr/2)UV+vV^2).                          (31)
```

Equations `(6),(7),(31)` show `V|Phi`, contradicting generic irreducibility.
This closes the last constant-row loophole.

## 6. Normal `S3` consequence, design boundary, and replay

Under the theorem's irreducibility hypotheses, the associated
Delone--Faddeev algebra is finite free.  The irreducible discriminant has
height-one valuation one, so the DVR index formula makes the order maximal
in codimension one; finite freeness gives `S2`, hence normality.  The generic
cubic is irreducible and its discriminant is nonsquare, so its Galois group
is `S3`.  Thus the obstruction is genuinely upstream of the normality and
connected `C3`-resolvent invoices: quadratic depth fails at place
confluence, reducibility, or monogenicity.

The first coefficient depth at which a moving triple root is possible is
three, because the twisted-cubic pullback degree can then equal three.
Therefore the next economical search should start with cubic coefficient
rows, not another quadratic perturbation.  This theorem does not exclude
unit-ideal discriminants of degree at most five, cubic or higher coefficient
depth, a different branch degree, Keller realization, or `JC(2)`.

Reproduce the exact packet with

```bash
python3 04-computation/jc2_quadratic_depth_one_point_sextic_obstruction_thm3908.py
python3 -O 04-computation/jc2_quadratic_depth_one_point_sextic_obstruction_thm3908.py
```

Both streams must byte-match the frozen output named in the metadata.
