---
id: THM-3582
title: "Mixed-cubic target-graph Euler no-go"
status: >
  PROVED + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.  Every
  collision-compatible target graph of exact total degree three and
  a-degree one for the fixed THM-1300 Keller map has irreducible complete
  pullback X with compactly supported Euler characteristic at least two.
  Hence X is not A2 and its defining polynomial has no source-coordinate
  factor.  The proof gives the exact cubic-fibre Euler formula and reduces
  the bound to rho(Delta)+rho(g3)>=5.  Three coefficient chambers are closed
  by finite saturated-ideal and resultant atlases over QQ.  Irreducibility
  follows from THM-3573's complete polynomial Pell-descent classification.
source: kps-s188 / delegated mixed-cubic target-graph attack, 2026-08-21
depends_on:
  - THM-3560-jelonek-euler-gate-monomial-target-shear-no-go
  - THM-3571-quadratic-target-graph-euler-no-go
  - THM-3573-polynomial-target-graph-pell-parameter-descent-classification
related:
  - THM-3570-universal-pell-conic-target-graph-factor-compiler
companion: 04-computation/jacobian_mixed_cubic_target_graph_euler_no_go_thm3582.py
output: 05-knowledge/results/jacobian_mixed_cubic_target_graph_euler_no_go_thm3582.out
script_sha256: 8f6f1329cbb55ba468984f2f5e5a53b2755b0294c357e0a80faa93a527894c84
output_sha256: 55399e31acc692068953b3b2749062337fe23ca2e78ae2bb0dc3aabaedadc00a
hash_basis: LF-normalized bytes
---

# THM-3582 -- mixed-cubic target-graph Euler no-go

**PROVED + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**  The next target-graph
cell after THM-3571 is closed.  The result is a computer-assisted theorem:
the geometric and Euler identities below are displayed algebra, while the
finite exceptional atlases are exact Groebner/resultant computations over
`QQ`, replayed with ordinary and optimized Python.  No independent audit is
claimed here.

All varieties are over `C`; `chi` means compactly supported topological Euler
characteristic.  The fixed Keller map is the map `F:A3->A3` of THM-1300.

## 1. Statement and the complete mixed-cubic cell

Every polynomial of exact total degree three, of degree one in `a`, whose
graph contains the collision value `(-1/4,0,0)`, is uniquely

```text
phi(a,b)=a(Ab^2+Bb+C)+Db^3+Eb^2+Fb+C/4,                (1)
```

with

```text
(A,B,C)!=(0,0,0),                 (A,D)!=(0,0).         (2)
```

The first condition is `deg_a(phi)=1`; the second says that the total degree
is exactly three.  The constant `C/4` is forced by
`phi(-1/4,0)=0`.  Put

```text
T_phi=V(c+phi(a,b)) ~= A2,
X_phi=F^(-1)(T_phi)=V(F3+phi(F1,F2)).                   (3)
```

Then

```text
chi(X_phi)>=2,                                          (4)
```

and `X_phi` is irreducible.  Consequently `X_phi` is not `A2`, the polynomial
`F3+phi(F1,F2)` is not a source coordinate, and it has no source-coordinate
factor.

## 2. Polynomial strict transform

THM-3571 resolves the quadratic discriminant of the Jelonek hypersurface by

```text
a=s(2b-s)/12,
c=4(-b+2s)/(3s^2).                                     (5)
```

On this sheet define

```text
G(s,b)=12s^2(c+phi(a,b))
      =g3(s)b^3+g2(s)b^2+g1(s)b+g0(s),                 (6)

g3=2s^2(As+6D),
g2=s^2(-As^2+2Bs+12E),
g1=-Bs^4+2Cs^3+12Fs^2-16,
g0=s(-Cs^3+3Cs+32).                                    (7)
```

The finite map

```text
(s,b) |-> (s(2b-s)/12,b)                               (8)
```

restricts to a constructible bijection

```text
V(G)  ->  D_phi=(T_phi intersect V(L))_red.             (9)
```

Indeed, the two roots of `s^2-2bs+12a` select the two exchanged factors in
THM-3571's sheet factorization.  Off the ramification diagonal exactly one
factor vanishes; on the diagonal the two roots represent one point.  The
apparent denominator boundary introduces no extra point because

```text
G(0,b)=-16b,                 partial_b G(0,0)=-16.      (10)
```

Thus

```text
chi(D_phi)=chi(V(G)).                                  (11)
```

## 3. Exact cubic-fibre Euler formula

Homogenize the fibre polynomial as

```text
Ghat_s(B,T)=g3 B^3+g2 B^2T+g1 BT^2+g0 T^3.            (12)
```

Its discriminant is

```text
Delta=g2^2 g1^2-4g3 g1^3-4g2^3 g0-27g3^2 g0^2
      +18g3 g2 g1 g0.                                  (13)
```

For a nonzero polynomial `P in C[s]`, let `rho(P)` denote its number of
distinct complex roots.  A nonzero binary cubic is a cube of one linear form
exactly when

```text
K1=g2^2-3g3g1=0,
K2=g1^2-3g2g0=0,
K3=g2g1-9g3g0=0.                                      (14)
```

Let `tau` be the number of `s` satisfying `(14)` for which not all four
coefficients `g3,g2,g1,g0` vanish.  Then

```text
chi(D_phi)=3-rho(Delta)-rho(g3)-tau.                   (15)
```

Here is the complete local proof.  A nonzero projective binary cubic has
three, two, or one distinct roots according as its multiplicity partition is
`1+1+1`, `2+1`, or `3`.  Therefore its projective root count is

```text
3 - 1_{Delta=0} - 1_{binary cubic is a cube}.           (16)
```

The point at infinity is a root exactly when `g3=0`, so subtract
`1_{g3=0}` to get the affine fibre count.  If all four coefficients vanish,
the fibre is a whole `A1`, of Euler characteristic one; `Delta=0` and
`g3=0` subtract two from three, while the definition of `tau` excludes this
zero binary cubic.  Thus `(16)` also gives one in the whole-fibre case.
Euler integration over `A1_s` proves `(15)` without a genericity assumption.

The omitted target curve has coordinates

```text
a=b^2/12,                         c=4/(3b),             (17)
```

with `b!=0`.  Multiplying its graph equation by `12b` gives

```text
f(b)=Ab^5+(B+12D)b^4+(C+12E)b^3+12Fb^2+3Cb+16.        (18)
```

Since `f(0)=16`, its roots are exactly the reduced omitted support, and

```text
chi(T_phi intersect E)=rho(f)<=5.                      (19)
```

Apply THM-3560's exact `3/1/0` Euler integration to `T_phi~=A2`.  Equations
`(15)` and `(19)` give

```text
chi(X_phi)
 =3-2chi(D_phi)-rho(f)
 =-3+2(rho(Delta)+rho(g3)+tau)-rho(f).                 (20)
```

It therefore suffices to prove the stronger root-support inequality

```text
rho(Delta)+rho(g3)>=5.                                 (21)
```

## 4. Discriminant packet

Exact expansion of `(13)` gives

```text
Delta=s^2 R(s),                    deg R<=14.           (22)
```

Write `r_j=[s^j]R`.  The low and high rows used by the proof are

```text
r0 =196608D,
r1 =32768A,
r2 =36864(E^2-12DF),
r3 =-12288(6AF-BE+6CD+108DE),                          (23)

r14=A^2(B^2-4AC),
r13=4AB(B^2-4AC),                                     (24)

r12=4(3A^3C-6A^2BF-8A^2C^2+36A^2CE
        -2AB^2C-6AB^2E-54ABCD+B^4+12B^3D),            (25)

r11=4(32A^3+9A^2BC+120A^2CF-48AB^2F
        -16ABC^2-12ABCE-216AC^2D
        +4B^3C+12B^3E+36B^2CD).                       (26)
```

The remaining exact rows are generated directly by `(7)` and `(13)` in the
companion.  We now prove `(21)` in the three exhaustive chambers determined
by `(A,D)`.

## 5. Case I: `A!=0`, `D!=0`

Here

```text
rho(g3)=rho(2s^2(As+6D))=2.                            (27)
```

Since `R(0)!=0`, `rho(Delta)=1+rho(R)`.  It is enough to show
`rho(R)>=2`.

The high packet leaves only

```text
deg R in {14,12,11}.                                  (28)
```

Indeed, if `B^2-4AC!=0`, then `r14!=0`.  On the tangent
`C=B^2/(4A)`,

```text
r12=3B(A^2B-8A^2F+4ABE-2B^2D).                        (29)
```

If `(29)` vanishes, either `B=0`, or it solves `F`; in both branches `(26)`
reduces to

```text
r11=128A^3!=0.                                         (30)
```

Suppose `R` had one distinct root.  Since `R(0)=196608D`,

```text
R=196608D(1+ps)^M,
p!=0,                    M in {14,12,11}.              (31)
```

The `s^1` coefficient gives `A=6MDp`.  The weight change

```text
D=dp^4, B=up^4, C=vp^3, E=ep^3, F=fp^2,
x=ps                                                       (32)
```

turns `(31)` into the same identity with `p=1`, `s=x`, `D=d`, and
`A=6Md`: every term `[s^j]R` has parameter weight `p^(j+4)`.

The normalized finite atlas is:

| `M` | exact elimination certificate |
|---:|---|
| 14 | The top ratio gives `B=294D`; `r14` and `r12` solve `C,F`.  The `r3` row is linear in `E`; its leading/constant gcd in `D` is one.  After solving it, the `r4,r5` gcd in `D` is one. |
| 12 | Tangency gives `C=B^2/(4A)` and `r3` solves `F`.  With `w=B/D`, away from `w=0,-54`, `r12` solves `E`; the three resultants of `r4` against `r5,r6,r7` have radical supported on `w(w+54)`.  The direct `w=0` row and saturated `w=-54` ideal are empty. |
| 11 | Again `C=B^2/(4A)` and `r3` solves `F`.  The `B=0` ideal is unit.  For `B=Dw`, away from `w=-54`, three resultants have radical supported on `w=0`; at `w=-54`, `r12` forces `D=880/891` and the remaining gcd in `E` is one. |

All gcds, resultants, and exceptional ideals are over `QQ`; a saturated
unit ideal means a literal Nullstellensatz contradiction after inverting the
listed nonzero chamber factors.  Hence `(31)` is impossible,
`rho(R)>=2`, `rho(Delta)>=3`, and `(21)` follows from `(27)`.

## 6. Case II: `A=0`, `D!=0`

Now `rho(g3)=1`, `R(0)=196608D`, and `R'(0)=0`.  The exact high-degree
stratification is

```text
r12=4B^3(B+12D),

B=-12D:       r11=-6912D^3(C+12E),
C=-12E:       r10=-82944D^3F,
F=0:          r9 = 248832D^3E,
E=0:          r8 =-110592D^3,

B=0:          r10=-3888C^2D^2.                        (33)
```

Because `(B,C)!=(0,0)`, this proves

```text
deg R in {8,9,10,11,12}.                              (34)
```

Thus `R` is nonconstant.  It cannot have one root: a one-root polynomial
with nonzero constant term is `c(1+ps)^M`, whose derivative at zero is
nonzero, contrary to `R'(0)=0`.

If `R` had exactly two roots, order their multiplicities `m<=n`.  The zero
derivative forces the normalized form

```text
R=196608D(1+ps)^m(1-(m/n)ps)^n,
p!=0,                    8<=m+n<=12.                  (35)
```

The weights `(32)` set `p=1`.  There are exactly

```text
4+4+5+5+6=24.                                          (36)
```

multiplicity pairs.  The full finite atlas is:

| chamber | multiplicity cells | certificate |
|---|---:|---|
| `E=0` | 24 | after `r2,r3` solve `F,C`, every saturated residual ideal is unit |
| `E!=0`, degree 8 | 0 | `(33)` would force `E=0` |
| `E!=0`, degree 9 | 4 | `B=-12D,C=-12E,F=0`; four unit ideals |
| `E!=0`, degree 10, `B=0` | 5 | five unit ideals |
| `E!=0`, degree 10, tangent | 5 | `B=-12D,C=-12E`; five unit ideals |
| `E!=0`, degree 11 | 5 | `B=-12D`; five unit ideals after inverting `C+12E` |
| `E!=0`, degree 12, `m<n` | 5 | put `w=B/D`; generic resultants plus one quadratic exception |
| `E!=0`, degree 12, `(m,n)=(6,6)` | 1 | direct symmetric exceptional ideal |

For completeness, the degree-twelve exceptional polynomial is

```text
q(w)=w^2-81w-972.                                     (37)
```

For each of the five asymmetric pairs, `r11` solves `E` away from `q=0`;
the resultants of `r12` against `r4,r5,r10` have gcd

```text
nonzero rational * w^15.                              (38)
```

Degree twelve excludes `w=0`.  On `q=0`, the saturated ideal generated by
`r11,r12,q` is unit.  In the symmetric `(6,6)` cell, `r11` is a nonzero
rational multiple of

```text
D^3 E w^2 q(w),                                       (39)
```

and the ideal from `r9,r10,r12,q`, saturated by `DEw`, is unit.  This closes
all 24 pairs in `(35)`.  Therefore `rho(R)>=3`,
`rho(Delta)=1+rho(R)>=4`, and `(21)` follows because `rho(g3)=1`.

## 7. Case III: `A!=0`, `D=0`

Here `rho(g3)=1`, and

```text
R=sS,                     S(0)=32768A!=0.              (40)
```

The high packet gives

```text
deg S in {13,11,10}.                                  (41)
```

It remains to exclude one- and two-root `S`.

### 7.1 One root

Write

```text
S=32768A(1+ps)^M,
M in {13,11,10}.                                      (42)
```

The first coefficient gives `A=9E^2/(8Mp)`, so `E!=0`; the weighted change
sets `p=1`.  The second coefficient solves `B` in terms of `F`.  The exact
atlas is:

| `M` | certificate |
|---:|---|
| 13 | the top ratio and leading row solve `B,C,F`; the gcd of the remaining coefficient rows is one |
| 11 | `C=B^2/(4A)`; the `R`-rows `r12,r4,r5` generate a saturated unit ideal |
| 10 | `C=B^2/(4A)`; the `R`-rows `r11,r12,r4` generate a saturated unit ideal |

Thus `S` cannot have one root.

### 7.2 Two roots with `E=0`

Normalize two distinct roots as

```text
S=32768A(1+s)^m(1+rs)^n,
r!=0,1.                                                (43)
```

Put `h_j=[s^j](1+s)^m(1+rs)^n`.  If `E=0`, the first coefficient forces

```text
h1=m+nr=0.                                             (44)
```

For `m+n in {10,11,13}` there are `5+5+6=16` pairs.  Substitute
`r=-m/n`; the second row solves `F`, and all sixteen residual ideals,
saturated by `A`, are unit.

### 7.3 Two roots with `E!=0`

Now `h1!=0`.  Normalize by one nonzero scale `u`:

```text
E=(8h1/9)u,                    A=(8h1/9)u^2.           (45)
```

For `M=13`, the high ratio and leading row solve

```text
B=A(mr+n)/(4r),
F=(B/u-(8/3)h2)/6,
C=(B^2-32768r^n/A)/(4A).                              (46)
```

For each of the six multiplicity pairs, resultants of the first residual row
against the next four have radical supported only on

```text
r(m+nr)(mr+n).                                        (47)
```

The factors `r` and `m+nr` are forbidden by `(43)` and `(45)`.  Directly
substituting the apparent branch `mr+n=0` gives coefficient gcd one in every
one of the six cells.

For `M=10`, the tangent splits into `B=0` and a generic branch.  The `B=0`
branch gives five unit ideals.  In the generic branch,

```text
F=-(8h2/3)(4+u)/(16+6u),
B=-8u(8h2/3)/(16+6u),
C=B^2/(4A).                                            (48)
```

For all five multiplicity pairs, resultants of the leading row against four
other rows have radical supported only on `r(m+nr)`, hence only on forbidden
loci.  The omitted denominator `u=-8/3` forces `h2=0`; there the leading
row is

```text
128A^3=32768Ar^n,                                     (49)
```

and the gcd of its difference with `h2` is one in each of the five cells.

For `M=11`, put

```text
k=8h1/9,
V=6F+(8/3)h2,
A=ku^2, E=ku, B=uV, C=V^2/(4k).                       (50)
```

After removing the universal factor `u^2`, the `R`-row `r4` is linear in
`u`; its coefficient is a nonzero rational multiple of

```text
-2368h1^2.                                             (51)
```

Solve it for `u`.  For each of the five multiplicity pairs, the resultants
of `r12` against `r5,r6,r10,r11` have radical supported only on
`r(m+nr)`.  These loci are forbidden, so all five cells are empty.

This exhausts one and two roots for `S`.  Hence `rho(S)>=3`, so by `(40)`

```text
rho(Delta)=1+rho(S)>=4.                                (52)
```

Together with `rho(g3)=1`, this proves `(21)` in Case III.

## 8. Euler conclusion and irreducibility

All three cases give `(21)`.  Since `tau>=0` and `rho(f)<=5`, equation `(20)`
now yields

```text
chi(X_phi)>=-3+2*5-5=2.                               (53)
```

It remains to justify that the complete pullback is irreducible.  THM-3573
classifies every reducible nonzero polynomial target graph: reducibility of
the core cubic is equivalent to

```text
phi=4H(1+bH+4aH^2),                 0!=H in C[a,b].    (54)
```

Moreover,

```text
deg_a(phi)=3deg_a(H)+1.                               (55)
```

Our `deg_a(phi)=1`, so `(55)` forces `H in C[b]`.  If `H` is constant,
`(54)` has total degree at most one.  If `deg_b(H)=r>=1`, the term

```text
16aH^3                                                   (56)
```

has total degree `3r+1>=4` and cannot cancel against the two terms without
`a`.  Thus no exact total-degree-three row `(1)` is reducible.

The core cubic over `C(a,b)` is therefore irreducible.  Because `F` is
quasi-finite, every irreducible component of the two-dimensional pullback
would dominate `T_phi`; generic-core irreducibility leaves exactly one
component.  Hence `X_phi` is irreducible.  If its defining polynomial had a
source-coordinate factor, irreducibility would make that factor the entire
hypersurface, contradicting `(53)` and `chi(A2)=1`.  This proves the theorem.

## 9. Hostile controls and sharp boundaries

The exact companion includes four controls chosen to attack a different
load-bearing part of the proof.

1. For `phi=b^3-12ab`,

   ```text
   Delta=-12288s^2(9s^8+132s^4-16),       f=16,
   (rhoDelta,rhoG3,tau,rhoF,chiD,chiX)=(9,1,0,0,-7,17).
   ```

2. For `phi=ab^2`,

   ```text
   Delta=128s^3(s^10-718s^5+256),
   (rhoDelta,rhoG3,tau,rhoF,chiD,chiX)=(11,1,0,5,-9,16).
   ```

   This row attains the omitted-support cap `rho(f)=5`.

3. For

   ```text
   (A,B,C,D,E,F)=(1/2,1/4,-16,0,0,193/48),
   ```

   one has `G(1,b)=b^3`.  Here `tau=1` and

   ```text
   (rhoDelta,rhoG3,tau,rhoF,chiD,chiX)=(12,1,1,5,-11,20).
   ```

   Thus the cube-fibre correction in `(15)` is genuinely necessary.

4. The scope boundary is exact.  Pure `b^3` has `deg_a=0` and `chi(X)=5`.
   On the reducible side, THM-3573 with `H=-b/2` gives

   ```text
   phi=-2ab^3+b^3-2b,                                  (57)
   ```

   of total degree four.  Reducible rows therefore appear immediately after
   the cubic theorem's degree boundary, but not inside it.

These controls do not assert that the numerical lower bound two is attained;
they show that the omitted-support cap, the `tau` correction, the pure-`b`
exclusion, and the irreducibility degree boundary are all real rather than
decorative.

## 10. Exact verification and scope

Run

```bash
python3 04-computation/jacobian_mixed_cubic_target_graph_euler_no_go_thm3582.py
python3 -O 04-computation/jacobian_mixed_cubic_target_graph_euler_no_go_thm3582.py
```

The ordinary and optimized transcripts agree.  The companion verifies
`(6)--(7)`, the discriminant formula `(13)`, every displayed coefficient in
`(23)--(26)`, the omitted polynomial `(18)`, the eight cubic fibre types,
all three root-support atlases, every exceptional ideal/resultant support,
and all four hostiles.  Truth gates use explicit exceptions rather than
Python `assert`, so optimization removes no check.

The theorem closes only collision-compatible target graphs of exact total
degree three with `deg_a(phi)=1` for the fixed THM-1300 map.  It does not
close arbitrary cubic target coordinates, nongraph targets, the quartic
Pell family, `JC(2)`, or `DC(2)`.  Independent proof and hostile audit remain
pending.

**QED.**
