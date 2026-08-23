---
id: THM-3721
title: "Automorphic Cohn single-right-shear nonentry"
status: >
  PROVED + VERIFIED-EXACT.  A polynomial-ring automorphic image of Cohn's
  non-elementary matrix has an infinite-dimensional family of elementary
  decorations in which one row really is a gradient.  Nevertheless, after
  either one arbitrary upper or one arbitrary lower right shear and either
  order of two alternating left shears, the matrix is never a Jacobian
  matrix.  Nontriangular right data prevent the first closure.  Every
  surviving triangular parameter produces a Broughton component, possibly
  plus a charge-separated pure-coordinate debt, and THM-3716 excludes the
  second closure.  This closes the complete single-right-shear
  automorphic-Cohn cell, not two right factors, general non-elementary
  classes, or JC(2).
source: root / 2026-08-22
audit: >
  INDEPENDENTLY HOSTILE-AUDITED.  The exact companion
  reconstructs the Cohn normalization, all four closure gates, both full
  triangular survivor formulas, the charge-coordinate homogeneous kernels,
  the forbidden Laurent layer, and the inherited Broughton debt.  An
  independent bounded scout checks 224 nontriangular cells, and normal and
  optimized runs byte-match the frozen transcript.
depends_on:
  - THM-3652-wright-elementary-jacobian-criterion-reduced-word-reproof
  - THM-3653-cohn-factorial-repair-and-weighted-rectangle-holonomy
  - THM-3716-monomial-broughton-hamiltonian-obstruction-family
related:
  - THM-3709-cohn-alternating-two-by-two-elementary-decoration-nonentry
  - THM-3719-cohn-aligned-identity-shortening-laurent-nonentry
script: 04-computation/jc2_automorphic_cohn_one_right_shear_thm3721.py
output: 05-knowledge/results/jc2_automorphic_cohn_one_right_shear_thm3721.out
script_sha256: e925b4772405df0c8e499158df5b9653defcc65e1e784d07029b9e7110717d0b
output_sha256: ecb76b1f84d43a201eb8e76023d25d79dc6a58dd278f6dda2d59b6db787d6811
semantic_sha256: 8933ed6981486aa546d29f3b63b997a18452446cb296343c837a78e1c7fe9852
hash_basis: raw LF bytes
---

# THM-3721 -- no single right shear repairs the automorphic Cohn core

**PROVED + VERIFIED-EXACT.**  This theorem records a positive first-stage
construction, not just another failed coefficient search.  The construction
puts an honest closed row into Cohn's nontrivial `SL_2/E_2` class for every
polynomial in either source variable.  The reason it still fails is exact:
the closed rows retain an isolated Broughton Hamiltonian-cokernel sector.

Work over a characteristic-zero field `k`, and put

```text
M_0 = [ 4T^2       2XT-1 ]
      [ 1+2XT      X^2   ].                              (1)
```

If `C(X,Y)` is Cohn's matrix, then

```text
M_0 = [0 -1;1 0] C(X,2T),                 det M_0=1.    (2)
```

The substitution `Y -> 2T` is a ring automorphism, and the constant rotation
in `(2)` is elementary.  Cohn's nonentry therefore gives
`M_0 notin E_2(k[X,T])`.

For arbitrary `v,f,g,h in k[X,T]`, set

```text
M(v)=M_0 E_+(v).                                         (3)
```

Neither of the alternating two-left decorations

```text
E_+(f)E_-(h)M(v),             E_-(g)E_+(h)M(v)          (4)
```

is a Jacobian matrix.

The same conclusion holds for the lower right orientation.  For arbitrary
`u,f,g,h in k[X,T]`, neither

```text
E_+(f)E_-(h)M_0E_-(u),       E_-(g)E_+(h)M_0E_-(u)    (4-)
```

is a Jacobian matrix.  Sections 1--4 prove the upper orientation; Section 5
proves the lower orientation.

## 1. The two exposed-row equations

Write the rows of `M(v)` as `alpha_v,beta_v`.  Explicitly,

```text
alpha_v=(4T^2, 2XT-1+4T^2v),
beta_v =(1+2XT, X^2+(1+2XT)v).                          (5)
```

For `curl(P,R)=P_T-R_X`,

```text
curl(alpha_v)=6T-4T^2v_X,
curl(beta_v) =-(1+2XT)v_X-2Tv.                          (6)
```

The first word in `(4)` exposes `beta_v+h alpha_v`.  Its closure is
equivalent to

```text
4T^2 h_T-(2XT-1+4T^2v)h_X +(6T-4T^2v_X)h
       =2Tv+(1+2XT)v_X.                                (7)
```

The second word exposes `alpha_v+h beta_v`.  Its closure is equivalent to

```text
(1+2XT)h_T-[X^2+(1+2XT)v]h_X
 -[(1+2XT)v_X+2Tv]h +6T-4T^2v_X=0.                    (8)
```

These are necessary first-stage equations before the other left parameter
can act.

## 2. Nontriangular right data cannot close either row

Suppose `p=deg_X v>=1`, and write

```text
v=a(T)X^p+lower,                  a!=0.                 (9)
```

First consider `(7)`.  If `h!=0`, put
`m=deg_X h` and `h=c(T)X^m+lower`.

- If `p>=2,m>=2`, the unique top term is

  ```text
  -4T^2(p+m)a c X^(p+m-1),                              (10)
  ```

  above both the right side and every term not containing `v`.
- If `p>=2,m=1`, comparison in degree `p` forces
  `c=-1/(2T)`, which is not a polynomial.  If `m=0`, the right side has
  larger `X`-degree than the left side.
- Let `p=1`.  For `m>=2`, the degree-`m` coefficient is

  ```text
  4T^2c'+(6-2m)Tc-4T^2(m+1)ac=0.                      (11)
  ```

  The last summand has strictly larger `T`-degree than the first two.  For
  `m=1`, the same coefficient equation, now including the right side, is

  ```text
  Tc'+c-2Tac=a.                                        (12)
  ```

  Thus `z=1+2Tc` would satisfy `z'=2az`; no polynomial with `z(0)=1`
  satisfies that equation for nonzero polynomial `a`.  The case `m=0`
  already fails in the `X` coefficient of the right side.

The case `h=0` is included in the last degree comparison.  Hence `(7)` has
no solution when `v_X!=0`.

For `(8)`, the product terms in top `X`-degree combine to

```text
-2T(p+m+1)ac X^(p+m).                                  (13)
```

If `p>=2`, `(13)` is uniquely highest.  If `p=1,m>=1`, its collision with
the unperturbed terms gives

```text
2Tc'-mc=2T(m+2)ac,                                     (14)
```

which is impossible by `T`-degree.  If `p=1,m=0`, the top equation is
`c'=2ac`, again impossible for nonzero polynomial `a`; `c=0` would require
`a=3/(2T)`.  Thus `(8)` also has no polynomial solution when `v_X!=0`.

The important point is structural: allowing arbitrary `X`-dependence does
not deform the Broughton debt.  It destroys the first integrability condition
one stage earlier.

## 3. The complete lower-row survivor is one Broughton orbit

It remains to put `v=b(T)`.  There is a unique polynomial `A(T)`, divisible
by `T`, satisfying

```text
A+2TA'=2Tb.                                             (15)
```

Indeed, for `b=sum b_jT^j`,

```text
A=sum [2b_j/(2j+3)]T^(j+1).                            (16)
```

Set

```text
H=A/(2T),                 L=X+A(T),
Q=L+TL^2.                                               (17)
```

A direct differentiation gives the positive construction

```text
dQ=beta_b+H alpha_b.                                   (18)
```

This is the only solution of `(7)`.  To see uniqueness, subtract `(18)` and
write the difference as `K(L,T)`.  Equation `(7)` becomes

```text
[4T^2 partial_T +(1-2LT)partial_L +6T]K=0.             (19)
```

Decompose by charge `q=i-j` on monomials `L^iT^j`.  A charge-`q` component
has the form `L^q F(z)`, `z=LT`, with `F` divisible by `z^(-q)` when `q<0`.
The left side of `(19)` is

```text
L^(q-1){z(1+2z)F'+[q+(6-2q)z]F}.                      (20)
```

If `F!=0`, with valuation `nu` and degree `d`, its lowest and highest terms
in `(20)` force simultaneously

```text
nu=-q,                         d=q-3.                  (21)
```

For `q>0` the first equality is impossible; for `q<=0` the second is
impossible.  Thus `K=0`.

After `(18)`, the other row can close only if some `f` satisfies

```text
J(f,Q)=curl(alpha_b)=6T.                               (22)
```

The triangular source automorphism `(X,T)->(L,T)` has Jacobian one, and
turns `Q` into `L+TL^2`.  THM-3716 says exactly that `(22)` has no polynomial
solution.  Therefore the first word in `(4)` always fails at its second
stage—even though its first stage succeeds for every `b(T)`.

## 4. The opposite survivor is the same obstruction

For `(8)`, continue with `v=b(T)`.  If `deg b>=1` and `h` has top `T`-term
`h_N(X)T^N`, then the unique coefficient in degree `N+deg(b)+1` is

```text
-2 b_top [h_N+X h_N'].                                 (23)
```

It is nonzero in characteristic zero.  If `b=0`, equation `(8)` is

```text
J(h,X+TX^2)=6T,                                        (24)
```

already impossible by THM-3716.  Finally let `b=c in k*`.  The same
top-`T` comparison forces `h in k[X]`; separating the constant and linear
`T` coefficients then gives uniquely

```text
h=3/c.                                                 (25)
```

With

```text
L=X+(2c/3)T,
Q=(3/c)(L+TL^2),                                       (26)
```

one checks `dQ=alpha_c+(3/c)beta_c`.  The remaining row has debt

```text
curl(beta_c)=-2cT.                                     (27)
```

Its last correction would require

```text
J(g,L+TL^2)=-(2c^2/3)T,                                (28)
```

again impossible by the nonzero scalar form of THM-3716.  This proves the
second word in `(4)` and completes the theorem.

## 5. The lower right shear has an asymmetric Broughton survivor

Now put `N(u)=M_0E_-(u)` and write its rows as `alpha_u,beta_u`:

```text
alpha_u=(4T^2+(2XT-1)u, 2XT-1),
beta_u =(1+2XT+X^2u, X^2).                             (30)
```

Their curls are

```text
curl(alpha_u)=6T+2Xu+(2XT-1)u_T,
curl(beta_u) =X^2u_T.                                  (31)
```

### 5.1 The lower exposed row

Closure of `beta_u+h alpha_u` is exactly

```text
[4T^2+u(2XT-1)]h_T-(2XT-1)h_X
 +[6T+2Xu+(2XT-1)u_T]h+X^2u_T=0.                      (32)
```

There is also a useful determinant form.  If this row is `dQ`, then
`det(alpha_u,dQ)=1`, so

```text
[4T^2+u(2XT-1)]Q_T-(2XT-1)Q_X=1.                      (33)
```

Let `d=deg_T Q>=1`.  The top `T`-degree in `(33)` first gives
`deg_T u<=1`.  Write `u=u_0(X)+T u_1(X)` and let `q_d(X)` be the top
coefficient of `Q`.  The degree-`d+1` equation is

```text
X q_d'=d(2+Xu_1)q_d.                                   (34)
```

Its `X`-valuation gives `q_d=X^(2d)R`; then `(34)` reduces to
`R'=d u_1R`.  Polynomial degree forces `u_1=0`.  Thus every survivor has

```text
u=u(X).                                                 (35)
```

Conversely, every `u(X)` really survives the first gate.  Let

```text
F'(X)=X^2u(X),
Q=X+TX^2+F(X).                                         (36)
```

Then `dQ=beta_u`, so `(32)` has `h=0`.  It is the unique solution.  If
`u!=0` has top term `u_sX^s` and a nonzero homogeneous correction has top
`X` coefficient `c_m(T)`, the coefficient of `X^(s+m+1)` in `(32)` is

```text
2u_s[Tc_m'+c_m],                                       (37)
```

forcing `c_m=C/T`.  If `u=0`, uniqueness is exactly the charge-kernel
calculation `(19)--(21)`.

Because `F` is divisible by `X^3`, write `F=X^2V(X)` and make the
Jacobian-one source shear

```text
S=T+V(X).                                               (38)
```

Then

```text
Q=X+X^2S,
curl(alpha_u)=6S+2[XV'-V].                             (39)
```

For the Broughton Hamiltonian operator in `(X,S)`,

```text
J(X^iS^j,X+X^2S)
 =-jX^iS^(j-1)+(i-2j)X^(i+1)S^j.                      (40)
```

Both output terms have charge `i-j+1`.  The target `6S` has charge `-1`,
whereas the pure-`X` correction `2(XV'-V)` has only nonnegative charges.
The isolated nonterminating THM-3716 chain in charge `-1` therefore cannot
be cancelled by the correction.  The second left shear is impossible.

### 5.2 The upper exposed row

If `alpha_u+h beta_u=dQ`, then `det(dQ,beta_u)=1`, hence

```text
[1+2XT+X^2u]Q_T-X^2Q_X=-1.                            (41)
```

The same top-`T` argument as `(33)--(35)` forces `u=u(X)`.  Put `z=2XT`
and write `Q_tilde(X,z)=Q(X,z/(2X))`.  Equation `(41)` becomes

```text
[2+z+2X^2u] partial_z Q_tilde
       -X partial_X Q_tilde=-X^(-1).                  (42)
```

Take the lowest Laurent exponent `r` in `X`.  The `2X^2u` term cannot
contribute at that layer.  If its coefficient is `F(z)`, the layer equation
is

```text
(2+z)F'-rF=0                         if r<-1,           (43)
(2+z)F'+F=-1                         if r=-1.           (44)
```

If `r>-1`, no term can produce the `X^(-1)` right side.  If `r<-1`, the
highest coefficient in `(43)` is multiplied by
`deg(F)-r`, so `F=0`.  Equation `(44)` has the unique polynomial solution
`F=-1`.  But an honest source monomial with Laurent exponent
`i-j=-1` has `j=i+1>=1`; its contribution is divisible by `z`.  Thus the
`r=-1` layer can never be `-1`.  No upper exposed row exists, completing the
lower-right proof.

## 6. Counterexample-search meaning and boundary

Every matrix in `(4)` or `(4-)` remains outside `E_2`, because it differs
from `M_0` only by elementary factors.  If one were a Jacobian matrix, Wright's
criterion in THM-3652 would turn it into a planar Jacobian counterexample.
The theorem therefore removes a genuine counterexample cell, while retaining
the positive information that an entire row can be integrated exactly.

The decisive triangularity coordinates and survivors are

```text
upper shear: v_X!=0 dies; v=v(T) gives a Broughton orbit,
lower shear: u_T!=0 dies; u=u(X) gives Broughton plus pure-X debt. (45)
```

This does **not** treat two nontrivial right factors, arbitrary
ring-automorphic images of Cohn, or arbitrary non-elementary matrices.  A
successful successor must couple the two source directions strongly enough
to alter the isolated Broughton charge sector; this points to a genuinely
reduced second right move rather than another triangular conjugacy.

Reproduce the exact checks with

```bash
python3 -B 04-computation/jc2_automorphic_cohn_one_right_shear_thm3721.py
python3 -B -O 04-computation/jc2_automorphic_cohn_one_right_shear_thm3721.py
```

The bounded scout covers 224 nontriangular cells with `deg h<=6`; it guards
the formulas, while the degree and charge arguments above are universal.
**QED.**
