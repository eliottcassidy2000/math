# Polynomial graphs close for every linear target plane

**Status: PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; unnumbered pending a
fresh identifier scan and independent hostile audit.**

THM-2705 classifies every linear target plane containing `d` on a polynomial
graph for the THM-2696 reflection quotient.  Its stated next boundary was a
target plane not containing `d`.  That boundary also closes, and the two
pieces together classify **every** affine-linear rank-two target projection on
every global graph `z=f(x,y)`.  The new mechanism is a projective-normal split:
one normal coordinate gives an invertible Euler operator, while the other
forces a rank-one zero-Jacobian boundary.

This is a proof candidate about the fixed THM-2696 quotient family.  It is not
a proof of `JC(2)`, `DC(2)`, or a statement about arbitrary polynomial source
surfaces or nonlinear target projections.

## 1. Projective target normals

On the graph `z=f(x,y)`, write

```text
A=x^2-2y,          B=y^2-2xf,          d=f.              (1)
```

A rank-two linear target projection is a two-plane in the row space with a
projective normal

```text
ell=(lambda,mu,nu) in P^2.                               (2)
```

If `nu=0`, the row plane contains `d`; this is exactly THM-2705.  Suppose
`nu!=0` and normalize

```text
ell=(a,b,1).                                             (3)
```

An explicit row basis for its orthogonal plane is

```text
U=A-a d,             V=B-b d,                            (4)
```

because `(1,0,-a) cross (0,1,-b)=(a,b,1)`.  Any other row
basis changes the planar Jacobian only by a nonzero scalar.  Thus `(3)--(4)`
lose no target plane or Keller map.

Direct differentiation gives

```text
Jac(U,V)=-2 P_f,                                         (5)

P_f=(a y+b+2x)f_x+(a f+b x+2x^2)f_y+2f-2xy.            (6)
```

Put

```text
P_f=c,                    c=-kappa/2.                    (7)
```

The Keller case is exactly `c!=0`.

## 2. The `a!=0` half is empty at nonzero Jacobian

Let `r=deg_y(f)` and let `p_r(x)` be the leading coefficient.

- For `r>=3`, the term `a f f_y` has unique top `y`-degree `2r-1`,
  with coefficient `a r p_r^2`.  Every other term in `(6)` has degree at
  most `r+1`, and `2r-1>r+1`.  This is impossible in characteristic zero.
- For `r=2`, the coefficient of `y^3` is

  ```text
  a(p_2'+2p_2^2).                                       (8)
  ```

  A nonzero polynomial cannot satisfy `p_2'+2p_2^2=0`: at positive
  degree the two terms have degrees `deg(p_2)-1` and `2deg(p_2)`, while at
  degree zero the square survives.

Hence `r<=1`, including the `r=0` boundary.  Write

```text
f=p(x)y+q(x).                                           (9)
```

The coefficients of `y^2,y,1` in `P_f-c` are

```text
a p',                                                   (10)
a p^2+a q'+(b+2x)p'-2x+2p,                             (11)
a p q+(b x+2x^2)p+(b+2x)q'+2q-c.                       (12)
```

Since `a!=0`, `(10)` makes `p` constant and `(11)` gives

```text
q'=2x/a-2p/a-p^2.                                      (13)
```

After integrating `(13)`, the `x^2` coefficient of `(12)` is

```text
3(ap+2)/a.                                             (14)
```

Thus `p=-2/a`.  At this value `(13)` becomes `q'=2x/a`, and every
term in `(12)` cancels except `-c`.  Therefore

```text
a!=0 and P_f=c  implies  c=0.                           (15)
```

This proves emptiness for `kappa=-2c!=0`.  It also identifies the exact
equality boundary.  When `c=0`, all survivors are

```text
f=-2(y-x^2/2)/a+q_0.                                   (16)
```

But then `U=A-a f=-a q_0` is constant, so `(U,V)` has rank at most one.
The cancellation in `(15)` is geometric rank collapse, not a hidden Keller
family.

## 3. The `a=0` half has one triangular automorphism

Now the normal is `(0,b,1)` and the target pair is

```text
(U,V)=(A,B-bd).                                        (17)
```

Set

```text
t=y-x^2/2.                                             (18)
```

In `(x,t)` coordinates, equation `(7)` becomes

```text
(2x+b)f_x+2f=2xt+x^3+c.                                (19)
```

After `u=x+b/2`, the left operator is

```text
2(u partial_u+1),                                      (20)
```

which acts on `u^i t^j` by the nonzero scalar `2(i+1)`.  It is therefore
invertible on the whole polynomial ring, not merely through a bounded degree.
The unique solution is

```text
f=((x-b/2)t)/2
  +((x-b/2)(x^2+b^2/4))/8
  +c/2.                                                (21)
```

For `b=0` and `c=-kappa/2`, `(21)` is precisely THM-2702's unique
`(A,B)` graph

```text
f=xt/2+x^3/8-kappa/4
 =xy/2-x^3/8-kappa/4.                                  (22)
```

More importantly, every shifted member is triangular:

```text
U=-2t,
V=-c x+(t+b^2/8)^2-bc/2.                              (23)
```

For `c!=0`, its polynomial inverse is

```text
t=-U/2,
x=((t+b^2/8)^2-bc/2-V)/c,
y=t+x^2/2.                                             (24)
```

The Jacobian is `-2c=kappa` exactly.

## 4. Combined classification and exact boundary

Together with THM-2705, the projective normal `ell` has no remaining cell:

```text
nu=0:       plane contains d; classified by THM-2705;
nu!=0,a=0:  one explicit triangular graph for each (b,kappa);
nu!=0,a!=0: empty when kappa!=0; rank-drop family when kappa=0.  (25)
```

Thus every nonzero-constant-Jacobian linear target projection of the fixed
THM-2696 quotient on a polynomial graph is an explicit polynomial
automorphism.  This strictly completes the target-plane boundary left open by
THM-2705.

The restrictions that remain are substantive: a source surface not globally
represented as `z=f(x,y)`, a nonlinear target projection, an arbitrary `S4`
resolvent, `JC(2)`, and `DC(2)` are not covered.

## 5. Exact companion and hostile controls

Run

```text
python 04-computation/jacobian_s4_polynomial_graph_all_linear_target_planes_referee_20260728.py
python -O 04-computation/jacobian_s4_polynomial_graph_all_linear_target_planes_referee_20260728.py
```

Both modes must byte-match

```text
05-knowledge/results/jacobian_s4_polynomial_graph_all_linear_target_planes_referee_20260728.out.
```

The companion independently checks the projective row normal, `(5)--(6)`,
the `r=2` Riccati coefficient, all three `r<=1` coefficient equations, the
forced slope and final `-c` collapse, the entire `c=0` rank-drop family,
the Euler inverse, `(21)--(24)`, the THM-2702 `b=0` boundary, and a shifted
`b=3,c=5` hostile.  Exact hashes should be frozen only after the fresh
identifier/audit pass.

