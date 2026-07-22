---
id: THM-2134
title: "Preterminal factor edges are coarsened powers or are terminally short"
status: >
  PROVED. For a reduced proper-power planar Keller pair, fix a factor of the
  primitive top root and take the globally exposed lower Newton edge from the
  top vertex in the corresponding factor valuation. If its slope is strictly
  beyond the radial slope, then either its residue-field edge polynomial is
  the exact m/gcd(m,n)-th power predicted by a coarsened approximate root, or
  the edge ends within an explicit terminal window. The floor-exact bound is
  D <= L+r0*floor(en/(tau*r0)). The proof combines THM-2127's full formal
  train with a sharp zero-block recurrence for fractional powers. This is a
  factorwise residue-field statement: it does not make the roots for different
  factors compatible or lift them to a global polynomial approximate root.
source: codex-2026-07-22-JC2-preterminal-factor-edge
depends_on:
  - THM-2127
related:
  - THM-2102
  - THM-2132
---

# THM-2134 -- the preterminal factor-edge power dichotomy

Let `w=(w_1,w_2)` be a positive integral weight, put

```text
W=w_1+w_2,
```

and let `h in C[x,y]` be nonconstant, power-free, and `w`-homogeneous of
degree `d`. Suppose `f,g in C[x,y]` and

```text
{f,g}=1,
in_w(f)=h^m,                    m>=2,
in_w(g)=c h^n,                  c!=0,                 (1)
```

where the mate has been reduced so that

```text
n=q m+r,              q>=0,       0<r<m,
D=(m+n)d-W>0.                                         (2)
```

Write the full weighted-face decompositions as Rees polynomials

```text
F(z)=sum_(delta>=0) f_delta z^delta,       f_0=h^m,
G(z)=sum_(t>=0) g_t z^t,                  g_0=c h^n.  (3)
```

Thus `deg_w(f_delta)=md-delta`, with a missing face interpreted as zero,
and

```text
{F(z),G(z)}=z^D.                                      (4)
```

## 1. The globally exposed factor edge

Fix a homogeneous irreducible factor `pi` of `h`, and set

```text
e=v_pi(h).
```

Use `v_pi(0)=+infinity` and define the steepest slope from the top factor
vertex by

```text
tau=max_(delta>0, f_delta!=0)
        [me-v_pi(f_delta)]/delta.                     (5)
```

Assume

```text
tau d>e.                                              (6)
```

This is exactly the nonradial case. In particular, every factor selected by
THM-2132 satisfies (6): its first-face slope `mu` obeys `mu d>e`, and the
global maximum `tau` is at least `mu`.

Let

```text
E={delta>=0 : f_delta!=0 and
                  v_pi(f_delta)+tau delta=me}.        (7)
```

The set contains zero and at least one positive defect. Put

```text
r0=gcd(E minus {0}),          L=max E,          s=L/r0. (8)
```

Every `tau delta` with `delta in E` is an integer. Bezout applied to the
positive elements of `E` therefore shows that

```text
tau r0 is a positive integer.                         (9)
```

Localize `C[x,y]` at the height-one prime `(pi)` and let

```text
kappa=Frac(C[x,y]/(pi))
```

be its residue field. For `delta in E`, reduce

```text
pi^[-(me-tau delta)] f_delta
```

modulo `pi`, and call the resulting nonzero element `p_delta in kappa`.
The factor-edge polynomial is

```text
P(T)=sum_(delta in E) p_delta T^(delta/r0)
       in kappa[T].                                  (10)
```

It has degree `s`, and, if `u=h/pi^e` in the local ring, then

```text
P(0)=bar(u)^m!=0.                                    (11)
```

Geometrically, (10) is obtained from the edge initial form by the integral
substitution

```text
T=z^r0/pi^(tau r0).
```

Put

```text
M=m/gcd(m,n).                                        (12)
```

Then at least one of the following conclusions is forced:

```text
P(T)=Q(T)^M                 for some Q in kappa[T],   (13)
```

or

```text
D<=L+r0 floor(en/(tau r0)).                           (14)
```

In the power branch, `Q` can be normalized by

```text
Q(0)=bar(u)^gcd(m,n),                                 (15)
```

after multiplying it by an `M`th root of unity. In particular,

```text
M divides L/r0.                                      (16)
```

Thus a noncentral arbitrary tail has only two possibilities at its globally
first factor edge: a factorwise coarsened approximate root, or an explicitly
short preterminal edge.

## 2. Why polynomiality creates a zero block

Give a coefficient of `z^t` the auxiliary value

```text
nu_tau(a z^t)=v_pi(a)+tau t.                         (17)
```

The definition of `tau` says that every term of `F` has value at least `me`,
and equality holds exactly on the edge `E`.

THM-2127 gives, through every defect `T_0<D`, the full formal-train identity

```text
G(z)=sum_(ell d<=T_0) lambda_ell z^(ell d)
      h^(n-ell)[F(z)/h^m]^((n-ell)/m)
      modulo z^(T_0+1),                              (18)
```

with `lambda_0=c`. The `ell`th train in (18) has auxiliary value at least

```text
en+ell(tau d-e).                                     (19)
```

Consequently, at a fixed defect `t`, every `ell>0` train has factor valuation
strictly greater than

```text
en-tau t,                                            (20)
```

by (6). Inside the seed-zero train, every selection involving a face off
`E` also has strictly greater valuation. The sum of all edge selections at
defect `t=k r0` is the coefficient of the distinguished formal branch

```text
Y(T)=P(T)^(n/m) in kappa[[T]],
Y(0)=bar(u)^n.                                       (21)
```

Here (21) means the unique binomial branch with the displayed constant term.
It satisfies `Y^m=P^n`.

Let `y_k=[T^k]Y`. If

```text
k r0<D                 and                 en-tau k r0<0, (22)
```

then

```text
y_k=0.                                                (23)
```

Indeed, if `y_k` were nonzero, (20) would be the unique least
`pi`-valuation in the exact formula for the polynomial face `g_(k r0)`.
That least valuation is negative by (22), while a polynomial has
nonnegative `pi`-valuation. No higher seed or off-edge selection can cancel
it. This contradiction proves (23).

The first integer `k` satisfying the second inequality in (22) is

```text
k_0=floor(en/(tau r0))+1,                             (24)
```

and the last one satisfying the first is

```text
k_1=floor((D-1)/r0).                                 (25)
```

Thus polynomiality supplies the exact consecutive zero block

```text
y_k=0        for k_0<=k<=k_1.                        (26)
```

There is no ceiling convention hidden here: defects and `D` are integers,
so `k r0<D` is precisely `k r0<=D-1`.

## 3. A sharp zero-block lemma for fractional powers

Write

```text
P(T)=sum_(i=0)^s p_i T^i,       p_0 p_s!=0,
Y(T)=sum_(k>=0) y_k T^k.                              (27)
```

Differentiating `Y^m=P^n` and using that `P,Y` are units in the formal power
series ring gives

```text
m P Y'=n P'Y.                                        (28)
```

For `k>=1`, the coefficient of `T^(k-1)` is

```text
m k p_0 y_k
 +sum_(i=1)^min(s,k) [m(k-i)-n i]p_i y_(k-i)=0.      (29)
```

Because the field has characteristic zero and `p_0!=0`, equation (29)
determines `y_k` from the preceding `s` coefficients. Therefore any block of
`s` consecutive zero coefficients propagates forever. In that event `Y` is
a polynomial.

If (14) fails, put

```text
A_0=floor(en/(tau r0)).
```

Since `L=s r0`, the strict reverse inequality is

```text
D>r0(A_0+s).
```

Both sides are integers, so `D-1>=r0(A_0+s)`. Hence

```text
k_1=floor((D-1)/r0)>=A_0+s,
```

and (26) contains the `s` consecutive indices

```text
A_0+1,...,A_0+s.                                     (30)
```

The recurrence proves `Y in kappa[T]`.

Now set

```text
g_0=gcd(m,n),             m=g_0 M,             n=g_0 N,
gcd(M,N)=1.                                            (31)
```

Unique factorization in `kappa[T]`, applied to `Y^m=P^n`, shows for every
irreducible `q` that

```text
M v_q(Y)=N v_q(P).
```

Thus `M` divides every irreducible multiplicity of `P`. We have
`P=a Q_0^M` for a unit `a in kappa^*`. But (11) makes the constant term of
`P` an `M`th power, since `m=g_0 M`; as `Q_0(0)!=0`, the unit `a` is also an
`M`th power. Absorbing it into `Q_0` proves (13), and normalization gives
(15). This completes the proof.

The recurrence length is optimal as a formal-series statement. For `s>=2`
and a nonintegral positive rational `alpha`, the nonpower

```text
P(T)=1+T^s                                             (32)
```

has

```text
P(T)^alpha=sum_(j>=0) binom(alpha,j)T^(sj),            (33)
```

so it has arbitrarily many zero blocks of length exactly `s-1`, but none of
length `s`. Therefore the proof cannot replace the `s` consecutive positions
in (30) by `s-1` without additional Keller information.

## 4. The terminally short alternative and two tame controls

The endpoint belongs to the edge, so

```text
0<=v_pi(f_L)=me-tau L,
L<=me/tau.                                            (34)
```

Consequently the short alternative (14) implies the useful slope invoice

```text
tau D<=e(m+n),                                        (35)
```

or, after substituting `D=(m+n)d-W`,

```text
(m+n)(tau d-e)<=tau W.                                (36)
```

Thus a nonpower edge is not merely finite: its slope excess above the radial
line is quantitatively trapped by the terminal bracket weight.

The short bound is attained by a genuine Keller pair. With ordinary weights,

```text
f=y^4+y^2+x,                   g=y.                   (37)
```

Here

```text
h=y,       m=4,        n=1,        D=3,
pi=y,      e=1,        tau=4/3,
E={0,3},   r0=L=3,     s=1.                            (38)
```

The defect-two face `y^2` lies strictly above the globally steepest edge,
while the terminal face `x` is its other endpoint. Thus

```text
P(T)=1+xT in C(x)[T],              M=4,               (39)
```

which is not a fourth power, and (14) is the equality `3<=3`. This also
shows why a central **first** lower face need not make the globally exposed
factor edge central: the terminal tail can select a steeper edge.

The standard infinite tame family realizes the power branch exactly. Let

```text
f=x+(x^2+y)^k,             g=x^2+y,             k>=2. (40)
```

Then, in the notation above,

```text
h=x,       m=2k,       n=2,       D=2k,
pi=x,      e=1,        tau=2,
E={0,1,...,k},         r0=1,       L=s=k.             (41)
```

Indeed,

```text
F(z)=(x^2+yz)^k+x z^(2k-1),
G(z)=x^2+yz.                                          (42)
```

The extra `x z^(2k-1)` term lies strictly above the factor edge. In the
residue field `kappa=C(y)`,

```text
P(T)=(1+yT)^k,
M=(2k)/gcd(2k,2)=k.                                  (43)
```

Moreover the short bound would read `2k<=k+1`, which fails for every
`k>=2`; the theorem therefore selects exactly the coarsened-power edge that
actually occurs.

## 5. Exact scope and remaining compatibility problem

The theorem uses the **globally steepest** edge for one fixed factor `pi`.
This packages every arbitrary later face into one lower-Newton-polygon object,
rather than iterating the first-tooth argument face by face. It strengthens
THM-2132 as follows:

```text
noncentral factor
  -> globally exposed pi-edge
  -> factorwise coarsened power or terminally short edge.                (44)
```

The first conclusion is deliberately local. The polynomial `Q` lives in the
residue field of the divisor `pi=0`. Different factors of `h` can have
different slopes, edge lengths, residue fields, and roots. Equation (13) does
not by itself:

1. make those factorwise roots compatible;
2. lift `Q` from `kappa[T]` to `C[x,y,z]`;
3. produce a common polynomial approximate root for both `F` and `G`; or
4. exclude the terminally short inequality (14).

The lost coordinate is precisely cross-factor root compatibility. A future
global descent needs either a simultaneous Hensel/Chinese-remainder sidecar
for these edge roots or a terminal bracket obstruction to (14). Central
slopes `tau d=e` also lie outside the theorem. Thus THM-2134 narrows the
arbitrary-tail factor-initial branch but does not prove planar JC.
