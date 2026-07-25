---
id: THM-2217
title: "Square-prefix pole alternative and odd-leading-degree terminal wall"
status: >
  PROVED + INDEPENDENTLY AUDITED. In THM-2181's exact-square-prefix Rees
  setup, but without assuming that the first remainder face avoids rad(h),
  polynomiality forces the sharp alternative j_0 rho>=D or h|A^j_0. Thus
  the first rational tooth is not automatically fatal, but its entire
  denominator must be absorbed by the first remainder face. For a polynomial
  nonmonic quartic P=H^2+qz+r with H=Vz^2+pz+s, odd positive v=deg V, and
  genuinely reduced mate degree 4R-2, every survivor has q!=0, R>=2, and
  R deg(q)>=v. Odd degree makes V nonsquare over C(x); combining with
  THM-2214 strengthens the terminal survivor condition to R>=4, equivalently
  reduced degree at least fourteen. This is a necessary wall, not terminal
  descent or a proof of JC(2).
source: codex-2026-07-24-square-prefix-factor-support
depends_on:
  - THM-2127-pre-resonance-ladder-and-affine-root-residue-obstruction
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2181-exact-square-prefix-compression-and-monic-depressed-quartic-closure
related:
  - THM-2071-planar-keller-quadratic-pencil-normal-form
  - THM-2202-uniform-all-degree-quartic-pole-closure
  - THM-2206-quadratic-deck-augmentation-and-hamiltonian-characteristic-incompatibility
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
---

# THM-2217 -- the factor support of the first square-prefix tooth

The radical-avoidance hypothesis in THM-2181 is used only after the first
rational tooth has been isolated. Retaining that tooth instead of declaring
it impossible gives a useful alternative on the remaining nonmonic quartic
branch.

## 1. Factor-sensitive square-prefix alternative

Fix positive integral weights `w=(w_1,w_2)` and put

```text
W=w_1+w_2.                                            (1)
```

Let `K,B in C[x,y]`, set

```text
f=K^2+B,                                              (2)
```

and write `h=in_w(K)`. Assume that `h` is nonconstant, not a proper power,
and `w`-homogeneous of degree `d`. Encode the lower faces by

```text
K(tau)=sum_(eta>=0) K_eta tau^eta,
deg_w(K_eta)=d-eta,                 K_0=h,

B(tau)=sum_(delta>=rho) B_delta tau^delta,
deg_w(B_delta)=2d-delta,            B_rho=A!=0.       (3)
```

Assume

```text
rho>0,                 alpha=deg_w(A)=2d-rho>0.       (4)
```

No divisibility condition on `A` is imposed. Suppose `{f,g}=1`, and reduce
the mate by polynomial target shears until

```text
in_w(g)=c h^s,                c!=0,             s odd. (5)
```

Put

```text
D=(s+2)d-W,                 j_0=(s+1)/2.              (6)
```

As in THM-2181, `D>0` is the Rees defect of the constant bracket. Then

> **Factor-sensitive pole alternative.**
>
> ```text
> j_0 rho>=D                 or                 h divides A^j_0.
>                                                               (7)
> ```

The first branch is equivalently

```text
d+j_0 alpha<=W.                                      (8)
```

Both branches may hold; (7) is a necessary disjunction.

### Proof

Put

```text
F(tau)=K(tau)^2+B(tau),
G(tau)=sum_(t>=0)g_t tau^t,
deg_w(g_t)=sd-t.                                      (9)
```

The Keller equation is

```text
{F(tau),G(tau)}=tau^D.                               (10)
```

For `ell>=0`, take the constant-one branch

```text
S_ell(tau)
 =K(tau)^(s-ell)
  (1+B(tau)/K(tau)^2)^((s-ell)/2)
 =F(tau)^((s-ell)/2).                                (11)
```

THM-2127's centralizer induction, exactly as in THM-2181, gives scalars
`lambda_ell`, with `lambda_0=c`, such that for every `T<D`,

```text
G(tau)=sum_(ell d<=T)
 lambda_ell tau^(ell d)S_ell(tau)
 modulo tau^(T+1).                                   (12)
```

Assume `j_0 rho<D`. The first seed expands as

```text
S_0=sum_(j>=0) binom(s/2,j)K^(s-2j)B^j.              (13)
```

Terms with `j<j_0` are polynomial, and terms with `j>j_0` have larger
defect. At defect `j_0 rho`, the unique nonpolynomial contribution from
this seed is

```text
c binom(s/2,j_0) A^j_0/h.                            (14)
```

Its coefficient is nonzero because `s` is odd. No later seed reaches a
pole at this defect. If `s-ell` is even and nonnegative, `S_ell` is a
polynomial power of `F`. If `s-ell` is odd, its first possible pole has
defect

```text
ell d+((s-ell+1)/2)rho
 =j_0 rho+ell(d-rho/2)
 =j_0 rho+ell alpha/2
 >j_0 rho                                             (15)
```

for `ell>0`. Since `rho<2d`, every seed relevant at this defect has
`ell<=s`, so this parity split is exhaustive.

The coefficient of `G` at defect `j_0 rho` is polynomial. Equation (14)
therefore forces

```text
h divides A^j_0.                                     (16)
```

This proves (7). Substituting `rho=2d-alpha` and `s=2j_0-1` makes its first
branch equivalent to (8). Notice what is not proved: without the
radical-avoidance hypothesis, THM-2181's later deductions
`d+alpha>=W` and `s=1` do not follow.

## 2. Odd-leading-degree nonmonic quartics

Let a planar polynomial Keller pair have an exact square-prefix component

```text
P=H^2+L,

H=V(x)z^2+p(x)z+s(x),
L=q(x)z+r(x),                         V!=0,           (17)
```

and reduce its mate by polynomial target shears. Suppose its remaining
fibre degree is

```text
n=4R-2,                         R>=1.                 (18)
```

Assume

```text
v=deg(V)>=1                       is odd.             (19)
```

Then every such survivor satisfies

```text
q!=0,                 R>=2,                 R deg(q)>=v.
                                                               (20)
```

### Choice of face

Give `x` weight one and `z` an integer weight `M`. Choose `M` sufficiently
large that:

1. the unique top face of `H` is the leading monomial of `Vz^2`;
2. `H^2` lies strictly above every face of `L`;
3. the reduced mate's top fibre term lies strictly above its lower fibre
   terms; and
4. when `q!=0`, the face of `qz` lies strictly above that of `r`.

There are only finitely many strict linear inequalities on `M`, so such a
positive integer exists. With `v_0=lc(V)`,

```text
h=v_0 x^v z^2,                   d=v+2M.              (21)
```

Because `gcd(v,2)=1`, `h` is not a proper power. The leading coefficient
of the reduced degree-`4R-2` Faber term is a nonzero scalar multiple of
`V^(2R-1)`. Hence the reduced mate has top face

```text
c h^(2R-1),                                         (22)
```

so Section 1 applies with

```text
s=2R-1,                         j_0=R.               (23)
```

If `q!=0`, write `m=deg(q)`. The first remainder face is

```text
A=q_0 x^m z,                      alpha=m+M>0.        (24)
```

Here

```text
d>W=1+M.
```

Thus the first branch of (8),

```text
d+R alpha<=W,
```

is impossible. Equations (7) and (24) force

```text
x^v z^2 divides (x^m z)^R.                           (25)
```

Comparison of the `z`- and `x`-exponents gives

```text
R>=2,                         Rm>=v.                 (26)
```

It remains to justify `q!=0`. If `q=0` and `r` is nonconstant, then the
first remainder face is a nonzero monomial in `x` alone. The first branch
is still impossible, while `h|A^R` fails because `h` contains `z`.
If `q=0` and `r` is constant, then

```text
{P,g}=2H{H,g}
```

is divisible by the nonconstant polynomial `H` and cannot be a nonzero
constant. These contradictions prove `q!=0`, and hence (20).

## 3. Combination with the spectral-curve closure

Odd `deg(V)` makes `V` nonsquare in `C(x)`: its valuation at the place at
infinity is odd. THM-2214 closes nonsplit terminal reduced degrees

```text
2,6,10.
```

These correspond to `R=1,2,3`. Combining that theorem with (20), every
noninvertible odd-leading-degree terminal survivor must satisfy

```text
R>=4,
n=4R-2>=14,
q!=0,
R deg(q)>=deg(V).                                    (27)
```

This is the first factor-support wall on the terminal nonmonic branch. It
does not classify the remaining spectral curves or turn `H` into a
quadratic coordinate.

## 4. Boundary and next tooth

The oddness in (19) is load-bearing for this proof. If `deg(V)` is even,
the monomial `x^v z^2` is a proper square and THM-2127's power-free
centralizer coordinate must first be changed; Section 1 cannot simply be
applied with that `h`.

The exact-square-prefix hypothesis is also load-bearing. A merely leading
square does not identify the full seed (11), and treating lower faces of
`H^2` as independent remainder faces recreates the error repaired by
THM-2181.

When `h|A^R`, the next term of the first seed has schematic denominator

```text
A^(R+1)/h^3.
```

It is tempting to iterate (7). Later seeds and lower faces can enter at
subsequent defects, so no higher-tooth divisibility cascade is asserted
without a new separation theorem. The present result is exactly the first
tooth alternative and its rigorously forced support inequality. It does
not close terminal nonmonic quartics, `JC(2)`, or `DC(2)`. QED.
