# A rational-primitive dichotomy: multi-root coefficients force a torsion-free unit response

**Status:** PROVED + EXACT HOSTILE-AUDITED; canonized in
[THM-3326](../01-canon/theorems/THM-3326-linear-in-z-unit-response-trichotomy-and-jet-torsion.md).

## Statement

Let `K` be a characteristic-zero field, let `f,g in K[x]` with `g` nonzero,
and put

```text
R=K[x,z],                 P=f(x)+g(x)*z,
D_P=P_x*partial_z-P_z*partial_x,
C_P=R/D_P(R),             theta=[1] in C_P.
```

If `g` is nonconstant and is **not** of the form

```text
g=c*(x-a)^d,              c in K*, a in K, d>=2,        (1)
```

then

```text
Ann_K[P](theta)=(0).                                     (2)
```

Thus no nonzero polynomial response `F(P)` is a Hamiltonian derivative.  In
particular there is no polynomial mate `Q` with `Jac(P,Q)=1`.

The rational-integration input is field-stable and exact:

> For every characteristic-zero extension `L/K`, the rational function
> `1/g` has a primitive in `L(x)` if and only if `(1)` holds.  In that case a
> primitive is `-1/(c*(d-1)*(x-a)^(d-1))`.

The pure-power case is necessary here, not automatically sufficient for a
polynomial response.  When additionally `f'(a)!=0`, translating `x-a` and
`P-f(a)` puts it under the simple-zero theorem in
`jacobian_deformed_torsion_agent_theorem_ready.md`, which gives

```text
Ann_K[P](theta)=((P-f(a))^(d-1)).                        (3)
```

### Complete classification under a unimodular gradient

For this entire linear-in-`z` family, the preceding results give an exact
trichotomy whenever `(P_x,P_z)=R`:

```text
g in K*:                         theta=0;
g=c*(x-a)^d, d>=2:               Ann(theta)=((P-f(a))^(d-1));
g has at least two geometric roots: Ann(theta)=(0).      (T)
```

Indeed, over an algebraic closure, a root `a_i` of `g` would give a common
gradient zero unless

```text
g'(a_i)=0 and f'(a_i)!=0.                                (G)
```

If `g'(a_i)!=0`, choose `z=-f'(a_i)/g'(a_i)`; if `g'(a_i)=0=f'(a_i)`, every
`z` is bad.  Conversely `(G)` at every root leaves no common zero, hence the
gradient ideal is the unit ideal (and this descends faithfully to `K`).  Thus
unimodularity forces every root of a nonconstant `g` to be repeated.  One
geometric root is precisely the pure-power lane covered by `(3)`; two or more
roots are covered by `(2)`.  The constant lane has
`D_P(-x/g)=1`.  This proves `(T)`.

## Proof

Work in `K(P)(x)`, where

```text
z=(P-f(x))/g(x),                D_P=-g(x)*partial_x      (4)
```

with `P` held fixed.  Suppose a nonzero `F(T) in K[T]` and `Q in R` satisfy

```text
D_P(Q)=F(P).
```

Since `P` is transcendental over `K`, `F(P)` is nonzero.  Equation `(4)`
then gives

```text
d/dx(-Q/F(P))=1/g(x) in K(P)(x).                         (5)
```

It remains to classify when `1/g` is rationally exact.  Extend constants to
an algebraic closure and factor

```text
g=c*product_i (x-a_i)^(m_i),
```

where there are `k` distinct roots and `d=sum_i m_i=deg g`.

1. A derivative of a rational function has zero residue at every finite
   pole.  Hence `g` cannot have a simple root; every `m_i>=2`.
2. If `H'=1/g`, subtract the constant value of `H` at infinity.  At each
   `a_i`, `H` has exact pole order `m_i-1`.  Its reduced denominator therefore
   has degree

   ```text
   sum_i(m_i-1)=d-k.                                    (6)
   ```

3. At infinity, integrating `1/g=c^(-1)x^(-d)+O(x^(-d-1))` shows that `H`
   has a zero of exact order `d-1`.  A nonzero proper rational function whose
   denominator has degree `d-k` can have order at infinity at most `d-k`.
   Hence `d-1<=d-k`, so `k<=1`.

Therefore `k=1` and `g=c(x-a)^d` with `d>=2`.  This conclusion descends
without a specialization argument: `c` is the leading coefficient in `K`,
and the coefficient of `x^(d-1)` is `-c*d*a`, so `a` lies in `K` because the
characteristic is zero.  The displayed primitive proves the converse over
`K` itself.  Applying this field-stable classification to `(5)` proves `(2)`.

## Smooth-gradient two-root family

The obstruction is not merely caused by a singular gradient.  Take

```text
f=x,
g=lambda*x^r*(1+a*x)^s,
lambda,a in K*,                 r,s>=2.                  (7)
```

Then `g` has two distinct roots, each of multiplicity at least two, so `(2)`
applies.  Nevertheless the gradient is unimodular.  Since `g` divides
`(g')^2`, the explicit polynomial row

```text
A=1-g'*z,
B=((g')^2/g)*z^2                                         (8)
```

satisfies

```text
A*P_x+B*P_z=(1-g'z)(1+g'z)+(g')^2*z^2=1.                (9)
```

Here

```text
(g')^2/g
 =lambda*x^(r-2)*(1+a*x)^(s-2)*(r+a*(r+s)*x)^2,
```

so polynomiality in `(8)` is explicit.  The canonical divergence
representative attached to this row is

```text
m=A_x+B_z=z*(2*(g')^2/g-g'').                           (10)
```

For the smallest hostile `lambda=a=1`, `r=s=2`,

```text
g=x^2*(1+x)^2,
A=1-2*x*(1+x)*(1+2*x)*z,
B=4*(1+2*x)^2*z^2,
m=(6+20*x+20*x^2)*z.                                   (11)
```

The rational obstruction can be seen directly:

```text
1/g = 1/x^2 - 2/x + 1/(x+1)^2 + 2/(x+1).               (12)
```

The two nonzero simple-pole residues force logarithms in every primitive.
Thus `(11)` is a smooth, gradient-unimodular first coordinate whose canonical
unit response is not merely nonzero: its cyclic `K[P]`-submodule is a free
copy of `K[P]`.

## Boundary and connection to the torsion ladder

- If `g` is constant, `theta=0` because `D_P(-x/g)=1`; this is outside the
  nonconstant dichotomy.
- If `g=c(x-a)` has degree one, `1/g` needs a logarithm, so `(2)` still holds;
  the exceptional rationally exact powers start at degree two.
- The pure power `g=lambda*x^r` has one deleted divisor and yields the finite
  ladder of THM-3318 and its simple-zero deformation.  Adding a second
  repeated root preserves gradient smoothness but destroys every possible
  polynomial-in-`P` killing response.  The root-support count, not just the
  pole multiplicities or gradient ideal, is therefore the missing sidecar.
