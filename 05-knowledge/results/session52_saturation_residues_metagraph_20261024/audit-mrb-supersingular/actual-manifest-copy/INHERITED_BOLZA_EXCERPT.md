# Pinned inherited Bolza quotient and CM proof excerpt

Source: THM-4012, `01-canon/theorems/THM-4012-weighted-leading-face-good-elliptic-factor-observer.md`, Sections 5.1–5.2, equations (35)–(46).
Repository commit: f50a8c2f3a307869ec51650851dbc43eaedfc176.
Complete source SHA-256: 769305f5ddccb41366088d3ffa6440b7f536491849e1cf25e399a627dc8c41da.
The source's status is PROVED + VERIFIED-EXACT + INDEPENDENTLY GEOMETRIC-AUDITED.

This is a verbatim bounded excerpt, not a new theorem or a substitute for the complete source. Its equation (33) is C_B: y^2=x(x^4+1). It records the inherited elliptic quotient and CM algebra only. The surrounding weighted-face/planar-Jacobian hypotheses and conclusions are outside the present MRB transfer. The original notation E is the current packet's E_minus, not E_plus; their -1 twist is proved in RESEARCH.md.

### 5.1 Two explicit elliptic quotients

On `(33)`, consider

```text
iota_1(x,y)=(1/x,y/x^3).                                (35)
```

The invariants

```text
u=x+1/x,                 v=y(x+1)/x^2                   (36)
```

give the degree-two quotient

```text
E:v^2=(u+2)(u^2-2)=u^3+2u^2-2u-4.                      (37)
```

For this integral Weierstrass equation,

```text
c_4=160,                 Delta=512,              j=8000. (38)
```

A conjugate involution is

```text
iota_2(x,y)=(-1/x, i*y/x^3).                            (39)
```

Indeed, conjugate `(35)` by `(x,y)->(-i*x,b*y)` with `b^2=-i`.
For the basis

```text
omega_0=dx/y,                    omega_1=x dx/y
```

of `H^0(C_B,Omega^1)`, one computes

```text
iota_1^*(omega_0)=-omega_1,     iota_1^*(omega_1)=-omega_0,
iota_2^*(omega_0)=-i*omega_1,   iota_2^*(omega_1)= i*omega_0.   (40)
```

The two invariant quotient differentials are

```text
omega_0-omega_1,                 omega_0-i*omega_1,      (41)
```

and are linearly independent. Therefore the product of the two quotient
pushforwards has full-rank differential and is an isogeny

```text
Jac(C_B)  ~  E^2.                                        (42)
```

### 5.2 Self-contained CM mismatch

Shift `x=u+2` in `(37)` to write

```text
E:y^2=x^3-4x^2+2x.                                      (43)
```

Quotienting by `(0,0)` gives the explicit degree-two isogeny

```text
X=x-4+2/x,
Y=y(2-x^2)/x^2,
E':Y^2=X^3+8X^2+8X.                                    (44)
```

The change `x=-X/2,y=tau*Y`, where `tau^2=-1/8`, identifies `E'` with
`E`. If `alpha` is the resulting degree-two endomorphism of `E`, direct
substitution into the duplication law gives

```text
alpha^2=[-2].                                            (45)
```

Since an elliptic endomorphism of degree two cannot be an integer
multiplication, `(45)` and the characteristic-zero classification of
elliptic endomorphism algebras give

```text
End^0(E)=Q(sqrt(-2)).                                    (46)
```

