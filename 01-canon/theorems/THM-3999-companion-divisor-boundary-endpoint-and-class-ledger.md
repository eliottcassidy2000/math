---
id: THM-3999
title: "The reduced 2:3 companion has an exact boundary endpoint and class ledger"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. On THM-3973's
  smooth completion X2, THM-3992's mixed residual factors as G=tQ with Q in
  B2. The boundary order is exactly ord_D(Q)=2, while the total strict
  companion meets D in Spec k[y]/(gamma-R(0,y)). It is disjoint from D iff
  R(0,y)=0 iff R lies in (p^2,py). The divisor identities are
  div(t)=L-2D, div(Q)=sum C_i+2D, and div(G)=L+sum C_i, hence
  [L]=2[D] and [sum C_i]=-2[D]. These are total-divisor statements: they do
  not prove Q irreducible, complete the node-address census, or turn the
  two-clutch graph critical group into a local class-group obstruction. On
  the THM-3997 live seam the normalized endpoint is the residual-mu5
  invariant E(y)=1-(R/gamma)(0,y), with E'(0)=-b. Thus b!=0 forces an
  endpoint, while boundary disjointness forces b=0 and R in (p^2,py) but
  retains the mandatory interior coefficient
  [p^2](R/gamma)=-16/(3A5^2)!=0.
source: root + conductor_incidence / planar Jacobian continuation, 2026-08-24
audit: >
  PASS (conductor_incidence + laurent_rows, 2026-08-24). Q membership,
  boundary saturation, the scheme-theoretic endpoint iff, reduced divisor
  multiplicities, class signs, the R=0 G_m hostile, irreducibility caveats,
  and the THM-3994/3996 firewall were independently reconstructed. A repeated
  endpoint gives the expected length-two local algebra. The provisional
  same-direction incidence matrix was corrected to the oppositely oriented
  two-cycle. The live-seam normalization, residual-mu5 invariance, endpoint
  derivative, and boundary-blind mandatory p^2 coefficient were replayed
  against THM-3997. Normal, optimized, and stored certificate outputs match.
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
related:
  - THM-3955-node-cotangent-normalization-kernel-and-conductor-torsion
  - THM-3994-double-resultant-collision-separates-two-address-and-length-two-seams
  - THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy
script: 04-computation/jc2_companion_divisor_endpoint_thm3999.py
output: 05-knowledge/results/jc2_companion_divisor_endpoint_thm3999.out
script_sha256: 835d0b663e8d8cc3251b30b109d039039d681651dac034174180f4d024c63214
output_sha256: fefdf2f5f29f7eddb01c96ae8dd3357e15f969174c5421d20b162db78fd6cdf3
hash_basis: raw LF bytes
---

# THM-3999 -- the companion divisor has a boundary ledger

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. On THM-3973's
completion put

```text
z=1+x^2t,                 p=zt,                 y=xtp,
B_2=k[x,z,p,y],           X_2=Spec(B_2),
D=V(x,z)=A1_y,            Cl(X_2)=Z[D].                  (1)
```

The open `D(x) union D(z)` is the Keller source `U=A2_(x,t)`. Begin in
THM-3992's normalized reduced `(2,3)` cell:

```text
G=C^2-A^3+(3a^2/4)A+a^3/4
 =gamma*u+alpha*p+R(p,y),

u=z-1=x^2t,               alpha=3a/(2gamma),
R in (p^2,y),             a*gamma!=0.                   (2)
```

Choose `H,K in k[p,y]` with

```text
R=p^2H+yK,                R_0(y)=R(0,y)=yK(0,y).         (3)
```

The theorem records exactly what `(2)` determines about the companion
divisor, and exactly what it forgets.

## 1. The companion quotient extends and has boundary order two

In `k(x,t)`, equation `(2)` factors as

```text
G=tQ,

Q=gamma*x^2+alpha*z+z*p*H(p,y)+x*p*K(p,y).              (4)
```

Every term on the right lies in `B_2`, so `Q` extends regularly across
`D`. The equality follows from

```text
u/t=x^2,             p/t=z,             p^2/t=zp,
y/t=xp.                                                    (5)
```

Near `D`, the chart relations from THM-3973 are

```text
z(z-1)^2=x^3y,                    p(z-1)=xy.              (6)
```

Here `x` is a uniformizer and `z-1` is a unit. Because `u=x^2t`,

```text
Q/x^2=G/u

 =gamma+alpha*xy/(z-1)^2+xypH/(z-1)^2+yK/(z-1).          (7)
```

This is regular near `D`, and its restriction is

```text
(Q/x^2)|D=gamma-yK(0,y)=gamma-R_0(y).                    (8)
```

The right side is not the zero polynomial: at `y=0` it equals
`gamma!=0`. Therefore

```text
ord_D(Q)=2,                         ord_D(G)=0.            (9)
```

The second identity also follows from `ord_D(t)=-2`. This is a valuation
repair: division by the source equation `t` produces a regular function,
but leaves a fixed order-two boundary component in its divisor.

## 2. The exact endpoint scheme

Let `C_str` denote the total strict closure in `X_2` of the divisor
`V(Q) subset U`; no irreducibility is assumed. Since `Q/x^2` is its local
equation after removing the fixed boundary multiplicity, equation `(8)`
gives the scheme-theoretic intersection

```text
C_str intersect D = Spec k[y]/(gamma-R_0(y)).             (10)
```

Because `k` is algebraically closed and `R_0(0)=0`, the endpoint scheme
is empty exactly when `R_0` is zero. The monomial-ideal identity

```text
(p) intersect (p^2,y)=(p^2,py)                           (11)
```

then gives the complete iff

```text
C_str intersect D empty
 iff R_0=0
 iff R in (p^2,py).                                      (12)
```

If `R_0!=0`, the endpoint scheme has total length `deg R_0`, counted
with multiplicity, and no endpoint is supported at `y=0`. A simple root
of `gamma-R_0` is a transverse endpoint. A repeated root records only
intersection multiplicity; it does not prove that a component is smooth or
that `Q` is irreducible.

### 2.1 The live seam has an invariant endpoint gate and a blind interior coefficient

Now impose the live continued seam of THM-3997. Put

```text
A5=a^5,               Rtilde=R/gamma,
b=[y]Rtilde=beta/gamma.                                  (12a)
```

The residual fifth-root action scales `R` and `gamma` by the same character,
fixes `p,y`, and fixes `A5`. Since `gamma` is a unit, multiplying the endpoint
equation by `gamma^-1` does not change its closed subscheme. Thus `(10)` has
the invariant normalized form

```text
C_str intersect D=Spec k[y]/(E(y)),
E(y)=1-Rtilde(0,y),          E(0)=1,       E'(0)=-b.      (12b)
```

In particular,

```text
b!=0  implies  C_str intersect D is nonempty.             (12c)
```

Indeed, `b!=0` makes `Rtilde(0,y)` nonconstant, and a nonconstant polynomial
over the algebraically closed field `k` takes the value one. The converse is
false: if `b=0`, higher pure-`y` terms may still produce endpoints.

The disjoint lane is sharper. Equations `(11)--(12)` force

```text
C_str intersect D empty
 implies R in (p^2,py),              b=0.                 (12d)
```

But THM-3997 simultaneously fixes

```text
[p^2]Rtilde=-16/(3A5^2)!=0.                              (12e)
```

Hence on this lane one can write

```text
Rtilde=p^2 Htilde+py Ktilde,
Htilde(0,0)=-16/(3A5^2).                                 (12f)
```

The strict companion can therefore be boundary-disjoint while its residual
is provably nonzero. The endpoint map sees `b` and the remaining pure-`y`
coefficients but kills the first mandatory live coefficient `(12e)`. This is
an exact information-loss statement, not a factorization, owner, node-address,
properness, or Jelonek conclusion.

## 3. Principal divisors and class signs

Let `L` be the strict closure of `V(t)` and let `C_i` be the prime
components of `C_str`. The target nodal cubic in THM-3992 is reduced.
An etale pullback of a reduced scheme is reduced, so `V(G) subset U` is
reduced. Moreover,

```text
Q(x,0)=gamma*x^2+alpha
```

is not identically zero, hence `t` and `Q` have no common component.
Thus every `L,C_i` occurs with multiplicity one. Together with `(9)`
and the fact that `D` is the only boundary prime, this proves

```text
div_X2(t)=L-2D,
div_X2(Q)=sum_i C_i+2D,
div_X2(G)=L+sum_i C_i.                                  (13)
```

Since `Cl(X_2)=Z[D]`, the principal-divisor relations give the exact signs

```text
[L]=2[D],                    [sum_i C_i]=-2[D].           (14)
```

Only if `Q` is independently proved irreducible may the second row be
sharpened to `[C]=-2[D]`. The endpoint polynomial alone does not supply
that proof.

## 4. The zero-residual hostile separates valuation from incidence

Although THM-3997 proves that a Keller pair cannot have `R=0`, that
specialization is a sharp ambient-completion control. In this case

```text
Q_0=gamma*x^2+alpha*z
   =alpha+x^2(gamma+alpha*t).                            (15)
```

It is primitive and linear in `t`, hence irreducible in `k[x,t]`. Its
coordinate ring is

```text
k[x,t]/(Q_0)=k[x,x^-1],
t=-gamma/alpha-x^-2.                                    (16)
```

Thus `C_0=V(Q_0)` is `G_m`. It meets `L` transversely at the two roots

```text
x^2=-alpha/gamma=-3a/(2gamma^2),                         (17)
```

but `C_0` is disjoint from `D`, while

```text
div_X2(Q_0)=C_0+2D.                                     (18)
```

Hence positive boundary order of a regular equation does not force its
strict component to meet the boundary. The fixed boundary component and the
strict endpoint scheme are different coordinates.

## 5. Node-address and class-group firewall

THM-3996 supplies the missing ownership sidecar only after the node-address
graph is complete and the appropriate properness hypothesis is imposed. If
the two known clutches form its full proper connected packet, their directed
edges have opposite orientations. With vertices `W_L,W_C`, the incidence
matrix is

```text
B=[[-1, 1],
   [ 1,-1]],                    BB^T=[[ 2,-2],[-2, 2]].  (19)
```

Its kernel is `Z(1,1)`, and the underlying two-edge multigraph has critical
group `Z/2`. This statement is conditional on the complete-packet
hypothesis. The two roots in `(17)`, or in `Q(x,0)` generally, are not a
census of all source addresses over the target node. In fact THM-3996's
Keller degree gap proves that this two-cycle cannot exhaust a proper finite
fibre: there is an additional disjoint address packet, or the node is in the
nonproperness locus.

The graph `Z/2` in `(19)` is unrelated to THM-3994's local class group

```text
Cl(k[L,V,Z]_(L,V,Z)/(V^2-LZ))=Z/2.                       (20)
```

That group belongs to one normal `A_1` surface singularity created by a
length-two Rees center. The two THM-3992 clutch points lie in the smooth
surface `X_2`, so their local class groups vanish; globally
`Cl(X_2)=Z[D]` is torsion-free. The coincidence of the two abstract
`Z/2` groups carries no canonical map and no Keller obstruction.

## 6. Preserved data, destroyed data, and next tests

The map from the mixed residual to the completion sidecar is

```text
R(p,y)
  |-> Q=G/t
  |-> (ord_D Q, endpoint polynomial, total divisor class)

  = (2, gamma-R(0,y), -2[D]).                            (21)
```

It preserves boundary valuation, endpoint length, and the total strict
class. It destroys factor ownership, normalized-address orientation,
completeness of the node fibre, and properness/Jelonek information.
THM-3996 restores the address and properness coordinates only after a
complete packet is specified.

On the THM-3997 live seam, division by `gamma` refines the endpoint coordinate
to the invariant polynomial `E(y)` in `(12b)`. Its first derivative records
`-b`, but the projection `p=0` erases the forced nonzero coefficient `(12e)`.

The cheapest next tests are therefore precise:

1. factor `Q` or control the normalization owners of its two known germs;
2. decide whether `R(0,y)` vanishes, is squarefree, or has repeated roots;
3. census all source addresses over the forced target node;
4. locate that node relative to the nonproperness locus.

None of `(4)--(21)` proves `Q` irreducible, that the two known clutches
exhaust the fibre, that the target node is proper, that the companion cycle
has an obstructing class, the full reduced `(2,3)` cell, or `JC(2)`.
**QED.**

## Reproduction

```bash
python3 04-computation/jc2_companion_divisor_endpoint_thm3999.py
python3 -O 04-computation/jc2_companion_divisor_endpoint_thm3999.py
sha256sum 04-computation/jc2_companion_divisor_endpoint_thm3999.py \
  05-knowledge/results/jc2_companion_divisor_endpoint_thm3999.out
python3 agents/check_docs.py
```
