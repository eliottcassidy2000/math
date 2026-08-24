---
id: THM-3934
title: "Infinity-component sextic has a unique reduced splitting conic"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For
  the THM-3932 rational one-place sextic
  4X^3Z^3-(Y^3-X^2Z)^2=0, the only reduced projective conic that splits in
  the quadratic double plane is XZ=0, up to scalar. All other conics whose
  pullback to the branch normalization is a square are nonreduced double
  lines. Consequently the genuine order-three resolvent class found in
  THM-3932 has no second representative supported over a plane curve of
  the splitting-curve construction in degree at most two. This does not
  compute the full resolvent class group:
  higher-degree splitting curves and abstract divisor classes remain open.
source: root / post-THM-3932 low-degree global three-class lane, 2026-08-23
depends_on:
  - THM-3932-infinity-component-linear-conic-torus-sextic-fold-classification
related:
  - THM-3879-rational-torus-sextic-c3-packet-one-place-tradeoff
  - THM-2387-degree-eighteen-h4-elliptic-three-isogeny-atlas
  - THM-3935-linear-conic-resolvent-class-group-unique-cubic-character
script: 04-computation/jc2_infinity_component_unique_splitting_conic_thm3934.py
output: 05-knowledge/results/jc2_infinity_component_unique_splitting_conic_thm3934.out
script_sha256: b658f2314f27ff391051d0886598b67c807b3af505cb0dc2bfde6f61c71aca18
output_sha256: c8660794b4b9d51cc1ba08da35109d218a7e79c417bb365b182797507135f6ee
semantic_sha256: 3325fe5ff3ae2e63a5bd1ed65076301e2eb706fd36864f219ff295ec5dedf989
hash_basis: raw LF bytes
---

# THM-3934 -- the first low-degree twist search is rigid

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero. Let

```text
Gamma: F=4X^3Z^3-(Y^3-X^2Z)^2=0                         (1)
```

be the THM-3932 sextic, and let the associated projective quadratic double
plane have function field

```text
k(P2)(sqrt(-F)).                                           (2)
```

A reduced conic `Q=0` is called **splitting** when its inverse image in this
double plane is generically the union of two components (on each normalized
component if `Q` is reducible). Then, up to a nonzero scalar,

```text
Q=XZ                                                        (3)
```

is the unique reduced splitting conic. More precisely, all conic forms whose
pullback to the normalization of `Gamma` is a square are exactly

```text
Q=L^2,                 L any projective line,
Q=lambda XZ,           lambda in k.                         (4)
```

The first family consists of nonreduced double lines. The second has one
reduced member up to scale, namely `(3)`.

Affinely `Z=1`, the component `X=0` is exactly the support that splits into
the two primes `D^+,D^-` of THM-3932. Thus `(4)` proves a sharply scoped
negative result:

```text
the splitting-curve construction produces no second resolvent
three-class from a curve of degree at most two.             (5)
```

It does not imply that the complete resolvent class group has only one
three-torsion direction.

## 1. Splitting forces a square on the branch normalization

THM-3932 gives the projective normalization

```text
nu([S:T])=[H^2:S^2 T H:S^6],          H=T^3-2S^3.         (6)
```

The branch has only two singular supports: a three-branch point
`O=[0:0:1]` and a unibranch point `I=[1:0:0]`. Away from these points, a
splitting conic has even intersection multiplicity with the unique smooth
branch of `Gamma`. At `I`, splitness forces even total contact, and there is
only one normalization branch on which that contact can lie. For a reducible
conic, apply splitness on each normalized line component and then add their
local intersection numbers. The only possible loss of information is
therefore at `O`: even *total* contact on the conic might conceivably hide odd
contact on two of the three branch addresses. It does not.

In the affine chart `Z=1`, use `t=H` as a local parameter at one of the three
addresses `alpha^3=2`. Equations `(6)` give

```text
x=t^2,                     y=alpha t+O(t^2).               (7)
```

Write a conic through `O` as

```text
a x^2+bxy+cx+dy^2+ey.                                    (8)
```

If `e!=0`, its order is one on all three branches, so the total order is
odd. If `e=0`, the order-two coefficient at the `alpha` branch is

```text
c+d alpha^2.                                               (9)
```

The three values `alpha^2` are distinct. Unless `c=d=0`, `(9)` can vanish
on at most one branch. The other two orders are two, so total evenness forces
the exceptional branch order to be even as well. If `c=d=0` but `b!=0`,
all three orders are three and their sum is odd. Finally `b=0` leaves the
even order-four row `a x^2`. Hence even total contact at `O` always upgrades
to even contact on every normalization branch.

It follows that for every reduced splitting conic,

```text
nu^*Q is an even effective divisor on P1.                  (10)
```

The pulled-back section has degree twelve. Since `k` is algebraically
closed, `(10)` is equivalent to

```text
Q o nu=G(S,T)^2                                            (11)
```

for a binary sextic `G`. This is the bridge that makes the finite exact
classification below complete even at the three-branch singularity.

## 2. Exact square-root classification

Set `S=1`, write `s=T/S` and `h=s^3-2`. The affine normalization is

```text
(x,y)=(h^2,sh).                                            (12)
```

The six conic monomials pull back to the linearly independent polynomials

```text
h^4,       s h^3,       h^2,       s^2h^2,       sh,       1.  (13)
```

Their coefficient matrix in degrees zero through twelve has rank six. Thus
restriction of conics to `Gamma` is injective.

Let

```text
G(s)=sum_(i=0)^6 g_i s^i,             G(s)^2=sum_(j=0)^12 c_j s^j.  (14)
```

Membership of `G^2` in the span `(13)` is equivalent to the seven independent
coefficient equations

```text
c_2+c_5=0,
3c_1+6c_4+8c_7=0,
-c_2+4c_8=0,
c_3+4c_6+8c_9=0,
-c_1+16c_10-2c_4=0,
c_11=0,
64c_12-c_3-4c_6=0.                                      (15)
```

Solving `(15)` by the highest nonzero coefficient of `G` gives the complete
projective chart ledger:

| `deg G` | solutions |
|---:|---|
| `6` | `G=A h^2+Bsh+C`, with `A!=0` |
| `5` | none |
| `4` | `G=Bsh+C`, with `B!=0` |
| `3` | `G=A h`, with `A!=0` |
| `2` | none |
| `1` | none |
| `0` | nonzero constants |

For example, after normalizing `g_6=1`, the exact Groebner chart has

```text
g_5=0,       g_3=-4,       g_1=-2g_4,       g_2^2=0,     (16)
```

and hence `g_2=0`. After normalizing `g_3=1`, the only extra chart has

```text
g_0=-2,       g_1^2=g_2=g_4=g_5=g_6=0,                   (17)
```

which is `G=h`. The remaining rows in the table follow from the same seven
quadrics; the exact companion records every chart.

The first, third, and last rows of the table combine as

```text
G=A h^2+Bsh+C=L o nu,                                    (18)
```

the pullback of an arbitrary line `L=AX+BY+CZ`. Therefore `G^2` is the
pullback of `L^2`. Injectivity of `(13)` forces `Q=L^2`. The only extra root
is `G=h`, and

```text
h^2=(XZ) o nu.                                             (19)
```

Injectivity again forces `Q` to be a scalar multiple of `XZ`. This proves
`(4)`.

## 3. The unique reduced conic really splits

The conic `XZ=0` is the reduced union of two distinct lines. On them the
branch polynomial restricts as

```text
F|_(X=0)=-Y^6,                    F|_(Z=0)=-Y^6.           (20)
```

Both are squares up to a scalar in the algebraically closed field. Thus
each normalized line splits in the quadratic double plane, proving the
converse and `(3)`.

On the affine chart, the line `X=0` gives

```text
(Y^3-X^2)^2-4X^3=Y^6,                                    (21)
```

and its two lifts are exactly the primes `D^+,D^-` used in

```text
div(X)=D^++D^-,                 div(Y^3-X^2+w)=3D^+.       (22)
```

Thus the unique reduced degree-two support recovers the already known
Cardano class and cannot supply an independent twist.

## 4. Why the local elliptic directions are still worth tracking

The weighted exceptional curve suggested by THM-3932 is

```text
E: w^2=y^6-4x^3 subset P(2,1,3).                          (23)
```

In the chart `y=1`, it is the `j=0` elliptic curve

```text
w^2=1-4x^3.                                                (24)
```

The natural divisors meet the flexes `(x,w)=(0,+/-1)`. After putting
`u=-x,v=w/2`, the three-division polynomial of
`v^2=u^3+1/4` is a nonzero scalar times

```text
u(u^3+1).                                                  (25)
```

Accordingly, the other six nonzero three-torsion points of `(24)` have

```text
x^3=1,                         w^2=-3.                     (26)
```

They suggest weighted initial ideals

```text
(x-lambda y^2, w-mu y^3),       lambda^3=1, mu^2=-3.      (27)
```

The naive polynomial graph lift already fails: substituting
`x=lambda y^2` into the affine resolvent branch gives

```text
y^6((1-lambda^2 y)^2-4),                                 (28)
```

whose residual quadratic has two distinct roots and is not a square. This
does not prove that the other elliptic directions never globalize. It says
exactly where they must hide: in splitting curves of degree at least three,
in non-plane-supported divisor classes, or not at all. The full class group
and involution-compatible nonmonogenic twist remain open.

## Reproduction

```bash
python3 04-computation/jc2_infinity_component_unique_splitting_conic_thm3934.py
python3 -O 04-computation/jc2_infinity_component_unique_splitting_conic_thm3934.py
```

Both runs must byte-match the frozen output. The assertion-free companion
checks the rank-six restriction map, all seven coefficient equations, every
highest-degree Groebner chart, the line and `h` components, splitting of
`XZ`, the finite parity bridge, and the exceptional-elliptic hostile lift.
It reports 26 exact gates. **QED.**
