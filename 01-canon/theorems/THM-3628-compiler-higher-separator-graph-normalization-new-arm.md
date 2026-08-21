---
id: THM-3628
title: "Compiler higher-separator graph normalization and new arm"
status: >
  PROVISIONAL PROVED + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.
  For every higher THM-3618 separator H=e^(m-1)z, m>=2, the graph
  closure has a uniform saturated two-relation presentation and image the
  complement of one new affine line. Its normalization is obtained by
  adjoining d and W=xe, not by reusing the THM-3622 normalization. The
  source omits the old d=0,-2 arms and a new d=-3 arm. The normalization is
  smooth for m=2 and has one A_(m-2) point for m>=3.
source: root / audit_thm3622 higher-separator continuation, 2026-08-21
audit: PENDING -- independent hostile reconstruction has not yet been completed.
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3618-compiler-one-graph-observable-fibre-separator-no-embedding
  - THM-3622-compiler-one-observable-graph-closure-normalization-arm-debt
related:
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
script: 04-computation/jc2_compiler_higher_separator_graph_normalization_thm3628.py
output: 05-knowledge/results/jc2_compiler_higher_separator_graph_normalization_thm3628.out
script_sha256: 574d75d7b5fb6afb72088abce687fc2b0d4d0d91fb4aa419fdc1e89a56927647
output_sha256: 112c4cfc7040b15cd2aa883a272570cdd54d05cf78902873e414e5cf3bd94b3a
hash_basis: raw LF bytes
---

# THM-3628 -- compiler higher-separator graph normalization and new arm

**PROVISIONAL PROVED + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.**

All rings and closed points are over `C`. Retain the compiler coordinates

```text
d=1+x^2q,
b=(d-1)(d+2)^2,
c=xd(d+2),
e=q(d+3),
z=xq,
A=C[x,q],                         R=C[b,c,e].          (1)
```

Fix an integer `m>=2`, put

```text
r=m-1,                    n=2m-1,
h=e^r z,                  S_m=R[h],
X_m=Spec(S_m).                                            (2)
```

Uppercase `B,C,E,H` denote the target coordinates. Addition of an element
of `R` and nonzero rescaling do not change the graph ring, so `(2)` covers
every higher separator in the THM-3618 classification.

## 0. Statement

The exact graph kernel has the uniform saturated presentation in section 1.
The source map

```text
psi_m:A2_(x,q) -> X_m                                  (3)
```

is injective, birational, and quasi-finite. Its exact set-theoretic image is

```text
psi_m(A2)=X_m minus V_infinity,

V_infinity={B=-4,C=E=0},                    H free.     (4)
```

Thus the visible complement remains one affine line, but for `m>=2` it is a
new vertical `H`-line rather than the `C`-line omitted when `m=1`.

Put

```text
W=xe=(d+3)z.                                           (5)
```

The normalization is

```text
N_m=S_m[d,W],                         Xtilde_m=Spec(N_m). (6)
```

It is finite and birational over `S_m`. The source is exactly

```text
A2 isomorphic to Xtilde_m minus (N_0 union N_-2 union V_-3), (7)
```

where the three omitted affine lines are

```text
N_0   ={d=0,  E=H=W=0},                 C free,
N_-2  ={d=-2, E=H=W=0},                 C free,
V_-3  ={d=-3, C=E=W=0},                 H free.        (8)
```

They are pairwise disjoint. The normalization is smooth when `m=2`. For
`m>=3` its only singular point is

```text
p_-3={d=-3,C=E=H=W=0},                              (9)
```

and its completed local equation is the rational double point

```text
CH=E^(m-1),                         type A_(m-2).       (10)
```

The old omitted arms acquire transverse numerical semigroup

```text
<2,2m-1>,                         delta=m-1,           (11)
```

while `V_-3` is a genuinely new `q`-infinity arm. In the common function
field,

```text
polar(x)=N_0+N_-2,
polar(z)=(m-1)V_-3,
polar(q)=2(m-1)V_-3.                                  (12)
```

## 1. Exact uniform saturated graph presentation

Let

```text
Delta=E^(2m-1)-H^2.                                   (13)
```

In `P=C[B,C,E,H]`, define

```text
F_B=B Delta^3-4H^2(3E^(2m-1)+H^2)^2,

F_C=C Delta^3
    -4H E^(m-1)(E^(2m-1)+3H^2)(3E^(2m-1)+H^2).        (14)
```

Then the exact kernel of

```text
P -> A,              (B,C,E,H) |-> (b,c,e,e^(m-1)z)  (15)
```

is the prime ideal

```text
J_m=(F_B,F_C):(E Delta)^infinity,
S_m=P/J_m.                                             (16)
```

The product in `(16)` selects a convenient dense inverse chart; no minimality
claim about the saturation multiplier is needed. On `E Delta!=0`, one has

```text
q=Delta/[4E^(2m-2)],
x=4H E^(m-1)/Delta,
d=[E^(2m-1)+3H^2]/Delta.                               (17)
```

Conversely `(17)` gives

```text
e=q(d+3)=E,                    e^(m-1)xq=H,            (18)
```

and substitution into the compiler gives exactly `(14)`. Therefore

```text
(P/J_m)_(E Delta)=C[E,H,(E Delta)^(-1)].               (19)
```

The open source locus in `(17)` is dense. Contracting its kernel is exactly
the saturation `(16)`, proving the scheme-theoretic kernel statement and
primality. No bounded-degree elimination assumption is used.

For `m=2` and `m=3`, independent direct elimination through the already
audited `z`-graph agrees exactly with `(16)`. In both cases the reduced lex
basis after imposing `E=0` is

```text
B(B+4),          (B+4)H^3,          CH^3,          E. (20)
```

The uniform reduced support is

```text
sqrt(J_m+(E))
 =(E,B,H) intersect (E,B+4,H) intersect (E,B+4,C).     (21)
```

Thus the `E=0` boundary consists of three lines:

```text
L_1     ={B=0, E=H=0},                    C free,
L_-4    ={B=-4,E=H=0},                    C free,
V_inf   ={B=-4,C=E=0},                    H free.       (22)
```

Only the reduced support is asserted uniformly in `(21)`; the exact
nonreduced fibre `(20)` is separately verified for the two hostile controls.

## 2. Exact image and the new visible line

On `E!=0`, relation `F_B=0` rules out `Delta=0`: otherwise
`H^2=E^(2m-1)` and its remaining term is nonzero. Thus equations `(17)`
recover the unique source point. The boundary
source inventory is equally explicit. If `e=0`, then either

```text
q=0:
  d=1,              (B,C,E,H)=(0,3x,0,0),

d=-3:
  x!=0, q=-4/x^2,   (B,C,E,H)=(-4,3x,0,0).             (23)
```

The first family covers all of `L_1`. The second covers
`L_-4 minus {C=0}`. No source point has `E=0,H!=0`, and the common endpoint

```text
p_infinity=(-4,0,0,0)                                  (24)
```

is also absent. Hence the missing set is exactly `V_infinity`, including
its origin, and `(4)` follows. Formula `(17)` and the two inverses in `(23)`
also prove injectivity without appealing only to point separation.
Equality of function fields proves birationality. Every closed fibre has at
most one point, so the finite-type morphism is quasi-finite.

The image is dense and open but not closed. It is not an open immersion:
sections 5 and 6 exhibit two normalization branches over each of the covered
lines `L_1` and `L_-4`, only one of which lies in the source.

## 3. Why the THM-3622 normalization cannot be reused

Let `X_1` be the `z`-graph of THM-3622. Its `E=0` boundary contains the
closed hyperbola

```text
D={B=-4,E=0,CZ=-12}.                                  (25)
```

Under the map `X_1 -> X_m`, `H=E^(m-1)Z`, this hyperbola maps to

```text
{B=-4,E=H=0,C!=0}=L_-4 minus {p_infinity}.             (26)
```

The image `(26)` is not closed. Consequently neither `S_1` nor its finite
normalization `T` from THM-3622 is finite over `S_m`. Equivalently, section
7 shows that

```text
z=H/E^(m-1)                                           (27)
```

has a pole on the new normalization divisor `V_-3`; it is not integral over
`S_m`. This is the smallest hostile obstruction to the tempting but false
claim that every separator graph has the same normalization.

## 4. Exact five-relation normalization

Put

```text
a=d(d+2),                    g=(d-1)(d+3).             (28)
```

In `C[d,C,E,H,W]`, let `K_m` be generated by the five relations

```text
CW=ag,
CE=aW,
W^2=gE,
(d+3)H=E^(m-1)W,
CH=a(d-1)E^(m-1).                                    (29)
```

Then

```text
N_m=C[d,C,E,H,W]/K_m=S_m[d,W].                         (30)
```

The source identities

```text
C=ax,              E=g/x^2,             W=g/x,
H=(d-1)g^(m-1)/x^(2m-1)                               (31)
```

prove that `(29)` maps to zero. Exactness is not inferred from these
identities alone. It follows from the three localization charts in section
5, which cover the quotient of `(29)`, are domains, and embed in the common
function field. Direct parametric elimination independently returns the
same ideal for `m=2` and `m=3`.

The map to the graph is

```text
B=(d-1)(d+2)^2,              (C,E,H)=(C,E,H).          (32)
```

It is finite because the two adjoined elements satisfy the monic equations

```text
d^3+3d^2-(B+4)=0,
W^2-gE=0.                                                (33)
```

It is birational by `(17)`. Section 5 proves that `(30)` is normal, so the
finite birational ring `(30)` is exactly the integral closure of `S_m` in
`C(x,q)`.

## 5. Three-open normality proof

The following principal opens cover `Spec(N_m)` because the three
polynomials

```text
g,                   a(d+3),                   a(d-1)  (34)
```

have no common zero in the `d`-line.

### 5.1 The `g!=0` chart

Equations `(29)` eliminate `E,H` and give

```text
(N_m)_g
 =C[d,C,W,g^(-1)]/(CW-ag).                             (35)
```

This hypersurface is smooth. A singular point would have `C=W=0`, hence
`a=0`; at `d=0,-2`, the derivative of `ag` is nonzero because `g` is a unit.

### 5.2 The `a(d+3)!=0` chart

Here `W=CE/a` and `H=CE^m/[a(d+3)]`, so

```text
(N_m)_(a(d+3))
 =C[d,C,E,(a(d+3))^(-1)]/(C^2E-a^2g).                 (36)
```

If `C!=0`, the derivative with respect to `E` is nonzero. If `C=0`, the
equation forces `d=1`, and the `d`-derivative is nonzero. Thus this chart is
smooth.

### 5.3 The remaining `d=-3` chart

On `a(d-1)!=0`, eliminate `W=CE/a`. Near the only locus not already covered,
`d=-3`, the two remaining equations are

```text
C^2E=a^2(d-1)(d+3),
CH=a(d-1)E^(m-1).                                     (37)
```

The derivative of the first right-hand side with respect to `d` is nonzero
at `d=-3`. Formal implicit elimination of `d+3`, followed by rescaling `H`
by a unit whose residue is `-12`, gives

```text
completed O_(Xtilde_m,p_-3)
 =C[[C,E,H]]/(CH-E^(m-1)).                             (38)
```

For `m=2`, `(38)` is smooth. For `m>=3`, it is the normal `A_(m-2)`
singularity and is singular only at its origin. Away from that origin one
of `C,H` is a unit and `(37)` is smooth. This proves normality of `(30)` and
the sharp singularity statement `(9)--(10)`.

The exact Jacobian controls agree: for `m=2`, the ideal of `(29)` plus all
`3x3` Jacobian minors is the unit ideal; for `m=3`, its reduced support is
exactly `(d+3,C,E,H,W)`.

## 6. Five boundary lines and the normalization-open source

Set-theoretically, imposing `E=0` in `(29)` first gives `W=0` and

```text
a g=0,               (d+3)H=0,               CH=0.   (39)
```

The normalized boundary therefore has five lines:

```text
d=1:       N_1,       E=H=W=0,                 C free,
d=-2:      N_-2,      E=H=W=0,                 C free,
d=0:       N_0,       E=H=W=0,                 C free,
d=-3:      D_-3,      E=H=W=0,                 C free,
d=-3:      V_-3,      C=E=W=0,                 H free. (40)
```

The last two meet at `p_-3`. Their images under normalization are

```text
N_1 and N_-2 -> L_1,
D_-3 and N_0 -> L_-4,
V_-3          -> V_infinity.                           (41)
```

The source covers all of `N_1` and `D_-3 minus {p_-3}`. It omits the three
pairwise-disjoint lines in `(8)`. This is an exact open statement. On the
complement of `N_0 union N_-2`, define

```text
x=C/a                    on a!=0,
x=g/W                    on W!=0.                      (42)
```

These opens cover, and the formulas agree because `CW=ag`. After also
removing `V_-3`, define

```text
q=E/(d+3)                on d+3!=0,
q=(d-1)a^2/C^2           on C!=0.                      (43)
```

These opens cover because a point with `d=-3,C=0` lies on `V_-3`. The
formulas agree by the derived relation

```text
C^2E=a^2g.                                             (44)
```

Equations `(42)--(44)` give

```text
d=1+x^2q,       C=xd(d+2),       E=q(d+3),
W=xE,           H=E^(m-1)xq.                           (45)
```

They construct the inverse to the source map and prove `(7)`.

## 7. Valuations and multiplied arm debts

At the generic points of `N_0` and `N_-2`, `C` is a unit and `W` is a
uniformizer. The leading packets are

| arm | leading terms in `W` |
|---|---|
| `N_0` | `d~-CW/6`, `E~-W^2/3`, `H~(-1)^(m-1)W^(2m-1)/3^m`, `x~-3/W`, `q~-W^2/9`, `z~W/3` |
| `N_-2` | `d+2~CW/6`, `E~-W^2/3`, `H~(-1)^(m-1)W^(2m-1)/3^(m-1)`, `x~-3/W`, `q~-W^2/3`, `z~W` |

Thus the target transverse ring on each omitted branch has numerical
semigroup

```text
C[[W^2,W^(2m-1)]],                                   (46)
```

whose delta invariant is `m-1`. Multiplication by `e^(m-1)` has changed the
old smooth transverse coordinate `z~W` into a cusp coordinate of order
`2m-1`.

At the generic point of `V_-3`, `H` is a unit and `E` is a uniformizer.
Writing `r=m-1`, equations `(29)` and `(37)` give

```text
C      ~ -12 E^r/H,
W      ~ -4 E^(r+1)/H,
d+3    ~ -4 E^(2r+1)/H^2,
x      ~ -4 E^r/H,
q      ~ -H^2/[4E^(2r)],
z       = H/E^r.                                      (47)
```

Hence `x` vanishes on the new arm, whereas `z` and `q` have pole orders `r`
and `2r`. The source-open complement `(7)` contains every other prime
divisor, proving the exact polar-divisor identities `(12)`.

The debt ledger is therefore:

| debt | normalization arm | target effect | visibility |
|---|---|---|---|
| old hidden collision debt | `N_-2` | shares `L_1` with covered `N_1` | invisible |
| transferred old closure debt | `N_0` | shares generic `L_-4` with covered `D_-3` | generically invisible |
| new higher-separator debt | `V_-3` | maps to omitted `V_infinity` | visible |

At `p_infinity`, both the `N_0` point and `p_-3` over `V_-3` are omitted.

## 8. Sharp comparison with `m=1`

For `m=1`, `H=z`; the relation at `d=-3` becomes

```text
CH=-12.                                                (48)
```

There is no finite `H`-axis and no `V_-3`. The normalization is the smooth
THM-3622 ring, the source omits only `N_0,N_-2`, and the visibly absent line
is the old `C`-line.

For `m=2`, the local relation is `CH=-12E`; the new arm appears, but the
normalization remains smooth.

For `m>=3`, the same new arm persists and its endpoint carries the exact
`A_(m-2)` singularity `(38)`. Thus the higher powers do not multiply the
number of missing divisors beyond three, but they increase both the cusp
delta `(11)` and the boundary singularity type.

## 9. Scope and exact companion contract

What is proved:

* the exact uniform saturated graph kernel `(16)` and reduced boundary
  support `(21)`;
* injectivity, birationality, quasi-finiteness, and exact image `(4)`;
* the exact five-relation normalization, its finiteness and normality;
* all boundary lines, the exact normalization-open source `(7)`, and inverse
  charts `(42)--(43)`;
* the three omitted arms, their valuations, semigroups, and polar divisors;
* the sharp `m=1`, `m=2`, and `m>=3` transitions;
* the explicit obstruction to reusing the THM-3622 normalization.

What is not claimed:

* that the normalization is smooth for `m>=3`;
* an exact uniform primary decomposition of the nonreduced `E=0` fibre;
* a classification for sums of multiple active separators;
* any counterexample to, proof of, or reduction of `JC(2)`.

The companion verifies, without truth-bypassing assertion statements:

* compiler identities and the generic inverse `(17)`;
* the saturation `(16)` against independent direct elimination for `m=2,3`;
* the exact fixed-case boundary basis `(20)` and all three boundary lines;
* the five normalization relations against direct parametric elimination;
* finite integrality, graph compatibility, chart identities, and Jacobian
  singular ideals;
* all five normalized boundary lines and the source inverse gluing;
* semigroup delta invariants and all leading valuation packets;
* the `m=1` hyperbola hostile witness and syntax hygiene.

The all-`m` saturation, chart normality, source-open identification, and
divisorial conclusions are proof-driven from uniform identities; the two
fixed cases are hostile controls, not finite substitutes for the theorem.
