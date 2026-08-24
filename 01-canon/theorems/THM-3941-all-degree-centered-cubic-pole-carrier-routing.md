---
id: THM-3941
title: "All-degree centered cubic pole carriers have three inertia grammars"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED AFTER MISTAKE-470
  QUANTIFIER REPAIR. Let A(u) be cubic and t(u) a nonconstant trace-zero
  rational function of arbitrary degree N in the centered linear-color
  repeated-root grammar. Its pole support obeys a trichotomy: no finite poles
  gives the root-regular exit; at least two finite poles share one A-value; or
  the nonempty collision-free support has exactly one of three carriers -- one
  pole of order nonzero mod 3 over a C3 point, one odd pole over a selected C2
  point, or odd poles over both C2 points. Only under
  the added genuine exact-double, reduced tame (2,1), maximal-normal Keller
  hypotheses does the shared case force a non-unibranch source arm. The
  infinity order is zero or nonzero mod 3, and all pole orders sum to N. This
  is carrier exhaustion, not all-degree polynomial-color closure. At N=5 it
  gives seven coarse non-root-regular signatures plus the root-regular exit.
source: jc_zero_debt_lift / recursive synthesis of THM-3933, THM-3936, and THM-3938, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS AFTER REPAIR (root audit and independent
  quantifier probe, 2026-08-24). The audit reconstructed the completed-local
  trace sieve, cubic critical packets, all-degree ledger and generating
  function, but found that the original opening silently imported the extra
  maximal-ramification hypotheses from Section 2. MISTAKE-470 records an
  infinite shared-pole counterfamily. The repaired theorem keeps that case as
  one branch of the pole-support trichotomy and invokes source non-unibranchness only with the
  stated Keller/maximal hypotheses. The assertion-free 2853-gate companion
  includes the hostile family; normal and optimized runs byte-match the frozen
  LF transcript, and raw and semantic hashes pass.
depends_on:
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
  - THM-3929-root-regular-one-place-linear-color-cubic-is-monogenic
  - THM-3936-centered-degree-three-infinite-root-value-nonentry
related:
  - THM-3933-centered-degree-three-root-map-pole-partition-octic-nonentry
  - THM-3938-centered-degree-four-root-map-nonentry
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
script: 04-computation/jc2_all_degree_centered_pole_carrier_routing_thm3941.py
output: 05-knowledge/results/jc2_all_degree_centered_pole_carrier_routing_thm3941.out
script_sha256: d3fc62f474fa9d0d615ca8202a4796de15b6dd50f1dbf5e3ca8909e426685243
output_sha256: 771726627e9f674d265e457e9c138aa1707802f549d481d8c4ede98f813d354e
semantic_sha256: 35b3b0eaffde82ab38c59a9b5d743a574a644b4e0ff938181f5dacaaef7ba9ef
hash_basis: raw LF bytes
---

# THM-3941 -- every depth uses the same two local characters

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED AFTER MISTAKE-470
QUANTIFIER REPAIR.** Work over an algebraically closed field `k` of
characteristic zero. Let an irreducible one-place component in the
linear-color binary-cubic grammar have normalization `A1_u`, with

```text
Phi=a(A)U^3+C U^2V+c(A)UV^2+d(A)V^3,                    (1)
A=A(u), C=C(u) in k[u],       deg A=3,                  (2)
k(A,C)=k(u).                                               (3)
```

Let its repeated-root map be the nonconstant rational function

```text
t=U/V in k(u),       deg(t:P1_u -> P1_[U:V])=N>=1.        (4)
```

We use the fixed centered root gauge. Assume that `t notin k(A)` and that
eliminating `C` from the two repeated-root equations gives the nonzero
primitive incidence

```text
a(A)t^3-c(A)t-2d(A)=0,                                   (5)
```

as the minimal polynomial of `t` over `k(A)`, up to a unit of `k(A)`. (Since
`[k(u):k(A)]=3` is prime, `t notin k(A)` is exactly the generation gate.) Its
missing quadratic coefficient gives

```text
Tr_(k(u)/k(A))(t)=0.                                     (6)
```

Under `(1)--(6)`, the pole support satisfies exactly one of the following
alternatives:

```text
(R) there are no finite poles (the root-regular exit);
(S) two finite poles share one A-value;
(C) the nonempty finite support is collision-free, and its carrier is one
    of C3, C2, or C2 x C2 from Section 3.                (6a)
```

In `(S)`, the repeated-root equations alone force a shared target address
`a(A0)=C=0`. They do **not** by themselves prove that the corresponding order
component remains genuine ramification in the maximal normalization. The
source non-unibranch exclusion needs the additional hypotheses stated in
Section 2. This is the MISTAKE-470 repair.

An arbitrary root translation need not preserve the linear `C U^2V` slot,
so `(1)` and `(5)` are hypotheses rather than a gauge available for every
binary cubic.

The theorem classifies **pole carriers** for every `N`. It does not assert
that the color reconstructed from

```text
C(u)=-(3a(A(u))t(u)^2+c(A(u)))/(2t(u))                   (7)
```

is polynomial, that the coefficient ideal is the unit ideal, or that any
listed carrier exists as a Keller completion. Those are the subsequent color
and surface tests.

## 1. The completed-local trace sieve

Let `p` lie over the finite value `A0`, with ramification index `e`. In tame
completed coordinates one may write

```text
k((z)) subset k((s)),              z=s^e.                (8)
```

For a Laurent series `f(s)=sum_j f_j s^j`, local trace is

```text
Tr(f)=sum_(xi^e=1) f(xi s)
     =e sum_l f_(el) z^l.                                (9)
```

Thus trace keeps exactly the trivial inertia character. Suppose a pole of
order `r` at `p` is the only pole over `A0`. If `e` divides `r`, its nonzero
leading coefficient survives in `(9)` and cannot be cancelled. Therefore

```text
e=1: no isolated pole;
e=2: an isolated pole has r odd;
e=3: an isolated pole has r not congruent to 0 mod 3.     (10)
```

If a pole is unramified, `(10)` says that trace zero requires at least one
other pole in the same `A`-fibre. More generally, every pole not obeying the
isolated rule must share its `A`-value with another pole. These are precisely
the shared-address alternative `(S)`; Section 2 explains when the Keller
application can route them to a source obstruction.

The polynomial map `A:P1_u -> P1_A` is totally ramified of index three at
infinity. Write

```text
m=max(0,-ord_infinity(t)) >= 0.                           (11)
```

If `m>0`, the same trace calculation gives `3 does not divide m`. If `m=0`,
the constant term in `(9)` gives

```text
t(infinity)=0.                                            (12)
```

This last address is often useful in color division, but it is not an extra
carrier.

## 2. A shared pole value becomes a source obstruction only in the maximal Keller application

Suppose distinct finite pole supports `p,q` satisfy `A(p)=A(q)=A0`. At both
supports the limiting homogeneous repeated root is `[U:V]=[1:0]`. In the
chart `U=1`, with `w=V/U`,

```text
Phi=a(A)+Cw+c(A)w^2+d(A)w^3.                             (13)
```

Multiplicity at `w=0` forces

```text
a(A0)=C=0.                                                (14)
```

Hence `p` and `q` are distinct normalization addresses of the same target
point `(A,C)=(A0,0)`.

For the Keller-completion application, assume the component is generically
exact-double and is a genuine reduced tame `(2,1)` ramification component of
the maximal cubic normalization. Let `B` be that normal finite
`k[A,C]`-algebra and `E` its residue-degree-one ramification prime. A normal
surface is Cohen--Macaulay, so miracle flatness makes `B` finite flat of rank
three over the smooth target. Every non-etale point of a geometric fibre has
local fibre length at least two. A rank-three fibre therefore cannot contain
two distinct points of `E`.

The two addresses `p,q` consequently map to one point of `E`. Since `E^nu`
is finite birational to the target component normalization, they are two
branches of `E` at that point. THM-3920 says that an irreducible boundary
curve of a normal affine surface containing an `A2` open is unibranch.
Deleting this ramification component can therefore never produce a Keller
`A2` open.

The qualifiers in the preceding paragraph are essential. Under `(1)--(6)`
alone, for every `q>=1` there is a shared-pole family

```text
A=u^3,                 t=u/(A-1)^q,
a=(A-1)^(3q),          c=0,                    d=A/2,
C=-(3/2)u(A-1)^(2q).                                  (14a)
```

It has trace zero, primitive minimal incidence

```text
(A-1)^(3q)T^3-A=0,
```

and `k(A,C)=k(u)`, but its three distinct unramified poles all lie over
`A=1`. Its one-place target component is

```text
C^3=-(27/8)A(A-1)^(6q).                                (14b)
```

Thus the opening grammar gives only the shared target address. A squared
order-discriminant factor which disappears in the maximal order cannot be
silently promoted to source ramification. MISTAKE-470 records the failed
implication and the repaired trichotomy.

## 3. Three collision-free carrier gauges

Restrict now to the complementary collision-free case `(C)`. Every finite
pole is isolated in its `A`-fibre. By `(10)` it
must lie at a finite ramification point of `A`. Riemann--Hurwitz gives total
ramification four for the cubic map. Infinity already contributes two, so
the finite ramification packet is exactly one of

```text
one point of index 3;                 or
two points of index 2.                                     (15)
```

Affine changes of `u` and `A` give the following complete list.

### Carrier C3

If the derivative has a double root, normalize it to zero:

```text
A=u^3.                                                     (16)
```

There is one finite pole support, at `u=0`, of order

```text
r>=1,                 r not congruent to 0 mod 3.          (17)
```

### Carrier C2

If only one of two simple critical points supports a pole, move the selected
point to zero and scale:

```text
A=u^3+u^2.                                                 (18)
```

The selected point is `u=0`; the unused critical point is `u=-2/3`, and their
critical values `0,4/27` are distinct. There is one finite pole order

```text
r>=1,                 r odd.                              (19)
```

Selecting the other critical point gives the same affine carrier orbit.

### Carrier C2 x C2

If both simple critical points support poles, normalize them to `u=1,-1`:

```text
A=u^3-3u.                                                  (20)
```

Their critical values are `-2,2`, so the supports do not share a target
`A`-value. The finite pole orders are

```text
r,s>=1,              r and s odd,                         (21)
```

and the involution `u -> -u, A -> -A` exchanges them. Thus only the unordered
pair `{r,s}` is carrier data.

No fourth carrier exists: an unramified pole was routed in Section 2, while
the derivative of a cubic has only the packets `(15)`.

## 4. The all-degree ledger

The degree of a rational map is the degree of its pole divisor. Consequently
every collision-free non-root-regular row is exactly one of

```text
C3:       N=m+r,       3 does not divide r;
C2:       N=m+r,       r odd;
C2 x C2: N=m+r+s,      r<=s, r and s odd,                 (22)
```

where in every row

```text
m=0 or (m>=1 and 3 does not divide m).                    (23)
```

If there are no finite poles, `t` is a polynomial, `m=N`, and `(23)` is the
root-regular exit routed to THM-3929. In particular no trace-zero
root-regular row exists when `3|N`.

For a compact all-degree count, put

```text
R3=(x+x^2)/(1-x^3),        R2=x/(1-x^2),
M=1+R3,
P2=x^2/[(1-x^2)(1-x^4)].                                 (24)
```

Here `M` records allowed infinity orders, `R3` and `R2` record one finite
pole, and `P2` records an unordered pair of positive odd orders. The number
of non-root-regular carrier rows of degree `N` is

```text
[x^N] M(R3+R2+P2).                                       (25)
```

The separate root-regular generating function is `R3`.

This formula shows precisely how the recurring binary and ternary grammars
meet: they are the nontrivial local characters of the `C2` and `C3` inertia
groups. It is not an analogy imposed after the calculation; equations
`(9)`, `(17)`, `(19)`, and `(23)` are the character projections themselves.

## 5. The inherited signatures and the new degree-five frontier

At `N=3`, `(22)` gives the five coarse nonregular carrier signatures used by
THM-3933/3936:

```text
m=0: C2(3);
m=1: C3(2), C2xC2(1,1);
m=2: C3(1), C2(1).                                      (26)
```

At `N=4`, it gives the five coarse nonregular signatures of THM-3938 and one
root-regular exit:

```text
m=0: C3(4), C2xC2(1,3);
m=1: C2(3);
m=2: C3(2), C2xC2(1,1);
m=4: root-regular.                                       (27)
```

At `N=5`, there are exactly seven coarse allowed nonregular signatures plus
the root-regular exit:

```text
m=0: C3(5), C2(5);
m=1: C3(4), C2xC2(1,3);
m=2: C2(3);
m=4: C3(1), C2(1);
m=5: root-regular.                                       (28)
```

There is no `m=3` signature. These are carrier possibilities, not existence
claims. The next finite mathematical problem is now sharply typed: impose the
centered incidence and polynomial color `(7)` on the seven signatures in
`(28)`. THM-3941 neither performs nor prejudges that color closure.

## 6. Exact reproducibility and scope boundary

Run

```bash
python3 04-computation/jc2_all_degree_centered_pole_carrier_routing_thm3941.py
python3 -O 04-computation/jc2_all_degree_centered_pole_carrier_routing_thm3941.py
```

The assertion-free companion verifies the three canonical critical packets,
the local `C2/C3` residue rules, the homogeneous shared-root equations
`a=C=0`, the MISTAKE-470 family `(14a)--(14b)`, direct row enumeration, the
independent generating function `(25)`, and hostile degrees through 30. It
prints the exact `N=3,4,5` signature counts.

The theorem closes only the carrier layer under `(1)--(6)`. Arbitrary root
gauges, nonlinear color slots, polynomial-color division for unbounded `N`,
maximal-order verification for a new survivor, construction of a Keller
atlas, and JC(2) remain **OPEN**.
