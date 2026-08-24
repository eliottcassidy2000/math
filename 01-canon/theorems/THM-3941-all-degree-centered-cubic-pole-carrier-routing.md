---
id: THM-3941
title: "All-degree centered cubic pole carriers have three inertia grammars"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Let A(u) be a
  cubic polynomial and let t(u) be a
  nonconstant trace-zero rational function of arbitrary degree N. In the
  centered linear-color repeated-root grammar, any two finite poles in one
  A-fibre land at the same target address. Under the genuine reduced tame
  (2,1) maximal-ramification hypothesis used for a cubic Keller completion,
  they route to the non-unibranch source-ramification obstruction.
  After deleting that routed case, every finite pole lies at a critical point
  of A. Up to the stated affine gauges, there are exactly three carriers: one
  pole of order nonzero mod 3 over a C3 point, one odd pole over a selected C2
  point, or odd poles over both C2 points. The infinity order is zero or
  nonzero mod 3, and the pole orders sum to N. This is carrier exhaustion,
  not all-degree polynomial-color closure. At N=5 it gives exactly seven
  non-root-regular rows plus the root-regular exit. The carrier count has an
  exact period-twelve quadratic formula, and the deterministic carrier order
  gives every fixed-degree task a zero-based natural ordinal.
source: jc_zero_debt_lift / recursive synthesis of THM-3933, THM-3936, and THM-3938, 2026-08-24
audit: >
  FOUR INDEPENDENT HOSTILE AUDITS PASS (root, jc3913_lattice_referee,
  jc_decic_lattice, and jc_tournament_response, 2026-08-24). They reconstructed
  completed-local trace sieve, the explicitly qualified shared-fibre source
  obstruction, the three affine cubic critical packets, the unordered
  C2xC2 partition series, and the N=3,4,5 ledgers. The root audit repaired
  the frontmatter so that common target address is unconditional while source
  non-unibranchness retains its genuine maximal-ramification hypothesis. The
  other audits verified that the
  result counts trace-realizable support/order carriers rather than affine
  orbits or polynomial-color survivors, and that the root-regular handoff is
  one routing bucket retaining the unit-ideal hypothesis. The third audit
  supplied the all-order triangular trace-block realization lemma. The
  strengthened assertion-free 2858-gate companion proves the exact rational
  generating-function identity, the period-twelve quadratic count, and a
  carrier-to-natural ordinal bijection; an independent 368-gate path agrees
  through N=360. Normal and optimized runs byte-match the frozen output;
  raw and semantic hashes and documentation checks pass.
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
script_sha256: f04c9c8d8e6c2456ec04df32a4a1c10dd9a777fa98f1005676d1e0a8be84a1bc
output_sha256: 7bf0c868f35a38cea0d7ed01877b91c2abc001795f4ebbd921e20a0a41487d87
semantic_sha256: 29439314946684930dff9f379b2a048cbf808df926fc7992df0af2b210a75bf1
hash_basis: raw LF bytes
---

# THM-3941 -- every depth uses the same two local characters

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of
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

Here a carrier is only the projection type together with the support and
orders of the pole divisor after all lower Laurent jets have been forgotten.
It is not an affine orbit of the full pair `(A,t)`. This distinction matters
both for the exact count and for the root-regular exit below.

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
isolated rule must share its `A`-value with another pole. Section 2 identifies
their common target address unconditionally and routes them under the explicit
maximal-ramification hypothesis before any affine normalization.

The polynomial map `A:P1_u -> P1_A` is totally ramified of index three at
infinity. Write

```text
m=ord_infinity(t)_infinity >= 0.                          (11)
```

If `m>0`, the same trace calculation gives `3 does not divide m`. If `m=0`,
the constant term in `(9)` gives

```text
t(infinity)=0.                                            (12)
```

This last address is often useful in color division, but it is not an extra
carrier.

## 2. A shared pole value is a source, not merely target, obstruction

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

The qualifiers in the preceding paragraph are essential. A squared order
discriminant factor which disappears in the maximal order is not silently
promoted to a source ramification component here.

## 3. Three collision-free carrier gauges

After Section 2, every finite pole is isolated in its `A`-fibre. By `(10)` it
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

These admissible tuples are all nonempty at the trace layer. For `A=u^3`,
the finite block `u^(-r)` has trace zero whenever `3` does not divide `r`.
For the selected-C2 gauge `A=u^3+u^2`, put `v=1/u`. Then

```text
Av^3-v-1=0.
```

If `S_j=Tr(v^j)`, Newton recurrence gives

```text
S_1=0,  S_2=2/A,  S_j=(S_(j-2)+S_(j-3))/A.              (23a)
```

Inductively, the leading coefficients of `S_(2n)` and `S_(2n+1)` at `A=0`
are `2A^(-n)` and `(2n+1)A^(-n)`. Hence
`S_2,S_4,...,S_(2n)` form a triangular cancellation system with determinant
`2^n`. There is a unique trace-zero block

```text
v^(2n+1)+sum_(j=1)^n b_j v^(2j)                          (23b)
```

having exactly one pole, of order `2n+1`, at the selected C2 point. The
affine involution supplies the other C2 point, and blocks over the two
distinct critical values may be added. Finally, for any allowed positive
infinity order `m`,

```text
u^m-Tr(u^m)/3                                             (23c)
```

is a trace-zero polynomial of degree `m`: because `3` does not divide `m`,
the leading cubic-inertia character cancels in the trace correction. Adding
the relevant finite and infinity blocks realizes every tuple in `(22)`
without changing its leading pole orders. Thus the generating function below
counts actual trace carriers, though color may kill every one.

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

The separate root-regular generating function is `R3`. It records one
**routing bucket** in each allowed degree, not one affine orbit: with empty
finite support the inequivalent projections `A=u^3` and `A=u^3-3u`, as well
as polynomial lower-jet moduli, all enter the same THM-3929 exit.

### 4.1 Exact natural ordinals for the generated tasks

Let `c_N` denote the coefficient in `(25)`, so `c_N` counts non-root-regular
carrier rows, not color survivors. For `N=12j+r`, with `j>=0` and
`0<=r<12`, exact coefficient extraction gives

```text
c_(12j+r)=6j^2+b_r j+d_r,                                (25a)

r:    0  1  2  3  4  5  6  7  8  9 10 11
b_r: 16 10 14 16 16 14 22 16 20 22 22 20
d_r:  0  2  4  5  5  7 10  9 12 14 15 16.
```

This is an identity, not a sampled fit. Put `y=x^12`. The companion verifies
the rational-function equality

```text
M(R3+R2+P2)
 = sum_(r=0)^11 x^r [6y(1+y)/(1-y)^3
                     +b_r y/(1-y)^2+d_r/(1-y)].          (25b)
```

The standard generating series for `1,j,j^2` prove `(25a)` for every
`j>=0`; the harmless extension `j=r=0` says `c_0=0`.

There is also a literal procedural map to natural numbers. At fixed `N`,
order the nonregular triples in `(22)` first by increasing infinity order
`m`, then by carrier order

```text
C3 < C2 < C2 x C2,
```

and finally, in the two-pole carrier, by the smaller odd part. Their
zero-based ranks are exactly

```text
0,1,...,c_N-1.                                            (25c)
```

This ordinal is a task address, not a new invariant of the underlying
curve. It makes the recursive frontier precise: degree five has tasks
`0,...,6`, each asking for its lower-principal-part trace solution, exact
color division, coefficient ideal, maximal-order factor, and address fibre.

This formula shows precisely how the recurring binary and ternary grammars
meet: they are the nontrivial local characters of the `C2` and `C3` inertia
groups. It is not an analogy imposed after the calculation; equations
`(9)`, `(17)`, `(19)`, and `(23)` are the character projections themselves.

## 5. The inherited rows and the new degree-five frontier

At `N=3`, `(22)` gives the five nonregular rows used by THM-3933/3936:

```text
m=0: C2(3);
m=1: C3(2), C2xC2(1,1);
m=2: C3(1), C2(1).                                      (26)
```

At `N=4`, it gives the five nonregular rows of THM-3938 and one root-regular
exit:

```text
m=0: C3(4), C2xC2(1,3);
m=1: C2(3);
m=2: C3(2), C2xC2(1,1);
m=4: root-regular.                                       (27)
```

At `N=5`, there are exactly seven nonregular rows plus the root-regular exit:

```text
m=0: C3(5), C2(5);
m=1: C3(4), C2xC2(1,3);
m=2: C2(3);
m=4: C3(1), C2(1);
m=5: root-regular.                                       (28)
```

There is no `m=3` row. The next finite mathematical problem is now sharply
typed: impose the centered incidence and polynomial color `(7)` on the seven
rows in `(28)`. THM-3941 neither performs nor prejudges that color closure.

## 6. Exact reproducibility and scope boundary

Run

```bash
python3 04-computation/jc2_all_degree_centered_pole_carrier_routing_thm3941.py
python3 -O 04-computation/jc2_all_degree_centered_pole_carrier_routing_thm3941.py
```

The assertion-free companion verifies the three canonical critical packets,
the local `C2/C3` residue rules, the homogeneous shared-root equations
`a=C=0`, direct row enumeration, the independent generating function `(25)`,
and hostile degrees through 30. It prints the exact `N=3,4,5` row counts.

The theorem closes only the carrier layer under `(1)--(6)`. Arbitrary root
gauges, nonlinear color slots, polynomial-color division for unbounded `N`,
maximal-order verification for a new survivor, construction of a Keller
atlas, and JC(2) remain **OPEN**.
