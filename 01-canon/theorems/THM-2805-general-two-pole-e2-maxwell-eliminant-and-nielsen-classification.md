---
id: THM-2805
title: "General two-pole e=2 Maxwell eliminant and complete Nielsen classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  two-pole partition (d,N-d) with 2<=d<=N-d, the
  balanced e=2 response chamber is the common-bitangent locus of
  x^d(x-1)^(N-d).  After removing the two diagonal inflection collisions,
  an explicit Maxwell polynomial has degree N-4 and is saturated by the
  ordered Nielsen atlas.  Unequal pole parts give N-4 maps; equal parts give
  N/2-2 unmarked maps and N-4 ordered charts.  Together with the d=1 edge
  this yields binomial(N-2,2) two-pole e=2 maps.  This is a response-layer
  theorem, not Keller-chart entry, JC(2), or DC(2).
source: root/jc-e2-general-two-pole-maxwell-2026-07-28
depends_on:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
  - THM-2800-two-pole-two-double-zero-stieltjes-recurrence-and-first-nielsen-pair
related:
  - THM-2799-one-pole-two-double-zero-chebyshev-response-classification
script: 04-computation/jc2_general_two_pole_e2_maxwell_thm2805.py
output: 05-knowledge/results/jc2_general_two_pole_e2_maxwell_thm2805.out
script_sha256: bf1b1dcee6ef07920a80b16dcf2beb4be7c995b02eb341a27bda85171f36c9ab
output_sha256: 10908b5aae6b3e4b746f58240b15f9909b583efd378ab738c48bf1f673bec798
hash_basis: LF-normalized bytes
---

# THM-2805 -- general two-pole e=2 Maxwell eliminant and complete Nielsen classification

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2800 treats the edge partition `(N-1,1)`.  Nothing essential about the
bitangent mechanism is confined to that edge.  For every two-pole
partition, equality of tangent intercepts still eliminates the product of
the two double zeros linearly, and equality of slopes leaves a one-variable
Maxwell polynomial.  Its degree agrees exactly with the Nielsen count.

## 1. Arbitrary two-pole bitangent equations

Fix

```text
N>=4,                    2<=d<=N-d=:b,
e=2,                     h=2,
pole partition=(d,b).                                  (1)
```

Normalize the pole of order `d` to `0`, the pole of order `b` to `1`, the
high point over the third branch value to infinity, and
`F(infinity)=1`.  Put

```text
D=x^d(x-1)^b,                    T=x(x-1),
E=(x-gamma)(x-delta)=x^2-sx+p.                        (2)
```

As in THM-2800, the third partition `(N-1,1)` gives

```text
F=B/D,                    B=S E^2,
B-D=A(x-z),               A!=0.                       (3)
```

The two double-zero conditions say that the line `-A(x-z)` is tangent to
`y=D(x)` at both `gamma` and `delta`.  Since

```text
phi(x):=x-D(x)/D'(x)
 =x ((N-1)x-(d-1))/(Nx-d),                            (4)
```

equality `phi(gamma)=phi(delta)` is, after division by
`gamma-delta`,

```text
N(N-1)p-d(N-1)s+d(d-1)=0.                             (5)
```

Hence

```text
p=p_(N,d)(s)
 =[d(N-1)s-d(d-1)]/[N(N-1)],                         (6)

z=[(N-1)s-(d-1)]/N=(N-1)p/d.                         (7)
```

Define

```text
h_0=1,                  h_1=s,
h_k=s h_(k-1)-p_(N,d)(s) h_(k-2).                    (8)
```

Then `x^k=h_(k-1)x-p h_(k-2)` modulo `E`.

Expand

```text
D'(x)
 =sum_(j=0)^b c_j x^(d+j-1),
c_j=(d+j)(-1)^(b-j) binom(b,j).                      (9)
```

The common-slope equation is

```text
H_(N,d)(s)
 :=[D'(gamma)-D'(delta)]/(gamma-delta)
 =sum_(j=0)^b c_j h_(d+j-2)=0.                       (10)
```

Here `H_(N,d)` has degree `N-2`.

## 2. Remove exactly the two diagonal inflections

On the diagonal `gamma=delta=u`, equations `(5)--(6)` become

```text
J_(N,d)(s)
 :=N(N-1)s^2-4d(N-1)s+4d(d-1)=0,
s=2u.                                                (11)
```

For `2<=d,b`, both roots are non-pole roots of `D''`.  Indeed,
after the pole factors are removed, `D''(u)=0` is precisely `(11)` with
`s=2u`.  Thus both simple roots of `J_(N,d)` are roots of
`H_(N,d)`, and

```text
Q_(N,d)(s):=H_(N,d)(s)/J_(N,d)(s)                    (12)
```

is a polynomial.  Since the leading term in `(10)` is
`N s^(N-2)`,

```text
deg Q_(N,d)=N-4.                                     (13)
```

The case `N=4,d=b=2` gives a nonzero constant and hence no map, exactly as
the Nielsen count below predicts.

## 3. Algebraic converse

For a root of `Q_(N,d)`, initially assume the standard squarefree and
disjointness gates.  Define

```text
A=p sum_(j=0)^b c_j h_(d+j-3),                       (14)
B=D+A(x-z),
C=-(N-1)A,
S=B/E^2.                                             (15)
```

The convention is `h_(-1)=0`.  Equations `(8)--(10)` give

```text
D'=-A,                    D=-A(x-z)          modulo E,
```

so `E^2` divides `B`.  Equations `(6)--(7)` also give the exact identity

```text
D-(x-z)D'
 =-(N-1)x^(d-1)(x-1)^(b-1)E.                        (16)
```

Therefore

```text
F=B/D,
G=C E/(2DT),
V=4 S D T^2/C^2                                     (17)
```

satisfy

```text
F'=2G,                    F=VG^2,
2VG'+V'G=2.                                         (18)
```

As before, the provisional gate is

```text
Disc(E) A Disc(S)
Res(S,E x(x-1)) Res(E,x(x-1)) !=0.                  (19)
```

The all-degree Nielsen saturation proves that every root passes `(19)`.

## 4. Nielsen saturation

Put `L=N-1`, fix the third inertia `(0 1 ... L-1)` with one fixed sheet,
and normalize one zero transposition to `(*,0)`.  Write the other as
`(r,r+k)`.  The two pole-cycle lengths are

```text
k,                         N-k.                     (20)
```

For unequal parts `d<b`, condition `(20)=(d,b)` allows

```text
k=d:       N-2-d placements,
k=b:       d-2 placements.                          (21)
```

Their sum is

```text
N-4.                                                 (22)
```

For equal parts `d=b=N/2`, there are

```text
N/2-2                                                 (23)
```

unmarked Nielsen classes.  Marking which equal pole is sent to `0` doubles
this to `N-4` ordered charts.

Each ordered chart has a normalized rational realization by the Riemann
existence correspondence and hence supplies a root of `Q_(N,d)`.  The
algebraic formula `(14)--(18)` makes the map unique once that root is
fixed.  In the equal-pole case the two markings send

```text
s -> 2-s.                                            (24)
```

They cannot coincide: if `s=1`, the double zeros are exchanged by
`x->1-x`; symmetry gives `D'(1-x)=-D'(x)`, so equal slopes force both
points to be the single non-pole critical point `x=1/2`, the excluded
diagonal.

Thus the ordered Nielsen count equals the degree `N-4` in every case.
Consequently:

1. `Q_(N,d)` is squarefree;
2. all of its roots are admissible and pass `(19)`;
3. unequal pole parts give exactly `N-4` affine/target classes; and
4. equal pole parts give `N/2-2` unmarked classes, represented by `N-4`
   ordered roots paired by `(24)`.

This is the complete two-pole classification for `2<=d,b`.

## 5. Total `e=2,h=2` count

Add the edge partition `(1,N-1)` from THM-2800, which has `N-3` classes.
Summing `(22)--(23)` over unordered two-part partitions gives

```text
# {all e=2,h=2 response maps of degree N}
 =binom(N-2,2).                                      (25)
```

Thus the first passport multiplicity is part of a closed quadratic-size
atlas, not an isolated pair.

Together with THM-2799, only the three-pole chamber `e=2,h=3` remains in
the full `e=2` response layer.

## 6. Exact controls

The companion checks, with exact rational polynomial arithmetic:

1. every middle pole partition through `N=9`;
2. the intercept relation, recurrence, two-factor diagonal quotient, and
   degree `N-4`;
3. exact `E^2` division, the response identity, and all gates `(19)`;
4. complementary-pole symmetry under `x->1-x`, `s->2-s`;
5. unequal, equal, marked, and unmarked Nielsen counts;
6. the total formula `(25)` through `N=30`; and
7. normal, optimized, and stored transcript identity.

The script has no Python `assert` node.  Run:

```text
python 04-computation/jc2_general_two_pole_e2_maxwell_thm2805.py
python -O 04-computation/jc2_general_two_pole_e2_maxwell_thm2805.py
```

The finite controls support but do not replace the all-degree proofs.

## 7. Scope and next exact target

This classification remains at the abstract response layer.  It does not
establish Faber-flux compatibility, Keller-chart entry, JC(2), or DC(2).
Its concrete contribution is to remove all two-pole ambiguity from the
`e=2` program and to identify the missing integer passport coordinate as a
finite bitangent/Nielsen address.

The next sharp target is `e=2,h=3`.  Its three pole positions leave a true
cross-ratio before normalization.  The cheapest decisive question is
whether the two-double-zero conditions make that cross-ratio algebraic of
bounded degree, or whether this is the first response chamber with a
positive-dimensional accessory parameter.
