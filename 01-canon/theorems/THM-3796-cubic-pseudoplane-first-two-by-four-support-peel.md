---
id: THM-3796
title: "Cubic pseudo-plane first two-by-four support peel"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every exact
  2x4 or 4x2 Euler-support cell on the THM-3785 cubic pseudo-plane with
  exactly one contribution collision is empty, uniformly in the support
  gaps.  Every common-arithmetic-progression 2x4 or 4x2 cell of step one or
  two is also empty.  Step three and the first two-disjoint-pair cell survive
  the necessary sign gate and remain open controls; arbitrary
  multiple-collision cells, arbitrary Darboux pairs, and JC(2) remain open.
source: root / cubic_pseudoplane_low_support first-2x4 lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS by thm2x4_hostile_audit, 2026-08-23.  The
  audit rederived the unbounded collision-path grammar, every zero seam in
  the uniform one-collision contradiction, the exhaustive step-one and
  step-two endpoint rows, every ODE primitive and Delta-divisibility factor,
  output-swap symmetry, and the positive step-three/disjoint-pair controls.
  Its independent universe has 445,275 active gates.  The canonical
  assertion-free companion has 658,881 active gates; normal and optimized
  executions LF-normalize exactly to the frozen transcript.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
related:
  - THM-3787-cubic-pseudoplane-complete-low-support-darboux-nonentry
  - THM-3790-cubic-pseudoplane-arm-nodal-immersion-gate
  - THM-3792-pure-first-normal-nodal-carriers-have-critical-points
script: 04-computation/jc2_cubic_pseudoplane_first_2x4_peel_thm3796.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_first_2x4_peel_thm3796.out
script_sha256: 7b7b8aab503ddf41e656d65c19a5bf390d60849d736d982545871ff418836e6e
output_sha256: ba24fe2232b757c1912c0f5f7f1f3d4de19cef905342df1a24aa05f8d992d4c5
semantic_sha256: 8a23e174422c538f1d5553336e2d7b5b585aa547739c656470cf72b7915c6087
hash_basis: raw LF bytes
---

# THM-3796 -- the first cubic-pseudoplane two-by-four support peel

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let `k` be algebraically closed of characteristic zero, let `c in k*`, and
let

```text
B=k[R,Z,E]/(R^2E-Z^3+c^3R)
```

have the THM-3785 Poisson bracket and Euler weights

```text
wt(R,Z,E)=(3,1,-3),                 bracket shift=+2.
```

Suppose `F,G in B` have exact active supports

```text
supp(F)={a,a+r},
supp(G)={b,b+s,b+t,b+u},       r>0, 0<s<t<u.          (1)
```

An active weight-zero scalar may be translated away and is not counted;
an active nonconstant weight-zero profile remains in the support.

Then:

1. if the eight contribution addresses in `(1)` have exactly one collision,
   `{F,G}` cannot be a nonzero scalar;
2. if the supports have one common step

   ```text
   r=d,                    (s,t,u)=(d,2d,3d),
   ```

   and `d=1` or `d=2`, then `{F,G}` cannot be a nonzero scalar.

The same statements hold for exact `4 x 2` support after swapping the
outputs.  No statement is made for the remaining multiple-collision
`2 x 4` cells, arbitrary `3 x 3`, arbitrary Darboux pairs, or `JC(2)`.

## 1. Exact collision grammar

Put `V={0,s,t,u}`.  The contribution addresses, before the bracket shift,
are

```text
V union (r+V).                                         (2)
```

A collision at `q` is therefore exactly an edge

```text
v -> v+r                         inside V,             (3)
```

and its address is `q=v+r`.  The directed graph `(3)` is a disjoint union
of paths.  On four vertices it has at most three edges.  Two edges give
either a three-vertex path plus an isolated vertex or two disjoint pairs;
three edges give exactly the four-term AP

```text
V={0,r,2r,3r}.                                        (4)
```

This classifies the support grammar without a bounded search.

## 2. Uniform one-collision no-go

Assume `(3)` has exactly one edge `v -> v+r`, and let `q=v+r` be the unique
collision address.  A scalar must occupy this collision: a scalar in a
singleton bucket would be a homogeneous Darboux pair, excluded by THM-3785.
Thus

```text
a+b=-q-2.                                             (5)
```

The bottom address is singleton.  Two nonzero-weight homogeneous profiles
that commute have equal strict signs.  If one weight is zero and the other
nonzero, the zero-weight profile is scalar and removable; `(0,0)` remains a
lawful seam.  Since the sum in `(5)` is negative, exact support forces

```text
a<0,                         b<0.                     (6)
```

There are two `G` components outside the collided pair.  Every edge from
either `F` component to either outside component is singleton.  Commutation
with the lower `F` component makes both outside `G` weights strictly
negative; commutation with either outside component then makes the upper
`F` weight strictly negative.  The two cross-edges from the collided
components which do not lie at `q` are also singleton, so all four `G`
weights are strictly negative.  In particular the top singleton contains
two strictly negative integers, but its sum is

```text
(a+r)+(b+u)=u-v-2.                                    (7)
```

The left side is at most `-2`, so `(7)` would give `u<=v`.  This contradicts
`v+r in V` with `r>0`.  Hence every exactly-one-collision cell is empty.

This proof includes every zero seam: `(6)` rules out bottom zeros; either
outside hub rules out an upper-`F` zero; the cross-singletons then rule out
zeros in the collided `G` pair.

## 3. Complete endpoint sign rows for common step `d<=2`

For common step the collision addresses are `d,2d,3d`; only the bottom and
top addresses are singleton.  For a scalar address `q`, `(5)` holds and the
bottom weights are strictly negative.  Requiring the top singleton to have
equal strict signs or to be exactly `(0,0)` gives the complete table:

```text
d=1:
  q=1: no row
  q=2: (a,b)=(-1,-3)
  q=3: no row

d=2:
  q=2: (a,b)=(-1,-3)
  q=4: (a,b)=(-1,-5)
  q=6: (a,b)=(-2,-6).                                (8)
```

There are no omitted zero-weight cases: the only endpoint-zero rows are
displayed in `(8)`, and nonconstant weight-zero profiles are retained below.

Use the THM-3785 chart

```text
w=c+x^3y,                    Delta=w^3-c^3.
```

Every negative profile is divisible by `Delta`.

## 4. Common step `d=1`

The unique row in `(8)` is

```text
F: {-1,0},                    G: {-3,-2,-1,0}.        (9)
```

Write `F_-1=p`, `F_0=q`.  Bottom endpoint commutation gives

```text
G_-3=lambda p^3.                                      (10)
```

Writing `G_-2=m`, the first nonscalar collision is

```text
-p m'+2p'm+3lambda q'p^3=0,
```

so, in the rational function field,

```text
(m/p^2)'=3lambda q',
m=p^2(3lambda q+mu).                                  (11)
```

Let `G_-1=n`.  The scalar bucket is

```text
2q'm+p'n-pn'.                                         (12)
```

Both `p` and `n` are negative profiles.  Write
`p=Delta P`, `n=Delta N`.  Equation `(11)` makes the first term of `(12)`
divisible by `Delta^2`, while

```text
p'n-pn'=Delta^2(P'N-PN').                             (13)
```

Thus `(12)` is divisible by `Delta^2` and is not a nonzero scalar.

## 5. Common step `d=2`

### 5.1 Lower scalar row

Here

```text
F=(-1,1),                   G=(-3,-1,1,3).            (14)
```

Write the `F` profiles as `p,q`.  Endpoint commutation gives
`G_-3=lambda p^3`, `G_3=mu q^3`.  The upper collision integrates to

```text
G_1=q(3mu pq+nu).                                     (15)
```

The middle collision says `qG_-1-pG_1` is constant.  Since the valid
weight-`1` profile `q` is nonconstant, polynomiality forces that constant
to vanish, so

```text
G_-1=p(3mu pq+nu).                                    (16)
```

Substitution into the scalar bucket gives exactly

```text
3(lambda-mu)p^2(pq)',                                 (17)
```

which is zero or divisible by `Delta^2`.

### 5.2 Central scalar row

Here

```text
F=(-1,1),                   G=(-5,-3,-1,1).           (18)
```

Endpoint commutation gives `G_-5=mu p^5`, `G_1=lambda q`.  The top
collision forces

```text
q(G_-1-lambda p)=constant.
```

Again `q` is nonconstant, so the constant is zero and
`G_-1=lambda p`.  The bottom collision integrates to

```text
G_-3=p^3(5mu pq+nu).                                  (19)
```

The first scalar summand vanishes and the second has a factor `p^2`, hence
is divisible by `Delta^2`.

### 5.3 Upper scalar row

Here

```text
F=(-2,0),                   G=(-6,-4,-2,0).           (20)
```

Write the `F` profiles as `p,q`.  Bottom endpoint commutation gives
`G_-6=lambda p^3`.  The two lower collisions integrate successively to

```text
G_-4=p^2(3lambda q+mu),
G_-2=p(3lambda q^2+2mu q+nu).                         (21)
```

If the other weight-zero endpoint profile is `h`, the scalar bucket is

```text
2(q'G_-2-ph'),                                        (22)
```

and retains the factor `p`, hence the factor `Delta`.  This closes the
last row in `(8)`.

## 6. Output swap, controls, and boundary

For homogeneous profiles,

```text
[v,g;u,f]=-[u,f;v,g].                                 (23)
```

Therefore swapping outputs changes only the bracket sign and proves the
`4 x 2` statements with no new weight or zero seam.

The deterministic companion retains two positive open controls after every
necessary sign/zero filter:

```text
common AP d=3:
  q=3: (a,b)=(-2,-3),(-1,-4)
  q=6: (a,b)=(-2,-6),(-1,-7)

two disjoint r-pairs:
  r=4, V=(0,1,4,5), q=4: (a,b)=(-3,-3).              (24)
```

The proof does not close `(24)`.  These are the first coefficient systems
on which to combine multiple collision ODEs with the arm normal/residual
sidecar.

The exact companion is

```text
04-computation/jc2_cubic_pseudoplane_first_2x4_peel_thm3796.py
```

with frozen output beside it.  It uses no Python `assert`; normal and `-O`
executions LF-normalize exactly to the frozen transcript.

```text
script_sha256 = 7b7b8aab503ddf41e656d65c19a5bf390d60849d736d982545871ff418836e6e
output_sha256 = ba24fe2232b757c1912c0f5f7f1f3d4de19cef905342df1a24aa05f8d992d4c5
semantic_sha256 = 8a23e174422c538f1d5553336e2d7b5b585aa547739c656470cf72b7915c6087
active_checks = 658881
```

**QED.**
