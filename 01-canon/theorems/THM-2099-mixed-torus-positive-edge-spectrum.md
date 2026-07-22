---
id: THM-2099
title: "The primitive relation formula and a dyadic no-go for the rank-eight tree gate"
status: >
  PROVED. Every restricted pair weight in THM-2098 is an exact rational
  functional of the primitive three-character relation. For independent
  terminal characters, if A is the order of the guard modulo their character
  lattice, the weight is a finite average of explicit trapezoid-CDF increments;
  every nonfactor layer A>=2 has weight at least 1/98. Dependent terminal
  pairs reduce to THM-1166's one-dimensional overlap times 5/7. The proposed
  universal rank-eight closure tau>5/49 is false. For the dyadic affine
  pencils c_0=(1,0), c_k=(a,2^k), the sign of tau-5/49 obeys an exact
  seven-residue law: it is negative for a=1,3 mod 7 and zero for a=0 mod 7.
  At a=1, tau=70541889/691400192=5/49-9151/691400192. None is a mixed cover:
  every nondegenerate common affine pencil has a one-dimensional safe slice. Thus pair
  weights alone cannot close rank eight; relation/pencil or
  higher-intersection data is indispensable. This is not LRC(14).
source: codex-2026-07-22-LRC-mixed-torus-positive-spectrum
depends_on:
  - THM-1166
  - THM-2097
  - THM-2098
related:
  - THM-1221
  - THM-2096
  - THM-2103
  - THM-2120
script: 04-computation/lrc14_mixed_torus_edge_spectrum_scout_codex_20260722.py
output: 05-knowledge/results/lrc14_mixed_torus_edge_spectrum_scout_codex_20260722.out
script_sha256: cba376cfd8ef7297aa9dba2bd8f65bf8fdc5908687f0736a2fc1227871377ce2
output_sha256: 863552279335b1c41f75f73bd05391a3fe09f3bf9e1b9ede59f33336da20003b
hash_basis: working-tree files with LF line endings
---

# THM-2099 -- primitive relation spectrum and the dyadic no-go

Retain THM-2098's guard character `g` and two terminal characters `u,v`,
each transverse to `g`. Put

```text
w(u,v;g)=measure({X in T^2:
  ||u.X||<1/14, ||v.X||<1/14, ||g.X||>1/7}).           (1)
```

This theorem computes (1) from the primitive relation and then uses the
formula to challenge the hoped-for strict rank-eight tree gate.

## 1. Independent terminal characters

Assume `u,v` are linearly independent and let

```text
Lambda=Z u+Z v.
```

Let `A` be the order of `g` in the finite quotient of the ambient character
lattice by `Lambda`. Equivalently, there is a unique primitive relation up to
simultaneous sign

```text
A g=B u+C v,                  A>0,                     (2)
```

where `gcd(A,B,C)=1`. Transversality makes `B,C` nonzero.

For `b=|B|`, `c=|C|`, let `U,V` be independent uniform variables on
`[-1,1]`, and put `S=bU+cV`. Its exact CDF is

```text
F_(b,c)(s)=1/(8bc) [
 (s+b+c)_+^2-(s-b+c)_+^2
 -(s+b-c)_+^2+(s-b-c)_+^2].                           (3)
```

For `0<=k<A`, define the separate fiber-indexed periodic guard intervals

```text
I_(A,k)=union_(m in Z)
 [14Am-14k-2A, 14Am-14k+2A].                          (4)
```

Then

```text
w(u,v;g)=1/(49A) sum_(k=0)^(A-1) P(S notin I_(A,k)).  (5)
```

In (5), each probability is a finite sum of increments of (3), since
`S in [-b-c,b+c]`. In particular every weight is rational and can be
evaluated without torus meshing.

### Proof

The torus homomorphism

```text
Phi:X |->(u.X,v.X)
```

is surjective with finite fibers. The image of `g` on `ker(Phi)` is the
cyclic subgroup of order `A`. Hence, conditionally on

```text
u.X=x,                    v.X=y,
```

the `g`-values are

```text
(Bx+Cy+k)/A,                  0<=k<A,                  (6)
```

with equal multiplicity. On the two terminal danger bands write
`x=U/14`, `y=V/14`. The condition that (6) belongs to the radius-`1/7`
guard is precisely

```text
BU+CV in I_(A,k).
```

Signs of `B,C` do not change the distribution. The terminal box has Haar
mass `1/49`, giving (5). Formula (3) is the standard rectangle-slice CDF,
obtained by translating `[-b,b]x[-c,c]` to a positive rectangle and applying
two-variable inclusion-exclusion. QED.

### The nonfactor floor

If `A>=2`, the `A` values in (6) form a translated uniform `A`-grid. Away
from a null set of endpoint alignments, an arc of length `2/7` contains at
most `ceil(2A/7)` grid points. For every `A>=2`,

```text
ceil(2A/7)<=A/2.                                      (7)
```

The cases `A=2,3,4` are immediate; for `A>=5`, use
`ceil(2A/7)<=(2A+6)/7<=A/2`. Thus at least half of every fiber lies outside
the guard, and

```text
A>=2  implies  w(u,v;g)>=1/98.                        (8)
```

For `A=1`, equation (5) is the exact integral-factor spectrum. In the
unwrapped range `b+c<=7`, it reduces to the elementary trapezoid tail

```text
P(|bU+cV|>2)=
  0                              if b+c<=2,
  (b+c-2)^2/(4bc)                if |b-c|<=2<=b+c,
  1-2/max(b,c)                   if 2<=|b-c|.          (9)
```

The zero row `b=c=1` is THM-2097 pair rigidity. The first positive row
`{b,c}={1,2}` has conditional fraction `1/8` and weight `1/392`.
No global claim that the bounded rows of (9) exhaust the positive spectrum
is needed below.

## 2. Dependent terminal characters

If `u=a p` and `v=b p` for a primitive character `p`, then `p` and `g` are
independent. Their torus map is surjective, so the terminal-pair event and
the guard event are independent. Therefore

```text
w(u,v;g)=(5/7) rho(|a|,|b|),                          (10)
```

where `rho` is THM-1166's exact one-dimensional pair-overlap formula. This
completes the pair spectrum in every rank configuration.

## 3. The dyadic affine pencils defeat the strict tree gate

First take

```text
g=(1,0),             c_k=(1,2^k),        0<=k<=7.     (11)
```

All terminals are transverse to `g`; specializing along `(1,D)` for large
positive `D` gives an odd positive guard and distinct positive terminal
speeds. For a pair at gap `r`,

```text
(2^r-1)g=2^r c_k-c_(k+r).                              (12)
```

Formula (5) gives the exact gap spectrum

```text
r       1        2        3         4          5          6          7
w_r   1/392    2/147    5/343    227/15680  177/12152   5/343   46441/3186176.
                                                               (13)
```

The largest value is `5/343`, attained exactly at gaps three and six. Those
edges connect the vertices with equal residue modulo three. Their graph has
three components and a maximum forest with five edges. The next edge is the
unique gap-seven edge; the next weight is the gap-five layer, any one of
whose available edges joins the last component. Kruskal's algorithm therefore
gives

```text
tau=5(5/343)+46441/3186176+177/12152
   =70541889/691400192
   =5/49-9151/691400192
   <5/49.                                               (14)
```

Consequently the tempting strengthening

```text
every eight transverse characters have tau>5/49
```

is false. THM-2098's pair/tree budget cannot by itself exclude the dyadic
coefficient plane.

This is a certificate counterexample, not a mixed-cover counterexample. At

```text
X=(1/3,0)
```

one has `g.X=c_k.X=1/3` for all `k`, so the guard and all eight terminals are
strictly safe.

### The exact seven-residue law

The failure is not confined to signed copies of the guard. For any positive
integer `a`, put

```text
c_k(a)=(a,2^k),                    0<=k<=7,             (15)
```

and write `tau_a` for the corresponding maximum tree weight. At gap `r`,
the primitive relation is

```text
a(2^r-1)g=2^r c_k(a)-c_(k+r)(a).                       (16)
```

In formula (5), the number of `A`-grid points in the guard is the number of
integers in an interval of length `2A/7`, shifted by the fixed variable `S`.
Replacing `A` by `A+7d` increases this count by exactly `2d` away from null
endpoints. Hence, for fixed `r`,

```text
w_r(a)-5/343=kappa_(r,a mod 7)/a.                      (17)
```

Every spanning tree has seven edges. Taking the maximum over the finite tree
set preserves the scaling law: `a(tau_a-5/49)` depends only on `a mod 7`.
Exact evaluation of the seven residue representatives by (5) gives

```text
a mod 7       0              1                 2
a(tau_a-5/49) 0        -9151/691400192      213/24304

a mod 7       3              4              5          6
a(tau_a-5/49) -1031/41818560  1/49           3/49       33/392.   (18)
```

Thus the tree gate is below budget on two infinite residue classes and
exactly saturates on multiples of seven. The order-three doubling orbit
modulo seven is visible in the edge spectrum, but it does not imply a mixed
cover.

## 4. The signed affine-pencil escape

The last observation has a more general structural form. Suppose there are
signs `sigma_i in {+1,-1}` and characters `p,q` such that `p` is independent
of both `g` and `q`, and

```text
sigma_i c_i-q in Q p                     for every i. (19)
```

After replacing `p` by the primitive generator of the rank-one lattice in
(19), restrict to the circle `K={X:p.X=0}`. The restrictions of `g` and `q`
are two nonzero integer characters on `K`. On that circle the guard danger
has measure `2/7` and the common terminal danger has measure `1/7`; their
union has measure at most `3/7`. Choose `X in K` outside both. Then

```text
||g.X||>1/7,             ||c_i.X||=||q.X||>1/14.       (20)
```

Thus (19) always supplies a strict mixed-threshold escape. Every family
(15) has `p=(0,1)`, `q=(a,0)`, and all signs positive, so the complete
negative/equality residue law (18) lies in this geometrically safe class.

Hence any rank-eight mixed cover must satisfy the new sidecar condition

```text
rank_Q{sigma_i c_i-sigma_j c_j:1<=i,j<=8}=2
                 for every sign vector sigma,          (21)
```

unless the common difference line is guard-proportional or the common
quotient character vanishes on its rank-one kernel. This signed affine-pencil
rank is invisible to the edge-weight multiset and to the maximum tree value.

## 5. Assumption challenge and Tournament Analysis

The challenged assumption was that quantifying positive pair weights would
automatically repair the binary zero graph. Equations (14) and (18) disprove
that assumption at the exact THM-2098 budget. What these dyadic rows hide is
not a slightly better pair inequality but a common signed affine fiber,
exposed by (19).

The pairwise observable is (1). Thresholds `0`, `1/392`, `1/98`, and
`5/343` give nested undirected relations; formula (2) is their relation-height
sidecar. Orienting by weight, then breaking ties by label, gives a tournament
search scheduler. In the dyadic row it has many ties (gaps three and six
already coincide), and its score histogram, cycles, SCCs, edge flips, and
Hamiltonian paths depend on the tie rule. None detects (19). The faithful
carrier is the weighted graphic matroid together with the primitive relation
labels and signed-pencil rank. QED.

## 6. Exact referee and scope

The companion evaluates (3)--(5) with `Fraction`, independently checks the
small unwrapped rows, verifies every relation and weight in (12)--(13),
checks (18) through five full residue periods, and replays Kruskal's exact
sum (14). It also reports structured triadic,
consecutive, Fibonacci, and alternating controls and a seeded random scout;
those bounded scans are diagnostics, not proof of a tail statement.

The theorem removes a universal pair-tree route and adds the lawful rank
sidecar (21). It does not decide whether any signed-pencil-free rank-eight
plane covers, does not discharge THM-2098's vertical branch, and does not
prove LRC(14).
