---
id: THM-2073
title: "LRC(14) dyadic deletion tower and exact Z/4 seam"
status: >
  PROVED REDUCTION; NOT LRC(14). In a strict THM-2061 seam, every
  imprimitive deletion of the primitive eleven-core has gcd two and deletes
  its unique odd speed. The next four lifts carry a literal 2+1+1 partition
  owned by that guard and the two original odd tails. The construction
  iterates: every quotient remains primitive and divisor-complete through 14,
  exactly one safe child survives at each binary lift, and the chain ends at
  a hereditarily primitive core after at most eight levels by the theorem's
  internal cardinality argument; THM-2080 subsequently sharpens the live
  depth bound to four. The resulting
  2-adic normal form has one speed of each valuation below the terminal
  scale. Exact capacity, ownership, endpoint, and denominator-26 residue
  shells pass normally and under Python optimization.
source: codex-2026-07-21-LRC14-dyadic-deletion-tower
depends_on:
  - THM-2061
  - THM-765
  - LRCUpTo13
related:
  - THM-775
  - THM-2066
  - THM-2068
  - THM-2072
  - THM-2076
  - THM-2080
script: 04-computation/lrc14_dyadic_deletion_tower_referee_codex_20260721.py
output: 05-knowledge/results/lrc14_dyadic_deletion_tower_referee_codex_20260721.out
script_sha256: a65a72a538897d3d9b2f6a25ce8c9099e88679fe8292b2e7213d7816e5b095e4
output_sha256: c99d6dcfb2d0aae805492e2a8349ef6a1793eda4494cf54fb440aedac4d79e24
hash_basis: normalized repository blobs (LF)
---

# THM-2073 -- LRC(14) dyadic deletion tower

Put `delta=1/14` and, for a finite nonempty set `Q`, write

```text
phi_Q(t)=min_(q in Q)||qt||,
G_Q={t in R/Z:phi_Q(t)>=delta}.                          (1)
```

Assume throughout that

```text
S=2C union {x,y},    |C|=11,    gcd(C)=1,
x,y distinct positive odd integers,    M(S)<delta.      (2)
```

This is exactly the strict dyadic seam isolated by THM-2060/2061. In
particular, THM-2061 says that `C` contains a multiple of every integer from
`2` through `14`.

## 1. Every imprimitive core deletion is dyadic

First note that the thirteen-speed row `S` is hereditarily primitive. If a
twelve-speed deletion `P=S\{w}` had gcd `d>1`, primitivity of `S` would give

```text
D=d/gcd(d,w)>=2.
```

The settled twelve-speed LRC gives `M(P)>=1/13`. Part B of THM-765 would then
give

```text
M(S)>=min(1/13,(D-1)/(2D))=1/13>1/14,                  (3)
```

contrary to (2).

Now fix `u in C`, put `P=C\{u}`, and suppose

```text
d=gcd(P)>1.                                             (4)
```

Then

```text
d=2,    u is odd,    P=2V with gcd(V)=1.                (5)
```

Consequently either every one-element deletion of `C` is primitive, or

```text
C=2V union {u},    |V|=10,    gcd(V)=1,                 (6)
```

where `u` is the unique odd member and deleting `u` is the sole imprimitive
deletion.

### Proof

Primitivity of `C` gives `gcd(d,u)=1`. The deletion

```text
S\{2u}=2P union {x,y}
```

is primitive by (3), so `d` cannot divide both odd tails. Choose
`w in {x,y}` with `d` not dividing `w`, and put

```text
D=d/gcd(d,w)>=2.                                       (7)
```

Let `t_0` maximize `P`. Since `P` has ten speeds, settled LRC gives
`M(P)>=1/11>delta`. All `d` points

```text
t_j=t_0+j/d,    j in Z/dZ,                              (8)
```

are strictly `P`-safe. At each point either `u` is strictly
`delta`-dangerous, or `t_j in G_C` and THM-2061 forces the chosen odd tail
`w` to be strictly `1/7`-eligible. Thus the deck is covered by

```text
||u t_j||<1/14       or       ||w t_j||<1/7.            (9)
```

Because `u` is a unit modulo `d`, the first condition occupies at most the
fraction

```text
alpha(d)=(floor(d/7)+1)/d.                              (10)
```

The second orbit has order `D`, each orbit point repeated equally often, and
occupies at most

```text
beta(D)=(floor(2D/7)+1)/D.                              (11)
```

For every integer at least two,

```text
alpha(d)<=1/2,    equality iff d=2,
beta(D)<=1/2,     equality iff D in {2,4}.               (12)
```

For `d,D>=5` this follows immediately from `floor r<=r`; the cases `2,3,4`
are direct. Since (9) covers the entire deck, the two capacities sum to at
least one. Both must be equalities, and the first forces `d=2`. Hence `u` is
odd and `P=2V`; equation `gcd(P)=2` makes `V` primitive. Applying the same
argument to any imprimitive deletion says its deleted member must be the
unique odd member, proving uniqueness. QED.

The open inequalities in (9) are correct for a strict counterexample. The
closed interval counts (10)--(11) are upper bounds for those open teeth, so
endpoint coincidences can only decrease capacity and cannot invalidate the
argument.

## 2. The first Z/4 seam is a literal 2+1+1 partition

Assume (6), so

```text
S=4V union {2u,x,y}.                                    (13)
```

For every `sigma in G_V`, the following strict eligibility bounds hold:

```text
||u sigma||<1/7,
||x sigma||<2/7,    ||y sigma||<2/7.                    (14)
```

Let `N_z(sigma)` denote the unique nearest integer to `z sigma`. On the four
lifts

```text
t_j=(sigma+j)/4,    j in Z/4Z,                          (15)
```

the strict-danger ownership sets are

```text
2u: {j:j=N_u(sigma) (mod 2)},
x : {j_x},    j_x=-x^(-1)N_x(sigma) (mod 4),
y : {j_y},    j_y=-y^(-1)N_y(sigma) (mod 4).            (16)
```

They are a disjoint partition of `Z/4Z` of sizes `2+1+1`. The singleton
owners are the two distinct residue classes of the parity opposite to the
guard class, so

```text
j_x-j_y=2 (mod 4),
N_x(sigma)y-N_y(sigma)x=2 (mod 4).                      (17)
```

### Proof

Every speed in `4V` is weakly `delta`-safe on all four lifts. The speed `2u`
can kill at most two lifts, while each odd speed `x,y` can kill at most one.
Strict failure of `S` requires all four lifts to be killed. The capacities
therefore all saturate and their danger sets are disjoint. Solving the
nearest-integer congruences gives (14)--(16). The complement of one parity
class consists of the other two singleton classes, which differ by two.
Multiplying their congruence difference by the odd unit `xy` gives (17).
QED.

Exactly one of the two lifts

```text
(sigma+k)/2,    k in Z/2Z,                              (18)
```

belongs to `G_C`; the other is blocked by `u`. Thus `u` is a persistent
strict `1/7`-eligible guard, and the map from `sigma` to its unique safe lift
is the first **safe-child map**.

There is also a useful metric sidecar. Put

```text
B=max(V),    mu=M(V),    rho=(mu-delta)/B.               (19)
```

The `B`-Lipschitz property gives an interval of radius `rho` in `G_V` about a
maximizer. On its real lift, (14) and (17) give

```text
2/(xy)+2rho < 2/(7x)+2/(7y),
u<1/(7rho),    x,y<2/(7rho).                            (20)
```

Since `V` has ten speeds, `mu>=1/11`, so

```text
rho>=3/(154B),
u<22B/3,    x,y<44B/3.                                 (21)
```

These bounds are not used in the tower induction, but they give a finite
relative-height sidecar at its first seam.

## 3. Divisor completeness crosses the first seam

The quotient `V` contains a multiple of every

```text
m=2,3,...,14.                                           (22)
```

### Proof

If `V` missed `m`, every unit fraction `a/m` would lie in the closed set
`G_V`: each nonzero least residue is at least `1/m>=1/14`. The guard and
singleton bounds (14), required at all these fractions, have the exact
all-unit residue shells

```text
odd z mod 2m with 7|za|_m<m for all a in U_m
  = empty                       if m is even,
  = {m}                         if m is odd;             (23)

odd z mod 2m with 7|za|_m<2m for all a in U_m
  = {m}    for m=3,5,7,9,11,13.                         (24)
```

Equation (23) excludes even `m`. For odd `m`, (23)--(24) force

```text
u=m r_u,    x=m r_x,    y=m r_y
```

with all three quotients odd. At `sigma=1/m`, formula (16) gives

```text
j_x=j_y=-m^(-1) (mod 4),                                (25)
```

so `j_x=j_y`, contradicting the disjoint singleton owners. More explicitly,
`x^(-1)N_x=(m r_x)^(-1)r_x=m^(-1) mod 4`, and identically for `y`.
This proves (22). QED.

## 4. The full dyadic safe-child tower

There is a finite chain

```text
C=Q_0,
Q_i=2Q_(i+1) union {h_i},    h_i odd,    0<=i<r,        (26)
```

ending at a quotient `Q_r` whose one-element deletions are all primitive.
Every `Q_i` is primitive and contains a multiple of every `m=2,...,14`. At
each displayed seam:

1. deleting `h_i` is the sole imprimitive deletion of `Q_i`;
2. `h_i` is strictly `1/7`-eligible at every point of `G_(Q_(i+1))`;
3. exactly one of the two lifts of each such point lies in `G_(Q_i)`, while
   the other is blocked by `h_i`.

If `C` is already hereditarily primitive, take `r=0`. Otherwise Sections
1--3 construct `Q_1=V` and all three invariants.

### Inductive construction

Suppose `i>=1`, the invariants hold through

```text
Q_(i-1)=2Q_i union {h_(i-1)},                            (27)
```

and `P=Q_i\{q}` has gcd `d>1`. Primitivity gives `gcd(d,q)=1`. Moreover
`d` cannot divide `h_(i-1)`: if it did, then

```text
Q_(i-1)\{2q}=2P union {h_(i-1)}
```

would be imprimitive, contradicting the uniqueness of the deletion of
`h_(i-1)` one level above.

At a maximizer of `P`, all `d` gcd-deck translates are strictly
`P`-safe by settled lower-dimensional LRC. Each is covered either by the
`delta`-danger tooth of `q`, or, when it lies in `G_(Q_i)`, by the strict
`1/7`-eligibility tooth of `h_(i-1)`. The same capacities (10)--(12) force

```text
d=2,    Q_i=2Q_(i+1) union {q},    gcd(Q_(i+1))=1.       (28)
```

The same argument makes `q` the unique odd member and its deletion the sole
imprimitive deletion. Set `h_i=q`.

For `sigma in G_(Q_(i+1))`, the odd guard `h_i` can block at most one of its
two lifts. It must block one. Otherwise both lifts would lie in `G_(Q_i)` and
the preceding odd guard `h_(i-1)` would have to satisfy `||h_(i-1)t||<1/7`
at two antipodal phases. For odd `h_(i-1)`, those two distances sum to `1/2`,
so this is impossible. Hence exactly one safe child survives, regenerating
the eligibility and safe-child invariants.

### Divisor transfer at later seams

It remains to show that `Q_(i+1)` is divisor-complete. If it missed
`m in {2,...,14}`, every unit `a/m` would belong to `G_(Q_(i+1))`. The guard
shell (23) excludes even `m`; for odd `m`, it forces `h_i=m r` with `r` odd.
The lift blocked by `h_i` has even numerator over `2m`, while the unique safe
lift has odd numerator. As `a` ranges over `U_m`, the safe children are
exactly all unit fractions of denominator `2m`, and they lie in `G_(Q_i)`.
The preceding guard would therefore be `1/7`-eligible on every unit fraction
of denominator

```text
2m in {6,10,14,18,22,26}.                               (29)
```

No odd residue has that property. An elementary gcd-orbit check gives the
possible maximum least residues

```text
6 : {1,3},       10:{3,5},       14:{5,7},
18: {7,3,9},     22:{9,11},      26:{11,13};            (30)
```

each is at least the corresponding denominator divided by seven. Thus the
strict universal `1/7` inequality fails. This proves divisor transfer and
completes the induction.

Cardinality decreases at each step, so the process terminates. Unwinding
(26) gives the exact 2-adic normal form

```text
C=2^r Q_r union {2^i h_i:0<=i<r}.                       (31)
```

A depth-`r` tower has exactly one speed of each 2-adic valuation
`0,1,...,r-1`; all remaining `11-r` speeds are divisible by `2^r`. A
hereditarily primitive terminal core has at least three elements: a singleton
cannot be both primitive and divisor-complete, while a hereditarily primitive
two-set would require both singleton deletions to be `{1}`, contradicting
distinctness. Hence

```text
0<=r<=8.                                                (32)
```

THM-2076 subsequently sharpens (32) to `r<=5` by the strict Haar-capacity
tax of the terminal guard, and THM-2080 sharpens it again to `r<=4` by an
unequal-comb Hunter bound. The present theorem retains (32) as the elementary
cardinality bound proved internally here.

## 5. Frontier effect and guardrails

THM-2072 proves that a fixed owner-clock bank cannot close all unbounded
cores and that an antipodal pair in `G_C` closes the seam. The present theorem
supplies the structural replacement on the non-hereditarily-primitive lane:
the safe set is pulled through a finite sequence of binary fibers, with one
and only one safe child at each level. The terminal state is a primitive,
divisor-complete, hereditarily primitive core plus the ordered odd guards.

This does **not** prove LRC(14). It leaves two honest residuals:

1. hereditarily primitive quotient cores, including the terminal `Q_r`;
2. proving that no such terminal can support the inherited safe-child and
   original folded-owner constraints.

Nor does (31) bound the absolute speeds. It is a valuation normal form and a
finite-depth proof carrier, not a height bound. THM-2066/2068 closes bounded
terminal cores; an unbounded completion needs an adaptive clock, an antipodal
safe pair, or a terminal-core theorem.

## 6. Assumption challenge and Tournament Analysis

The challenged assumption is that the recursive vertices should be runners
or arcs. The faithful vertices are binary sheet obligations. At the first
seam, four sheet vertices carry the literal `2+1+1` partition (16); at later
seams, the safe-child map selects one vertex in each binary fiber. Quotienting
to 2-adic valuations preserves the normal form (31) but destroys which child
is safe, so the nearest-integer parity labels remain load-bearing.

One may orient guards by which labelled child they block, using reflection
of the two-sheet fiber as the switch/gauge and depth order as the tie
Hamiltonian path. The resulting tournament is transitive, with no directed
cycles or nontrivial SCCs. It discards the capacity equality, the all-unit
shells, and the safe-child function. The correct carrier is therefore the
rooted binary sheet tree with owner labels, not its tournament shadow.

## 7. Exact referee

The companion verifies (12) through `100,000`, checks the complete shells
(23)--(24) and (29), exhausts the four-lift congruence formulas through
denominator `80` and odd speed `49`, and checks over one million exact
antipodal exclusions. It uses explicit `require` calls, passes under
`python -O`, and the stored output ends in `PASS`. QED.
