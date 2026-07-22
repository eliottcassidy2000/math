---
id: THM-2114
title: "Finite-ring Kakeya needles, the primitive toothpick graph, and rank-eight tree equality"
status: >
  PROVED. The exact moving-root law for a transverse terminal pair gives the
  sharp positive floor 1/392 and, under primitive guard/terminal hypotheses,
  identifies every edge below 1/147 with one of two self-similar toothpick
  moves. The resulting low graph has maximum degree four and is C4-free, so
  its rank-eight complement is connected and tau>=1/21. Independently, a
  missing common-nonzero vector modulo N=2,...,6 is an immediate finite-ring
  Kakeya escape; N=7 also escapes after a one-sided guard drift. Prime rows
  force a guard-or-terminal 13-content blocker through terminal rank 12 and
  an 11-content blocker through rank 10. The odd quotient-height branch has a
  uniform quarter-torsion witness. An exact eight-character row has tau=5/49
  but no signed affine-pencil gauge, independently strengthening THM-2103's
  no-go classification. It is nevertheless excluded twice: by a mod-5
  needle and by equality rigidity across all maximum graphic-matroid bases.
  These results replace the false tree-or-pencil binary target by a strict-
  tree / maximum-basis-equality / finite-torsion routing. They do not exclude
  every rank-eight coefficient plane or prove LRC(14).
source: codex-2026-07-22-LRC-rank-eight-needles-and-tree-equality
depends_on:
  - THM-2097
  - THM-2098
  - THM-2099
related:
  - THM-1166
  - THM-1221
  - THM-2096
  - THM-2103
  - THM-2104
  - THM-2105
  - THM-2116
  - THM-2120
script: 04-computation/lrc14_quarter_torsion_toothpick_codex_20260722.py
output: 05-knowledge/results/lrc14_quarter_torsion_toothpick_codex_20260722.out
script_sha256: 2293c15edee4a81a5f0178857903c613fc932b30536cfce3fb6013a2d73da182
output_sha256: 72bec20c3f07c04a5f02875851aa19c48cb0f6ee34daa95d8e71f88346087ddb
hash_basis: working-tree files with LF line endings
---

# THM-2114 -- finite-ring needles and the rank-eight equality carrier

Retain THM-2098's connected two-torus.  Write

```text
I={u in R/Z:||u||<1/14},
J={u in R/Z:||u||>1/7},
w(f,f';g)=measure{f in I, f' in I, g in J}.            (1)
```

All boundaries below are Haar-null.  This theorem first sharpens the local
pair spectrum, then records two pieces of structure that pair weights forget:
small torsion needles and simultaneous equality over all maximum trees.

## 1. The moving-root law

Let `f,f'` be independent characters, each transverse to `g`.  Put

```text
L=Z f+Z f'
```

in the rank-two character lattice, and let `m` be the order of `g+L` in the
finite quotient.  There is a primitive relation

```text
m g=a f+b f',                 m>0, a*b!=0.             (2)
```

Conditionally on `f=u,f'=v`, the values of `g` are, with equal multiplicity,
the `m` roots of

```text
m z=a u+b v.                                            (3)
```

Define

```text
N_m(s)=#{z in R/Z:mz=s and ||z||<1/7}.
```

If

```text
2m=7q+r,                    0<=r<7,                    (4)
```

then, almost everywhere in `s`,

```text
N_m(s)=q+1_{||s-q/2||<r/14}.                            (5)
```

Here `q/2` is read modulo one.  Indeed, the roots in (3) are an evenly spaced
`m`-grid, and counting them in the moving interval of length `2/7` is the same
as counting integers in an interval of length `2m/7`.  Such a count is `q` or
`q+1`; its excess interval has length `r/7`, and endpoint parity places its
center at `q/2`.  This also proves (5) when `r=0`, away from the null endpoint
alignment.

Haar disintegration over `(f,f')` now gives the exact formula

```text
w=(m-q)/(49m)
  -(1/m) measure{(u,v) in I^2:
                 ||a u+b v-q/2||<r/14}.                (6)
```

Thus the cardinality of the full kernel of `(f,f')` is irrelevant.  The only
finite-fiber datum seen by the guard is its monodromy order `m`.

The correction set in (6) has area at most `1/49`, with zero correction when
`r=0`.  Therefore

```text
w >= [m-ceil(2m/7)]/(49m).                              (7)
```

For `m>=2`, `ceil(2m/7)<=m/2`: check `m=2,3,4`, and for `m>=5` use
`2m/7+1<=m/2`.  Hence

```text
m>=2  implies  w>=1/98.                                (8)
```

This floor is sharp at `2g=f+f'`.  It is a floor statement, not pointwise
monotonicity in `m`: root averaging can lower a larger coefficient-layer
weight.

### The primitive factor layer

Assume now that `g` is primitive and `m=1`.  Then `g=a f+b f'` forces
`gcd(a,b)=1`.  Put `1<=A=|a|<=B=|b|`.  On scaling the terminal box by `14`,

```text
49w=Pr{dist(AU+BV,14Z)>2},       U,V uniform on [-1,1]. (9)
```

If `B=1`, coprimality gives `(A,B)=(1,1)` and (9) is zero.  If `B=2`, it gives
`(A,B)=(1,2)`; the two corner triangles in (9) have total relative area
`1/8`, so `w=1/392`.

Suppose `B>=3` and fix `U`.  As `V` varies, `AU+BV` traverses an interval of
length `2B`.  For `3<=B<=6`, this interval meets the periodic guard-danger
intervals, each of length four and separated by ten, in total length at most
four.  Its outside fraction is at least `1-2/B>=1/3`.  For `B=7k+j>=7`,
`0<=j<=6`, the first `14k` units contribute exactly `10k` outside units, so

```text
outside fraction >=10k/(14k+2j)>=1/3.                  (10)
```

Consequently, in the primitive factor layer,

```text
w=0       iff {A,B}={1,1},
w=1/392   iff {A,B}={1,2},
otherwise w>=1/147.                                    (11)
```

Together with (8), this gives the global positive floor `1/392` for an
independent primitive triple.

## 2. The primitive toothpick threshold graph

Assume `g` and eight distinct terminal bands are primitive and transverse.
After `GL_2(Z)` normalization and sign choices, write

```text
g=(1,0),                    f_i=(x_i,d_i),
d_i>0,                      gcd(x_i,d_i)=1.             (12)
```

For a pair put `s=gcd(d_i,d_j)` and

```text
Delta=d_j x_i-d_i x_j.
```

Its primitive relation is

```text
(Delta/s)g=(d_j/s)f_i-(d_i/s)f_j.                      (13)
```

By (8)--(11), the graph `N` of edges with `w_ij<1/147` has exactly the moves

```text
(x,d) -- (x+/-1,d),
(x,d) -- (2x+/-1,2d).                                  (14)
```

The first is the zero layer and the second is the `1/392` layer.  Equation
(14) is a literal self-similar toothpick graph: every level sprouts two
dyadic children.

Primitivity improves the naïve degree-six count to

```text
Delta(N)<=4.                                           (15)
```

If `d` is odd, a vertex has no integral parent and has at most two same-level
neighbors and two children.  If `d` is even, primitivity makes `x` odd, so
both same-level candidates are nonprimitive; only two parents and two
children remain.

The graph is also `C_4`-free.  More strongly, two primitive vertices have at
most one common neighbor.  If their heights are `d<=e`, a common-neighbor
height belongs to

```text
{d/2,d,2d} intersect {e/2,e,2e},
```

so `e/d` is `1,2`, or `4`.

- For `e=4d`, the only height is `2d`; the coordinate sets
  `{2x+/-1}` and `{(y+/-1)/2}` have spacings two and one and meet at most once.
- For `e=2d`, height `d` and height `2d` cannot both contribute: their
  coordinate equations make `y-2x` respectively odd and even.  Two common
  neighbors at height `2d` would force `(y,e)=(2x,2d)`, nonprimitive.
- For `e=d`, common neighbors at heights `d,2d,d/2` require respectively
  `|x-y|=2,1,2`.  In the only apparent double case, `d` is even and `x,y`
  are odd; the same-level midpoint has both coordinates even and is not a
  vertex.

Thus `N` has no `C_4`.  Its complement `H` on eight labels is connected.  If
it were disconnected, (15) gives every `H`-component at least two vertices;
two vertices in one component and two in another would form a `K_(2,2)` in
`N`, hence a `C_4`.  Every `H`-edge has weight at least `1/147`, so

```text
tau>=7/147=1/21.                                       (16)
```

The rank-eight covering cap is `5/49=15/147`; (16) leaves the honest gap
`8/147`.  It is a structural floor, not rank-eight closure.

## 3. Finite-ring Kakeya needles

Let `c_0=g,c_1,...,c_8 in Z^2` be the guard and terminal characters.  Fix
`N in {2,3,4,5,6}`.  If there is

```text
z in (Z/NZ)^2  with  c_i.z !=0 mod N for every i,       (17)
```

then any integer lift gives a torus point `X=z/N`.  Every character value is
nonzero `N`-torsion, hence has circle distance at least

```text
1/N>1/7.                                               (18)
```

It is therefore simultaneously guard-safe and terminal-safe.  In particular,
any mixed cover must satisfy the finite blocking conditions

```text
(Z/NZ)^2=union_(i=0)^8 ker(c_i mod N),
                         N=2,3,4,5,6.                  (19)
```

For prime `N`, (19) is a projective-line Kakeya condition: unless one
character vanishes modulo `N`, the nine kernel directions must exhaust
`P^1(F_N)`.

There is a useful row-count form beyond the all-nonzero test.  Work in the
actual saturated character lattice of the connected coefficient torus, and
choose a unimodular `Z`-basis `(h,p)`, with `p` the primitive guard-ray
generator, such that

```text
g=B p,                    c_i=d_i h+k_i p.             (19a)
```

Let `2<=Q<=13` and choose `a mod Q` with
`||Ba/Q||>1/7`.  On the torsion row

```text
(h.X,p.X)=(b/Q,a/Q),
```

terminal `i` is dangerous exactly when

```text
d_i b+k_i a=0 mod Q.                                  (19b)
```

Every nonzero `Q`-residue has distance at least `1/Q>1/14`.  The number
`N_i(Q,a)` of solutions in `b` is either zero or `gcd(d_i,Q)` (with the
usual `Q` solutions when both coefficients in (19b) vanish).  Therefore

```text
sum_i N_i(Q,a)<Q  implies a strict mixed escape.       (19c)
```

For a prime `p in {11,13}`, if `p` does not divide `B` and no terminal column
`(d_i,k_i)` vanishes modulo `p`, choose `a` so that `Ba` has a residue of
distance greater than `p/7`.  Each terminal forbids at most one value of `b`.
Thus any terminal rank `n<p` satisfies (19c).  In particular, every rank-eight
cover must obey the divisibility handoff (whose counting cutoff `n<p` is sharp)

```text
p divides B  or  p divides both d_i,k_i for some i,
                                      p=11 and p=13.   (19d)
```

This routes failures directly into the guard-proportional and p-adic ledgers.
Composite values `Q=9,10,12` retain the more flexible gcd invoice (19c).
The coordinate-free consequence, with `cont(c)` denoting content in that
saturated lattice (and hence invariant under `GL_2(Z)`), is

```text
13 divides cont(g) product_i cont(c_i)   for every rank n<=12 cover,
11 divides cont(g) product_i cont(c_i)   for every rank n<=10 cover.        (19e)
```

Hence a rank-two cover of terminal rank at most twelve cannot have all guard
and terminal characters primitive.  This applies to the transverse and
guard-proportional branches alike and is a new arithmetic filter on all of
THM-2098's ranks `8..11`, not only on the example below.

If an integer direction `d` specializes the characters to the guard `h=g.d`
and terminal speeds `q_i=c_i.d`, divisibility of content implies divisibility
of the specialization.  Thus (19e) has the directly usable scalar shadow

```text
13 divides h product_i q_i       for terminal rank n<=12,
11 divides h product_i q_i       for terminal rank n<=10.          (19f)
```

The converse need not hold, so (19e) remains the stronger invariant.  But in
the depth-four LRC terminal sizes `7..10`, every surviving rank-two cover must
display both an `11`-divisible and a `13`-divisible guard/terminal entry.  This
is an exact two-prime invoice before any tree computation.

The same missing-kernel implication holds at `N=7` after a one-sided drift.
At `X_0=z/7`, every terminal has distance at least `1/7>1/14`.  If the guard
has distance at least `2/7`, it is already safe.  If its distance is exactly
`1/7`, choose a real direction which moves the local signed lift away from
zero (its derivative has the same sign as that lift), and move to
`X_0+epsilon H`.  For sufficiently small positive `epsilon`, the
guard becomes strictly safe while every terminal retains its `1/14` margin.
Thus a cover must satisfy (19) for `N=7` as well.  When no character is zero
modulo seven, its nine normals must realize all eight directions of
`P^1(F_7)`.  This is the exact finite-field needle/Fano-style sidecar that a
pair-weight multiset erases.

### Uniform quarter torsion

There is a useful composite-torsion specialization.  Let `p,f` be independent,

```text
g=c p,                         c odd,
c_i=b_i p+a_i f,               every a_i odd.           (20)
```

Surjectivity of `(p,f)` supplies a point with

```text
p.X=1/2,                       f.X=1/4.                 (21)
```

Then

```text
||g.X||=1/2,
||c_i.X||=||b_i/2+a_i/4||=1/4.                          (22)
```

This is a strict mixed escape independent of every `b_i`.  In an LRC
coefficient plane, odd specialization of the guard forces its content `c` to
be odd.  Hence a surviving cover cannot have all transverse quotient heights
odd.  This branch is broader than one affine pencil: only the common nonzero
class modulo two in the quotient by the guard ray is used.

## 4. A clean guard-pencil no-go for the strict tree target

Take

```text
g=(1,0),                     f_k=(k,1).                 (23)
```

For a pair at distance `D=|k-l|`, conditionally on the two terminal values
the guard lies on a shifted `D`-grid.  Write `D=7q+r`, `0<=r<7`.  If `r>0`,
the guard-danger count is constantly `2q+1`; if `r=0`, it is `2q` away from a
null alignment.  Therefore the normalized pair weight `P_D=49w_D` is

```text
P_D = 5/7                         if r=0,
P_D = (5q+r-1)/D                  if 1<=r<=6.           (24)
```

For

```text
K={0,7,14,15,16,17,24,31},                               (25)
```

the complete normalized spectrum is

```text
5/7:6, 22/31:1, 17/24:2, 12/17:3, 7/10:2,
11/16:2, 2/3:5, 5/8:2, 1/2:2, 0:3.                    (26)
```

Kruskal first takes four `5/7` edges spanning the two triples
`{0,7,14}` and `{17,24,31}`, then the gap-31 edge, then the two gap-16 edges
that attach `15,16`.  Hence

```text
49 tau=4(5/7)+22/31+2(11/16)=8579/1736,
tau=8579/85064=5/49-101/85064.                          (27)
```

The direction `(1,2)` specializes the guard to one and the terminals to

```text
2,9,16,17,18,19,26,33,                                 (28)
```

so all LRC positivity/distinctness conditions survive.  Nevertheless (20)
applies with every quotient height equal to one; `(1/2,1/4)` is explicitly
safe.  Thus a tree deficit can encode a torsion hole rather than a near-cover.

## 5. A cap-equality refuter and maximum-tree equality rigidity

The binary target already refuted by THM-2103 also fails at the exact covering
cap.
Let `g=(1,0)` and, in label order, put

```text
V=((0,7),(-4,10),(11,2),(-18,10),
   (1,8),(-19,9),(-3,18),(9,2)).                        (29)
```

The direction `(1,3)` gives guard one and distinct positive terminal speeds

```text
21,26,17,12,25,8,51,15.                                (30)
```

Put `P_ij=49w_ij` and `S={0,1,3,4,5,6,7}`.  Exact evaluation by (6) gives

```text
P_ij=5/7  for all ij in K_S and for ij=02,              (31)

P_12=42/59,     P_23=52/73,       P_24=61/86,
P_25=2347/3288, P_26=655/918,     P_27=1/2.             (32)
```

Thus every edge is at most `5/7`, while the equality graph is `K_7` on `S`
plus the pendant edge `02`.  It is connected, and

```text
tau=7(5/343)=5/49.                                     (33)
```

No sign gauge puts the eight vectors in an affine line.  Since global sign is
irrelevant, fix `sigma_0=1`; for each of the `2^7=128` remaining gauges the
exact referee checks the rank-one criterion

```text
det(sigma_i V_i-V_0, sigma_j V_j-V_0)=0 for all i,j.    (34)
```

The number of passing gauges is zero.  Equations (33)--(34) give an independent
refuter of the proposed `tau>5/49` or signed-affine-pencil dichotomy.

This is not a cover.  Modulo five, the vector `z=(1,-1)` pairs nontrivially
with `g,V_0,...,V_7`; the residues are

```text
1,3,1,4,2,3,2,4,2.                                    (35)
```

The point `(1,-1)/5` is the finite-ring needle from Section 3.

There is also a purely graphic-matroid exclusion which exposes the next
rank-eight gate.  For any spanning tree `T` and active-label set

```text
A(X)={i:X in B_i},
```

the pointwise Hunter slack is

```text
|A|-e_T(A)-1_(A nonempty)>=0.                           (36)
```

Because `T[A]` is a forest, equality in (36) holds exactly when `A` is empty
or `T[A]` is connected.  In an eight-band transverse cover, if
`sum_(ij in T)w_ij=5/49`, the integrated Hunter inequality is equality.
Therefore almost every nonempty `A(X)` must induce a connected subgraph in
**every** maximum tree `T`.

This has a reusable higher-order form.  Let `H_M(i,j)` be the union of the
vertex sets of the unique `i`--`j` paths over all maximum spanning trees.  If
`i,j` are active, connectivity in each maximum tree forces every vertex on
each such path to be active.  Hence the **maximum-basis path hull** satisfies

```text
B_i intersect B_j intersect C
  subset intersection_(v in H_M(i,j)) B_v              (36a)
```

almost everywhere.  The left side is only pairwise data, while the right side
can be an arbitrarily higher collision.  This is the precise information
discarded by retaining only the scalar value `tau`.

For (31), every spanning tree of `K_7(S)` together with `02` is maximum.  A
nonempty label set connected in all of them is exactly one of

```text
a singleton,       {0,2},       S,       S union {2}.   (37)
```

Indeed, any proper subset of `S` of size between two and six is disconnected
by a star centered in its nonempty complement.  The pendant vertex `2` then
gives the remaining cases in (37).

Equivalently, `H_M(2,i)` is all eight labels for every
`i in S minus {0}`.  The event `B_2 intersect B_i intersect C` can therefore
occur only in the all-eight state.  Its measure would be independent of `i`,
forcing

```text
w_21=w_23=w_24=w_25=w_26=w_27,                          (38)
```

contrary to (32).  Thus maximum-basis exchange closes this cap-equality row
even without exhibiting (35).

## 6. The corrected carrier and remaining frontier

The rank-eight route now has three genuinely different layers:

```text
tau>5/49:   immediate Hunter contradiction;
tau=5/49:   intersect the connected-set systems of all maximum trees;
tau<5/49:   search finite-ring needles and then higher intersections.       (39)
```

Signed affine rank and quotient parity are useful sidecars, but neither can
replace the explicit torsion search.  Equation (39), not the binary conjecture
falsified in THM-2103, is the lawful next target.

### Assumption challenge and Tournament Analysis

Candidate vertices considered here are terminal characters, primitive
relation types, threshold edges, kernel submodules, torsion points, maximum
graphic-matroid bases, and Hunter connectivity obligations.  The pairwise
observable is (1); thresholding it produces the toothpick graph and its higher
layers.  One may orient edges by weight and use a label tie-break to obtain a
tournament scheduler, but score histograms, SCCs, cycles, edge flips, and
Hamiltonian paths do not recover either the missing kernel point (17) or the
all-bases condition (37).  The faithful carrier is

```text
exact relation-labelled edge weights
 + graphic-matroid maximum bases
 + finite-ring kernel cover
 + higher active-set intersections.                    (40)
```

The challenged assumption is that the continuous two-torus obstruction is
fundamentally continuous.  Sections 3 and 5 show that its first hidden layer
is often a finite Kakeya/blocking problem, while Section 2 shows that its
smallest pair relations form a dyadic self-similar graph.  These sharpen the
rank-eight object but do not yet discharge every signed-pencil-free,
finite-ring-blocked coefficient plane.  LRC(14) remains open. QED.
