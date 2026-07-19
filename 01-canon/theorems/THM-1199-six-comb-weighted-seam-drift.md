---
id: THM-1199
title: Six-comb weighted seam drift and unconditional last-tooth descent
status: PROVED analytic composition / COMPUTER-EXACT arithmetic audit.  Every six-comb cover pays a c/16 weighted quantum on every edge of every connected-nerve spanning tree, yielding a phase-free twelve-piece functional H-drift.  The five-prefix dual survivor has physical length at least 1/(49c), so the beat mesh gives the unconditional bound (d_6+d_5)/c<49C_bar with no H_5<2 hypothesis.  The valid local gcd refinement uses the actual endpoint owners; a naive gcd(d_i,d_j) endpoint factor is refuted.  This is not universal six-comb noncoverage or LRC(14)
source: codex-2026-07-18 weighted-functional-drift session
depends_on: [THM-1178, THM-1196, THM-1198]
related: [THM-1094, THM-1176, THM-1192]
script: 04-computation/lrc14_six_comb_weighted_seam_drift_thm1199.py
output: 05-knowledge/results/lrc14_six_comb_weighted_seam_drift_thm1199.out
---

# THM-1199 -- six-comb weighted seam drift and unconditional last-tooth descent

At radius `1/14`, let

```text
G=[(14k+1)/(14c),(14k+13)/(14c)]                    (1)
```

be a complete safe gap of the integer speed `c`, and suppose the six distinct
faster integer speeds

```text
c<d_1<...<d_6                                           (2)
```

have strict danger combs covering `G`.  THM-1198 supplies the symmetric
six-bin probability density

```text
f=(3/4,13/12,7/6,7/6,13/12,3/4)                    (3)
```

on the six equal subintervals of `[0,1]`.  In particular,

```text
integral f=1,              3/4<=f<=7/6.             (4)
```

This theorem extracts two consequences that neither the unweighted seam
surplus nor the one-comb operator norm sees alone.

## 1. Normalization and the strict/closed endpoint audit

Write the affine slow-gap coordinate as

```text
t=(14k+1)/(14c)+(6/(7c))x,             0<=x<=1,     (5)
L_i=6d_i/(7c),
phi_i=d_i(14k+1)/(14c) mod 1.                         (6)
```

Let `U_i` be the preimage in `[0,1]` of the **strict** danger comb
`||d_i t||<1/14`, and put

```text
A_i=integral_[0,1] f(x)1_(U_i)(x)dx.                 (7)
```

THM-1198 states its operator norm using closed teeth.  This causes no gap:
inside a compact interval, strict and closed teeth differ at finitely many
points, so their `f dx` loads are equal.  Using closed teeth only enlarges
the set in the upper bound, and therefore

```text
A_i<=Pbar(L_i)<=7/36,                                (8)
```

where `Pbar` is THM-1198's exact twelve-piece phase envelope on
`6/7<=L<=3` and its BV majorant `1/7+1/(7L)` for `L>3`.

Every `U_i` is nonempty.  Otherwise five strict combs cover `G`, hence their
closed enlargements cover it, contrary to THM-1198.  The intersection graph
of the six relatively open sets `U_i` is connected: a disconnected vertex
partition would split the connected interval `[0,1]` into two nonempty,
disjoint, relatively open unions.  Thus every spanning tree used below is a
genuine tree of positive-length handoffs.  This is why the proof takes the
nerve of the strict teeth before invoking the closed-tooth load bound.

## 2. Weighted nerve-tree seam theorem

Fix any spanning tree `T` of the connected intersection graph.  For every
edge `e={i,j}` choose one connected component `J_e` of `U_i intersect U_j`
in the physical `t` coordinate.  Its two endpoint owners

```text
u_e,v_e in {c,d_i,d_j}                               (9)
```

record whether the relevant endpoint is a slow-gap wall or a tooth wall.
The rational endpoint calculation from THM-1178 gives

```text
|J_e|>=gcd(u_e,v_e)/(14u_e v_e).                    (10)
```

Scaling (10) by `7c/6` and then using the lower density `f>=3/4` gives the
new weighted quantum

```text
integral_(J_e) f(x(t)) (7c/6)dt
 >=(7c/6)(3/4)gcd(u_e,v_e)/(14u_e v_e)
 =c gcd(u_e,v_e)/(16u_e v_e).                       (11)
```

At each covered point, if `N(x)=#{i:x in U_i}`, the tree induced by those
`N(x)` vertices has at most `N(x)-1` edges.  Multiply this pointwise fact by
`f(x)` and integrate.  Equations (7) and (11) prove:

> **Theorem 1 (weighted seam drift).**  For every chosen nerve spanning tree
> and every choice of one positive intersection component on each edge,
>
> ```text
> sum_(i=1)^6 A_i-1
>  >=(c/16)sum_(e in T) gcd(u_e,v_e)/(u_e v_e)
>  >=(c/16)sum_({i,j} in T)1/(d_i d_j).              (12)
> ```

The last inequality is phase-free.  Indeed
`gcd(u,v)/(uv)=1/lcm(u,v)`, and the lcm of any two actual owners in
`{c,d_i,d_j}` is at most `d_i d_j`.

The owner labels in (12) are essential.  The local endpoint calculation does
**not** permit replacing the edge term uniformly by
`gcd(d_i,d_j)/(d_i d_j)`.  For

```text
(c,d_i,d_j,k)=(7,10,34,2),                           (13)
```

the pair intersection has exactly one component.  Its endpoint owners are
`(7,34)` and its length is

```text
1/3332=gcd(7,34)/(14*7*34)
       <gcd(10,34)/(14*10*34)=1/2380.               (14)
```

Thus (10) is the sharp local gcd-refined statement; the unlabeled
`1/(d_i d_j)` floor is the correct universal projection.  A stronger global
cover theorem involving `gcd(d_i,d_j)` would require an additional argument
that forbids boundary-owned handoffs; (14) shows it cannot come from endpoint
arithmetic alone.

**Later chronological upgrade.**  THM-1253 supplies exactly that additional
argument for one deletion-minimal cover by *individual* teeth.  Its raw
consecutive handoffs have fast-owned endpoints and are all pairwise disjoint,
so the weighted coverage excess charges the entire owner word:

```text
sum_i Pbar(6d_i/(7c))-1
 >=(c/16)sum_(a=1)^(N-1)1/lcm(s_a,s_(a+1)).          (14a)
```

Thus (12) remains the phase-free arbitrary-nerve theorem, while (14a) is the
stronger phase-located full-word consumer now available on the canonical
chronological carrier.

## 3. The phase-free functional H-drift

Combining (8) and (12) gives the phase-free envelope

```text
F_6(d_1,...,d_6;c)
 :=sum_(i=1)^6 Pbar(6d_i/(7c))-1
 >=(c/16)sum_({i,j} in T)1/(d_i d_j).                (15)
```

This is a genuine functional form of the H-drift: unlike the scalar harmonic
pressure, its left side remembers twelve exact ratio chambers and the
oscillatory one-comb load.  Rooting `T` at the largest speed gives

```text
F_6>=c/(16d_6) sum_(i=1)^5 1/d_i.                   (16)
```

In the normalized reciprocal coordinates `x_i=c/d_i`, equations (15)--(16)
become

```text
sum_i Pbar(6/(7x_i))-1
 >=x_6/(16c) sum_(i=1)^5 x_i.                       (17)
```

For comparison, the simultaneous unweighted THM-1178 law is

```text
H_6-1>=7x_6/(12c) sum_(i=1)^5 x_i,
H_6=sum_i x_i.                                      (18)
```

The two inequalities constrain different functionals on the same ordered
reciprocal packet.  They must not be added or identified.  Also, since
`A_i<=7/36`, (12) supplies the finite seam budget

```text
c sum_({i,j} in T)1/(d_i d_j)<=8/3.                 (19)
```

Both (17) and (18) retain the unavoidable `1/c` toothpick drift.  Neither is
a scale-invariant contradiction.

## 4. A quantitative survivor after any five combs

There is a second, more consequential use of both extrema in (4).  Let `E_x`
be the complement in `[0,1]` of any five normalized strict danger combs.
Their union load is at most the sum of their loads, so THM-1198 gives

```text
integral_(E_x) f
 >=1-sum_(i=1)^5 A_i
 >=1-5(7/36)=1/36.                                  (20)
```

Since `f<=7/6`,

```text
|E_x|>=1/42.                                         (21)
```

Returning through (5) yields the uniform physical survivor law

```text
|E_t|>=1/(49c).                                      (22)
```

More sharply, define the positive phase-free deficit

```text
eta=1-sum_(i=1)^5 Pbar(6d_i/(7c))>=1/36.             (23)
```

Then the same proof gives

```text
|E_x|>=6eta/7,             |E_t|>=36eta/(49c).      (24)
```

This upgrades THM-1198's qualitative dual noncoverage to a uniform interval-
measure input for the recursive beat descent.

## 5. Unconditional last-tooth descent

Now take the first five combs in the putative six-cover as the prefix and let
`E_t` be their physical survivor in `G`.  For `i<=5`, let `n_i` count the
open `d_i` teeth meeting `G`, and put

```text
C=1+sum_(i=1)^5 n_i,
C_bar=1+sum_(i=1)^5 ceil((6d_i+c)/(7c)).             (25)
```

THM-1196 proves that `E_t` has at most `C` components and avoids every sum
and difference beat lattice joining `d_6` to the prefix.  The lattice of
denominator `d_6+d_5` is among them, so every survivor component has length
strictly less than `1/(d_6+d_5)`.  Hence

```text
|E_t|<C/(d_6+d_5).                                   (26)
```

Combining (22) or (24) with (26) proves:

> **Theorem 2 (unconditional last-tooth bound).**  Every six-comb cover of a
> complete `c`-slow gap satisfies
>
> ```text
> (d_6+d_5)/c<49C<=49C_bar,                          (27)
> ```
>
> and the phase-free functional refinement
>
> ```text
> (d_6+d_5)/c<49C/(36eta)<=49C_bar/(36eta).          (28)
> ```

Unlike THM-1196's earlier harmonic-to-Farey estimate, (27) has no
`H_5<2` hypothesis and no denominator that degenerates on the prefix-critical
face.  It proves that **every fixed five-speed prefix has a finite final-speed
range**.  The bound still grows with the prefix tooth count; it does not make
the first five ratios compact.

## 6. Kakeya carrier and Tournament Analysis audit

The protected slow gap is the Kakeya needle, but one vertex set does not
faithfully carry both new arguments.  The weighted drift uses

```text
endpoint-labelled positive handoff components on a connected nerve tree,
```

while the descent uses

```text
phase-local Farey cells separated by forbidden beat points.
```

We explicitly challenged runners, gaps, fixed sections, section boundaries,
wall-crossing events, residues, overlap components, beat points, and proof
obligations as possible vertices.  Speed labels preserve only the ordered
corollary (16) and the fact that `d_6+d_5` is the largest sum denominator.
They destroy the endpoint owners in (10), the gap address `k`, the beat
numerators, and the survivor-cell occupancy predicate.

For the mandatory pairwise audit, take observable `d_j-d_i`; its sign is the
switch/gauge, and increasing speed breaks no ties because (2) is strict.  The
resulting tournament is transitive, with score histogram
`0,1,2,3,4,5`, no directed cycle, six singleton SCCs, and the unique tie
Hamiltonian path

```text
d_1 -> d_2 -> d_3 -> d_4 -> d_5 -> d_6.             (29)
```

This tournament retains denominator order but not the predicates proved in
(12) and (26).  The correct structural object is therefore a two-layer
certificate: an endpoint-labelled nerve tree laid over a phase-local beat
mesh, not a tournament on runners alone.

## 7. Exact replay and honest frontier

The dependency-free `Fraction` referee checks the density extrema and every
new scale factor, then scans `44,655` phase-labelled pair rows containing
`29,750` positive overlap components.  Every component is an integer multiple
of its owner-labelled quantum, and the phase-free ratio attains exactly one.
It also freezes the counterexample (13)--(14), and checks the tooth-count
ceiling on `29,400` phase rows.  Normal and optimized runs are byte-identical
to the stored output.  Frozen SHA-256 hashes are

```text
source  ae58b171908ba0c83ade2c4c1a2c600b91162bcbba7d228d09ded7ed9527f0d5
output  e36312f013a66175a3f69ca53ab18bb8158b4dd70158f76a13f9261d8d12e6ea
```

Theorem 1 gives a new phase-sensitive functional drift and Theorem 2 removes
the last harmonic hypothesis from the final-tooth descent.  They do not rule
out the remaining compact or prefix-escaping five-speed families, produce a
uniform all-prefix bound, prove universal six-comb noncoverage, or prove
LRC(14).  Those are the honest residuals.
