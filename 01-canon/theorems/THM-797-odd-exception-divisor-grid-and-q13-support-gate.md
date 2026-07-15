---
id: THM-797
title: Odd-exception divisor grids, the q=13 folded-support gate, and the global-component correction
status: PROVED (general odd-grid residue-shell lemma and q=13 support corollaries) + VERIFIED (sharp aligned survivor, exact deep components, and global erosion witness)
source: codex-2026-07-14-S10 global-erosion selection analysis
depends_on:
  - THM-772   # two-sheet packet arithmetic and mandatory deep exception
  - THM-774   # folded-diamond equivalence
  - THM-789   # global erosion target
related: [THM-769, THM-775, HYP-6800, HYP-6820]
verification:
  - 04-computation/lrc13_odd_exception_divisor_grid_codex_S10.py
  - 05-knowledge/results/lrc13_odd_exception_divisor_grid_codex_S10.out
---

# THM-797 — Odd-exception divisor grids and the q=13 support gate

Use THM-772's two-sheet notation

```text
A=2U union {x,y},       |U|=10,       x,y odd,
E_U={t:phi_U(t)>=1/11},
R_U={d:max_(u in U)||ud||<2/143}.
```

Put

```text
a=(x+y)/2,              b=|x-y|/2,
H_(x,y)={t:||at||+||bt||>=11/13}.
```

Thus THM-774 identifies `H_(x,y)` with the two-exception eligibility and
opposite-sheet-colour predicate, and THM-789's residual target is

```text
E_U not subset H_(x,y) minus R_U.                         (0)
```

## 1. The general odd-grid theorem

Let `q>=3` be odd.  Write `[z]_q` for the unique balanced residue in

```text
{-(q-1)/2,...,(q-1)/2}.
```

For a unit numerator `p mod q`, considered modulo sign, define

```text
D_q(U)={p mod +/- :
        |[up]_q|>=ceil(q/11) for every u in U},           (1)

A_q(x,y)={p mod +/- :
          |[xp]_q|,|[yp]_q|<=floor(2q/13),
          [xp]_q != [yp]_q (mod 2)}.                      (2)
```

Both predicates are invariant under `p->-p`.

### Theorem 1 — deep grids and exception shells

For every unit `p mod q`,

```text
p/q in E_U       iff p in D_q(U),                         (3)
p/q in H_(x,y)   iff p in A_q(x,y).                       (4)
```

Consequently,

```text
D_q(U) not subset A_q(x,y)
  implies E_U not subset H_(x,y) minus R_U.               (5)
```

If `q|x`, then (2) simplifies to the explicit shell

```text
A_q(x,y)={p mod +/- :
          q does not divide y,
          1<=|[yp]_q|<=floor(2q/13),
          [yp]_q is odd}.                                 (6)
```

In particular, if `q` divides both exceptions, this shell is empty.

### Proof

For (3), simply observe that

```text
||up/q||=|[up]_q|/q.
```

Because the numerator is integral, the inequality `>=1/11` is exactly the
threshold in (1).

For (4), write

```text
xp=qN_x+r_x,       yp=qN_y+r_y,
r_x=[xp]_q,        r_y=[yp]_q.                            (7)
```

Since `q` is odd, the balanced residues have magnitude strictly less than
`q/2`, so `N_x,N_y` are the unique nearest integers.  Closed exception
eligibility is

```text
||xp/q||,||yp/q||<=2/13,
```

equivalently

```text
|r_x|,|r_y|<=floor(2q/13).                               (8)
```

Modulo two, `q,x,y` are all one.  Formula (7) therefore gives

```text
N_x=p-r_x,          N_y=p-r_y          (mod 2).           (9)
```

The two sheet colours are opposite exactly when `N_x,N_y` have opposite
parity, which by (9) is exactly the last condition in (2).  THM-774 identifies
this lifted predicate with the folded diamond, proving (4).

Now `0 in R_U`, so morphological erosion satisfies

```text
H_(x,y) minus R_U subset H_(x,y).                         (10)
```

A class in `D_q(U)\A_q(x,y)` supplies through (3)--(4) a point of
`E_U\H_(x,y)`, and (10) proves (5).

Finally suppose `q|x`.  Then `r_x=0`.  Also `x/q` is odd, and the exact
nearest integer `N_x=(x/q)p` has parity `p`.  The two colours are opposite
precisely when `r_y=[yp]_q` is odd.  Eligibility supplies the magnitude bound;
an odd residue is nonzero, so `q` cannot divide `y`.  This proves (6). ∎

## 2. The mandatory q=13 support gate

Now assume the arithmetic conclusions of THM-772 for a hypothetical primitive
tight two-sheet packet.  In particular, no member of `U` is divisible by 13,
and at least one deep-branch exception is.

Let

```text
C=(Z/13Z)^*/{+/-1},
S(U)={+/-u mod 13:u in U} subset C.                       (11)
```

Inversion acts on these six folded classes.  Since

```text
ceil(13/11)=2,
```

condition (1) says that no `up` is `+/-1`.  Hence

```text
D_13(U)=C minus S(U)^(-1).                                (12)
```

Also `floor(26/13)=2`.  If `13|x` but `13` does not divide `y`, the only odd
nonzero balanced residues in the shell (6) are `+/-1`, so

```text
A_13(x,y)={+/-y^(-1)}.                                   (13)
```

Therefore hypothetical containment in (0), or even the weaker containment
`E_U subset H_(x,y)`, forces

```text
C minus S(U)^(-1) subset {+/-y^(-1)},
```

and inversion gives the exact support gate

```text
C minus S(U) subset {+/-y}.                               (14)
```

Thus

```text
S(U)=C
or
S(U)=C minus {+/-y}.                                     (15)
```

If both `x` and `y` are divisible by 13, (6) is empty.  Containment forces
`D_13(U)` empty, hence

```text
S(U)=C.                                                   (16)
```

This proves the following uniform eliminations before any height search:

1. every core supported on at most four folded mod-13 classes;
2. every five-class core whose missing class is not the folded class of the
   off-divisor exception; and
3. every non-full support core when both exceptions are 13-divisible.

The only q=13-invisible cases are full six-class support and the **aligned
five-class** pattern in (15).

## 3. Sharpness of the residue gate and the global-component correction

The aligned pattern really can trap the entire prime grid.  Take

```text
U_0=(1,2,3,4,7,9,10,11,12,16),
(x,y)=(13,5).                                             (17)
```

This core is primitive, contains a multiple of every `m=2,...,12`, has no
13-multiple, and has `B=16`.  Exact breakpoint enumeration gives

```text
M(U_0)=2/13,
argmax phi_(U_0)={5/13,8/13}.                             (18)
```

The breakpoint certificate is finite and lossless: `phi_U` is piecewise
linear, and all its possible maxima occur at a runner cusp `k/(2u)` or an
intersection of two signed linear pieces, whose denominators divide `u+v` or
`|u-v|`.  The verifier enumerates exactly these rationals and evaluates them
with `Fraction` arithmetic.

The row also satisfies every earlier two-sheet arithmetic tax:

```text
x,y<=11B,             min(x,y)<=11B-36,
rho=(M(U_0)-1/13)/B=1/208,
1/(xy)+2rho=1/40 <= 36/845
             =2/(13x)+2/(13y).                           (19)
```

The folded support multiplicities on the classes `1,...,6` are

```text
(2,2,3,2,0,1).                                           (20)
```

Thus `S(U_0)=C\{5}` and the missing class is exactly the class of `y`.
The only deep thirteenth numerators are `p=5,8`; at both,

```text
phi_(U_0)(p/13)=2/13,
Q_(9,4)(p/13)=12/13.                                     (21)
```

So every point visible to the q=13 gate lies inside the folded diamond.

Nevertheless global escape is immediate at another component endpoint:

```text
t=7/33,
phi_(U_0)(t)=1/11,
Q_(9,4)(t)=8/33<11/13.                                   (22)
```

Thus `t in E_(U_0)\H_(13,5)`, and since `0 in R_(U_0)`, (22) proves the full
erosion failure (0).

The complete closed `1/11`-superlevel set has eight components:

```text
[12/77,7/44],       [23/110,7/33],
[67/176,43/110],    {9/22},
{13/22},             [67/110,109/176],
[26/33,87/110],      [37/44,65/77].                       (23)
```

Their endpoint-owner pairs are

```text
7->12, 10->9, 16->10, {10,12},
{10,12}, 10->16, 9->10, 12->7.                           (24)
```

They form four reflection orbits.  The orbit containing the global maxima
`5/13,8/13` is entirely trapped: its minimum `Q` is `159/176>11/13`.
The second orbit contains (22) and escapes with a large margin.  Therefore
**selecting only global maximizers is insufficient**.  The faithful global
state must retain every closed `1/11`-superlevel component, including
threshold endpoints and endpoint owners.

## 4. What remains after the support gate

THM-789 showed why refining one trapped deep anchor cannot prove global
noncontainment.  Theorem 1 supplies a different global selector whenever an
odd exception divisor has a deep unit class outside the other exception's
explicit shell.  The q=13 specialization eliminates most residue supports,
but (17)--(21) prove the gate is sharp.

The remaining two-sheet program should therefore branch on two axes:

1. **arithmetic divisor grids:** test every useful odd divisor of an exception
   through (1)--(6), with q=13 always first; and
2. **geometric deep components:** when a grid is silent, enumerate or bound
   the owner-labelled components of `E_U` and compare their minimum folded-
   diamond margins, rather than refining only an argmax component.

For the mandatory q=13 branch, only full six-class support and aligned
five-class support need this second stage.  No claim of complete two-sheet or
uniform `n=12` closure is made here.

## 5. Tournament Analysis and challenged assumptions

There are two useful but deliberately decorated tournament views.

### Folded residue obligations

Take the six folded classes in `C` as vertices.  For (17), orient by core
multiplicity, breaking ties by the numerical class order; use inversion

```text
(2 6)(3 4)
```

as the switch/gauge and the sorted order as the tie Hamiltonian path.  Before
and after the switch the paths are

```text
(5,6,1,2,4,3),          (5,2,1,3,6,4).
```

The gauges differ on five edges.  Both tournaments are transitive: score
histogram `0,...,5`, zero directed triangles, singleton SCCs, and one
Hamiltonian path.  The bare fingerprints do not expose the missing class or
the exception shell; the decorated incidence `(class,multiplicity,inverse,
accepted-by-y)` does.

### Deep-component alternatives

Take the four reflection orbits in (23) as vertices.  The first pairwise
observable is escape margin

```text
11/13-min_(t in C) Q_(9,4)(t),
```

and the switch/gauge replaces it by component width.  With index tie-breaking,
the two Hamiltonian paths are

```text
(2,0,3,1),             (3,1,0,2),
```

again with five edge flips.  Both tournaments are transitive with scores
`0,...,3`, zero cycles, singleton SCCs, and one Hamiltonian path.  Yet only
the decorated first gauge remembers the **sign** of escape and the endpoint
owners.

Runner vertices lose the divisor grid; residue vertices lose component
geometry; a single global-max vertex loses alternative components; bare
tournaments lose all metric signs.  The minimal carrier suggested by this
theorem is a bipartite incidence object between folded divisor obligations and
owner-labelled deep components, with exact escape margins as a sidecar.  This
is the assumption challenge: the underlying object is a family of global
alternatives, not one privileged deepest point.

## Exact replay

The verifier exhausts every odd `q=3,...,101`, every unit numerator, and every
odd pair of speed classes modulo `2q`: `10,971,770` exact tests of the direct
lifted predicate, the balanced-residue shell, and the folded diamond, with no
failures.  It separately checks `143,346` divisor specializations, all `832`
q=13 support/exception profiles, the exact breakpoint certificate (18), every
component in (23), the witness (22), and both tournament fingerprints.  No
floating point or sampled-circle verdict enters the proof.
