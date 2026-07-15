---
id: THM-789
title: Two-sheet erosion, global gap/Kneser budgets, and the symmetric return-packet upgrade; exact local and component-tournament obstructions
status: PROVED (pointwise and global metric erosion + Kneser/Macbeath) + VERIFIED (exact local trap, component-erosion tournament liar, and signed-slope liar)
source: codex-2026-07-14-S10 phase-noncontainment analysis
depends_on:
  - THM-774   # folded-diamond equivalence
  - THM-780   # quantitative phase pigeonhole
  - THM-782   # one-sided return packet, strengthened here
  - LRC(<=13)
related: [THM-772, THM-775, THM-776, THM-797, HYP-6820]
verification:
  - 04-computation/lrc14_two_sheet_erosion_trap_codex_S10.py
  - 05-knowledge/results/lrc14_two_sheet_erosion_trap_codex_S10.out
  - 04-computation/lrc14_two_sheet_erosion_budget_liar_codex_S10.py
  - 05-knowledge/results/lrc14_two_sheet_erosion_budget_liar_codex_S10.out
---

# THM-789 — erosion and the symmetric return packet

Let `U` be ten distinct positive integers and put

`B=max(U)`, `alpha=1/13`, `beta=1/11`, `gamma=beta-alpha=2/143`,

`phi_U(t)=min_(u in U)||ut||`,

`G_U={t:phi_U(t)>alpha}`,

`E_U={t:phi_U(t)>=beta}`,

`R_U={d:max_(u in U)||ud||<gamma}`.

For odd exceptions `x,y`, let `H_(x,y)` denote THM-774's folded diamond.

## 1. Bohr erosion and the pointwise thickness tax

### Theorem 1

`E_U+R_U subset G_U`.                                             (1)

Consequently, if `2U union {x,y}` is tight at `1/13`, then

`E_U subset H_(x,y) minus R_U`,                                   (2)

where the morphological erosion is

`H minus R={t:t+R subset H}`.

More strongly, tightness forces, for every `t in G_U` and each
`w in {x,y}`,

`||wt|| + (w/B)(phi_U(t)-1/13) <= 2/13`.                           (3)

At every `1/11`-deep time this specializes to

`||wt|| <= 2(11B-w)/(143B)`.                                      (4)

In particular `w<=11B`. If `w=11B`, then `wt` is integral at every
`t in E_U`, and necessarily `M(U)=1/11`.

### Proof

For `t in E_U` and `d in R_U`, the circle triangle inequality gives

`||u(t+d)|| >= ||ut||-||ud|| > beta-gamma=alpha`

for every `u`, proving (1). Tightness is exactly `G_U subset H_(x,y)` by
THM-774, hence (2).

For the pointwise statement, put

`r=(phi_U(t)-alpha)/B`.

The centered open interval of radius `r` about `t` lies in `G_U`: for
`|d|<r`,

`||u(t+d)|| >= phi_U(t)-u|d| > alpha`.

Under tightness it lies in `H_(x,y)`. Membership in the folded diamond implies
the individual eligibility condition `||ws||<=2/13`. A connected interval
contained in this union of disjoint w-teeth lies in the tooth containing `t`.
Its radius is at most the distance to that tooth's nearer boundary, so

`wr <= 2/13-||wt||`,

which is (3). If `t in E_U`, then `r>=gamma/B`, giving (4). Its right side
must be nonnegative, hence `w<=11B`. Equality forces `||wt||=0` for all deep
times. Settled `LRC(<=13)`, applied to these ten speeds, makes `E_U` nonempty.
If `M(U)>1/11`, continuity
would make `E_U` contain an interval, impossible inside the finite w-grid;
therefore `M(U)=1/11`. ∎

### Theorem 1A — the global gap and return-mass budgets

Put

```text
delta=gamma/B=2/(143B),       I_delta=(-delta,delta).
```

Let the cyclic complementary gaps between the connected components of
`E_U` have lengths `g_1,...,g_r`, including the wraparound gap.  Hypothetical
tightness forces

```text
mu(E_U)+sum_i min(g_i,2delta) <= mu(H_(x,y)) <= 8/117,   (4a)
mu(E_U)+mu(R_U)                <= mu(H_(x,y)),            (4b)
mu(R_U) >= max(2delta,2*72^(-10)).                        (4c)
```

In particular,

```text
mu(E_U) <= mu(H_(x,y))-max(4/(143B),2*72^(-10)).         (4d)
```

### Proof

If `|d|<delta`, then for every `u in U`,

```text
||ud|| <= u|d| < B delta=gamma,
```

so `I_delta subset R_U`.  By (1)--(2), tightness gives

```text
E_U+I_delta subset E_U+R_U subset H_(x,y).
```

Thickening a finite union of closed circle intervals by a centered interval
of radius `delta` removes exactly `min(g_i,2delta)` from each complementary
gap.  Hence

```text
mu(E_U+I_delta)=mu(E_U)+sum_i min(g_i,2delta),
```

which proves (4a).  Macbeath's circle sumset inequality gives

```text
mu(E_U+R_U)>=min(1,mu(E_U)+mu(R_U)).
```

Since this sumset lies in `H_(x,y)` and `mu(H_(x,y))<1`, (4b) follows.  The
interval inclusion gives `mu(R_U)>=2delta`.  The difference packet `D` in
Theorem 2 satisfies `||ud||<1/72<2/143=gamma` for every `d in D`, hence
`D subset R_U`; (5) below gives the second lower bound in (4c).  Finally use
THM-774's sharp diamond cap `mu(H)<=8/117`. ∎

## 2. Symmetric return-packet upgrade

### Theorem 2

Choose a `1/11`-deep time `t_0` and a heavy `72^10` joint phase-cell preimage
`A` as in THM-782. Put `D=A-A`. Then

`D=-D`, `0 in D`, `t_0+D subset G_U`,

`mu(D)>=2mu(A)>=2*72^(-10)`.                                    (5)

Consequently

`mu(G_U)>=2*72^(-10)`,                                          (6)

and `G_U` has a component of length at least

`2*72^(-10)/(sum_(u in U)u) >= 72^(-10)/(5B)`.                  (7)

This improves THM-782's recorded measure floor by a factor of two and its
component-width floor by a factor of four. Under hypothetical tightness,

`D subset (H_(x,y)-t_0) intersect (t_0-H_(x,y))`.                (8)

### Proof

Two points of `A` lie in the same length-`1/72` interval in every U-coordinate,
so every `d in D` obeys `||ud||<1/72`. The same clearance calculation as in
THM-782 gives `t_0+D subset G_U`; symmetry and `0 in D` are immediate.

The set `A` lies inside the preimage of any one coordinate cell, whose measure
is `1/72`, so `mu(A)<=1/72<1/2`. The Kneser-Macbeath difference-set inequality
on the circle therefore gives

`mu(A-A)>=min(1,2mu(A))=2mu(A)`,

proving (5)--(6).

For a fixed `u`, `{t:||ut||<=1/13}` is a union of `u` closed arcs. The union of
all ten danger sets has at most `sum u` components, so its open complement
`G_U` also has at most `sum u` components. Averaging (6), then using
`sum u<=10B`, proves (7). Under tightness, `t_0+D subset H`; applying both `d`
and `-d` gives (8). ∎

## 3. Exact obstruction: the full natural return set can be locally trapped

Take

`U_0={1,2,3,5,7,8,9,10,11,12}`, `t_0=4/17`, `(x,y)=(13,9)`.

This core is primitive, hereditarily primitive, divisor-complete through 12,
contains no multiple of 13, and has

`M(U_0)=2/17`.

At `t_0` its ten clearance numerators over 17 are

`(4,8,5,3,6,2,2,6,7,3)`,

so it is `1/11`-deep. The nearest integers to `13t_0` and `9t_0` are 3 and 2,
of opposite parity, and the determinant is already sharp:

`13*2-9*3=-1`.

Writing `a=(13+9)/2=11`, `b=(13-9)/2=2`, one has

`Q(t_0)=||11t_0||+||2t_0||=15/17`,

and the whole interval

`t_0+[-8/1989,8/1989] subset H_(13,9)`.                         (9)

On this interval the signed 11- and 2-errors keep their signs, so
`Q(t_0+d)=15/17-9d`; its minimum occurs at the positive endpoint and equals
`11/13`, proving (9) including the endpoints.

The full return Bohr set is exactly

`R_(U_0)=(-1/858,1/858)`.                                      (10)

Indeed coordinate 1 first restricts the circle representative to
`|d|<2/143`, and coordinate 12 then forces `|d|<(2/143)/12=1/858`; all other
coordinates satisfy the converse inclusion. Since `1/858<8/1989`, (9)--(10)
give

`t_0+R_(U_0) subset H_(13,9)`.                                 (11)

Likewise every `72`-grid U-cell `C` obeys

`C-C subset (-1/864,1/864)`,                                  (12)

because coordinates 1 and 12 force it. Thus every anchored same-cell packet,
its full symmetric difference packet, and every further literal refinement by
the `x,y,a,b` phases remains trapped at this base time. The ordinary Lipschitz
safe interval has radius

`(2/17-1/13)/12=3/884<8/1989`

and is trapped as well.

This is not a counterexample to global escape. For the same triple,

`tau=14/19`

is `2/19`-deep while

`||11tau||+||2tau||=11/19<11/13`,

so `tau` lies outside the diamond. It proves a precise method boundary: local
symmetrization, added phase coordinates, and return-width arguments at a fixed
deep time cannot establish structured noncontainment. The sharpened residual is

`E_U not subset H_(x,y) minus R_U`;                              (13)

the proof must globally select a different deep time or component.

## 4. Exact obstruction: the raw component tournament loses erosion incidence

The global budgets are stronger than a local return packet, but their scalar
and raw-ranking shadows still do not decide (13).  For the same core `U_0`,
the `1/11`-deep set has six components

```text
C1=[12/77,7/44]       C2=[23/99,21/88]
C3=[23/88,32/121]     C4=[89/121,65/88]
C5=[67/88,76/99]      C6=[37/44,65/77].
```

For an odd pair define the raw and eroded escape margins

```text
m(C)=sup_(t in C)(11/13-Q(t)),
m_R(C)=sup_(t in C+R_(U_0))(11/13-Q(t)).                       (14)
```

Positive sign means that the component, respectively its full return
thickening, escapes the folded diamond.  Exact evaluation for two odd pairs
gives:

| pair | raw margins on `C1,...,C6` | eroded margins on `C1,...,C6` |
|---|---|---|
| `(13,9)` | `(159/572,-7/1144,447/1573,447/1573,-7/1144,159/572)` | `(15/52,5/1144,2825/9438,2825/9438,5/1144,15/52)` |
| `(17,13)` | `(197/1001,-59/1144,538/1573,538/1573,-59/1144,197/1001)` | `(1301/6006,-125/3432,3415/9438,3415/9438,-125/3432,1301/6006)` |

Orient components toward larger raw margin and break ties by left endpoint.
Both rows then have the identical transitive tournament order

```text
C2 < C5 < C1 < C6 < C3 < C4,
```

the identical raw sign decoration (only `C2,C5` are negative), score
histogram `(0,1,2,3,4,5)`, no directed cycles, six singleton SCCs, and one
Hamiltonian path.  Both folded diamonds have four components and

```text
mu(H)=8/169,        mu(H minus R_(U_0))=212/5577.               (15)
```

Nevertheless their erosion incidence differs.  The closed thickening of
`C2` is

```text
[595/2574,823/3432].                                           (16)
```

For `(13,9)` the relevant diamond tooth is `[37/169,28/117]`, and (16)
overshoots its right endpoint by `5/10296`; hence `m_R(C2)=5/1144>0`.
For `(17,13)` the relevant tooth is `[50/221,41/169]`, and (16) lies inside
with left/right slacks `215/43758` and `125/44616`; hence
`m_R(C2)=-125/3432<0`.  Reflection gives the same discrepancy on `C5`.

There is a pointwise version of the same loss.  At `t_0=4/17`, the pairs
`(13,9)` and `(43,13)` have the same folded value `Q=15/17`, the same
unsigned odd-error multiset `{1/17,2/17}`, opposite parity, and the same sharp
determinant `-1`.  The first pair traps all of `t_0+R_(U_0)` by (11), whereas
for `d=1/1000 in R_(U_0)`,

```text
Q_(43,13)(t_0+d)-11/13=-1503/221000.
```

Thus the signed affine tooth slope, not the unsigned errors or determinant,
decides the return incidence.  The faithful component vertex must retain its
tooth address, signed slope, and return-set incidence—or directly `m_R(C)`.
Raw component order, raw signs, and the scalar eroded measure in (15) are
provably insufficient.

## 5. Tournament and assumption challenge

Runner vertices, gap vertices, and fixed-time phase-cell vertices all lose the
global choice exposed by (11). A more faithful candidate vertex set is the set
of connected components of `E_U`. Give component `C` the escape observable

`m(C)=sup_(t in C)(11/13-Q_(a,b)(t))`,

orient pairs toward larger `m`, and use left endpoint as a tie gauge. The tie
Hamiltonian path is the sorted component list. This tournament is transitive
(score histogram `{0,...,r-1}`, no directed cycles, singleton SCCs, one
Hamiltonian path); its value is not nontransitive structure but preservation of
the **global alternatives**. In the exact example the escaping `14/19` choice
beats the locally trapped `4/17` choice.  Section 4 now shows the exact limit
of that proposal: even the raw sign-decorated component tournament and the
total eroded-diamond measure can agree while componentwise erosion incidence
differs.

Even this quotient must retain component intervals, owner labels, signed
folded-tooth addresses, and the eroded pointwise margin; its bare isomorphism
class destroys the LRC predicate.  The challenged assumption is that refining
one successful local phase anchor adds information.  It is monotone trapping
here.  Recursion must branch across deep components, then transport the full
return incidence on each branch, not only descend inside one phase cell.

THM-797 implements the first arithmetic branch of this prescription.  An odd
exception-divisor grid class outside its explicit acceptance shell gives an
immediate global escape; at `q=13` only full and aligned-five folded supports
survive.  Its sharp aligned row traps every global maximizer but escapes at a
different threshold component, confirming that the branch set in this section
must contain all closed deep components and their endpoint owners.
