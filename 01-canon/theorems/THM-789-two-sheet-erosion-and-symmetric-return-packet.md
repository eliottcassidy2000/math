---
id: THM-789
title: Two-sheet erosion, pointwise thickness tax, and the symmetric return-packet upgrade; exact local-trapping obstruction
status: PROVED (elementary metric erosion + Kneser/Macbeath; exact obstruction independently refereed and replayed)
source: codex-2026-07-14-S10 phase-noncontainment analysis
depends_on:
  - THM-774   # folded-diamond equivalence
  - THM-780   # quantitative phase pigeonhole
  - THM-782   # one-sided return packet, strengthened here
  - LRC(<=13)
related: [THM-772, THM-775, THM-776, HYP-6820]
verification: 04-computation/lrc14_two_sheet_erosion_trap_codex_S10.py
  (+ 05-knowledge/results/lrc14_two_sheet_erosion_trap_codex_S10.out)
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

## 4. Tournament and assumption challenge

Runner vertices, gap vertices, and fixed-time phase-cell vertices all lose the
global choice exposed by (11). A more faithful candidate vertex set is the set
of connected components of `E_U`. Give component `C` the escape observable

`m(C)=sup_(t in C)(11/13-Q_(a,b)(t))`,

orient pairs toward larger `m`, and use left endpoint as a tie gauge. The tie
Hamiltonian path is the sorted component list. This tournament is transitive
(score histogram `{0,...,r-1}`, no directed cycles, singleton SCCs, one
Hamiltonian path); its value is not nontransitive structure but preservation of
the **global alternatives**. In the exact example the escaping `14/19` choice
beats the locally trapped `4/17` choice.

Even this quotient must retain component intervals, owner labels, and the
pointwise margin; its bare isomorphism class destroys the LRC predicate. The
challenged assumption is that refining one successful local phase anchor adds
information. It is monotone trapping here. Recursion must branch across deep
components, not only downward inside one phase cell.
