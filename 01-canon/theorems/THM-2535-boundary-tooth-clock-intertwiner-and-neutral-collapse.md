---
id: THM-2535
title: "Boundary-tooth clock intertwiner, neutral collapse, and clock-holonomy dipole"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The C_13 prime-necklace marker does not by
  itself select one of THM-2517's seven owner epochs: a pure septimal CRT
  translation fixes every root-mask datum and moves the clock freely.  The
  retained THM-2508 affine cut chart supplies exactly the missing C_7
  coordinate.  If the selected wall is s->t=s+tau and the cut is (a,b),
  its tooth-one source is kappa=a^(-1)(1-b); at output t the cut contains
  the distinguished summand d(s,kappa), and kappa transforms by
  kappa->B kappa+C under source-chart transport.  Reindexing the whole cut
  bank by kappa gives an exact C_13-by-C_7 root--clock table whose primitive
  Fourier transform is a nonzero geometric multiplier times the original
  mixed transform.  Thus all 5,184 primitive cut modes survive in a genuine
  clock-covariant local system.  On every target-anchored singleton/adjacent-
  pair deep comb, the marker-to-target chord supplies such a tooth even when
  the adjacent marker wall misses the target; this edge-derived cut slope
  need not be the physical guard slope.  Its neutral clock character vanishes
  identically, sharply: the uniform diagonal sum is zero for every row-zero
  defect.  A correlated selector/scheduler incidence cycle escapes that
  collapse and gives an exact positive marker-to-target root dipole on the
  full clock orbit, but its uncharted signed table cancels identically and it
  is an incidence boundary rather than a temporal transition.  An explicit
  nonnegative anchored-table hostile has all 72 mixed modes while the marker-
  selected tooth is zero.  The construction therefore identifies the correct
  clock cocycle but does not put the old signed cut defect on the Boolean
  selector/owner-path ancestry, make the empty boundary an arrival, emit an
  owner loop, exclude a row, or prove LRC(14).
source: codex-2026-07-27-boundary-tooth-clock-intertwiner
depends_on:
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2517-cubic-anchored-spectrum-flt3-and-three-time-boolean-lift
  - THM-2527-owner-weighted-all-mode-odd-bank-and-boolean-cut-coordinate
  - THM-2531-prime-necklace-guard-boundary-selector
related:
  - THM-2512-lawful-interaction-cut-bundle-transplant-and-replica-dichotomy
  - THM-2528-intrinsic-four-arm-boolean-path-and-joint-autocorrelation-scalarization
  - THM-2532-cyclic-tournament-cayley-algebra-and-chi7-even-quotient
  - THM-2533-owner-weighted-phase-and-mixed-gain-radon-ladders
script: 04-computation/lrc14_boundary_tooth_clock_intertwiner_thm2535.py
output: 05-knowledge/results/lrc14_boundary_tooth_clock_intertwiner_thm2535.out
script_sha256: 5b439d24cff33e57c626f0f8bff5f105b4b1aefdb75896132c905576b37d9073
output_sha256: 51ec45b3ed6c179f5e38b4dbc132bf3821c06815dbc85dd3029af019e41f0232
hash_basis: working-tree bytes (LF)
---

# THM-2535 -- the cut intercept exposes the clock, and correlated holonomy survives

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2517 leaves the phase-zero owner at a free seven-epoch torsor.  THM-2531
now gives every nonconstant thirteen-root mask a canonical guard-directed
occupied-to-empty wall.  Those facts do not compose directly: the marker
lives over `C_13`, while the owner clock lives over `C_7`.

The full affine cut chart supplies the missing bridge, but with a sharp
qualification:

```text
root marker alone                         -> no C_7 section;
root boundary + retained affine cut       -> exact C_7 clock coordinate;
Fourier transform in that clock           -> all primitive mixed modes;
uniform/source-neutral clock contraction  -> identically zero.             (1)
correlated full-clock incidence lift       -> nonzero root dipole only.
```

Thus the algebraic clock-covariant intertwiner exists.  What remains is to
put it on the same positive ancestry object as the selector and late owner.

## 1. The two torsors are independent before a cut is retained

Let

```text
e:F_13->{0,1},                    e nonconstant,
tau in F_13^*.                                                        (2)
```

Use THM-2531's selector notation

```text
alpha=alpha_tau(e),
q=q_tau(e),
s=s_tau(e),                       t=t_tau(e)=s+tau.             (3)
```

Then `e_s=1` and `e_t=0`.  Under a root affine relabelling

```text
h -> U h+H,
```

the oriented slope and selector transform as

```text
tau -> U tau,
(alpha,q,s,t)->(U alpha+H,q,U s+H,U t+H).                      (4)
```

Suppress the large delay in THM-2517's full-support scheduler and write its
seven disjoint cells as

```text
P_kappa=product_(d in F_7) G_(kappa+d)(T^(D_d)x),
                                                   kappa in F_7.       (5)
```

The phase table `gamma=kappa+d` is the cyclic Latin square.  In the fixed
owner gauge, the actual phase-zero owner occurs in `P_kappa` at the unique
slot

```text
d=-kappa.                                                       (6)
```

A source translation by `C` sends the phase table to the row `kappa+C`.
It does not change the thirteen root labels.  Conversely, a pure root
translation changes (3) and not the scheduler row.  These are independent
translations because CRT gives

```text
Z/91Z = F_13 x F_7.                                            (7)
```

This independence gives a formal no-section theorem.

> **Root-only clock no-go.**  No function of the full root mask, its
> lexicographic marker, run length, selected wall, word tournament, root
> autocorrelation, or four-arm root addresses can be a translation-covariant
> choice of `kappa in F_7`.

Indeed, use the CRT translation `K=13`.  It has

```text
K=0 mod 13,                       K=6 mod 7.                    (8)
```

It fixes every root-only input to the proposed function but would have to
add six to its clock output, a contradiction.  Reducing the binary marker
score, the maximum run length, or a tournament rank modulo seven cannot
repair (8).  The problem is not insufficient ingenuity in choosing a root
statistic; it is a free action on a coordinate absent from the input type.

For completeness, an affine source relabelling `gamma -> B gamma+C` takes
the Latin rule to

```text
B gamma+C=(B kappa+C)+B d.                                    (9)
```

Thus `kappa->B kappa+C` only after the slot label is transported by
`d->B d`.  THM-2517's seven physical epochs `D_d` are not themselves fixed
by that permutation.  The existing scheduler has an exact translation law;
the full affine statement in (9) is a transported schedule, not an equality
of the original time-labelled cells.  This distinction will be retained.

## 2. Every affine cut is an ordered seven-tooth source chart

Let

```text
d:F_13 x F_7 -> Q,
sum_(r in F_7)d(h,r)=0                         for every h,      (10)
```

and recall THM-2508's cut operator

```text
R_(tau,a,b)(v)
 =sum_(r in F_7)d(v-tau rep(ar+b),r),

a in F_7^*,                    b in F_7.                       (11)
```

For a tooth label `j in F_7`, define

```text
kappa_j(a,b)=a^(-1)(j-b).                                     (12)
```

Changing variables from `r` to `j=ar+b` rewrites the whole cut as

```text
R_(tau,a,b)(v)
 =sum_(j=0)^6 d(v-tau j,kappa_j(a,b)).                         (13)
```

Thus `(a,b)` is not merely one coefficient pair.  It is an ordered matching
between all seven source rows and the seven-root ray

```text
v, v-tau, ..., v-6tau.                                        (14)
```

At the marker-selected output `v=t`, tooth one is exact:

```text
kappa:=kappa_1(a,b)=a^(-1)(1-b),

t-tau=s,

the j=1 summand in R_(tau,a,b)(t) is d(s,kappa).               (15)
```

This is the first exact incidence between the selected root boundary and a
source coordinate.  Each fixed `a` makes `b->kappa` a bijection, and the
full `42`-cut atlas contains exactly six charts over every `kappa`.

The covariance is equally exact.  Transport coordinates by

```text
h'=U h+H,                         r'=B r+C.                    (16)
```

The same geometric cut is represented in the new chart by

```text
tau'=U tau,                 v'=U v+H,
a'=a B^(-1),                b'=b-a B^(-1)C.                   (17)
```

Then `a'r'+b'=ar+b`, so every tooth label is fixed and

```text
kappa'_j=B kappa_j+C,
v'-tau'j=U(v-tau j)+H.                                      (18)
```

In particular, the distinguished coordinate in (15) obeys

```text
kappa' = B kappa+C.                                          (19)
```

Equation (19) is precisely the affine law of a source point.  It is not a
root-to-source homomorphism: the retained cut chart carries the missing
source origin and scale.

### Every deep-comb mask has a marker-to-target clock chord

Nothing in (12)--(19) requires the ordered edge to be adjacent in the guard
cycle.  It requires only two distinct marked roots `s,t` and the derived
nonzero slope

```text
tau_edge=t-s.                                                  (19a)
```

This removes one apparent endpoint restriction on THM-2531's actual deep
comb.  In its target-anchored chart, root `0` is excluded and every mask is

```text
{j},                  1<=j<=12,
{j,j+1},              1<=j<=11.                               (19b)
```

For any retained marker direction, let `alpha` be the occupied marker root.
Then

```text
s=alpha,                 t=0,                 tau_edge=-alpha (19c)
```

is an occupied-to-empty **selector-star chord** on every one of the `23`
masks.  At output zero, tooth one of every chart is

```text
d(alpha,kappa),             kappa=a^(-1)(1-b).                 (19d)
```

Under a root affine chart change both endpoints transport, so
`tau_edge->U tau_edge`; the chord is covariant.  This bypasses THM-2531's
adjacent-wall endpoint-diagonal restriction at the level of the full cut
atlas: one no longer needs the adjacent selected wall itself to end at root
zero.

The qualification is load-bearing.  The slope in (19c) is selected from the
marker-to-target chord.  It need not equal THM-2526/2527's independently
retained physical guard slope.  The full THM-2508 bank contains every
nonzero slope, so the signed algebra below applies; the guard-oriented
positive odd path does not automatically move to this slope.  Nor does the
empty target acquire arrival semantics.  The corollary supplies a clock
chart on every deep fibre, not a positive target current.

There is a faithful tournament interpretation, but it adds no shortcut.
THM-2531's thirteen binary words give a transitive tournament with source
`alpha`.  A cut `(a,b)` gives a separate transitive order of the seven
source rows through `j=0,...,6`.  Equation (15) joins one root wall to the
rank-one tooth of that order.  The object is a typed bipartite incidence,
not a tournament on twenty vertices.  Forgetting `(a,b)` forgets both the
source order and `kappa`; reversing `a` reverses its oriented ray.

## 3. Boundary-relative cut coordinates match the owner scheduler

For fixed `tau,a`, reindex the cut translation by its tooth-one source:

```text
Q_(tau,a)(kappa,t)
 :=R_(tau,a,1-a kappa)(t),

(kappa,t) in F_7 x F_13.                                    (20)
```

This is a bijective change of coordinates in the full cut bank, not a new
observer.  Its semantic advantage is that the `j=1` summand is always

```text
d(t-tau,kappa).                                              (21)
```

Pair the chart in (20) with THM-2517's cell `P_kappa`.  Under a pure source
translation `C`,

```text
kappa -> kappa+C,
b=1-a kappa -> b-aC=1-a(kappa+C),                             (22)
```

so the diagonal family of pairs is permuted.  In the fixed owner gauge its
actual phase-zero block lies at slot `-kappa`.  With the transported-slot
qualification in (9), the same statement obeys the full affine law (19).

Consequently (20)--(22) give a **clock-covariant section of the product
label bundle**.  They do not give one fixed physical epoch.  More
importantly, `R` is still a signed contraction of a response-table defect,
whereas `P_kappa` is a delayed Boolean scheduler cell.  Canon has not yet
constructed their product on one ancestry fibre.

THM-2527/2528 supply the complementary positive half.  Replacing their late
owner factor by a sufficiently delayed positive `P_kappa` can retain the
positive odd root path in a clock-labelled cell, and THM-2531 can retain its
selected wall.  That path has no source-row defect `d`.  Equation (20) gives
the source-row local system but has no positive selector/path semantics.
The missing operation is exactly a same-ancestry pairing of these two halves.

## 4. Fourier transform in the clock is lossless

Let `zeta=zeta_13`, `xi=zeta_7`, and define

```text
dtilde(alpha,beta)
 =sum_(h,r)d(h,r)zeta^(-alpha h)xi^(-beta r),                  (23)

Qhat_(tau,a)(alpha,beta)
 =sum_(t,kappa)Q_(tau,a)(kappa,t)
                   zeta^(-alpha t)xi^(-beta kappa).            (24)
```

In (13), put `b=1-a kappa`.  For fixed `r,j`,

```text
kappa=r+a^(-1)(1-j).                                         (25)
```

Also put `h=t-tau j`.  Direct reindexing gives

```text
Qhat_(tau,a)(alpha,beta)
 =xi^(-beta a^(-1)) K^partial_(tau,a)(alpha,beta)
    dtilde(alpha,beta),                                      (26)

K^partial_(tau,a)(alpha,beta)
 =sum_(j=0)^6
    (zeta^(-alpha tau)xi^(beta a^(-1)))^j.                    (27)
```

For `alpha,tau!=0`, let

```text
lambda=zeta^(-alpha tau)xi^(beta a^(-1)).                     (28)
```

Then

```text
lambda^7=zeta^(-7 alpha tau)!=1.                              (29)
```

Hence `lambda!=1` and the geometric sum in (27) is nonzero.  Therefore,
for every `alpha,beta!=0`,

```text
Qhat_(tau,a)(alpha,beta)!=0
 iff dtilde(alpha,beta)!=0.                                  (30)
```

This is exactly THM-2508's cut-character intertwiner in boundary-relative
coordinates.  If `Psi` denotes its cut-translation transform, then

```text
Qhat_(tau,a)(alpha,beta)
 =xi^(-beta a^(-1))
   Psi_(tau,a)(alpha,-beta a^(-1)).                           (31)
```

Thus no new nonvanishing premise has been smuggled in.  Rather, (31)
identifies the old nontrivial cut character with the nontrivial **owner-clock
character**.  Whenever the lawful defect has all `72` primitive mixed modes,
the reindexed bank has all

```text
12*6*12*6=5,184                                             (32)
```

primitive root--clock coefficients.

This is the exact algebraic clock-covariant intertwiner requested by the
full-support branch of THM-2517.  Its carrier is a signed local system, not
yet an emitted Boolean current.

## 5. The neutral clock character vanishes sharply

At `beta=0`, (23) vanishes by the row-zero law (10).  There is also a direct
pointwise identity.  For every fixed `tau,a,t`,

```text
sum_(kappa in F_7)Q_(tau,a)(kappa,t)

 =sum_(kappa,r)
   d(t-tau rep(1+a(r-kappa)),r)

 =sum_(j=0)^6 sum_r d(t-tau j,r)

 =0.                                                         (33)
```

The second equality uses that `kappa->1+a(r-kappa)` is a bijection for each
`r`.  Thus:

```text
source-neutral uniform clock sum = zero,
nontrivial clock character       = lossless mixed source phase.            (34)
```

This is the exact cocycle obstruction behind the scheduler tradeoff.  Fixing
one `kappa` retains a clock but breaks source neutrality.  Summing all clocks
with equal scalar weight restores neutrality and kills the row-zero signal.
If a delayed pairing asymptotically factorizes so that all seven scheduler
cells contribute one common scalar mean, its leading term is that common
mean times (33), hence zero.  A surviving physical pairing must retain a
nontrivial clock character or a finite-delay correlation with the response;
mere independent mixing cannot supply it.

Equation (33) does not kill THM-2527's positive odd path, because that path
has no row-zero `F_7` defect.  It says that attaching only a formal chart
label to that positive event does not transfer the mixed source charge.

## 6. A correlated selector/scheduler incidence cycle escapes the collapse

The uniform collapse is sharp, but it is not the end of the clock thread.
The selector and scheduler can be correlated before the chart is forgotten.
This produces a genuine clock-labelled signed incidence current.

Work in the full-support branch and in the target-anchored deep chart.  Let
`E_a` be any positive Boolean selector class whose occupied marker is
`a!=0`, and put

```text
m_a=measure(E_a)>0,
X_kappa^L=E_a intersection P_kappa^L,
x_kappa^L=measure(X_kappa^L).                                 (34a)
```

The seven `X_kappa^L` are disjoint because the scheduler cells are disjoint.
The higher-order mixing argument of THM-2517 remains valid with the fixed
finite-horizon factor `1_(E_a)` inserted.  If

```text
P=product_(gamma in F_7)q_gamma>0,
```

then, uniformly over the seven clocks,

```text
x_kappa^L -> m_a P,
M_a^L:=sum_kappa x_kappa^L -> 7m_aP>0.                        (34b)
```

Hence `M_a^L>0` for every sufficiently large admissible delay.
The classes `E_a`, `a in F_13^*`, partition the target-anchored mixed-mask
locus.  Retain the full `a`-bank: root affine changes permute it, every
positive class supplies the construction below, and zero classes contribute
zero.  Thus existence of a positive class never requires choosing a root
origin or breaking the marker gauge.

For each `kappa`, form the two-clock signed scalar table

```text
D_kappa^L(h,r)
 =1_(h=a)[x_kappa^L 1_(r=kappa)
          -x_(kappa+1)^L 1_(r=kappa+1)].                     (34c)
```

Its positive and negative coefficients are masses of the two disjoint
Boolean events `X_kappa^L` and `X_(kappa+1)^L`.  Use the marker-to-target
slope and the adjacent clock chart

```text
tau_edge=-a,                    cut scale=1,
cut intercept=-kappa.                                            (34d)
```

Source `kappa` is tooth zero and source `kappa+1` is tooth one.  Directly
from (11),

```text
R_kappa^L(v)
 =x_kappa^L 1_(v=a)-x_(kappa+1)^L 1_(v=0).                   (34e)
```

Summing the entire covariant clock orbit gives the exact root dipole

```text
J_a^L(v):=sum_kappa R_kappa^L(v)
          =M_a^L[1_(v=a)-1_(v=0)].                           (34f)
```

It is nonzero for all sufficiently large `L`.  Every one of its twelve
nontrivial root Fourier coefficients is nonzero:

```text
Jhat_a^L(alpha)=M_a^L(zeta^(-alpha a)-1)!=0,
                                             alpha in F_13^*. (34g)
```

This construction is exactly source-translation covariant.  Translation by
`C` sends `kappa->kappa+C`, the chart intercept to `-kappa-C`, and permutes
the seven summands in (34f).  Root affine transport moves the chord and its
slope as in (19c).  Thus (34f) is source-neutral without choosing one owner
epoch, while every labelled summand retains which epoch contains the actual
phase-zero owner.

The escape mechanism is a nontrivial chart cocycle.  Each local table has

```text
sum_r D_kappa^L(a,r)=x_kappa^L-x_(kappa+1)^L ->0,             (34h)
```

but need not be row-zero at finite delay.  More strikingly,

```text
sum_kappa D_kappa^L=0                                      (34i)
```

exactly: every event mass appears once positively and once negatively in
the same uncharted `(h,r)` cell.  The nonzero dipole (34f) survives only
because each local table is evaluated in its transported chart before the
charts are summed.  Forget the chart first and the whole construction is
zero.  This is why it does not contradict (33).

Equations (34a)--(34i) are a rigorous improvement over a purely formal pair
of source columns.  The two coefficients come from actual disjoint Boolean
selector/scheduler intersections, and positivity follows at a finite large
delay.  They still stop short of a semantic LRC current for three reasons.

1. The `F_7` coordinate is the scheduler clock, not the inherited
   THM-2449 source/deep residue.  Thus (34f) does not transplant the old
   mixed source charge.
2. Every `X_kappa^L` is used in two incidence components.  This is a lawful
   boundary in a labelled local system, not one observed temporal
   transition of that event.
3. Root zero is empty for the predecessor mask.  The negative endpoint in
   (34f) is an incidence boundary coefficient, not a positive arrival event.

The construction closes the **abstract clock-cocycle and positive-mass
incidence** problem.  It uses the whole torsor cycle and is explicitly not an
equivariant fixed-epoch section.  It leaves the empty-target-to-arrival and
inherited source/deep-semantic transfers open.

## 7. Full mixed spectrum does not force the selected tooth

The distinguished summand in (15) can vanish even when every primitive mode
of the cut bundle survives.  Define the integral defect

```text
d(h,r)
 =(1_(h=0)-1_(h=2))(1_(r=0)-1_(r=1)).                        (35)
```

Both margins vanish.  Moreover

```text
dtilde(alpha,beta)
 =(1-zeta^(-2alpha))(1-xi^(-beta))!=0                         (36)
```

for every `alpha,beta!=0`.  Hence (30) gives all `5,184` primitive
root--clock coefficients.

This defect occurs inside the nonnegative rational anchored-table axioms.
Let `I_(r,h)=d(h,r)` and define

```text
A_(r,0)=1_(r=0),

A_(r,h)=d(h,r)-d(0,r)+1_(r=0)+1,              h!=0.           (37)
```

Every entry is nonnegative, the column `h=0` is the positive delta anchor,
and the owner row is nonconstant.  The last three terms after `d` in (37)
are row and column main effects, so the ANOVA interaction of `A` is exactly
`I`.  This is an algebraic anchored response table; it is not asserted to
be a physically realized THM-2449 row.

Now take the singleton root mask

```text
e=1_{1},                         tau=1.                        (38)
```

Its THM-2531 selector is

```text
alpha=1,             q=1,             (s,t)=(1,2).            (39)
```

But (35) has

```text
d(s,kappa)=d(1,kappa)=0                 for every kappa.       (40)
```

In the target-anchored deep chart the same singleton has the universal
selector-star chord

```text
1 -> 0,                         tau_edge=-1.                   (40a)
```

Its distinguished tooth is still `d(1,kappa)=0` for every `kappa`, while
the full root--clock bank at `tau_edge=-1` retains every primitive mode.
Thus passing from the adjacent wall to the universal target chord removes
the incidence-domain obstruction but not the charge obstruction.

Thus all `42` boundary charts have zero distinguished tooth while their
complete primitive Fourier bank is nonzero.  The other six teeth carry the
spectrum.  Equations (38)--(40) are a finite hostile to the inference

```text
selected actual boundary + complete mixed spectrum
  => selected source tooth is charged.                                 (41)
```

The hostile is possible precisely because current canon has no joint law
coupling the full root mask to the signed ANOVA defect.  It does not refute
such a law once a common physical ancestry construction is supplied.

## 8. Exact gain and remaining proof obligation

The C_7 clock ledger is now sharper.

1. THM-2531's root marker is not and cannot be a clock section by itself.
2. The full affine cut chart carries an intrinsic boundary-relative clock
   `kappa=a^(-1)(1-b)`.
3. Every target-anchored deep-comb mask has a marker-to-target version of
   that chart, even when its adjacent marker wall misses the target.
4. In that coordinate, the cut character is exactly the owner-clock
   character, and all primitive modes survive losslessly.
5. The invariant neutral character is forced to zero.
6. Correlating the selector with the whole scheduler orbit before applying
   the transported charts gives an exact positive root dipole; forgetting
   those charts makes its signed incidence table zero.
7. A selected root wall or target chord need not meet the old charged tooth
   without an additional joint-ancestry theorem.

THM-2532/2533 sharpen, but do not erase, the seventh item.  They now supply
an actual late-owner-attached centred predecessor profile and `216`
slope-labelled nonzero Cayley/Radon currents, and one such oriented current
recovers that profile by a signed sawtooth inverse.  What they explicitly do
not supply is an equality between those currents and the ordered THM-2508
cut/carry defect, or a coupling to the positive selector incidence used
here.  Thus the common-profile problem is closed while the common-ancestry
and target/terminal intertwiner remain open.

The next lawful object is therefore not another marker, tournament, or
colour census.  It must be one of:

```text
a same-ancestry identification of the incidence cycle (34c) with the old
Q(kappa,t), THM-2528 path, or THM-2449 source/deep coordinate;

a proof that one marker-selected tooth d(s,kappa) is a positive semantic
source--arrival current on the live family; or

a different clock-covariant Boolean local system carrying both the
THM-2528 path imbalance and the THM-2508 source charge.                   (42)
```

Even then, THM-2531's empty endpoint `t` remains absence of the predecessor
event, not a semantic arrival.  No owner loop, terminal-current emission,
row exclusion, or proof of LRC(14) follows here.

## 9. Exact companion

Run

```bash
python3 04-computation/lrc14_boundary_tooth_clock_intertwiner_thm2535.py
python3 -O 04-computation/lrc14_boundary_tooth_clock_intertwiner_thm2535.py
```

The dependency-free referee works in `F_547`, which contains primitive
seventh and thirteenth roots.  It verifies:

- the `7 x 7` Latin assignment law, all owner slots, source translations,
  and the transported affine slot law;
- all `57,330` nonconstant-mask/owner-epoch states and the pure septimal
  CRT hostile `K=13`;
- all `6,552` boundary-tooth incidences and six-to-one chart fibres;
- all `966` marker-to-target chord incidences on the `23` target-anchored
  deep-comb masks;
- `12,348` source-chart and `170,352` root-chart covariance identities;
- `73,008` uniform-diagonal identities on a basis of every row-zero defect;
- all `6,048` primitive geometric kernels and `48,384` direct DFT
  factorizations on eight independent exact controls; and
- the `72`-mode hostile (35), its `42` zero selected teeth, and all `5,184`
  surviving root--clock modes; and
- `7,644` correlated selector/scheduler cut identities, `7,644` entries of
  its identically zero forgotten table, and `1,092` orbit-dipole identities.

Normal and optimized executions reproduce the stored transcript
byte-for-byte. **QED.**

The independent audit rederived the CRT no-section argument, source/root
chart covariance, the boundary-relative Fourier factorization and its
neutral-character collapse.  It separately checked the marker-to-target
edge slope, the nonnegative anchored-table hostile, and the full-clock
incidence cycle: its uncharted table telescopes to zero while its transported
cut sum is exactly `M_a^L(delta_a-delta_0)`.  Normal and optimized executions
reproduced the stored transcript byte-for-byte.
