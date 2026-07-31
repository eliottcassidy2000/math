---
id: THM-2670
title: "Sharp-graph clock-incidence atlas and physical-gluing boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On the canonical THM-2616 raw nonnegative carrier, impose the THM-2629
  graph q=h, r=-h-1 and retain the THM-2630 predecessor digit, half-tooth,
  and binary carry before integration.  The full delayed guard cospan gives
  an exact family of 13-by-13 common-x incidence matrices.  Every one of the
  72 fixed-source constant-clock-step sevenfold Boolean products is zero,
  but only because its clock cycle contains a universally zero matrix.
  This does not survive varying clock steps or clock marginalization.  The
  honest seven-clock graph has a rich 32-edge Hamiltonian generic class and
  two reflected 19-edge exceptional classes; its matrix labels are nowhere
  nontrivially translation-equivariant.  Formal generic safe and guard-free
  Hamiltonian products service all 13 state displacements but do not service
  all 169 endpoint pairs, while the danger labels obey a separate U^3=0
  state law.  These are exact formal label products, not a typed handoff
  between the disjoint owner-clock strata, chronology, a transition, or an
  LRC(14) conclusion.
source: root-2026-07-28-sharp-clock-incidence
depends_on:
  - THM-2616-cross-time-target-future-diagonal-and-principal-action-no-go
  - THM-2624-two-clock-root-tomography-and-disjoint-carrier-holotopy-boundary
  - THM-2629-fixed-deep-affine-graph-spectrum-and-puncture-cancellation-boundary
  - THM-2630-old-wall-affine-clutching-and-successor-sector-no-go
  - THM-2642-cyclic-difference-relation-saturation-and-thick-holotopy-no-go
related:
  - THM-2623-guard-safe-danger-cospan-and-residual-unit-wall
  - THM-2637-derangement-character-fixed-branch-holotopy-principle
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2680-dilation-reversed-two-edge-clock-fibre-products-and-source-drift-boundary
  - THM-2682-central-arrival-clock-trap-and-three-event-dilation-nilpotence
  - THM-2684-three-tooth-rail-envelope-diagonal-arrival-law-and-full-dilation-nilpotence
scripts:
  - 04-computation/lrc14_successor_private_sharp_graph_clock_collapse.py
  - 04-computation/lrc14_guard_cospan_successor_private_clock_collapse.py
  - 04-computation/lrc14_clock_graph_hamiltonian_audit.py
outputs:
  - 05-knowledge/results/lrc14_successor_private_sharp_graph_clock_collapse.out
  - 05-knowledge/results/lrc14_guard_cospan_successor_private_clock_collapse.out
  - 05-knowledge/results/lrc14_clock_graph_hamiltonian_audit.out
---

# THM-2670 -- the sharp graph has a formal clock atlas, not yet a physical one

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The common-`x` carrier has repeatedly produced dense local root support but
no lawful transition.  This theorem retains one more coordinate: the
numerical predecessor digit before integration.  The resulting object is
best viewed neither as a scalar support table nor as a tournament.  It is a
seven-vertex directed graph whose arrows carry Boolean relations on thirteen
digit states.  That reframing explains both the striking constant-step zeros
and why they do not close LRC(14).

## 1. The exact common-x incidence bank

Put

```text
p=13,                    clocks=F_7,
d in F_7^*,              s in F_13^*,
ell_5=ell_4+d.                                             (1)
```

Start with the `162` positive middle rails of THM-2616.  Impose the sharp
coefficient graph of THM-2629,

```text
q=h,                    r=-h-1 mod 13,                    (2)
```

and refine before integration by the predecessor digit `j`, the later-tooth
half `epsilon in {0,1}`, and binary carry `kappa in {0,1}`.  THM-2630's
pointwise identity becomes

```text
r=2j+kappa+epsilon,
h=-2j-kappa-epsilon-1                 mod 13.             (3)
```

For a delayed sector

```text
sigma in {safe,danger,free},
free = safe disjoint-union danger,                         (4)
```

For an inherited middle-rail occurrence `alpha`, let

```text
A^sigma_(alpha,d,s,ell;j,h,epsilon,kappa) >= 0            (5)
```

be the exact raw integer numerator of this common-`x` atom.  Define the
Boolean incidence relation by forgetting the rail and binary refinements:

```text
K^sigma_(d,s,ell)(j,h)=1
  iff some A^sigma_(alpha,d,s,ell;j,h,epsilon,kappa)>0.    (6)
```

The word `predecessor` here is numerical: `(6)` does **not** assert that `j`
is a semantic owner or that `h` is an adjacent-clock input state.  Counts
below retain the rail label `alpha`; they are labelled occurrences, not
counts of distinct Boolean matrix entries.

The exact universe and partition checks are:

```text
guard-safe initial atom candidates:              42,768
guard-cospan atom candidates:                     46,656
guard-cospan atom partitions:                     46,656
sharp-graph partitions:                           11,664
h=12 -> r=0 structural zero checks:                   972. (7)
```

The positive labelled atom-occurrence counts are

```text
safe=7,436,             danger=1,636,             free=8,360. (8)
```

These counts need not add: a single displayed atom label can have positive
mass in both disjoint physical sectors.  Exact nonnegative additivity gives

```text
A^free_(alpha,...)=A^safe_(alpha,...)+A^danger_(alpha,...),
K^free=K^safe OR K^danger                              (9)
```

entrywise.  The danger leg restores `h=0`; `h=12` would force the absent
sheet `r=0` and remains structurally zero.

## 2. The complete constant-step zero atlas

For each sector, the zero matrices among the `6*12*7=504` triples
`(d,s,ell)` are the same.  Their number at each step is

```text
d:          1   2   3   4   5   6
zero cells: 18  39  16  16  39  18.                     (10)
```

More decisively, the clock labels at which the matrix vanishes for **every**
source `s` are

```text
d=1: {5},          d=2: {0,4,6},      d=3: {2},
d=4: {5},          d=5: {0,1,3},      d=6: {2}.          (11)
```

Since every nonzero `d` traverses all of `F_7`, `(11)` proves

```text
K^sigma_(d,s,0) K^sigma_(d,s,d) ... K^sigma_(d,s,6d)=0  (12)
```

for every sector `sigma`, all six steps, and all twelve sources: `72/72`
fixed-source, constant-step candidate products vanish in each sector.

This zero is completely explained by an individually zero factor.  It is
not a holonomy cancellation.  In particular the danger leg fills none of
the eighteen `d=1` holes, so `(10)`--`(12)` are invariant under the delayed
guard split rather than caused by the guard-safe truncation.

## 3. Sharp hostiles to a broader clock no-go

The conclusion of `(12)` fails as soon as the clock itinerary is allowed to
vary.  For every generic source

```text
s in {1,2,3,4,5,8,9,10,11,12},                           (13)
```

the itinerary

```text
(0,1,2,3,6,4,5)                                          (14)
```

has safe/free formal support `70/104`, and

```text
(0,6,2,4,3,5,1)                                          (15)
```

has safe/free support `88/143`.  Both supports are zero at `s=6,7`.

Forgetting the initial clock label is even more destructive.  Put

```text
M^sigma_(d,s)=OR_(ell in F_7) K^sigma_(d,s,ell).          (16)
```

Then the guard-free seventh power has the exact shape

```text
(M^free_(d,s))^7 = F_13 x (F_13\{12})                    (17)
```

for every one of the `72` pairs `(d,s)`, hence support `156`.  The safe
seventh powers have support `110` in `71` cases and `75` at `(d,s)=(5,6)`;
the danger powers are zero.  Thus neither varying steps nor clock
marginalization inherits the constant-step zero.

Even the carrier-free affine grammar in `(3)`, together with the elementary
half compatibility

```text
kappa=0 => h<=6,             kappa=1 => h>=6,             (17a)
```

is unobstructed: its Boolean powers have supports

```text
27, 56, 113, 169.                                        (18)
```

It becomes the full relation after four steps.  The holes in `(10)` are
therefore genuine intersections with the canonical carrier, not a modular
artifact of the predecessor formula.

## 4. The honest seven-clock graph

For fixed `s`, use the seven clocks as vertices and draw

```text
ell -> ell'  iff  K^sigma_((ell'-ell),s,ell) is nonzero.  (19)
```

The scalar graph is the same in all three guard sectors.  The full labelled
matrix bank has exactly three source classes:

| source class | arrows | unordered-pair signature `(both,one,none)` | normalized Hamilton cycles | Hamilton paths |
|---|---:|---:|---:|---:|
| the ten sources in `(13)` | 32 | `(13,6,2)` | 92 | 876 |
| `{6}` | 19 | `(3,13,5)` | 0 | 0 |
| `{7}` | 19 | `(3,13,5)` | 0 | 0 |

The exceptional graphs are related by `ell -> -ell`.  At source `6`, clocks
`1,2` have indegree zero; at source `7`, clocks `5,6` do.  Since a directed
spanning path can start at only one vertex, this proves the absence of a
Hamilton path without enumeration.  Each exceptional graph nevertheless
has exactly fourteen directed paths through six clocks.  In the generic
class only `16/92` normalized cycles survive reversal, and none of the six
constant-step cycles occurs.

This is not a tournament: bidirected and absent pairs are intrinsic data.
Forcing an orientation would erase precisely the two phenomena that explain
the generic/exceptional split.

## 5. Matrix labels, displacements, and the danger nilpotent

The scalar graph alone is insufficient.  Compose its `13`-state edge labels
formally by Boolean matrix multiplication, which identifies a printed output
digit with the next printed input digit.

For a generic source, every one of the `92` scalar Hamilton cycles has a
positive safe product and a positive free product.  Their pair supports
satisfy

```text
safe: 49 <= |support| <= 88,    free: 87 <= |support| <= 143. (20)
```

and every product meets all thirteen cyclic diagonals

```text
h-j=c,                    c in F_13.                      (21)
```

None is the full `169`-pair relation.  Ordinary, rather than Boolean,
products give nonuniform state-path counts; no cycle has all thirteen formal
cyclic-displacement totals divisible by `13`.  Thus `(21)` is not
THM-2642's equivariant
convolution law in disguise.

The danger sector has the same scalar clock arrows, but every Hamiltonian
matrix product is zero.  This has a clock-independent explanation.  After
forgetting all clock and step labels its union relation is exactly

```text
U={(0,11),(5,1),(6,0),(12,0),(12,1)},
U^2={(6,11),(12,11)},              U^3=0.                (22)
```

Every danger edge is a subrelation of `U`, so **every formal product of any
three danger matrices is zero**, regardless of its clock labels.  The danger
obstruction is state nilpotence, not failure of the clock skeleton.

## 6. Why equivariant closure loses the information

Among the `504` edge matrices in each sector, `146` are zero and `358` are
positive.  No positive matrix is invariant under simultaneous translation

```text
(j,h) -> (j+a,h+a).                                      (23)
```

The numbers of occurring increments `h-j` among the positive edges have
histograms

```text
safe:   9^68 10^101 11^189,
danger: 1^24  3^57   4^44  5^233,
free:  10^24 11^101 12^233.                              (24)
```

The least translation-equivariant relation containing an edge with
increment set `D` is `F_13 x D` in difference coordinates.  Consequently,
THM-2642 implies that the composition of **any two** positive safe hulls
services every clutch class with at least

```text
13(9+9-13)=65                                            (25)
```

formal sections; for free hulls the floor is

```text
13(10+10-13)=91.                                         (26)
```

Thus equivariant closure immediately saturates the desired obstruction.  It
is a hostile quotient, not a repair of the missing transition.

## 7. The missing object is a typed chart handoff

> **Subsequent resolution (THM-2680/2682/2684).**  The candidate `D` below does
> produce positive physical two-edge fibre products, but the arrival-six
> return forces a repeated shallow clock and makes every three-event product
> empty.  The full inherited rail bank is a three-tooth envelope on which `D`
> has only diagonal arrival transitions; its endpoint three-returns also force
> a repeated clock.  Thus the formal/physical distinction established here
> was load-bearing; this handoff is closed on the entire parent rail bank.

An entry of `(6)` means

```text
there exists x in E_(ell,ell';j,h) subset O_ell,          (27)
```

where `O_ell` is the corresponding owner-clock stratum.  Boolean
multiplication replaces two instances of `(27)` by

```text
E_(ell0,ell1;j,h) is nonempty
and E_(ell1,ell2;h,k) is nonempty,                        (28)
```

possibly with two different witnesses in two different clock charts.  A
literal intersection is not the repair:

```text
E_(ell0,ell1;j,h) intersect E_(ell1,ell2;h,k) = empty      (29)
```

for distinct consecutive clocks, because THM-2624 proves the owner-clock
strata `O_ell` are physically disjoint.  A physical composition instead
requires additional typed data, for example either

```text
T_(ell0,ell1): O_ell0 -> O_ell1
with T_(ell0,ell1)(E_01) intersect E_12 nonempty,          (30)
```

or two maps from successive edge events to one shared boundary object and a
nonempty fibre product there.  Such a map must preserve the state label and
the target predicate while recording which component information it
forgets; the Boolean matrices supply none of this.  For a cycle the handoffs
must also compose coherently, not merely exist edge by edge.

Thus the forgotten datum is the witness/component **together with its typed
chart-handoff or descent map**.  The next exact experiment is to enumerate
candidate handoffs already supplied by the natural-extension/endpoint
sidecars, test their two-edge boundary fibre products, and then audit the
higher Čech nerve along the `92` generic Hamilton cycles.  Direct same-`x`
intersection and pairwise nonemptiness are both insufficient.

There is one canonical candidate that fixes the orientation of this test.
Let

```text
D(x)={13x}.                                                (31)
```

Away from half-open boundaries, the numerical digits satisfy

```text
j(Dx)=h(x),                                                (32)
```

and the shallow clock phase at `Dx` is the rail-owner phase at `x`:

```text
{13Dx}={13^2x}.                                           (33)
```

Thus the `D`-compatible candidate chronology traverses the clock arrow
**opposite** to the bookkeeping convention `(19)`.  If `a` is the shallow
clock and `b` the rail clock, its oriented label is

```text
B^sigma_(a,b)=K^sigma_(a-b,s,b).                          (34)
```

The first genuine two-edge test is therefore not `(29)`, but

```text
E^(s0)_(b,a;j,h) intersect
  D^(-1) E^(s1)_(c,b;h,k),                                (35)
```

with all source pairs `(s0,s1)` retained.  Equation `(32)` supplies exactly
the formal intermediate equality `h_current=j_next`; it does not prove that
the other carrier factors, source shift, or component labels survive `D`.
Reversing the scalar graph preserves its cycle counts, but the matrix labels
are noncommutative.  The exact fixed-source reversed audit gives

```text
                         positive cycles  support range  all differences
safe                           92             45..100          92
danger                          0               --              0
free                            92             82..144          92,

positive Hamilton paths:       safe 876 (30..106),
                                free 876 (65..154),
                                danger 0.                     (36)
```

None of the `184` positive reversed products is the full `169`-pair relation,
and all `72` reversed constant-step candidates still vanish in every sector.
The changed ranges relative to `(20)` prove that orientation is not cosmetic,
even though the qualitative hostile survives.  This exhausts the formal
fixed-source reversal only; the physical `D`-pullback sets and source
transport in `(35)` remain the next exact audit.

No common-witness gluing, adjacent-clock chronology, unit row, endpoint
owner, positive transition, holonomy trivialization, scalar row exclusion,
or LRC(14) conclusion follows from this theorem.

## 8. Reproduction and audit boundary

Run each companion normally and in optimized mode:

```bash
python3 04-computation/lrc14_successor_private_sharp_graph_clock_collapse.py
python3 -O 04-computation/lrc14_successor_private_sharp_graph_clock_collapse.py
python3 04-computation/lrc14_guard_cospan_successor_private_clock_collapse.py
python3 -O 04-computation/lrc14_guard_cospan_successor_private_clock_collapse.py
python3 04-computation/lrc14_clock_graph_hamiltonian_audit.py
python3 -O 04-computation/lrc14_clock_graph_hamiltonian_audit.py
```

Each execution must byte-match its stored output.  The LF-normalized SHA-256
hashes `(script,output)` are

```text
successor/private incidence:
  f49eb399af964dd0a6efab5ea5b48c485d9e3368d9a9fb24cbd0a1ddd586fdff
  b43715e745f01089d3c8e236af64121dd2e0b9ce51a89faf662a7d2d449453f3
guard cospan:
  aa279eba2512513e3b1d004bbaea739db080b1d3320f3f0ac2f83b2e4b74e29a
  9884d64dc337b948009c4b438d23851ce7f56994addcd456f09d5cb262cee94b
clock/Hamilton audit:
  bd8ceec21ff77623852a58de865cf152ab66730c9358a79ee69d76ef7777cb0b
  53a438917a2ab141358531b2418ac1c892b6874aa82327543a64ca0d83af4191
```

An independent hostile audit replayed all three original normal executions,
checked their hashes, audited the rail and sector types, and caught the
disjoint-chart and half-compatibility boundaries now explicit in Sections 3
and 7.  It then independently derived the reversed counts in `(36)` before
they were added to the clock companion.  Normal and optimized executions of
that augmented companion byte-match its updated stored output and current
hashes above.  The computations exhaust their stated finite carrier and
include positive and hostile controls.  They do not construct the handoff
required after `(29)`; that is deliberately left as the next mathematical
obligation.

QED.
