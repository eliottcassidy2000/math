---
id: THM-2656
title: "Endpoint-anchored rainbow selector moduli and reflection fixed-edge tax"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT AUDIT.  Once one
  absolute endpoint two-set is supplied, THM-2647 reconstructs the other and
  the THM-2648 rainbow selectors split into genuinely different moduli.  An
  orientation-independent atlas of two individually reflection-equivariant
  rainbow charts exists on every eleven-by-eleven relation; both charts must
  contain the canonical midpoint edge, and explicit matched-wall templates
  make this their only common edge, retaining 21 of 121 edges.  A separate
  twelve-template bank gives two edge-disjoint rainbow charts on every
  relation and retains 22 edges, but cannot keep both charts individually
  reflection-equivariant.  Both atlases cover all thirteen carries, have
  profile 1^4 2^9, retain all twelve charged characters, and have energy
  36/169.  Edge overlap is selector data, not a relation invariant; no
  physical endpoint sidecar or lawful LRC selector is constructed.
source: deep-energy-audit-2026-07-28-rainbow-selector-moduli
depends_on:
  - THM-2647-endpoint-anchored-two-point-deconvolution-and-thirteen-halves-signed-tax
  - THM-2648-two-rainbow-thinning-full-carry-cover-and-noncanonical-selector-boundary
related:
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - HYP-2233-missed-problem-frontier-carrier-atlas
script: 04-computation/lrc14_endpoint_anchored_rainbow_selector_moduli_thm2656.py
output: 05-knowledge/results/lrc14_endpoint_anchored_rainbow_selector_moduli_thm2656.out
script_sha256: 5f6c79afc3aab536e4c18f47cd6e5b9bf02256c2e6424a88ba567980137d42cc
output_sha256: 23294c0a5cedc9ab3a8de5178fea76367a3cc1ba4a6053489add8cf86af5facd
hash_basis: LF-normalized bytes
---

# THM-2656 -- the midpoint is a symmetry tax, not an overlap invariant

**PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT AUDIT.**

THM-2648 produced two rainbow charts from a saturated eleven-by-eleven
relation, but its one/three-edge overlap census belonged to the displayed
selector atlas.  It was not intrinsic to the relation.  THM-2647 now exposes
the correct organizing coordinate: after one absolute endpoint two-set is
supplied, the other endpoint is reconstructed, and selector space separates
into a reflection-fixed lane and a symmetry-breaking lane.

## 1. Endpoint reconstruction and the intrinsic reflection

Let

```text
G=F_13,        A={a_0,a_1},        B={b_0,b_1},
S=G\A,         T=G\B,
m=1_A*1_B.                                                (1)
```

Suppose the full absolute left endpoint `A` is supplied in the common-origin
gauge.  THM-2647 reconstructs `B` exactly from `m`; no ordering of the two
members is required.  Put

```text
mu_A=(a_0+a_1)/2,       mu_B=(b_0+b_1)/2,                (2)
rho_A(x)=2mu_A-x,       rho_B(y)=2mu_B-y.                (3)
```

These definitions depend only on the unordered sets.  The reflections
preserve `S,T`, and `mu_A in S`, `mu_B in T` are their unique fixed points.

A matching `f:S->T` is **individually reflection-equivariant** when

```text
f rho_A = rho_B f.                                       (4)
```

Equation (4) forces

```text
f(mu_A)=mu_B.                                             (5)
```

Consequently any two individually equivariant matchings share the canonical
midpoint edge `(mu_A,mu_B)`.  In particular their union has at most `21`
edges.  This fixed-point argument is the exact symmetry tax.

## 2. An intrinsic symmetric atlas attaining the tax

Choose temporary orderings and normalize by

```text
x=a_0+(a_1-a_0)u,
y=b_0+(a_1-a_0)v,
t=(b_1-b_0)/(a_1-a_0).                                  (6)
```

Then the missing sets are `{0,1}` and `{0,t}`, and both reflections are

```text
u -> 1-u,             v -> t-v.                          (7)
```

Away from `t=+-1`, take the two affine charts

```text
f_+(u)=tu,             f_-(u)=t(1-u).                    (8)
```

They are rainbow, their hole pairs are respectively

```text
{0,1+t},               {t,1},                            (9)
```

and these are disjoint.  Both charts obey (7), and solving `f_+=f_-` leaves
only `u=1/2`; hence their unique common edge is the midpoint edge.

On the two matched walls, retain the unique rainbow affine orientation and
use the following nonlinear second rows.  Entries are the targets for
`u=2,3,...,12`:

```text
t= 1: (3,4,2,6,9,7,5,8,12,10,11),   holes {3,12};       (10)
t=-1: (2,3,1,5,8,6,4,7,11,9,10),    holes {2,11}.       (11)
```

For `t=1` the first chart is `v=u` with holes `{0,2}`.  For `t=-1` it is
`v=u-1` with holes `{1,12}`.  Each row in (10)--(11) is a target bijection,
is rainbow, commutes with (7), has the displayed disjoint holes, and meets
the affine chart only at `u=1/2`.

The construction is independent of the temporary orderings.  Off the walls,
changing an endpoint orientation merely permutes the two affine charts.  On
the walls, if `g_+,g_-` denote the rows (10)--(11), then

```text
g_-(u)=g_+(u)-1=-g_+(1-u),                               (12)
```

which is exactly the source/target orientation-change law.  Thus the
unordered pair of charts descends to the unordered endpoint sets recovered
by THM-2647.

Transporting (8), (10), and (11) through (6) proves:

> Every endpoint-anchored relation has an orientation-independent pair of
> individually reflection-equivariant rainbow charts with disjoint colour
> holes and exactly the midpoint edge in common.

The atlas retains exactly `21` of the `121` relation edges on every one of
the `6,084` endpoint pairs.  Its chart census is the familiar

```text
11,154 affine,              1,014 nonlinear.             (13)
```

## 3. A maximal edge-disjoint atlas

The midpoint cost disappears as soon as individual chart equivariance is
dropped.  Keep the first normalized chart

```text
f_1(u)=tu       (t!=-1),              f_1(u)=u-1 (t=-1), (14)
```

and, for each `t=1,...,12`, use the corresponding second row below.  Again
the entries are targets for `u=2,...,12`:

```text
 1: (11,12, 5, 2, 4, 9, 6, 3, 8,10, 7)
 2: ( 1, 4, 9, 3, 5, 7,11, 6, 8,12,10)
 3: ( 2,10, 5,11, 1, 7, 4, 6, 8,12, 9)
 4: ( 3,10,11, 9, 1, 5, 8, 2,12, 6, 7)
 5: ( 4,10, 1, 7, 8, 2, 3,12, 6, 9,11)
 6: ( 5, 8, 9, 1,11, 7, 2, 3,12,10, 4)
 7: ( 6,10, 1,11, 9, 3,12, 2, 4, 8, 5)
 8: ( 7,10, 1, 9, 2, 3, 4,11, 5, 6,12)
 9: ( 8,10,11, 4, 1, 5, 3,12, 6, 7, 2)
10: ( 9, 1, 3, 8, 2, 7, 4, 6,12, 5,11)
11: (10, 1, 6, 8, 9, 2, 3,12, 4, 5, 7)
12: (10,11, 4, 1, 3, 8, 5, 2, 7, 9, 6).               (15)
```

Each row is a bijection from `G\{0,1}` to `G\{0,t}`, has eleven distinct
carry colours, has holes disjoint from those of (14), and uses no edge of
(14).  Affine transport therefore gives two edge-disjoint rainbow charts on
every endpoint pair.  Their union has `22` edges, the maximum possible for
two eleven-edge matchings.

Both charts cannot satisfy (4), because (5) would give them a common edge.
The explicit bank (15) makes a temporary endpoint ordering; the theorem does
not claim that every setwise-equivariant, chart-swapping atlas has been
classified.  The proved contrast concerns **individual** chart symmetry.

## 4. The carry spectrum does not see the selector component

For either atlas, two disjoint two-point hole sets give the same carry
incidence profile

```text
w=2-1_(H_0 union H_1),             w=1^4 2^9.            (16)
```

Thus all twelve charged characters survive by the `Phi_13` rational
criterion, and normalized Parseval gives

```text
sum_(k!=0)|what(k)|^2=36/169,
max_(k!=0)|what(k)|^2 >= 3/169.                           (17)
```

The symmetry-fixed atlas and the edge-disjoint atlas are indistinguishable
by (16)--(17), even though their physical edge unions have sizes `21` and
`22`.  Edge overlap is therefore selector data, not an invariant of the
underlying relation or its carry spectrum.  In particular, THM-2648's
one/three-edge census remains correct for its displayed atlas but is not an
optimal or intrinsic overlap law.

## 5. Holotopy interpretation and physical boundary

The reflection-fixed selector lane behaves like a fixed stratum of a section
space: the unique fixed source and target force an intersection number one.
After forgetting individual equivariance, the selector moves to a free
stratum where the two sections can be disjoint.  The carry projection
contracts these strata to the same `1^4 2^9` spectral point.  This is a small
but exact holotopy model of information lost by quotienting a relation to its
character profile.

THM-2647 supplies the minimal **algebraic consumer** once an endpoint anchor
exists.  It does not supply that endpoint on the physical LRC packet.
Likewise, neither atlas turns the static eleven-sheet coefficient fibre into
a same-base positive transition table, and neither provides a lawful
measurable selector preserving the owner, word, clock, root, and endpoint
current.  The extra ordered normalization used by (15) is also not a
canonical physical datum.

Therefore THM-2656 proves no row exclusion.  The LRC ledger remains `165`.

## 6. Exact reproduction

Run

```bash
python 04-computation/lrc14_endpoint_anchored_rainbow_selector_moduli_thm2656.py
python -O 04-computation/lrc14_endpoint_anchored_rainbow_selector_moduli_thm2656.py
```

The dependency-free exact companion checks all `6,084` endpoint pairs and
all four temporary endpoint orientations.  It independently replays the
THM-2647 deconvolution, verifies every symmetric chart, reflection identity,
midpoint intersection, nonlinear wall transport, every row of the disjoint
bank, both edge-union counts, all `146,016` charged-character guards, and both
copies of (16)--(17).  Normal and optimized runs must byte-match the stored
transcript and end in `PASS`.

QED (candidate; independent audit pending).
