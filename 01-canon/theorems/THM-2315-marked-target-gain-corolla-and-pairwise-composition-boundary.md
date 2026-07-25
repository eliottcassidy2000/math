---
id: THM-2315
title: "The marked target-gain corolla and the pairwise composition boundary"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. THM-2309's
  two-target quotient is the projective line over F_13 with two labelled
  boundary axes. Its twelve mixed fork gains form one orbit under the
  labelled-boundary torus, so the terminal-word quotient has no equivariant
  fork section. Target swap sends gain g to g^-1; the fixed gains +/-1 force
  ties in every swap-equivariant binary head selector. The exact analytic
  object is a clocked measured word-corolla, or at full witness level a span,
  decorated by THM-2312's cubic current. Its pair marginals have a
  one-dimensional fork-versus-two-pure ambiguity; total source mass repairs
  masses but not witnesses or phase. Actual span pullback composes
  associatively, while positive prescribed-clock spans are not closed under
  composition. Thus neither a tournament, ordinary gain graph, signed
  hypergraph, underlying matroid, bare hafnian kernel, nor 1-skeleton
  2-Segal summary is faithful to every retained predicate. No scalar row is
  excluded and LRC(14) remains open.
source: codex-2026-07-25-target-gain-corolla
depends_on:
  - THM-2290-context-selected-colored-pair-kernel-is-hafnian-complete
  - THM-2294-anchored-plucker-tournament-and-kakeya-address-bank
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2310-quantitative-beta-shift-gluing-of-positive-handoffs
  - THM-2312-sparse-root-bispectrum-positive-word-current
related:
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2307-dual-rank-six-reconstruction-spectrum-and-selector-no-go
script: 04-computation/lrc14_target_gain_corolla_thm2315.py
output: 05-knowledge/results/lrc14_target_gain_corolla_thm2315.out
script_sha256: 99c398529307e9463b7197b62f2440ec3677c0b130278f0077f73452b74323ca
output_sha256: 82ef8f3b156ba52804690d7a1835c31a8649f2fe7422b3237df0c104163b528f
hash_basis: working-tree bytes (LF)
---

# THM-2315 -- the owner relation is a marked gain corolla, not a tournament

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2305 and THM-2309 expose two objects with the same three support types:

```text
analytic terminal words:         {a}, {b}, {a,b};
algebraic target-plane classes:  a-axis, b-axis, mixed support.          (1)
```

They are not yet identified. This theorem classifies each object, the exact
quotient relating their types, and the operations under which familiar
pairwise summaries fail.

The closest positive mechanisms are:

- THM-2305's canonical word partition and measured pure/fork alternative;
- THM-2309's exact target-plane quotient;
- THM-2312's nonzero cubic current on the same complete word;
- THM-2290's theorem that a pair kernel is exact only after all selector
  indices and its contraction contexts are retained.

The canonical hostile controls are a simultaneous fork versus two separate
pure events, target gains fixed by reversal, and disjoint arrival/departure
subsets with identical edge data.

## 1. The marked projective target line

Work over

```text
k=F_13,
V=k e_a direct_sum k e_b.                            (2)
```

THM-2309 identifies the owner-normalized quotient `K/L` with `V`. Its
nonzero projective classes form

```text
G=P(V)=P^1(k).
```

Write

```text
A=[1:0],                  the pure a-axis,
B=[0:1],                  the pure b-axis,
G_g=[1:g], g in k^*,      the mixed gains.           (3)
```

Define the support-word map

```text
s:G->{ {a},{b},{a,b} }

s(A)={a},       s(B)={b},       s(G_g)={a,b}.        (4)
```

The fibres have projective sizes

```text
1, 1, 12,                                             (5)
```

or vector sizes `12,12,144` before projectivization.

### Exact automorphism group

The projective-linear transformations preserving `A` and `B` individually
are exactly the diagonal transformations modulo common scalar. Thus

```text
T=PGL(V)_(A,B)=k^*,
c.G_g=G_(cg).                                        (6)
```

The mixed locus is one free transitive `T`-orbit. If target labels may be
swapped, the full boundary normalizer is

```text
N=N_(PGL(V))({A,B})=k^* semidirect C_2,              (7)
```

where the nontrivial element interchanges `A,B` and sends

```text
g -> g^(-1).                                         (8)
```

Indeed, a projective matrix fixing both coordinate axes is diagonal, while
one interchanging them is anti-diagonal. This proves (6)--(8).

The support map (4) is the orbit quotient for the labelled-boundary torus:

```text
G/T={A},{B},{mixed locus}.                            (9)
```

Consequently there is no `T`-equivariant section of `s` on the fork word.
Such a section would have to choose a mixed point fixed by every
`c in k^*`, but the action in (6) is free. Thus a terminal fork determines
the full twelve-point gain fibre and no canonical gain inside it.

This gives the sharp form of THM-2309's landing problem:

```text
word support can demand an address in the correct fibre of s;
word support alone cannot select its exact gain.                    (10)
```

Any exact-gain selection must use additional data that break this torus
symmetry. A root character, endpoint, shell address, or chosen coordinate
scale is only a candidate sidecar until an equivariant landing map using it
is proved.

## 2. Why the gain does not produce a tournament

A tempting repair is to use `g` to choose one head of a fork. Let

```text
h:k^*->{a,b,tie}
```

be compatible with target swap:

```text
h(g^(-1))=swap(h(g)),       swap(tie)=tie.            (11)
```

The inversion-fixed gains in `F_13^*` are exactly

```text
g=1,-1.                                               (12)
```

Equation (11) forces `h(1)=h(-1)=tie`. Therefore:

> **Head-selector no-go.** Every target-swap-equivariant binary summary of
> the mixed gains has at least two ties. No tie-free tournament head
> selector exists.

Away from the two fixed gains, the other ten gains form five inverse pairs.
There are exactly `2^5=32` equivariant selectors attaining the minimum of
two ties; each makes five arbitrary binary choices. None is canonical under
the full torus (6).

Quadratic character does not repair this:

```text
chi_13(g^(-1))=chi_13(g).                            (13)
```

It gives a symmetric six-square/six-nonsquare coloring of the fork gains.
Likewise, representing `G_g` by `(1,g)`, the alternating determinant of two
gains is

```text
det((1,g),(1,h))=h-g.
```

Reversal negates this determinant. Since `chi_13(-1)=1`, THM-2294's exact
boundary applies: the quadratic-character shadow is an undirected coloring,
not an orientation.

A signed hypergraph retaining `chi_13(g)` is therefore stronger than the
bare word but weaker than the gain line. For example `g=1` and `g=4` are
distinct mixed points with the same sign. The linear selector

```text
y-x
```

vanishes at `G_1` and not at `G_4`, so the lost gain is observable under
linear algebraic contexts.

The simple rank-two matroid on the fourteen projective points is just

```text
U_(2,14).                                             (14)
```

Its abstract automorphism group is the full symmetric group on fourteen
points; it can exchange a pure axis with a mixed gain and forgets every
cross-ratio. An oriented matroid would require additional lift/sign data
not supplied by the finite-field target plane. The faithful algebraic
carrier is therefore

```text
(P^1(F_13); labelled boundary A,B; exact mixed gain g),             (15)
```

not its underlying matroid or a chirotope tournament.

## 3. The analytic object is a measured word-corolla

Fix the THM-2305 source owner `j`, prescribed clock `k_j`, and targets
`a,b`. For a terminal word

```text
sigma in W_j={{a},{b},{a,b}},
```

let

```text
F_(j,sigma)
 =E_j intersection T^(-k_j)Q_(j,sigma).              (16)
```

Modulo null sets, these are disjoint transition witnesses and

```text
measure(F_(j,sigma))=rho_(j,sigma).                  (17)
```

Their disjoint union `F_j` carries the literal word map

```text
omega_j:F_j->W_j.                                    (18)
```

The exact support-level object is the clocked corolla

```text
j -> sigma,
```

with a pure edge for a singleton word and one irreducible two-head
hyperedge for `{a,b}`. Its measured quotient is

```text
mu_j=(p_a,p_b,d)
    =(rho_(j,{a}),rho_(j,{b}),rho_(j,{a,b}))
    in R_(>=0)^3.                                    (19)
```

THM-2305 retains the source, clock, word, and mass. THM-2312 further gives,
on every positive word, a nonzero gauge-invariant cubic phase-closure current

```text
C_(r,s)(Q_(j,sigma))                                 (20)
```

for some nonzero root-character pair. On one-sheet fibres this current can
be pure mass cubed, so (20) is not a theorem about an ordinary terminal
component phase. Thus the best current analytic summary is a **clocked
measured word-corolla decorated by a cubic closure current**. It still
forgets the individual transition witness and does not land the current on
a relation address in the gain fibre (10).

The canonical relation between the analytic and algebraic summaries is
therefore only the incidence correspondence

```text
I_j={
 (x,G) in F_j cross P^1(F_13):
 omega_j(x)=s(G)
}.                                                   (21)
```

No map `F_j->P^1(F_13)` is presently proved. Equation (21) records exactly
what a future landing theorem must refine.

## 4. Pair marginals lose the fork

The pairwise marginal map of the measured corolla is

```text
M:R^3->R^2,
M(p_a,p_b,d)=(p_a+d,p_b+d).                          (22)
```

It records how much mass hits each target, regardless of whether the two
hits occur simultaneously. Its kernel is

```text
ker M=span((-1,-1,1)).                               (23)
```

The minimal collision is

```text
(p_a,p_b,d)=(1,1,0)        two separate pure events,
(p_a,p_b,d)=(0,0,1)        one simultaneous fork,   (24)
```

both of which have pair marginals `(1,1)`.

Hence every ordinary directed graph, binary relation, or tournament built
only from target incidences identifies (24). The Boolean relation image of
one hyperedge `j->{a,b}` is literally the same pair set

```text
{(j,a),(j,b)}
```

as the image of the two pure edges. Relational or powerset composition reads
the fork as nondeterministic “a or b,” while the scalar word means the
conjunction “a and b at the same terminal point.”

There is one exact mass repair. If

```text
rho=p_a+p_b+d                                        (25)
```

is retained with the two marginals `m_a,m_b`, then

```text
p_a=rho-m_b,
p_b=rho-m_a,
d=m_a+m_b-rho.                                      (26)
```

Thus a weighted pair graph plus total source mass reconstructs the three
word **masses**. It still does not reconstruct the sets `F_(j,sigma)`, their
terminal components, or the cubic current (20).

## 5. Exact comparison with the Krenn pair kernel

Map the measured corolla to the symmetric star kernel on vertices
`{j,a,b}`:

```text
A_(j,a)=p_a+d,
A_(j,b)=p_b+d,
A_(a,b)=0.                                          (27)
```

This is an explicit instance of THM-2290's selected pair primitive. Add a
fourth vertex `z`. A context containing only the unit forced edge `{b,z}`
outside the star has hafnian response `A_(j,a)`; the analogous context with
forced edge `{a,z}` has response `A_(j,b)`. Hence forced-leaf contexts
reconstruct exactly the two marginals in (27).

They cannot recover more: the two corollas in (24) give the identical
kernel (27), so THM-2290's contextual completeness says every hafnian
extension depending only on that kernel agrees. The loss occurred before
hafnian contraction, in the map (22).

If the terminal word `sigma` is retained as an endpoint-color selector,
one may place `p_a,p_b,d` in three separate selected-kernel entries and
recover them by forced contexts. This is faithful because it manually keeps
the hyperedge selector; it is not a derivation of simultaneity from an
uncolored pair graph. No matching-family or hafnian identity is being
asserted for the LRC dynamics.

A still coarser attempt sends a mixed vector `(x_a,x_b)` to the single
target-pair weight `x_ax_b`. It identifies, for example,

```text
(1,2) and (3,5) in F_13^2,
```

because both products equal `2`, while their gains are `2` and `6`.
It also sends both pure axes to zero. Thus a bare pair weight cannot replace
the marked gain line.

## 6. Span composition and the association boundary

At witness level, a pure handoff is the span

```text
E_j <- F_(j,{t}) -> E_t,                             (28)
```

where the right map is `T^(k_j)`. Two such spans compose by fibre product
over the intermediate owner set. Span composition is associative up to the
canonical bijection of iterated fibre products. This is the faithful
set-level association law.

Positive prescribed-clock handoffs are not closed under this composition.
The obstruction is already finite. Let an intermediate probability space
have four atoms of mass `1/4`, and take all arrival and departure maps below
to be subset inclusions into it. In one model, incoming arrival and outgoing
departure supports are both `{1,2}`; in another they are `{1,2}` and
`{3,4}`. Decorate both models with the same labels, word types, and any
identical allowed gain data. Their incoming/outgoing masses are all `1/2`,
whereas their prescribed fibre products have masses

```text
1/2 and 0.                                           (29)
```

Thus no composition rule determined only by the hypergraph, gains, and edge
masses can preserve positive composability. Actual witness spans, or at
least their intermediate intersection, are necessary.

THM-2310 repairs support composition by inserting a sufficiently long
mixing wait. In span language it inserts a positive delayed correspondence
between the arrival support and the next departure support. The wait is
chosen from masses and variations, changes the prescribed clock, and is not
a canonical identity correspondence. It therefore does not turn the
original prescribed handoffs into a strict path category.

This is also the exact **one-skeleton reconstruction boundary** for a
2-Segal interpretation. A Hall/2-Segal multiplication determined from the
one-step corollas would have to assign the same two-step structure constant
to the two models in (29), while their actual pullback witnesses differ.
Retaining the full spans gives ordinary bicategorical associativity;
forgetting them to the one-skeleton does not determine the 2-simplex.
This is not a no-go for a witness-enriched 2-Segal object. THM-2310 supplies
existence after a noncanonical delay, not a 2-Segal identification of
prescribed decompositions.

Forks require one more enlargement. Their terminal set
`Q_(j,{a,b})` is not either exclusive-owner state `E_a` or `E_b`. To compose
a fork without turning conjunction into nondeterminism, the object set must
include the simultaneous state `{a,b}` and retain its common witness.
No outgoing transition theorem from that simultaneous state is presently
available.

## 7. Exact classification and stopping boundary

The faithful hierarchy is:

| retained object | exact predicate retained | first information lost |
|---|---|---|
| directed graph / binary relation | which target labels occur | fork versus two pure events |
| measured pair graph plus `rho` | three word masses | event witnesses and phase |
| directed measured hypergraph | pure/fork word and mass | intermediate intersection |
| marked gain line `(15)` | exact algebraic target residue | analytic word/current landing |
| clocked witness span | prescribed set-level composition | Fourier component phase |
| word-span plus cubic current | owner, clock, word, witness, gauge-invariant cubic phase closure | terminal-component and shell/address landing |

Accordingly:

1. an ordinary gain graph is too small: the pure axes `0,infinity` are not
   elements of the gain group `F_13^*`, and a fork is one hyperedge rather
   than a pair edge;
2. a signed hypergraph retains the word and perhaps square class, but
   collapses the twelve gains to two six-element classes and loses the
   complex current;
3. the underlying or oriented-matroid shadow needs boundary, field, and
   phase sidecars;
4. the marked projective **gain corolla** is faithful for THM-2309's
   algebraic quotient;
5. the clocked measured **word span with cubic current** is faithful for the
   analytic data now proved by THM-2305/2312.

The exact unresolved map is

```text
positive word-current
  -> one actual relation/Fourier address in the matching support fibre
  -> a shell edge with the required unit colour and component phase.    (30)
```

Neither tournament orientation nor abstract association can manufacture
this landing. No scalar row is excluded and LRC(14) remains open.

## 8. Exact companion

The companion enumerates the marked projective line, its torus and boundary
normalizer, the three support orbits, inversion-fixed gains, all minimal-tie
head selectors, quadratic-character colors, the marginal kernel and repair,
the forced-context hafnian probes, the pair-product collision, and the
composition witness. Reproduce with

```bash
python3 04-computation/lrc14_target_gain_corolla_thm2315.py
python3 -O 04-computation/lrc14_target_gain_corolla_thm2315.py
```

Both transcripts must match the stored output byte-for-byte after LF
normalization. QED.
