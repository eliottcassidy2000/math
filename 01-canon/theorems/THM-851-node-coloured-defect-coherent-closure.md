---
id: THM-851
title: One factorization deck closes the node-coloured defect algebra at n=3 and n=5..7
status: PROVED ABSTRACT CLOSURE LEMMA + FINITE-EXACT n=3..7
source: codex-2026-07-15-S15
depends_on: [THM-549, THM-781, THM-796, THM-818, THM-830]
related: [THM-813, THM-840, THM-843, HYP-6825, HYP-6880]
verification:
  - 04-computation/node_coloured_defect_factorization_codex_S15.py
  - 05-knowledge/results/node_coloured_defect_factorization_codex_S15.out
  - 05-knowledge/results/node_coloured_defect_factorization_codex_S15.json
---

# THM-851 - one factorization deck closes the node-coloured defect algebra

THM-830 identifies every fixed-path tiling with an arrow in a disjoint union
of elementary abelian pair groupoids.  This theorem gives the exact missing
refinement after colouring those arrows by tournament nodes.  At `n=3` and
for each `5<=n<=7`, one two-arrow factorization deck recovers every oriented
tiling; its reversal-compatible version recovers exactly every
grid-reflection orbit.  Pairing the resulting colours under all-tile
complementation then recovers the line quotient, including its node endpoints
and blue/black colour.  The analogous reconstruction fails at `n=4`.

## 1. The coloured pair groupoid

The computation starts on the unquotiented tiling groupoid `X_n`; it does not
start on THM-830's already complement-paired line quotient.  Write

```text
f=f_n=floor((n-1)/2),       r=h_(n-1)=floor((n-2)^2/4),
h=f+r,                      m=h+r=binom(n-1,2).
```

In each trace sheet, let `x in F_2^(f+r)` be an object and let
`delta in D=F_2^r` be a defect.  If `iota:D->F_2^(f+r)` pads the defect by
zero on the fixed trace, the arrow law is

```text
(x,alpha);(x+iota(alpha),beta)=(x,alpha+beta).             (1)
```

Let `tau(x,delta)` be the corresponding tiling and let `q(tau)` be its
converse-merged tournament node.  The initial node-coloured defect is

```text
C_0(x,delta)=(delta,q(tau(x,delta))).                      (2)
```

Grid reflection `sigma` reverses the arrow:

```text
sigma(x,delta)=(x+iota(delta),delta),
C_0(sigma(x,delta))=C_0(x,delta).                          (3)
```

Thus (2) preserves the complete defect word, not merely Hamming weight or the
binary blue/black predicate.  It still forgets where a coloured arrow sits in
its trace component.

For an output `a=(x,gamma)`, every `alpha in D` gives the unique
factorization

```text
a=(x,alpha);(x+iota(alpha),alpha+gamma).                   (4)
```

Define the ordered and reversal-compatible factor decks

```text
F_C^->(a)=multiset_(alpha in D)
 (C(x,alpha),C(x+iota(alpha),alpha+gamma)),                (5)

F_C^pm(a)=multiset_(alpha in D)
 {C(x,alpha),C(x+iota(alpha),alpha+gamma)}.                (6)
```

The braces in (6) are an unordered pair, but multiplicity and both complete
factor colours are retained.  Put

```text
R^->(C)=(C,F_C^->),             R^pm(C)=(C,F_C^pm).        (7)
```

## 2. General closure lemma

For any finite coloured groupoid, iteration of `R^->` stabilizes at the
coarsest refinement of `C` having constant ordered composition numbers.  If
the seed colour is inversion-invariant, `C(g)=C(g^(-1))`, iteration of
`R^pm` stabilizes at the coarsest inversion-invariant refinement having
constant symmetrized, or Jordan, composition numbers.

Indeed, (7) only splits colour cells, so it terminates.  A fixed point says
precisely that the number of factorizations with each ordered colour pair is
constant on every output colour; this is the intersection-number condition.
If an equitable colour `E` refines `C`, its constant `E`-pair counts can be
summed over the fibres `E->C`, so `E` refines `R^->(C)`.  Induction proves that
`E` refines every iterate and hence the stable partition.  The same argument,
with `(a,b)` and `(b,a)` summed, proves the symmetric statement.  For
off-diagonal colours the unordered count is `p_(ab)+p_(ba)`, i.e. the Jordan
coefficient before the conventional factor `1/2`.  This is the canonical
recursive closure requested in THM-830, not a chosen feature list.

## 3. Exact closure through n=7

The complete committed tiling-to-node atlases give:

| `n` | `(f,r)` | arrows/tilings | `C_0` cells | `R^->(C_0)` | `R^pm(C_0)` | `X_n/<sigma>` |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | `(1,0)` | 2 | 2 | 2 | 2 | 2 |
| 4 | `(1,1)` | 8 | 4 | 6 | 5 | 6 |
| 5 | `(2,2)` | 64 | 23 | **64** | **40** | 40 |
| 6 | `(2,4)` | 1,024 | 255 | **1,024** | **544** | 544 |
| 7 | `(3,6)` | 32,768 | 7,926 | **32,768** | **16,640** | 16,640 |

Consequently, for `5<=n<=7`, equality of ordered one-step colours is equality
of literal tilings, while equality of symmetric one-step colours is exactly

```text
t'=t or t'=sigma(t).                                      (8)
```

No second refinement round is needed.  At `n=4`, the stable ordered and
symmetric closures have respectively six and five cells, so the complete
separation statement genuinely starts at five.

The unrefined node colour is not composition-equitable.  It already fails in
the symmetric algebra at `n=4`.  At `n=5`, masks

```text
8  = {(5,2)},             24 = {(5,2),(4,2)}              (9)
```

have the same defect `01` and the same merged node `n5-a03`, but the factor
pair

```text
{(0,n5-a00),(01,n5-a03)}                                  (10)
```

occurs once for the first arrow and zero times for the second.  Thus exact
defect plus output node is not enough even to multiply colours consistently.

## 4. Precisely what must be preserved

When `gamma!=0`, the two distinct `alpha=0,gamma` terms in (5) first attach
the identity colours at the source and target objects.  For `gamma=0` they
coincide.  Those endpoint nodes explain much, but not all, of the refinement.
The controlled-forgetting census is:

| `n` | symmetric `C_0` | plus endpoint nodes | unlabelled midpoint-node deck | full defect-labelled deck | missing orbits before full deck |
|---:|---:|---:|---:|---:|---:|
| 4 | 4 | 4 | 5 | 5 | 1 |
| 5 | 23 | 29 | 39 | 40 | 1 |
| 6 | 255 | 401 | 543 | 544 | 1 |
| 7 | 7,926 | 14,562 | 16,604 | 16,640 | 36 |

Here the unlabelled deck retains the multiset of midpoint **node pairs** but
forgets which factor defect produced each pair.  At `n=5` it merges the two
blue masks `52` and `62`; at `n=6` it merges blue masks `32` and `352`.
Therefore the last necessary field is not another scalar node statistic.  It
is the incidence relation

```text
(factor defect) <-> (ordered or unordered pair of factor nodes). (11)
```

This is a small coloured triangle deck.  It is the recursive information that
node histograms, endpoint nodes, and unlabelled midpoint decks destroy.

## 5. Recovering tilings, nodes, and lines

Let `kappa` complement every tile and define, from the symmetric one-step
colour,

```text
Lambda([t])={R^pm(C_0)(t),R^pm(C_0)(kappa t)}.             (12)
```

Equation (12) is invariant under swapping the complement endpoints and under
grid reflection.  Exhaustively for `3<=n<=7`, it is injective on reflection
orbits of complement lines.  Their counts are

```text
n                         3    4    5     6      7
line orbits                1    3   20   272   8320.       (13)
```

There are `2^(m-1)` literal complement lines and `2^(h-1)` blue ones.  A black
line cannot be fixed setwise by reflection because a fixed grid coordinate is
unchanged by `sigma` and toggled by `kappa`.  Burnside therefore gives

```text
|L_n/<sigma>|=(2^(m-1)+2^(h-1))/2
             =2^(h-2)(2^r+1)                              (14)
```

when the exponent form is applicable; direct counting covers the tiny base
sizes.  Equations (8), (12), and (14) show exactly how the three tracked sorts
fit together:

```text
factor colour -> tiling modulo converse/grid reflection -> merged node,
complement pair of factor colours -> line modulo reflection -> blue/black. (15)
```

The literal line multiplicity is still a separate field.  Projecting (12) to
an unordered node pair can merge many line instances, as THM-843 shows at
`n=8`.

## 6. Tournament Analysis and scope

At `n=7`, use eight information carriers as Tournament Analysis vertices:
blue/black, merged node, node plus defect weight, node plus exact defect,
endpoint-node refinement, unlabelled midpoint deck, symmetric factor colour,
and ordered factor colour.  The pairwise observable is the number of unordered
arrow pairs separated.  The switch divides by the ceiling of the binary
address length; carrier name is the tie Hamiltonian path.

Both gauges are transitive, with score histogram `0,1,...,7`, no directed
triangle, eight singleton SCCs, and one Hamiltonian path.  They differ on 18
edges.  Thus information volume and information economy order these exact
carriers differently even though both comparisons are acyclic.

The challenged assumption is that vertices must be runners, tournament
vertices, or original arcs.  Here the faithful vertices are coloured arrows
and the proof obligations are their factorization triangles.  The theorem
preserves fixed-path tilings, merged tournament nodes, full defects,
composition multiplicity, reflection, complement incidence, and line colour.
It destroys runner speeds, metric scale, gap geometry, owners, wall chronology,
continued-fraction carry, and the LRC loneliness predicate.  It proves a
finite metagraph reconstruction theorem, not LRC(14).  Whether the one-step
separation in (8) survives at `n>=8` remains open. ∎
