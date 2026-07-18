# The current LRC14 object is a light-relation representation with an owner-obligation nerve

*codex-2026-07-17-S66.  Structural synthesis after THM-963/969/970/971/972/974/976/977/978/980/981.
This is a research map, not a new global theorem.*

## 1. Two exact objects have finally met

The recent calculations expose two objects that had mostly been studied in
separate lanes.

For a speed packet `v=(v_i)`, define its weight-fourteen relation system

```text
L_14(v)={alpha in Z^I : sum_i alpha_i v_i=0,
                         sum_i |alpha_i|<=14}.
```

THM-972 and `LRCRealRelationLock.lean` say that at a common bad phase, after
choosing the corresponding integer witnesses `w_i`, every relation in
`L_14(v)` is inherited exactly:

```text
sum_i alpha_i w_i=0.
```

This includes the rational pair locks, all sum triples, and the first exact
mediant triple counts.  The witness vector is therefore not merely close to
the speed ray.  It is a representation of the same light integer-relation
system.

For one fixed sheet context `sigma`, let `U_sigma` be the finite set of global
unit words and let

```text
O_o(sigma)={u in U_sigma : u covers every sheet at owner o}.
```

The correct terminal object is the simplicial intersection nerve

```text
N_sigma={F subset Owners : intersection_(o in F) O_o(sigma) is nonempty}.
```

A packet survives exactly when the full owner set is a face.  Minimal
nonfaces are human-sized unsatisfiability certificates even when the raw bank
has hundreds of millions of labelled rows.

The emerging underlying object is thus the pair

```text
(light relation system L_14(v), owner-obligation nerve N_sigma),
```

with CRT masks giving a representation functor from the first object to the
second.  The computational bottleneck is increasingly this functor, not the
terminal contradiction.

## 2. Standardizing the graph language

The repository has sometimes called both a pair-intersection graph and its
complement a “nerve.”  The distinction matters.  In this note:

- the **intersection graph** is the one-skeleton of `N_sigma`;
- the **conflict graph** joins owners whose obligation sets are disjoint;
- the full nerve also records triple and higher intersections.

With that convention, the exact scale sequence is:

| scale | first decisive carrier | terminal obstruction |
|---|---|---|
| `5` | owner obligations | eight minimal triple nonfaces; pair data alone is insufficient |
| `6` | affine four-flats over `F_2` | three antipodal conflict edges plus no triple faces |
| `7` | row-character matrix over `F_7` | the two quadratic cosets die before an owner nerve is needed |
| `8` | signed-cycle owner obligations | distance-two conflict graph `K3,3` |
| `9` | divisor-coloured owner obligations | intersection graph `3K2`, or mixed induced graph `2K2`; every cross-pair is a conflict |
| `10` | projective owner classes `F_13^*/{+-1}` | conflict cycle `C6`; intersection graph `K6\C6`, the triangular prism |
| `11` | owner obligations on 66 all-prime-order supports | 64 supports have empty intersection graph; the two quadratic cosets have intersection graph `3K2` |
| `12` | owner obligations on 64 all-order-twelve sign transversals | complete owner orthogonality: the intersection graph is empty and every obligation has size 48 |
| `13` | no primitive common-scale face | THM-860: every retained and replacement speed would be divisible by 13 |
| `14` | fourteen local sheets at one owner | THM-977: scalar capacity leaves 576 rows, but every owner-local union misses at least two sheets; no global nerve remains |
| `15` | labelled owner-feasibility subset and maximum-union vector | THM-978: 2,184 scalar rows have 0, 1, 2, or 4 feasible owners, never all six |
| `16` | labelled owner-feasibility subset and maximum-union vector | THM-980: 2,540 scalar rows have 0, 1, or 2 feasible owners, so at least four owners fail |
| `17` | alternating signed-projective support classes, then owner-local masks | THM-981: scalar capacity leaves the QR/NQR rows; all twelve owner maxima are exactly sixteen of seventeen sheets |
| `18` | labelled owner-feasibility subset and maximum-union vector | THM-982: 13,098 scalar rows have 0, 1, 2, or 4 feasible owners, so at least two owners fail |
| prime `p>=19` | weighted ratio-capacity digraph before any sheet nerve | THM-983: residue-block recurrence plus two exceptional support scans kills scalar capacity uniformly |
| `20` | labelled capacity and owner-feasibility/max-union vectors | THM-986 claim: 12,584 scalar rows have 0, 1, or 2 feasible owners; independent replay is pending |

The carrier changes because each quotient is legal only after checking which
predicate it preserves.  At scale ten, projectivization is legal because the
owner-pair law factors through sign classes.  At scale eleven, the quadratic
cosets are exceptional precisely because negation pairs retain identical
obligations.

## 3. What the tournament viewpoints did and did not reveal

Several tournament constructions remain useful diagnostics:

1. **Runner vertices.** Orient by speed, residue, or boundary arrival.  This
   orders a search but forgets simultaneous sheet coverage.
2. **Provider or sheet vertices.** Orient by which mask covers more sheets.
   This retains local capacity but forgets that one unit word is shared.
3. **Owner-obligation vertices.** Use nonempty pair intersection as the
   observable, choose a sign/projective gauge when available, and complete
   ties along a deterministic Hamiltonian path.  This is closest to the proof
   object but still forgets intersection multiplicities and higher faces.
4. **Relation-circuit vertices.** Orient two light circuits by first forced
   witness coordinate or coefficient weight.  This sees algebraic locking,
   but without the mask representation it cannot decide coverage.
5. **Boundary-event vertices.** Orient by first wall crossing.  This retains
   continuation order but needs owner and component labels to be sound under
   recursion.

At scale eleven the completed tournaments for the generic empty graph and the
exceptional `3K2` graph can both be made transitive with the same score
sequence, SCC data, and Hamiltonian-path count.  The tournament is therefore
provably nonfaithful even at the coarsest distinction between the two terminal
species.  The missing data is not a more clever orientation: it is the edge
meaning and the higher nerve.

Scale fourteen sharpens this warning.  THM-977's independent certificates
show that the contradiction occurs before owner
obligations exist as global sets: the faithful vertices are the fourteen
sheets of one owner and the invariant is maximum union cardinality.  Any
tournament on runners or owners literally quotients away the deficit.

Scale fifteen is the transitional case.  Some owner obligations are locally
nonempty, but every context has at least two empty owners.  Ranking owners by
`(feasible,maximum union)` always produces a transitive tournament, while the
proof is the absolute threshold statement encoded by the labelled feasibility
subset.  This identifies a new pre-nerve layer between local sheet incidence
and the global obligation nerve.

Scale sixteen strengthens that layer: two independent exact reconstructions have
no row with more than two locally feasible owners.  Its completed owner
tournament is again always transitive.  Keeping the exact ordered pair
observable would let one reconstruct the labelled maximum-union vector, but
the orientation alone forgets the deficit magnitudes and the absolute
sixteen-sheet threshold.  A tournament can serialize the obstruction; it
cannot replace it.

Scale seventeen separates the two roles especially cleanly.  Signed
projective classes classify the scalar bank: the only multiplicity words are
the two alternating `2,0,2,0,2,0` color classes.  That quotient cannot close
the problem, because the terminal fact is one level lower in the mask
representation: every owner reaches sixteen sheets but not seventeen.  Since
all six owner summaries tie, even the completed tournament has literally no
edge signal; its deterministic Hamiltonian path is pure gauge.

Scale eighteen shows that the pre-nerve obstruction persists even when the
divisor lattice and scalar bank become much richer: 63 multiplicity patterns
and 13,098 scalar rows survive, yet no row has more than four nonempty owner
projections.  Scale twenty's primary census strengthens the same phenomenon
to at most two owners despite 830 of 924 supports surviving the scalar gate.
This is evidence that the robust quantity is **projection debt**, not scarcity
of scalar supports.  The scale-twenty statement remains a claim until its
independent referee lands.

Prime scales expose a different earliest carrier.  THM-983 never needs sheet
incidence: `a_(p+13)(r)=a_p(r)+2` turns the problem into weighted ratio
capacity.  The two numerical exceptions also collapse structurally: `p=23`
is killed by the symmetrized low-ratio Cayley spectrum (`9>6` internal-edge
debt), and `p=29` by multiplicative closure of `{1,4,5,8,9}` to all of
`F_13^*`.  Ordering weights gives a transitive tournament, whereas orienting
reciprocal capacities gives a strongly connected tournament with 40 triangles
at both exceptions.  Their disagreement is informative: the faithful object
is the weighted directed ratio graph plus labelled six-support and absolute
owner threshold, not either tournament completion.

## 4. The pre-nerve filtration is a projection/gluing problem

The earlier two-object description should be refined by one layer.  Fix a
scalar context `sigma`, and for each owner `o` let

```text
R_o(sigma)={ union_i M_{i,o}(e_i) : e_i ranges over the locally allowed
                                         unit states of provider i },
m_o(sigma)=max{|A|:A in R_o(sigma)},
F(sigma)={o:m_o(sigma)=c}.
```

Here the choices `e_i` are allowed to depend on `o`.  Thus `F(sigma)=Owners`
is necessary but not sufficient for one global unit word: it tests the six
owner projections separately and deliberately forgets covariance between
them.  If some owner is absent from `F(sigma)`, however, the global fibre is
immediately empty.  This is exactly what closes scales fourteen through
eighteen, and is the claimed scale-twenty obstruction.

Only after all projections are nonempty should one form the global product
`U_sigma` of covariant unit words and the inverse-image obligations

```text
O_o(sigma)={u in U_sigma : the o-projection of u is the full c-sheet mask}.
```

The resulting proof object is a filtration:

```text
scalar capacity
    -> owner-local reachable families (R_o,m_o,F)
    -> globally covariant unit fibre U_sigma
    -> obligation nerve N_sigma
    -> metric continuation.
```

Scale fourteen dies in the second layer at every owner; scales fifteen and
sixteen die there because `F(sigma)` is always a proper labelled subset.
Scales eight through twelve reach the nerve layer, where pairwise conflicts or
higher nonfaces close them.  This explains why one universal graph was never
quite the underlying object: the honest carrier is a finite projection/gluing
system, and the nerve is its derived compatibility invariant only after local
nonemptiness.

Categorically flavored language is useful here but should not be oversold.
The local reachable families are projections of a finite global word fibre;
the proof asks whether six prescribed full-mask sections glue to one word.
This wording preserves the crucial quantifier order and suggests reusable
algorithms: compute projection images, delete empty obligations, then build
only the surviving part of the nerve.  It also identifies the smallest
unsatisfiable certificate: an empty projection, a minimal nerve nonface, or a
later metric obstruction, in that order.

## 5. A relation-rich / relation-poor dichotomy

The light-relation theorem suggests a principled split before any large
metric recursion.

### Relation-rich packets

If the span of `L_14(v)` has high rank, then the witness vector has few free
coordinates.  Sum-triple and mediant circuits should be eliminated first,
followed by the induced CRT representation on the remaining parameters.  The
canonical pair layer is now closed in Lean (THM-964/967/971), and THM-972
opens the triple layer.  The natural next invariant is the circuit matroid of
`L_14(v)`, not a list of isolated pair ratios.

### Relation-poor packets

If `L_14(v)` has low rank, the packet is genuinely dissociated at the scale
seen by the band inequality.  This is exactly the regime where Hunter credit,
pair-overlap arcs, Bonferroni deviation, and density-floor estimates have a
chance to be uniform.  The negative canonical-family B5 limit is then no
paradox: that family is relation-rich and belongs in the locked/direct lane,
not in a dissociated live-floor theorem.

This dichotomy aligns two formerly competing proof programs.  Algebraic locks
handle resonance; analytic overlap handles the absence of light resonance.

## 6. A proposed recursive state

The global sporadic completeness bridge still needs a state that survives
legal insertions and deletions without losing labels.  A candidate state is

```text
S=(span L_14(v), sigma, (R_o,m_o,F)_o, U_sigma, N_sigma,
   labelled safe components and remaining rays,
   denominator debt n*m-q).
```

The pre-nerve coordinates retain the quantifier order that was invisible in
the earlier state: local unit choices may vary with the owner, whereas a point
of `U_sigma` is one shared choice.  The last coordinate is positive and
quantized on the sporadic branch;
`LRCSporadicDiscreteCap.lean` turns it into the improved top-speed cap.  The
nerve records exact simultaneous compatibility.  The component/ray payload
records continuation semantics.  No single coordinate is currently known to
be monotone, so the next experiment should log transitions of all five under
the existing Hamming insertion trees and test lexicographic or multiset
potentials.

The most promising well-founded alternatives are:

- rank defect of the inherited light-relation system;
- number and dimension of unresolved nerve faces;
- number of labelled continuation components;
- denominator debt `n*m-q` together with the exact completion cap.

A valid completeness theorem would show that every transition decreases one
coordinate without increasing any earlier coordinate, or lands in a closed
scale/ramification language.  This is a concrete falsifiable program, not a
claim that such an order has already been proved.

## 7. Immediate proof targets

1. Formalize the mask-to-nerve bridges at scales eight through twelve.  The
   terminal `c=8` and `c=9` nerve contradictions are already kernel-pure;
   scale-ten and scale-twelve modules are in progress, while the large native
   reductions remain the honest trust surface.
2. Extend THM-972 from sum triples to a circuit-indexed counting lemma for
   all primitive weight-at-most-fourteen relations.
3. Compute `L_14(v)` and its circuit incidence on the eliminated
   `c=14,15,16,17,18` scalar banks and the provisional `c=20` bank; use circuit rank to
   choose algebraic versus analytic dispatch before any metric recursion.
4. Instrument the existing metric insertion trees with the proposed state and
   search for an actual decreasing transition statistic.
5. Keep the global quantifier boundary explicit: closing common scales does
   not classify arbitrary non-AP-centred or deep sporadic packets.

The conceptual compression is now clearer: the problem is not fundamentally
a tournament on runners.  It is a finite projection/gluing problem for
representations of light integer relations.  Its earliest obstruction is an
empty owner projection; if those projections survive, the obstruction is read
in an owner-labelled intersection nerve; tournaments supply useful but lossy
telemetry at either layer.
