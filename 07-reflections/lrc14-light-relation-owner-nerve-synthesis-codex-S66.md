# The current LRC14 object is a light-relation representation with an owner-obligation nerve

*codex-2026-07-17-S66.  Structural synthesis after THM-963/969/970/971/972/974/976/977/978.
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

## 4. A relation-rich / relation-poor dichotomy

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

## 5. A proposed recursive state

The global sporadic completeness bridge still needs a state that survives
legal insertions and deletions without losing labels.  A candidate state is

```text
S=(span L_14(v), sigma, N_sigma,
   labelled safe components and remaining rays,
   denominator debt n*m-q).
```

The last coordinate is positive and quantized on the sporadic branch;
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

## 6. Immediate proof targets

1. Formalize the mask-to-nerve bridges at scales eight through twelve.  The
   terminal `c=8` and `c=9` nerve contradictions are already kernel-pure;
   scale-ten and scale-twelve modules are in progress, while the large native
   reductions remain the honest trust surface.
2. Extend THM-972 from sum triples to a circuit-indexed counting lemma for
   all primitive weight-at-most-fourteen relations.
3. Compute `L_14(v)` and its circuit incidence on the eliminated `c=14,15`
   scalar banks and on the live `c=16` scale; use circuit rank to choose algebraic
   versus analytic dispatch before any metric recursion.
4. Instrument the existing metric insertion trees with the proposed state and
   search for an actual decreasing transition statistic.
5. Keep the global quantifier boundary explicit: closing common scales does
   not classify arbitrary non-AP-centred or deep sporadic packets.

The conceptual compression is now clear: the problem is not fundamentally a
tournament on runners.  It is a finite representation problem for light
integer relations whose obstruction is read in an owner-labelled
intersection nerve, with tournaments supplying useful but lossy telemetry.
