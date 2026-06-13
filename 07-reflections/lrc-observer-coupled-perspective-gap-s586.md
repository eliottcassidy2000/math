---
source: codex-2026-06-03-S586
status: repo archaeology + exact tournament audit + LRC proof-stack synthesis
tags: [LRC, tournaments, perspectives, observer-coupled, A000568, HYP-2121]
---

# Observer-coupled perspective gap

The user brought back an old tournament curiosity:

```text
n=3 has 2 unlabelled tournament structures, but 4 perspectives:
  transitive contributes 3, cyclic contributes 1.

n=4 has 4 structures, and their perspective counts are 4+4+2+2=12.

12 is U(5), the number of unlabelled tournaments on 5 vertices.
```

The new prompt adds the missing interpretive key:

```text
observer-blind versus observer-coupled.
```

That distinction changes how I now read the whole thread.

## Repo Occurrence Map

I searched for literal phrases and nearby aliases:

```text
perspective(s), unique perspectives, distinct perspectives,
rooted/pointed tournament, vertex orbit, observer-pointed,
P(3)=4, P(4)=12, P(5)=48, 2+2+4+4,
56-48, gap=8, perspective conjecture.
```

The important occurrences are:

1. `03-artifacts/drafts/simplex-cube-perspectives-S25g.md`
   records the original wording almost exactly: transitive `3`, cyclic `1`,
   total `4`; four-vertex classes contribute `2+2+4+4=12`.  It defines
   `P(n)` as the sum of vertex-orbit counts.

2. `04-computation/simplex_perspective_counting.py` is the first exhaustive
   script for the small counts and the simplex/cube view.

3. `00-navigation/TANGENTS.md` T075 records the human Tournament Tiling
   Explorer conjecture: vertex-orbit counts at `n` predict unlabelled classes
   at `n+1`; it works for `3->4` and `4->5`, fails for `5->6`.

4. `03-artifacts/visualizations/tournament-tiling-explorer.html` embeds the
   same UI note: perspectives at `n=5` predict `48`, while actual `n=6` is
   `56`.

5. `03-artifacts/drafts/tiling-class-structure.md` leaves the failure as an
   open question: why does the vertex-orbit conjecture fail at `5->6`?

6. `04-computation/perspectives_investigation.py`,
   `04-computation/perspectives_deep_s92t.py`, and their result file
   `05-knowledge/results/perspectives_deep_s92t.out` establish:

   ```text
   P(2)=2=U(3)
   P(3)=4=U(4)
   P(4)=12=U(5)
   P(5)=48<U(6)=56
   P(6)=296<U(7)=456
   ```

7. `04-computation/perspective_class_map.py` and
   `05-knowledge/results/perspective_class_map.out` tried to map five-rooted
   perspectives to six-classes by adding a vertex with all incident words.
   The result file already contains the key correction, though it did not name
   it this way: all `56` six-classes have parents; there are `0` orphan
   classes.

8. `04-computation/perspective_A093934_check.py` notices an OEIS-name/offset
   confusion.  The count itself is robust: `P(n)` is pointed/rooted tournament
   count.  The external label is not what matters here.

9. `01-canon/theorems/THM-260-rooted-tournament-layer-decomposition.md` gives
   the exact Burnside layer:

   ```text
   P(n) = n*U(n) - D(n),
   ```

   where `D(n)` is an odd-partition automorphism correction.  This explains
   ordinary symmetry tax; it does not explain the LRC-safe predicate.

10. `07-reflections/lrc-tournament-56-bridge-s369.md` and
    `07-reflections/chirality-perspective-atlas-s370.md` tied the first
    perspective gap `56-48=8` to the eight LRC alpha stencils, the `12+44`
    self-converse/chiral split, and the `14+42` outer/interior LRC split.

11. `05-knowledge/hypotheses/HYP-1977-lrc-a000568-projection-defect.md` and
    `07-reflections/lrc-a000568-isoclass-analogy-s509.md` show the LRC version
    directly: unmarked classes and even observer-pointed classes can mix safe
    and unsafe cells.

12. Recent T600/T629 notes in `00-navigation/TANGENTS.md` already point toward
    the right quotient stack: raw circular body, source-deleted target,
    observer-source, threshold/gap/apex labels, kinetic flow, and
    blocker/pressure labels.

## What I Think We Misread

The old phrasing sometimes says:

```text
48 pointed five-tournaments but 56 six-tournaments,
so eight six-tournaments are not reachable.
```

That is false under the extension operation in `perspective_class_map.py`.
S586 recomputed the same thing:

```text
n=5->6:
  rooted=48
  target=56
  covered targets=56
  orphan targets=0
```

Even more important:

```text
root-sensitive parent classes=0.
```

The extension operation fixes an unrooted parent tournament, then ranges over
all `2^n` incident words for a new vertex.  The root variable is present in the
rooted list but not used by the extension rule.  So the operation is
observer-blind while being counted with observer-coupled objects.

That is exactly the kind of mistake the LRC quotient can make.

The small equality `P(n)=U(n+1)` is real, but it is not a natural bijection.
It is a rank-shift count:

```text
one level of pointing at n
looks numerically like
one level of unpointing at n+1
```

for a few small `n`.  At `5->6`, the equality fails by `8`; at `6->7`, it fails
by `160`.  The lesson is not "there are orphan next-level classes."  The lesson
is "pointing is not the extension payload."

## The Better Dictionary

For tournaments:

```text
unmarked class          = observer-blind structure
rooted/pointed class    = first observer-coupled lift
incident word           = how a new observer/vertex couples to every old vertex
threshold-coloured word = LRC-relevant incident word
endpoint-owner labels   = arithmetic payload over the incident word
```

For LRC:

```text
unmarked A000568 class       : useful cache, not exact safety
observer-pointed class       : necessary, still not enough
observer-source threshold    : exact one-observer loneliness predicate
endpoint/gap owner sheaf     : exact safe-box persistence predicate
proof automaton              : efficient "hit once" checker
```

The S509 projection-defect examples are now easier to understand.  If a safe
and unsafe cell share an observer-pointed half-turn class, the missing datum is
not another root.  It is the threshold/endpoint fiber above the root.

After rebasing, Opus S582/HYP-2120 landed with the sharper exact LRC slice:
source-perspectives satisfy `source(n)=T(n-1)` with no gap, because deleting a
source is canonical.  This reflection complements that result.  HYP-2120 says
which slice LRC should use for the clean skeleton; S586/HYP-2121 says why the
full rooted-perspective gap and the old add-a-vertex extension map were
misleading unless incident threshold fibers are retained.

## Which Clocks Matter

Safe to ignore for exactness, after using them as filters:

```text
raw unmarked tournament class,
raw Hamiltonian-path count H,
score sequence,
plain c3,
plain rooted class when the root does not change the switch,
the small P(n)=U(n+1) coincidence after n=4.
```

These are not useless.  They are cache keys, metagraph vertices, or entropy
barriers.  They just cannot certify "the orbit hit the safe box."

Must keep:

```text
which runner is observer,
observer incident threshold arcs,
whether observer is source/almost-source,
endpoint and gap labels,
owner residues,
pressure SCCs/good cuts,
middle automaton state,
Cprime/Phi/proof-obligation gates.
```

THM-381 and THM-385 are the reason observer-source data matters: loneliness of
the observer is not a metaphor there; it is exactly the source predicate in the
threshold tournament.

## A More Efficient Safe-Box Model

Instead of asking:

```text
sample time until a runner is lonely
```

ask:

```text
does the periodic/equidistributed orbit hit an accepting state in a labelled
quotient automaton?
```

State stack:

```text
base:       coarse quotient class / wall-crossing chamber
root:       observer-pointed version of that class
switches:   thresholded observer incident arcs
fiber:      endpoint/gap/owner labels and pressure SCCs
automaton:  L/M/R or proof-obligation state
accepting:  observer-source with compatible endpoint labels
```

This makes the question finite in the rational-ratio case and measure/discrepancy
driven in the irrational/equidistributed case.  It also clarifies what should
be cached: base classes and rooted classes are cache layers, while threshold
and endpoint labels are the accepting-state payload.

## Connection Back to the Eight

S369/S370 were right to care about:

```text
56 = 48 + 8,
56 = 12 + 44,
56 = 14 + 42.
```

But I would now phrase the shared `8` as:

```text
the first observer-coupled extension residue.
```

In tournaments, five-rooted perspectives do not supply enough data to account
for all six-unrooted structures by a one-per-perspective count.

In LRC, the eight alpha stencils are not just "missing classes"; they are the
places where a coarse quotient needs endpoint/threshold payload.  Multiplication
by the seven odd shifts is then the arithmetic orbit layer over that stencil
fiber.

This is also compatible with the recent automaton work.  A binary edge is often
the terminal shadow of a three-state process.  The perspective gap is another
place where the middle state is hidden: the incident word has been summed out,
leaving only a root count.

## Tournament Analysis

S586 used proof-state models as tournament vertices, explicitly challenging
the assumption that tournament vertices must be runners:

```text
unmarked_A000568_class
observer_pointed_class
root_plus_incident_word
observer_source_threshold
endpoint_owner_sheaf
proof_obligation_automaton
```

Pairwise observable:

```text
predicate preservation,
compression,
computability,
proof power,
extension payload,
projection risk.
```

The resulting tournament is transitive:

```text
endpoint_owner_sheaf
  > observer_source_threshold
  > proof_obligation_automaton
  > root_plus_incident_word
  > observer_pointed_class
  > unmarked_A000568_class
```

Metrics:

```text
score_histogram = {0:1,1:1,2:1,3:1,4:1,5:1}
directed_3_cycles = 0
scc_sizes = [1,1,1,1,1,1]
Hamiltonian paths = 1
```

This is not meant as a surprising tournament.  It is a clean priority order:
keep the endpoint-owner sheaf and observer-source threshold state; use the
unmarked class only as a base cache.

## New Hypothesis

This reflection creates HYP-2121:

```text
the perspective gap is an observer-coupling defect, not an orphan-class defect.
```

The next concrete test is to revisit the S509 mixed pointed buckets and add
threshold incident words first, endpoint-owner labels second.  If this framing
is right, most safe/unsafe mixing should split at one of those two lifts.  The
remaining mixed cases are the genuine middle automaton/proof-gate residuals.
