# HYP-2570 - LRC section check-off should reduce exotic cases to wall-switch Hall packets

**Status:** OPEN proof program  
**Source:** codex-2026-06-17-S3  
**Computation:** `04-computation/lrc_section_checkoff_switches_codex.py`  
**Stored output:** `05-knowledge/results/lrc_section_checkoff_switches_codex.out`

## Claim

Fix the slowest-runner gauge for an `n`-runner LRC instance by writing the
speeds as

`V=(0,v_1,...,v_{n-1})`.

Divide the circle into the `n` fixed sections

`[s/n,(s+1)/n)`.

For each runner `r`, record the sections in which `r` has a lonely witness,
meaning

`||(v_r-v_j)t|| >= 1/n` for every `j != r`.

This gives a bipartite graph

`runners <-> sections`.

The user's check-off dream is exactly the existence of a perfect matching in
this graph.  The proof program is:

1. Work with the compactified graph, where section boundary witnesses count for
   both adjacent sections.
2. Prove that any failure of strict-open matching is represented by a small
   Hall-deficit packet plus a finite wall-debt list.
3. Prove a wall-switch lemma: every open Hall packet either switches across a
   boundary to enlarge its section union, or descends to an endpoint-mouth /
   observer-source case already handled by the LRC14 dihedral machinery.

In this form, exotic cases are not arbitrary.  They compress to:

- a runner subset `R`,
- its available open section union `S`,
- the deficit `|R|-|S|`,
- the wall-debt sections reachable only at equality,
- a section-support tournament on the contested sections.

## Evidence

The script performs exact rational open-cell plus wall-event scans in the
canonical gauge.

Small primitive scans:

| total runners `n` | speed max | primitive rows | compact perfect matching | strict-open matching |
|---:|---:|---:|---:|---:|
| 4 | 14 | 325 | 325/325 | 0/325 |
| 5 | 10 | 205 | 205/205 | 0/205 |
| 6 | 8 | 56 | 56/56 | 0/56 |

Thus the compactified check-off dream is ubiquitous in the tested range, while
strict-open check-off always exposes wall debt.  That is not a failure of the
section idea; it identifies the missing lemma.

Named examples:

- LRC4 AP `(0,1,2,3)` has compact matching and open cover `2/4`; runners
  `(0,3)` are wall-only, giving the open Hall packet
  `deficit=2, R=(0,3), S=()`.
- LRC4 uneven `(0,1,2,4)` has compact matching and open cover `3/4`; the open
  Hall packet is `deficit=1, R=(1,2,3), S=(1,2)`, with wall debt in section `3`.
- LRC14 AP and the Goddyn-Wong row both have compact matching and open cover
  `12/14`; endpoint sections `0` and `13` carry the visible wall debt.

## Tournament Analysis

This hypothesis deliberately challenges the default vertex choice.

- Vertices: fixed loop sections, not runners.
- Pairwise observable: section support vector
  `(witnessed runner count, total open witness measure, private runners)`.
- Switch/gauge: orient `a -> b` when section `a` has larger support vector.
- Tie Hamiltonian path: increasing section index.
- Fingerprints: score histograms, directed 3-cycles, SCCs, Hamiltonian-path
  counts, and edge flips under the alternate private-first switch.

The tested examples are mostly transitive section ledgers, which means raw
section support is too scalar by itself.  The proof-bearing tournament is likely
the smaller tournament induced on Hall packets, wall-debt sections, or boundary
switch events.

## Assumption Challenge

Considered vertex sets:

- runners,
- gaps,
- fixed circle sections,
- section boundaries,
- wall-crossing events,
- residues,
- cover arcs,
- Fourier modes,
- matroid-like circuits,
- proof obligations.

Chosen quotient: fixed sections plus Hall packets.  It preserves the LRC
predicate "which regions can host a runner's lonely moment" and the matching
obstruction to checking off all regions.  It destroys timing, most pairwise
runner identity, and Galilean invariance unless a gauge is fixed.  The
observer-relative S539 section-boundary functor is the invariant version.

## Proof Target

Prove a compactified section-checkoff lemma and a wall-switch lemma:

> If the compactified runner-section graph has a perfect matching but the
> strict-open graph has a Hall deficit, then some endpoint or boundary switch
> either opens a new section for the deficient packet or reduces to an
> endpoint-mouth / observer-source certificate.

For LRC14, this should interface with HYP-2569: the endpoint-mouth determinant
is exactly the local mechanism by which wall-only section debt becomes positive
open interval mass.

Cross-links: HYP-2024, S539, HYP-2486, HYP-2568, HYP-2569, OPEN-Q-108, THM-523,
T838.
