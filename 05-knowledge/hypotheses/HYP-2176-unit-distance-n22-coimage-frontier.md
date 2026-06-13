---
id: HYP-2176
status: SUPPORTED synthesis from S614 scout and current literature; exact value of u(22) open
source: user-2026-06-03; codex-2026-06-03-S614
tags: [unit-distance, n22, coimage, tournament-analysis, LRC, Moser-lattice, grid-conjecture]
---

# HYP-2176: unit-distance n=22 is a coimage-with-side-channels frontier

The planar unit-distance problem at `n=22` is the first small case where the
published graph-only upper coimage is visibly too coarse while the exact
geometric optimum remains one bit open:

```text
60 <= u(22) <= 61.
```

The claim is not that S614 determines `u(22)`.  The claim is that the productive
proof object for `n=22` is a retained-side-channel object, analogous to the
recent LRC `n=14` quotient work: graph coimages, Moser coordinates, dense
deletion cores, and unfaithful-subgraph certificates must be kept together.

## Evidence

Alexeev-Tikhonov prove exact values through `n=21`, ending with `u(21)=57`,
and show `u(22)<=61`: the `62`-edge `F`-free frontier has exactly two
graph-only candidates, and both contain totally unfaithful subgraphs.  The
known lower bound is `60`, realized by Moser-ring database examples from
Engel-Hammond-Lee-Su-Varga-Zsamboki.

Rebase integration note: incoming HYP-2170/S599 records `49` unit distances
inside the triangular/Harborth lattice search family.  That is a useful
structured baseline and a warning about lattice-specific coimages, but it is
not the planar `n=22` optimum frontier once the published Moser-ring examples
with `60` edges are included.

S614 decodes the five exact `n=21`, `57`-edge graph6 cores from
Alexeev-Tikhonov and records their degree and cycle profiles.  This gives the
immediate deletion lemma:

1. If a `61`-edge UDG on `22` vertices exists, its average degree is `122/22`,
   so it has a vertex of degree at most `5`.
2. Deleting such a vertex leaves at least `56` edges on `21` vertices.
3. Since `u(21)=57`, a `61`-edge UDG cannot have a vertex of degree `0..3`.
4. Therefore every `61`-edge solution has minimum degree `4` or `5` and lies
   over a `56`- or `57`-edge `21`-vertex UDG core with a degree-5 or degree-4
   embeddable extension vertex.

For a `60`-edge witness, the same argument only forces minimum degree at least
`3`, and dense `21`-decks in the `55..57` range.

## Tournament reading

Tournament vertices should be proof routes or proof obligations, not merely
points.  S614 orients the route tournament by a pairwise observable combining
proof burden, side information retained, relation rank, `n=22` specificity,
coimage loss, and geometric contact.  The high-score routes are:

- Moser-ring beam/search toward `61`;
- extending exact or near-exact `21`-vertex cores;
- mining totally unfaithful obstructions.

The graph-only `62` transfer and plain triangular/square grid variants score
low because they forget the side information that already proved decisive.

## Relation to LRC

The analogy with recent LRC `n=14` work is literal at the proof-object level.
In LRC, the `Res_27` quotient is useful only when owner labels, carry cocycles,
and endpoint/pinch side channels are reattached.  In unit distances, the
`F`-free graph quotient is useful only when geometric embeddability, totally
unfaithful subgraphs, and Moser-coordinate data are reattached.

The recent asymptotic disproof of the Erdos unit-distance conjecture reinforces
the same moral.  The old grid intuition fixes a low-degree visible lattice
(`Q(i)` for the classical construction).  The counterexample fixes small
splitting primes and varies the CM field degree through a Golod-Shafarevich
class-field tower, then projects one coordinate to the plane.  The visible
planar picture is a coimage of a higher-dimensional arithmetic carrier.

## Next targets

1. Obtain or sample the Engel et al. Moser-ring `22`-vertex, `60`-edge
   examples and compute their `21`-decks, degree profiles, and unfaithful
   subgraph residues.
2. Enumerate or collect `21`-vertex `56`-edge near-max UDG cores; a `61`
   witness must extend one of the `56/57` cores by degree `5/4`.
3. Build a geometry-aware extension solver: given a dense embedded `21`-core,
   enumerate possible points at unit distance from `4` or `5` core vertices
   and test whether any new edges force `61`.
4. Turn the totally-unfaithful filter into a reusable obstruction library for
   candidate `61` graphs.

## Artifacts

- `04-computation/unit_distance_n22_tournament_lrc_s614.py`
- `05-knowledge/results/unit_distance_n22_tournament_lrc_s614.out`
- `07-reflections/unit-distance-n22-tournament-lrc-grid-disproof-s614.md`
