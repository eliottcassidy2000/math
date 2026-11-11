# Reverse target search exposes a neutral-padding boundary

**FINITE-EXACT construction probe.** This complements the strong-source
two-birth certificate in [ALTERNATE141](../history/ALTERNATE141.md). It does
not produce a new strong H141 source with a one-birth425 extension.

Take all distinct literal graphs of count425 returned by Session251's
normalized inverse for the five cores, all90 raw cuts and both flank
addresses. There are42 such target graphs. Inspecting every one-vertex
deletion gives six labelled hits of count141, all with11-vertex targets
and10-vertex parents. Each target uses C3 with one flank address72 and
the other zero; the complete cuts and label maps are in
[reverse425.json](reverse425.json).

Every resulting H141 parent is WEAK: it has exactly one global source or
sink. Removing this neutral singleton leaves a strong9-vertex H141 source.
Thus the actual construction has count trace

    strong141 -> weak141 -> strong425.

This is a one-birth extension of a weak H141 tournament and a two-birth
extension of its strong ancestor. It does not contradict the independent
search finding no one-birth425 extension among its nineteen strong bank
classes, nor the all-order obstruction for the two specified C4 sources.
The underlying source is different.

For example the target rows are

    (1018,2044,2041,1008,992,960,1920,768,512,1024,441).

Delete vertex1. The H141 parent is

    (508,1021,504,496,480,960,384,256,512,221).

Its vertex1 is a global source. Adding a new vertex with outmask1022
recovers the target up to the explicit label permutation in the receipt.
The receipt separately records the strong ancestor and every old/new label.

The inverse normalization preserves count support, not all graph embeddings.
Accordingly this is an exhaustive deletion check only on the42 returned
literal targets. It is not an all-order classification of H141->425 edges.
`python -I -B reverse425.py` reproduces the exact receipt with literal subset
Hamilton counts and strongness checks.
