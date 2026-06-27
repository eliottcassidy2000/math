# LRC14 Automatic Fiber Zipper - S187

S187 turns the automaton thread into a proof obligation.  The earlier
Fermat-Catalan, Moser-de Bruijn, fibbinary, Ostrowski-Hadamard, and Hurwitz
ideas are still useful, but the full HYP-2963 bank says they are warning
labels unless the magnitude coordinate is retained.

The computational signal is sharp:

```text
automatic_word           225 fibers, 143 mixed-route fibers
residue_terminal_fiber 16555 fibers, 265 mixed-route fibers
magnitude_cocycle      21909 fibers, 0 mixed-route fibers
barcode_shadow         21913 fibers, 0 mixed-route fibers
packet_zipper          21913 fibers, 0 mixed-route fibers
```

So the next proof lemma should not be "finite automata classify LRC14."  It
should be:

```text
inside a fixed automatic/residue fiber, the exact magnitude cocycle is forced
to one of the known route mechanisms.
```

The user's outside cues fit this form.  The 2-adic Littlewood/Hurwitz-doubling
paper is about controlling continued-fraction states under repeated doubling;
here that suggests a native-transition state, not a terminal-state quotient.
The Ostrowski-Hadamard gap theorem warns that lacunary support can carry a
hard boundary, so gap labels should stay attached to boundary/certificate
data.  Moser and fibbinary remain the right finite automata for base-4
even-position support and Zeckendorf/no-adjacent carries, but the route-purity
test says they are sidecars.  Fermat-Catalan power data is a no-lift ledger
for perfect-power payloads, not an LRC certificate by itself.

Source check: arXiv:2506.04110 is the 2025 Vitorino-Vukusic paper on bounds
for the 2-adic Littlewood conjecture and explicitly uses Hurwitz doubling of
continued fractions, so it supports a native doubling-transition sidecar:
https://arxiv.org/abs/2506.04110.  The Ostrowski-Hadamard page records the
lacunary ratio-gap natural-boundary theorem, so the LRC use is only a gap
boundary warning, not an analytic-continuation proof:
https://en.wikipedia.org/wiki/Ostrowski%E2%80%93Hadamard_gap_theorem.  The ACM
PDF endpoint was access-blocked during this session, but DOI metadata identifies
the linked paper as Hall-Holt, Katz, Kumar, Mitchell, and Sityon's "Finding
large sticks and potatoes in polygons"; I am using it only as the convex-subset
"potato" analogy for packet lenses, not as a theorem input:
https://doi.org/10.1145/1109557.1109610.

The concrete target is still `MFCMMCCFFFCCC`.  It contains AP/GW equality,
q-witnesses, one petal row, open-Haar loose rows, and covering rows.  Residue
terminal fibers reduce the mixing but keep AP in a 30-row mixed fiber with
covering and q-witness rows.  Exact magnitude splits it.

That gives the next manual proof route:

```text
1. Partition MFCMMCCFFFCCC by exact M.
2. For each of its 33 M-values, identify the family formula or certificate.
3. Show no family crosses route after magnitude is fixed.
4. Promote that local statement to every automatic/residue fiber in the bank.
```

This is the cleanest creative angle from the automaton stack so far: prove the
zipper reconstructs the labelled packet route, rather than trying to make a
single sequence shadow carry the theorem.
