# LRC14 Lonely-Profile Persistence Barcode

codex-2026-06-26-S179

## What Changed

The creative turn is to treat the lonely profile

```text
m_S(t)=min_v ||v t||
```

as a one-dimensional persistence object.  At threshold `1/14`, every positive
safe component is a bar with an exact length, peak time, peak height, and
height margin.  This keeps the proof predicate, but it avoids collapsing all
open witnesses to a single scalar `M` or a single safe mass.

The exact scout
`04-computation/lrc14_lonely_profile_persistence_barcode_codex_s179.py`
confirms the expected boundary split:

```text
AP13 and GW 12->24: zero bars, M=1/14
K33 12->36:         two bars, M=3/41
petals:             two bars, M=2/27
covering 12->84:    eight bars, M=7/89
fibbinary first13:  38 bars, M=3/25
Moser first13:      64 bars, M=1/6
```

So automatic/lacunary rows are not merely positive; they have many bars.  The
hard near-boundary rows are thin-bar rows, especially K33 and some covering
components.  That is useful routing information for Fejer, endpoint-owner, and
state-lift certificates.

## External Prompt Fit

The ACM link in the prompt is the "large sticks and potatoes in polygons"
paper, so the useful import is not Fermat-Catalan arithmetic but an extremal
subset viewpoint: find a large certified subset inside a constrained geometry.
The barcode is the exact LRC version of that move, with each safe component a
candidate "stick" to certify.

The arXiv 2-adic Littlewood prompt and the Ostrowski-Hadamard gap prompt remain
guardrails from HYP-3009..HYP-3012: automatic and lacunary shadows can suggest
where positive mass lives, but the proof carrier has to remember exact endpoint
geometry.  The barcode does that by keeping the bar boundary, peak, and
persistence data attached to the fibbinary and Moser rows.

## Why It Is Useful

Raw `M` answers "how high does the best witness get?"  The barcode also answers:

```text
how many strict witness components exist
how wide the strongest component is
where a rational certificate should be centered
how far the row survives a threshold perturbation
whether a quotient erased safe topology
```

This makes the barcode a bridge between exact Haar/Baire safe components and
dual certificates.  A Fejer manifest can choose the peak or midpoint of a
specific bar; a boundary-moment route can focus on the lowest-persistence bars;
an automaton or divisor side channel can be tested for mixed barcode fibers.

## Assumption Challenge

Do not take runners, arcs, or sequence names as the default vertices here.

Alternate vertices considered:

```text
runners, gaps, fixed circle sections, section boundaries, wall-crossing events,
residues, cover arcs, Fourier modes, matroid circuits, persistence bars,
proof obligations.
```

The chosen vertices for Tournament Analysis are proof carriers.  The quotient
preserves the LRC predicate "there is an open bar above threshold" and destroys
runner identity only after exact height, topology, and peak geometry have been
stored.

## Next Pull

Run the barcode sidecar over the HYP-2963 bank.  The key finite target is not
"largest M"; it is the list of positive rows with the smallest persistence
margin or shortest longest-bar.  Those are the rows where the proof stack most
needs endpoint-owner and interval-Fejer precision.
