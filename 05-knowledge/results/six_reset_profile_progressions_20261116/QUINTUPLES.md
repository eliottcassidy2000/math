# Two all-height, full-interface five-digit windows

In each line the five prescribed values are q*2^t+I+2j, j=0,1,2,3,4.
Each has a strong tournament with a normalized root, with common
(g,E,O)=(g,q*2^(t-1),q*2^(t-1)+g):

| q | I | g | heights |
|---:|---:|---:|---:|
|167|−75|15|t≥3|
|601|−315|33|t≥2|

The first family starts at 1261,1263,1265,1267,1269.
The check_quintuples.py data are deletion tournament outgoing masks, root's
outgoing mask, and the two seed invariants used in the search. It
makes each seed with the root as vertex 8 and recalculates every
marked statistic by subset DP; it does not implement the search.
Its output also gives raw full rows and raw profiles.

Write Z=E+O-g, I'=H-2O+2g. The inherited near-sink transfer
(H,E,O,g) -> (H+2E,E+O-g,E+O,g) doubles Z and preserves
I' and g, while after the first step O=E+g.
The literal seeds give Z=q*2^e and I'=I+2j, so append
t-e near-sinks. The respective exponents e are
(2,0,1,2,0) and (0,0,1,0,0), ensuring t-e≥1.
The literal verifier checks the first two allowed heights for each
family. In particular the first allowed height is strong; every
subsequent near-sink has incoming old nonroot arcs and an outgoing
arc to the root, so it preserves strongness.

The common full continuation distinguishes this from an old ray on
which only E matches. Any common fixed-root near-sink/ordinary
sink word preserves equality of (g,E,O) and shifts all five counts
equally, because the scalar updates add expressions independent of
H while the last three coordinates evolve independently of H.
This is a five-way selectable core on explicitly anchored rays and
their shared-word images; it is not a claim that these absolute
windows cover all targets, or that width five is maximal. A word
ending in an ordinary sink need not itself give strong witnesses.
