# LRC14 design/Hodge carriers

S105 tested the user's unital/Clebsch/truncated-octahedral prompts against the
current HYP-2889/HYP-2890 residual-leak state.

The productive interpretation is not that any one graph is the LRC proof.  The
productive interpretation is a carrier split:

```text
q=3 unital and Clebsch = pair-uniform tight frames;
truncated octahedral and octahedral L(K4) = signed curl complexes.
```

Post-rebase namespace and proof guardrail: HYP-2891/T1005 is the primary
Clebsch/Bruhat labelled-residual quotient from the incoming S105 commit.  This
note is the HYP-2892/T1006 refinement focused on the `q=3` unital pair frame
and Hodge split.  KPS S31m also refuted the sparse-design exact-tiler reading:
the exact threshold tiler is the full residue system `Z_14 \ {0}` and its
dilates.  Thus the unital is not a tiling model; it is a candidate averaging
frame on the k=8 pair-slot residual after HYP-2890's same-frequency packet.

The `q=3` unital is the strongest surprise.  The Hermitian construction over
`GF(9)` gives `28` points and `63` four-point blocks, with every pair of points
in exactly one block.  Since `28=C(8,2)`, this is exactly sized for the k=8 AP
pair slots.  Its incidence identity `N^T N=8I+J` is a ready-made frame for
averaging additive quadruple residuals.

The Clebsch graph contributes a parallel but different frame.  As the folded
5-cube it has SRG parameters `(16,5,0,2)`, and closed neighborhoods form a
`2-(16,6,2)` design with `N^T N=4I+2J`.  That looks better suited to folded
parity/sign packets than to literal AP pair slots.

The truncated octahedral graph is even more directly connected to the failed
compression attempts.  It is the Cayley graph of `S4` on adjacent
transpositions: `24` vertices, `36` edges, `6` square commutation faces and `8`
hex braid faces.  A one-step compression is an edge.  HYP-2889 already showed
edge monotonicity is false.  The Coxeter complex says the next invariant should
be curl around square/hex faces.

This dovetails with HYP-2887: the support-six repeated-packet graph is
`L(K4)`, the smaller octahedral sphere with cycle rank `7`.  So the project now
has two Hodge carriers: truncated octahedral for local AP/Freiman
permutations, and octahedral `L(K4)` for repeated residue packets.

No LRC14 proof is claimed.  The new concrete next step is to put the S104
residual-leak functional on Bruhat edges/faces for low-depth near-AP rows.  If
positive leaks are face boundaries and the remaining block averages are
nonpositive in the unital/Clebsch frames, the HYP-2890 inequality becomes a
finite Hodge certificate plus HYP-2636 tail cancellation.

External prompts checked:

- https://mathworld.wolfram.com/Unital.html
- https://mathworld.wolfram.com/ClebschGraph.html
- https://mathworld.wolfram.com/TruncatedOctahedralGraph.html
