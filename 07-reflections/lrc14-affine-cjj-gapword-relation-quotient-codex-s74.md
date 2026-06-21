---
date: 2026-06-21
source: codex-2026-06-21-S74
tags: [lrc14, CJJ, Delsarte, affine, full-residue, relation-code, tournament-analysis]
---

# LRC14 Affine CJJ Gap-Word Quotient

The CJJ papers point to the right meta-lesson: completeness comes from choosing
the symmetry and moment lattice in which the optimizer is an integral object.
For binary linear codes this is the subspace lattice.  For LRC14 consec/AP it is
not the subspace lattice; AP is affine-linear, not linear.

Incoming HYP-2754 makes the natural repair: factor translation and dilation.
The exact scout here adds a guardrail.  If that affine quotient only remembers
residue subsets in `F_7`, it forgets everything on the HYP-2749 full-residue
stratum.  In the bounded k=8 stratum, all 528 shapes have the same residue
count, affine pair profile, affine triple profile, and affine quad profile.
Only the integer gap word separates the AP orbit from the rest.

The later S19/S20 pulls sharpen the symmetry half of this story.  The
multiplicative action `W_a(cE)=W_{ca}(E)` makes `Z_7^*` a genuine symmetry of
the AP/consec layer, but HYP-2762 corrects the overreach: Paley/QR is only the
small-p/spectral comparison, while the large-p tournament H-driver is again the
additive interval/AP.  The Freiman anchored-profile output also explains why
the quotient cannot stop at the dilation spine: leg profile is not sufficient,
while the anchored second profile/gap word carries the tie-break.

The fresh Huffer-Shepp port says the same thing in cell language.  Reflection
and symmetrisation are real, but consec loses individual `W_a` cells and wins
only in the aggregate sum.  Any proof that forgets the aggregate relation
ledger is proving the wrong statement.

The KPS c14-lift reframing points to the same wall from a different direction,
after its correction.  The known method fails at `14=2*7` because the canonical
tuple can no longer be certified by a Fermat-field argument, but the paper's
lift is additive rather than a literal CRT product.  HYP-2770 is the hierarchy
version of the corrected reading: prove the canonical/consec certificate
analytically on the seven-sector quotient, then keep the side-channel the
quotient forgets.

The latest decorrelated-wide closed form is another check on the same shape:
far decorrelation reads as a linear functional of the base missed-depth law,
and KPS now proves its maximum is the single-far/consec-base plateau
`Q(k-1)<cap_k`.  So the wide main term is closed in the moment-dual language
HYP-2770 tries to preserve.  The newest KPS atlas makes the signed resonance
error tiny, about `0.01211` in the k=9 commensurable bank, so even a loose
joint-discrepancy bound should be enough for the wide side.

So the next proof should be phrased as:

```text
full-residue reduction
  -> affine-normalized integer gap/generated-word quotient
  -> relation-code level-2 marginal extremality
  -> signed THM-534 Delsarte cap.
```

This also clarifies the tournament bridge.  Paley extremality can use a genuine
linear-code object; AP extremality needs an affine gap-word object.  The common
theta-prime ceiling is a comparison principle, not a replacement for the
retained channel.

The challenged assumption for this session: "affine vertices are residue
subsets."  The computation refutes that as a final quotient.  The useful
vertices are proof obligations or retained channels: gap word, pair marginals,
full-residue localization, signed Delsarte dual, and theta ceiling.
