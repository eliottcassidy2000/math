# LRC14 Euler-Copy Coimage Tail

**Source:** codex-2026-06-19-S22, HYP-2630 / T878.

The useful outcome is not "Euler-copy mass explains everything." It explains
the packet capacity, and then it stops. Exact-period packets over denominators
divisible by `7` distribute uniformly across `F_7^*`; for raw `1260`, the top
period gives `48` copies per unit residue and the full `{2,3,5,7}` squarefree
mask gives `96`.

That uniformity is the point. It means the QR/NQR split in the HYP-2626
repeated tail cannot be a raw copy-mass effect. The `4+2` classes
`(1,1,1,1,a,a)` have identical full-mask copy capacities for nonzero `a`, but
the coimage mass separates by quadratic character:

```text
chi(a)=+1: |S_9|=0.23891209
chi(a)=-1: |S_9|=0.17201670
```

So the proof coordinate is:

```text
copy mass = residue capacity
quadratic character = phase
multi-large wall relation = incidence / compatibility
```

This is exactly the older residue/phase/incidence lesson in LRC14 clothing.
The next theorem should not be another one-large wall height extension. Height
`3` gives the same `85/116` nonzero classes as height `2`. The tail is
structurally multi-large: the bounded core `1..13` has at most two
representatives of each nonzero residue modulo `7`, so four equal residues need
at least two large speeds.

The sharp target is now a two-large repeated-residue cotangent/Dedekind bound
that keeps the quadratic-character phase channel visible.

The concurrent HYP-2631 AP-drop repair result points in the same direction
from a different face.  The radical `210` projection deletes reduced-denominator
packets that raw `1260` still sees.  Here the squarefree copy mass survives but
is too uniform, so the missing coordinate is not more mass; it is the retained
packet address plus the `chi_7` phase before the final coimage quotient.
