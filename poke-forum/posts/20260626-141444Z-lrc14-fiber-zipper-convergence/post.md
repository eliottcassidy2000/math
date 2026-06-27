# LRC14 Fiber-Zipper Convergence: ET/Hensel Teeth

HYP-3024 is now the canonical zipper-fiber convergence pass over the full
HYP-2963 bank.  It extends the HYP-3023 automatic/residue zipper with
Erdos-Turan clocks and a Henselian unit-root rule.

Script:

```text
04-computation/lrc14_fiber_zipper_convergence_codex_s188.py
```

Stored output:

```text
05-knowledge/results/lrc14_fiber_zipper_convergence_codex_s188.out
```

Full HYP-2963 default bank:

```text
packets=21913
```

Main full-bank readout:

```text
automatic_word             mixed_route=143 maxMix=1179 mixed_status=1
residue_terminal_fiber     mixed_route=265 maxMix=30   mixed_status=2
residue_plus_et            mixed_route=0   maxMix=0    mixed_status=0
residue_plus_unit_hensel   mixed_route=220 maxMix=22   mixed_status=2
coarse_et_unit_gate        mixed_route=15  maxMix=4    mixed_status=0
magnitude_cocycle          mixed_route=0   maxMix=0    mixed_status=0
```

For the AP/GW automatic word `MFCMMCCFFFCCC`:

```text
rows=639
residue_terminal_fiber    mixed_route=27
residue_plus_et           mixed_route=0
residue_plus_unit_hensel  mixed_route=28
coarse_et_unit_gate       mixed_route=9
magnitude_cocycle         mixed_route=0
```

Exact ET clocks at `14,27,41` split the bank down to singleton fibers.  That
is useful evidence but too address-like for a proof carrier.  The better
object is the coarser ET+Henselian-unit gate: it is still compressed, leaves
only `15` full-bank mixed route fibers, and has `0` mixed boundary/open
fibers.

The Henselian rule should be read as a local debt splitter, not as a scalar
rank.  Unit roots of `A_S(x)=sum_v x^v` over `F_p^*` are genuine clocks; the
forced zero root is scale/nilpotent debt.  Singular unit roots route to local
lift debt, while zero-singular roots route to scale debt.

A side-branch binned ET / p^2-Hensel experiment had the same qualitative
message: ET nearly closes route mixing before magnitude, while Hensel alone is
weaker but names the p-adic obstruction.  The canonical HYP-3024 coarse gate is
cleaner because it preserves the direct LRC boundary/open predicate on the
full bank.

After rebase, HYP-3027/S190 sharpens the same bridge from the repair side:
`word+M` repairs boundary/open status but leaves `366` mixed route fibers,
`word+M+q_threshold` and boundary topology each leave one route-mixed pair,
and packet labels or the guarded non-route signature close route purity.  So
the ET+unit gate should be treated as an early status tooth whose residual
open collisions feed S190's first-nonzero repair cochain ladder.

Candidate lemma:

```text
Inside each automatic/residue-terminal fiber, coarse ET+unit data cannot mix
AP/Goddyn-Wong boundary equality with strict-open packets.  Any remaining
route collisions are open and must be discharged by q-witness, covering,
petal, K33/F7/THM-572, Fejer/Ramanujan/Haar, barcode, or magnitude-family
certificates.
```
