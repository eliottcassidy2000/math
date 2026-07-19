# The proposed gap theorem is false—and would have been stronger than INV

> **SECOND CORRECTION (THM-1153 / MISTAKE-170):** the comparison target
> `INVcov` used below is itself false on the doubled AP `2*{1,...,13}`.  The
> proposed `1/12` gap remains independently refuted by THM-1131.  Thus the old
> logical implication points to a refuted premise and supplies no proof route;
> a live inverse target needs primitive normalization plus a re-Covering bridge.

*boxeph-2026-07-18-S111; scope correction codex-2026-07-18 (MISTAKE-167). Owner: prove the gap theorem
"non-AP core ⟹ M ≥ 1/12". The script evaluated five named families; it did **not** enumerate covering
families and therefore did not verify a universal non-AP floor.  THM-1131 now gives the primitive fully
Covering non-AP row `W={2,3,5,8,9,11,12,13,14,15,17,20,23}` with
`M(W)=3/37` strictly between `1/13` and `1/12`, so the proposed gap theorem and empty window are false.
What survives is only the logical implication that this false strengthening would have implied the
corrected covering inverse target `INVcov`.  The one genuinely proved piece is the non-covering sieve
branch. LRC(14) is not closed.*

Throughout this historical reflection, “INV” should be read as the corrected target `INVcov` on the
fully covering stratum (MISTAKE-166).  Statements that the `1/12` floor or empty window is “verified” are
finite evidence from the five displayed rows only and are refuted by THM-1131.

## What the five-row probe shows

For **covering** 13-families (the INV regime), computing `M`:

| family | covering? | AP core? | `M` |
|---|---|---|---|
| deep well `{1..12}∪{182}` | yes | yes (AP) | `14/183 = 0.0765` |
| `{2,4,…,24}∪{182}` | yes | yes (dilated AP) | `7/92 = 0.0761` |
| `{2,4,…,22,26}∪{780}` | yes | **no** | `1/12 = 0.0833` |
| `{…,14}∪{182}` (repl 7→14) | yes | **no** | `1/10 = 0.100` |
| `{…,24}∪{182}` (repl 5→24) | yes | **no** | `1/9 = 0.111` |

In this five-family probe, AP cores dip below `1/13` and the smallest displayed non-AP value is exactly
`1/12`.  The probe simply missed THM-1131's fully Covering non-AP row at `3/37`; it supplies no live
conjecture about an empty window.

## The correction: the gap theorem is ≥ INV, not < INV

S110 proposed the gap theorem as a *more tractable* route to INV, on the "stability with slack is easier"
intuition. **That is wrong here**, and the logic is clean:

> `INV` (contrapositive): non-AP core ⟹ `M ≥ 1/13`.
> gap theorem:           non-AP core ⟹ `M ≥ 1/12`.
> Same hypothesis; `1/12 > 1/13`, so the gap theorem's **conclusion is strictly stronger**. Hence
> **gap theorem ⟹ INV**. Proving the gap theorem proves LRC(14)'s covering case.

The implication is sound, but its premise is false: THM-1131 gives
`1/13<3/37<1/12` on the exact fully Covering non-AP row `W`.  Thus the displayed `M=1/12`
row is not a sharpness witness, and there is no conjectural empty gap to exploit.  S110's route is
retracted rather than merely left open.

The direction that *does* have room is `INVcov` itself: to prove `non-AP ⟹ M ≥ 1/13` on the fully
covering stratum must work at the actual `1/13` boundary and allow non-AP values immediately above it.
The content remains irreducibly INV.

## What IS provable (the non-covering side)

The gap theorem does split, and one half is clean:

> **Provable (sieve):** if some modulus `q' ∈ {2,…,12}` divides no speed, then `t = 1/q'` gives
> `M ≥ 1/q' ≥ 1/12` (`lonely_of_no_multiple`, S106). So **non-covering-at-`{2..12}` ⟹ M ≥ 1/12**,
> for any core, AP or not.

This reduces the gap theorem to families **covering all of `{2,…,12}`**. Among those, the ones also
covering `13` are fully covering (the INV regime, open); the ones missing only `13` get `M ≥ 1/13` from
the sieve at `13` (not necessarily `≥ 1/12`). So the honest residue of the gap theorem is again the
covering / INV core — nothing new is bought.

## Net (honest)

- **Finite evidence only:** five named covering families were evaluated; among the three displayed
  non-AP examples the minimum is `1/12`.  The sample missed a fully Covering row at `3/37`.
- **Refuted (THM-1131):** the gap theorem is false, although it is logically stronger than `INVcov`.
  The supposed slack does not exist.
- **Provable fragment:** non-covering-at-`{2..12}` ⟹ `M ≥ 1/12` (sieve). The remainder is INV.
- **Right target:** `INVcov` at the `1/13` boundary.  Relating it to the global Tao `n=12` inverse still
  needs the separate structural bridge recorded by the corrected formalization.

The PG(2,13) / metagraph picture (S110) remains a parameter analogy, not a proved structural connection:
the residue quotient can identify two integer lifts with different global maxima.  Honest terminus:
LRC(14) rests on the correctly scoped covering inverse target, with lift data retained.

Cross-links:
[[lrc14-lives-at-PG-2-13-and-INV-is-the-transitive-class-isolation-a-doubling-gap-boxeph-S110]]
(route corrected),
[[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]],
THM-724 (covering-min), [[the-non-covering-sieve-dispatch-is-now-kernel-pure-lean-boxeph-S106]].
