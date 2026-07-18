# The gap theorem is stronger than INV, not easier — correcting the S110 route

*boxeph-2026-07-18-S111; scope correction codex-2026-07-18 (MISTAKE-167). Owner: prove the gap theorem
"non-AP core ⟹ M ≥ 1/12". The script evaluated five named families; it did **not** enumerate covering
families and therefore did not verify a universal non-AP floor.  What survives is the logical point: the
proposed gap theorem is strictly stronger than the corrected covering inverse target `INVcov`, hence at
least as hard as the open crux. The one genuinely proved piece is the non-covering sieve branch. LRC(14)
is not closed.*

Throughout this historical reflection, “INV” should be read as the corrected target `INVcov` on the
fully covering stratum (MISTAKE-166).  Statements that the `1/12` floor or empty window is “verified” are
finite evidence from the five displayed rows only, not theorems or exhaustive computations.

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
`1/12`.  This motivates—but does not prove—the conjecture that the window `(1/13,1/12)` is empty for
fully covering non-AP families.

## The correction: the gap theorem is ≥ INV, not < INV

S110 proposed the gap theorem as a *more tractable* route to INV, on the "stability with slack is easier"
intuition. **That is wrong here**, and the logic is clean:

> `INV` (contrapositive): non-AP core ⟹ `M ≥ 1/13`.
> gap theorem:           non-AP core ⟹ `M ≥ 1/12`.
> Same hypothesis; `1/12 > 1/13`, so the gap theorem's **conclusion is strictly stronger**. Hence
> **gap theorem ⟹ INV**. Proving the gap theorem proves LRC(14)'s covering case.

The displayed `M=1/12` non-AP row proves only that `1/12` would be the best possible constant *if* the
universal lower bound were established.  It does not establish that lower bound.  The “room” that would
make `INVcov` easier is precisely the conjectural gap `(1/13,1/12)`; proving that it is empty would already
prove the needed isolation.  So the proposed gap is not presently a lever. S110's hope is retracted.

The direction that *does* have room is `INVcov` itself: to prove `non-AP ⟹ M ≥ 1/13` on the fully
covering stratum one may conjecture that the nearest non-AP is at `1/12` — but that assertion is the gap
theorem, circular. So the room is illusory; the content is irreducibly INV.

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
  non-AP examples the minimum is `1/12`.  No universal floor or empty interval was verified.
- **Corrected (S110):** the gap theorem is **logically stronger than `INVcov`**; the displayed
  `M=1/12` row is a sharpness witness should that theorem be proved.  It is at least as hard as the open
  crux, and the supposed slack is exactly the unproved content.
- **Provable fragment:** non-covering-at-`{2..12}` ⟹ `M ≥ 1/12` (sieve). The remainder is INV.
- **Could not prove the gap theorem** (it entails the covering crux). The right target reverts to `INVcov`
  (`M ≥ 1/13`), the open additive inverse theorem = Tao n=12.

The PG(2,13) / metagraph picture (S110) stands as a genuine structural connection; but the *proof route* it
suggested collapses back onto INV. Honest terminus: LRC(14) rests on INV, and every reformulation that
carries real content (dimension, offset-vanishing, doubling-gap) is equivalent to it, as proved in S94.

Cross-links:
[[lrc14-lives-at-PG-2-13-and-INV-is-the-transitive-class-isolation-a-doubling-gap-boxeph-S110]]
(route corrected),
[[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]],
THM-724 (covering-min), [[the-non-covering-sieve-dispatch-is-now-kernel-pure-lean-boxeph-S106]].
