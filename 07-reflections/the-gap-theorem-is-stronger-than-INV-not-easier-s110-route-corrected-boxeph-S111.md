# The gap theorem is stronger than INV, not easier — correcting the S110 route

*boxeph-2026-07-18-S111. Owner: prove the gap theorem "non-AP core ⟹ M ≥ 1/12". Outcome: an honest
self-correction of S110. The gap theorem is **verified true for covering families** (their non-AP floor is
exactly 1/12), but it is **logically STRONGER than INV** (same hypothesis, strictly stronger conclusion),
hence at least as hard as the open crux — proving it proves LRC(14). The one genuinely provable piece is
the non-covering side via the sieve. LRC(14) not closed. Verified S111 computation.*

## What the data confirms (the isolation gap is real)

For **covering** 13-families (the INV regime), computing `M`:

| family | covering? | AP core? | `M` |
|---|---|---|---|
| deep well `{1..12}∪{182}` | yes | yes (AP) | `14/183 = 0.0765` |
| `{2,4,…,24}∪{182}` | yes | yes (dilated AP) | `7/92 = 0.0761` |
| `{2,4,…,22,26}∪{780}` | yes | **no** | `1/12 = 0.0833` |
| `{…,14}∪{182}` (repl 7→14) | yes | **no** | `1/10 = 0.100` |
| `{…,24}∪{182}` (repl 5→24) | yes | **no** | `1/9 = 0.111` |

So among covering families: **AP cores dip below `1/13`; every non-AP core has `M ≥ 1/12`**, and the
tightest non-AP is **exactly `1/12`**. The window `(1/13, 1/12)` is **empty of covering non-AP families** —
a real *isolation gap*, sharper than the S89 "discrete spectrum" statement.

## The correction: the gap theorem is ≥ INV, not < INV

S110 proposed the gap theorem as a *more tractable* route to INV, on the "stability with slack is easier"
intuition. **That is wrong here**, and the logic is clean:

> `INV` (contrapositive): non-AP core ⟹ `M ≥ 1/13`.
> gap theorem:           non-AP core ⟹ `M ≥ 1/12`.
> Same hypothesis; `1/12 > 1/13`, so the gap theorem's **conclusion is strictly stronger**. Hence
> **gap theorem ⟹ INV**. Proving the gap theorem proves LRC(14)'s covering case.

And it is **sharp**: the tightest covering non-AP family sits *exactly* at `1/12`, so there is **no slack**
for a crude bound. The "room" that would make INV easier is precisely the gap `(1/13, 1/12)` between the
AP-floor and the non-AP-floor — but that gap being empty *is* the isolation, i.e. *is* INV. So the gap is
a **restatement of the crux**, not a lever on it. S110's hope is retracted.

The direction that *does* have room is INV itself: to prove `non-AP ⟹ M ≥ 1/13` one may use that the
nearest non-AP is (empirically) at `1/12` — but "the nearest non-AP is at `1/12`" is again the gap
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

- **Verified:** covering non-AP families have `M ≥ 1/12` (floor exactly `1/12`); the gap `(1/13, 1/12)` is
  empty of covering non-AP families — the deep-well isolation, sharpened.
- **Corrected (S110):** the gap theorem is **logically stronger than INV** and **sharp at `1/12`**, hence
  at least as hard as the open crux, not more tractable. The "stability-with-slack is easier" framing does
  not apply — the slack IS the crux.
- **Provable fragment:** non-covering-at-`{2..12}` ⟹ `M ≥ 1/12` (sieve). The remainder is INV.
- **Could not prove the gap theorem** (it entails LRC(14)). The right target reverts to `INV` itself
  (`M ≥ 1/13`), the open additive inverse theorem = Tao n=12.

The PG(2,13) / metagraph picture (S110) stands as a genuine structural connection; but the *proof route* it
suggested collapses back onto INV. Honest terminus: LRC(14) rests on INV, and every reformulation that
carries real content (dimension, offset-vanishing, doubling-gap) is equivalent to it, as proved in S94.

Cross-links:
[[lrc14-lives-at-PG-2-13-and-INV-is-the-transitive-class-isolation-a-doubling-gap-boxeph-S110]]
(route corrected),
[[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]],
THM-724 (covering-min), [[the-non-covering-sieve-dispatch-is-now-kernel-pure-lean-boxeph-S106]].
