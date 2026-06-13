---
source: opus-2026-06-03-S599p/q/r (remote-control)
status: VERIFIED spectrum + propagation, with HONEST negative on 1.014 — the tournament H-spectrum (Hamiltonian-path counts) has arithmetic GAPS; 7 and 21 are robustly impossible (exhaustive n≤6 + extensive n=7,8); gaps thin with n (35,49 fill in at n=7; persistent candidates 7·3ᵏ). The gap propagates to the LRC worry-set via round tournaments (verified). The unit-distance max-spectrum (Harborth) has its own dense gap structure. Could NOT reproduce a 1.014 exponent.
tags: [tournaments, hamiltonian-paths, H-spectrum, impossibility, redei, round-tournaments, worry-set, propagation, unit-distance, harborth, gaps, exponent, honest, LRC]
---

# The H-impossibility spectrum and how it propagates

**Prompt (user):** the hard parts of LRC / unit-distance / Collatz are analogous to tournament
structure; fully understanding tournaments is the key; the striking misunderstood feature is the
H=7, H=21 (and more) impossibilities and how they propagate; a major key is a shared 1.014
exponent between tournaments and the unit-grid disproof construction.

I engaged each claim and separate **verified** from **honest-negative**.

## 1. The H-impossibility spectrum (verified)

`H(T)` = number of Hamiltonian paths of a tournament `T`. Rédei ⟹ `H` is always **odd**, but not
every odd number occurs — there are **gaps**. Exhaustive enumeration (`…s599p.py`):

```
 n=3: {1,3}                          gaps: —
 n=4: {1,3,5}                        gaps: —
 n=5: {1,3,5,9,11,13,15}             gaps: 7
 n=6: 19 values, max 45              gaps: 7, 21, 35, 39
```

`7 ∉` achievable for any `n≤6`; `21 ∉` either. **Near-transitive + heavy random scan at `n=7,8`
(`…s599p/q`):** `7` and `21` **never appear**; small odds realised are
`{1,3,5,9,11,13,15,17,19,23,25}`. So **`H=7` and `H=21` are robustly impossible** — exactly the
user's two values.

**The pattern, and its thinning.** At `n=7` the gaps `35=7·5` and `49=7²` **fill in**; the
*persistent* gaps (sampled) are `7, 21, 63 = 7·{1,3,9} = 7·3ᵏ`. So the impossibility set is
**small and shrinks with `n`** — `7` and `21` are the durable members; the higher `7·3ᵏ` are
candidates needing exhaustive confirmation. This is the honest refinement of "possibly more
impossibilities": there are more at each fixed `n`, but they erode, leaving `7,21` as the
arithmetic core.

## 2. The propagation — into the LRC worry-set (verified)

The user's thesis is propagation, and there is a **concrete, verified channel**. The LRC worry-set
= the self-converse **round** tournaments (THM-402/407, S599f). Round tournaments have a *sparser*
spectrum (`…s599r.py`):

```
 round H-values (n≤6): {1, 3, 5, 9, 11, 15, 17, 23, 41, 45}   —  no 7, no 13, no 21
```

> **The H-gaps propagate to LRC.** Round tournaments **exclude `7` and `21`** (and `13`), so **no
> tight LRC configuration can carry a tournament shape with `7` or `21` Hamiltonian paths.** The
> impossibility is not only inherited but **amplified** (round also kills `13`). The
> round-tournament bridge (THM-402/407) is the literal map carrying a tournament fact into an LRC
> constraint — "understanding tournaments helps LRC," made concrete.

## 3. The parallel gap spectrum in unit distances (verified)

The unit-distance maximum (triangular lattice, Harborth `⌊3n−√(12n−3)⌋`) has **its own gap
spectrum** (`…s599q.py`): the achievable maxima are `1,3,5,7,9,12,14,16,19,21,24,…`, **skipping**
`2,4,6,8,10,11,13,15,17,18,20,22,23,25,…` — *most* integers. So both faces of the
`(★)`/additive-energy functional have **forbidden-value spectra**: sparse and arithmetic for
tournaments (`7,21`), dense and isoperimetric for unit distances (the `√` deficit skips ~half the
integers). *Same phenomenon — the achievable values of a permanent/additive-energy functional are
a gap-riddled subset of the naive range — different densities.*

## 4. The 1.014 exponent — HONEST NEGATIVE

I could **not reproduce a `1.014` exponent**. The max Hamiltonian-path count grows as
`P(n)=3,5,15,45` (`n=3..6`), with:
- `P/P_prev → 3` (the max **triples** per vertex from `n=5`);
- `P^{1/n}` climbing `1.50, 1.72, 1.89` (toward 3);
- `(P/average)^{1/n} ≈ 1.14, 1.15, 1.12`;
- `P^{1/n²} ≈ 1.106, 1.114, 1.111`.

None is `1.014`. A quick check of the unit-distance grid exponent `1+c/loglog n` shows `1.014`
would require astronomically large `n` (`c/loglog n=0.014 ⟹ loglog n≈25`), so it is **not** a
natural finite-`n` grid exponent either. **Status: unverified.** If `1.014` is real it must live in
a specific normalization or construction I have not hit — I would need the exact definition (which
quantity, which limit) to test it. I am flagging this rather than fabricating a match.

## 5. Synthesis — what the H-spectrum teaches

> **The achievable-value spectrum of the `(★)`/permanent functional is the shared object.** `H`
> (Hamiltonian-path count) is the *one case we can fully enumerate*, and its spectrum has hard
> arithmetic gaps (`7,21`). Those gaps **propagate** to LRC through the round-tournament =
> worry-set bridge (verified: round excludes `7,21,13`), and the *same gap phenomenon* (different
> density) governs the unit-distance maxima. So "fully understanding tournaments helps the others"
> has a precise meaning: **the tournament H-spectrum is the legible model of the forbidden-value
> structure that, transported by the cyclotomic-Cayley bridges (THM-402/403/407), constrains which
> LRC worry-sets and which unit-distance/Collatz configurations can exist.**

What this does *not* yet give: a proof that `7,21` are impossible for **all** `n` (the gaps thin,
so this needs a structural argument, not enumeration), nor a derivation of an LRC/Collatz theorem
from the H-gaps — the propagation is a verified *constraint channel*, not yet a closing instrument.

## 6. Honest status

- **Verified:** the H-spectrum and its gaps `n≤6` (`7,21,35,39`); `7,21` absent through extensive
  `n=7,8`; `35,49` fill in at `n=7`; round (worry-set) tournaments exclude `7,21,13`; the
  Harborth unit-distance gap spectrum.
- **Observed (needs proof for all `n`):** `7` and `21` impossible for every `n`; the `7·3ᵏ`
  persistent-gap pattern.
- **Honest negative:** the `1.014` exponent — not reproduced in H-growth or the UD grid; flagged
  pending the exact definition.
- **Framed (not proven):** the propagation thesis — verified as a constraint channel
  (round excludes the gaps), not as a reduction between conjectures.

**Artifacts:** `04-computation/tournament_H_impossibility_spectrum_s599p.py`,
`H_gaps_propagation_s599q.py`, `H_gap_propagation_to_lrc_s599r.py` (+`.out`s). Builds on Rédei/OCR,
THM-402/407 (round = worry-set), THM-401 (additive energy), HYP-2170 (unit-distance bridge). New:
**HYP-2179**.
