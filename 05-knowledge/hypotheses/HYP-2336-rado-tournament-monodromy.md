---
id: HYP-2336
title: The Rado graph as a tournament = the random/universal homogeneous tournament; the non-residue "loop" swaps the two outputs (monodromy); Paley → Rado via Weil/√p
status: OPEN (synthesis); the non-residue swap + Paley k-EC + Weil rate VERIFIED
source: claudebox-2026-06-06-S641
related:
  - HYP-2321  # Paley = QR-trienement; quadratic reciprocity / the √-monodromy
  - HYP-2311  # character-ratio spectrum (= the Weil/Gauss bound = √p)
  - THM-258   # Paley subtournament uniformity (quasi-randomness)
---

# HYP-2336 — the Rado-tournament, and the loop that swaps the outputs

"Consider the Rado graph as a tournament" = the **unique countable homogeneous (random) tournament**
— the Fraïssé/`n→∞` limit of all finite tournaments, the infinite object our finite program (H,
delta field, glass) approximates. "A loop of the input swaps the two outputs" = **monodromy**: a loop
around the `√` branch point swaps two sheets (`Z/2`).

## Verified

1. **The loop = a quadratic non-residue.** In the Paley tournament `T_p`, `x ↦ g·x` for a non-residue
   `g` sends `T_p → T_pᵒᵖ` (reverses every arc = self-complementation), because `χ(g·d)=−χ(d)`. And
   `g² ∈ QR`, so `g` generates the `Z/2` of the QR double cover `(ℤ/p)* → (ℤ/p)*/(squares)`. So **the
   non-residue "loop" swaps the two outputs (the two orientations / the two `√` branches)** — the
   monodromy is literally self-complementation. (Verified p=7,11,19,23.)

2. **Paley → Rado via Weil/√p.** Paley tournaments satisfy the extension property (k-existential
   closure): every disjoint `U,V`, `|U|+|V|≤k`, has a common "dominator" `x` (`x→U`, `V→x`). Count
   `≈ p/2^{|U|+|V|} ± (|U|+|V|)√p` (Weil). Verified k-EC: 2-EC at p=7,11; 3-EC at p≥19. As `p→∞`,
   `T_p →` the random tournament (k-EC for **all** k). **The `√p` error is the Gauss sum = the
   character-ratio spectrum (S638) = the monodromy `Z/2`** — the quadratic "loop" is exactly the rate
   of approach to the universal limit.

## The synthesis
The finite tournament program (H, the delta field, the glass transition, the character-ratio
spectrum) has its `n→∞` limit in the **Rado/random tournament**; the **Paley tournaments are the
explicit quasi-random approximants**; and the **quadratic monodromy (the non-residue loop swapping
`T↔Tᵒᵖ`, i.e. the `√` `Z/2`) is the engine** — both the self-complement symmetry of each finite Paley
tournament and, via Weil/√p, the convergence to the universal limit. The "loop swaps two outputs" is
the same `±1`/`√`/out-of-phase structure as quadratic reciprocity (HYP-2321) and the FTA branch
(HYP-2326), now seen as the deck transformation taking the finite program to its Rado limit.
