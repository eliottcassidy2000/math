---
source: mac-mini-2026-07-21-S160 (cont.)
status: THM-1966 — the signed Rédei count |R| is spectral (n≤5), (spectrum,H)-measurable (n=6), and genuinely independent of (spectrum,H) from n=7; answers the S160 highest-leverage handoff
tags:
  - redei
  - signed-count
  - invariant-ladder
  - co-spectral
  - hamiltonian-paths
  - THM-1966
---

# The signed Rédei count is a genuinely new invariant — from n=7

The S160 handoff asked the highest-leverage question about `R(T)=Σ_{Ham paths}sgn(π)`: **is `|R|`
a genuinely new invariant, or is it derivable from the Rédei count `H` and the spectrum?** The
answer places `|R|` precisely on THM-1780's moment ladder — and it is a *staged* answer.

## The three rungs

- **`n ≤ 5`: `|R|` is spectral.** Constant on every co-spectral class — a trace moment, like `H`.
- **`n = 6`: `|R|` leaves the ladder *with* `H`, but adds nothing beyond it.** It splits exactly
  the same 3 co-spectral classes `H` splits (perfect coupling — 0 classes split one but not the
  other), and `(spectrum,H,|R|)` = `(spectrum,H)` distinguish the *same* 32 of 56 iso classes.
  So at `n=6`, `|R| = f(spectrum, H)`. It looked derivable.
- **`n = 7`: `|R` DECOUPLES.** Explicit verified witness — two non-isomorphic 7-tournaments with
  identical char poly `x⁷−9x⁴−12x³−16x²−8x−1` and identical `H = 81`, but `|R| = 1` vs `17`. So
  `(spectrum,H,|R|)` *strictly refines* `(spectrum,H)`: **`|R|` carries information that neither
  the spectrum nor `H` captures.** The `n=6` coupling was a small-`n` coincidence.

The lesson is a familiar one in this project (per-path identity `n≥6`, `H`-split `n=6`): **do not
read a small-`n` coincidence as a law.** `|R|=f(spectrum,H)` held exhaustively at `n≤6` and *broke*
at the very next `n`. Only the explicit `n=7` witness settles it.

## Why `|R|` and `H` are the signed/unsigned pair

Both are Hamiltonian-path statistics; both leave the spectral floor at exactly `n=6` (the #P
threshold of THM-1780/1870); at `n=6` they are coupled, from `n=7` independent. Equivalently the
2-D invariant `(#even-sign, #odd-sign Ham paths) = ((H+|R|)/2, (H−|R|)/2)` is the natural object,
and it is strictly finer than `H`. Concretely `(H,|R|)` distinguishes **31 of 56** iso classes at
`n=6` — beating `H` alone (**19**) *and the full spectrum* (**28**). Two combinatorial counts beat
the linear-algebra spectrum.

## The extremes

- **Maximizer (handoff 2):** `max|R|` at `n=7` is `147`, achieved by the **QR(7) Paley (regular)
  tournament** (`H=189`). The sequence `3,3,15,15,147` is *not* the double factorial (`7!!=105`).
- **Strong-atom spectrum (handoff 1):** strong tournaments realize `|R|` in
  `{3},{1},{3,5,7,11,15},{1,3,5,7,9,11,13}` for `n=3..6`; strong-`6` caps at `13` while
  decomposable `5▷1` reaches `15` (THM-1936 factorization). Characterizing the strong-atom
  spectrum is the residual open thread.

## Takeaway

`|R|` earns its place in the invariant zoo: the **signed Rédei count is a genuinely new tournament
invariant from `n=7`**, the signed partner of `H` on THM-1780's ladder, and `(H,|R|)` is a
strictly finer universal code than the spectrum. It should be added to the WOWII/zoo invariant set
(klein-S399) and to the H-spectrum "universal code" fingerprint (engineering mandate).

→ THM-1966; builds on THM-1936 (R join-multiplicative), THM-1780 (H-ladder), THM-1870 (cycle
counts spectral); descends from Rédei.
