# Character vs interval: one Fourier duality on Z_p behind both max-H and LRC, graded by χ(−1)

*opus-2026-06-29. Owner asked to prove the Paley coincidence ends at 19, find the maximizer for
larger n, and think Fourier-abstractly toward LRC. The proof is a clean exact computation; the
"why" is a duality between the two canonical functions on Z_p — the multiplicative character
(Gauss-flat) and the interval/sign (Dirichlet-concentrated) — and it unifies this whole arc.*

## PROOF: Paley is not the H-maximizer at n=19 (exact, two routes)
Both via independent Held–Karp routes (paths-from-0 ×n, and paths-from-all-starts):
```
H(Paley_19)      = 1,172,695,746,915   (matches the dihedral-session count)
H(half-turn_19)  = 1,184,212,824,763
difference       = +11,517,077,848  > 0.
```
A second tournament beats Paley, so **Paley is not the maximizer at n=19. QED.** The "coincidence"
(Paley = maximizer) holds for n=3,7,11 and **fails at n=19** — even though Paley is genuinely
doubly-regular there. (Independently corroborates LEM-004's "p=19 boundary.")

## The Fourier structure of a circulant tournament (clean, proved)
A circulant tournament ⇔ an **odd sign function** `s: Z_n→{±1}` (`s(−d)=−s(d)`), `C={d:s(d)=+1}`.
Its adjacency eigenvalues are the Fourier transform of the connection set, and they are **rigidly
constrained**:
> `μ_0=(n−1)/2`, and for `j≠0`, `μ_j = −1/2 + (i/2)·ŝ(j)` — **all nonzero eigenvalues lie on the
> vertical line `Re=−1/2`** (verified n=7,11,19,23), with imaginary part `= ½ ŝ(j)` (`ŝ` purely
> imaginary, since `s` is odd-real).

So a circulant tournament *is* its odd function `s`, and its "shape" is `ŝ`. Two canonical choices:

| tournament | `s` | `ŝ` | spectrum `|Im μ_j|` (j≠0) |
|---|---|---|---|
| **Paley** | Legendre character `χ` | **Gauss sum** `g=√(χ(−1)p)` | **FLAT** `=√n/2` (Ramanujan-extremal) |
| **half-turn** | sign / sawtooth on `{1..(n−1)/2}` | **Dirichlet / cotangent kernel** | **SPREAD**, large at low freq (n=19: 0.08–6.03) |

`χ` is a Fourier **eigenfunction** (`χ̂ = g·χ̄`) — maximally pseudorandom, energy spread flat. The
interval/sign is its opposite — localized in the group, hence spread in frequency (Dirichlet).

## Why the half-turn wins big (the mechanism, honest)
Both regular ⇒ same `t₃`. Power sums `Σμ_j^k = tr(A^k)` give the cycle counts: the **flat** Gauss
spectrum maximizes the **short** cycles (`t₅,t₇` — Paley wins those), which dominate `H` at small n.
But `H=1+2Σαₖ` is ruled by **disjoint odd-cycle families** for large n, and the half-turn's
**spread** (large low-frequency) spectrum — the antithesis of quasi-random — boosts `H` above the
quasi-random value that Paley (flat/Ramanujan) sits near. So: **max-H selects spectral SPREAD
(localized/Dirichlet), not flatness (Gauss/quasi-random).** Pseudorandomness is the *wrong* target.
The precise spectral law for the disjoint-family term is the open piece; the crossover (Paley→half-turn)
is at n=13–19.

## Larger n (status)
Half-turn is the **full circulant maximizer** at n=13,15,17 (exhaustive) and beats Paley + samples at
n=19. Conjecturally it continues (the spread-spectrum argument), but neither "half-turn is the global
maximizer" nor "Paley never returns for n>19" is proven — both are open. `R_n` (half-turn H) is the
natural sequence: `…,3711175,198464295,13689269499,1184212824763,125547534942879` (n=13..21).

## The grand Fourier unification (the abstract picture)
Everything in this arc is the **Fourier transform of a ±1 function on `Z_p`**, and the **same
discriminant `χ(−1)=(−1)^{(p−1)/2}=ε`** governs it all:
- `ε` = whether the **Gauss sum** is real (`√p`, p≡1) or imaginary (`i√p`, p≡3) = the master
  discriminant `(−1)^{C(p,2)}` of the reversal/palindrome session.
- `ε` = Paley reflection automorphism (p≡1) vs **anti**-automorphism (p≡3) = the dihedral session.
- **Two poles of the duality:** the **character/Gauss** pole (flat, pseudorandom, imaginary, the LRC
  **sign-obstruction**, Paley) vs the **interval/Dirichlet** pole (concentrated, ordered, the LRC
  **unsafe band `φ̂`**, the half-turn, the SOS/provable side).
- **Max-H** picks the interval/Dirichlet pole. **LRC's hardness** is the character/Gauss pole. They
  are the *same* duality, optimized in opposite directions.

So "why Paley/circulant/nothing-else" resolves into: the maximizer is the **Dirichlet/interval**
circulant; Paley (the **Gauss/character** circulant) only wins in the small-n short-cycle regime; and
the whole question is one face of the character-vs-interval Fourier duality on `Z_p` that also carries
the LRC obstruction on its other face — both keyed to `χ(−1)=ε`.

## Status
- **Proved:** Paley not max at 19 (exact); eigenvalues on `Re=−1/2`, `Im=½ŝ`; Paley flat (Gauss),
  half-turn spread (Dirichlet).
- **Surmise/open:** max-H ⇔ spectral spread (disjoint-family law); half-turn = global/asymptotic
  maximizer; Paley never returns.

Artifacts: `04-computation/maxH_fourier_spectrum_opus_20260629.py` (+ `.out`),
`maxH_why_mechanism_opus_20260629.py`. Related: THM-027 (H=trM), THM-374 (half-turn), THM-088/127
(ε, χ(−1)), the half-principle + free-monoid reflections, OPEN-Q-108.
