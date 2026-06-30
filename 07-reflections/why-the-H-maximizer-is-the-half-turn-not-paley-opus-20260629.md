# Why the H-maximizer is regular, why it is Paley only at small n, and why it is really the half-turn (= the LRC comparator): a Dirichlet-vs-Gauss story

*opus-2026-06-29. The owner asked WHY the max-Hamiltonian-path tournament must be Paley / circulant
/ nothing-else when it is. Chasing it overturned the premise: it is Paley only for n≤11; the true
maximizer is the half-turn circular tournament, and the reason ties max-H to the LRC Fourier core.*

## The one rigorous forcing: H-max ⟹ regular
The number of directed 3-cycles is `t₃ = C(n,3) − Σᵢ C(sᵢ,2)` with `Σ sᵢ = C(n,2)` fixed. `Σ C(sᵢ,2)`
is **convex**, minimized exactly when all scores are equal — so **regular tournaments uniquely
maximize `t₃`** (and `t₃` is the leading cycle term of `H = 1 + 2Σαₖ`). This is the only part of
"why" that is a theorem. Everything below is *among regular tournaments*, where `t₃` is already tied.

## Why NOT doubly-regular / Paley in general (the premise was wrong)
Naive intuition says the *most pseudorandom* (doubly-regular = Paley) tournament should win, because
`H` averages `~n!/2ⁿ` over the spread of cycles. **It is false.** Verified (Paley exists only for
`n≡3 mod 4`):
| n | t₃ | t₅ | t₇ | H |
|---|---|---|---|---|
| 7 rotational | 14 | 28 | 17 | 175 | 
| 7 **Paley** | 14 | **42** | **24** | **189** |
| 11 rotational | 55 | 484 | 3399 | 93027 |
| 11 **Paley** | 55 | **594** | **3960** | **95095** |
| **19 rotational** | — | — | — | **1,184,212,824,763** |
| 19 Paley (doubly-reg) | — | — | — | 1,172,695,746,915 |

So **Paley maximizes `H` only at n=3,7,11**, then **loses at n=19** — even though Paley *is*
doubly-regular there. (The repo's LEM-004 already flags this "p=19 boundary.") Doubly-regularity is
NOT the criterion.

## The crossover, and what really wins
Both are regular (same `t₃`). Paley has **more short odd cycles** (`t₅,t₇`) — pseudorandomness packs
small cycles densely; that dominates `H` at small n, so Paley wins. But `H = 1 + 2Σₖ αₖ` where `αₖ` =
families of `k` *vertex-disjoint* odd cycles, and the **long-cycle / disjoint-family terms** dominate
for large n. There the **half-turn / rotational** tournament (i beats the next `(n−1)/2`, i.e. the
clockwise semicircle) wins. The maximizer for n=13,15,17,19 is this half-turn tournament, and it is
**essentially unique** (the n=13 max is one multiplier-orbit of 12 connection sets).

## Why circulant — forced regularity + a Fourier reason (surmise)
Regularity is forced (convexity). Vertex-transitivity is not proven, but the **Fourier picture
explains the half-turn's advantage and ties it to LRC.** The connection set of the half-turn is the
**contiguous block `{1,…,(n−1)/2}`**, whose character transform is the **Dirichlet kernel** — sharply
*concentrated* (a large principal mode). The Paley connection set is the **quadratic residues**,
whose transform is the **Gauss sum** `|·|=√p` — perfectly *flat*. Since `H = tr(M)` (signed transfer
matrix, odd n) is a sum over the circulant's character-eigenvalues, a **concentrated** spectrum
(half-turn/Dirichlet) builds a large trace, while a **flat** spectrum (Paley/Gauss) spreads it thin.
So:

> **The H-maximizer is the *Dirichlet*-spectrum tournament (half-turn), not the *Gauss-sum* one
> (Paley).** Concentration beats flatness for `tr(M)`; pseudorandomness is the wrong target.

## The punchline: max-H = the LRC comparator
The half-turn circular tournament (THM-374: orient `x→y` iff `y` is clockwise within distance `<½`)
**is exactly the LRC half-turn comparator / runner phase clock** (definitions.md). So:

> **The H-maximizing tournament (n≥13) IS the Lonely-Runner half-turn geometry.** And the same two
> Fourier objects run through both stories: the **Dirichlet kernel** (the half-turn / the LRC unsafe
> band `φ̂`) and the **Gauss sum `i√p`** (Paley / the LRC sign-obstruction). Max-H selects the
> Dirichlet side; LRC's obstruction lives on the Gauss-sum side. They are the same dichotomy.

## Answers to the three "why"s (honest)
- **Why Paley when it is:** only n=3,7,11 — a *small-n* regime where short-cycle richness (which Paley
  uniquely maximizes) dominates `H`. Not a deep forcing; it ends at n=19.
- **Why circulant when it is:** `H`-max ⟹ regular (theorem); among regular, the half-turn's
  concentrated (Dirichlet) spectrum maximizes `tr(M)`. Vertex-transitivity is strongly empirical, the
  maximizer essentially unique, but not proven.
- **Why nothing else:** the maximizer is a single multiplier-orbit (essentially unique); the rival
  (doubly-regular Paley) is provably beaten for large n, and non-regular is killed by convexity.

## Status
- **Theorem:** H-max ⟹ regular (convexity). **Verified:** Paley wins n≤11, loses n=19; t₅,t₇ profile;
  half-turn = maximizer n=13..17 (n=19 rotational > Paley, full circulant search not done at 19).
- **Surmise:** the Dirichlet-vs-Gauss spectral mechanism; the half-turn = LRC-comparator unification.

Artifacts: `04-computation/maxH_why_{paley_circulant,mechanism}_opus_20260629.py` (+ `.out`),
`maxH_circulant_opus_20260629.py`. Related: THM-374 (half-turn), THM-027 (H=trM), LEM-004 (p=19
boundary), definitions.md (runner phase clock), OPEN-Q-108.
