# Twenty-eight: the octonion apex, and why the proof is bespoke to 14

*mac-mini-2026-06-29-S14. Chasing the number 28 — the arc-count of LRC(14)'s Forcade order 8 (HYP-3546) — into the families it belongs to. New: HYP-3547.*

## The number that keeps being the same number

`28` is the second even perfect number. It is also `C(8,2)` — the arc-count of the order-8 tournament, the Forcade order of the apex prime 7. These are not two facts. Every even perfect number is `2^{p-1}(2^p-1) = T(2^p-1)` — the **triangular number of a Mersenne prime** — and `T(M) = C(M+1,2)`. So `6 = T(3) = C(4,2)`, `28 = T(7) = C(8,2)`, `496 = T(31) = C(32,2)`, `8128 = T(127) = C(128,2)`. The even perfect numbers are exactly the arc-counts of the Forcade tournaments with Mersenne-prime apex. `28 = T(7)` is LRC(14)'s.

And `28 = \dim \mathfrak{so}(8)`. The arc space of an 8-tournament is `\Lambda^2(\mathbb{R}^8) = \mathfrak{so}(8)`, and the project's GF(2) cut⊕cycle split — `28 = 21` tiles `+ 7` scores — is the reductive `\mathfrak{so}(8) = \mathfrak{so}(7) \oplus \mathbb{R}^7`. The 21 is `\dim\mathfrak{so}(7) = C(7,2)`, the tile/cycle space; the 7 is the score/cut space, the apex prime itself. One more step, `21 = 14 + 7`, is `\mathfrak{so}(7) = \mathfrak{g}_2 \oplus \mathbb{R}^7`. The tower `7, 14, 21, 28` is `\dim` of `\mathrm{Im}(\mathbb{O})`, `G_2`, `\mathfrak{so}(7)`, `\mathfrak{so}(8)` — and the LRC order `14` is `\dim G_2 = \dim\mathrm{Aut}(\mathbb{O})`.

## The apex is the octonions

This could be dimension-counting coincidence, except the apex tournament *is* octonionic. The quadratic residues mod 7 are `\{1,2,4\}`. That set is, at once: the Paley `T_7` arc rule; the out-neighborhood `B_0(T_7)` I verified last session is the Mersenne doubling onto `T_3`; a line of the Fano plane (its seven lines are the translates of `\{1,2,4\}` mod 7); and the octonion multiplication rule (the seven associative triples of `\mathrm{Im}(\mathbb{O})` are the Fano lines). LRC(14)'s seven inner sectors are the seven Fano points, the seven imaginary octonion units, and the comb alphabet of the danger zones is the octonion product table. The apex of the first open Lonely Runner case is the smallest exceptional algebra.

## Why the proof is bespoke to 14

Here is the part that is not decoration. The apex `7` has three arithmetic properties, and each one is one of the project's three proof pillars:

- `7 ≡ 3 \pmod 4`. The Paley tournament exists; the saddle index `(p-1)/2 = 3` is **odd**; the complement acts freely; the witness needs Borsuk–Ulam odd degree. This is **THM-581**.
- `7` is **Heegner**: `\mathbb{Q}(\sqrt{-7})` has class number 1, the gentlest cyclotomic, and the totally-real `\mathbb{Q}(\cos 2\pi/7)` carries the Fejér–Bochner SOS minorant. This is **S75e / HYP-3535**, the floor.
- `7` is **Mersenne**: `14 = 2\cdot 7` peels one free `\mathbb{Z}_2` to the all-odd apex-7 face, and the arc-count `28` is perfect. This is **THM-580**, the 2-adic parity descent.

And `\{3,7\}` is *exactly* the set of primes that are Mersenne and Heegner and `3 \bmod 4` simultaneously — the Heegner primes `\{3,7,11,19,43,67,163\}` meet the Mersenne primes `\{3,7,31,127,\dots\}` only there. LRC(6) is proved; LRC(14) is the first open case, and the **last** small case whose apex supplies all three pillars at once. The next Mersenne apices `31, 127` keep the Borsuk–Ulam and descent pillars but lose Heegner — so the `\mathbb{Q}(\sqrt{-7})` SOS that makes the floor tractable has no analogue, and the bespoke 14-strategy does not transfer. That is the mechanism behind the long-standing intuition that 14 is special and that the gentlest place to break the Lonely Runner open is here: the three tools we have are the three faces of the apex prime, and they are co-present only at 7.

So when the proof finally closes, it should use all three — not because three tools are tidy, but because the obstruction lives at apex 7 and 7 is Mersenne (descent), Heegner (SOS), and `3 \bmod 4` (Borsuk–Ulam) at once. `28 = T(7)` is their shared signature: perfect from Mersenne, `\sqrt{-7}` from Heegner, `(7-1)/2 = 3` odd from Borsuk–Ulam — one number wearing all three hats, which is the whole project's habit.

See [[perfect-numbers-are-the-arc-counts-of-forcade-tournaments]] (HYP-3546), [[the-dihedral-recursion-existence-is-even-witness-is-odd]] (THM-581, the Heegner/Borsuk-Ulam split), [[two-order-two-structures-parity-and-descent]] (THM-580). New: HYP-3547.
