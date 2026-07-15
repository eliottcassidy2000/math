---
id: THM-871
title: FERMAT-RUNG RIGIDITY AND THE SEDENION CENSUS — for prime n with n−1 a 2-power (Fermat primes: the CD rungs 3, 5, 17, 257, …), the multiplier group Z_{n−1} is a 2-group, so (tournament Aut being odd) every rotational tournament has Aut = Z_n exactly; at n = 17 the Z₁₆ action on the 256 connection sets is FREE (−1 lies in every nontrivial subgroup): exactly 16 rotational classes, all rigid, NO canonical representative — the sedenion rung loses the distinguished unit (rungs 3, 5 have a UNIQUE rotational class; the octonion rung 9 is anti-rigid: C₉{1,3,4,7} = C₃[C₃] with Aut = Z₃≀Z₃, order 81 = n²)
status: PROVED (free action + census: machine-exact; Aut = Z_n via Alspach's theorem for prime-order vertex-transitive tournaments [classical] + the odd-Aut constraint; wreath identification at n = 9: |Aut| = 81 computed by full 9!-sweep, block partition preserved by all 81)
source: klein-2026-07-15-S313 (cont.3); answers THM-868 "Named next steps" item 2 (sedenion rung degeneration)
depends_on: [THM-868 §5 (Feit–Thompson face), Alspach 1970 (Aut of prime-order circulant tournaments = Z_p ⋊ M), Elspas–Turner (Ádám conjecture at prime order)]
verification: 04-computation/icosian_kakeya_milgram_sedenion_klein_S313.py (sections B; 24/24); c4-invariants distinguish ≥ 5 of the 16 classes directly
---

# THM-871 — Fermat-rung rigidity; the sedenion rung census

## 1. The rigidity mechanism (Feit–Thompson leverage)

For a rotational (circulant) tournament on prime p, Aut = Z_p ⋊ M with M ≤ Z_{p−1}
the multiplier stabilizer (Alspach). Tournament automorphism groups are odd, so M
lies in the odd part of Z_{p−1}. **If p − 1 is a 2-power — p a Fermat prime, the CD
rungs 3, 5, 17, 257 — the odd part is trivial: Aut = Z_p exactly. Every rotational
tournament at a Fermat rung is rotation-rigid.**

## 2. The sedenion census (n = 17)

The 256 connection sets (one choice per ±pair) carry the Z₁₆ multiplier action.
Z₁₆ has a UNIQUE involution (−1), contained in every nontrivial subgroup; a set
fixed by any m ≠ 1 would satisfy S = −S — impossible. So the action is **FREE**:

> exactly **16 rotational classes at n = 17**, a free Z₁₆-torsor — no fixed points,
> no distinguished class. (By Ádám-at-primes these are honestly 16 iso classes;
> the c4 invariant already separates values 1105…1428 across them.)

Contrast along the CD tower:
| rung | 3 (C) | 5 (H) | 9 (O) | 17 (S) |
|------|-------|-------|-------|--------|
| rotational classes | 1 (C₃) | 1 (C₅ round) | many | 16 |
| max rotational Aut | Z₃ | Z₅ | **Z₃≀Z₃ (81 = n²)** | Z₁₇ (all rigid) |

At rungs 3 and 5 there is a UNIQUE rotational tournament (a canonical unit). At the
octonion rung 9 (composite: 9 − 1 = 8 but 9 = 3² is not prime) rigidity fails
maximally: C₉{1,3,4,7} **is** the wreath C₃[C₃] (blocks {0,3,6},{1,4,7},{2,5,8}),
|Aut| = 81 = Z₃≀Z₃, the largest of any 9-tournament type — extra symmetry, not less.
At the sedenion rung 17 rigidity returns but uniqueness is destroyed: 16 anonymous
rigid classes.

> **The sedenion degeneration, tournament-side: the loss is not of symmetry but of
> CANONICITY — the first CD rung with no distinguished rotational unit.** (Echo of
> the algebra: sedenions keep the involution structure but lose the norm that
> singles units out; flagged as reading, the census is the theorem.) At 257 the
> same mechanism gives 2^128/256 = 2^120 anonymous rigid classes.

Note the QR/Paley thread: Paley tournaments need q ≡ 3 (mod 4); Fermat primes ≥ 5
are ≡ 1 (mod 4), so no Paley anchor exists at rungs 5, 17, 257 — consistent with
anonymity. (At rung 3, QR₃ = C₃ is the unique class.)
