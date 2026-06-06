# HYP-2275 — Goldbach/Lemoine: even=unordered (σ-symmetric), odd=ordered (σ-broken), and 6=2·3 the meeting point

**Session:** claudebox-2026-06-03-S630. **Prompt (user):** even = sum of two primes (unordered pairs); odd = prime +
doubled prime (ordered pairs); which evens/odds map to the same pairs; the duplicates (p,p) where a prime doubles to
an even and triples to an odd. **Threads:** HYP-2270 (add/mult exclusion), HYP-2185 (apex=σ-fixed), HYP-2265 (π/3=6),
HYP-2175 (Collatz 2-vs-3).

## The two faces (formalized)
- **Even = Goldbach `p+q` = UNORDERED pair** `{p,q}`: the swap `σ:(p,q)↦(q,p)` FIXES `p+q` (symmetric) — the
  **additive, σ-symmetric** face.
- **Odd = Lemoine/Levy `p+2q` = ORDERED pair** `(p,q)`: the doubling `2q` BREAKS σ (`p+2q ≠ q+2p` off-diagonal) — the
  **multiplicative, σ-broken** face. Formalized: `lemoine_swap_symm_iff_diagonal` (`p+2q = q+2p ↔ p=q`),
  `lemoine_parity` (`(p+2q) % 2 = p % 2` — the doubling is invisible to parity, only the first prime sets it).

This is exactly the S629 additive/multiplicative complementarity and the perspective key: unordered/ordered =
σ-symmetric/σ-broken; the doubling (×2, the 2-adic) is the symmetry-breaker (as in LRC/Collatz).

## The diagonal = the apex (σ-fixed): 2p and 3p
The duplicate pairs `(p,p)` are the σ-FIXED points (the apex of HYP-2185): `p` doubles to `2p` (even, the Goldbach
`p+p`) and triples to `3p` (odd, the Lemoine `p+2p`). The 2-face and the 3-face, on the diagonal — the same 2 vs 3 as
Collatz (×3 vs ÷2) and the π/3 cube root.

## THE MEETING POINT: 6 = 2·3 = 3·2 (formalized)
**`2p = 3q` for primes `p,q` ⟹ `p=3, q=2`**, value **6** — the UNIQUE number that is both a doubled prime and a
tripled prime (verified + formalized `double_eq_triple_unique`). `6 = 2·3 = 3·2` is the single point where the
doubling-face and the tripling-face COMMUTE — and it is the SAME 6 as the hexagonal/Eisenstein lattice, the LRC gap
`δ=1/6`, the chord `dZ=1/6` (S623), the cube-root angle `π/3`, and `Φ₃` (S628). The whole arc's "2 and 3" converge on
`6 = lcm(2,3) = 2·3`, and here it is the prime `(3,2)` where doubling meets tripling.

## The collision structure (which evens/odds share pairs)
An even `E = p+q` shares its pair `{p,q}` with the two odds `E+p = 2p+q` and `E+q = p+2q` (the Lemoine reps of the
same pair). So `O − E ∈ {p,q}` (a Goldbach prime of `E`). The diagonal even `6=(3+3)` links to the single odd
`9 = 3·3` (the triple); `p=2` is special (`3·2=6` even — the only even triple, the meeting point again).

## Formalized (math-lean, sorry-free) — `Math/NumberTheory/GoldbachLemoine.lean`
`double_eq_triple_unique` (`2p=3q ⟹ p=3∧q=2`, the 6 meeting), `six_double_three_eq_triple_two`,
`lemoine_swap_symm_iff_diagonal` (doubling breaks swap-symmetry except on the apex), `lemoine_parity`.

## The synthesis
Goldbach/Lemoine is another face of the arc's spine: even/odd = additive/multiplicative (S629) = σ-symmetric/σ-broken
(the perspective key) = unordered/ordered; the doubling (2-adic) breaks the symmetry; the diagonal (p,p) is the apex;
and `6 = 2·3` is the unique commutation point of the 2-face and the 3-face — the hexagonal/π-3/dZ=1/6 of the unit
distance, LRC, and tournament threads. Goldbach and Lemoine are the additive (even) and multiplicative (odd)
conjectures, and the diagonal/apex/6 is where they meet.

## Open
- Is the shared-pair "collision graph" (evens linked to odds via `O−E` prime) a known structure? Its components.
- The Goldbach (σ-symmetric, even) vs Lemoine (σ-broken, odd) as the additive/multiplicative conjecture pair — a
  formal duality? (Both are "every n>bound has a representation"; the σ-orbit structure.)
