# The n=4 Klein-four is the shadow of a free monoid: tournament iso-classes factor into strongly-connected primes, the complement is word-reversal, and H is multiplicative — a recursive compression that names the incompressible cores

*opus-2026-06-29. The owner gave the n=4 tiling/half-tiling Klein-four `{T,+,-,S}` and asked to
extend the compression hard and recursively. It extends all the way: the Klein-four is the
weight-4 slice of a universal **free monoid on the strongly-connected tournaments**, and that
monoid IS the recursive compression. Everything connects back to the complement/reversal
involution `τ` of the last three sessions.*

## The owner's two structures (verified)
- **Flip-composition group** (scheme 1→2): the 4 iso-classes at n=4 form the **Klein-four**
  `(Z/2)²`, `T`=identity, `S=+·−`; scheme-2's 2-bit cross-section bijects with them; scheme-1's
  3rd bit (`c≈ab`) is redundant (the 4 extra `S`-witnesses). **Complement = the `(+↔−)`
  automorphism**, fixing `{T,S}`. So `{T,S}` (= the SC self-converse classes) is the `+1`
  eigenspace of complement and `{+,−}` (the NS pair) is the `−1` eigenspace: **the Klein-four is
  the complement-involution's eigenspace split made into a group.**

## The deep extension (verified): a FREE MONOID on strongly-connected primes
Every tournament decomposes **uniquely** as a transitive chain of its strongly-connected
components (the condensation). So under the composition `X⊕Y` ("X dominates Y"):

> **Tournament iso-classes = the free monoid on the strongly-connected "primes."**
> `#iso(n) = Σ_{compositions of n} ∏ SC(partᵢ)`, `SC(n)=1,0,1,1,6,35,353,…` (the irreducibles).
> Verified exactly: free-monoid count `= 1,1,2,4,12,56` `= A000568(1..6)`.

The Klein-four is just the weight-4 words in primes `{a(1), C₃(3), S(4)}`:
`T=[1,1,1,1]`, `+=[3,1]`, `−=[1,3]`, `S=[4]`. So **`+↔−` is literally word-reversal**, and
`(X⊕Y)^op = Y^op⊕X^op` ⟹ **complement `τ` = reverse the word + complement each letter**;
**SC (self-converse) iso-classes = the `τ`-palindromic words.** This is exactly the
reversal/palindrome involution of the [[the-reversal-involution-and-the-half-principle]] session:
the free-monoid word is the tiling; the palindromic word is the blue/grid-symmetric/SC tiling;
the "half" that determines a palindrome is the compression.

## The recursive compression and its practical payoff (verified)
`H` (Hamiltonian-path count) is **multiplicative over the condensation**:
> `H(X⊕Y) = H(X)·H(Y)` (once a Ham path enters the dominated block it cannot return).
> Verified for all 56 classes at n=6. So **`H(T) = ∏ H(its SC primes)`** — a *prime
> factorization of the invariant*. The "H-atoms" are the SC primes; their H-values are all
> **odd** (3, 5, then {9,11,13,15} at n=5, {15..45} at n=6).

**Practical recursive applications:**
1. **Compression / canonical code.** Encode a tournament by its condensation word + each
   prime's code, recursively. Reducible structure costs ~0 bits; only the SC primes carry
   information. The transitive tournament is maximally compressible (all trivial primes); a
   single large SC tournament (e.g. the regular/resonant one) is **incompressible** — it *is* a
   prime. This is the project's "universal tournament code," made into a factorization.
2. **Fast invariants.** Any condensation-multiplicative invariant (`H`, and — known — the full
   independence/forward polynomial up to a shift) is computed by factoring + multiplying atoms.
   This is the right home for the engineering "H-spectrum fingerprint."
3. **The hardness gradient is the prime size.** `SC(n)/A000568(n) → 1`: almost all tournaments
   are single primes, so compression helps the reducible tail; the irreducibles are the cost.

## Surmise: the same compression is the LRC proof shape (free, unproven)
The whole project's complement/reversal `τ` (= `T^op` = path reversal = `t↦−t` on the LRC
circle = `+↔−` here) and its sign `ε=(−1)^{C(n,2)}` are one object seen many ways. The
compression principle then reads, on the LRC side:

- **Factor the cover/cap into "transitive (reducible)" ⊕ "strongly-connected (irreducible)"
  parts;** the reducible part multiplies out (provable, even/`τ`-symmetric), and the
  **obstruction lives in the irreducible `τ`-odd residual** — the resonant/AP core, which is the
  LRC analogue of a large SC prime. This is exactly what every route this week found (the
  bounded resonant core, the sign-isotypic obstruction). The "half"/score compression keeps the
  `τ`-odd residual and discards the rest — and *that residual is the thing to certify.*
- Concretely testable next step: define a condensation-style decomposition of a covering speed
  set and check whether the cap is multiplicative over it, isolating an irreducible core whose
  size is bounded (the THM-527 bounded-core, now as a "prime").

## Status
- **Verified:** free-monoid `=A000568`; Klein-four = weight-4 slice; `+↔−`=word-reversal;
  complement = reverse+complement-letters; SC = palindromic words; `H` multiplicative over
  primes (all n≤6); H-atoms odd.
- **Synthesis:** the Klein-four is one slice of the universal free monoid; the complement/
  reversal `τ` and `ε` unify this with the dihedral, palindromic, and LRC threads; compression =
  factor-into-primes = quotient-by-`τ`-keeping-the-odd-residual.
- **Surmise / open:** the LRC "cap is multiplicative over a condensation, obstruction = the
  irreducible `τ`-odd prime." A concrete, falsifiable program.

Artifacts: `04-computation/{iso_compression,condensation_monoid,H_factorization}_opus_20260629.py`,
`05-knowledge/results/iso_class_condensation_monoid_opus_20260629.out`. Related: THM-084/088/051
(reversal), THM-127/131 (dihedral), the half-principle reflection, h21-s617 (H over strong
components, Moon), OPEN-Q-108.
