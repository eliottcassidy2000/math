# HYP-2327 — LRC at 19 (prime) vs 14 (composite): the fiber toolkit fails, the provable corner generalizes, and 19's leverage is Heegner √−19

**Session:** S649
**Status:** CONFIRMED (consecutive config formalized for every n incl. 19; the prime/Heegner structure mapped)
**Provenance forward:** math-lean `Math/LonelyRunner/LonelyNineteen.lean` (sorry-free)
**Prompt:** work on the full LRC, but 19 not 14.

---

## 0. Why 19 is a different beast from 14

`14 = 2·7` is **composite**, and the whole S640 attack was to fiber it via CRT `ℤ/14 ≅ ℤ/2 × ℤ/7` over
the 7-clock (the doubling `⟨2⟩` is `ord₇(2) = 3` = two 3-cycles — small, nested, *reducible*). **`19` is
prime**, so:

```
  divisors(19) = {1, 19}        → NO nontrivial CRT fiber.
  ord₁₉(2) = 18 = 19 − 1        → 2 is a PRIMITIVE ROOT; doubling = a single 18-cycle (maximal mixing).
  2·19 − 1 = 37 (prime).
```

So **the n=14 divisor/fiber toolkit does not transfer** — prime `n` is exactly the hard case for
reduction methods (no sub-shell, no fiber-bundle, the doubling tower is one long cycle). This is the
honest headline: *19 is harder than 14 for the structural attack, not easier.*

---

## 1. The provable corner (formalized for every n, incl. 19)

What survives is the **consecutive / tight-extremal configuration**, and it generalizes cleanly. At
`t = 1/n`, runner `k ∈ {1,…,n−1}` has `k/n ∈ [1/n, 1 − 1/n]`, so `dZ(k/n) ≥ 1/n` — the config is lonely.

**Formalized (math-lean, sorry-free): `Math/LonelyRunner/LonelyNineteen.lean`**
- `consecutive_lonely : 1 ≤ k → k + 1 ≤ n → 1/n ≤ dZ(k/n)` — LRC for `{1,…,n−1}` at `t = 1/n`, **every
  `n`** (generalizing S648's concrete `n = 14`; `dZ_ge_of_mem` + `gcongr` + `div_le_one`).
- `lonely_nineteen : 1 ≤ k → k ≤ 18 → 1/19 ≤ dZ(k/19)` — **a machine-checked LRC proof for the canonical
  19-runner config `{1,…,18}`**, witnessed by `t = 1/19`.

Verified (`lrc19_prime_structure_s649.py`): `{1,…,18}` at `t = 1/19` has min clock distance **exactly
1/19** (equality on the unit runners `k = 1, 18`) — the bound is *achieved, not beaten*: the tight
extremal. (And the friendliest config, S647: lonely set = the single instant `1/19`, measure zero.)

**HONEST SCOPE.** This is the one tight config. The **full** LRC for all 18-speed sets is open (proven
only up to 7 runners, Barajas–Serra). A random 18-speed search easily finds configs reaching gap `0.216`
≫ `1/19`; the conjecture is the *uniform lower bound* `1/19` for **every** config — untouched here.

---

## 2. Where 19's real leverage lives: CM / Heegner, not fibers

`19` is special on the **cube-root / CM side** the arc converges on, not the divisor side:

- **Heegner number.** `19 ∈ {1,2,3,7,11,19,43,67,163}`; `ℚ(√−19)` has class number 1; `19 = 4·5 − 1` is
  the rotation field for Eisenstein norm `N = 5` (HYP-2277). It is the conjectural **`χ = 5` chromatic
  step** in the Moser/Heegner tower `√−3 → √−11 → √−19` (S687/S641). So `n = 19` sits exactly where the
  chromatic-number-of-the-plane tower would take its next rung.
- **Centered hexagon.** `19 = 1 + 6 + 12 = hex(2)` (Eisenstein lattice, radius 2); `2n−1 = 37 = hex(3)`.
- **Paley-19 exists** (`19 ≡ 3 mod 4`, self-converse tournament, S638), `|QR| = 9`. But `2` is a
  *non*-residue and `⟨2⟩ = (ℤ/19)*` (primitive root), so `⟨2⟩ ≠ QR`: **19 is NOT in the `p ≡ 7 mod 8`
  "doubling = Paley" family** (S640, where `7` lived). `19 ≡ 3 mod 8`.

> So `n = 14`'s leverage was its composite divisor `7` (the fiber); `n = 19`'s leverage is its CM field
> `√−19` (the Heegner/χ=5 frontier). The two frontier cases sit on the two halves of the arc — the
> 2-adic/divisor seam (14) and the cube-root/CM seam (19).

---

## 3. The right attack on full LRC(19) (handoff)

Since fibers fail, the prime case wants the **cyclotomic-depth `q*`** approach (opus S704/THM-439): a
rational witness `t = j/q` clearing the floor `1/n`, with the natural moduli `q = 19` and `q = 37 = 2n−1`
(the shell modulus). The witness lives in the abelian `ℚ(ζ_q)`; LRC(19) ⟺ the depth `q*(V) < ∞` uniformly
— and the cube-root/Heegner structure of `√−19` is where any genuine progress for prime 19 should come
from, not the (absent) divisor fiber.

## 4. New threads / handoffs
- **General `{1,…,n−1}` is now formalized for all `n`** (`consecutive_lonely`) — the tight family is done.
- **Full LRC(19):** the `q*` cyclotomic-depth attack at `q ∈ {19, 37}` (S704); does `√−19` (the χ=5 step)
  give a CM handle the way `√−11` did for the Moser spindle?
- **Prime-vs-composite dichotomy:** is LRC genuinely *harder* at primes (no fiber) and the composite
  cases (2p, the S640 family) the tractable frontier? Catalogue which `n` admit a fiber reduction.
