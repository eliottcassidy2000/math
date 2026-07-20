---
id: THM-1735
title: "THE FINITE-PLACE HALF OF TNC IS CLOSED (and it CORRECTS THM-1730): for a tunable trinomial, CT(m0) and CT(2m0) are coprime mod EVERY prime p not dividing Res_a(CT(m0),CT(2m0)) -- so, since Res != 0 (THM-1680/1710), ALL BUT FINITELY MANY finite places separate, with an EXPLICIT smallest good prime. The correct criterion is gcd(CT(m0),CT(2m0)) = 1 in F_p[a] (Lucas/carry-CA reductions), NOT the coarse Newton-polygon root-valuation used in THM-1730 s1: that coarse criterion falsely reported {-2,1,4} as needing the archimedean place, but gcd(3(a^2+1), 15(a^4+4a^2+1)) = 1 mod 7 (Res = 72900 = 2^2 3^6 5^2, bad primes {2,3,5}, so p=7 separates). Verified: {-2,1,4} good at p=7, {-2,3,6} Res=1 (coprime mod EVERY p), {-3,-1,3} p=3, {-3,1,5} p=3, {-3,2,7} p=11, {-4,1,6} p=7. THE CELLULAR-AUTOMATON CONTENT: CT(m) mod p is a LUCAS product of base-p digit multinomials = the Sierpinski/Pascal carry CA (HYP-2491), and the BAD primes are exactly the prime factors of the resultant. So the finite-place half reduces to: bound the resultant's prime factors by the charge data (a carry bound), making the separating prime explicit. INTEGRATES kind-pasteur THM-1720: the recurrence's structural leading-coefficient roots are roots-of-unity (monodromy) exponents, never positive integers, so the recurrence is nonsingular at every m>=1 -- the same nondegeneracy on the recurrence side"
status: PROVED (finite-place separation at every p not dividing Res, given Res != 0 from THM-1680/1710) + VERIFIED (6/6 tunable trinomials, explicit good primes). Corrects THM-1730 s1.
author: opus-2026-07-20-S426
depends_on: [THM-1680 (trinomial gcd, Res != 0), THM-1710 (resultant), THM-1730 (adelic synthesis -- s1 corrected here), kind-pasteur THM-1720 (recurrence structural roots), HYP-2491 (Sierpinski carry CA)]
corrects: THM-1730 §1 (the {-2,1,4} archimedean-only claim)
---

# THM-1735 — The finite-place half is closed

## 0. Correction to THM-1730 §1

THM-1730 §1 used a **coarse criterion** — disjointness of the `p`-adic Newton-polygon *root
valuations* — and reported `{−2,1,4}` as separating "only at the archimedean place." **That
was wrong.** The Newton-polygon valuations of both `CT(m_0)` and `CT(2m_0)` are `0` at every
`p` for `{−2,1,4}` (their leading and trailing coefficients are equal in magnitude), so the
coarse test never fires — but it is only *sufficient*, not necessary. The **exact** criterion
is gcd:

```
gcd( 3(a²+1),  15(a⁴+4a²+1) )  =  1   in 𝔽_7[a] .
```

`{−2,1,4}` separates at the **finite** place `p = 7`. No archimedean place is needed.

## 1. The finite-place theorem (PROVED)

> **Theorem.** For a tunable trinomial, `CT(Λ^{m_0})` and `CT(Λ^{2m_0})` (a-contents) are
> **coprime in `𝔽_p[a]` for every prime `p ∤ Res_a(CT(m_0), CT(2m_0))`.**

*Proof.* Reduction mod `p` commutes with the resultant: `Res(\bar f, \bar g) ≡ Res(f,g)
\pmod p` (when leading coefficients survive). If `p ∤ Res(f,g)` then `Res(\bar f, \bar g) ≠
0` in `𝔽_p`, so `\bar f, \bar g` share no root in `\overline{𝔽_p}` — coprime. ∎

Since **`Res(CT(m_0), CT(2m_0)) ≠ 0`** (THM-1680/1710), its integer value has finitely many
prime factors, so **all but finitely many primes separate the two levels** — the finite-place
half of HYP-8530 is complete, with a concrete certificate.

**Verified**, with the smallest good prime:

| pattern | `|Res|` | bad primes `p\|Res` | smallest good `p` |
|---|---|---|---|
| `{−2,1,4}` | `72900 = 2²3⁶5²` | `2,3,5` | **7** |
| `{−2,3,6}` | `1` | — | any (coprime mod every `p`) |
| `{−3,−1,3}` | `2·7·101·…` | `2,7,101` | **3** |
| `{−3,1,5}` | `2⁷5²7·…` | `2,5,7` | **3** |
| `{−3,2,7}` | `2…3…5…` | `2,3,5` | **11** |
| `{−4,1,6}` | `2…3…5…2287` | `2,3,5,2287` | **7** |

## 2. The cellular-automaton content (Lucas = Sierpiński carry)

The reductions are not opaque: **`CT(Λ^m) \bmod p` is a Lucas product**. By Lucas' theorem a
multinomial `\binom{m}{x,y,z} \bmod p` factors over the base-`p` digits of `m, x, y, z`, and
this digit product is exactly the **Sierpiński/Pascal carry cellular automaton** (Rule 90 at
`p = 2`, HYP-2491/2497). So:

> The gcd `CT(m_0) \bmod p`, `CT(2m_0) \bmod p` is computed by the carry CA, and **the bad
> primes are precisely the prime factors of the resultant** — the places where the two CA
> reductions acquire a common root.

This turns the remaining question into a **carry bound**: *bound the prime factors of
`Res(CT(m_0),CT(2m_0))` by the charge data `(N; j, d)`.* The resultant is a product of root
differences `∏(α_i − β_j)` (roots of the two levels), each an algebraic multinomial-ratio; a
Kummer/Lucas estimate of its height would bound the bad primes and make the smallest good
prime **explicit in the pattern** — closing the finite-place half *uniformly*, not just
per-pattern.

## 3. Integration with kind-pasteur THM-1720 (the recurrence side)

kind-pasteur's THM-1720 attacks the same object from the recurrence: the order-`D` toral
recurrence's leading coefficient factors `P_D(m) = STRUCTURAL(m)·APPARENT_R(m)`, and the
**structural roots are root-of-unity (monodromy) exponents** — negative rationals with
denominators dividing `M, N`, hence **never positive integers**, so `STRUCTURAL(m) ≠ 0` for
all `m ≥ 1` unconditionally. That is the *recurrence-side* nondegeneracy; §1 here is the
*coefficient-side* nondegeneracy (coprimality mod `p`). Both say the two levels cannot align:
one via the recurrence's leading coefficient, one via the resultant's `p`-adic reduction. The
shared object is the **roots-of-unity branch monodromy** — kind-pasteur's structural roots
and my Kummer/Lucas reductions are two readouts of the same `μ_M × μ_N` covering symmetry.

## 4. Status of the adelic completion

| half | status |
|---|---|
| **finite places** | **CLOSED** (§1): coprime mod every `p ∤ Res`, given `Res ≠ 0` |
| archimedean | not needed (§0 correction); the amoeba is one alternative certificate, superseded by finite places |
| uniform bound | open: bound `p \| Res` by `(N;j,d)` via a carry/height estimate (§2) |

**Net:** the finite-place half of the adelic picture closes the trinomial TNC with an
explicit finite certificate, resting only on `Res ≠ 0` (THM-1680) and Lucas reduction. The
one remaining refinement is making the good prime uniform in the pattern.

## 5. Next

1. **Carry bound on the bad primes.** `p \| Res` ⟹ `CT(m_0), CT(2m_0)` share a root mod `p`;
   bound such `p` by the multinomial heights (Kummer). A clean `p ≤ f(j,d,N)` closes the
   uniform finite-place statement.
2. **k-nomial.** The same reduction: the `(k−1)` levels are coprime (empty variety) mod every
   `p` not dividing the elimination resultant; the Lucas/carry CA computes the reductions.
3. **Merge with kind-pasteur's APPARENT factor.** Their residual (apparent factor has no
   positive integer root) and my resultant bad primes are the same obstruction from two
   sides; a joint desingularization statement would close both.

## Verification

`04-computation/tnc_finite_place_modp_opus_S426.py` (resultant bad-prime table, smallest good
prime per pattern), `04-computation/tnc_modp_correction_opus_S426.py` (the `{−2,1,4}` gcd-mod-7
correction). Outputs in `05-knowledge/results/`.
