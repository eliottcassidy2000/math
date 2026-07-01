# "Large primes forced" made RIGOROUS for the single-swap case: for a prime q with (n−1)/2 < q ≤ n−1, NO single-swap set S = ({1,…,n−1}∖{q}) ∪ {g} is tight — proved by two witnesses, (i) if g is not a multiple of q then S has no multiple of q so the q-witness at t=1/q gives M≥1/q>1/n, and (ii) if g=mq (m≥2, the only covering option) then S misses the ±pair {q, n} mod (n+q), so the radius-1 k-witness at rotation a=q⁻¹ mod (n+q), t=a/(n+q), gives M≥2/(n+q)>1/n (every runner v∈{1,…,n−1}∖{q} has va∉{0,q,n}=−{−1,0,1}, and runner mq lands at residue m∉{0,1,n+q−1}) — verified exactly n=14,18 all large primes, m=2,3; the FULL "large primes forced" (multi-swap / general tight S) is exactly the residual case n∈S, where the (n+q)-witness is defeated by runner n hitting residue −1 — that residual is the repo's HYP-3749 CRT-linkage crux, not closed here

*opus-2026-06-30. Owner: comprehensively catalog the repo's number theory + make "large primes forced"
rigorous. The single-swap case is now a clean elementary theorem; the general case is pinned to one residual =
the repo's open crux. Plus: a 58-thread number-theory landscape (search integrated below).*

## Theorem (single-swap large-prime-forced) — RIGOROUS
Let `q` be prime with `(n−1)/2 < q ≤ n−1`, and `S = ({1,…,n−1} ∖ {q}) ∪ {g}` for a single `g`. Then
**`M(S) > 1/n`** (so `S` is not tight). Two cases:

**(i) `g` not a multiple of `q`.** Since `q > (n−1)/2`, `q` is the UNIQUE multiple of `q` in `{1,…,n−1}`
(`2q > n−1`); with `q` removed and `g` a non-multiple, `S` has **no multiple of `q`**. By the q-witness
(THM-523) at `t = 1/q`: every `v∈S` has `v ≢ 0 (mod q)`, so `‖v/q‖ ≥ 1/q`, giving `M(S) ≥ 1/q > 1/n`
(as `q < n`). ∎

**(ii) `g = mq`, `m ≥ 2`** (the only way to cover `q`). Let `p = n+q`. Since `q` is prime, `q > n/2`, and
`q ≠ n`, `gcd(q, p) = gcd(q, n) = 1`, so `a := q⁻¹ (mod p)` exists. Take `t = a/p`. Compute `v·a (mod p)`:
- `v ∈ {1,…,n−1}∖{q}`: `va ∈ {0, 1, p−1}` ⟺ `v ≡ 0, q, −q≡n (mod p)`; but `v≠0`, `v≠q` (removed),
  `v≠n` (`n>n−1`). So `va ∉ {0,1,p−1}`.
- `v = mq`: `va = mq·q⁻¹ ≡ m (mod p)`; for `2 ≤ m ≤ p−2` (i.e. `g = mq < 2n`), `m ∉ {0,1,p−1}`.
So **no runner has `va ∈ {−1,0,1} (mod p)`**. By the radius-1 k-witness (HYP-3741): `M(S) ≥ 2/p = 2/(n+q)`,
and `2/(n+q) > 2/(2n) = 1/n` since `q < n`. ∎

*Verified exactly:* `M(S) = 2/(n+q)` at `t = q⁻¹/(n+q)` for n=14 (q=11,13) and n=18 (q=11,13,17), m=2,3;
all runners avoid `{0,1,n+q−1}` (0 exceptions).

## The witness is fresh: the missing pair {q, n} mod (n+q)
Removing `q` and adding `mq` leaves the **±pair `{q, n}` uncovered mod `n+q`** (`q ≡ −n`, and neither `q` nor
`n` is a runner). This is a NEW modulus — not THM-401's `2n−1` pair-sum modulus, but a `q`-specific
`n+q`. It says: *the large prime `q` and its complement `n` form a resonant pair that only `q` itself (or `n`)
can break; substituting a multiple `mq` shifts the runner to residue `m`, leaving the pair open.* This is the
mechanism the repo's HYP-3749 ("punctured-core wide hole") sees measure-theoretically; here it is an exact
finite witness for the single swap.

## The residual = the crux (honest)
The `(n+q)`-witness needs **`n ∉ S`** — if `n ∈ S` (as a large runner), then `n·a = n·q⁻¹ ≡ −1 (mod p)`, so
runner `n` sits at residue `p−1` and **covers the witness**. So the FULL "large primes forced" (general tight
`S` missing `q`, possibly with `n ∈ S` and multi-element substitutes) is NOT closed: it reduces to ruling out
this residual `n∈S` case, which is precisely the **CRT-linkage** of HYP-3749/3745 (filling one modulus shifts
the deficit to another). Single-swap: closed. General: the repo's open crux.

## Number-theory landscape (58-thread catalog, integrated)
The repo's number theory clusters into five anchors, all bearing on the tight-skeleton:
1. **Continued fractions / Farey / Ostrowski** (THM-523 q-witness; HYP-3738/3739 covering-min = CF
   `[0;n−1,k]` semiconvergent; HYP-3732 Stern-Brocot; three-gap HYP-3717) — the combinatorial covering
   backbone.
2. **Cyclotomic / p-adic** (Φ₆=(n−1)n+1; THM-580 2-adic parity descent; THM-590 apex-7 cyclotomic gap;
   CRT escape HYP-3745/3751) — the modulus arithmetic behind covering & the residual.
3. **Additive / Fourier** (HYP-+2873 additive-energy = Fejér, refuted as a tightness lens; HYP-2974
   Toeplitz-PSD certificates; THM-515 singular series = theta form; THM-503/504 convergence) — the analytic
   route that localizes to covering.
4. **Primes / Jacobsthal** (HYP-2893 doubling = Jacobsthal-gated; THM-590 apex-7; the doubling operad =
   covering-preserving `k→2k`) — the sporadic gating and THIS large-prime thread.
5. **Burnside / Heegner / elliptic** (A000568 tournament classes; HYP-3730 class-number-1 at n=14=2·7;
   X₀(14) cusps) — the deeper arithmetic of the `n=14` arena.

**Niche ideas flagged for the general "large prime forced" residual** (the `n∈S` / CRT-linkage case):
- **Ostrowski k-witness beyond three-gap** — express the residual deficit as a higher Ostrowski digit (ties
  lowness to covering; the repo's HYP-3741/3739 direction).
- **Minimal-transversal / matroid parity** — is `{1,…,n−1}` the unique minimal transversal of the covering
  hyperplanes `{v ≡ 0 mod q}`? A matroid-rank argument could force the skeleton without case analysis.
- **Baker linear-forms-in-logs (p-adic)** — bound the Ostrowski height of a `q`-missing tight set below the
  reachable range, ruling out the residual.
- **Ihara/Ramanujan spectral gap** of the resonance graph — a spectral lower bound could replace the
  witness-by-witness CRT chase.
- **Pentagonal / Φ₆ partition support** — Φ₆ sits in Euler's pentagonal identity; the covering lattice's
  degeneracies may forbid the residual mediants.

## Status
- **PROVED (opus, elementary):** single-swap large-prime-forced — `q` prime in `((n−1)/2, n−1]` cannot be the
  lone swapped-out element of a tight set (`g` non-multiple ⇒ q-witness `≥1/q`; `g=mq<2n` ⇒ `(n+q)`-witness
  `≥2/(n+q)`; both `>1/n`). Verified n=14,18.
- **New object:** the `(n+q)`-witness from the missing pair `{q,n}` — an exact finite instance of HYP-3749's
  wide hole.
- **Residual (crux):** general `S` with `n∈S` — the CRT-linkage of HYP-3749/3745, not closed.
- **Catalog:** 58 number-theory threads across 5 anchors; niche leads flagged for the residual.

Related: BRAVE-lowness-attack-… (the structure), the-tight-skeleton-is-COMBINATORIAL-COVERING (the covering
frame), PATCH-TUNING-… (bounded g), THM-523 (q-witness), HYP-3741 (k-witness, used here), HYP-3749/3745
(the CRT residual), HYP-2893 (Jacobsthal), THM-401 (the OTHER pair modulus 2n−1); OPEN-Q-108.
