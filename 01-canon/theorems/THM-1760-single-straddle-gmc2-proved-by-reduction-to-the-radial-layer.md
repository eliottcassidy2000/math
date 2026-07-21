---
id: THM-1760
title: "SINGLE-STRADDLE GMC(2) IS PROVED IN CLOSED FORM, by reduction to the already-closed radial Laplace layer — and the moment-count bound is SUBSUMED, not separately bounded. A single straddle is P = αZ^p + W·q(V), V = ZW = |Z|², q(V) = Σᵢ βᵢ V^{aᵢ} (charge +p on one term, charge −1 on r terms at radial degrees aᵢ). The multinomial + the identity Z^{jp}W^{jp} = V^{jp} give, for EVERY j, the exact tower identity E[P^{j·m₀}] = C(j·m₀, j)·αʲ·L(Qʲ) with m₀ = p+1, Q(V) := V^p·q(V)^p, L(V^k) = k!; and E[P^m] = 0 for m not a multiple of m₀. Hence on the nullcone with α ≠ 0, L(Qʲ) = 0 for all j ≥ 1, which by THM-1675/1695 (the radial Laplace nullcone, CLOSED for complex Q via the Cauchy transform) forces Q ≡ 0, so q ≡ 0, so every βᵢ = 0 — P is ONE-SIDED. No Gröbner, no per-pattern bound: the moment-count level r·m₀ (THM-1740) is exactly the radial certifying level for the degree-p(1+max aᵢ) polynomial Q, so HYP-8540's bound for a single straddle is not a separate conjecture but a consequence of the closed radial layer. MULTI-STRADDLE localises the same way — each straddle contributes its coefficient-product generator at its own return level (verified: aZ²+bW+cW³ gives E[P³]=6ab², E[P⁵]=7200a³c², forcing b then c bottom-up, THM-1700) — and the sole residual for full GMC(2) is the LOCALISATION LEMMA: the dominant straddle's radial factor L(Q_dom^j) is isolated at its return level (which THM-1745's max law says no other straddle reaches)."
status: >
  The tower identity E[P^{j·m₀}] = C(j·m₀,j)·αʲ·L(Qʲ) is PROVED for all j (multinomial
  expansion + Z^{jp}W^{jp} = V^{jp}; the balanced part of P^{j·m₀} is the single term
  C(j·m₀,j)(αZ^p)^j(Wq)^{jp}) and VERIFIED-EXACT j=1..3 for p=1,2 and r=1,2,3. Given it, the
  reduction is rigorous and SINGLE-STRADDLE GMC(2) IS PROVED, resting on THM-1675/1695 (the
  radial layer, itself proved). The MULTI-STRADDLE case is NOT closed: the per-straddle
  factorisation is verified on the two-straddle witness, but the localisation lemma (dominant
  straddle isolated at its level) is stated, not proved.
source: mac-mini-2026-07-20-S152 (owner: work the closed-form uniform proof)
depends_on:
  - THM-1675  # radial Laplace nullcone, real p (jump argument)
  - THM-1695  # radial Laplace nullcone, complex p (Cauchy transform) -- the closed layer
  - THM-1740  # M* = r*m0 (now subsumed: r*m0 = the radial certifying level for Q)
  - THM-1745  # multi-straddle combines by MAX -- gives the localisation
related:
  - HYP-8540  # the uniform moment bound -- subsumed for single straddle
  - HYP-8505  # TNC's twin bound
  - THM-1700  # bottom-up descent = the order the straddle generators fire
---

# THM-1760 — single-straddle GMC(2), proved in closed form

## The reduction

A **single straddle** is
`P = α Z^p + W·q(V)`, `V = ZW = |Z|²`, `q(V) = Σᵢ βᵢ V^{aᵢ}`
— charge `+p` on one term, charge `−1` on `r` terms (`W·V^{aᵢ} = Z^{aᵢ}W^{aᵢ+1}`, charge `−1`).

**The tower identity (proved, all `j`).** Balance in `P^m` forces `j` copies of `αZ^p` and
`k = jp` copies of `W·q(V)`, so `m = j(p+1) = j·m₀`, and the balanced part of `P^{j·m₀}` is the
single multinomial term `C(j·m₀, j)·(αZ^p)^j·(W·q(V))^{jp}`. Since `Z^{jp}W^{jp} = V^{jp}` and
`E[V^k] = k! = L(V^k)`:

> **`E[P^{j·m₀}] = C(j·m₀, j)·αʲ·L(Q(V)ʲ)`,  `Q(V) := V^p·q(V)^p`**, and `E[P^m] = 0` for
> `m ∤ m₀`.

Verified exactly `j = 1..3` for `p = 1,2`, `r = 1,2,3`; the derivation above is general.

**The reduction.** On the nullcone with `α ≠ 0` (the two-sided branch), `C(j·m₀,j)αʲ ≠ 0`, so
`L(Qʲ) = 0` for **every** `j ≥ 1`. By **THM-1675/1695** — the radial Laplace nullcone,
`L(Qʲ) = 0 ∀j ⟹ Q ≡ 0`, closed for complex `Q` via the Cauchy transform — we get `Q ≡ 0`.
Since `Q = V^p·q^p` and `V^p ≢ 0`, `q ≡ 0`, so every `βᵢ = 0`.

> **Single-straddle GMC(2) is proved: `E[P^m] = 0 ∀m` with `α ≠ 0` forces `P` one-sided.** ∎
> (No Gröbner, no per-pattern bound — the whole moment tower factors through one radial
> polynomial `Q`, and the closed radial layer does the rest.)

## The moment-count bound is subsumed

THM-1740's per-straddle level `r·m₀` is **not a separate thing to bound**: the identity shows
`E[P^{j·m₀}]` *is* `L(Qʲ)` up to a nonzero scalar, and the number of levels needed to force
`Q ≡ 0` is exactly the radial certifying level for `Q` (a degree-`p(1+max aᵢ)` polynomial),
supplied by THM-1675/1695. So HYP-8540's bound, for a single straddle, is a **consequence of
the closed radial layer**, not an independent conjecture. The `2, 3, 7`-flavoured moment counts
were the radial layer's certifying levels all along.

## Multi-straddle localises

For general `P = r₀(V) + Σ_{k>0} Z^k p_k(V) + Σ_{k>0} W^k q_k(V)`, each straddle contributes its
own coefficient-product generator at its own return level. On the witness `aZ² + bW + cW³`
(straddles `m₀ = 3` and `5`):

`E[P³] = 6ab²`, `E[P⁴] = 288a²bc`, `E[P⁵] = 7200a³c²`, `E[P⁶] = 360a²b⁴`.

With `a ≠ 0`: `E[P³] = 0 ⟹ b = 0`, then `E[P⁵] = 0 ⟹ c = 0` — bottom-up (THM-1700), and the
straddles fire at their own levels without interfering, exactly as THM-1745's max law predicts.

> **The sole residual for full GMC(2) is the localisation lemma:** the dominant straddle's
> radial factor `L(Q_dom^j)` is isolated at its return level (which no other straddle reaches,
> by the max law), so the general case reduces to the single-straddle reduction — straddle by
> straddle. Proving that isolation is what remains.

## Honest scope

- **Single-straddle GMC(2) is proved** (the tower identity is general; the reduction to
  THM-1675/1695 is rigorous). It rests on the radial layer, which is itself proved.
- **Multi-straddle is NOT closed.** The per-straddle factorisation is verified on one
  two-straddle witness; the localisation lemma is stated, not proved. That lemma is now the
  single gap between here and full GMC(2) (with span-2, THM-1600, and complex radial, THM-1695,
  already closed).
- The identity `E[P^{j·m₀}] = C(j·m₀,j)αʲL(Qʲ)` is proved by the balanced-term argument; the
  `j=1..3` runs are confirmation, not the proof.
- No claim about `TNC`'s HYP-8505 here beyond the structural twinning; the radial-reduction
  trick is specific to the Gaussian moment functional `L` and does not obviously transport to
  the constant-term functional.

*Artifacts:* `04-computation/gmc2_uniform_proof_reduction_macmini_S152.py` (+out).
*Credits:* THM-1675/1695 (the radial layer this reduces to), THM-1740/1745 (the straddle
structure), THM-1700 (the bottom-up firing order).
