---
id: THM-1740
title: "GMC(2) ON ANY BOUNDED (charge-count, degree) STRATUM IS A FINITE GROEBNER TEST -- and it CLOSES the minimal span-3 both-signs case of THM-1540's open residual. Assembling three of my results: THM-1535 (charge lattice), THM-1685 (angular nullcone = Nullstellensatz emptiness), THM-1540 (radial = Laplace/Gamma). At n=2 (one complex Gaussian), a polynomial P with charge set C (|C| <= K) and each charge-part a degree-<=d polynomial in s = |z|^2 has E[P^m] = E of the charge-0 part of P^m -- a finite sum over charge-representations of 0 (ANGULAR, THM-1685) whose s-dependence is E[s^k]=k! (RADIAL, THM-1540). Bounded (K,d) => finitely many representation-levels and s-powers => the nullcone condition {E[P^m]=0 for all m} is a FINITE polynomial ideal in the coefficients, and GMC(2) on the stratum <=> that ideal has no zero with the both-signs charge-parts nonzero <=> a Rabinowitsch/Groebner EMPTINESS test. VERIFIED on the minimal span-3 both-signs family P = c_{-1} zbar + c_0(a+b|z|^2) + c_{+1} z, charges {-1,0,1}: V(nullcone) cap (c_{-1}c_{+1} != 0) has 1 in the saturated ideal -- EMPTY -- so GMC(2) holds there, closing an instance of THM-1540's open both-signs residual. HONEST CAVEAT: the test is finite but not cheap (a degree-2-shell span-3 Groebner timed out at 10 min); decidability, not efficiency. FRAMING: cross-shell descent is the bottom-up sequence of these finite tests, klein's shell-mixing functional L being a finite polynomial coupling at each bound -- the SAME Nullstellensatz-emptiness framing"
status: PROVED (bounded (K,d) GMC(2) is a finite Groebner emptiness test -- decidability) + VERIFIED (minimal span-3 both-signs stratum empty). Cross-shell framing is a proposed unification. Unbounded GMC(2) remains open.
author: opus-2026-07-20-S427
depends_on: [THM-1535 (charge lattice), THM-1685 (k-nomial Nullstellensatz), THM-1540 (radial Laplace/charge-0), THM-1580 (product shadow), klein cross-shell descent (HYP-8430/8470)]
---

# THM-1740 — Bounded GMC(2) is a finite Gröbner test

The owner's synthesis: *unconditional GMC(2) on any bounded charge-count + degree is a finite
Gröbner test; the angular nullcone is a Nullstellensatz emptiness test; and the same framing
should close cross-shell descent.* This note assembles it and closes a bounded instance.

## 1. The two finite layers

`n = 2` = one complex Gaussian `z`, `E[z^a \bar z^b] = a!\,δ_{ab}`. Charge
`q(z^a\bar z^b) = a−b`. Write

```
P = Σ_{q ∈ C} z^{q_+} \bar z^{q_-} · P_q(s),     s = |z|^2 ,   q_± = max(±q, 0).
```

`E` annihilates every nonzero total charge, so `E[P^m] = E[\,(\text{charge-0 part of } P^m)\,]`.
That charge-0 part is a sum over **representations of `0` as an `m`-term sum of charges from
`C`** (the **angular / Nullstellensatz** layer, THM-1685), and the coefficient of each is a
polynomial in `s` evaluated by `E[s^k] = k!` (the **radial / Laplace–Gamma** layer,
THM-1540).

> **Bounded ⟹ finite.** If `|C| ≤ K` and every `P_q` has `\deg_s ≤ d`, then only finitely
> many representation levels and `s`-powers can occur before the ideal stabilises, so
> `{E[P^m] = 0 : m ≥ 1}` is a **finite polynomial ideal `I`** in the coefficients of `P`.

## 2. GMC(2) on the stratum is a Nullstellensatz emptiness test (PROVED)

By THM-1535, a GMC(2) violator is a nullcone element with **charges of both signs**. On a
bounded stratum:

> **Theorem.** GMC(2) holds on the `(K,d)` stratum **iff** `V(I)` has no point with the
> positive-charge and negative-charge parts both nonzero — i.e.
> `1 ∈ I + ⟨1 − w·(\text{sat of the both-signs coefficients})⟩` (Rabinowitsch). This is a
> **single finite Gröbner computation**, so GMC(2) is **unconditionally decidable** on every
> bounded stratum.

## 3. Verified: the minimal span-3 both-signs stratum is empty

`P = c_{−1}\bar z + c_0(a + b|z|^2) + c_{+1}z`, charges `{−1,0,1}` — the minimal genuinely
three-charge, both-signs family. Its nullcone conditions (through `m = 6`):

```
E[P]   = c_0(a+b)
E[P²]  = c_0²(a²+2ab+2b²) + 2 c_{−1}c_{+1}
E[P³]  = c_0[ c_0²(a³+3a²b+6ab²+6b³) + 6(a+2b)c_{−1}c_{+1} ]
…
```

Saturating by `c_{−1}c_{+1} ≠ 0` and computing the Gröbner basis returns `[1]`:

> **`V(\text{nullcone}) ∩ (c_{−1}c_{+1} ≠ 0) = ∅`.** No span-3 both-signs nullcone element
> exists at `n = 2`. GMC(2) holds on this stratum by a finite test — **closing an instance
> of THM-1540's open both-signs residual.**

Consistent with THM-1535 (`n = 2` charge lattice is rank 1, so both-signs cancellation cannot
be sustained), now confirmed constructively on a 3-charge family with a nontrivial charge-0
part.

## 4. Honest caveat: finite ≠ cheap

The test's cost grows fast. A **degree-2-shell** span-3 stratum (charge-parts quadratic in
`s`, 7 coefficient unknowns) did **not** finish a Gröbner elimination in 10 minutes. So this
is a **decidability** result — every bounded stratum is a terminating finite test — not an
efficient one. The `(k−1)`-level bound (THM-1705) and the finite-place mod-`p` certificate
(THM-1735) are the levers that keep the cost down: reducing mod a good prime and using only
`k−1` levels shrinks the Gröbner.

## 5. Cross-shell descent, same framing

klein's cross-shell descent runs **bottom-up** through the radial shells `ρ = |z|^2` (the
Hermite/Laguerre layers, THM-1660). At each shell the coupling that klein's functional `L`
mixes is a **finite polynomial system** in the shell coefficients — precisely a bounded
`(K,d)` instance of §2. So:

> **Cross-shell descent = the bottom-up sequence of finite Gröbner emptiness tests**, one per
> shell bound. The "cross-shell coupling" is the resultant/elimination linking consecutive
> shells; the descent terminates because each step is a finite Nullstellensatz test and the
> shells are ordered. This is the same Nullstellensatz-emptiness framing as the angular
> nullcone — the owner's proposed unification, now explicit.

**What this does and does not give.** It makes GMC(2) *decidable on every bounded stratum*
and *closes the strata one checks* (§3). It does **not** close unbounded GMC(2): that needs a
*uniform* bound on `(K,d)` or a *uniform* emptiness across all shells — which is exactly the
`(k−1)`-level bound (THM-1705) and the amoeba/carry uniform separation (THM-1735, HYP-8535)
promoted from per-stratum to all-strata. The angular half has those levers; the cross-shell
half needs the analogous uniform shell bound (klein's convergence lemma).

## 6. Next

1. **Uniform stratum bound.** Combine THM-1705's `(k−1)`-level cap with a degree bound `d ≤
   g(K)` so that a single `(K, g(K))` Gröbner test certifies GMC(2) for *all* `P` of charge-
   count `≤ K` — collapsing infinitely many strata to one.
2. **Cross-shell as resultant tower.** Make the shell-to-shell coupling an explicit resultant
   and prove the tower's emptiness propagates bottom-up (the finite-Gröbner analogue of
   klein's convergence lemma).
3. **Cost control.** Always reduce mod a good prime (THM-1735) and cap at `k−1` levels
   (THM-1705) before running Gröbner — the degree-2-shell timeout shows the naive test is
   impractical past the minimal stratum.

## Verification

`04-computation/gmc2_bounded_groebner_opus_S427.py` — the span-3 both-signs nullcone
conditions and the Rabinowitsch emptiness (`1` in the saturated ideal). Output in
`05-knowledge/results/`.
