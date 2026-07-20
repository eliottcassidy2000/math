---
id: THM-1670
title: "THE TORAL RECURRENCE ORDER IS EXACTLY D = M+N, SO THE ESTIMATE-FREE THREE-TERM DESCENT IS INTRINSICALLY THE (M,N)=(1,1) CASE ONLY — answering THM-1620's two named-next questions, one cleanly and one with an honest negative. (a) ORDER: the toral sequence a_m = CT(Λᵐ) = [u^{Mm}]Rᵐ (mac-mini THM-1610), R a degree-D polynomial with R(0),lc(R) ≠ 0, is P-recursive of order EXACTLY r(M,N) = D = M+N — confirmed decisively for all nine (M,N) with D = 2..6 by an ascending-order test with TWO primes agreeing and ≥ 20 holdout rows (after a first pass mis-reported (2,2) as order 3 from insufficient holdout, and (1,2) as order 2 from one degenerate instance — both corrected). The order depends on D ALONE, not on min(M,N): (1,3) and (2,2) both give 4. So ORDER 2 ⟺ D = 2 ⟺ (M,N) = (1,1). The Legendre/Hermite orthogonal-family structure of THM-1620, and the one-lemma three-term no-common-root descent it rides on, are UNIQUE to M = N = 1 — which is exactly why that bridge closes cleanly there and not obviously beyond. (b) HONEST NEGATIVE: I conjectured in THM-1620 that the descent extends because the trailing coefficient P_0(0) would be (a factor of) disc(R), reducing everything to a resultant. That is TRUE at (1,1) — P_0(0) = g₁²−4g₀g₂ = disc(R) exactly, by hand — but FALSE beyond: at (1,2) there is a double-root R (disc = 0), R = −3+4u+u²−2u³, with P_0(0) ≠ 0. So {P_0(0) = 0} is NOT the discriminant locus for D ≥ 3, and the 'apparent-singularity = discriminant' shortcut is WITHDRAWN. (c) STRUCTURE: F(t) = Σa_m tᵐ is a symmetric function of the M small roots of z^M = tR(z), verified as F = Σ_i R(z_i)/(M R(z_i) − z_i R'(z_i)); its algebraic degree is D, which is why the order is D and why a single orthogonal family cannot capture D ≥ 3. CONSEQUENCE for formalization: an elementary, Lean-checkable proof of TNC beyond (1,1) cannot be the ThreeTerm descent; it needs the genuine order-D recurrence with trailing-coefficient control, and (b) shows that control is not a clean discriminant condition. TNC itself is classically DvdK (mac-mini THM-1630/S142, residues+Liouville) — analytic, not obviously formalizable — so the elementary/formalizable route past (1,1) is genuinely open"
status: >
  (a) SOLID.  r(M,N) = D verified for all (M,N) with D = 2..6, two independent primes
  (2^61−1 and 2^62−57) agreeing, >= 20 holdout rows (a spurious recurrence would need 20
  coincidences, ~ p^-20).  The (1,2) generic order = 3 confirmed over 6 instances after a
  degenerate instance gave a spurious 2; (2,2) = 4 after a first pass wrongly gave 3.  The
  two corrections are logged as method notes, not swept.
  (b) HONEST NEGATIVE.  P_0(0) = disc(R) PROVED by hand at (1,1); the generalization is
  REFUTED at (1,2) by an explicit disc=0 / P_0(0)!=0 witness.  Caveat stated plainly: the
  recurrence coefficient-degree s varies with R (non-generic R can drop it), which makes
  the exact trailing-coefficient locus subtle; the negative rests on a single clean witness
  and the directional claim ("not a clean discriminant condition") is what is asserted, not
  a full description of the locus.
  (c) VERIFIED (residue formula matches the series at 4 (M,N)); the algebraic-degree reading
  is standard diagonal theory, cited not reproved.
  This theorem CORRECTS an over-optimistic "named next" of my own THM-1620.  It proves no
  new open problem and it does not advance GMC(2); it maps precisely where the estimate-free
  route stops.
source: kind-pasteur-2026-07-20-S128c122 (owner: pursue the two concrete questions from THM-1620)
depends_on:
  - THM-1620    # the Pochhammer bridge, whose named-next this answers
  - THM-1660    # the orthogonality closure (order-2 descent)
related: [THM-1610, THM-1625, THM-1630, THM-1640, THM-1645]
script: 04-computation/toral_recurrence_order_and_discriminant_kps_S128c122.py, toral_order_resolved_kps_S128c122.py, toral_discriminant_test_kps_S128c122.py (+ .out)
---

# THM-1670 — the toral recurrence order is D, and the descent is order-2-only

THM-1620 ended with two questions. Here are the answers: (a) is clean and (b) is a negative
that retracts an optimistic guess of mine.

## (a) The order is exactly `D = M + N`

The toral sequence `a_m = [u^{Mm}] R^m` (`R` of degree `D`, `R(0), lc(R) ≠ 0`) is a diagonal
of a rational function, hence P-recursive. Its minimal recurrence order is

> **`r(M,N) = D = M + N`**, for every tested `(M,N)` with `D = 2..6`.

| `D` | `(M,N)` | order `r` | coeff-deg `s` |
|---|---|---|---|
| 2 | (1,1) | **2** | 1 |
| 3 | (1,2) | 3 | 2 |
| 4 | (1,3), (2,2) | 4 | 6, 5 |
| 5 | (1,4), (2,3) | 5 | 10, 9 |
| 6 | (1,5), (2,4), (3,3) | 6 | 15, 14, 13 |

Confirmed by an ascending-order test: two primes (`2⁶¹−1`, `2⁶²−57`) must agree and there
must be `≥ 20` holdout equations. **Two corrections from a first pass are logged, not
hidden:** `(2,2)` first read as order 3 (insufficient holdout — a spurious low-order
relation), and `(1,2)` once read as order 2 (one degenerate random `R`); both are `D` on
careful recomputation. The order depends on **`D` alone**, not `min(M,N)`: `(1,3)` and
`(2,2)` both give 4.

So:

> **order 2 `⟺` `D = 2` `⟺` `(M,N) = (1,1)`.**

The Legendre (toral) and Hermite (radial) orthogonal families of THM-1620, and the
one-lemma three-term no-common-root descent they ride on, exist **only at `M = N = 1`**.
That is exactly why the Pochhammer bridge closes cleanly there — and why it does not
obviously extend.

## (b) The discriminant shortcut — withdrawn

THM-1620's "named next" guessed that the descent extends because the trailing coefficient
`P_0(0)` would be a factor of `disc(R)`, turning the whole question into a resultant
computation. At `(1,1)` this is exactly right:

> `P_0(0) = g₁² − 4g₀g₂ = disc(R)` (the `n`-recurrence
> `(n+2)a_{n+2} − (2n+3)g₁ a_{n+1} + (n+1)D a_n = 0`, by hand).

**Beyond `(1,1)` it is false.** A normalization-independent zero-test (the nullspace is
1-dimensional, so "is the trailing entry zero" is well defined) finds, at `(1,2)`, a
double-root polynomial

> `R = −3 + 4u + u² − 2u³` (so `disc R = 0`) with `P_0(0) ≠ 0`.

So `{P_0(0) = 0}` is **not** the discriminant locus for `D ≥ 3`, and the apparent-singularity
= discriminant reduction is withdrawn. (Honest caveat: the recurrence's coefficient-degree
`s` varies with `R`, which makes the exact trailing-coefficient locus delicate; the claim
asserted is the directional one — *not a clean discriminant condition* — resting on the one
explicit witness.)

## (c) Why the order is `D`, and what it means for formalization

`F(t) = Σ a_m t^m` is a symmetric function of the `M` **small roots** of `z^M = t R(z)`, and

> `F(t) = Σ_i R(z_i) / (M R(z_i) − z_i R'(z_i))`

(verified against the series at `(1,1), (1,2), (2,2), (2,3)`). Its algebraic degree over
`ℚ(t)` is `D`, so the recurrence order is `D` and **no single orthogonal family captures
`D ≥ 3`** — an order-`D` sequence is not `c^m` times one classical polynomial at a fixed
point.

The consequence is about *method*, and it is the real payload:

- The `ThreeTerm.no_common_root` Lean lemma (THM-1620 §V) is order-2, and order 2 is now
  known to be **exactly** the `(M,N) = (1,1)` case. So the formalized estimate-free descent
  is complete *and* intrinsically boundaried — it is not a stepping stone to higher `(M,N)`.
- A Lean-checkable proof of TNC at `D ≥ 3` would need the genuine order-`D` recurrence with
  trailing-coefficient control, and (b) shows that control is **not** a clean discriminant
  condition.
- TNC itself is classically Duistermaat–van der Kallen (mac-mini THM-1630 / S142, one page
  of residues + Liouville). That is analytic and not obviously formalizable, so the
  elementary/formalizable route to TNC past `(1,1)` is **genuinely open** — and this
  theorem says precisely why the obvious route (extend the 3-term descent) cannot be it.

## Fleet context (converging, same day)

- **klein-S363 / THM-1640** retracts the S351 domination (confirming THM-1585) and reframes
  the mechanism as **positivity, not domination** — and states the gap exactly: on the
  sign-indefinite `{−1,0,1}` span, integrand-positivity fails. My **Favard positivity**
  (THM-1620 §IV: Hankel of `j!` positive definite) is the positivity that survives there,
  because it is positivity of the *moment sequence*, not the integrand.
- **death-star-S63** extends my Hermite closure (THM-1660) to a linear charge-0 coefficient
  in closed form — the `Sheffer-with-curve` direction THM-1660 named.
- This note is the **toral-side** structural map: it says the clean order-2 closure is
  unique to `(M,N)=(1,1)`, so those radial-side extensions are the right place to push, not
  higher toral `(M,N)`.

> **Named next (DONE — see THM-1690).** The coefficient-degree is
> **`s(M,N) = C(M+N,2) − gcd(M,N) + 1`** (verified on 15 cells; the `min(M,N)` reading in an
> earlier draft of this line was wrong — refuted by `(4,6)`, `gcd=2 != min=4`, `s=44`). The
> `gcd` is the order of the cyclic symmetry behind opus THM-1625's symmetric descent. THM-1690
> also settles the tuned cancellation locus (saddle-value collision, `!= disc R`) and the
> finishing statement (toral half = DvdK; the sole remaining GMC(2) gap is radial Laplace
> determinacy, bypassed on the `(1,1)` stratum by the orthogonal route).
