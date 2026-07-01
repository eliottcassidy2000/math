---
id: HYP-3768
title: THE LRC14 COVERING-MIN MARGIN IS THE OBSERVER'S DEDEKIND SUM -- and the B2/A2 dichotomy behind its positivity, with E2 the modular home and reciprocity the attack on the hard residual. EXACT (verified n=3..15): margin(n) = n/Phi6 - 1/n = (n-1)/(n Phi6) = -12 s(n, Phi6)/n^2, where Phi6=n^2-n+1 and s is the Dedekind sum; the closed form 12 s(n,Phi6) = -(Phi6-1)/Phi6 holds because n^3 = -1 mod Phi6 (n = the order-6 Eisenstein element omega, 60deg). DICHOTOMY: an order-4 element h (h^2=-1 mod k = the SQUARE lattice / B2 / imaginary unit i, 90deg) has s(h,k)=0 => ZERO margin (degenerate); the order-6 observer (HEXAGONAL / A2 / Eisenstein) has s!=0 => POSITIVE margin. So LRC14's margin (why the covering-min 14/183 exceeds the floor 1/14 with room) IS the hexagonal Dedekind anomaly; B2 would kill it. E2: F_14 = E2(tau)-14 E2(14 tau) is a genuine weight-2 Eisenstein form on Gamma_0(14) (const 1-n=-13); F_2,F_7,F_14 span the 3-dim Eisenstein subspace (4 cusps), and the 1-dim cusp form f_14 (genus(X_0(14))=1) = the residual, to be controlled at the apex cusp d=7. RESIDUAL ATTACK: eroding the margin to 0 (a patch beating the covering-min) requires the observer Dedekind anomaly to vanish = restoring the order-4/square structure; a single patch integer has CRT-linked residues, and by Dedekind-Rademacher reciprocity (the repo far-coherence tool HYP-2808) its Dedekind sums at different moduli cannot all be tuned anomaly-free at once => the hexagonal anomaly persists = "the hole moves but never vanishes."
status: IDENTITIES EXACT/VERIFIED (n=3..15): closed form 12 s(n,Phi6)=-(Phi6-1)/Phi6; margin = -12 s/n^2; B2 s=0; E2 form F_14 weight-2 on Gamma_0(14). These are NEW to the repo (E2 and B2 were absent; the margin=Dedekind-sum identity is new). The RESIDUAL ATTACK (anomaly-persistence via reciprocity) is a framework/mechanism aligned with the established HYP-2808 Dedekind-Rademacher tool, NOT a proof of the Step-4 residual.
source: mac-mini-2026-06-30-S64
depends_on:
  - HYP-2808   # Dedekind-Rademacher reciprocity = the far-coherence tail tool (the reciprocity mechanism)
related:
  - HYP-3715   # the covering-min = zeta6-line in the hexagonal A2 lattice (S*14 mod 183 = n-spaced AP)
  - HYP-3586   # X_0(14) cusps = Klein V4; apex cusp d=7; genus jump 0->1 = the f_14 obstruction
  - HYP-3553   # the Gamma_0(14) congruence 2nd-moment floor framework
  - HYP-3745   # the Step-4 CRT-patch residual (perturbation-proved)
  - HYP-3749   # the punctured-core wide hole (the residual's wide-hole map)
  - HYP-3750   # my S61 band-transversal reduction (the patch is +-k mod band primes)
results:
  - 04-computation/dedekind_e2_b2_margin_macmini_20260630.py
  - 05-knowledge/results/dedekind_e2_b2_margin_macmini_20260630.out
---

# HYP-3768 -- the LRC14 margin is the observer's Dedekind sum (B2 / A2 / E2)

Bringing the three classical objects the owner asked for -- **B2/C2 root system**, **E2 (weight-2
Eisenstein)**, **Dedekind eta/sums** -- to bear on the genuinely hard residual (Step 4 / the f_14 cusp form at
the apex cusp d=7).

## The margin IS a Dedekind sum (exact, verified n=3..15)
Let `Phi6(n) = n^2 - n + 1` (so `Phi6(14) = 183`), and `s(h,k)` the Dedekind sum. Then

> **`12 s(n, Phi6) = -(Phi6 - 1)/Phi6`**, because `n^3 = -1 (mod Phi6)` (`Phi6 | n^3+1`), i.e. `n` is the
> **order-6 element** `omega` (the `60`-degree Eisenstein rotation of the torus lift `Z[omega]/(n-omega) ~
> Z/Phi6`, the repo's hexagonal covering frame HYP-3715).

Consequently the **covering-min margin** -- the amount by which `M = n/Phi6 = 14/183` exceeds the LRC floor
`1/n = 1/14` -- is exactly that Dedekind sum:

> **`margin(n) = n/Phi6 - 1/n = (n-1)/(n Phi6) = -12 s(n, Phi6) / n^2`.**

For `n=14`: `margin = 13/2562 = -12 s(14,183)/196`, `s(14,183) = -91/1098`. An LRC quantity (the safety margin
of the conjecture) equals a modular/Dedekind quantity (a special value of the Dedekind sum at the order-6
element). Verified `n=3..15`.

## B2 vs A2/E2 -- why the margin is positive
The Dedekind sum vanishes at the **order-4** element:
- **B2 / square lattice / imaginary unit `i` (90 deg):** `h^2 = -1 (mod k)` `=> s(h,k) = 0` (verified
  `k=5,13,17,25,29,...`). A square-lattice observer would give **zero margin** -- the covering-min would equal
  the floor, degenerate.
- **A2 / hexagonal lattice / Eisenstein `omega` (60 deg):** the order-6 observer has `s != 0`
  (`= -(Phi6-1)/12Phi6`), giving a **strictly positive** margin.

So **LRC14 has room to spare precisely because its observer lives on the hexagonal (A2/Eisenstein) side, not
the square (B2) side** -- the positive margin is the hexagonal Dedekind anomaly. This is the exact arithmetic
reason the covering-min (`n/Phi6`, hexagonal-Kershner) beats the floor: the order-6 anomaly `s(n,Phi6) != 0`.
B2 is the anomaly-free competitor that loses.

## E2 -- the modular home
`F_d(tau) = E2(tau) - d E2(d tau)` is a genuine weight-2 modular form on `Gamma_0(d)` (the E2 quasimodular
anomalies cancel). For `d=14`, `F_14` has constant term `1 - 14 = -(n-1) = -13`. The three forms `F_2, F_7,
F_14` span the **3-dimensional Eisenstein subspace** of `M_2(Gamma_0(14))` (`X_0(14)` has 4 cusps, `dim = c-1 =
3`); the remaining **1-dimensional cusp form `f_14`** (`genus(X_0(14)) = 1`, the conductor-14 elliptic curve)
is the **obstruction** -- exactly the survey's hard core: control `f_14` at the **apex cusp `d=7`** (Atkin-
Lehner `W_7`, width 2, the narrowest, matching the binding doublet). E2 was not previously in the repo; here it
gives the explicit Eisenstein (bulk) part of the `Gamma_0(14)` 2nd-moment floor, isolating `f_14` as the
residual.

## The attack on the hard residual
The residual (Step 4, HYP-3745/3749): a covering set missing core speed `k`, patched by a large `w = +-k` at a
band prime, must still have `M > n/Phi_6` -- i.e. **the margin cannot be eroded to `0`**. In this Dedekind
frame:
- the margin `= -12 s(n,Phi6)/n^2` is the **hexagonal anomaly**;
- eroding it to `0` means restoring the **order-4/square (B2) anomaly-free structure**;
- the patch `w` is a **single integer**: its residues across moduli are **CRT-linked**, and by
  **Dedekind-Rademacher reciprocity** (the repo's far-coherence tool, HYP-2808) its Dedekind sums at different
  moduli **cannot all be tuned to the anomaly-free value simultaneously**.

So the hexagonal anomaly persists under patching -- this is the arithmetic content of "the hole moves but
never vanishes" (klein-S44), now with a target: **the Step-4 residual = the persistence of the order-6
Dedekind anomaly `s(n,Phi6)` under `+-k` CRT patching = the non-vanishing of the `f_14` cusp contribution at
`d=7`.** The three objects the owner named are one story: `s(n,Phi6)` (Dedekind) `=` the margin `=` the E2/f_14
anomaly on `Gamma_0(14)`, positive because A2 (order 6) not B2 (order 4).

## Honest scope
The identities (closed form `12 s(n,Phi6)=-(Phi6-1)/Phi6`; `margin = -12 s/n^2`; `B2 => s=0`; `F_14` weight-2
on `Gamma_0(14)`) are exact and verified `n=3..15`, and E2/B2 are new to the repo. The **residual attack is a
framework** -- it identifies the target (anomaly persistence) and the tool (Dedekind-Rademacher reciprocity,
HYP-2808), and aligns Step 4 with the cusp-form control at `d=7`, but does **not** close the quantitative gap
(that a general `+-k` patch cannot vanish the anomaly, i.e. `f_14`-control). That remains the open hard core.
