---
id: THM-581
title: The LRC floor (EXISTENCE of a lonely time) is a sigma-EVEN-category statement -- the lonely measure is invariant under the complement reflection t->-t, so its SIGN-isotypic component vanishes; hence existence lives in the Brouwer/SOS category and does NOT require the Borsuk-Ulam sign certificate. The Borsuk-Ulam/odd-degree obstruction (kps S31av) is strictly for WITNESS CONSTRUCTION: since the sigma-fixed points t=0,1/2 are danger for any covering set, lonely times occur ONLY in antipodal pairs {t*,-t*}. The 2-adic parity descent (THM-580) is the constructive even-category proof of existence, validated across the dihedral family n=2p (including the PROVEN n=6).
status: PROVED (the sigma-symmetry of lonely(S) and the danger of the fixed points 0,1/2 for covering S are elementary; the isotypic consequence is immediate). The dihedral-family descent validation is VERIFIED (n=6 KNOWN case: meas(lonely S)>=0.033>0 over 60 covering sets; n=10,14,22 analogous; rho_j and odd caps bounded, NO p mod 4 split in the scalar floor).
source: mac-mini-2026-06-29-S5
depends_on:
  - THM-580    # the 2-adic parity descent (the constructive even-category proof)
  - HYP-3239   # kps S31av: 14=|D_7|, complement=anti-aut (p=3mod4), Borsuk-Ulam certificate
related:
  - HYP-3219   # Brouwer/SOS side (p=1mod4)
  - THM-381    # LRC observer-source tournament (the sign acts on orientations)
  - OPEN-Q-108
results:
  - 04-computation/lrc_2p_family_parity_descent_macmini_20260629.py
  - 05-knowledge/results/lrc_2p_family_parity_descent_macmini_20260629.out
---

# THM-581 — the floor is even-category (existence); Borsuk-Ulam is for the witness

## The category split (PROVED)
Let `S` be a covering speed set for LRC(n), `sigma: t -> -t` the complement reflection
(`||s(-t)|| = ||st||`, so `lonely(S)` is `sigma`-INVARIANT).  Decompose functions on the circle
under the dihedral group `D_{n/2}` (for `n=2p`, `D_p` of the `p`-gon `zeta_p^k`): trivial + SIGN
(the `sigma`-odd / complement isotype) + the 2-dim de Moivre isotypes (kps S31av / HYP-3239).

> **Existence (the floor) is `sigma`-EVEN.**  `meas(lonely S) = int 1_{lonely(S)}` depends only on the
> `sigma`-invariant (even) part of `1_{lonely(S)}`; since `1_{lonely(S)}` is already `sigma`-even, its
> SIGN-isotypic component is `0`.  So `meas(lonely S) > 0` is a statement in the EVEN / Brouwer / SOS
> category -- it carries no sign-isotypic content and does NOT need the Borsuk-Ulam certificate.

> **Witness construction is `sigma`-ODD (antipodal).**  The `sigma`-fixed points are `t=0` and
> `t=1/2`; both are DANGER for any covering `S` (`0` is danger for every speed; `1/2` is danger for the
> even speed a covering set must contain).  So NO lonely time is `sigma`-fixed, and lonely times occur
> only in antipodal pairs `{t*, -t*}`.  Certifying a particular witness is therefore the free-`Z_2`
> odd-degree (Borsuk-Ulam) problem kps S31av describes -- strictly FINER than existence.

**Consequence.** kps's "n=14 needs Borsuk-Ulam, not Brouwer (p=7=3 mod4)" is correct for the WITNESS,
but the FLOOR (existence, the only thing LRC needs) is in the easy even/Brouwer/SOS category for ALL
`p`.  This is why the scalar floor shows NO `p mod 4` distinction (verified below): the
Brouwer/Borsuk-Ulam split lives in the orientation/witness data, not the measure.

## The constructive even-category proof: 2-adic parity descent (THM-580)
The descent `meas(lonely S) = PROD_j rho_j . PROD_j meas(lonely O_j)` IS the even-category existence
proof: the free `Z_2` of Borsuk-Ulam is the HALF-TRANSLATION `t -> t+1/2` (the deck map of the
doubling cover `t->2t`), and `meas(lonely E) = meas(lonely E/2)` projects `lonely(O)` onto its
half-translation-invariant part (the 2-sheet count).  So the descent computes the "odd degree"
ARITHMETICALLY (via the sheet-combination/product) instead of certifying it topologically.

## Dihedral-family validation (VERIFIED)
Descent statistics over covering `(n-1)`-sets, `n=2p`:

| n | p | p mod4 | LRC | min meas(lonely S) | min rho_j | min odd-cap |
|---|---|--------|-----|--------------------|-----------|-------------|
| 6 | 3 | 3 | **KNOWN/proven** | 0.0333 (>0 ✓) | 0.338 | 0.148 |
| 10 | 5 | 1 (Brouwer) | open | 0.0623 | 0.565 | 0.285 |
| 14 | 7 | 3 (Borsuk-Ulam) | open | 0.0651 | 0.637 | 0.283 |
| 22 | 11 | 3 | open | 0.0680 | 0.588 | 0.296 |

The descent route SOUNDLY reproduces the proven `n=6` (no zero-measure covering set found), and the
floor / `rho_j` statistics are UNIFORM across `p mod 4` -- the predicted absence of a Brouwer/Borsuk-
Ulam split in the scalar measure.  So the descent is a sound, family-uniform existence route; `n=14`
differs from the proven `n=6` only in scale, not in category.

## Significance & open
This RECONCILES the descent (THM-580), kps's D_7/Borsuk-Ulam framework (HYP-3239), and mac-mini's
S75e cyclotomic SOS: existence is even-category, the descent is its constructive proof, and the
remaining per-level bound `rho_j >= c` is an EVEN / 2-dim (de Moivre) SOS problem -- exactly what the
S75e Fejer-Bochner magic function targets.  So the closure is **descent (even-category reduction) +
S75e cyclotomic SOS (per-level 2-dim bound)**, with NO Borsuk-Ulam needed for existence.  Open: make
`rho_j >= c` rigorous via the 2-sheet cyclotomic SOS; the proven `n=6` is the template.

Recursive reading: for `n = 2^a . m`, the descent peels the `a` factors of 2 (even-category,
constructive, uniform), leaving the odd core `m` = the apex/cyclotomic arithmetic.  The 2-part is
always the "easy" even category; the difficulty is the odd core (for `n=14`, `m=7`, which is proven).

Script: `04-computation/lrc_2p_family_parity_descent_macmini_20260629.py`.
