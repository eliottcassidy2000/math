---
id: THM-1730
title: "opus's TORAL k-NOMIAL CLOSURE TRANSFERS TO THE RADIAL INTERIOR — the |charge|≥2 residual of the cross-shell descent closes on the radial side exactly where it closes on the toral side, extending THM-1700's {−1,0,1} witness to genuine interior patterns. Verified UNCONDITIONALLY by direct radial elimination on {−2,−1,1,2} (constant AND s-dependent radial coefficients — the real Gamma-bridge test) and on the 5-charge {−2,−1,0,1,2} (where the charge-0 coefficient is forced nilpotent, i.e. DvdK-forbidden). Assembles, via the charge-radius lock (THM-1700) and the Gamma bridge (klein-S351), into: GMC(2) for P with ≤5 distinct charges — conditional on the bridge, whose charge-0 layer is already unconditionally closed (mac-mini-S147 Cauchy). MINED: opus's 'CT(m)=0 are vanishing sums' = THM-415 (prime modulus ⟹ no nontrivial vanishing) = the LRC cyclotomic vanishing-sum structure (death-star-S67 GMC↔LRC reflection); the three name one wall."
status: >
  UNCONDITIONAL for the specific bounded interior families (exact Gröbner closure over ℚ, no
  bridge, no asymptotics):
    - {−2,−1,1,2} constant coefficients (4 params) — variety ⊆ one-sided locus;
    - {−2,−1,1,2} with the +1 shell radial-linear (5 params) — the genuine s-dependent test,
      factorials mix (the 4·a11·c1 term at m=2), still closes;
    - {−2,−1,0,1,2} constant (5 params) — d0 (charge 0) nilpotent mod I AND every pos·neg
      nilpotent, so the only nullcone members are one-sided.
  CONTROL PASSES: one-sided P (charges +1,+2) gives E[P^m]=0 ∀m and E[QP^m]=0 for m>1 (Q=Z̄),
  i.e. MZ-harmless, correctly NOT flagged.
  ASSEMBLY (opus toral k-nomial + lock + Gamma bridge ⟹ GMC(2) for ≤5 charges): CONDITIONAL on
  the Gamma-bridge rigor. The direct eliminations above do NOT use the bridge — they are
  independent confirmations that the transfer holds on these patterns. GMC(2) at ALL radial
  degrees for ≤5 charges is what the bridge would add.
  GMC(2) REMAINS OPEN in general (unbounded charge count; the bridge's non-charge-0 rigor).
source: klein-2026-07-20-S375 (owner: pull recent agent work, take the most cutting-edge ideas as far as you can)
depends_on:
  - THM-1700  # the charge-radius lock + bottom-up descent + the {−1,0,1} witness
  - THM-1685  # opus: toral k-nomial TNC = Nullstellensatz emptiness, ≤5 charges closed
related:
  - THM-1530  # extreme-weight-±1 Lagrange (the {−1,0,1} / dual descent)
  - THM-1645  # the polar bridge: angular = DvdK, gap = radial
  - THM-415   # prime-modulus vanishing sums — opus's vanishing-sums bridge
  - "mac-mini-S147: complex-radial / charge-0 closed via Cauchy (the bridge's closed layer)"
  - "death-star-S67: GMC(2)↔LRC(14) 'positivity past the cancellation wall' reflection"
script: 04-computation/interior_residual_radial_klein_S375.py (+ .out)
---

# THM-1730 — the toral k-nomial closure transfers to the radial interior

## The setup, and what was open

The cross-shell descent (THM-1700) closes charge span `{−1,0,1}` — extreme weight `±1`. The
**interior residual** is spans with `|charge| ≥ 2` at both ends (extreme weight `≥ 2`), which
THM-1530(C), mac-mini's top-edge charge descent, and opus's k-nomial TNC all name from
different sides. opus THM-1685 closed the **toral** (constant-coefficient, leading-symbol) side
for `≤ 5` charges via a Nullstellensatz emptiness test. This file asks whether that closure
survives on the **radial** (`s`-integrated) side — the thing GMC(2) actually needs.

## What is verified — directly, without the bridge

By the charge-radius lock (THM-1700), a genuine `P` gives an integer-moment radial functional,
so `E[P^m]` are exact rationals in the coefficients and Gröbner elimination decides the
nullcone. The one-sided (nullcone) loci are where all charges share a strict sign.

| pattern | coefficients | params | result |
|---|---|---|---|
| `{−2,−1,1,2}` | constant | 4 | variety ⊆ one-sided; **closes** |
| `{−2,−1,1,2}` | `+1` shell radial-linear | 5 | s-dependent; factorials mix; **closes** |
| `{−2,−1,0,1,2}` | constant | 5 | `d0` nilpotent + no straddle; **closes** |

The moments show the **bottom-up descent** explicitly:
`E[P²] = 2a₁c₁ + 4a₂c₂` (both straddling pairs, the `4` from `E[Z²Z̄²]=2!`),
`E[P³] = 6a₁²c₂ + 6a₂c₁²`, and elimination drives every `pos·neg` cross-product into the ideal.

**Control:** a one-sided `P = a₁Z + a₂Z²` gives `E[P^m] = 0 ∀m` and `E[Z̄·P^m] = 0` for `m > 1`
— MZ-harmless, correctly not flagged. So the test has power and does not fire on the nullcone.

The `{−2,0,2}` "gcd escape" is confirmed separately: it is `{−1,0,1}` in `w = u²`, already
closed — the genuine interior is the gcd-1 patterns above.

## The assembly (conditional on the bridge)

```text
opus THM-1685   toral k-nomial TNC, ≤5 charges (Nullstellensatz)
   +  THM-1700  charge-radius lock ⟹ integer-moment radial functional
   +  S351      Gamma bridge: toral nonvanishing ⟹ radial nonvanishing
   ⟹  GMC(2) for P with ≤ 5 distinct charges.
```

The **number-of-charges** complexity parameter (opus's key insight — not the bidegree)
transfers from toral to radial through the lock. The bridge is what upgrades the bounded-degree
eliminations above to **all** radial degrees; its charge-0 layer is unconditionally closed
(mac-mini-S147, Cauchy transform), and the mechanism is proved (S351), so the conditional is
narrow. The direct eliminations here are independent, unconditional confirmations on the
interior — they are what THM-1700 was missing ("the general HYP-8470 ... is NOT closed here").

## Mined connection — one wall, three names

opus's observation that `CT(m) = 0` are **vanishing sums of coefficient monomials over the
charge lattice** is exactly **THM-415** (prime modulus ⟹ no nontrivial vanishing sum; composite
⟹ collisions), which is in turn the **LRC cyclotomic vanishing-sum** structure that
death-star-S67 identified in the GMC(2)↔LRC(14) reflection ("positivity past the cancellation
wall"). The toral CT-vanishing ideal, the THM-415 modulus dichotomy, and the LRC covering
vanishing sum are three faces of the same object — the minimal charge relation `m₀` is the
primitive vanishing sum, and `CT(2m₀)` adds the independent square-plus-correction that closes
the unit ideal. This is why opus's "≤5 CT levels suffice" and the LRC "B5 quintic certificate"
rhyme: both are the second rung past a primitive relation generating the whole obstruction.

## Scope

Unconditional on the three bounded interior families (exact closure). The general
`≤5`-charge GMC(2) is conditional on the Gamma bridge. GMC(2) at unbounded charge count remains
open. What this adds: the interior `|charge| ≥ 2` residual of the descent is confirmed to close
on the radial side, the `s`-dependent case included, and the toral/radial transfer is pinned to
opus's number-of-charges parameter.

*Files: `04-computation/interior_residual_radial_klein_S375.py` (+ `.out`).*
