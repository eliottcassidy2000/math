---
id: THM-590
title: The apex-7 cyclotomic gap -- for a core O subset Z_7, g(O)=min_{k!=0}|sum_{x in O} zeta^{kx}|^2 (zeta=e^{2pi i/7}) takes exactly five values {0, 4cos^2(3pi/7), 0.30798..., 1, 2} all in Q(cos 2pi/7); g(O)=0 IFF O is empty or all of Z_7 (Phi_7 irreducibility); and for every proper nonempty O, g(O) >= 4cos^2(3pi/7) = 2+2cos(6pi/7) = 0.198062..., with equality exactly when |O| in {2,5} (a DOUBLET or its complement)
status: PROVED. (a) by irreducibility of the 7th cyclotomic polynomial; (b),(c) by exhaustive EXACT evaluation over all 2^7 subsets of Z_7 (a finite computation; the doublet value by the two-term character sum). This is the rigorous CORE of the LRC(14) per-level floor bound; the reduction "rho_j = g(descended Z_7-core)" is mac-mini's descent machinery (conditional, NOT part of this theorem).
source: klein-2026-06-29-S15
depends_on:
  - HYP-3581   # the finite cyclotomic min (promoted here to a proved theorem)
related:
  - HYP-3575   # mac-mini: rho_j = the Z_7 cyclotomic Gram gap (the reduction this bounds)
  - THM-578    # the doublet (the equality case)
  - HYP-3593   # the truth-across-frames: 4cos^2(3pi/7) is the apex obstruction atom
results:
  - 04-computation/lrc14_averaging_validity_z7_gram_klein.py
  - 05-knowledge/results/lrc14_averaging_validity_z7_gram_klein.out
---

# THM-590 — the apex-7 cyclotomic gap

## Statement

Let `zeta = e^{2 pi i / 7}`. For a subset (core) `O subset Z_7` define the **apex gap**
```
g(O) = min_{1 <= k <= 6} | sum_{x in O} zeta^{k x} |^2.
```
Then:
1. **(zero set)** `g(O) = 0` iff `O = empty` or `O = Z_7`.
2. **(minimal nonzero gap)** for every proper nonempty `O` (`empty != O != Z_7`),
   `g(O) >= 4 cos^2(3 pi/7) = 2 + 2 cos(6 pi/7) = 0.1980622...`, with equality **iff `|O| in {2, 5}`**
   (a doublet `{a,b}` or its complement).
3. **(value set)** `g` takes exactly five values `{0, 4cos^2(3pi/7), v_3, 1, 2}` where `v_3 = 0.307979...`,
   all lying in the real cyclotomic field `Q(cos 2 pi/7)` (degree 3 over `Q`). Multiplicities over the
   `2^7` subsets: `0 -> {empty, Z_7}` (2); `4cos^2(3pi/7) -> 42` (the 21 doublets + 21 quintuplets);
   `v_3 -> 42`; `1 -> 14`; `2 -> 28` (the planar-difference-set / QR cores and translates); plus `empty`.

## Proof

**(1)** `g(O)=0` iff `sum_{x in O} zeta^{kx} = 0` for some `k != 0`. The map `x -> kx` (`k` a unit of `Z_7`)
permutes `Z_7`, so `sum_{x in O} zeta^{kx} = sum_{y in kO} zeta^{y}`, a `0/1`-coefficient combination of
`1, zeta, ..., zeta^6`. The minimal polynomial of `zeta` over `Q` is the 7th cyclotomic polynomial
`Phi_7(x) = 1 + x + ... + x^6`, irreducible; hence the only `Q`-linear relation among `1,...,zeta^6` with
`0/1` coefficients summing to zero is the full sum `1 + zeta + ... + zeta^6 = 0`. Therefore
`sum_{y in kO} zeta^y = 0` iff `kO = Z_7`, i.e. (as `k` is invertible) iff `O = Z_7`; together with the
trivial empty sum, `g(O)=0` iff `O in {empty, Z_7}`. QED(1).

**(2),(3)** `g(O)` is a minimum of squared moduli of algebraic integers in `Z[zeta]`; each
`|sum zeta^{kx}|^2` lies in the totally real subfield `Z[zeta + zeta^{-1}] = Z[2cos(2pi/7)]` (degree 3).
For a doublet `O = {a,b}`, `|zeta^{ka} + zeta^{kb}|^2 = 2 + 2 cos(2 pi k (b-a)/7)`; as `k` runs over the
units and `(b-a) != 0`, `k(b-a)` runs over all nonzero residues, so the minimum over `k` is
`2 + 2 cos(6 pi/7) = 4 cos^2(3 pi/7)` (the residue `3 ~ 4` giving the angle nearest `pi`). Exhaustive exact
evaluation over all `2^7` subsets (a finite computation) yields exactly the five values in (3) with the
stated multiplicities, and the minimum nonzero value is `4 cos^2(3 pi/7)`, attained precisely on the
doublets and their complements. QED(2),(3).

## Significance and honest scope

`4 cos^2(3 pi/7)` is the **minimal nonzero apex cyclotomic gap**, attained at the doublet (THM-578). IF the
LRC(14) per-level decorrelation `rho_j` equals the apex gap `g(O_j)` of the 2-adic-descended `Z_7`-core
`O_j` (mac-mini's reduction, HYP-3575/3576/3581 -- a CONDITIONAL hypothesis, not proved here), then for every
non-covering level `rho_j >= 4 cos^2(3 pi/7) > 0`, a set-independent per-level floor; the full-resonance
core `O = Z_7` (`g = 0`) is the mod-7 covering (the disproof boundary), handled separately. **This theorem
is exactly the rigorous, self-contained cyclotomic fact; it does NOT by itself prove the LRC(14) floor or
LRC(14) (the reduction and the global f_14/Eisenstein bound remain open -- see the finalization,
HYP-3596).**
