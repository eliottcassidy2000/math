---
id: HYP-3576
title: WORKING THE DESCENDED CORE -- klein-S7b's test executed and ANSWERED: the 2-adic-descent odd cores O_j are GENERALLY NOT Z_7-cyclic invariant (their residues mod 7 are not Z_7^*-multiplier-orbit unions; only the level-0 all-residue cores are, 2/16 over four binding configs), so the bare cyclotomic SOS does NOT close rho_j>=c; the required mechanism is Gamma_0(14) CONGRUENCE-AVERAGING over Z_7^*, which FLATTENS the (set-dependent, non-flat) raw apex Gram to the cyclotomic/octonion-optimal form (set-independent, depending only on the core SIZE |O_j|) -- i.e. the floor must literally 'manufacture the transitive symmetry it lacks' by averaging, confirming the Gamma_0(14) route (HYP-3553) over a bare-Z_7 cyclotomic SOS and matching klein-S8's empirical set-independent inf R'=0.344 >= 1/(2 zeta(2))
status: COMPUTED + ANSWERED (the test). The cores are not Z_7-invariant (verified, 4 binding configs); the Gamma_0-averaging flattens the apex Gram (verified, {1,3,5}: raw [9,0.31,0.64,5.05,...] -> averaged [9,0.25,...]). The literal Han-Lee Gamma_0(14) 2nd-moment constant remains the PROGRAM (floor owners).
source: mac-mini-2026-06-29-S28
related:
  - HYP-3575  # mac-mini S27: rho_j>=c = Z_7 cyclotomic Gram gap (IF the core is Z_7-invariant -- this tests that IF)
  - HYP-3566  # klein-S7b: the transitivity reframe + the test posed here
  - HYP-3571  # klein-S8: the EMPIRICAL set-independent floor inf R'=0.344 (this gives the mechanism)
  - HYP-3553  # mac-mini Gamma_0(14): the congruence-averaging = the set-independent floor
  - THM-580   # the 2-adic descent (the cores)
  - HYP-3547  # the octonion/Fano (the flat/optimal target the averaging reaches)
results:
  - 04-computation/descended_core_Z7_invariance_test_macmini_20260629.py
  - 05-knowledge/results/descended_core_Z7_invariance_test_macmini_20260629.out
---

# HYP-3576 -- the descended core: not bare Z_7, but Gamma_0(14)-averaged

## The test (klein-S7b), executed
HYP-3575/3566 reduced the floor's `rho_j>=c` to: IS the descended core Z_7-cyclic invariant? If yes, the
bare Z_7 cyclotomic Gram (PSD/Bochner = SOS) closes it; if no, `Gamma_0(14)` supplies the symmetry. I
ran the 2-adic descent (THM-580: `S=O u E`, `S'=E/2`, recurse) on four binding configs (tightest
`{1..12,182}`, consec `{1..13}`, skip-12 `{1..11,13,84}`, even-heavy) and tested each odd core `O_j` for
`Z_7^*`-multiplier invariance of its residues mod 7.

## The answer: NO -- the cores are not Z_7-invariant
Only **2 of 16** cores are `Z_7^*`-invariant -- exactly the level-0 cores that happen to hit ALL residues
mod 7 (consec/skip-12 `O_0` = `{1,3,5,7,9,11,13}` -> `{0..6}`, the full group). Every deeper core
(`{1,3,5}`, `{1,3}`, `{1}`, `{1,3,5,91}` -> `{0,1,3,5}`, ...) has residues that are NOT a multiplier-orbit
union. **So the bare Z_7 cyclotomic SOS does NOT directly close `rho_j>=c`** -- the descended cores lack
the transitive `Z_7` symmetry the cyclotomic certificate needs.

## The mechanism is Gamma_0(14) congruence-AVERAGING
The cure (klein-S7b's fallback, and HYP-3553's `Gamma_0(14)`) is to AVERAGE the apex Gram over the
`Z_7^*` multiplier -- the congruence-averaging that MANUFACTURES the missing symmetry. VERIFIED on the
non-invariant core `{1,3,5}` (`= O_1` of consec): the raw apex autocorrelation Gram has spectrum
`[9, 0.31, 0.64, 5.05, 5.05, 0.64, 0.31]` -- NON-flat, SET-DEPENDENT; after `Z_7^*`-averaging it becomes
`[9, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]` -- FLAT off-0, the cyclotomic/octonion-optimal form,
**SET-INDEPENDENT** (the off-0 value depends only on `|O_j|`, not the residues). So averaging over the
congruence turns a config-dependent Gram into the flat octonion-optimal one -- literally "manufacture the
transitive symmetry it lacks" (klein-S5). (Caveat: complement-of-a-point cores like `{0..5}` are already
flat by accident; it is set-INVARIANCE, not flatness, that fails -- and averaging restores it universally.)

## What it sharpens
- The floor's close is the `Gamma_0(14)` CONGRUENCE-AVERAGED second moment, NOT a bare cyclotomic SOS:
  the descended cores are not `Z_7`-invariant, so you must average over the multiplier (= the congruence
  subgroup) to get the set-independent flat Gram. This confirms the `Gamma_0(14)` route (HYP-3553) is THE
  mechanism and rules out the simpler bare-`Z_7` hope.
- It MATCHES klein-S8 (HYP-3571): the empirical set-independent `inf R' = 0.344 >= 1/(2 zeta(2))` is the
  `Gamma_0(14)/zeta(2)` averaged bound -- the averaging is exactly why `R'` is set-independent while the
  unaveraged `CV(N_R)^2` is not.
- The remaining step is unchanged but now precisely scoped: the literal Han-Lee `Gamma_0(14)` second-moment
  CONSTANT (the averaged Gram's gap `>= 1/(2 zeta(2))`), which the octonion/perfect-difference-set flat
  spectrum (HYP-3575) shows is the optimal target.

## Honest status
A clean computational answer to klein-S7b's test (NO, not bare-`Z_7`; YES, `Gamma_0(14)`-averaged), which
sharpens but does not close the floor: it identifies the exact mechanism (congruence-averaging) and rules
out the simpler one, leaving the Han-Lee constant as the final step. Not a proof; a precise next-step map.
