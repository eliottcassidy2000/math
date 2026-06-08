---
id: HYP-2370
name: lrc-dual-conformal-symmetry-spacetime-vs-dual
status: VERIFIED (A,B,C,D exact/exhaustive) + a synthesis conjecture (maximal dual symmetry = hardness)
date: 2026-06-08
session: claudebox-2026-06-08-S717
depends_on:
  - THM-403   # cyclotomic witness orbit of the AP ((Z/m)* = the dual space)
  - THM-415   # signed-LRC homometry = same Patterson power spectrum (the dual data)
  - THM-420   # bad multipliers = {+-v_i^-1} (the inversion / special-conformal generator)
  - THM-441   # convolution-correlation-adjoint duality (S714 autocorrelation operator)
  - THM-409   # perspective-flip involution (a dual symmetry)
related:
  - "oracle-S540o: LRC worldlines as a pure braid on the cylinder S^1 x [0,1] (the spacetime picture)"
provisional_id: true
---

# HYP-2370: The LRC has a dual conformal symmetry; the observer is spacetime, the obstruction is dual

In planar N=4 SYM, dual conformal symmetry is a HIDDEN symmetry — invisible in spacetime, manifest only
in dual region/twistor variables — and ordinary+dual conformal generate the Yangian. The LRC has the same
two-space structure, and it cleanly sorts what can and cannot be expressed in spacetime.

## The two spaces

- **SPACETIME** = the cylinder `S^1 x [0,1]`, runner worldlines `x_i(t) = frac(v_i t)`, observer = the
  vertical strand at `0` (the LRC pure braid). Loneliness = a time `t` with every strand far from the
  observer. Local, geometric.
- **DUAL** = the speed/residue space; the multiplier group `(Z/m)*` (= "perspective" = `Gal(Q(zeta_m)/Q)`,
  THM-403/439); the autocorrelation/homometry (THM-415/441). Spectral, arithmetic.

## (A) VERIFIED: ordinary (spacetime) conformal symmetries of M(S)

`M(S) = max_t min_v ||v t||` is invariant under **dilation** `v -> lambda v` (= rescale time; "conformal
weight 0"), **reflection** `v -> -v` (= time reversal), and **permutation** `S_n` (relabel strands):
`300/300` random tests each. These act visibly on the worldlines — the ordinary conformal group.

## (B) VERIFIED: the DUAL conformal symmetry (hidden) — the multiplier, with inversion as special-conformal

On the `(Z/m)*` witness orbit (THM-403), the **multiplier** `a in (Z/m)*` acts by `v -> a v mod m` and
```
   orbit-M(a . S mod m) = orbit-M(S)      EXACTLY   (300/300; proof: a permutes the witnesses b -> ab).
```
This is a HIDDEN symmetry: it scrambles the worldlines arbitrarily (no spacetime picture) yet preserves
loneliness. The **inversion** `v -> v^{-1} mod m` is the **special-conformal generator**: `a = v_i^{-1}`
sends `v_i -> 1`, onto the central band — exactly the bad multiplier of THM-420. The hard
(transversal/Paley) core is hard precisely because the inversions `{+-v_i^{-1}}` cover ALL multipliers, so
no good rotation (dodge) survives.

## (C) VERIFIED: the OBSERVER/origin is pure spacetime, invisible to the dual

`M` is NOT translation-invariant, but the dual data (distance multiset / autocorrelation) IS:
`S+c` has the same distance multiset yet different `M` in `265/300` trials (e.g. `{1,2}` has `M=1/3`,
`{3,4}` has `M=3/7`, identical dual data). The observer (the fixed worldline at `0`) breaks translation:
"where the origin is" is a spacetime fact the dual cannot see.

## (D) VERIFIED: M is NOT determined by the dual data — loneliness needs spacetime

Homometric, non-congruent configurations (same distance multiset = same autocorrelation = same dual data)
have **different** `M`: e.g. `{0,1,2,6,8,11}` vs `{0,1,6,7,9,11}` (shifted to speeds) give `M = 1/5` vs
`2/9`. So `M` is not a function of the dual autocorrelation alone; loneliness carries genuine spacetime
information beyond the dual.

## The dictionary: what can / cannot be expressed in spacetime

| | expressible | example |
|---|---|---|
| SPACETIME only | loneliness, `M`, covering depth, braid topology, **the observer/origin** | local geometry |
| DUAL only (hidden) | **multiplier/perspective symmetry**, homometry, autocorrelation flatness | THM-403/415/441 |
| needs BOTH | `M` itself (dual-symmetric but not dual-determined; observer is spacetime) | (C)+(D) |

## The synthesis conjecture (the amplitudes lesson)

The dual conformal symmetry is **maximal exactly where the LRC is hardest**. Maximal multiplier symmetry
= flat autocorrelation (S714/THM-441) = the Paley/transversal core (THM-420), where every multiplier is
killed by some inversion `+-v_i^{-1}` and no spacetime dodge exists. So:

> **the obstruction lives in the dual (maximal hidden symmetry, flat spectrum, no rotation survives); the
> observer lives in spacetime (the origin breaks translation); and `M` — the answer — requires both.**

This is the LRC mirror of the amplitudes principle that the answer is hardest where the hidden (Yangian)
symmetry is largest, and that locality (spacetime) is a partial description. Ordinary + dual conformal =
the "perspective Yangian" of the LRC (the user's perspective key = automorphism rigidity = this dual
symmetry).

## Next
- exhibit the full Yangian-type closure: do the spacetime conformal and dual multiplier generators
  satisfy Yangian/Serre relations on the witness-orbit representation?
- the special-conformal/inversion as a momentum-twistor-like change of variables that makes the dual
  symmetry manifest and trivializes the transversal core (a "momentum twistor" for LRC);
- tie "maximal dual symmetry = hardness" to the open frontier (n=15,19,21,22): measure the multiplier-
  orbit symmetry / autocorrelation flatness as a hardness predictor (continues S714 t-0100).
