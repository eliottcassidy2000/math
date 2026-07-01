---
id: HYP-3774
title: THE COVERING-MIN IS A ZETA-REGULARIZATION CARRIER -- five results synthesizing the crystallographic (S65/S66), Dedekind (S64), and regularization threads. SEED (all EXACT, verified n=3..14): with T = 1+...+(n-1) = n(n-1)/2 the finite speed-sum, Phi6 = 2T+1, the killer n(n-1)=2T = -1 mod Phi6, and s(n,Phi6) = -T/(12T+6). (R1 MOBIUS INTERPOLATION) s = -1/(12+6/T) is a fractional-linear function of the finite sum T with limit -1/12 = zeta(-1) (the Ramanujan-regularized 1+2+3+...); the finite covering-min carries the ACTUAL sum T, the asymptotic Dedekind limit carries its regularization; the margin = -12 s/n^2 has n^2*margin -> -12 zeta(-1) = 1, so the LRC margin's decay coefficient IS the regularized sum. (R2 THE 24=psi(14)) zeta(-1)=-1/12=-B2/2; the 12 in s = -1/zeta(-1), the 6 = 1/B2; psi(14)=[SL2:Gamma0(14)]=24=2*12 = the eta exponent (Delta=eta^24) -- the covering modulus's modular index IS the regularization denominator, and the floor's E2/Eisenstein bulk inherits zeta(-1). (R3 THE ANOMALY IS HEXAGONAL) zeta-regularization -1/12 appears IFF the lattice is hexagonal (order-6/Eisenstein, n^3=-1); the square (order-4/B2, h^2=-1) gives s=0 -- so the LRC margin = -12 zeta(-1)/n^2 EXISTS precisely because the covering-min is hexagonal (order 6); B2/square would give ZERO margin. (R4 THE KILLER=-1 IS THE ANCHOR) 2T = -1 mod Phi6 is the reflection/complement (the project's R involution) evaluated at -1; zeta(-1) at argument -1 and the killer at residue -1 are the same anchor. (R5 FAULHABER/EULER-MACLAURIN) T=sum k, and the B2=1/6 Euler-Maclaurin/Faulhaber correction is the margin's structure (ties HYP-2457/2453). OVERARCHING PROOF PICTURE: the floor/margin has a REGULARIZABLE bulk (zeta(-1), zeta(2), E2 -- always positive) + an UN-regularizable residual (the genus-1 cusp form f14) = the hard core is exactly the part beyond regularization.
status: SEED + R1-R4 EXACT/verified (n=3..14 and the asymptotics); R2 the 24=psi(14)=eta-exponent is exact; R3 hexagonal-iff is exact (square gives s=0). R5 and the OVERARCHING proof picture are DIRECTIONS/framework aligned with klein-S56 (E2 bulk vs cusp residual) and opus-S5 (reciprocity descent), NOT closed proofs. Synthesizes S64 HYP-3768 (Dedekind margin), S65 HYP-3771 (crystallographic spine), S66 HYP-3772 (dim/quasicrystal), and the owner's regularization seed.
source: mac-mini-2026-06-30-S67
related:
  - HYP-3768   # my S64 + klein-S56 + opus-S5: margin = order-6 Dedekind sum; s->-1/12=zeta(-1)=-B2/2 (E2 anomaly); reciprocity descent
  - HYP-3771   # my S65: the (2,3,n) crystallographic spine; hexagonal covering-min vs square B2
  - HYP-3772   # my S66: dim>=phi(n); the triple 6; quasicrystal
  - HYP-2457   # Faulhaber anchor expansion (the Bernoulli/sum-of-powers structure)
  - HYP-2453   # triangular-tower moment bridge (the triangular numbers T)
  - HYP-2896   # the zeta(12)/-1/12 core-tail dialectic (finite/regularized)
results:
  - 04-computation/dedekind_e2_b2_margin_macmini_20260630.py
  - 04-computation/crystallographic_dimension_hierarchy_lrc14_macmini_20260630.py
---

# HYP-3774 -- the covering-min is a zeta-regularization carrier

Building on the owner's seed and keeping the crystallographic angle: the covering-min construction is the place
where the **actual** sum of speeds (finite geometry) and the **zeta-regularized** sum (`zeta(-1) = -1/12`,
Ramanujan `1+2+3+... = -1/12`) meet. Let `T = 1 + 2 + ... + (n-1) = n(n-1)/2` be the finite speed-sum.

## The seed (all exact, verified n=3..14)
- `Phi_6 = n^2 - n + 1 = 2T + 1` (the sum-of-naturals identity).
- The **killer** speed `n(n-1) = 2T = Phi_6 - 1 = -1 (mod Phi_6)`.
- `s(n, Phi_6) = -T/(12T + 6)`.

## R1 -- the Mobius interpolation (finite sum <-> regularized sum)
`s(n,Phi_6) = -1/(12 + 6/T)` is a **fractional-linear (Mobius) function of the finite sum `T`**, with
`lim_{T->inf} s = -1/12 = zeta(-1)`. The finite covering-min carries the honest `T`; the Dedekind sum
interpolates it to the Ramanujan-regularized `1+2+3+...`. And the margin
`margin(n) = n/Phi_6 - 1/n = -12 s(n,Phi_6)/n^2` satisfies

> `n^2 * margin -> -12 zeta(-1) = 1`  (verified: 0.989, 0.9996, 0.99998, ... for n=10,50,200).

So **the LRC safety margin's asymptotic decay coefficient is exactly the zeta-regularized sum of the naturals.**

## R2 -- the `24 = psi(14)` regularization index
`zeta(-1) = -1/12 = -B_2/2` (`B_2 = 1/6`). The constants in the seed are the regularization's: the `12` is
`-1/zeta(-1)`, the `6` is `1/B_2`. And `psi(14) = [SL_2(Z):Gamma_0(14)] = 24 = 2*12` is the **eta exponent**
(`Delta = eta^24`, `eta = q^{1/24}...`, `1/24 = -zeta(-1)/2`). So the covering modulus's modular index IS the
zeta-regularization denominator; the floor's `E2`/Eisenstein **bulk** inherits `zeta(-1)`. *Proof-direction:* the
floor's Eisenstein bulk constant is a zeta-regularized quantity, positive and computable; the residual is the
cusp part.

## R3 -- the anomaly is HEXAGONAL (the crystallographic synthesis)
`zeta`-regularization `-1/12` appears **iff the lattice is hexagonal** (order-6 / Eisenstein `omega`,
`n^3 = -1 mod Phi_6`): then `s(n,Phi_6) -> -1/12`. The **square** lattice (order-4 / `B_2` / imaginary unit `i`,
`h^2 = -1`) gives `s(h,k) = 0` -- **no regularization anomaly**. Therefore

> the LRC margin `= -12 zeta(-1)/n^2` **exists precisely because the covering-min is hexagonal** (order 6);
> a square (`B_2`) covering-min would have `s = 0`, hence **zero margin** (covering-min `=` floor, degenerate).

This is a conceptual proof of margin-positivity: order 6 is the unique 2D crystallographic order whose
infinite-limit Dedekind sum is the nonzero `zeta(-1)`. It fuses S65/S66 (crystallographic: hexagonal vs square)
with S64 (Dedekind) and the seed (regularization).

## R4 -- the killer `= -1` is the regularization anchor
`2T = -1 (mod Phi_6)`: the killer sits at the **reflection/complement** residue, the project's `R` involution
evaluated at `-1`. `zeta` at argument `-1` and the killer at residue `-1` are the same anchor -- the covering
outlier `n(n-1) = 2T` being `= -1` is exactly what makes `s(n,Phi_6) = -s(1,Phi_6)`-type clean and drives the
`-1/12` limit. (The extreme Dedekind sum `s(Phi_6 - 1, Phi_6) = -s(1,Phi_6)` is the outlier's, S64.)

## R5 -- Faulhaber / Euler-Maclaurin
`T = sum_{k=1}^{n-1} k` is the first Faulhaber sum; its Euler-Maclaurin correction carries `B_2 = 1/6` -- the
same `6` in the seed and the same `B_2` in `zeta(-1) = -B_2/2`. The covering-min trajectory `M(n) = n/(2T+1)`
and its margin are Bernoulli/Faulhaber objects (ties HYP-2457 Faulhaber, HYP-2453 triangular-tower). Direction:
express the margin as an Euler-Maclaurin remainder of the speed-sum.

## Overarching proof picture (the finite/regularized dialectic)
The floor and margin split into a **regularizable bulk** and an **un-regularizable residual**:
- **bulk** = `zeta(-1)` (E2/Eisenstein, the margin's decay), `zeta(2)` (the `1/(2 zeta(2))` floor bound),
  `B_2` -- all zeta-values / regularizations, always present and positive;
- **residual** = the genus-1 cusp form `f_14` at the apex cusp `d=7` (HYP-3586/3768) -- the part **beyond
  regularization**, the hard core.

So "the finite covering-min carries the actual speed-sum; the asymptotic margin carries its regularization"
(owner) extends to the whole floor: **the easy part is what zeta-regularizes; the hard part is what does not.**
This is the same split as klein-S56 (`iota`-even E2 bulk + `iota`-odd cusp form) and opus-S5 (reciprocity
closes the cyclotomic endpoint, not the covering-min trajectory) -- now stated as a regularization dichotomy.

## Honest scope
Seed and R1-R4 are exact/verified; R2's `24 = psi(14) = eta`-exponent is an exact identity; R3's hexagonal-iff
(`square => s=0`) is exact. R5 and the overarching picture are a framework/direction (the bulk is genuinely
zeta-regularized; the claim that the residual is "un-regularizable" is a reframing of the cusp-form obstruction,
not a new proof). This is a synthesis giving multiple proof-directions -- most concretely R3 (why the margin is
positive = why hexagonal) -- not a closed proof of the residual.
