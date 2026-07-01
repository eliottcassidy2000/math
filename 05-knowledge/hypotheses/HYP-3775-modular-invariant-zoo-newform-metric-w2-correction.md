---
id: HYP-3775
title: A ZOO of X_0(N) invariants ("metrics like genus") + a SHARPER frontier metric + the rho<->R lever DEGREE CORRECTION. (A) SHARPER METRIC: the dim of NEW cusp forms at level 2p, dim S_2^new(2p) = g(2p) - 2 g(p) = 0,0,1,0,2 for p=3,5,7,11,13 -- SHARPER than the raw genus (0,0,1,2,2, HYP-3587): n=14 is the FIRST 2p with a GENUINELY NEW obstruction (f_14 = curve 14a), while n=22 is ALL-OLD (g(22)=2 but 2g(11)=2, J_0(22)~J_0(11)^2, ZERO new forms). This pinpoints n=14 as the frontier better than the genus. (B) LEVER CORRECTION (honest, refines HYP-3773): THM-580's 2-adic descent is a DEGREE-2 measure-preserving cover (u=2t on the circle); my S58 used the DEGREE-3 modular degeneracy X_0(2p)->X_0(p) (Gamma_0(2) index=3). Different covers. The degree-2 candidate is the ATKIN-LEHNER involution W_2, with 4 fixed points on X_0(14) (RH, quotient X_0(14)/W_2 genus 0; CM points). So S58's R=6 is the (real) ramification of the degeneracy, NOT of THM-580's degree-2 descent; the descent = W_2. (C) THE ZOO: cusps nu_inf(2p)=4 CONSTANT (= the Eisenstein bulk, dim=cusps-1=3, HYP-3587); index psi(2p)=3(p+1); psi(14)=24=eta exponent; orbifold -chi_orb=psi/6; nu2,nu3 elliptic; Phi6(n) factorization (large prime = covering-min binding core); Dedekind s(n,Phi6). Diverse metrics, all computed exact
status: MIXED (exact zoo + one sharper metric + honest correction). PROVED/standard: all classical invariants (genus formula sanity-checked; cusps, elliptic, index); dim S_2^new(2p)=g(2p)-2g(p) (oldform multiplicity 2 from level p) = 0,0,1,0,2 (VERIFIED); W_2 fixed points = 4 on X_0(14) (Riemann-Hurwitz, quotient genus 0 known). CORRECTION to HYP-3773: THM-580 is degree 2 (u=2t), the degeneracy is degree 3 -- S58's R=6 is the degeneracy's ramification, not the descent's; the degree-matched descent object is W_2 (4 fixed points). SYNTHESIS: the newform-count-as-frontier-metric and the cusps=constant-bulk reading. The rho_j<->(W_2 fixed points) numerical correspondence remains OPEN (rho real in [.5,1] vs 4 integer).
source: klein-2026-06-30-S59
depends_on:
  - HYP-3773   # my S58 Riemann-Hurwitz descent (this CORRECTS its degree)
  - HYP-3587   # genus = obstruction dim (this SHARPENS it to newform count)
related:
  - THM-580    # the 2-adic descent (degree 2 = u=2t = W_2, not the degeneracy)
  - HYP-3768   # covering-min = E2/Eisenstein; residual = f_14 (the new form here)
  - HYP-3774   # mac-mini: covering-min = zeta-regularization carrier (bulk + residual)
  - HYP-3041   # cusps of X_0(14) = Klein-4 = Atkin-Lehner (2/Z2)^2 (the W_2)
  - HYP-3732   # covering-min Farey rung (Phi6 large-prime core)
results:
  - 04-computation/modular_invariants_zoo_klein.py
  - 05-knowledge/results/modular_invariants_zoo_klein.out
---

# HYP-3775 — a zoo of X_0(N) invariants, a sharper frontier metric, and the lever correction

## (A) A sharper "metric like genus": the NEW-form count pinpoints n=14
The raw genus of `X_0(2p)` is `0,0,1,2,2` (HYP-3587). But the LRC-relevant obstruction is the *new*
cusp forms at level `2p` (the old ones are inherited from level `p`):
```
    dim S_2^new(2p) = g(2p) - 2 g(p)   (oldform multiplicity 2 = #divisors of 2p/p)
      p     :  3   5   7   11  13
      g(2p) :  0   0   1   2   2
      2g(p) :  0   0   0   2   0
      NEW   :  0   0   1   0   2
```
**`n=14` is the FIRST `2p` with a genuinely NEW obstruction** (`f_14` = elliptic curve `14a`), while
**`n=22` is ALL-OLD** (`g(22)=2` but `J_0(22) ~ J_0(11)^2`, zero new forms -- its genus is inherited
from `p=11`). So the new-form count is a **sharper frontier metric than the genus**: it isolates `n=14`
as the first level carrying its own obstruction (matching mac-mini HYP-3771 "`p=7` is the frontier").
The residual `f_14` (S56/HYP-3768, the `iota`-odd degree S55, the ramification S58) is this new form.

## (B) The rho<->R lever: an honest DEGREE correction to HYP-3773
THM-580's 2-adic descent is the **degree-2** measure-preserving cover `u = 2t` on the circle (exact,
`meas(lonely E)=meas(lonely E/2)`). My S58 (HYP-3773) matched it to the modular degeneracy
`X_0(2p) -> X_0(p)`, which has **degree 3** (`= [Gamma_0(p):Gamma_0(2p)] = ` index of `Gamma_0(2)`).
**Different covers.** So S58's `R=6` is the (genuine) ramification of the *degree-3 degeneracy*, but NOT
of THM-580's *degree-2* descent.

The **degree-matched** object is the **Atkin-Lehner involution `W_2`** (degree 2). On `X_0(14)` (genus 1
= curve `14a`), `W_2` has **4 fixed points** (Riemann-Hurwitz: `2-2*1 = 2(2-2*0) - f => f=4`, since
`X_0(14)/W_2` has genus 0), which are CM points (discriminants among `-4,-8,-7,-56`, `h(-56)=4`). So the
corrected reading: **the 2-adic descent = the `W_2` involution (degree 2, 4 CM fixed points), quotient
`X_0(14)/W_2` genus 0.** The `4` is the ramification of the degree-matched cover (vs S58's `6` for the
degree-3 degeneracy). HYP-3041 confirms `W_2` lives in the Klein-4 Atkin-Lehner group of `X_0(14)`.
[Honest: the numerical `rho_j <-> 4` correspondence is still open -- `rho` is real in `[0.5,1]`, `4` is
an integer count; what is matched is the geometric object (degree-2 involution, CM fixed points).]

## (C) The invariant zoo (diverse, exact)
```
  N    g  cusps nu2 nu3 index psi  psi/6   note
  6    0    4    0   0    12      2      2p, p=3
 10    0    4    2   0    18      3      2p, p=5
 14    1    4    0   0    24      4      psi=24 = eta exponent (Delta=eta^24)
 22    2    4    0   0    36      6
 26    2    4    2   0    42      7
```
- **cusps `nu_inf(2p) = 4` CONSTANT** across the family = the **Eisenstein bulk** (dim `= cusps-1 = 3`,
  HYP-3587): the bulk is the *constant* invariant, the genus/newforms the *growing* obstruction -- the
  local-global split read directly off the invariants.
- **index `psi(2p) = 3(p+1)`**; `psi(14) = 24` = the exponent in `Delta = eta^24` (mac-mini HYP-3774's
  `eta` exponent) `= 2*12 = -2/zeta(-1)`.
- **orbifold `-chi_orb = psi/6`** (hyperbolic area `/pi = psi/3`); `nu2,nu3` elliptic (the `(2,3,p)`
  cone points, mac-mini HYP-3771).
- **`Phi_6(n)` factorization**: `183=3*61, 651=3*7*31, ...` -- the large prime factor is the
  covering-min binding's arithmetic core (Farey rung HYP-3732); **`Dedekind s(n,Phi6)`** as in HYP-3768.

## Net
Searching broadly for invariants "like genus": the sharpest LRC frontier metric is the **new-form count**
`g(2p)-2g(p) = 0,0,1,0,2`, isolating `n=14` as the first level with its own obstruction (`n=22` is
all-old). The invariant zoo shows the **cusps are constant (Eisenstein bulk) while the newforms grow
(the obstruction)** -- the local-global split as invariants. And the `rho<->R` lever is corrected:
THM-580's degree-2 descent is the **Atkin-Lehner `W_2`** (4 CM fixed points on `X_0(14)`), not S58's
degree-3 degeneracy (`R=6`) -- both real covers, but only `W_2` is degree-matched to the 2-adic peel.
