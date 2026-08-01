---
source: kind-pasteur-2026-07-24-S155 (Opus 4.8)
status: RESULT (concrete transfer test, honest NEGATIVE with mechanism). Ran the proposed test -- evaluate the
  repo's dimension-3 Special-Image-Conjecture witness (THM-2801 f_1) under the Factorial-Conjecture functional.
  REPRODUCED THM-2801 exactly, then SHOWED the witness does NOT transfer to FC(3): its factorial-moment vanishing
  requires OFF-DIAGONAL (angular/Hopf) terms, while FC lives on the DIAGONAL subalgebra (radial). The diagonal
  projection X1(X2-X3) has non-vanishing factorial moments. This explains, structurally, why FC(2) can be true
  (radial KZ-rigidity) while SIC(2)/GMC are false (angular escape) -- and refines kps-S154 lever 2.
tags: [factorial-conjecture, image-conjecture, mathieu-zhao, gaussian-moments, periods, transfer-test, series]
related: [kps-S154, THM-2801, THM-2022, THM-1790]
refines: [kps-S154 lever-2]
---

# The repo's dim-3 image-conjecture witnesses do NOT transfer to FC(3)

## 0. Setup: FC = the diagonal restriction of the SIC/Gaussian functional
FC(m): `f in C[X_1..X_m]`, `L(f^k)=0 for all k` `=>` `f=0`, `L(X^a)=prod a_i!`. Since
`a! = E[(g bar g)^a]` for a standard complex Gaussian, `L(f^k)=E[iota(f)^k]` where the **diagonal embedding**
`iota(f)(g,bar g)=f(g_1 bar g_1,...,g_m bar g_m)` uses only the products `xi_i z_i` (`xi=bar g`, `z=g`).
THM-2801's operator `E_n` on balanced polys is exactly this Gaussian/factorial contraction (`F_n(xi^a z^b)=a! d_{ab}`, eq 5).
> **So FC(m) is precisely `E_n`/GMC restricted to the DIAGONAL subalgebra `C[xi_1 z_1,...,xi_m z_m]`**
> (functions of the moduli `|g_i|^2` only -- the radial part). SIC/GMC allow off-diagonal `xi_i z_j` (`i!=j`).

## 1. Reproduced THM-2801 exactly (machinery validated)
Witness `f_1 = tau(t+z)(w t - v y)`, pairs `(tau,t),(w,z),(v,y)` (THM-2801 eq 29). Computed exactly (Python, integer):
- `E_3(f_1^m) = 0` for `m=1..6`. (THM-2801 eq 27.) [OK]
- `E_3(z f_1^m) = (m+1)! m! * t + [m=1] z`: got `2t+z, 12t, 144t` for `m=1,2,3`. (THM-2801 eq 30.) [OK]
So `f_1^m in ker E_3` for all `m` but `z f_1^m` escapes -- the Mathieu-Zhao failure, `not SIC(3)`.

## 2. The transfer FAILS -- and exactly why
Split `f_1` by the pairing (`alpha=beta` vs `alpha!=beta`):
```
DIAGONAL part  (alpha=beta):   tau*w*t*z - tau*v*t*y   = P1 P2 - P1 P3   (P_i = xi_i z_i)  = X1(X2 - X3)
OFF-DIAGONAL   (alpha!=beta):  tau*w*t^2 - tau*v*z*y
```
The FC functional `L` sees only the diagonal projection. Its factorial moments are
```
L( (X1(X2-X3))^k ) = 0, 4, 0, 576, 0, 518400, 0, ...   (NOT all zero; even moments nonzero)
```
(and `L((X1(X2+X3))^k) = (k+1)(k!)^2`, a clean check). **`f_1`'s moment-vanishing `E_3(f_1^m)=0` is destroyed by
dropping the off-diagonal terms** `tau*w*t^2, tau*v*z*y`. Those off-diagonal terms are the **angular / Hopf phase**
freedom that THM-2801 sec 2.1 integrates out (its `A_m(1+u)`, `(1-t^2)^m`, `u^m`-coefficient mechanism, eqs 14-20).
> **The SIC(3) witness is not a diagonal polynomial, so it is not an FC(3) object at all; and its breaking
> mechanism (off-diagonal angular cancellation) is exactly what FC forbids.** No FC(3) counterexample results.

## 3. Why this is the right answer (and it explains the friend's claims)
FC's diagonal/radial restriction removes the angular escape that breaks SIC/GMC. That is the structural reason:
- **FC(2) true (friend: implied by KZ).** On the diagonal, `L(f^k)=int f(Y)^k e^{-|Y|}dY` (`Y_i=|g_i|^2 ~ Exp`)
  is a purely radial exponential period; KZ-rigidity has no angular direction to hide a cancellation, so vanishing
  forces `f=0`.
- **SIC(2)/GMC false.** The bigraded algebra has the off-diagonal `xi_i z_j` (Hopf-angular) directions, and the
  witnesses live there (THM-2801 sec 7: "the first failure lives on the rank-one Segre cone ... a CP^1 coefficient
  problem" -- angular).
- **FC(3) "may still be false" (friend).** If so, it needs a **genuinely radial/diagonal** mechanism, NOT the
  repo's angular SIC/GMC one. A quick scan of small diagonal `f(X1,X2,X3)` (`X1-X2`, `X1X2-X3^2`, `X1-X2 X3`, ...)
  finds none with all-vanishing factorial moments -- FC(3) resists cheap counterexamples.

## 4. Consequence for the program (refines kps-S154 lever 2)
kps-S154 lever 2 optimistically listed the repo's dim-3 SIC/GMC counterexamples as "candidate sources for an
FC(3) counterexample." **Refined honest verdict: they do NOT transfer** -- they are off-diagonal, FC is diagonal.
The value is the sharpened target:
> An FC(3) counterexample must be a **diagonal** polynomial `f(X1,X2,X3)` (radial in `|g_i|^2`) whose exponential
> moments `L(f^k)=int f^k e^{-|Y|}` all vanish -- a purely radial/period phenomenon. The repo's angular machinery
> (Hopf reduction, THM-2801 sec 2.1) is the wrong tool; the right one is exponential-period rigidity (KZ) applied
> to a candidate radial `f`. This is exactly the JC<->period edge kps-S154 named.

MISTAKE-discipline note: this is a **negative** transfer, honestly recorded (cf. the repo's Redei-parity-Mathieu
negative). The bridge kps-S154 survives -- it predicted the *machinery* is shared (Gamma/exponential moments,
KZ), and it correctly located FC as the radial restriction; it did NOT promise the specific witnesses port, and
they don't.

Files: `/tmp/{fc3,fc3b}.py`. Reproduces THM-2801 eqs 27/30 exactly; refines kps-S154.
