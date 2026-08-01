---
source: kind-pasteur-2026-07-24-S166 (Opus 4.8)
status: SYNTHESIS + reframe (identifies the missing tool). The FC(3) residual my methods left open -- an ISOLATED
  non-composition point where an exponential period vanishes (the "KZ coincidence", kps-S160) -- is PRECISELY an
  E-function / exponential-period nonvanishing, the exact object the external `int_0^1 e^Q != 0` proof addresses
  via Siegel-Shidlovskii / Beukers. So E-function transcendence theory is the missing piece that, combined with
  opus's saddle-cap (Morse case) and my transversality (finite degree), would COMPLETE FC(3). Same family closes
  the series non-elementarity (G-functions). Separately: AMM's golden gamma* is an LP-bound optimum in the exact
  sense of Cohn-Elkies sphere-packing (ten-proofs #1), with the artanh certificate a weak dual.
tags: [factorial-conjecture, E-functions, transcendence, siegel-shidlovskii, sphere-packing, LP-bounds, synthesis, reframe]
related: [kps-S154, kps-S159, kps-S160, kps-S165, HYP-9078]
---

# E-functions as the missing piece for FC(3); the LP-bound analogy for AMM

The owner flagged the OpenAI "ten advances" (Lean-formalized, prose-provenance to be checked -- HYP-9078 sec 4
already caught one abstract/arXiv mismatch) and a friend's `int_0^1 e^{Q(t)}dt != 0` (transcendental) result.
Here is what that machinery unlocks in this repo.

## 1. The FC(3) residual IS an E-function nonvanishing
My FC(3) program (kps-S156-S162) reduced "no counterexample" to a single residual: an **isolated, non-composition,
non-transversal** `f` where the exponential period vanishes -- the "KZ coincidence" no algebraic, analytic, or
finite-field method reached. Name it exactly:
> A FC(m) counterexample is a nonzero polynomial `f` with `Phi_f(t) := int_{R_+^m} e^{t f - |x|}dx == const`, i.e.
> all `d^k/dt^k` at 0 vanish. `Phi_f(t)` is a **period of an exponential motive, holonomic in `t`** -- for the
> isolated point it is exactly an `E`-function rigidity, the SAME shape as `int_0^1 e^{Q}` (both: an exponential
> integral of algebraic data, forced constant/zero, must have transcendental "frequency" `2 pi i`).
The external proof's tools -- **Siegel-Shidlovskii, Andre, Beukers lifting** -- are precisely for E-function
nonvanishing. So:
> **CLAIM (reframe): E-function transcendence theory is the missing tool that closes the FC(3) isolated point.**
> Combined with (a) opus's saddle-cap -- Morse/isolated-max-modulus closed analytically at all degrees, on the
> `{J=0}` fold = my rank-drop locus -- and (b) my finite-degree transversality (kps-S159/S161, exact to D=5) and
> the good-prime search (kps-S162), the E-function layer handles the remaining non-Morse isolated locus. That is
> the completion path for FC(3), and it is the same machinery the friend's `int e^Q` proof runs.
Honest: this identifies the tool and the path; it does not execute Siegel-Shidlovskii here. But it converts
"an unreachable KZ coincidence" into "an E-function nonvanishing," which is a solved *kind* of problem.

## 2. Same family, the series non-elementarity
My series values `S_a(k)=int_0^1 2F1(a,1-a;1;x^k)dx` (kps-S146-S153) are **G-function** periods (the algebraic/
classical twin of E-functions). Their conjectured non-elementarity for `k>=P(n)` (which I could only evidence,
never prove -- kps-S153) is a **G-function transcendence** statement, in Andre's / Beukers-Brownawell-Heckman
theory. The `int e^Q` result is the E-side demonstration that this family of machinery is now AI-tractable; the
G-side is the exact tool for the series locus. **One transcendence toolbox, two of the repo's open non-elementary
walls (FC(3) isolated point; series `C_{1/n}` upper bound).**

## 3. AMM's golden is a Cohn-Elkies-type LP optimum
Ten-proofs #1 (sphere packing to the Cohn-Elkies threshold) is an **LP duality**: a dual "magic function"
(`f-hat >= 0`, `f <= 0` past radius 1) certifies the packing bound, tight in dims 8, 24. AMM 12592's
`C* = 1 + gamma*` is the *same shape* (THM-2966): find nonneg box sequences `W_m, V_m` meeting the spine identity;
the minimal growth `gamma*` is an LP optimum. Its value is the **golden** `gamma* = log_5(phi^2) = 0.598`
(opus/klein THM-3027, min-max `x*=1/phi^2`) -- the "magic function" optimum. My **artanh two-bias log-odds
certificate is a weak feasible dual** (`gamma >= 2457/6592 = 0.3727`, kps-S164); the *tight* dual is the golden
Lagrangian. So the open AMM gap `[0.3727, 0.598]` is a duality-gap-closing problem of exactly Cohn-Elkies type --
find the extremal dual functional that matches the golden primal. **Transfer: the sphere-packing LP-bound
technology (positive-definite dual construction) is the natural tool to upgrade the artanh certificate to the
golden lower bound.**

## 4. Discipline (the ten-proofs lesson, owner-flagged)
The advances are *curated, formally-checked-but-prose-provenance-uncertain*; "fetch the primary, the capstone is
not the source" (opus/klein, MSG-1414). Same lesson as my S151 false negative and my kps-S163 golden-flag error
(now corrected, kps-S165). FC(2)-homogeneous is **Liu-Sun 2020 Thm 2.6** (cite, don't reprove); the `int e^Q =>
FC(2)` general implication is UNVERIFIED (HYP-9078, weights `j!` vs `1/(j+1)`).

## 5. Net (three concrete transfers)
1. **FC(3):** the isolated KZ point = an E-function nonvanishing -> Siegel-Shidlovskii/Beukers is the tool; with
   opus's saddle-cap + my transversality this is the completion path.
2. **Series `C_{1/n}`:** non-elementarity = G-function transcendence (Andre) -- same toolbox.
3. **AMM `[0.3727,0.598]`:** duality-gap closing of Cohn-Elkies type -- positive-definite dual to match the golden.
All three are "algebraic data cannot realize a transcendental (`2 pi i` / period) resonance" -- the continuous
Conway-Jones governor (kps-S165).

Files: verify inline. Builds on kps-S154/S156/S160/S165; engages HYP-9078, THM-2966/3027, Liu-Sun 2020, ten-proofs.
