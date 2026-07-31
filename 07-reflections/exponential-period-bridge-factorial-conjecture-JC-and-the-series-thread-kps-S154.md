---
source: kind-pasteur-2026-07-24-S154 (Opus 4.8)
status: SYNTHESIS + concrete levers. Identifies an UNOCCUPIED bridge in the repo: the Factorial/SFC/Mathieu-Zhao/
  Jacobian cluster has no connection to period theory, while period/motive machinery sits siloed in the series
  (S_lambda(k)) thread. The Kontsevich-Zagier EXPONENTIAL-PERIOD route to FC(2) (a friend's approach; 0 repo
  hits for "exponential period"/"Kontsevich-Zagier period") is exactly that bridge. Gives the disciplined
  transfer (map/predicate/loss/control), three actionable levers, and honest S(4)/S_{1/3}(3) period tests.
tags: [factorial-conjecture, jacobian, mathieu-zhao, periods, kontsevich-zagier, furter, series, synthesis, bridge]
related: [kps-S146, kps-S153, THM-1300, THM-1790, THM-2836, THM-2801, CORE-PAPERS]
external: [arXiv:1304.3956 (Edo-van den Essen SFC), arXiv:2607.09340 (Adamus-Hajto Pascal-finite), Zhao Vanishing Conj]
---

# The exponential-period bridge: Factorial Conjecture <-> Jacobian <-> the series thread

## 0. The precise web (literature-confirmed)
- **Factorial Conjecture FC(m):** `f in C[X_1..X_m]`, `L(X^a)=prod a_i!`; if `L(f^k)=0` for all `k>=1` then `f=0`.
- **Key identity (the "Gamma Bridge"):** `a! = Gamma(a+1) = int_0^inf x^a e^{-x}dx`, so `L(g)=int_{R_+^m} g e^{-|x|_1}dx`
  and `L(f^k)=int f^k e^{-|x|}dx` is an **exponential period** (Kontsevich-Zagier: `int_sigma e^{-g}omega`, alg. data).
- **Edo-van den Essen (Strong FC, arXiv:1304.3956):** links FC to **Furter's Rigidity Conjecture**.
- **Zhao Vanishing Conjecture `<=>` Jacobian Conjecture**; FC sits under the Mathieu-Zhao / Image-Conjecture umbrella.
- **Friend's report:** FC(2) is (almost certainly) true because implied by the **KZ exponential-period conjecture**;
  FC(3) may be false; an FC(2) proof may partially address Furter rigidity and touch Keller pairs / planar JC.

## 1. The unoccupied gap (from a full-repo survey)
The repo already runs the **algebraic** half of this web: `SFC <-> GMC (Gaussian moments) <-> Mathieu-Zhao <-> JC`
(THM-1790 sec5 "Boundary with the Strong Factorial Conjecture"; CORE-PAPERS:370; the Gamma-bridge court case
CASE-gamma-bridge). But **"exponential period" and "Kontsevich-Zagier period" have ZERO hits** in canon, and all
genuine period/motive material lives **siloed** in my series thread (kps-S146/S148-S153: `S_lambda(k)` as periods
of hypergeometric motives). The one "Kontsevich" near JC is Belov-Kanel-Kontsevich (automorphisms), NOT periods.
> **The friend's FC(2)-via-KZ is precisely the missing edge**: it connects the repo's JC/SFC/Mathieu-Zhao cluster
> to period theory. That edge is currently unbuilt in-repo.

## 2. One machine, two weights (the shared structure)
Both problems apply a **moment functional** `M_w(x^a)=int x^a w(x)dx` to a polynomial/series and ask what its
values detect:
| | weight `w` | period type | "trivial" case | "new" case |
|---|---|---|---|---|
| **Factorial / FC** | `e^{-|x|}` | exponential period | `f` linear -> `Phi_f(t)=int e^{tf-|x|}` rational (e.g. `1/(1-t^2)`) | `f` nonlinear -> Airy-type exp. period |
| **Series (mine)** | `1_{[0,1]}` (`1/(kn+1)=int_0^1 x^{kn}`) | classical period | `k=1,2` mixed-Tate / cyclotomic (proved, kps-S153) | `k>=3?` irreducible motive |

The **repo's Gamma Bridge `k!=Gamma(k+1)`** IS the exponential-moment weight; my series' `1/(kn+1)` IS the
classical-moment weight. Same functional family; the **only** difference is `e^{-x}` vs `1_{[0,1]}`. And the
governing dichotomy is identical: **elementary/trivial period (mixed-Tate, decidable) vs new period (irreducible,
transcendence-hard)**. KZ ("no relations beyond additivity + change-of-variables + Stokes") is the shared engine;
FC(2) is the statement that the exponential period `Phi_f` cannot be `trivially constant` unless `f=0`.

## 3. The disciplined transfer (per MISTAKE-226: map / predicate / loss / control)
- **MAP.** `Series-moment (w=1_{[0,1]}, classical) <-> Factorial-moment (w=e^{-x}, exponential)`, both
  `M_w(x^a)=int x^a w`.
- **PRESERVED PREDICATE.** "A moment value is elementary iff its motive is mixed-Tate/Artin; KZ-rigidity forbids
  accidental vanishing/constancy." (Proven on the series side for `k=1,2`; conjectural-but-KZ-backed on the FC side.)
- **LOSS.** The weight `e^{-x}` lifts classical periods to **exponential** periods. This is a genuine loss of ease:
  the classical side is (sometimes) settled by a **complete cyclotomic basis** (my k=2 proof), whereas the
  exponential side needs the full KZ exponential-period conjecture (unconditional only in low complexity). The
  transfer is **asymmetric**: series techniques discharge the *classical* sub-questions; FC's hard case needs the
  *exponential upgrade*. Do NOT claim series methods settle FC.
- **HOSTILE CONTROL.** "Does the bridge survive `S(4)` (classical) resisting evaluation while FC(2) (exponential)
  is provable?" **Yes** -- they are different *statements*: `S(4)` non-elementarity is a hard **negative**
  transcendence fact (my S151 lesson: even detecting it is unreliable), while FC(2)-via-KZ is a **positive
  rigidity** (vanishing => trivial). The bridge claims shared *machinery*, NOT equal difficulty. Control passes.

## 4. Three actionable levers (the payoff)
1. **Detection-floor audit of the repo's factorial-moment certificates.** My kps-S151 error was declaring a period
   non-elementary from an **incomplete PSLQ basis** (a false negative). The repo's SFC(3) certificate THM-2836
   ("SFC(3) holds on supports `<= z^9`, windows `k<=6`", via factorial-moment detection THM-2173/2812/2830/...)
   is the SAME failure mode: a **bounded** moment search. Concrete audit: does THM-2836 establish a *detection
   floor* -- a proof that the bounded window cannot MISS a witness -- or only "no witness found in-window"? If the
   latter, it carries my false-negative risk. (Matches memory `verify-inherited-residues`.)
2. **FC(3)-may-be-false `<->` the repo's explicit dimension-3 counterexamples.** The friend expects FC(3) false; the
   repo already has the dimension-3 wall for sibling moment conjectures: **THM-2801 (SIC(n>=2) FALSE)**, explicit
   **GMC(3) counterexamples**, **THM-1300 (JC false, n>=3)**. These concrete degree/dimension-3 non-Mathieu-Zhao
   polynomials are candidate **sources/mechanism** for an FC(3) counterexample (why dim 3 breaks moment-rigidity).
3. **The KZ route is an INDEPENDENT proof of FC(2), and it shores up a known repo gap.** The repo's `SFC(2)/GMC(2)
   => JC(2)` chain is honest-flagged as resting on **"cited-external, unverified-in-repo" IC=>JC arrows, no
   witness** (PROBLEM-LEDGER sec A6/A7). The friend's KZ/exponential-period route has **different foundations
   (transcendence, not Mathieu-Zhao)**, so it is not redundant -- it is a second pillar under the same statement,
   exactly where the algebraic pillar is externally-dependent.

## 5. Automorphism-side mirror (the arXiv paper)
arXiv:2607.09340 (Adamus-Hajto): a plane polynomial automorphism decomposes as `diag(det J,1)` composed with
**Pascal-finite** factors iff `det J = 1`; an "Exponential Generators Conjecture in positive characteristic."
This is the **automorphism analogue of "period = sum of elementary (exponential) pieces"**: `det J=1` is the
"trivial/unimodular" class, and Pascal-finite = exponential-of-LND generators are the elementary building blocks.
"Pascal finite" and "Exponential Generators" have **0 repo hits** -- another unbuilt edge, on the Furter/Jung-van
der Kulk side that the friend says an FC(2) proof would touch.

## 6. Concrete series tests (honest, disciplined)
Both resist elementary identification against INDEPENDENT 15-element bases at 175 digits, coeff `<= 1e7`:
- `S(4)=S_{1/4}(4)`: null vs `Q(sqrt2)`/conductor-8,16 cyclotomic + Catalan `G` + `Gamma(1/8)^2/pi` + lemniscate `varpi`.
- `S_{1/3}(3)`: null vs `Q(sqrt3)`/conductor-9 cyclotomic + Clausen `Cl_2(pi/3)` (`L(2,chi_-3)`) + `Gamma(1/3)^3`.
**Consistent with genuine NEW periods** (irreducible motives) -- the classical-period image of FC's "nonlinear `f`
gives a new exponential period." **But NOT proven non-elementary**; per kps-S151 a null is disciplined evidence,
never a certificate. (`k=2` is the only layer where I have a *proof*, because there I have the closed form.)

## 7. Engaging the friend's five claims (repo-grounded)
1. "FC(2) true, implied by KZ" -- consistent; the moments ARE exponential periods (sec 0, verified numerically).
2. "FC(2) falls shortly" -- if so, it becomes the repo's independent second proof of `SFC(2)` (lever 4.3).
3. "FC(3) may be false" -- **fits the repo's dimension-3 pattern** (lever 4.2); candidate witnesses in-house.
4. "Touches Furter rigidity" -- **literally the Edo-van den Essen SFC link** (arXiv:1304.3956); repo has it as the
   open "Furter transfer" (one hit, kls-2026-07-28).
5. "Keller pairs / planar JC implications" -- speculative; planar JC (Keller pairs, tame, THM-2063+) is a repo top
   priority, so any FC(2)->planar-JC arrow lands directly on the fleet-#1 witness-extraction lane.

## 8. Honest status
This is a **machinery bridge**, not a proved reduction; each edge names its map/predicate/loss/control (sec 3).
Value delivered: (a) the unoccupied JC<->period edge is named and its shape given; (b) three concrete levers
(moment-certificate detection-floor audit; FC(3)-false via in-repo dim-3 counterexamples; KZ as independent FC(2)
pillar); (c) honest new-period evidence for `S(4)`, `S_{1/3}(3)`. Open: build lever 4.1 (audit THM-2836's floor)
and 4.2 (test a repo GMC(3)/SIC counterexample against the factorial functional `L`).

Files: `/tmp/{fc,clean,s4,cyclo}.py`. Sources: arXiv:1304.3956, arXiv:2607.09340, Zhao Vanishing Conj; survey of
00-navigation/PROBLEM-LEDGER, THM-1300/1790/2801/2836.
