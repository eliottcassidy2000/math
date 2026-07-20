# An atlas of lenses: how many famous problems the GMC/LRC machinery links

**death-star-2026-07-20-S68** (HYP-8520). Owner: see how many famous problems these lenses link,
or define new ones. This is a **navigational map**, not a claim of reductions: each lens is a
structural pattern; the famous-problem placements are established facts (Weil positivity,
Marcus–Spielman–Srivastava, Sendov–Tao, Littlewood–Konyagin, …); the *links to our threads* are the
contribution. Two cross-links are machine-anchored below. Strength is flagged **[textbook] /
[strong] / [speculative]**. Builds on the S67 GMC↔LRC bridge.

---

## Lens 1 — Positivity past the Cancellation Wall (the "circle-method skeleton")

**Pattern.** The target = a **positive** main term (singular series / positive-definite form /
autocorrelation energy) minus a **signed** error; the goal is to show it is nonzero *without*
absolute-bounding the error, which over-counts or diverges. **Escape: reformulate positive-definite.**

| problem | positive term | signed wall | strength |
|---|---|---|---|
| **Riemann Hypothesis** | Weil explicit formula as a PD quadratic form; Li's criterion `λ_n≥0` | prime oscillation `Σ_ρ x^ρ/ρ` | [textbook] RH ⟺ positivity |
| **Goldbach / twin primes / Waring** (circle method) | singular series `𝔖(N)>0` (local densities) | minor-arc cancellation | [textbook] |
| **Erdős discrepancy** (Tao 2015) | logarithmic-Chowla / entropy positivity | the signed discrepancy sum | [strong] |
| **Chowla / Sarnak** (Möbius) | Gowers-norm / entropy PD methods | `Σ μ(n)e(nθ)` cancellation | [strong] |
| **LRC(14)** (ours) | autocorrelation `disc_v` (klein-S287) | covering `corrsum`, S266 divergence | our work |
| **GMC(2)** (ours) | Hankel-PD `α=r|c|²≥0` (klein-S363) | sign-indefinite `β`-coupling; domination dead (MISTAKE-202) | our work |

**Transfer.** Our `MISTAKE-202` ("never absolute-bound the signed part; square it") is a *local
instance of the universal law* that already governs RH (Weil), the circle method (singular series),
and Erdős discrepancy (entropy). The reflection: these are one method, discovered independently many
times. When *our* signed estimate resists, the historical answer is Parseval/PD, not a sharper bound.

---

## Lens 2 — L^p-energy of a Structured Exponential Sum (the "jump-sum")  ⟨anchored⟩

**Pattern.** An exponential sum `S(θ)=Σ_{x∈A} c_x e(θx)` over structured support `A`; the object is an
`L^1/L^2/L^4/L^∞` norm or an **autocorrelation energy** `Σ_k|C_k|²` (= a Parseval energy of `|S|²`).

- **Littlewood's L^1 problem** — `‖S‖_1 ≥ c·log|A|` (Konyagin; McGehee–Pigno–Smith). [textbook]
- **Merit factor / Golay / flat polynomials** — minimize the L^4 autocorrelation energy of `±1`
  sequences. [textbook]
- **LRC(14) `disc_v`** and **GMC reconstruction jumps** — ours (S67): `disc_v = Σ|S_{mv}|²/(2πmv)²`,
  `S_m=Σ_j sign_j e(mx_j)` over arc endpoints — *the same jump-sum*.

**Anchored** (`04-computation/lens_atlas_crosslinks_deathstar_S68.py`): the merit-factor autocorrelation
energy equals its Parseval form `(∫|S|⁴−N²)/2` (verified `3=3.000`), and this is the same functional
shape as `disc_v`. So **merit factor, Littlewood, LRC `disc_v`, GMC reconstruction are one object** —
an `L^p`-energy of a structured spectrum. **Transfer:** the three-gap arc endpoints give the LRC
energy in closed form (S67 backlog lead); the merit-factor literature's autocorrelation bounds are
directly importable.

---

## Lens 3 — Folds, Critical Points, and the Loop (the "geometry of roots")

**Pattern.** The critical points / folds of a function and their non-degeneracy / non-cancellation,
controlled by **topological (loop / winding) or interlacing** invariants rather than estimates.

- **Sendov's conjecture** (Tao, high degree) — every root has a critical point within distance 1.
  [strong]
- **Smale's mean value conjecture** — critical values vs the polynomial. [textbook open]
- **Gauss–Lucas / Ilieff–Sendov** — critical points inside the root hull. [textbook]
- **GMC(2) pinch/loop** (ours, S66) — the arc between coincident folds is a *loop* ⟹ opposite `g'`
  signs ⟹ real ⊥ imaginary ⟹ non-cancellation.

**Transfer.** My S66 loop argument (a self-intersection forces opposite derivative signs, so the two
folds cannot cancel) is a topological cousin of the root/critical-point geometry in Sendov and
Gauss–Lucas — a *proof shape* for "a degenerate stratum survives by winding, not by an estimate."

---

## Lens 4 — No Common Root / Real-Rootedness / Interlacing (moment rigidity)  ⟨anchored⟩

**Pattern.** A moment sequence whose orthogonal polynomials are **real-rooted and interlace** (share
no root), forbidding cancellation of a perturbation.

- **Kadison–Singer** (Marcus–Spielman–Srivastava 2013) — interlacing families / the barrier method.
  [textbook, resolved]
- **Lee–Yang / de Bruijn–Newman** (`Λ≥0`, Rodgers–Tao 2018) — RH as the survival of real-rootedness
  of `ξ`-approximants under heat flow. [textbook]
- **Hamburger/Stieltjes moment determinacy**. [textbook]
- **GMC(2) Hermite/Legendre** (ours) — the constant-coefficient closure is exactly "consecutive
  Hermite share no root."

**Anchored**: consecutive probabilists' Hermite polynomials interlace and share no root (verified,
`He_3..He_5`). So **the GMC nullcone "no common root," Lee–Yang real-rootedness (RH), and MSS
interlacing (Kadison–Singer) are one rigidity.** **Transfer:** the GMC failure mode — the recognition
is exact but the *recurrence opens a hierarchy* for non-constant coefficients (S64) — is the same
obstruction MSS overcame with interlacing *families* (not a single sequence). The MSS
"expected-characteristic-polynomial + barrier" is a candidate for the GMC `β`-coupling hierarchy.

---

## Lens 5 — Charge Grading / the Killed Kernel (this is our *core* project too)

**Pattern.** A functional (expectation, constant term, trace) **kills a graded piece** (a charge, a
Fourier mode, an odd-cycle parity); the kernel is conjecturally the "one-sided / trivial" part;
Mathieu–Zhao spaces formalize it.

- **Jacobian Conjecture** — via GMC / Mathieu–Zhao spaces (`GMC(n)⟹JC(n)`). [strong, our thread]
- **Mathieu conjecture / Duistermaat–van der Kallen** (tori) — the toral kernel is one-signed.
  [textbook, = our TNC]
- **Dixmier conjecture** (Weyl algebra endomorphisms). [textbook open]
- **Rédei's theorem / the Odd-Cycle Collection Formula** (the project's *heart*) — the tournament
  "charge" is **odd-cycle parity**, and the parity functional kills nonzero charge. Rédei (odd number
  of Hamiltonian paths) and the OCR are **charge-grading statements** in the same family as GMC's
  `E[Z^aW̄^b]=a!δ_{ab}` (kills nonzero Fourier charge). [our core]

**Transfer, the striking one.** The GMC "charge = Fourier mode of `arg Z`, `E` kills nonzero charge"
and the tournament "charge = odd-cycle parity, the OCR kills nonzero parity" are the **same U(1)/Z₂
grading-killed-by-a-functional structure**. The project's flagship object (Rédei/OCR) and the GMC
bridge to the Jacobian Conjecture sit in *one lens* — worth mining: is the tournament parity functional
a Mathieu–Zhao space? (That would tie Rédei directly to JC.)

---

## Lens 6 (meta) — Representation vs Function level (which proofs survive)

boxeph's MISTAKE-203: **representation-level** arguments (domination, the pinch-*sweep*) dodge and die
under review; **function-level** invariants (the Radial Lemma, the orbit-product TNC, my loop) survive.
This meta-lens applies across all the above: RH's failed "elementary" attacks were representation-level;
what advanced it (Weil positivity, Rodgers–Tao heat-flow) is function-level. **Rule for every problem
here:** exhibit an invariant separating the two sides at the *function* level, not the term level.

---

## The Rosetta (which problem sits in which lens)

| | L1 Positivity | L2 Jump-sum | L3 Fold/loop | L4 No-common-root | L5 Charge |
|---|:-:|:-:|:-:|:-:|:-:|
| Riemann Hypothesis | ● | | | ● (Lee–Yang) | |
| Goldbach / Waring | ● | ● (minor arcs) | | | |
| Erdős discrepancy | ● | ● (autocorr) | | | |
| Littlewood / merit factor | | ● | | | |
| Sendov / Smale | | | ● | ● | |
| Kadison–Singer | | | | ● | |
| Jacobian / Mathieu / Dixmier | | | | | ● |
| **Rédei / OCR (core)** | | | | | ● |
| **LRC(14) (ours)** | ● | ● | ○ (AP=resonance) | | |
| **GMC(2) (ours)** | ● | ● | ● (pinch/loop) | ● (Hermite) | ● |

**GMC(2) touches every lens** — which is *why* it has been such a productive testbed, and why its
lessons transfer so widely. **LRC and the core Rédei project each share ≥2 lenses with it.**

## Honest scope
A map + two anchored identities, not proofs or reductions. Its value is navigational: it says *which
famous machinery to import* for each of our open items (Weil/entropy positivity → LRC/GMC signed parts;
merit-factor bounds → `disc_v`; MSS interlacing → the GMC `β`-hierarchy; Mathieu–Zhao → is Rédei one?).
The single most actionable new question it raises: **is the tournament odd-cycle-parity functional a
Mathieu–Zhao space?** — which would place Rédei's theorem in the Jacobian ecosystem directly.

## Credit
death-star-S67 (the GMC↔LRC bridge this generalizes), klein-S363/S287 (the positivity move on both),
boxeph-S181/S182 (reconstruction, MISTAKE-203 meta-lesson), the project's Rédei/OCR canon; external:
Weil, Li, Konyagin, Marcus–Spielman–Srivastava, Rodgers–Tao, Tao (Sendov, discrepancy), van den Essen
/ Mathieu / Duistermaat–van der Kallen.

## Cross-links
S67 bridge · PROBLEM-LEDGER.md (JC ecosystem) · MISTAKE-202/203 · `04-computation/lens_atlas_crosslinks_deathstar_S68.py` · HYP-8520.
