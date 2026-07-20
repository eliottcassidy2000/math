# The Hurwitz principle: JC₂'s golden corner, the LRC ghost-convergent law, and the g-bonacci kernels — one shape, three problems

**Instance:** mac-mini-2026-07-20-S137 (owner: work JC₂ creatively; pull from past threads
unexpectedly; think Fibonacci). Data: HYP-8155, `jc2_golden_corner_macmini_S137.out`.

## 1. The principle

> **Obstruction objects live at worst-approximable data.**

Hurwitz's theorem names the extremal: the golden ratio is the number the rationals
approximate worst. Wherever a problem's reduction/attack calculus is driven by rational
approximation (continued-fraction convergents, Euclid steps, Farey mediants), the object
that survives longest sits at the noble/golden end. Three instances now live in this repo:

1. **JC₂ (this session).** The two-variable Jacobian conjecture's reduction theory acts on
   candidate degree pairs like the Euclidean algorithm (each elementary automorphism eats
   the polygon along a convergent). Lamé: Euclid is slowest on consecutive Fibonacci pairs.
   Census under the conservative citable battery [Magnus gcd=1; divisibility reduction;
   Moh ≤ 100]: among all 14,661 surviving pairs ≤ 300, the maximal Euclid-proxy chain
   (length 9) is achieved by EXACTLY ONE pair — **(178, 288) = 2·(89, 144)**, slope
   144/89, the best φ-approximant in range. The reduction-resistant corner of the JC₂
   landscape is golden, on the nose. (Scope: the chain statistic is the Euclid proxy for
   the true Abhyankar–Moh polygon calculus — a frame, not a theorem about JC₂.)

2. **LRC(14) (opus-S407's ghost loop; THM-1291/1292).** The extremal family FORBIDS ITS OWN
   PENULTIMATE CONVERGENT — its excluded element is the best rational approximation of its
   maximizer, and D is the CF remainder of the convergent it keeps. Extremality =
   resisting your best approximant. Same shape, measured side.

3. **The (r,g) shadow lattice (klein-S313).** The g-bonacci kernels 1 − x − x^{g+1} govern
   every gap-g shadow (Fibonacci at g = 1); the missing-region law's deficit-1 is the
   kernel's boundary effect. The repo's native Fibonacci machinery, generating-function side.

The JC₃ counterexample fits the principle from the other direction: it fell where NO
approximation calculus constrained it (the Campbell-minimal S₃ stratum, escape at
infinity), while JC₂ — the survivor — is exactly the case where the Jung–van der Kulk
amalgam gives a full reduction calculus. **Conjecture-shaped reading: JC₂ survives because
in two variables every candidate is inside the reach of the convergent calculus; the golden
corner (g·Fibonacci degrees, φ-slope polygons) is where a counterexample would have to
live, and where the calculus is slowest but never absent.** A proof of JC₂ along this line
= a quantitative "the calculus always catches up" statement — a Lamé-style bound that the
polygon reduction terminates even at noble slopes. That is precisely a Zaremba/Hurwitz-type
uniform-CF problem, which is the LRC corpus's home turf — the unexpected pull the owner
asked for: **the repo's Farey/three-gap/convergent technology is native equipment for the
golden corner of JC₂.**

## 2. The engine-dimension lemma (why n = 3 fell first)

THM-1340's lens at n = 2: for y-affine F = (a(x)y+b(x), c(x)y+d(x)), the det's
y-coefficient is a′c − ac′ ≡ 0 ⟹ a = κc ⟹ triangular ⟹ invertible (verified symbolically;
likely classical). The "engine" P(a : c) is a POINT of ℙ¹ — no room to curve. The first
curved engine (the cuspidal cubic that powers the ℂ³ counterexample) needs ℙ², i.e. three
variables. Dimension 2 is the engine-starved case: another face of "JC₂ is different."

## 3. Filed follow-ups

- (i) Formalize the golden-corner statement inside the true polygon calculus (Abhyankar–Moh
  edges, not the Euclid proxy); the target lemma: reduction length ≤ C·log(min degree) with
  the constant achieved exactly on g·Fibonacci pairs — Lamé for polygons.
- (ii) The Zaremba bridge: bounded-partial-quotient slopes vs reduction reach — does the
  repo's CF-law machinery (opus-S403's D = convergent remainder) transfer a statement?
- (iii) Complete the exclusion battery with the full literature list (Appelgate–Onishi,
  Nagata, Heitmann gcd results — not implemented; census is conservative) and recompute the
  golden survivors above 100.

**Cross-links:** HYP-8155, THM-1340 (engine lens), opus-S407/THM-1291/1292 (ghost
convergents), opus-S403 (CF active-leg law), klein-S313 (g-bonacci), THM-1330 §4 (JC₂ as
the open floor), Jung–van der Kulk, Magnus, Moh, Lamé, Hurwitz.
