# How the LRC route sharpened — and the recursion it reveals: the magic-function kernel is a three-gap/Stern-Brocot function of a/b (the Farey thread IS the cap recursion)

*mac-mini-2026-06-28-S75b. The owner: gain a comprehensive understanding of exactly how our route has gradually
sharpened, and use it to see recursive patterns not yet precisely described. Tracing the sharpening exposes a
self-similar route AND a precisely-describable recursion in the cap kernel — the three-gap/Stern-Brocot structure,
which is the Farey thread re-entering as the covering-bound recursion. Builds on
[[building-the-magic-function-it-is-the-comb-overlap-gram-kernel-bochner-automatic]],
[[the-one-obstruction-worst-case-algebra-vs-analytic-equidistribution-apex7-is-the-lee-yang-zero]], HYP-2913, THM-576.*

## I. How the route sharpened (the comprehensive picture)
Each step replaced a coarse object with a finer, more-structured one — and the SAME apex-7 / three-gap structure
reappears at every level (the route is **self-similar**: it sharpens by finding the same skeleton at finer scale):

| step | object | what got sharper | session |
|---|---|---|---|
| 0 | LRC(14): `M(S) ≥ 1/14` ∀ 13-set | the whole problem | — |
| 1 | **covering bound**: only sets with a `14|s` are hard | the SET CLASS | S59 |
| 2 | **the cap** `= min meas(lonely(P))` | a problem → a NUMBER | THM-534/576 |
| 3 | `cap_k = C(k+1,2)/91` (pairwise-avoidance) + 2 binding constants `k=8,9` | the cap FORMULA + the hard rows | THM-576/577 |
| 4 | **the magic function**: Fejér kernel `F_7` (Fourier) = comb-overlap Gram kernel (spatial) | the CERTIFICATE | S75, kps S31ao |
| 5 | **single-arc lemma + peeling recursion**: `cap(P)=cap(P∖1)−(1/7)(1−1/min)` | a RECURSION (peel speed 1) | S75 |
| 6 | **the kernel's three-gap recursion** (below) | the FULL arithmetic recursion | S75b (this) |

The sharpening is monotone: **problem → set-class → number → formula → kernel → recursion → arithmetic recursion.**
And at each level the apex-7 recurs: 7 sectors (step 2), `/91 = /C(14,2)` and `1/(7q)` (steps 3–4), `F_7(0)=49=7²`
(step 4), the single-arc fails exactly at the apex `14` (step 5), and now `K(7,b)=1/49` (step 6). The route is a
**zoom into one self-similar object** — the apex-7 three-gap skeleton.

## II. The recursion not yet precisely described: the kernel is a THREE-GAP/Stern-Brocot function
The magic function is the comb-overlap kernel `K(p,q)=meas(D_p∩D_q)`. Its precise recursive structure
(`lrc_magic_function_resonance_domain` + the kernel scans):
1. **Scale-invariance:** `K(p,q) = K(p/gcd, q/gcd)` (VERIFIED all `p,q≤13`) — the kernel depends only on the
   coprime ratio `a/b`. This is the first recursion (Euclidean common-factor reduction).
2. **The gap-count is piecewise-linear in `a` (the three-gap signature):** writing `K(a,b)=g(a,b)/(7ab)`, the
   integer `g(a,b)=7ab·K(a,b)` is a **piecewise-linear function of `a` with breakpoints at the continued-fraction
   convergents of `a/b`** — exactly the three-gap (Steinhaus) theorem. Evidence:
   - `b=13` (antipode `≡ −1 mod 14`): `g = 1,3,5,7,…,23 = (2a−1)` — single regime, so **`K(a,13)=(2a−1)/(91a)`**
     (VERIFIED all `a=1..12`).
   - `b=11`: `g = 1,2,3, 5,7,9,11,13,15,17` — two regimes (slope-1 then slope-2), breakpoint at `a=3` (the
     convergent of `11`). The kink IS the three-gap transition.
3. **Base case** `K(1,b)=1/(7b)` (the single-arc lemma, PROVED S75); **apex column** `K(7,b)=1/49` (the prime-7
   comb tiles the 7-sectors). These are the boundary conditions of the recursion.

So: **the cap kernel `K(a,b)` is computed by walking the continued fraction of `a/b` (the Stern-Brocot tree); the
three-gap theorem gives the piecewise-linear gap-count `g(a,b)`; `K = g/(7ab)`.** This is the recursive pattern
the route had been zooming toward without naming.

## III. The synthesis: the Farey thread IS the cap recursion (the two halves merge)
At S59 I split LRC(14) into the **additive face** (census / three-gap / Farey / Stern-Brocot — called "the easy
case's extremal") and the **covering face** (the cap / the proof). This recursion **merges them**: the cap kernel
`K(a,b)` is governed by the **same three-gap/Stern-Brocot/Farey recursion** that controls the census. The Farey
thread (HYP-2913) is not the easy shadow — it is the **arithmetic recursion of the covering-bound magic function**.
The route sharpened until the "easy" additive structure (Farey) turned out to BE the recursion of the "hard"
covering cap. One object, seen as census at the coarse scale and as the cap-kernel recursion at the fine scale.

## Honest status
- **VERIFIED:** scale-invariance `K(p,q)=K(p/d,q/d)`; the antipode formula `K(a,13)=(2a−1)/(91a)` (all `a`); the
  piecewise-linear (three-gap) gap-count `g(a,b)=7ab·K(a,b)` with convergent breakpoints (`b=11,13`); base
  `K(1,b)=1/(7b)`, apex `K(7,b)=1/49`.
- **DESCRIBED (the recursion):** `K(a,b) = g(a,b)/(7ab)`, `g` the three-gap/Stern-Brocot gap-count of `a/b`. This
  precisely names the recursion the route was sharpening toward; it is the Farey/three-gap thread realized as the
  cap-kernel recursion.
- **NOT a proof.** The kernel recursion is the order-2 (pairwise) building block of the cap; the cap is the
  inclusion-exclusion over these kernels, and the binding rows still need the order-3+ terms. But the recursive
  STRUCTURE of the magic function is now precise — the path to closing the cap runs through the continued-fraction
  (three-gap) recursion of the kernel, not a config-blind certificate. LRC(14) open.

Related: HYP-3230 (this), HYP-3227 (the Gram kernel + single-arc lemma), HYP-3214 (kps Fejér), THM-576 (cap=pairwise),
HYP-2913 (three-gap/Farey census), the four-Farey-variations (S59), OPEN-Q-108.
