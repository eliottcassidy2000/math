# The Lee-Yang circle web: coverage = R⁶, the dip = the off-circle φ⁴, and the compression hierarchy

*mac-mini-2026-06-27-S72. Owner: merge a dense web — Lee-Yang circle polynomial = binomial coefficients =
"cap = pair-normalized Pascal mass" + de Moivre-Laplace, dip = the φ⁴ off-circle correction; full extremality =
circular zeros (low λ) + large radius R, `q₀ = q₆·R⁶`; ear decompositions (SC ⟺ ear, factor-critical ⟺ odd ear)
= the odd-cycle/Ω face; the k=8 dual `L_y = q₀+q₆+(1/10)q₃` = the bimodality functional; the flip-action is a
transformation MONOID (not V₄) — solvable via dual degree ≤4 ⟹ Galois ⊆ S₄; Newton/Maclaurin quartic moment
inequality extremal at the AP; and compression failures beyond commutativity. The threads form one coherent
object. Continues [[invariants-as-arc-cube-functions-the-compression-and-the-parity-split-of-k8]],
[[the-cap-is-a-phi4-field-theory-and-the-quartic-marks-the-hard-row]].*

## The Lee-Yang circle: coverage is a radius (VERIFIED)
The miss-count PGF `G_N(z) = Σ q_t z^t` (degree 6) has, by Vieta, `q₀ = q₆·∏|z_i|`, and its zeros lie **near a
circle of radius `R = (q₀/q₆)^{1/6}`** (`lrc_leeyang_circle_galois_newton_macmini_S72.py`):
```
 k:   8      9      10     11     12     13
 R:  1.59   1.69   1.78   1.85   1.91   1.96      (zero-radius ratio max/min = 1.14 .. 1.36)
 q₀ = q₆·R⁶ EXACTLY ;  coverage p0 = q₀ = q₆·R⁶.
```
So **coverage IS a radius**: `p0 = q₆·R⁶`, and the extremality has two knobs — **R (the radius = the coverage
scale)** and **λ (the off-circle deviation = how far the zeros leave the circle)**:
- **λ → 0 (zeros on the circle):** the **Lee-Yang circle** case = the **binomial** (Pascal/pair-Pascal)
  coefficients = the **de Moivre-Laplace Gaussian** limit. This is the **EVEN/symmetric face** (`a+b, a·b`):
  the `cap = C(k+1,2)/91` pair-normalized Pascal mass (THM-577), solvable.
- **λ > 0 (off-circle):** the **φ⁴ correction** = the **dip** = the **ODD/antisymmetric face** (`a^b, b^a`):
  the Worpitzky/ear/odd-cycle content (S71/codex HYP-3147).

So the owner's whole sentence is one statement: **cap = binomial/de-Moivre circle (low λ, even); dip = φ⁴
off-circle (odd); extremality = large R + low λ.** consec is the (large R, low λ) extremizer.

## The bimodality functional and Newton-Maclaurin (VERIFIED)
The k=8 dual `L_y = 10q₀ + q₃ + 10q₆ = 10·(q₀+q₆ + 0.1 q₃)` is **literally the bimodality functional** — the
extreme mass `q₀+q₆` (the two modes / the two poles of the circle) plus a touch of the middle. **consec
maximizes `L_y`** (`L_y(consec) ≤ 10·cap`, HYP-3085). And this is the same object as everything else:
- Newton-Maclaurin: the normalized moment means `p_k = S_k/C(6,k)` have **all Newton defects `< 0`** at consec
  (`p_k² − p_{k-1}p_{k+1} = −0.048, −0.004, …`) — **consec is the EXTREMAL of the moment-inequality VIOLATION**.
- max Newton-violation = max bimodality = max extreme-mass = **0 real roots** (S66) = the (large R, low λ)
  circular-zero extremizer. **One object, five names.** The "Newton/Maclaurin quartic moment inequality
  extremal at the AP" is the AP maximizing the *violation* (it is anti-log-concave / bimodal, not log-concave).

## The rigorous solvability: Galois ⊆ S₄ (correcting V₄)
The owner's correction: **the flip-action is a transformation MONOID, not the group `V₄`** — the apex arc `c`
is **absorbing** (homogenizes `T,+,− → S` and swaps `T↔S`), so there are no inverses. My S70/S71 "Klein-four"
was the 2-arc *slice*, not the full action. The rigorous solvable connection is instead:
> **dual degree ≤ 4 ⟹ Galois group ⊆ S₄ ⟹ solvable by radicals** (`S₄ ▷ A₄ ▷ V₄ ▷ 1`).
And the specific gK8 duals are even simpler: `(t−1)(t−2)(t−4)(t−5)` and `(t−2)(t−3)(t−6)` have **rational
roots** ⟹ **trivial Galois** (split over ℚ) — *why `cap`, `dip` are exact rationals* (THM-577). The solvability
is real but lives in "degree ≤ 4 ⟹ S₄-solvable", not in a `V₄` group law.

## Compression failures beyond commutativity (the hierarchy)
S71's "score determines iso (n≤4), fails (n=5)" is the **commutativity/linearity** rung. The owner's "monoid"
correction reveals the full hierarchy of what the iso-class-function (on the arc-cube) does/doesn't respect:
| property | structure | status | meaning |
|---|---|---|---|
| **linearity** (`a+b`, score, cut) | additive group | **exact n≤4, fails n=5** | the score determines the class up to n=4 |
| **associativity** (flip-composition = XOR) | a **monoid** | **always holds** | flips compose; the natural frame |
| **invertibility** (`V₄` / group) | a group | **FAILS** | the **absorbing apex arc** collapses `T,+,− → S` (info loss) |
| **multiplicativity** (`a·b`, H, cycle) | the OCF | **the irreducible nonlinear part** | the odd/off-circle/dip content |
The pattern: the **compressible** part is the linear/commutative/even (score/cut/binomial/biquadratic), exact
while degree ≤4; the **incompressible** part is the absorbing-apex / multiplicative / odd-cycle (the φ⁴
off-circle dip). **The "absorbing apex arc" that collapses information is the tournament avatar of the apex
prime 7** — the element whose flip homogenizes, the source of the irreducible (odd, off-circle) content.

## Ear = odd-cycle/Ω = the odd/off-circle face
SC ⟺ ear decomposition; **factor-critical ⟺ ODD ear decomposition** — the odd ears are the odd cycles, i.e.
`Ω` / `H = I(Ω,2)`, i.e. the **odd cumulant κ₃** = the **off-circle φ⁴** = the **Worpitzky** content. So the ear
decomposition is the combinatorial certificate of the off-circle (incompressible) part, exactly the antisymmetric
face that the score-compression and the circle (even) face cannot see.

## The unified web (toward the proof)
```
 EXTREMALITY = (large R, low λ) on the Lee-Yang circle, coverage p0 = q₆·R⁶
   EVEN / on-circle (λ=0):   binomial = pair-Pascal cap = de Moivre-Laplace = biquadratic = Galois-solvable
   ODD  / off-circle (λ>0):  the dip = φ⁴ = odd cumulant κ₃ = odd ear = Ω = Worpitzky = the absorbing apex
 consec extremizes ALL: max R, max bimodality L_y, max Newton-violation, 0 real roots.
 PROOF = bound the dip = [even biquadratic, solvable degree-2, S70] + [odd Worpitzky/ear/κ₃, codex HYP-3147],
         with Galois⊆S₄ guaranteeing the even part is explicit and the Lee-Yang circle keeping λ small.
```

## Honest status
VERIFIED: `q₀=q₆R⁶`; zeros near a circle of radius R (coverage=R⁶); Galois trivial/⊆S₄ (solvable); the
bimodality functional; Newton defects all negative at consec (the violation extremal). CORRECTED: the flip-action
is a monoid (absorbing apex), not V₄; solvability is via degree≤4 ⟹ S₄, not a group law. SYNTHESIS/BOLD: the
λ (off-circle) ↔ dip ↔ odd-cumulant ↔ ear ↔ Worpitzky identification as a *single* incompressible quantity; that
bounding the dip = bounding λ = the even(biquadratic)+odd(Worpitzky) split. The reusable lesson: **coverage is a
radius, the dip is the off-circle, and the off-circle is the apex's irreducible odd content.**

Related: HYP-3152 (this), HYP-3150 (parity split), HYP-3147 (Worpitzky/ear = odd face), HYP-3132 (biquadratic =
even face), HYP-3122 (φ⁴), HYP-3103 (PGF zeros = the circle), THM-577 (cap = binomial), OPEN-Q-108.
