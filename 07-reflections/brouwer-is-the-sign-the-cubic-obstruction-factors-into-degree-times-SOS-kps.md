# Brouwer is the SIGN: the non-SOS cubic obstruction factors into (topological degree) × (SOS magnitude)

*kind-pasteur-2026-06-28-S31at. The owner asked to merge in the Brouwer fixed-point theorem, think
abstractly about what it and the other concepts represent, and keep working toward a proof. Brouwer is the
prototype of a **topological** certificate — existence/degree/sign — for facts that **algebra cannot see**.
The project's "one obstruction" (mac-mini HYP-3221) is precisely a fact algebra cannot see: the **non-SOS**
odd/cubic half, where every algebraic certificate (SOS, Lee-Yang, moment-LP, Bonferroni) stalls. This
reflection merges Brouwer in by showing the cubic obstruction **factors into a sign (Brouwer/degree) times
an SOS magnitude (algebra)** — so it is certifiable after all, by topology × algebra, and `14 = 2·7` is
exactly that product.*

## What each thing represents (the abstract frame)
| concept | represents | LRC role |
|---|---|---|
| **Brouwer / degree theory** | **existence via topology; the SIGN/orientation algebra can't compute** | the trace-sign of the obstruction; the witness as a minimax equilibrium |
| **SOS / Positivstellensatz** | algebraic positivity (the MAGNITUDE) | mac-mini S75e: the Fejér = cyclotomic square, even half SOS |
| **Equidistribution (Erdős–Turán)** | the analytic fixing of the leading term | mac-mini S74: what rescues the apex comb |
| **Cyclotomic ℚ(ζ₇)** | the arithmetic of the fixed point | the value (de Moivre/Fejér), HYP-3216/3217 |
| **Fixed point (Brouwer)** | the AP / the de Moivre–Fejér RG attractor | HYP-3218 self-dual config |

## The computation: the cubic obstruction FACTORS (verified)
mac-mini named the obstruction "non-SOS (Motzkin / AM-GM-of-3)." On the **de Moivre periods**
`a,b,c = 2cos(2πj/7)` (the cubic Gaussian periods, roots of `x³+x²−2x−1`, HYP-3217), the AM-GM-of-3 form is
```
a³+b³+c³ − 3abc  =  (a+b+c)·(a²+b²+c²−ab−bc−ca)  =  (a+b+c)·½[(a−b)²+(b−c)²+(c−a)²]
                 =  (e₁) · (e₁²−3e₂)  =  (−1)·(7)  =  −7 .
```
- the **first factor `e₁ = a+b+c = −1`** is the **trace** (the Newton trace, the `x²`-coefficient) — a pure
  **sign**, `±1`;
- the **second factor `½Σ(a−b)² = 7`** is a **sum of squares**, and it equals `|g(χ₇)|² = √disc = 7` — the
  Gauss-sum modulus / the apex prime (mac-mini's even-half SOS, S75e);
- so the AM-GM defect `= −7 = −(apex)`; and `disc = 49 = 7² = (the SOS factor)²`.

**The obstruction is non-SOS for one reason only: the trace sign `e₁=−1<0`.** Its magnitude is SOS times
the discriminant. A definite-sign × SOS form is exactly the kind a single SOS certificate misses (an SOS
form is `≥0` everywhere; this one is `≤0` because of the sign) — which is **why algebra stalled** — and
exactly the kind a **degree/sign argument (Brouwer)** certifies.

## The merge: `14 = 2·7` = (Brouwer SIGN, ℤ/2) × (SOS DISCRIMINANT, 7)
> **The cubic obstruction = (the trace sign, a Brouwer/degree datum in `ℤ/2 = {±1}`) × (an SOS magnitude =
> the discriminant 7).** The two factors are the two prime factors of `14`: the `2` is the **sign/orientation
> (Brouwer, the complement `ℤ/2`)**, the `7` is the **SOS magnitude (the cyclotomic discriminant)**.

This sharpens mac-mini's diagnosis. "No algebraic certificate" is right *for the sign* — SOS cannot produce
the negative trace. But the sign is a **single bit**, fixed by the cyclotomic (`e₁ = −1` is the Newton trace
of the de Moivre cubic, an integer), and **Brouwer/degree is exactly the tool that supplies a sign**. So:
```
obstruction certified  ⇐  Brouwer (the trace sign e₁ = −1, fixed) × SOS (the magnitude ½Σ(a−b)² = 7).
```
The de Moivre cubic **provides its own Reznick/Artin denominator** `½Σ(a−b)²` for free — the Motzkin-type
form is SOS-up-to-its-own-fixed-trace-sign. No external multiplier search needed.

## Where the three certificates meet (analysis = topology = the trace)
The trace `e₁ = a+b+c` is simultaneously:
- the **Brouwer/degree** datum (the sign/orientation/index);
- the **Frobenius / equidistribution** leading term (the `n=0` Weyl coefficient — mac-mini's analytic finish
  *fixes this sign*: the integer speeds' equidistribution makes the leading exponential sum the trace);
- the **Möbius (degree-1) mode** of HYP-3217 (the trivial-character / trace face).

So the three "rescues" the project has invoked — **analysis** (equidistribution), **topology** (Brouwer
degree), and the **degree-1 cyclotomic mode** (trace) — are the **same factor `e₁`**, and the **algebra**
(SOS) is the **other factor `7`**. The obstruction is their product. The fight (HYP-3218) dissolves again:
analysis and topology are not competing routes to the sign — they are the same sign.

## Honest status
- **VERIFIED:** the AM-GM-of-3 form on the de Moivre periods factors as `(e₁)(½Σ(a−b)²) = (−1)(7) = −7`; the
  second factor is SOS and equals `√disc = |g(χ₇)|² = 7`; `disc = 49 = 7²`.
- **MERGE (frame):** Brouwer = the sign/degree certificate for the non-SOS half; the obstruction = (Brouwer
  sign, fixed by the cyclotomic trace) × (SOS magnitude = discriminant); `14 = 2·7` = sign × magnitude. The
  trace `e₁` unifies analysis (equidistribution), topology (degree), and the degree-1 mode.
- **NOT a proof.** This is the cyclotomic *core* of the obstruction (at the de Moivre fixed point); the full
  lonely measure carries the same factorization only modulo the equidistribution tail (HYP-3218). Value: it
  tells the team to **stop trying to make the odd half SOS** (impossible — it has the wrong sign) and instead
  **certify its single fixed sign topologically (Brouwer/degree = the trace) and its magnitude by SOS** — a
  concrete two-factor certificate, not an analytic black box.

→ HYP-3219 (this), HYP-3221 (one obstruction = non-SOS), HYP-3235 (mac-mini: Fejér = cyclotomic square, even
half SOS), HYP-3217 (de Moivre = cubic mode), HYP-3218 (equidistribution = Fejér, the fight), Reznick/Artin
(SOS denominators), Brouwer/degree theory, OPEN-Q-108.
