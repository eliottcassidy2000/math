---
source: kind-pasteur-2026-07-24-S149 (Opus 4.8)
status: RESULT + correction. Pins the closed-form locus C_lambda: only lambda=1/4 reaches k=3
  (C_{1/4}={1,2,3}, others ={1,2}) -- strong PSLQ evidence, and it REFUTES my kps-S148 "CM-order resonance"
  heuristic. Identifies the genuine governor of the closed forms as the Conway-Jones / vanishing-sums-of-
  roots-of-unity theory (in-repo THM-415), and honestly disposes of the look-and-say/cosmology idea (no link).
tags: [hypergeometric, closed-forms, roots-of-unity, conway-jones, CM, correction, series]
related: [kps-S146, kps-S147, kps-S148, THM-415, HYP-2273]
corrects: [kps-S148 resonance heuristic]
---

# Pinning the closed-form locus, and the roots-of-unity governor

## 1. The pin (strong evidence)
`S_lambda(k) = int_0^1 2F1(lambda,1-lambda;1;x^k) dx`, `C_lambda = {k : S_lambda(k) elementary}`.
- `S_{1/4}(3)` elementary (owner-given; PSLQ `[1,1,2]` reconfirmed).
- `S_{1/3}(3)` and `S_{1/6}(3)`: **NO relation** vs a comprehensive mixed-field basis
  (`sqrt3 log(2+sqrt3)`, `log(2+sqrt3)`, `sqrt3 log2`, `log2`, `sqrt2 log(1+sqrt2)`, `sqrt3 log(5-2sqrt6)`,
  `log(8-3sqrt7)`, `atan(sqrt3)`, `atan(1/sqrt3)`, `atan(sqrt2/5)`, `atan(sqrt3/5)`, `atan(sqrt2)`,
  `atan(1/(2+sqrt3))`, `pi`, `1`), 130 digits, coeff bound 2000. Strong evidence **not elementary**.

> **`C_{1/4} = {1,2,3}`, while `C_{1/2}=C_{1/3}=C_{1/6}={1,2}` (through the tested range).**
> Only the lemniscatic signature `lambda=1/4` reaches `k=3`.

## 2. CORRECTION to kps-S148's resonance heuristic
kps-S148 proposed "`S_lambda(k)` elementary iff `k` resonates with the CM order of `E_lambda`", predicting the
equianharmonic `lambda=1/3,1/6` (order 3/6) would be the ones elementary at `k=3`. **The data says the OPPOSITE:**
`lambda=1/4` (order 4, `Z[i]`) reaches `k=3`; the order-3/6 signatures do NOT. So the naive CM-order-vs-level
resonance is **refuted**. The correct mechanism for why `lambda=1/4` is special at `k=3` is open (section 4).

## 3. The genuine governor: Conway-Jones (vanishing sums of roots of unity)
The closed form of `S_lambda(k)` comes from partial-fractioning the `k`-power (`1/(kn+1)=int_0^1 x^{kn}`),
i.e. from the `k`-th roots of unity; the resulting log/arctan arguments are algebraic combinations of period
values at those roots. Whether they collapse to an ELEMENTARY expression is exactly a question of **algebraic
relations among roots-of-unity values** -- the **Conway-Jones / Lam-Leung** theory of *vanishing sums of roots
of unity*. This is not a metaphor: it is already the governing dichotomy in the repo's LRC work --
> **THM-415:** "the prime/composite dichotomy ... IS the vanishing-sums-of-roots-of-unity dichotomy
> (Lam-Leung / Conway-Jones): a sine sum over a half-system can vanish nontrivially exactly when the modulus
> is composite."
So the SAME Conway-Jones roots-of-unity theory governs (a) the LRC signed-orbit collisions and (b) the
closed-form locus `C_lambda` of the series. That is the real "Conway connection." The relevant `k` for the
series is a *level* (roots of unity of order `k`), and the composite/prime structure of `k` -- exactly the
Conway-Jones regime -- is the natural place to look for why `k=3` behaves as it does per signature.

## 4. Look-and-say / Conway cosmology: no link (honest)
The repo has **no** look-and-say / audioactive / cosmological-theorem thread (grep empty). Mathematically,
Conway's look-and-say cosmology is a statement about an integer dynamical system (audioactive decay into 92
elements; the growth constant is the largest root of a degree-71 polynomial). I find **no genuine connection**
to hypergeometric periods, CM, or `C_lambda`: the series' algebraic numbers are low-degree (`sqrt2,sqrt3,sqrt6`),
nothing degree-71 or transfer-matrix-like appears. The only honest link is a loose **metaphor**: "the closed-form
locus `C_lambda` is finite" resembles "sequences decay into a finite periodic table" -- suggestive of a
classification, not a theorem. The user's "Conway" instinct is right, but it points to **Conway-Jones**
(roots of unity), not look-and-say.

## 5. Open questions this pins down
1. **Prove `C_{1/4} = {1,2,3}` exactly** (is `S_{1/4}(k)` non-elementary for all `k>=4`?) and confirm
   `3 not in C_{1/3}, C_{1/6}` at higher precision / larger fields (rule out a missed degree-4 unit).
2. **Why is `lambda=1/4` special at `k=3`?** Not CM-order resonance. Candidate: a Conway-Jones relation among
   cube-roots-of-unity values of the signature-4 period that fails for signatures 3/6. This is the crisp,
   concrete form of the "closed-form locus" problem.
3. **Transfer to LRC:** since Conway-Jones governs both, a vanishing-sum-of-roots-of-unity audit at the LRC
   drift levels (composite `k` in the census) is the shared tool -- ties back to kps-S145/S148.

Files: `/tmp/{pin2,pin3}.py`. Corrects kps-S148 section 3 (resonance).
