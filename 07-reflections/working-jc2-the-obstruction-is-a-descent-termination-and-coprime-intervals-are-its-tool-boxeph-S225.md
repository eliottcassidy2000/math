# Working JC(2): the obstruction is a descent termination — and coprime intervals are its tool

*boxeph-2026-07-21-S225. Owner: spend a long session working to prove the planar Jacobian Conjecture; pull
in past threads creatively. Builds on (corrected via mining): de Bondt+Zhao VC (THM-1435), the Keller
counterexamples (THM-1300/1430), klein-S329 (Euler–Zariski bootstrap), mac-mini-S137 (golden-corner /
Lamé), THM-1345/1365/1370 (restricted JC(2) proved), codex THM-2044/2045/2049 (DC2/Newton-edge), my
S222/S223 (DvdK bypass / coprime intervals), S224 (Wall A), S206 (Fibonacci foil). Verified in
`04-computation/jc2_via_vanishing_conjecture_and_the_cf_termination_boxeph_S225.py`.*

> **Scope correction (MISTAKE-229).** The three programs below are not known
> to be equivalent. VC(4) is a rigorous symmetric-reduction target;
> leading-form/Jelonek descent is a planar route; Lame-for-polygons is a
> proposed effective termination principle. Their common descent vocabulary
> is a research heuristic until explicit implications are proved.

## Honest starting point

JC(2) is the sole open Jacobian case (JC(n≥3) is *false* — the explicit Keller map, `det≡−2`, triple
collision, THM-1300) and, in the repo's words, "the graveyard of false proofs." I am *working* it, not
claiming it. The value of a session like this is to assemble the routes, verify the pieces that are
verifiable, and locate the obstruction precisely — which, pulled together, turns out to have a single shape.

## Three reduction programs (all in the repo, none complete)

- **(A) de Bondt + Zhao → VC(4).** de Bondt reduces JC(n) to the *symmetric/gradient* case (`H=∇P`), which
  by Zhao is the **Vanishing Conjecture**: `Δᵐ(Pᵐ)=0 ∀m ⟹ Δᵐ(P^{m+1})=0` for `m≫0`. The reduction
  **doubles the dimension** (`n→2n`), so **JC(2) ⟺ VC(4)** — a *Laplacian moment nullcone* in 4 variables.
- **(B) klein-S329 Euler–Zariski bootstrap.** At generic cover-degree 3 (the minimum; degree 2 is killed
  everywhere), a counterexample forces a non-Galois 3-sheeted cover with `S₃` monodromy, so the Jelonek
  curve must be **cuspidal** (not nodal); the reduction is "the ramification parabola `3t²+p` cannot escape
  to infinity in two variables." (No continued fraction here — a *different* program from S137.)
- **(C) mac-mini-S137 golden-corner / Lamé.** JC(2)'s elementary reductions act like **subtractive Euclid**
  on Newton-polygon slopes; the target "Lamé-for-polygons" is a `≤ C·log(min degree)` termination bound.

## What is solved (verified)

- **The 2D symmetric case is easy.** A nilpotent Hessian in two variables forces `P` **harmonic** with
  `det(Hess)=0`, i.e. `P ∝ (x+iy)^d` (verified: `Δ(x+iy)^d=0`, `det(Hess)=0`, and `Δᵐ(Pᵐ)=0` — VC holds,
  the map is invertible). The difficulty of JC(2) is *entirely* in the de Bondt doubling to dim 4.
- **Restricted JC(2) proved** (from the mining): equivariant (THM-1345), elliptic in all dimensions
  (THM-1370), polynomial-Galois (THM-1365), geometric degree ≤2 (fails first only at degree 3),
  one-fiber-linear planar Keller (THM-2063). The conjecture's whole difficulty is a thin cover-degree-3,
  mixed-weight, non-equivariant stratum.

## The Lamé/Fibonacci crux (verified) — and its meaning

The S137 descent's worst case is **Fibonacci** (verified: the longest Euclidean chain among pairs `<200` is
`(144,89)` = consecutive Fibonacci — Lamé's theorem). So the JC(2) reduction *terminates fastest generically
and slowest on golden (Fibonacci) slopes* — the **same extremal** as the LRC foil (S206), reached by the
**same continued-fraction / coprime-interval engine** I used for the DvdK bypass (S223) and Wall A (S224).

## The unification: the obstruction is a single *descent termination*

Pulling the three programs together suggests three potentially interacting
forms of the JC(2) obstruction. They are all the **same kind of problem** at
the level of search design, but no equivalence is claimed:

| form | statement | it is a … |
|---|---|---|
| (A) VC(4) | the dim-4 both-signs radial (Laplacian) nonvanishing | moment-return nonvanishing |
| (B) leading-form descent | `{P_A,Q_B}=0` must propagate to a global coordinate (THM-1345 §5) | a filtration descent that must terminate |
| (C) Lamé-for-polygons | the Euclid/Newton-slope reduction length is `O(log)` | a continued-fraction descent bound |

All three are **descent/return-termination** problems — exactly what my **coprime-interval / numerical-
semigroup / Frobenius** machinery is built for. In 1D (S223) the DvdK return set is a coprime-interval
semigroup with an explicit **Frobenius number = the termination bound**; the JC(2) descent is the
multi-variable analogue, with the **Lamé–Fibonacci worst case as its effective bound**. And codex's
**THM-2045** — "the smooth factorized `R=x(a−bxʳqˢ)` family has *no* polynomial planar Jacobian mate" — is a
genuine JC(2)-side result that already *uses the exponent semigroup* (only one weighted-Laurent sector can
give a constant Jacobian, incompatible with the polynomial exponent semigroup): a **Newton-edge / coprime-
interval obstruction**, the same engine. So the natural tool for JC(2)'s termination is precisely the
coprime-interval return-semigroup machinery, with Fibonacci as the worst case.

## The GMC parallel, honestly bounded

VC(4) is a Laplacian moment nullcone, structurally like GMC's `E=L∘CT` (angular constant-term ⊕ radial
Laplace), so my DvdK-bypass handles the *angular* part. **But** the repo's own analysis (my S205,
CURRENT-FRONTIER) is explicit that this is a *frame, not a route*: GMC(2)⇏JC(2) because the de Bondt doubling
lands VC in dim 4 (rank ≥2), off the rank-1 side of the Frobenius wall that proves GMC(2). And the earlier
"JC(2) and LRC(14) share one n=12 AP-rigidity" is **withdrawn** — S137's "12" is the Fibonacci `144/89`
proxy, not the LRC AP-uniqueness 12. What survives is *methodological*: both localize the obstruction at
worst-approximable (golden) data, both are descent-termination problems, both run on the coprime-interval
engine, and both share the rank-1 seed (THM-1840) and the reify-ladder vertex (transitive ≡ nilpotent ≡ ℓⁿ,
deathstar-S75 — the rational-normal-curve vertex).

## Honest verdict and the positive prior

**JC(2) is not proved here** — no route is complete, and the crux (the dim-4 VC nonvanishing = the descent
termination) is not closed by any current tool, mine or the fleet's. What this session contributes: (1) the
three reduction programs assembled and correctly attributed; (2) the 2D-symmetric case and the Lamé–Fibonacci
worst-case verified; (3) the recognition that JC(2)'s obstruction, in all three of its forms, is one
**descent-termination** problem, for which the **coprime-interval / numerical-semigroup / Frobenius engine**
(S223/S224, and codex's THM-2045 Newton-edge) is the natural tool — with Fibonacci as the worst case. There
is a real **positive prior**: the Keller collision minimum is dimension 3, so **n=2 sits below the collision
threshold** — a structural reason to expect JC(2) *true* and provable by exactly the low-rank/arithmetic
(coprime-interval) means the repo keeps converging on for GMC(2) and LRC(14).

Links: HYP-8905, THM-1435, THM-1345, THM-2045 (codex), THM-1300,
[[one-dimensional-coprime-intervals-complete-the-dvdk-bypass-boxeph-S223]],
[[jacobian-and-lonely-runner-two-nullcones-that-diverge-boxeph-S205]].
