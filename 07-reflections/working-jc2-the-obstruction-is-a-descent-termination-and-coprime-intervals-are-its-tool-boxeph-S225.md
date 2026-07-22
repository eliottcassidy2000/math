# Working JC(2): the obstruction is a descent termination — and coprime intervals are its tool

*boxeph-2026-07-21-S225. Owner: spend a long session working to prove the planar Jacobian Conjecture; pull
in past threads creatively. Builds on (corrected via mining): de Bondt+Zhao VC (THM-1435), the Keller
counterexamples (THM-1300/1430), klein-S329 (Euler–Zariski bootstrap), mac-mini-S137 (golden-corner /
Lamé), THM-1345/1365/1370 (restricted JC(2) proved), codex THM-2044/2045/2049 (DC2/Newton-edge), my
S222/S223 (DvdK bypass / coprime intervals), S224 (Wall A), S206 (Fibonacci foil). Verified in
`04-computation/jc2_via_vanishing_conjecture_and_the_cf_termination_boxeph_S225.py`.*

> **Scope correction (MISTAKE-237).** The three programs below are not known
> to be equivalent. VC(4) is a rigorous symmetric-reduction target;
> leading-form/Jelonek descent is a planar route; Lame-for-polygons is a
> proposed effective termination principle. Their common descent vocabulary
> is a research heuristic until explicit implications are proved.

## Honest starting point

JC(2) is the sole open Jacobian case (JC(n≥3) is *false* — the explicit Keller map, `det≡−2`, triple
collision, THM-1300) and, in the repo's words, "the graveyard of false proofs." I am *working* it, not
claiming it. The value of a session like this is to assemble the routes and
verify individual pieces. The routes below do not yet have a single proved
shape or a transfer theorem between them.

## Three reduction programs (all in the repo, none complete)

- **(A) de Bondt + Zhao → a four-variable symmetric/VC target.** de Bondt reduces JC(n) to the *symmetric/gradient* case (`H=∇P`), which
  by Zhao is the **Vanishing Conjecture**: `Δᵐ(Pᵐ)=0 ∀m ⟹ Δᵐ(P^{m+1})=0` for `m≫0`. The reduction
  **doubles the dimension** (`n→2n`), so **JC(2) ⟺ VC(4)** — a *Laplacian moment nullcone* in 4 variables.
- **(B) leading-form / Jelonek geometry.** Degree two is excluded and degree
  three is the first live geometric stratum. It is not the whole residual:
  higher generic degrees remain live, and the proposed cusp/ramification
  bootstrap is not a classification.
- **(C) mac-mini-S137 golden-corner / Lamé.** JC(2)'s elementary reductions act like **subtractive Euclid**
  on Newton-polygon slopes; the target "Lamé-for-polygons" is a `≤ C·log(min degree)` termination bound.

## What is solved (verified)

- **One binary homogeneous symmetric family is exact.** For
  `P=A(x+iy)^d+B(x-iy)^d`, `d>=2`, the determinant calculation forces `AB=0`;
  THM-2063 then makes the gradient map tame. The stored script checks the
  one-sided monomials, not the full symmetric target.
- **Restricted JC(2) proved** (from the mining): equivariant (THM-1345), elliptic in all dimensions
  (THM-1370), polynomial-Galois (THM-1365), geometric degree ≤2 (fails first only at degree 3),
  one-fiber-linear planar Keller (THM-2063). These are necessary exclusions,
  not a claim that only cover-degree three survives.

## The Lamé/Fibonacci crux (verified) — and its meaning

The stored computation verifies only that the ordinary Euclidean algorithm's
longest chain among pairs `<200` occurs at `(144,89)`, as Lamé's theorem
predicts. No decreasing invariant has been proved for the proposed JC(2)
Newton reduction, so Fibonacci is a test family, not a verified Keller
worst-case or a transfer from LRC/GMC.

## Three descent-shaped programs, without an equivalence

Pulling the three programs together suggests three potentially interacting
forms of the JC(2) obstruction. They are all the **same kind of problem** at
the level of search design, but no equivalence is claimed:

| form | statement | it is a … |
|---|---|---|
| (A) VC(4) | the dim-4 both-signs radial (Laplacian) nonvanishing | moment-return nonvanishing |
| (B) leading-form descent | `{P_A,Q_B}=0` must propagate to a global coordinate (THM-1345 §5) | a filtration descent that must terminate |
| (C) Lamé-for-polygons | the Euclid/Newton-slope reduction length is `O(log)` | a continued-fraction descent bound |

The common vocabulary suggests experiments, but MISTAKE-234 retracts the
mixed-sign semigroup noncancellation claim and MISTAKE-237 retracts the
equivalence. THM-2045 is a genuine exponent-sector obstruction in one
factorized family; it does not install a coprime-interval termination engine
for arbitrary Keller pairs. Any future transfer must exhibit the decreasing
quantity and prove that coordinate changes cannot reset it.

## The GMC parallel, honestly bounded

VC(4) is a Laplacian moment nullcone and thus has moment-like syntax, but no
DvdK bypass has been mapped to an angular part of the VC target. The repo's own analysis (my S205,
CURRENT-FRONTIER) is explicit that this is a *frame, not a route*: GMC(2)⇏JC(2) because the de Bondt doubling
lands VC in dim 4 (rank ≥2), off the rank-1 side of the Frobenius wall that proves GMC(2). And the earlier
"JC(2) and LRC(14) share one n=12 AP-rigidity" is **withdrawn** — S137's "12" is the Fibonacci `144/89`
proxy, not the LRC AP-uniqueness 12. What survives is a research prompt: test
whether a lawful decreasing quantity sees badly approximable slopes. No
common coprime-interval engine or rank-one transfer is proved.

## Honest verdict and the positive prior

**JC(2) is not proved here.** What survives is (1) a list of three separate
programs; (2) the exact binary homogeneous subcase; and (3) the independent
Euclidean/Fibonacci toy calculation. The explicit dimension-three
counterexample explains why JC(2) is now the sole live dimension, but it is not
evidence that any one of these methods must close it.

Links: HYP-8905, THM-1435, THM-1345, THM-2045 (codex), THM-1300,
[[one-dimensional-coprime-intervals-complete-the-dvdk-bypass-boxeph-S223]],
[[jacobian-and-lonely-runner-two-nullcones-that-diverge-boxeph-S205]].
