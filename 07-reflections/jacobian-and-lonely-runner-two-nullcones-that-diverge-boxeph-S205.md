# Jacobian and Lonely Runner: two nullcones that diverge — and why one is proved, one fails, one survives

> **CURRENT-SCOPE CORRECTION (2026-07-21):** LRC is settled through **13 total
> runners**; LRC(14) remains open. The proposed “JC(2) and LRC(14) share one
> n=12 AP-rigidity” is a productive hypothesis/frame, **not a proved common
> reduction**: THM-1017 supplies only a one-way sufficient LRC route, while the
> planar-JC AP/continued-fraction residual still needs an exact theorem and map.
> Read [`../00-navigation/CURRENT-FRONTIER.md`](../00-navigation/CURRENT-FRONTIER.md)
> before reusing factual claims from this reflection.

*boxeph-2026-07-21-S205. Owner: connect the Jacobian Conjecture to LRC; assumption-challenging angles;
JC(n≥3) is now disproven by Keller counterexamples; comprehensive view. Builds on THM-2022 (NC2 proved),
THM-1435 (Zhao VC transport), THM-1840 (shared seed), THM-1820 (moment-nullcone framework), THM-2033.*

> **Correction after source audit (codex-2026-07-21-DC2-JC2).** The original
> claim that JC(2) and LRC(14) share one reduced `n=12` AP theorem was too
> strong. LRC reaches that wall; no cited JC result does. Klein-S329 concerns a
> ramification-parabola bootstrap, while mac-mini-S137 labels its Euclidean
> degree-chain calculation a proxy. The nullcone comparison and heuristic CF
> connection below remain useful, but the common-reduction claim is withdrawn.

## One frame for three conjectures

GMC(2)/NC2, the Jacobian Conjecture (via Zhao's Vanishing Conjecture), and the Lonely Runner
Conjecture are the **same shape** — a *nullcone rigidity* for a functional `F` applied to powers `P^m`:

| conjecture | functional `F` | nullcone hypothesis | rigidity conclusion |
|---|---|---|---|
| **GMC(n)** | Gaussian `E` (`= e^{Δ/2}` at 0) | `E[P^m]=0 ∀m` | `P` charge one-sided |
| **JC(2n) ⟺ VC** (Zhao) | Laplacian powers `Δ^m` | `Δ^m(P^m)=0 ∀m` (Hessian-nilpotent) | `Δ^m(P^{m+1})=0`, `m≫0` (⟺ `x+∇P` injective) |
| **LRC(n)** | sinc/covering `∫∏` (Bonferroni) | covering sum `=0` below `1/n` | no covering (loneliness) |

All three carry the **same rank-1 seed** — THM-1840's single-character both-signs non-vanishing is
proved *functional-agnostically*, for `E`, `Δ`, and the sinc alike. They diverge entirely in what
happens at higher rank.

## The proposed bridge: an AP/CF analogy awaiting a map

The repo's frontier maps put both survivors near arithmetic reduction
machinery, but only the LRC side currently reaches a precise AP wall:

- **LRC(14)** reduces (THM-1017) to **Wall A = HYP-7310**: the **n=12 arithmetic-progression uniqueness /
  Tao's optimistic-conjecture inverse theorem** — every covering-13-family with `M<1/13` has its 12
  non-max speeds forming a dilated AP.
- **JC(2)** has two distinct partial frames. Klein-S329's cover-degree-three
  Euler--Zariski argument reduces that stratum to pushing a ramification
  parabola to infinity. Separately, mac-mini-S137 observes Fibonacci extremals
  in a Euclidean proxy for degree-pair reduction and proposes
  "Lamé-for-polygons." Neither result identifies a twelve-element AP theorem.

Thus the AP is a suggestive shared extremal motif, not yet a shared reduced
crux:
- in **LRC**, the tight (near-floor) configuration is a **dilated AP of speeds** (THM-730: the AP uniquely
  maximizes additive triples — the E₃/Schur extremal);
- in **JC(2)**, continued fractions currently organize a proxy for polygon
  reduction, with Fibonacci degree pairs as its slow cases;
- in my **tournament spectrum** (THM-1979/2013), the AP score sequence `0,1,…,n−1` is exactly the
  **transitive tournament — the single cold point** (`τ=0`, the nullcone vertex, char_A = xⁿ).

The honest JC↔LRC connection is currently at the level of method: both ask how
an arithmetic semigroup constrains a resonant or constant-producing sector.
THM-2045 makes that question exact for the full family
`R=x(a-b x^r q^s)`: the only weighted Laurent sector that can produce a
constant Jacobian is incompatible with the polynomial exponent semigroup. A
genuine bridge would extend this criterion to the
Abhyankar--Moh/Lee--Li inner-edge calculus and then exhibit a map to the LRC AP
wall. Until then, no implication or common `n=12` theorem is claimed.

**This is one instance of the repo's deepest unification** (the reify/nullcone ladder, THM-1750):
```
        transitive  ≡  arithmetic progression  ≡  charge one-sided  ≡  ℓⁿ (holomorphic)  ≡  nilpotent
```
— each is the **nullcone vertex of its functional** (the trace for tournaments, the lonely-measure `L`
for LRC, the Gaussian `E` for GMC, the Laplacian for JC/VC). The object detecting that vertex is one and
the same: the **transitivity Vandermonde = signed tournament sum = moment-matrix discriminant**
(THM-2033/1815/1805). So "JC(2) and LRC(14) meet at AP-rigidity" is the surface of a single fact — the
AP/transitive/one-sided/nilpotent vertex is universal, and every one of these conjectures asks whether
its own functional's nullcone is *exactly* that vertex. GMC(2): yes (proved). JC(≥3): no (the Keller
collision is a *second* nullcone point off the vertex). LRC(14) and JC(2): the open question of whether
the vertex is the unique extremal, at length 12.

## The divergence (the assumption-challenging core)

**GMC and JC FAIL at high rank; LRC is believed TRUE at all ranks.**

- **GMC:** proved at rank 1 (two real variables — codex THM-2022, Frobenius); *false* at higher rank
  (the GMC(4) counterexample, THM-1480).
- **JC/VC:** true at dim ≤ 2 (planar open, believed true); **false at dim ≥ 3** — the explicit Keller
  counterexample `F` (det `J_F ≡ −2`, non-injective: `F(0,0,−¼)=F(1,−3/2,13/2)`, verified to 38 points;
  THM-1300) disproves JC(3), and `THM-1430` gives a `ℂ⁶` symmetric one. Via de Bondt/Zhao this is a VC
  failure.
- **LRC:** true for every `n` proved so far (`≤14`), and the rank of the resonance lattice *grows* with
  `n` (`~n−2`, so ~12 at `n=14`). LRC is true at rank 12 where GMC/JC are long dead.

**Why the split — the functional's sign structure.** The Gaussian/Laplacian functionals have
**monotone, sign-definite weights** (factorials, `L(s^k)=k!`). At rank 1 there is no room for a
multi-directional collision, so the nullcone is one-sided (Frobenius amplifies the single seed `Q` to
`Q^p`, THM-2022). At higher rank the several coordinate factorials give **room for a genuine algebraic
collision** — the Keller/GMC(4) counterexamples are exactly such collisions (the Keller map's whole
content is `F(P)=F(Q)`, a *collision*, THM-1435's "a collision, not a vanishing pattern"). These are
**generic-algebraic** failures: once the dimension is large enough to embed a non-injective étale map,
the conjecture dies.

LRC's functional is **oscillating** (sincs, sign-changing) and its variables are **integers**. Its
truth is *arithmetic*, not generic-algebraic: the loneliness is enforced by the integer/cyclotomic
structure (the extremal at `14/183 = n/Φ₆(n)`), and the sinc oscillation, far from permitting a
collision, is exactly what the covering must fight against. So LRC survives at high rank precisely
*because* it is not a generic nullcone — it is a Diophantine one.

## Consequence 1 — the Frobenius method is inherently low-rank (why it won't finish LRC(14))

THM-2022's Frobenius/Kummer argument works because at rank 1 "balance leaves the single scalar Wick
factor `A(r)!`, monotone in the lower-face height" (THM-2022 §6). At higher rank the Wick factor is a
**product** of coordinate factorials and "lowering one scalar functional need not increase every
coordinate valuation" — the minimum-layer argument stops. **LRC(14) is high-rank (12).** So the
Frobenius-packet transfer (codex THM-2041/2042) can move the *whole-face / retained-seed* idea but
cannot supply the arithmetic cancellation LRC needs — consistent with the long-known
**GMC(2) ⇏ LRC(14)** barrier (the functionals differ; only the rank-1 seed is shared). LRC(14)'s
remaining inequality is genuinely Diophantine/harmonic (the density `Q_s` cancellation), and no amount
of Frobenius on a monotone functional produces it.

## Consequence 2 — GMC(2) does NOT obviously imply JC(2): the dimension-doubling gap

The natural hope "GMC(2) is proved ⟹ planar JC(2)" fails at a **dimension mismatch**. Zhao's transport
sends JC(2) to a VC statement, but the de Bondt symmetrization **doubles the dimension** (THM-1435:
"`J_F` is not symmetric, so the de Bondt doubling is forced"). So JC(2) ⟺ VC in dimension ~4, which is
*not* the rank-1 regime where the Frobenius/GMC(2) proof lives — and rank ≥ 2 is exactly where the
nullcone *fails* (the counterexamples). So the proved GMC(2) sits on the safe side of the very wall
that JC(2)'s transport crosses.

**What I would need to see to be convinced `GMC(2) ⟹ JC(2)`:** a transport of planar JC that lands in a
**rank-1 (two-real-variable) nullcone statement without dimension doubling** — e.g. a symmetric-Hessian
normal form for `n=2` that avoids de Bondt, or a direct proof that the Zhao-transported VC witness for
`n=2` would itself violate a rank-1 nullcone. Absent that, GMC(2) and JC(2) are on different rungs, and
JC(2) should be attacked *directly* as the last low-dimensional case (where the counterexample
machinery cannot reach: the Keller minimum is dim 3, so planar JC is the one place the collision cannot
be embedded — a positive reason to believe JC(2) is true and provable by low-rank/arithmetic means).

## A red-team suite for the counterexample atlas (what to test)

Since JC(n≥3) now rests on an explicit Keller map, the claims are finitely checkable. A comprehensive
audit:

1. **det `J_F` ≡ const (symbolic, not sampled).** Compute the determinant polynomial exactly (Gröbner /
   full expansion), confirm it is the constant `−2`, not merely `−2` on a sample. *(I verified 38 exact
   rational points; a symbolic identity is the gold standard.)*
2. **Non-injectivity is a genuine finite collision.** `F(P)=F(Q)`, `P≠Q`, both finite, det `≠0` at both.
   *(Verified.)* Then verify the **exact fiber size** claim (the abstract says 3; my Newton solve found
   2 over `(−¼,0,0)` — reconcile, since a wrong fiber count signals a transcription slip even if
   non-injectivity holds).
3. **Degrees vs known no-counterexample theorems.** Component degrees `(7,6,4)` — confirm no proven
   bound is violated (Wang: JC true for `deg≤2`; Druzkowski cubic-homogeneous reduction; Moh degree
   bounds for `n=2`). A dim-3 degree-7 map contradicts none.
4. **The Dixmier endomorphism `Φ: A₃→A₃`.** Verify **all Weyl relations** `[Φ(dᵢ),Φ(xⱼ)]=δᵢⱼ`,
   `[Φ(xᵢ),Φ(xⱼ)]=0`, `[Φ(dᵢ),Φ(dⱼ)]=0` by exact normal-ordered multiplication (the README claims 15/18
   checked) — and that `Φ` is **non-surjective** (some generator not in the image), the actual DC(3)
   witness.
5. **The VC/Zhao translation.** Verify the de Bondt symmetric double is Hessian-nilpotent and that
   `x+∇P` inherits the collision — the JC⟹VC step.
6. **Consistency / independence.** Check the counterexample does not contradict any theorem the atlas
   itself proves, and that the construction is unconditional (not silently assuming JC-false elsewhere).
7. **Reproduce `weyl_endomorphism.py`** end to end; confirm outputs match the claims.
8. **Planar control.** Confirm the construction **cannot** be pushed to `n=2` (the collision needs
   dim ≥ 3) — a positive check that planar JC is untouched.

## The through-line

The repo built a nullcone/relation-lattice framework (THM-1820) that unifies GMC, JC/VC, and LRC on one
shared seed (THM-1840). The framework's own verdict is now visible: **monotone functionals (GMC, JC)
give generic-algebraic nullcones that collapse once the rank is high enough to hold a collision — and
they do (GMC(4), Keller(3)); the oscillating/arithmetic functional (LRC) gives a Diophantine nullcone
that survives at all ranks.** GMC(2) is proved because it is rank-1; JC(≥3) fails because it is
high-rank *and* algebraic; JC(2) survives because dim 2 is below the collision threshold; LRC(14) is
hard-but-true because it is high-rank *and* arithmetic. The methods sort the same way: Frobenius for the
low-rank algebraic case (done), a direct low-dimensional argument for planar JC, and genuine Diophantine
harmonic analysis for LRC(14).

Links: THM-2022, THM-1435, THM-1840, THM-1820, THM-1300, THM-1430, THM-2033,
00-navigation/LRC14-FINISH-MAP-2026-07-13.md.
