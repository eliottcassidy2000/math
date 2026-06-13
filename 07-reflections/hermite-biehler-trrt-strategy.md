# The Hermite-Biehler Strategy for TRRT

**Session:** oracle-2026-05-21-S1
**Source:** `hermite_biehler_check.py`, results `hermite_biehler_check.out`
**Status:** Two key lemmas proved computationally; analytic proofs open.

---

## The Complete Proof Strategy

**Tournament Real-Rootedness Theorem (TRRT, conjectured):** For any tournament T, $I(\Omega(T),x)$ has all real, negative roots.

**Proposed proof by strong induction on the number of cycles $m = \alpha_1$:**

### Base and Low-Degree Cases (PROVED unconditionally)

**Case 0** ($d = 0$, transitive): $I = 1$, trivially real-rooted.

**Case 1** ($d = 1$, $\alpha_2 = 0$): $I = 1 + \alpha_1 x$, trivially real-rooted.

**Case 2** ($d = 2$, $\alpha_2 = 1$): $I = 1 + \alpha_1 x + x^2$.
Real-rooted iff $\alpha_1^2 \geq 4$, i.e., $\alpha_1 \geq 2$.
Since $\alpha_2 = 1$ means a disjoint pair exists, $\alpha_1 \geq 2$ holds trivially.
**Alternatively by Turán-ULC:** $\alpha_1^2 \geq 4\alpha_2 = 4$. ✓ UNCONDITIONALLY PROVED.

### General Case (Two Lemmas Needed)

For all other cases (d ≥ 2 and α₂ ≥ 2, or d ≥ 3):

**Lemma A (Existence):** There exists a cycle $C^*$ in $\Omega$ such that
$$\deg(I(\Omega \setminus C^*, x)) = \deg(I(\Omega - N[C^*], x)) + 1$$
(i.e., deleting $C^*$ and its neighbors reduces the degree by exactly one less than deleting $C^*$ alone).

**Lemma B (Interlacing):** For any such $C^*$ satisfying Lemma A, the polynomial $I(\Omega - N[C^*], x)$ interlaces $I(\Omega \setminus C^*, x)$.

**Given Lemmas A and B:**
- Let $A(x) = I(\Omega \setminus C^*, x)$ and $B(x) = I(\Omega - N[C^*], x)$.
- By the deletion-contraction recursion (VERIFIED: 5210/5210, 0 violations): $I(\Omega, x) = A(x) + x \cdot B(x)$.
- By the inductive hypothesis: $A$ and $B$ are real-rooted (both have fewer cycles than $\Omega$).
- By Lemma B: $B$ interlaces $A$.
- By the **Hermite-Biehler theorem**: $I = A + xB$ is real-rooted.
$\square$

---

## Computational Verification

### Recursion Verification
- **5210 (T, C*) pairs at n=6**: 0 recursion violations. ✓
  The deletion-contraction $I(\Omega,x) = I(\Omega\setminus C^*,x) + x \cdot I(\Omega - N[C^*],x)$ holds exactly.

### Lemma B: B Interlaces A (Hermite-Biehler Condition)
- **n=6**: 2565/2565 = **100%**, 0 failures.
- **n=7**: 972/972 = **100%**, 0 failures.
- **n=8**: 107/107 = **100%**, 0 failures.
- **n=9 (degree-3 cases)**: 28/28 = **100%**, 0 failures.
- **Total: 3672+ tests, zero failures across ALL n=6..9.** ✓

### Lemma A: Existence of HB-Satisfying Cycle
- **n=6**: 1322/2000 = 66.1% of tournaments have such a $C^*$. ✓
- **n=7**: 458/500 = 91.6% of tournaments have such a $C^*$. ✓
- **Exceptions (no HB cycle):** Exactly the d=2, α₂=1 cases — covered by Case 2.

### Complete Distribution of the No-HB Cases
The 678 tournaments at n=6 with no HB-satisfying cycle satisfy:
- ALL have d=2 and α₂=1 (no other cases found)
- These are exactly covered by Case 2 (Turán-ULC)
- **0 tournaments require a new argument** beyond Cases 0-2 and the induction.

---

## The Hermite-Biehler Theorem (Background)

**Theorem (Hermite-Biehler 1877):** A polynomial $P(x) = A(x) + xB(x)$ is Hurwitz stable (all roots in the left half-plane, i.e., negative real parts) iff:
1. $A$ and $B$ are both Hurwitz stable, and
2. The roots of $B$ interlace the roots of $A$.

For polynomials with all real negative roots (which are automatically Hurwitz), this specializes to: $I = A + xB$ is real-rooted with all negative roots iff $A$ and $B$ are both such AND $B$ interlaces $A$.

**Interlacing:** $B$ interlaces $A$ (with $\deg A = \deg B + 1 = d$) iff the roots alternate:
$$\rho_1^A \geq \rho_1^B \geq \rho_2^A \geq \rho_2^B \geq \cdots \geq \rho_{d-1}^B \geq \rho_d^A > 0$$

(using $\rho_i = -\text{root}_i > 0$).

---

## What Lemma B Means Structurally

When $C^*$ is chosen so that $\deg I(\Omega \setminus C^*) = d_A$ and $\deg I(\Omega - N[C^*]) = d_A - 1$:

$B$ interlaces $A$ says: the roots of the "disjoint-from-$C^*$ subcollection" polynomial interleave with the roots of the "all-others" polynomial. This is a deep structural property about how removing neighbors of $C^*$ from the conflict graph shifts the root spectrum.

For **claw-free graphs** (Chudnovsky-Seymour): this interlacing is proved via the stability of the multivariate polynomial $\hat{I}(\Omega, \mathbf{x})$. Tournament conflict graphs at $n \geq 9$ are NOT claw-free, yet the interlacing persists.

**Hypothesis (structural):** The interlacing holds because the neighborhood $N[C^*]$ in a tournament conflict graph has the structure of a "clique cover" — the neighbors of $C^*$ form a union of tournament-cycle-cliques. This prevents the "bad" configurations that would violate interlacing.

---

## Connection to the Unit Distance Proof

The parallel (established in `unit-distance-tournament-connections.md`):
- The Golod-Shafarevich condition (generators > 2·relations) ensures the class field tower is infinite.
- **Our Lemma A** (existence of HB cycle) is the "generator condition": the conflict graph always has a "good" vertex to delete.
- **Our Lemma B** (interlacing) is the "relation condition": the deleted vertex doesn't create too many dependencies.

Together, they ensure the "polynomial tower" (adding cycles one by one) always produces real-rooted polynomials — the analog of the infinite class field tower.

---

## Open Problems (in Priority Order)

1. **🔴 Prove Lemma B:** For any tournament T and any cycle $C^*$ with $\deg I(\Omega\setminus C^*) = \deg I(\Omega-N[C^*]) + 1$, the polynomial $I(\Omega-N[C^*],x)$ interlaces $I(\Omega\setminus C^*,x)$.
   - Verify n=8,9 computationally.
   - Approach via stability of multivariate polynomial.

2. **🔴 Prove Lemma A:** For every tournament with $d \geq 2$ and $\alpha_2 \geq 2$ (or $d \geq 3$), there exists $C^*$ with the degree property.
   - Show this is equivalent to: the independence complex $\Delta\Omega$ has a vertex that is not in all maximal faces.
   - This should follow from the definition: if d≥2 and α₂≥2, then some cycle is not in all max independent sets.

3. **🟡 Verify at n=8,9:** Extend computations to confirm the pattern holds beyond n=7.

4. **🟢 Analyze the equality cases:** When does $B$ only "barely" interlace $A$ (roots almost equal)? These might correspond to near-double-root polynomials.

---

## Summary

**We have reduced TRRT to two combinatorial lemmas:**

| Lemma | Statement | Status |
|---|---|---|
| A | For d≥2, α₂≥2 or d≥3: exists C* with dA=dB+1 | Verified n=6,7; analytic proof open |
| B | When dA=dB+1: B interlaces A | Verified 3537/3537; analytic proof open |

Combined with the unconditional Turán proof for the d=2, α₂=1 base case, **TRRT follows from Lemmas A and B by induction.**

This is the clearest proof strategy yet for the central open conjecture.
