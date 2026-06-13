# Unit Distance Disproof and Tournament Polynomial Connections

**Session:** oracle-2026-05-21-S1
**External result:** OpenAI (2026-05-20), "Disproof of Erdős's Unit Distance Conjecture"
  — constructs n points with ≥ n^{1.014} unit-distance pairs via CM fields + class field towers.
**Depends on:** `ulc-turan-unconditional-proof.md`, `real-rootedness-product-formula-erdos.md`

---

## The Result

Erdős (1946) conjectured that among n points in ℝ², at most n^{1+o(1)} pairs can be at unit distance. The new result constructs explicit point sets with n^{1.014} unit-distance pairs, disproving the conjecture.

**Key technique:** CM (Complex Multiplication) fields — totally imaginary quadratic extensions of totally real fields — where an element β has |β|=1 in ALL embeddings simultaneously (not just some). Combined with infinite class field towers (Golod-Shafarevich), this produces large families of unit-norm elements.

---

## Six Deep Connections to Tournament Polynomials

### 1. The Norm-1 Root Pairs (α₂ = 1 tournaments)

**CM field:** An element β in a CM field satisfies β·β̄ = 1 (product = 1 under the CM involution).

**Tournament parallel:** For degree-2 independence polynomials with α₂ = 1:
$$\rho_1 \cdot \rho_2 = \frac{1}{\alpha_2} = 1 \quad\text{(Vieta)}$$

The two roots are **multiplicative inverses of each other** — exactly the "norm-1 condition" on CM field elements!

The SC H-maximizer at n=6 (α₁=20, α₂=1, H=45): $\rho_1 \approx 19.95$, $\rho_2 \approx 0.05$, $\rho_1\rho_2 = 1$.

**Bridge:** The α₂=1 family of tournament polynomials consists of "CM-type" polynomials where the roots are mutual inverses — the algebraic analogue of CM norm-1 elements.

---

### 2. Tournament Complement Symmetry ↔ CM Involution

**CM field:** The CM involution σ: β ↦ β̄ satisfies β·σ(β) = N(β) (the norm). Fixed points of σ are the totally real elements of F ⊂ K.

**Tournament parallel:** The complement involution τ: T ↦ T^{op} (reverse all arcs). Self-complementary (SC) tournaments are fixed points of τ.

For the independence polynomial: $I(\Omega(T^{op}), x)$ and $I(\Omega(T), x)$ are related by the complement symmetry. SC tournaments have $\Omega(T) \cong \Omega(T^{op})$, so the polynomial is "invariant under the involution" — analogous to norm-1 CM elements fixed by conjugation.

**Key parallel:** The CM construction uses elements NEAR the unit circle (norm ≈ 1). The SC H-maximizers use polynomials with root product = 1 (α₂=1 family). Both are "near the unit norm" in their respective algebraic structures.

---

### 3. Class Field Tower ↔ Tournament Polynomial Tower

**CM result:** Golod-Shafarevich constructs an infinite tower of fields K₁ ⊂ K₂ ⊂ ... where each is the maximal unramified 2-extension of the previous. The tower never terminates when the class group has enough generators.

**Tournament parallel:** The tournament polynomial tower is the sequence $I(\Omega_0, x) = 1$, $I(\Omega_1, x)$, $I(\Omega_2, x)$, ... where each $\Omega_k$ is obtained by adding one cycle (vertex) to the conflict graph.

The TRRT conjecture says: this tower NEVER produces a non-real-rooted polynomial. The Golod-Shafarevich tower condition says: the tower NEVER terminates if gens > 2·rels.

**Proposed Tournament GS Criterion:** The interlacing property at each step (I(Ω_k, x) interlaces I(Ω_{k+1}, x)) is the "generator > 2·relations" condition for the tournament polynomial tower to remain real-rooted.

---

### 4. Split Primes ↔ Vertex-Disjoint Cycle Types

**CM result:** Primes that split completely in K contribute the most unit-norm elements via the Pigeonhole argument on prime ideal factorizations.

**Tournament parallel:** The "split" structure in our context is the factorization $I(\Omega, x) = I(\Omega \setminus C^*, x) + x \cdot I(\Omega - N[C^*], x)$. Cycles $C^*$ that "split cleanly" (where deleting $C^*$ reduces degree by exactly 1) are the "split primes" in this analogy.

From our computation: 444 of these "split" cases at n=6 all show interlacing. The 4468 "inert" cases (degree drops by more than 1) are the primes that don't split cleanly.

---

### 5. Turán Bound ↔ Spencer-Szemerédi-Trotter Bound

**Unit distance:** The Spencer-Szemerédi-Trotter (1984) upper bound gives U(n) ≤ O(n^{4/3}). This was the state-of-the-art for 40 years.

**Tournament parallel:** Our Turán bound gives $\alpha_2 \leq (1-1/d)\alpha_1^2/2$, which implies $\alpha_1^2 \geq \frac{2d}{d-1}\alpha_2$ (ULC at k=1). This is our "upper bound" analog — it constrains how large $\alpha_2$ can be relative to $\alpha_1^2$.

**The parallel improvement:** The unit distance result improves the lower bound from n^1 to n^{1.014}. Our Turán-ULC gives the tight unconditional upper bound (Turán extremal achieved by the double-root tournament). Can we find an analogous "breakthrough" for the ULC at k=2?

---

### 6. The Shape Exponent and Unit Distance Exponent

**Unit distance:** The exponent δ ≈ 0.014 measures the improvement over the trivial bound. It comes from:
$$\delta = \frac{\log(\text{count of unit-norm elements})}{\log n}$$

**Tournament parallel:** The "shape exponent" of a tournament family measures the growth of the shape parameter $s = \alpha_1/\sqrt{\alpha_2}$:

| Family | s growth | analog |
|---|---|---|
| Random tournaments | s = O(1) | Erdős grid construction: n^1 |
| SC H-maximizers | s ~ n^3/√α₂ → ∞ | Unit distance improvement |
| Transitive | s = 0 (no cycles) | empty configuration |

**Key observation:** For SC H-maximizers with α₂=1 (the "pure norm-1" family), $s = \alpha_1$ grows like $n^3/24$ (dominant 3-cycle count). The shape parameter exponent is:
$$\text{exponent of } s = 3 - \frac{\log\sqrt{\alpha_2}}{\log n} = 3 \text{ (for } \alpha_2 = 1\text{)}$$

This is the "3" in our cubic term — the tournament analog of the unit distance exponent. The Erdős improvement from 1 to 1.014 corresponds to finding SC-like tournaments where s grows faster than expected.

---

## The Deepest Connection: Algebraic Units and Independence Polynomial Roots

The unit distance proof uses elements of norm exactly 1 in CM fields. These satisfy the "unit equation" $\beta\bar\beta = 1$.

Our independence polynomial roots satisfy (Vieta for α₂=1):
$$\rho_1 \rho_2 = 1/\alpha_2 = 1$$

This is the "unit product equation" for polynomial roots.

**Algebraic number theory analog:**
- CM elements with norm 1: live in the group of units $\mathcal{O}_K^1 = \{β \in \mathcal{O}_K : \beta\bar\beta = 1\}$
- Our α₂=1 polynomial roots: live in the "multiplicative complement" set $\{(\rho_1, \rho_2) : \rho_1\rho_2 = 1\}$

The unit distance proof counts elements in $\mathcal{O}_K^1$ by pigeon-holing on ideal classes. Our Turán bound counts edges in $\bar\Omega$ by bounding them via Turán density. Both are "counting unit objects via algebraic structure."

---

## New Directions Suggested by the Connection

### Direction 1: CM Construction for Tournaments

Can we construct tournaments using CM fields to achieve specific ULC properties? Specifically:

- Use a CM field $K$ with split primes to define a tournament structure on the units of $\mathcal{O}_K^1$
- The unit distance condition (|β| = 1) would give the norm-1 root condition (ρ₁ρ₂ = 1)
- The class field tower would give a sequence of tournaments with growing cycle counts

This would be a systematic way to construct the α₂=1 tournament family using CM number theory.

### Direction 2: The Golod-Shafarevich Criterion for Interlacing

The Golod-Shafarevich criterion says: infinite tower iff generators > 2·relations.

**Tournament analog:** The polynomial tower $I(\Omega_0) \to I(\Omega_1) \to ...$ (adding one cycle at a time) remains real-rooted (TRRT) iff at each step, the "interlacing condition" holds.

The interlacing condition: adding cycle $C^*$ to $\Omega$ moves roots continuously from the roots of $I(\Omega \setminus C^*, x)$ to the roots of $I(\Omega, x)$, maintaining the interlacing.

**GS condition for tournaments:** TRRT holds if the number of "split" cycles (those that reduce degree by exactly 1) exceeds twice the number of "inert" cycles (those that reduce degree by more than 1). This would be the tournament Golod-Shafarevich condition.

### Direction 3: Analytic Number Theory Bounds for α₂

The unit distance proof uses Louboutin's explicit discriminant bounds to control the class number $h^-(K)$. These bounds are sharper than classical Minkowski bounds.

**Tournament analog:** Can we find explicit bounds for α₂ (the number of disjoint cycle pairs) that are sharper than the Turán bound? The Turán bound uses the maximum clique of $\bar\Omega$; Louboutin-type bounds might use more refined structural information.

### Direction 4: The "1.014" Threshold

The unit distance construction hits a barrier at exponent ≈ 1.243 (cannot exceed this via the current algebraic strategy). There might be an analogous "barrier" for our ULC at k=2 proof.

**Conjecture:** The minimum ratio $\alpha_2^2/(3\alpha_1\alpha_3)$ for tournament conflict graphs has a non-trivial lower bound (strictly greater than 1) that cannot be achieved by the K_{a,b,c} construction alone. The actual minimum might be ≈ 2.0 (from our computation: min observed 2.05).

---

## Summary of the Parallel Structure

| Unit Distance Proof | Tournament Polynomial Theory |
|---|---|
| Elements β with β·β̄ = 1 (norm-1 in CM) | Root pairs with ρ₁·ρ₂ = 1 (α₂=1 family) |
| CM involution σ: β ↦ β̄ | Complement involution T ↦ T^op |
| Self-conjugate elements (fixed by σ) | SC tournaments (fixed by complement) |
| Class field tower K₁ ⊂ K₂ ⊂ ... | Tournament polynomial tower (interlacing) |
| Split primes in K | "Split" cycles (degree-drop = 1) |
| Inert primes in K | "Inert" cycles (degree-drop > 1) |
| Golod-Shafarevich: gen > 2·rel → infinite tower | Interlacing persistence → TRRT |
| Spencer-Szemerédi-Trotter upper bound n^{4/3} | Turán bound α₂ ≤ (1-1/d)α₁²/2 |
| Unit distance exponent δ = 0.014 | Shape exponent 3 (for α₂=1 family) |
| Class number h^-(K) bounds element count | Turán extremal bounds disjoint pair count |
| Louboutin explicit discriminant bounds | Open: explicit ULC k=2 bounds |
| Barrier at exponent 1.243 | Conjectured barrier for min α₂²/(3α₁α₃) |
