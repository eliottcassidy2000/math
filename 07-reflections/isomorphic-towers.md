# Isomorphic Towers: What Shares the Same Shape

*opus-2026-03-23-S266*

## The Three Towers We Know

### Tower 1: Cayley-Dickson (algebraic doubling)
```
R(1) → C(2) → H(4) → O(8) → S(16) → ...
dim:  1    2     4     8      16
lost: —  commut  assoc  alter  power-assoc
n:    2    3     5     9      17
```
Each step doubles dimension and loses one algebraic property.

### Tower 2: Burnside (symmetry quotients)
```
2^m (labeled) → V_n (iso classes) → V_merged → V_even → V_score
exponent:  m    arc_orbits      arc_orbits-1   cycle_null   k-1
quotient: —        S_n              Z_2         cycle proj   cut proj
```
Each step quotients by a symmetry, reducing the Burnside exponent.

### Tower 3: Cartan (Lie algebra decomposition)
```
gl(n) = so(n) ⊕ p ⊕ R·I
dim:     n²    C(n,2)  C(n+1,2)-1   1
type: antisymm  symm    scalar
tour: tournament  cooperation  identity
```
The antisymmetric (directed) + symmetric (undirected) + scalar.

## The Common Pattern

All three towers share **five structural features**:

1. **Progressive decomposition**: A complex object splits into layers
2. **Property loss at each level**: Each descent loses information
3. **Factored degrees of freedom**: The exponent/dimension splits additively
4. **Asymptotic independence**: Layers decouple at large n
5. **Controlled by 2 and 3**: The generators have orders 2 and 3

## Seven More Towers with the Same Shape

### Tower 4: Postnikov (homotopy decomposition)
```
X → ... → X_{n} → X_{n-1} → ... → X_0 → *
```
A topological space decomposes into its homotopy groups, one level at a time. Each Postnikov section X_n captures π_k for k ≤ n and kills all higher homotopy. The k-invariants (obstructions to extending sections) are the "residuals" at each level.

**Isomorphism**: The Burnside tower quotients tournament space level by level. The Postnikov tower quotients homotopy type level by level. Both factor a complex object into a sequence of simpler pieces, with the obstruction to reconstruction living in a specific cohomology group.

**Tournament analog**: V_n is the "π_0" (connected components = iso classes). The clique complex Betti numbers β_1, β_2 are the "higher homotopy groups" of tournament space.

### Tower 5: Chromatic (stable homotopy)
```
S^0 → ... → L_n S^0 → L_{n-1} S^0 → ... → L_0 S^0 = S^0_Q
```
The sphere spectrum localizes at successive chromatic heights. Each height n captures phenomena controlled by the nth Morava K-theory K(n). The chromatic convergence theorem says the tower is exhaustive.

**Isomorphism**: The Burnside exponent factorization arc_orbits = cycle_null + (k-1) is a "chromatic decomposition" — the cycle part (height 1) and the cut part (height 0) are the two chromatic levels of tournament theory.

**Why only 2 levels**: Tournaments live in the "p=2 chromatic world." The cycle space is GF(2)-linear. The cut space is GF(2)^{n-1}. Higher chromatic heights would require p=3 structure, which appears only in oriented graphs (3 states per cell).

### Tower 6: Central Series (group theory)
```
G = G_0 ⊃ G_1 ⊃ G_2 ⊃ ... ⊃ G_c = {e}
```
The lower central series of a group, where G_{i+1} = [G, G_i]. Each quotient G_i/G_{i+1} is abelian. The nilpotency class c measures "how far from abelian" the group is.

**Isomorphism**: The tournament's "nilpotency" is measured by the residual sequence. residual(n)/T → 0 means the cut-cycle interaction becomes nilpotent. At n=3: residual=0 (abelian = cut and cycle commute). At n≥4: residual>0 (non-abelian = cut and cycle don't commute).

### Tower 7: Spectral Sequence (homological algebra)
```
E_2^{p,q} ⟹ H^{p+q}(X)
```
A spectral sequence converges to a target through successive pages E_r. Each page kills some differentials. The E_∞ page gives the associated graded of a filtration on the target.

**Isomorphism**: The Burnside sum is a "spectral sequence" converging to V_n. The "E_2 page" is the identity term (2^m/n!). The "differentials" are the non-identity corrections. The spectral sequence degenerates rapidly (identity dominates at 99.9% by n=11).

### Tower 8: Renormalization Group (physics)
```
Λ → Λ/b → Λ/b² → ... → 0 (IR fixed point)
```
In quantum field theory, the renormalization group flows from high energy (UV) to low energy (IR), integrating out degrees of freedom at each scale. Relevant operators survive; irrelevant operators vanish.

**Isomorphism**: The Burnside tower IS a renormalization flow. "High energy" = labeled tournaments (2^m degrees of freedom). "Low energy" = iso classes (V_n). The S_n action "integrates out" the labeling degrees of freedom. The twin_SL is the "relevant operator" (survives at all scales). The residual is the "irrelevant operator" (vanishes at large n).

The master formula E ≈ T × (2^{n-1}-2)/2^n is the "beta function" — it describes how E scales with n. The fixed point is E/T → 1/2 (every arc reversal changes class).

### Tower 9: Adams Resolution (homological)
```
X → K_0 → K_1 → K_2 → ...
```
An Adams resolution replaces a spectrum X by a sequence of Eilenberg-MacLane spectra, with each K_i capturing one layer of information.

**Isomorphism**: The tournament decomposes as T = score (cut) + cycles (even graph) + interaction (residual). Each layer is "resolved" by its own Burnside sum. The score layer uses the cut space K_0 = GF(2)^{n-1}. The cycle layer uses the cycle space K_1 = GF(2)^{C(n-1,2)}. The interaction is the d_1 differential connecting them.

### Tower 10: Hodge Decomposition (differential geometry)
```
Ω^k = H^k ⊕ dΩ^{k-1} ⊕ d*Ω^{k+1}
```
Every differential form decomposes into harmonic + exact + coexact. The three pieces are orthogonal.

**Isomorphism**: The tournament adjacency matrix decomposes as:
- **Antisymmetric** (so(n)): the tournament direction information
- **Symmetric** (p): the cooperation/conflict structure (Ω)
- **Scalar** (trace): the identity/Rédei quantum

This IS the Cartan decomposition. The Hodge star * corresponds to the complement operation T → T^op. The harmonic forms correspond to SC tournaments (self-complementary = self-dual).

## The Meta-Pattern

All ten towers share:

| Feature | Description | Tournament Version |
|---------|-------------|-------------------|
| **Filtration** | Object splits into graded layers | arc_orbits = cycle_null + cut_free |
| **Quotient** | Each layer is simpler than the previous | S_n → Z_2 → cycle proj → cut proj |
| **Obstruction** | Residual between adjacent layers | residual = cut-cycle interaction |
| **Convergence** | Tower stabilizes / obstructions vanish | identity dominance, res/T → 0 |
| **Generator count** | Controlled by 2 and/or 3 | PSL(2,Z) = ⟨S(2), ST(3)⟩ |

## Why 2 and 3

The deepest reason all these towers are isomorphic: they are all **manifestations of PSL(2,Z)** acting on different mathematical spaces.

- Cayley-Dickson: the doubling (order 2) and the 3-cycle (order 3) structure of hypercomplex multiplication
- Burnside: odd-order permutations (no even cycles = no 2-torsion) with 3-cycle corrections
- Cartan: so(n) vs p split by the involution (order 2), with the 3-fold structure of sl(2) ↪ gl(n)
- Postnikov: π_n kills even information at each step, 3-torsion in homotopy groups controls obstructions
- Chromatic: height n controlled by formal group law, heights 0 and 1 = rational and K-theory = orders ∞ and 2
- RG flow: relevant operators (order 2, survive) vs irrelevant (higher order, die)

The fraction 2/3 = order(S)/order(ST) appears everywhere because it IS the fundamental ratio of the modular group, and the modular group IS the symmetry of the hyperbolic plane, and the hyperbolic plane IS the universal cover of every Riemann surface of genus ≥ 2.

Tournament theory sits at genus 0 (the modular curve), which is the simplest nontrivial case. This is why tournament theory is "the simplest thing that isn't simple" — it has exactly the complexity of PSL(2,Z), no more and no less.
