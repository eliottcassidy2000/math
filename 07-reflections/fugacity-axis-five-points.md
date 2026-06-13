# The Five-Point Fugacity Axis

**Session:** oracle-2026-05-17-S1
**Source:** `fugacity_axis_full.py`, `05-knowledge/results/fugacity_axis_full.out`
**Depends on:** `fugacity-axis-and-vanishing-theorem.md`, `cubic-family-x6-coloring.md`

---

## The Five Special Fugacities

The independence polynomial $I(\Omega(T), x)$ evaluated at five canonical points:

| Fugacity | Evaluation | Interpretation |
|---|---|---|
| $x=-1$ | $I(-1) = 1-\alpha_1+\alpha_2-\cdots$ | Negative reduced Euler characteristic: $I(-1) = -\tilde{\chi}(\Delta\Omega)$ |
| $x=0$ | $I(0) = 1$ | Empty collection (trivial) |
| $x=1$ | $I(1) = 1+\alpha_1+\alpha_2+\cdots$ | Total # vertex-disjoint odd cycle collections |
| $x=2$ | $I(2) = H(T)$ | Hamiltonian paths (the OCF) |
| $x=6$ | $I(6)$ | $S_3$-decorated cycle collections (triangle colorings) |

**All five are linearly determined by $(\alpha_1, \alpha_2)$ for $n \leq 8$** (where $\alpha_k=0$ for $k \geq 3$).

---

## The Topology Point: $I(\Omega, -1)$

**Theorem.** $I(\Omega(T), -1) = -\tilde{\chi}(\Delta\Omega)$, where $\tilde{\chi}$ is the reduced Euler characteristic of the **independence complex** $\Delta\Omega$ (the simplicial complex whose $k$-simplices are $(k+1)$-element independent sets of the conflict graph $\Omega(T)$).

**Proof.** The simplicial complex $\Delta\Omega$ has exactly $\alpha_{k+1}$ faces of dimension $k$. Its Euler characteristic is:
$$\chi(\Delta\Omega) = \sum_{k\geq 0} (-1)^k \alpha_{k+1} = -\sum_{j\geq 1}(-1)^j \alpha_j = -(I(\Omega,-1)-1) = 1-I(\Omega,-1)$$
Since $\tilde{\chi} = \chi - 1$, we get $\tilde{\chi} = -I(\Omega,-1)$. $\square$

**Consequence:** The reduced Euler characteristic of the tournament independence complex is:
$$\tilde{\chi}(\Delta\Omega(T)) = \alpha_1 - \alpha_2 + \alpha_3 - \cdots$$

For $n \leq 8$ (degree $\leq 2$): $\tilde{\chi} = \alpha_1 - \alpha_2 - 1 + 1 = \alpha_1 - \alpha_2$.

Wait, let me restate: $\tilde\chi = -I(\Omega,-1) = -(1-\alpha_1+\alpha_2) = \alpha_1-\alpha_2-1$.

**Verified exhaustively** at $n=6$ (47 iso classes, zero errors).

### Topological spectrum (n=6):

| Class | $\alpha_1$ | $\alpha_2$ | $I(-1)$ | $\tilde{\chi}$ | topology |
|---|---|---|---|---|---|
| Transitive | 0 | 0 | 1 | -1 | $\tilde{\chi}=-1$ (empty complex) |
| Any 1-cycle | 1 | 0 | 0 | 0 | $\tilde{\chi}=0$ (contractible) |
| Double root | 2 | 1 | 0 | 0 | $\tilde{\chi}=0$ (path = contractible) |
| SC max | 20 | 1 | -18 | 18 | 18 "holes" in topology |

**The most SC tournament (H=45, $\alpha_1=20$, $\alpha_2=1$) has the most topologically complex independence complex**, with $\tilde{\chi}=18$.

---

## Extraction Formulas via $I(-1)$ and $I(1)$

The pair $(I(-1), I(1))$ alone determines all of $(\alpha_1, \alpha_2)$:

$$\boxed{\alpha_1 = \frac{I(\Omega,1) - I(\Omega,-1)}{2}}$$

$$\boxed{\alpha_2 = \frac{I(\Omega,1) + I(\Omega,-1)}{2} - 1}$$

**Proof:** $I(1)+I(-1) = (1+\alpha_1+\alpha_2)+(1-\alpha_1+\alpha_2) = 2+2\alpha_2$ and $I(1)-I(-1)=2\alpha_1$. $\square$

**Verified:** Exact (error $< 10^{-9}$) for all 47 n=6 iso classes.

### The H formula from I(-1) and I(1):

For $n \leq 8$ (degree $\leq 2$):
$$\boxed{H = 3 \cdot I(\Omega,1) + I(\Omega,-1) - 3}$$

**Proof:** $3I(1)+I(-1)-3 = 3(1+\alpha_1+\alpha_2)+(1-\alpha_1+\alpha_2)-3 = 1+2\alpha_1+4\alpha_2 = H$. $\square$

This says: **H equals 3 times the "counting" evaluation plus the "topology" evaluation, minus 3.**

For $n \geq 9$ (where $\alpha_3 > 0$): the formula has correction term $+6\alpha_3+12\alpha_4+\cdots = \sum_{k\geq 3}(2^k-2)\alpha_k$.

---

## The "Odd-Even" Symmetry

The independence polynomial has an odd-even parity structure:

| | Formula | Measures |
|---|---|---|
| $I(1)+I(-1)$ | $2+2\alpha_2$ | $2 \times$ (disjoint cycle pairs + 1) |
| $I(1)-I(-1)$ | $2\alpha_1$ | $2 \times$ (odd cycle count) |

The **sum** depends only on $\alpha_2$ (the "pairing" content) and the **difference** only on $\alpha_1$ (the "cycle" content). This is a parity decomposition of the tournament structure.

---

## The Root Ratio Formula (Exact)

For degree-2 polynomials with roots $\rho_1 \geq \rho_2 > 0$:
$$r = \frac{\rho_2}{\rho_1} = \frac{\alpha_1 - \sqrt{\alpha_1^2-4\alpha_2}}{\alpha_1 + \sqrt{\alpha_1^2-4\alpha_2}}$$

**Asymptotic expansion** for large $\alpha_1/\alpha_2$:
$$r \approx \frac{\alpha_2}{\alpha_1^2} \quad \text{(leading term)}$$

For $\alpha_2=1$: $r \approx 1/\alpha_1^2$. At $\alpha_1=20$ (SC H-max): $r \approx 1/400 = 0.0025$. Exact: $0.00251$ (0.4% error). ✓

**The H-ratio formula** ($\alpha_2=1$):
$$r \approx \frac{1}{\alpha_1^2} = \frac{4}{(H-5)^2} \quad \text{for degree-2, } \alpha_2=1$$

Since $H = 5+2\alpha_1$ for $\alpha_2=1$ classes, the most asymmetric roots correspond to **largest $H$** within the $\alpha_2=1$ family, i.e., the SC H-maximizer.

---

## Alpha-1 Gap at n=7

From 5000 random samples, **degree-1 ($\alpha_2=\alpha_3=0$) achievable $\alpha_1$ values at n=7:**
$$\{0,1,2,4,5,6,7,8,9,11,12,13,15,16,18,19,22,23,24,25,34\}$$

**Gaps at n=7 (degree-1):** $\{3, 10, 14, 17, 20, 21, 26, \ldots\}$

The **permanent gap at $\alpha_1=3$** persists from n=6 to n=7. All other gaps are likely n-specific. The key theorem: $\alpha_1=3$ with $\alpha_2=0$ is IMPOSSIBLE for ANY $n$ (this is the forbidden root $r=-1/3$).

Note: $\alpha_1=3$ CAN occur with $\alpha_2 \geq 2$ at $n=7$ (from prior computations: 30 tournaments at n=7 with $\alpha_1=3$, $\alpha_2=2$). The permanent forbidden case is specifically ($\alpha_1=3$, $\alpha_2=0$).

---

## Summary Table: All n=6 Iso Classes (47 found by key)

The full table $((\alpha_1, \alpha_2), H, I(-1), I(1), I(6), \text{SC})$ shows:

**SC tournaments have:**
- Most negative $I(-1)$ (most complex independence complex topology)
- Largest $I(1)$ at each H level (most total collections)
- Smallest root ratio $\rho_2/\rho_1$ (most asymmetric)

**The transitive tournament has $I(-1)=1$** (the unique class with positive Euler characteristic among degree-0 classes). All others have $I(-1) \leq 0$ (non-positive Euler characteristic).

**The double root class ($\alpha_1=2$, $\alpha_2=1$, $H=9$, SC) is topologically special**: $I(-1)=0$ means the independence complex is homotopy equivalent to a point or odd-dimensional sphere.

---

## Open Questions

1. **At n=9 (degree-3 polynomials):** What is the topological interpretation of the degree-3 classes? The formula $H = 3I(1)+I(-1)-3+6\alpha_3$ needs more evaluations.

2. **Higher independence complex topology:** Is $\Delta\Omega(T)$ always homotopy equivalent to a wedge of spheres? If so, $H_k(\Delta\Omega) = 0$ for all $k$ except one, and $\tilde{\chi} = \pm \text{rank}(H_k)$.

3. **The most topologically complex tournament:** As $n \to \infty$, what is the maximum $\tilde{\chi}(\Delta\Omega) = \alpha_1-\alpha_2$ over all $n$-vertex tournaments? This is the maximum of (# odd cycles) - (# disjoint cycle pairs), which has a concrete combinatorial meaning.

4. **Shelling of $\Delta\Omega$:** Is the independence complex of tournament conflict graphs always shellable? Shellable complexes have $\tilde{H}_k = 0$ for $k < \dim$, simplifying the topology.
