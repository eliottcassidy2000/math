# The Three Functors: H as a Composition

*opus-2026-04-04-S7*

## The Deepest Simplification

Everything about the Hamiltonian path polynomial H(t) follows from a **three-step composition**:

```
{0,1}^m  ---(1)--->  Tournaments  ---(2)--->  Graphs  ---(3)--->  Z
 tiling       linear        conflict graph    independence poly
   t     ------>  T(t)   ------>  Omega(T)  ------>  I(Omega, 2) = H
```

### Functor 1: Tiling → Tournament (linear, invertible)

Each bit t_k controls one arc. The base path is fixed. This is a coordinate chart on tournament space — a bijection between {0,1}^m and the 2^m tournaments containing a fixed Hamiltonian path.

**This functor is trivial.** It just relabels binary strings as tournaments.

### Functor 2: Tournament → Conflict Graph (nonlinear, local)

T maps to Omega(T), the graph whose vertices are directed odd cycles and whose edges connect vertex-sharing cycles. This is where all the nonlinearity lives.

**Properties of this functor:**
- **Local**: each arc participates in finitely many cycles. Flipping one tile adds/removes a bounded number of Omega-vertices.
- **Monotone at the origin**: from transitive (no cycles), each tile flip creates cycles but never destroys them (since the transitive has none to destroy).
- **Non-monotone in general**: at dense tilings, a tile flip can destroy more cycles than it creates.
- **Source property**: when vertex v is a source (beats all others), no odd cycle uses v as an intermediate vertex with a "receiving" arc. This explains recursive preservation (THM-299).

### Functor 3: Conflict Graph → Integer (polynomial evaluation)

Omega maps to I(Omega, 2) = Σ_k α_k · 2^k where α_k counts independent k-sets.

**This functor is algebraic but global.** It depends on the entire graph structure of Omega, not just local neighborhoods.

## Why Each Theorem Follows

| Theorem | Which functor explains it |
|---------|--------------------------|
| Redei (H odd) | Functor 3: α_0 = 1 always → I ≡ 1 (mod 2) |
| Recursive preservation (THM-299) | Functor 2: source vertex means Omega unchanged |
| Sign rule — same-end (THM-301a) | Functor 2: shared vertex forces opposite cycle directions |
| Sign rule — cross-end (THM-301b) | Functor 2: relay vertex forces same directions |
| Degree cap 2⌊(n-1)/2⌋ | Functor 3: max independent set size in a tournament's cycle graph |
| c_3 predicts 91% of H | Functor 3: I ≈ 1 + 2α_1 ≈ 1 + 2·c_3 + 2·c_5 + ... |
| Antiferromagnetic log-coupling (THM-290) | Functor 2: cycles compete for vertices |
| Grid reflection symmetry (THM-303) | Functor 1: vertex relabeling v↔n+1-v |

## The Tropical Perspective

The real tropicalization trop(H)(t) = max_S log|c_S| reveals the **dominant monomial** at each tiling:

- At small n (n=5): **linear terms dominate** (the apex tile c_(n,1) = 2^{n-2} wins for ~50% of tilings)
- At large n (n=7): **quadratic terms dominate** (disjoint pairs like c_{(6,1),(7,2)} = 52 win for 25%)
- The **tropical shift** from degree 1 to degree 2 reflects the growing importance of cycle interactions as the conflict graph becomes denser.

The 2-adic tropicalization is trivial (Redei: all H are odd), but the real tropicalization reveals the dominant interaction scale.

## The Frustration Threshold

H is NOT monotone under tile flips. Starting from a tiling with many backward tiles, flipping one more backward can DECREASE H. This happens when the negative quadratic interactions (same-end competition) overwhelm the positive linear contribution.

**Critical threshold for tile k**: gradient flips sign when approximately 2^{skip(k)-2} of tile k's same-end neighbors are backward. This is the **frustration threshold** — the point where adding another backward arc hurts rather than helps.

## What Remains Unknown

1. **Why n-2 negative eigenvalues?** The Q_same decomposition accounts for the negative directions, but the exact count n-2 (not n, not n-1) is unexplained.

2. **The disjoint sign rule.** Same-direction cycles outnumber opposite-direction cycles for disjoint pairs, but no clean proof exists.

3. **The gap structure.** H values form an odd arithmetic-like sequence with gaps of 2 and occasional gaps of 4. Why 4?

4. **The independence polynomial at x=2 specifically.** Why is x=2 the right evaluation point? (It comes from the GF(2) nature of arc orientations, but the geometric reason is deeper.)
