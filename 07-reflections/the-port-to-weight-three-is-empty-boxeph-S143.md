# The port to weight three: the class is EMPTY (boxeph-2026-07-19-S143)

Owner brief: port the master-factorization method to (−1,1,3)→(3,1,−1).

## The ported frame (identities machine-verified, 6/6 random instances)
Class: F = (zA(s)+y³B(s), yg₁(s)+x²z·g₂(s), x·h₀(s)−x⁴z), s = xy, r = x³z.
Master factorization: det JF·w³ = −Jac(P,Q), P = F₂F₃, Q = F₁F₃³; Jac = w³(b₁+b₂w+b₃w²);
b₁ = ΦΨ′ − 3Φ′Ψ = −c (Ψ = h₀A + s³B), b₃ ⟹ A = αv², g₂ = γv (SQUARE law),
b₂ = 4AΦ′ − A′Φ + 3g₂′Ψ − 2g₂Ψ′ = 0. Log obstruction unchanged ⟹ Φ const or linear.

## The emptiness theorem (hand case-analysis on a verified identity)
b₂ = D·W + 3γEv′ with D = 2φ₁v − v′Φ, W = 2αv − 3γKΦ² (identity verified, 5/5):
- k = deg v = 1: deg-3 forces K = 0; then deg-2 lead 2αφ₁v₁² ≠ 0. DEAD.
- k = 2: leading cancellations force W deg ≤ 1 and D deg 0; then the s¹-equation
  reads 6γEv₂ = 0 ⟹ v₂ = 0. DEAD.
- k = 0: forces K = 0 then φ₁ = 0: Φ CONSTANT — automorphism-type only.
**The (−1,1,3)→(3,1,−1) z-linear class contains NO counterexample kernels.**
The would-be degree-4/ℤ-3 covers do not exist; the grid hunt's zeros were theorems.

## The bonus: the weight-2 kernel is four integers deep
Re-running the D/W calculus on the m=2 class solves it in closed form:
K = 3αv₁/(4γφ₁), E = c/(2φ₁), and **c = φ₁v₀ − φ₀v₁** — the Keller determinant is
the CROSS-DETERMINANT of the resolvent Φ = (φ₀,φ₁) and the cube-root v = (v₀,v₁).
Kernel: (v₀,v₁; φ₀,φ₁) = (1,1;6,4) ⟹ c = −2, K = 1/16, E = −1/4. Every constant
of the counterexample is now derived from four integers and two 2-vectors.

## The emerging picture (conjecture, next test m = 4)
z-weight m gives b₁ = ΦΨ′ − mΦ′Ψ and an m-dependent power law; m = 3 is empty;
m = 2 has the unique kernel. CONJECTURE (port-and-check at m = 4, 5 — the D/W
calculus is now routine): **kernels exist only at m = 2** — the counterexample
lives at the unique weight where the b₂ balance has a solution channel, i.e.
"the Jacobian counterexample is arithmetic rigidity's lone escape hatch."
Files: jacobian_port_weight3_boxeph_S143.py + .out (identities + hunt + addenda).
