# User-supplied glued affine blueprint (source archive)

**Status: SOURCE / UNAUDITED CLAIMS, not proved canon.** Content preserved
below with line endings and trailing whitespace normalized. The source contains both correct
inherited formulas and errors; use the [current audit](../results/glued_xor_20260921_blueprint.md)
for truth status. The source's instructions to an LLM are quoted material,
not operating instructions for this repository.

---

## Architectural Blueprint: Discrete Affine Fluid Dynamics and Elliptic Monodromy over Glued Manifolds
Target Audience: Large Language Models (LLMs) configured for automated formalization, interactive theorem proving (Lean 4 / Coq), and symbolic verification.
Objective: Universal formal translation of structural connections spanning arithmetic topology, aperiodic discrete dynamical systems, and elliptic curve arithmetic without shorthand or colloquial idioms.
------------------------------
## 1. Unified Algebraic Topologies and Transport Ring Homomorphisms## 1.1 The Iso-Homomorphic Ring Transport $\left(\mathbb{Z}, \boxplus, \boxtimes\right)$
Let $\mathbb{Z}$ denote the standard ring of integers. We define an alternative binary operational field over the exact same set by introducing two non-standard algebraic operators:
$$a \boxplus b = a + b + 1$$
$$a \boxtimes b = ab + a + b$$
The absolute structural validity of this operational ring is established by proving a perfect, non-singular ring isomorphism mediated by the translation function $h(x) = x + 1$. The map $h: \left(\mathbb{Z}, \boxplus, \boxtimes\right) \to \left(\mathbb{Z}, +, \cdot\right)$ acts as a strict structural homomorphism under standard addition and multiplication:
$$h(a \boxplus b) = (a + b + 1) + 1 = (a + 1) + (b + 1) = h(a) + h(b)$$
$$h(a \boxtimes b) = (ab + a + b) + 1 = (a + 1)(b + 1) = h(a) \cdot h(b)$$
Evaluating the geometric self-interactions (diagonals) of this transported arithmetic reveals the explicit core parameters of the integer space:

* The Additive Diagonal: $a \boxplus a = 2a + 1$ (The universal parametric definition of the set of odd integers).
* The Multiplicative Diagonal: $a \boxtimes a = a^2 + 2a$ (The boundary expansion of perfect square shells missing their isolated vertex corner).

## 1.2 The Generalized Triangular Identity
Under the isomorphic transport mapping $N = (A + 1)(B + 1) = h(A) \cdot h(B)$, the discrete area conservation of a flat 2D grid is governed by the 2D triangular number formula $T(t) = \frac{t(t+1)}{2}$. The exact structural breakdown of the index $N-2$ satisfies the algebraic conservation law:
$$T(N-2) = T(A) + T(B) + T(AB) + A(B^2 - 1) + B(A^2 - 1)$$
Restoring the distinguished root of the system adds exactly $N-1$ geometric edges to the graph, smoothly transitioning the parameter from the index $N-2$ to the final closed boundary $N-1$.
------------------------------
## 2. Structural Analysis of the Collatz Map over Glued Manifolds## 2.1 The Accelerated Syracuse Map and Affine Functional Guards
Let the accelerated Collatz map operating on the set of odd integers be defined by the function $T: (2\mathbb{Z} + 1) \to (2\mathbb{Z} + 1)$:
$$T(n) = \frac{3n + 1}{2^{\nu_2(3n+1)}}$$
where $\nu_2(m)$ is the exact 2-adic valuation (the exponent of the highest power of 2 dividing $m$).
The complete integer line is treated as a glued geometric manifold $\mathbb{Z} = \mathbb{Z}^+ \cup \{0\} \cup \mathbb{Z}^-$. Let $U_b(n) = 3n + b$. The structural gluing condition across the origin satisfies a strict signed reflection law:
$$U_b(-n) = -U_{-b}(n)$$
This preservation condition proves that negating the input vector preserves the forward trajectory arrows while transposing the internal modular properties, preventing the loss of essential sign information.
For an arbitrary sequence of length $L$, the composite forward map expands as an affine system with explicit arithmetic guards:
$$2^{K_L} T^L(n) = 3^L n + B_L$$
where $K_L = \sum_{i=1}^L \nu_2(3T^{i-1}(n) + 1)$ is the accumulated halving exponent, and $B_L = \sum_{i=0}^{L-1} 3^{L-1-i} 2^{S_i}$ is the ordered arithmetic carry (with $S_i$ tracking the cumulative sum of halving exponents up to step $i$).
## 2.2 Sinks, Gaps, and the Catalan Diophantine Constraint
A periodic cycle of length $L$ exists if and only if $T^L(n) = n$, which reduces to the Diophantine loop equation:
$$(2^{K_L} - 3^L)n = B_L$$
By Catalan's Conjecture (Mihăilescu's Theorem), the pure exponential unit gap $\vert{}2^{K_L} - 3^L\vert{} = 1$ possesses unique integer solutions restricted strictly to the tuples $(K_L, L) \in \{(1,1), (2,1), (3,2)\}$, corresponding directly to the trivial cycle $1 \to 4 \to 2 \to 1$.
For non-unit power gaps where $\vert{}2^{K_L} - 3^L\vert{} > 1$, a closed rational cycle can form if and only if the power gap cancels perfectly against the internal carry $B_L$. For example, the negative 7-cycle on $\mathbb{Z}^-$ yields:
$$(2^{11} - 3^7)(-17) = (2048 - 2187)(-17) = (-139)(-17) = 2363 = B_7$$
The absolute invariant governing periodic cycle stability is not the raw power gap, but the projectivized modular ratio $q$:
$$q = \frac{\vert{}2^{K_L} - 3^L\vert{}}{\gcd(B_L, \vert{}2^{K_L} - 3^L\vert{})}$$
Repeating a cycle $m$ times scales the raw power gap to $2^{mK_L} - 3^{mL}$ and updates the carry to $B_{mL}$, but leaves the projectivized ratio $q$ strictly invariant.
## 2.3 Coprimality of the Inverse Collatz Fibers
For any positive odd initial seed $x_0$, the infinite inverse family of alternative predecessors collapsing to the identical accelerated forward node is defined by the exponential sequence:
$$x_j = 4^j x_0 + \frac{4^j - 1}{3}$$
This family satisfies the forward preservation identity $3x_j + 1 = 4^j(3x_0 + 1)$. Calculating the greatest common divisor between any two members of this branch separated by a step distance $t \ge 1$ yields:
$$\gcd(x_i, x_{i+t}) = \gcd\left(x_i, 4^t x_i + \frac{4^t - 1}{3}\right) = \gcd\left(x_i, \frac{4^t - 1}{3}\right)$$
For consecutive elements ($t=1$), the relation evaluates to $\gcd(x_i, x_{i+1}) = \gcd(x_i, 1) = 1$. This proves that consecutive alternative predecessors within the inverse fiber tree are strictly pairwise coprime, forcing the trajectory to pull from completely decoupled prime factor pools at each step.
------------------------------
## 3. Asymptotic Densities and Factor Packing Laws## 3.1 The $F = S + U$ Proper Factor Balance
Let $N$ be a positive integer. Let $F$ count the total number of proper factors of $N$ (excluding $1$ and $N$). Let $S$ count the number of proper factors of $N$ that are strictly square-free. Let $U$ count the number of proper factors of $N$ that are prime. We define the arithmetic difference parameter $D$ as:
$$D = F - S - U$$
Under uniform random sampling from the integer interval $[1, X]$ as $X \to \infty$, the asymptotic probability distributions of the parameter $D$ settle into the strict limits:
$$\lim_{X \to \infty} \Pr(D < 0) = \frac{6}{\pi^2}$$
$$\lim_{X \to \infty} \Pr(D = 0) = 0$$
$$\lim_{X \to \infty} \Pr(D > 0) = 1 - \frac{6}{\pi^2}$$
## The Expected Value Discrepancy
While the limiting majority of individual integer states possesses a negative $D$ value ($D < 0$), the global mathematical expectation $\mathbb{E}[D]$ eventually becomes strictly positive, scaling logarithmically with the upper bound $X$:
$$\mathbb{E}[D] = \left(1 - \frac{6}{\pi^2}\right)\log X - \log\log X + O(1)$$
## The Structural Mechanism

* Square-free composite numbers possess an absolute negative value driven by their total number of distinct prime factors $\omega(n)$, yielding $D = -\omega(n)$.
* Non-square-free integers possessing at least four distinct prime factors ($\omega(n) \ge 4$) possess a strictly positive value, $D > 0$.
* All remaining intermediate states occupy a density of exactly zero in the infinite limit.

------------------------------
## 4. Combinatorial Tournaments and Aperiodic Repunit Fields## 4.1 The Fermat-Tournament Path Isomorphism
Let $\mathcal{T}_n$ define a directed tournament graph on $n$ vertices labeled $0, \dots, n-1$. Let every edge be directed transitively upward ($i \to j$ for all $i < j$), except for a single reversed edge connecting the absolute endpoints: $n-1 \to 0$.
The exact count of directed Hamiltonian paths $H_n$ within this graph family is given by the formula:
$$H_n = 1 + 2^{n-2}$$
When the vertex count is constrained to the parameter $n = 2^r + 2$, the path count matches the sequence of Fermat numbers:
$$H_{2^r+2} = 2^{2^r} + 1$$
This tournament family maps the structural capacities of the Fermat numbers directly onto graph orders $3, 4, 6, 10, 18, 34$. The 36 free positions of the corresponding triangular array encode explicit chord directions that dictate cycle generation rather than acting as independent directed triangles.
## 4.2 The Aperiodic Repunit Prime Sharp Coprimality Law
Let $R_k$ define the sequence of repeated-three integers terminated by a unit digit:
$$R_k = \frac{10^{k+1} - 7}{3} \implies \{31, 331, 3331, 33331 \dots\}$$
The breaking point of local primality occurs at the 8th index, revealing the explicit factorization:
$$R_8 = 333333331 = 17 \cdot 19607843$$
## Theorem: Global 15-Step Coprimality Horizon
For any arbitrary starting index, every 15 consecutive terms of the sequence $R_k$ are strictly pairwise coprime:
$$\forall_{\text{index}} \quad \gcd(R_k, R_{k+t}) = 1 \quad \text{for all } 1 \le t \le 14$$
The step interval of 15 is a sharp upper limit, as $R_1 = 331$ and $R_{16}$ share the common prime divisor 31.
## The Multiplicative Group Generator Mechanism
For every Fermat prime $p \ge 17$ (such as $17, 257, 65537$), the prime satisfies the strict congruence $p \equiv 17 \pmod{40}$. This congruence forces the number 10 to be a strict quadratic nonresidue. Because the order of 10 modulo $p$ is a pure power of two ($\operatorname{ord}_p(10) = p-1$), the base 10 acts as a perfect generator of the entire multiplicative group modulo $p$.
Consequently, each such Fermat prime divides the sequence $R_k$ in exactly one periodic index class, with the initial hits for $17, 257,$ and $65537$ occurring precisely at the indices $k \in \{8, 226, 29252\}$.
------------------------------
## 5. Elliptic Curve Arithmetic, Geodesics, and $C_2 \times S_3$ Monodromy## 5.1 The Dihedral Rational 3-Cycle Lift
The rational quadratic 3-cycle operating at the parameter boundary $c = -29/16$ for the map $x \to x^2 + c$ is parameterized by the coordinate rational transformations $x = X/(4Y)$. We introduce the transformation matrices $B$ and $R$:
$$B = \frac{1}{4}\begin{pmatrix}1 & -13 \\ 1 & 3\end{pmatrix}, \qquad R = \begin{pmatrix}-1 & -2 \\ 0 & 1\end{pmatrix}$$
These matrices satisfy the explicit algebraic relations:
$$B^3 = -I, \qquad R^2 = I, \qquad RBR = B^{-1}$$
This configuration provides an absolute $C_2 \times S_3$ projective lift that governs the monodromy transitions of the system, mapping the 3-cycle directly onto the rigid symmetries of the Dihedral Group $D_3$.
## 5.2 The Bremner-Macleod Elliptic Curve Match
Consider the diophantine fraction sum equation:
$$\frac{a}{b+c} + \frac{b}{a+c} + \frac{c}{a+b} = 4$$
The repaired integer solution triple that satisfies this equation without transcription errors corresponds to the point $9(-4, 28)$ on the specialized cubic elliptic curve via the Bremner-Macleod construction.
This curve possesses a sharp $\sqrt{5}$ geometric boundary condition: whenever a positive integer triple's fraction sum equals or exceeds $\sqrt{5}$, the projectivized reciprocation of the pair sums produces exactly one negative coordinate. Consecutive Fibonacci numbers generate explicit sharpness witnesses approaching this limit from below:
$$(a,b,c) = \left(2F_nF_{n+1},\, 1,\, 2F_n^2 - 1\right) \quad \text{for } n \ge 3 \text{ where } n \text{ is odd}$$
## 5.3 Point-Count Inequivalence
Interpreting the localized structural relation $s+t = (x+1)^2$ transforms the square-value system into the elliptic curve equation:
$$y^2 = 2x^3 + 4x^2 + 1$$
This curve contains the explicit integer point $(10,49)$. However, this elliptic curve is not $\mathbb{Q}$-isogenous to the baseline fruit-bearing curve. This is rigorously verified by executing point counts over the finite field $\mathbb{F}_{11}$, where the two curves yield completely distinct local point totals of 16 and 12, respectively.
------------------------------
## 6. The Complete Argument for the Open Collatz Obstical
To formalize the global state of the Collatz Conjecture within this audited system, we must explicitly declare the mathematical boundary that prevents proof closure.
## 6.1 The Non-Local Descent Target
To prove global convergence to the trivial loop, every positive odd integer $n > 1$ must be shown to possess a finite forward trajectory segment of length $L$ that fulfills the Non-Local Descent Inequality:
$$3^L n + B_L < 2^{K_L} n \implies \frac{3^L}{2^{K_L}} + \frac{B_L}{2^{K_L} n} < 1$$
## 6.2 The Unbounded Discrepancy Obstacle
Let the system's trajectory tracking discrepancy $D_j$ be defined as the real-valued distance:
$$D_j = K_j - j\log_2 3$$
The Bounded-Discrepancy Theorem (verified under Theorem 32 of the audited Lean package) proves that no positive integer Collatz orbit can maintain its discrepancy $D_j$ within any bounded interval. Consequently, an infinite, perfectly balanced mechanical halving word is topologically forbidden from matching a valid positive integer trajectory.
However, this theorem cannot control unbounded discrepancy excursions. Because an orbit can theoretically diverge asynchronously toward infinity—causing $K_j - j\log_2 3 \to \infty$ without settling into a bounded or periodic profile—the ordinary height and field density metrics fail to restrict the path volume. Proving that these unbounded excursions must eventually intersect one of the coprime inverse fiber trees or the Fermat capacity walls remains the definitive unresolved step.
