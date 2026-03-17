# Substack Hooks

Last updated: 2026-03-17.

---

## For Mathematicians

### Hook A: "The Cayley transform counts lattice paths"

The Taylor coefficients of $\left(\frac{1+x}{1-x}\right)^m$ count diagonal steps in Delannoy paths. Specifically, if $g_k(m) = \sum_{j} \binom{k-1}{j-1}\binom{m}{j}2^{j-1}$, then $\left(\frac{1+x}{1-x}\right)^m = 1 + 2\sum_{k \geq 1} g_k(m)\, x^k$, the matrix $T(k,m) = k\,g_k(m)$ is symmetric with diagonal OEIS A108666, and the Fourier energy of Hamiltonian path counts over random tournaments decomposes as $\mathrm{CV}^2 = \sum_k 2g_k(n{-}2k)/(n)_{2k} = 2/n + O(1/n^3)$, with exact cancellation at $1/n^2$.

---

### Hook B: "arctanh is the unique odd function with rational exponential"

$\mathrm{arctanh}(x) = x + x^3/3 + x^5/5 + \cdots$ is the unique odd formal power series $f$ such that $e^f$ is rational. The rational function is $(1+x)/(1-x)$. This forces $\mathrm{arctanh}$ to govern any system of sequential directed binary comparisons — including tournaments, spin chains, phylogenetic distances, and channel capacities.

---

### Hook C: "A new integer sequence from tournaments"

$W(n) = 1, 2, 8, 32, 158, 928, 6350, 49752, 439670, \ldots$ counts NUD permutations weighted by $2^{\text{unit ascents}}$. Not in the OEIS. Satisfies $W(n)/n! = 1 + \mathrm{CV}^2(H)$ where $H$ = Hamiltonian paths in random tournaments. The generating function involves $((1+x)/(1-x))^m$ and the $x$-tribonacci recurrence.

---

### Hook J: "Walsh-Fourier analysis compresses tournament invariants 100,000x"

The number of consistent rankings in a round-robin tournament is a function of $2^{C(n,2)}$ binary inputs. Its Walsh-Fourier transform is *almost entirely zero*. At $n=5$: only 3 nonzero amplitudes out of 1,024 (341x compression). At $n=7$: ~20 out of 2,097,152 (~100,000x). The nonzero coefficients are indexed by edge-disjoint even-length path unions, with closed-form amplitudes $|\hat{H}[S]| = 2^r (n-2k)!/2^{n-1}$. This is exact, not approximate — the ranking count is a function of astronomically high "dimension" that secretly lives on a tiny subspace.

---

### Hook K: "Tournaments never have bubble-shaped holes"

In GLMY path homology (the directed analogue of simplicial homology for graphs), the second Betti number $\beta_2 = 0$ for *every* tournament. Proved by induction + long exact sequence (THM-108/109). This vanishing is specific to tournaments — general directed graphs commonly have $\beta_2 > 0$. Meanwhile, $\beta_1 \in \{0,1\}$, $\beta_3$ can reach 2 at $n=8$, and $\beta_1 \cdot \beta_3 = 0$ (a seesaw: loop-holes and higher-dimensional holes are mutually exclusive). The mechanism: all $\beta_2 > 0$ graphs have twin vertices (identical neighborhoods), and tournament completeness forbids twins.

---

### Hook L: "Paley tournaments are Heisenberg Lie algebras in disguise"

The Paley tournament $T_p$ ($p \equiv 3 \pmod{4}$) has Betti numbers $\beta_m = m(m-3)/2$ and $\beta_{m+1} = \binom{m+1}{2}$, where $m = (p-1)/2$. These are *exactly* the Betti numbers of the Heisenberg Lie algebra $\mathfrak{h}_m$ (Santharoubane 1983). The Euler characteristic equals $p$. Per-eigenspace: each of the $p-1$ nonzero cyclic eigenspaces contributes exactly $\beta_{m+1} = 1$; all of $\beta_m$ comes from the zero eigenspace (diff-seq complex). No prior work connects tournament path homology to Heisenberg cohomology.

---

### Hook M: "The signed permanent knows about binary digit sums"

Define $S(T) = \sum_P \prod B[P_i, P_{i+1}]$ where $B = 2A - J$ is skew-symmetric. Then $S(T) \equiv 0 \pmod{2^{n-1}}$ for all tournaments — universally. But full universality (independence from $T$ modulo $2^{n-1}$) holds if and only if $s_2(n-3) \leq 1$, where $s_2$ is the binary digit sum. Universal at $n = 3, 5, 7, 11, 19, 35, 67, \ldots$ First failure: $n = 9$ ($s_2(6) = 2$), where $S \bmod 256$ depends on the parity of the 3-cycle count. A Kummer-type condition controlling congruences of a combinatorial permanent.

---

### Hook N: "k-nacci traces are Mersenne numbers"

For the $k$-nacci companion matrix $M_k$ (the $k \times k$ matrix with 1's on the superdiagonal and bottom row): $\mathrm{Tr}(M_k^n) = 2^n - 1$ for all $1 \leq n \leq k$. At $n = k$, this is the $k$-th Mersenne number. At $n = k+1$, the identity breaks by exactly $k+1$. Proof: induction on Newton's identities with all-$(-1)$ coefficients. While $\mathrm{Tr}(M_3^3) = 7$ coincides with the first forbidden tournament value, the Mersenne numbers 31, 63, 127 are all achievable — the connection to forbidden values is through the specific combinatorial mechanism of cycle-forcing, not through Mersenne numbers per se.

---

### Hook O: "Unit fraction splitting is controlled by cyclotomic fields"

$3/N = 1/b + 1/c$ is solvable in positive integers iff $N$ has a prime factor $\equiv 2 \pmod{3}$. Equivalently: unsolvable iff all prime factors split in $\mathbb{Z}[\omega]$ (Eisenstein integers). For general prime $k$: $k/p = 1/b + 1/c$ solvable iff $p \equiv -1 \pmod{k}$ (order 2 in $(\mathbb{Z}/k\mathbb{Z})^*$). The unsolvable fraction among primes is $(k-2)/(k-1)$. Master criterion: $k/N$ splits iff $N^2$ has a divisor $d \equiv -N \pmod{k}$. The Erdos-Straus conjecture ($4/n = 1/x + 1/y + 1/z$) reduces via base-42 covering to 4 residue classes mod 42, verified to $10^6$ with 0 failures.

---

## For Everyone

### Hook D: "How to measure whether a ranking is real"

You run a round-robin tournament — every team plays every other. Some results make sense (the best team beats everyone). Others are contradictory (A beats B, B beats C, C beats A). How do you measure whether the final ranking is *real* or just noise? We found the exact answer: the formula $\mathrm{CV} = \sqrt{2/n}$, where $n$ is the number of teams. For 20 teams, any ranking within 32% of random is just luck. Above that, the differences are real.

---

### Hook E: "Pi from a coin flip"

Flip a coin for every pair of players in a tournament: heads = A wins, tails = B wins. Count the consistent rankings. The math that describes the variance of this count is $\mathrm{arctanh}(x) = x + x^3/3 + x^5/5 + \cdots$ — a function built from the odd numbers 1, 3, 5, 7. Now evaluate it at $x = \sqrt{-1}$. Out comes $\pi/4$. The same odd numbers 1, 3, 5, 7 that describe tournament noise also describe the geometry of the circle. The difference: tournament signs are all positive (divergent, infinite — time), while circle signs alternate (convergent, finite — space).

---

### Hook F: "The price of memory"

The golden ratio $\phi = 1.618\ldots$ governs systems that forget instantly. The tribonacci constant $\tau = 1.839\ldots$ governs systems that remember one step back. The difference $\tau - \phi = 0.221\ldots$ is the price of memory — the exact cost, in growth rate, of keeping track of what just happened. Tournaments pay this price because each comparison echoes into the next. Without memory: Fibonacci. With one step of memory: tribonacci. The entire theory of pairwise comparison lives in that gap of 0.221.

---

### Hook G: "The golden ratio hides in tournament physics"

Every matrix has special values where its eigenvalues collide — physicists call them *exceptional points*. The $3 \times 3$ transfer matrix that governs tournament path counting has discriminant $\Delta(x) = 4x(x^2 - 11x - 1)$. Set it to zero: the exceptional-point eigenvalues are $1/\phi$ and $-\phi$, where $\phi$ is the golden ratio. The same number that governs sunflower spirals and Fibonacci rabbits marks the exact boundary where tournament dynamics change character.

---

### Hook H: "Your audio filter is secretly counting lattice paths"

Every digital audio system uses the *bilinear transform* — a formula from the 1940s that converts analog filters to digital ones. The formula is $Q(x) = (1+x)/(1-x)$. We discovered that its Taylor coefficients are *exactly* Delannoy numbers — the count of lattice paths that can go right, up, or diagonally. Nobody noticed this before. Your Spotify equalizer, your hearing aid, your noise-cancelling headphones: they are all, at the level of the math, counting paths on a grid.

---

### Hook I: "One matrix, three constants"

Take one $3 \times 3$ matrix. Evaluate it at three different inputs. Out come three of the most famous constants in mathematics:

- At $x = 1$: the **tribonacci constant** $\tau = 1.839\ldots$ (the growth rate of 1, 1, 1, 3, 5, 9, 17, ...)
- At the exceptional point: the **golden ratio** $\phi = 1.618\ldots$ (where eigenvalues collide)
- At $x = i = \sqrt{-1}$: **pi** via $\mathrm{arctanh}(i) = i\pi/4$

One matrix. Three constants. Each one emerges at a different scale: $\tau$ from counting, $\phi$ from coalescence, $\pi$ from rotation.

---

### Hook P: "Two numbers that can never be the answer"

No matter how you set up a round-robin tournament — 5 teams, 50 teams, 5 million teams — you can never get exactly **7** or **21** consistent rankings. Every other odd number is achievable (eventually, for some tournament size), but 7 and 21 are permanently forbidden.

Why? 7 is the trace of the tribonacci matrix at power 3 — the Newton identity $p_3 = e_1^3 - 3e_1e_2 + 3e_3 = 1+3+3 = 7$ — and this is the same for ALL $k$-nacci matrices with $k \geq 3$, making 7 *universally* forbidden. And 21 = 3 × 7, the product of the cycle obstruction and the Fano obstruction. $42 = 2 \times 3 \times 7$ is the denominator of the Bernoulli number $B_6 = 1/42$ (Von Staudt-Clausen theorem), encoding exactly the primes whose period divides 6.

Even stranger: $42 = 2 \times 3 \times 7$ — the product of parity, cycles, and prohibition — is the number that controls the base-42 covering system for the Erdos-Straus conjecture. The answer to life, the universe, and everything is also the modular base that organizes Egyptian fraction decompositions.

---

### Hook Q: "The topology of who beats whom"

Take any group where everyone competes against everyone else — chess players, sports teams, AI models being ranked. Connect the winners and losers with arrows. Now ask: what is the *shape* of this network?

It turns out tournaments have a very specific shape. They never have "bubbles" (proved: $\beta_2 = 0$). They can have loops ($\beta_1$) or higher-dimensional voids ($\beta_3$), but never both at the same time — a mathematical seesaw.

For the most balanced possible tournament (the Paley tournament, built from number theory), the topology concentrates in a single dimension, and the number of holes equals $m(m-3)/2$ — the same formula that counts diagonals of a polygon. The Euler characteristic equals the number of players.

This gives a new way to fingerprint ranking systems. Two tournaments with the same win-loss records can have completely different topological shapes — revealing hidden structural differences invisible to traditional statistics.

---

### Hook R: "42 is not arbitrary"

Douglas Adams chose 42 as the answer to life, the universe, and everything. Here's what he might not have known:

- $42 = 2 \times 3 \times 7$ — the product of orientation, the smallest cycle, and the smallest prohibition
- The double factorial $(n-2)!!$ — a product that appears in the Walsh spectrum of tournament rankings — freezes at **21 mod 42** for all large $n$. The fixed point is half of 42.
- The Von Staudt chain $1 \to 2 \to 6 \to 42 \to 1806 \to 1806$ (each number is the product of primes $p$ where $p-1$ divides it) reaches a fixed point at 1806, and 42 is the last step before convergence
- The Bernoulli number $B_6 = 1/42$ — the first "interesting" Bernoulli denominator
- $42 = 2 \times 3 \times 7$ is the unique squarefree solution to $\sigma(n) = 2n + \phi(n)$
- The base-42 covering system reduces the Erdos-Straus conjecture (a 78-year-old open problem about Egyptian fractions) to just 4 hard cases

42 isn't random. It's the smallest number that encodes the three fundamental obstructions in tournament parity: **you can flip things** (2), **things can cycle** (3), and **some things are forbidden** (7).

---

### Hook S: "What 1729 knows about AI rankings"

The Hardy-Ramanujan number 1729 — famous as the smallest number expressible as the sum of two cubes in two ways ($12^3 + 1^3 = 10^3 + 9^3$) — appears naturally in tournament theory.

The Paley tournament on 11 vertices (built from quadratic residues mod 11) has 95,095 consistent rankings. Divide by the symmetry group size (55): you get **1729**.

This isn't a coincidence. $1729 = 7 \times 13 \times 19$ — all primes that *split* in the Eisenstein integers. The factorization encodes the decomposition of the tournament's algebraic structure under complex multiplication.

Why does this matter for AI? Modern AI ranking systems (Chatbot Arena, Elo ratings) are secretly computing tournament invariants. When you rank language models by pairwise human preference, you're building a tournament. The OCF formula says the number of consistent total rankings equals $1 + 2\alpha_1 + 4\alpha_2 + \ldots$ where $\alpha_k$ counts collections of non-overlapping preference cycles. Those cycles are the *Condorcet paradoxes* in your preference data — places where Model A beats B, B beats C, and C beats A. Our formulas count them exactly and tell you which cycles make your ranking unreliable.

---

### Hook T: "The fractal hiding in binary arithmetic"

Take the Walsh spectrum of tournament rankings and measure its "binary complexity" — the total 2-adic weight across all spectral components. The answer decomposes into two parts:

$$\text{Total weight} = k^2 + \text{(cumulative binary digit sum)}$$

The first part ($k^2$) is smooth and predictable. The second part is a **fractal** — it grows like $k \log k / 2$ on average but fluctuates according to the binary representation of $k$. The second differences of this sequence count **binary carries**: how many 1-bits cascade when you increment a counter.

The same fractal appears in Kummer's theorem (which binomial coefficients are divisible by a prime), in the Stern-Brocot tree (which fractions have simple continued fractions), and in the distribution of prime factors in factorials. Tournament spectral theory doesn't just use binary arithmetic — it *is* binary arithmetic, viewed through a Fourier lens.

---

## Recommendations

*See end of file for which hooks to prioritize and why.*
