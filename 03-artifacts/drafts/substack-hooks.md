# Substack Hooks

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

## For Everyone

### Hook D: "How to measure whether a ranking is real"

You run a round-robin tournament — every team plays every other. Some results make sense (the best team beats everyone). Others are contradictory (A beats B, B beats C, C beats A). How do you measure whether the final ranking is *real* or just noise? We found the exact answer: the formula $\mathrm{CV} = \sqrt{2/n}$, where $n$ is the number of teams. For 20 teams, any ranking within 32% of random is just luck. Above that, the differences are real.

---

### Hook E: "Pi from a coin flip"

Flip a coin for every pair of players in a tournament: heads = A wins, tails = B wins. Count the consistent rankings. The math that describes the variance of this count is $\mathrm{arctanh}(x) = x + x^3/3 + x^5/5 + \cdots$ — a function built from the odd numbers 1, 3, 5, 7. Now evaluate it at $x = \sqrt{-1}$. Out comes $\pi/4$. The same odd numbers 1, 3, 5, 7 that describe tournament noise also describe the geometry of the circle. The difference: tournament signs are all positive (divergent, infinite — time), while circle signs alternate (convergent, finite — space).

---

### Hook F: "The price of memory"

The golden ratio $\phi = 1.618\ldots$ governs systems that forget instantly. The tribonacci constant $\tau = 1.839\ldots$ governs systems that remember one step back. The difference $\tau - \phi = 0.221\ldots$ is the price of memory — the exact cost, in growth rate, of keeping track of what just happened. Tournaments pay this price because each comparison echoes into the next. Without memory: Fibonacci. With one step of memory: tribonacci. The entire theory of pairwise comparison lives in that gap of 0.221.
