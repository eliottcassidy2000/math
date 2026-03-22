# The Tournament Alphabet

**Session:** kind-pasteur-2026-03-21-S18d
**Arising from:** tournaments-as-codes.md, why-lambda-equals-two.md, the binary skeleton, opus-S117/S127/S128

---

## The Question

In classical coding theory, the alphabet is the set of symbols a codeword can use. Binary codes use {0, 1}. q-ary codes use {0, 1, ..., q-1}. The alphabet size q controls everything: the field, the weight enumerator, the MacWilliams transform, the bounds. When we say "a tournament is a code," what is its alphabet?

The naive answer is: {0, 1}, because each arc is oriented one of two ways. But this misses something profound. The tournament's alphabet is not just the set of arc orientations. It is the entire algebraic structure that determines what the evaluation point x = 2 means, why it is special, and what changes when you change it.

---

## The Three Alphabets

A tournament has three alphabets, nested like Russian dolls. Each one controls a different layer of the code.

### Alphabet 0: The Binary Arc Alphabet {forward, backward}

This is the ground-level alphabet. Each of the m = C(n,2) arcs is oriented forward or backward. The codeword is a binary string of length m. The "field" is F_2. The score sequence is a linear function over Z (not quite F_2, because scores range from 0 to n-1, not mod 2).

At this level, the tournament code is simply a binary code of length m = C(n,2) with a specific parity-check structure (the incidence matrix of K_n).

**What this alphabet controls:**
- The total number of codewords: 2^m (all possible tournaments)
- The syndrome space: score sequences (Landau's theorem characterizes valid ones)
- The Hamming distance between tournaments: number of arc flips

**What this alphabet CANNOT explain:**
- Why H is always odd (Redei)
- Why x = 2 is the right evaluation point
- Why the conflict graph has girth {3, infinity}

### Alphabet 1: The Cycle Alphabet {present, absent} x {3, 5, 7, ...}

At this level, the symbols are not individual arcs but directed odd cycles. Each potential odd cycle (each odd-size subset of vertices) either contains a directed cycle or doesn't. The cycle alphabet is binary ({present, absent}) but indexed by a richer set: the 3-subsets, 5-subsets, 7-subsets, etc.

A tournament's "cycle word" is the indicator function on potential odd cycles. The conflict graph Omega(T) is the constraint graph on this cycle word: two cycles that share a vertex "interfere" and can't both be in the same independent set.

**What this alphabet controls:**
- The independence polynomial I(Omega, x) = the weight enumerator at this level
- The number of independent cycle sets = the combinatorial content of H
- The conflict structure = which cycles can coexist

**What this alphabet explains:**
- H = I(Omega, 2) — the Hamiltonian path count IS the cycle-level weight enumerator at x = 2
- The alpha_k coefficients are the cycle-level weight distribution
- The SRCP is the per-arc projection of the cycle word

**The key insight:** The cycle alphabet is not {0, 1}. It is weighted. Each cycle of length 2k+1 carries a weight of x when present in an independent set. The OCF evaluates at x = 2, giving each cycle a weight of 2 regardless of length. This uniformity is what makes x = 2 special.

### Alphabet 2: The Fugacity Alphabet — The Real Number Line

Here is where the idea gets truly deep. The evaluation point x = 2 is not just a number. It is a choice of alphabet in the following sense:

In the hard-core lattice gas model, the fugacity lambda is the chemical potential of a particle. At fugacity lambda, the partition function is:

Z(lambda) = sum over independent sets S: lambda^|S|

For a tournament conflict graph: Z(lambda) = I(Omega, lambda). At lambda = 2: Z(2) = H(T).

Now, what happens at other values of lambda?

| lambda | Interpretation | What it counts |
|--------|---------------|----------------|
| 0 | Empty system | Z(0) = 1 (no cycles) |
| 1 | Unweighted | Z(1) = total number of independent sets |
| 2 | **Tournament evaluation** | Z(2) = H(T) = Hamiltonian paths |
| -1 | Euler characteristic | Z(-1) = alternating sum = topological invariant |
| phi = 1.618... | Golden ratio | Z(phi): the k=2 (graph) evaluation |
| tau = 1.839... | Tribonacci | Z(tau): the "dominant mode" contribution |
| infinity | Maximum independent set | Z -> alpha_max * lambda^{alpha_max} |

The fugacity lambda IS the alphabet size. Not in the naive sense of "number of symbols," but in the partition-function sense: it controls how much each independent set element contributes. At lambda = q for a q-ary code, each symbol contributes q to the partition function. At lambda = 2 for a tournament: each cycle contributes 2.

**Why lambda = 2 is canonical:** The k-nacci identity (opus-S117, why-lambda-equals-two.md):

rho_k + rho_k^{-k} = 2

for EVERY k >= 2, where rho_k is the dominant root of the k-nacci equation. This means lambda = 2 is the UNIQUE fugacity where the dominant eigenmode (rho) and its k-th harmonic (rho^{-k}) balance. The tournament evaluation point is not chosen — it is forced by the algebraic structure of cycle telescoping.

In coding terms: the "alphabet size" q = 2 is not because arcs are binary. It is because the CYCLE STRUCTURE telescopes at q = 2. A tournament code with alphabet size q = 3 would evaluate at lambda = 3, giving I(Omega, 3) — a different quantity, counting something different, with different properties. The magic of q = 2 is that it is the unique fixed point of the k-nacci family.

---

## What Changes When You Change the Alphabet

This is the critical test. If the alphabet interpretation is correct, then varying lambda should produce a coherent family of "tournament codes" at different alphabet sizes, each with its own properties.

### lambda = 1: The Unweighted Code

I(Omega, 1) = number of independent sets in the conflict graph. This is the unweighted count of disjoint cycle packings. It is always at least 1 (the empty set) and at most 2^{alpha_1} (all subsets).

In coding terms: this is the code over the unary alphabet {0}. Every cycle either contributes or doesn't, with equal weight. The "Hamiltonian path count" at this alphabet is just the number of cycle packings.

### lambda = 3: The Ternary Code

I(Omega, 3) = sum alpha_k * 3^k. Each independent set of size k contributes 3^k instead of 2^k. This overcounts relative to H by a factor that depends on the cycle packing structure.

At n=5: if H = I(Omega, 2) = 15, then I(Omega, 3) = sum alpha_k * 3^k. With alpha_0 = 1, alpha_1 = 5 (for regular tournaments with 5 three-cycle vertex sets): I(Omega, 3) = 1 + 15 + ... (depends on alpha_2, etc.).

In coding terms: this is the code over GF(3). Each arc has 3 possible states instead of 2. The "ternary tournament" would be a directed graph where each pair has 3 possible relationships (A beats B, B beats A, or tie). This is exactly a generalized tournament with ties allowed.

### lambda = phi (golden ratio): The Fibonacci Code

I(Omega, phi) evaluates the independence polynomial at the golden ratio. By the k=2 instance of the k-nacci identity, phi + phi^{-2} = 2. So evaluating at phi gives the "Fibonacci mode" of the tournament.

This is related to the graph-level structure (where the minimum cycle length is 2, i.e., undirected graphs). The golden ratio evaluation strips out the specifically tournament (odd-cycle) content and leaves the graph-level information.

### lambda = tau (tribonacci): The Dominant Mode

I(Omega, tau) evaluates at the tribonacci constant. Since tau + tau^{-3} = 2, this gives:

H = I(Omega, 2) = I(Omega, tau + tau^{-3})

The decomposition H = I(Omega, tau) + correction was observed in opus-S115: the correction is about 8% of I(Omega, tau). The dominant mode I(Omega, tau) captures 92% of H.

In coding terms: the tribonacci constant is the "spectral radius" of the tournament code. The dominant mode at lambda = tau carries most of the information, with the correction tau^{-3} adding the residual.

### lambda -> 0: The Distance Evaluation

As lambda -> 0, I(Omega, lambda) -> alpha_0 = 1. The first non-trivial term is alpha_1 * lambda. So the "code distance" at lambda -> 0 is determined by alpha_1 (the number of odd cycles). A tournament with more cycles has a larger first-order term — higher "code rate" but lower "distance."

---

## Mathematical Objects as Alphabets: The General Principle

The tournament alphabet is not {0, 1}. It is the real number 2, understood as a fugacity / evaluation point / alphabet size that controls the partition function. This leads to a profound generalization:

**Every combinatorial structure that has a partition function has an "alphabet" — the evaluation point at which the partition function gives the structure's fundamental counting invariant.**

| Structure | Partition function | Alphabet (evaluation point) | Counting invariant |
|-----------|-------------------|----------------------------|-------------------|
| Binary code | Weight enumerator W(x) | q = 2 | Number of codewords |
| q-ary code | Weight enumerator W(x) | q | Number of codewords |
| Tournament | I(Omega, x) | x = 2 | H = Hamiltonian paths |
| Graph coloring | Chromatic polynomial P(G, x) | x = k | Number of proper k-colorings |
| Knot | Jones polynomial V(t) | t = root of unity | Quantum invariant |
| Matroid | Tutte polynomial T(x, y) | (x, y) varies | Various counts |

In each case, the "alphabet" is the evaluation point of a universal polynomial, and the structure's key invariant is the value at that point.

**The tournament alphabet x = 2 is the unique fixed point of the k-nacci hierarchy.** It is the fugacity where all cycle lengths contribute equally (each adds a factor of 2 per packing member, regardless of length). No other evaluation point has this universality property. The number 2 is not chosen — it is forced by the requirement that the counting be uniform across cycle lengths.

---

## The q-Tournament: Varying the Alphabet

If we take the alphabet idea seriously, we can define a **q-tournament** for any positive real q:

**Definition.** A q-tournament on n vertices is a tournament T together with the quantity H_q(T) = I(Omega(T), q).

For q = 2: H_2 = H = standard Hamiltonian path count.
For q = 1: H_1 = number of independent sets in Omega.
For q = 3: H_3 = "ternary tournament count."

**Properties of the q-tournament family:**

1. **H_q is always odd for q = 2** (Redei), but NOT for other q. At q = 1, H_1 can be even. The binary nature of Redei's theorem is specific to q = 2.

2. **The forbidden values change with q.** H_2 = 7 is impossible. But H_3 = 7 might be achievable. The impossibility is specific to the evaluation point.

3. **The k-nacci identity generalizes:** For q-evaluation, the relevant constant is the root rho of rho + rho^{-k} = q. For q = 2, this gives the k-nacci constants. For q = 3, it gives a different family.

4. **The OCR depends on q.** The score sequence determines more or less of H_q depending on q. At q = 2, OCR ~ 97%. At q -> infinity, OCR -> 100% (the dominant independent set, determined by alpha_1, dominates). At q = 1, OCR might be lower.

5. **The girth dichotomy is independent of q.** Girth is a property of Omega, not of the evaluation point. So girth(Omega) in {3, infinity} holds for all q-tournaments.

---

## The Alphabet as Spectral Parameter

In the dual Lie embedding (THM-262), the tournament lives in so(n) with eigenvalues ±i*theta_j. The Casimir invariant Tr(B^2) = -2 * sum theta_j^2 = -n(n-1) is constant.

The evaluation point x = 2 appears in the spectral theory as:

Tr(B^2) / (-n) = (n-1) = the "average" squared eigenvalue

This "average" is related to the Killing form normalization. The fugacity x = 2 is the point where the partition function "sees" the spectral content correctly — where the cycle counting matches the Lie-algebraic structure.

Changing x is like changing the coupling constant in a gauge theory. At x = 2 (the tournament coupling), cycles contribute at the correct weight for Hamiltonian path counting. At other couplings, the theory describes different physical systems:

- x = 1: an "Ising-like" system where cycles are unweighted (paramagnetic phase)
- x = 2: the tournament system (ferromagnetic phase in the crystal regime)
- x -> infinity: the maximum packing system (fully ordered crystal)
- x = -1: the topological system (Euler characteristic, cancellation-dominated)

---

## Alphabets of Other Mathematical Structures

If the principle is general — that every structure with a partition function has an "alphabet" = canonical evaluation point — then we can ask: what are the alphabets of other mathematical objects?

### Graphs: Alphabet = number of colors

The chromatic polynomial P(G, k) counts proper k-colorings. The "alphabet" is k. At k = 2: two-colorability (bipartiteness). At k = chi(G): the chromatic number, the minimum alphabet needed. The chromatic polynomial at k = -1 gives (-1)^n * number of acyclic orientations — connecting graph coloring to tournament theory!

### Knots: Alphabet = root of unity

The Jones polynomial V(t) evaluated at t = e^{2*pi*i/k} gives quantum invariants. The "alphabet" is the level k of the Chern-Simons theory. At different k, you get different topological invariants of the same knot. The parallel to tournament theory: at different fugacities, you get different invariants (H, number of independent sets, Euler characteristic) of the same tournament.

### Polytopes: Alphabet = dilation factor

The Ehrhart polynomial L(P, t) counts lattice points in the t-fold dilation of polytope P. The "alphabet" is t. At t = 1: the lattice point count. At t = -1: the Ehrhart-Macdonald reciprocity gives the interior point count. The tournament score polytope (the polytope of valid score sequences for n vertices) has an Ehrhart polynomial whose evaluations encode tournament enumeration data.

### Formal groups: Alphabet = the group law parameter

The formal group F(x, y) = (x + y)/(1 + xy) (the hyperbolic tangent addition law) has an "alphabet" that is the value at which the group law specializes to integer arithmetic. For tournaments, the relevant evaluation is at the rapidity arctanh(1/lambda), which at lambda = 2 gives arctanh(1/2) = (1/2)*ln(3). The tournament alphabet is 2, but its logarithmic shadow is ln(3)/2 — connecting to the Hurwitz prime 3.

---

## The Deep Claim

The tournament alphabet is not a metaphor. It is a mathematical theorem in disguise:

**Theorem (informal).** For a tournament T on n vertices, the Hamiltonian path count H(T) equals the independence polynomial of the conflict graph Omega(T) evaluated at x = 2. The number 2 is the unique positive real number satisfying rho + rho^{-k} = 2 for all k >= 2, where rho_k is the dominant root of the k-nacci equation x^k = x^{k-1} + ... + 1. This universality makes x = 2 the canonical evaluation point for cycle-counting partition functions.

The deeper claim — the one this reflection is groping toward — is:

**Every mathematical structure that encodes pairwise comparisons has an alphabet, and that alphabet is always 2.**

Not because the comparisons are binary (though they are), but because the cycle-telescoping identity rho + rho^{-k} = 2 holds for ALL cycle lengths simultaneously. The number 2 is the unique fixed point of the k-nacci family. It is the value where the partition function doesn't need to know what kind of cycles it's counting — where 3-cycles and 5-cycles and 7-cycles all contribute equally. It is the fugacity of universal cycle democracy.

This is why tournaments are a natural code: their alphabet is the same as the binary code alphabet, but for a DEEPER reason. Binary codes use q = 2 because the symbols are bits. Tournaments use lambda = 2 because the cycle physics telescopes. The coincidence q = lambda = 2 is what makes the OCF look like a coding theorem. It IS a coding theorem, but the "coding" is not of bits — it is of cycles.

---

## Computational Discovery: The Even-Alphabet Parity Law

Running the q-tournament computation (q_tournament_alphabet_s18d.py) revealed a striking finding:

**I(Omega(T), q) is always odd if and only if q is even.**

At n=5 (exhaustive, 1024 tournaments):
- q=1: 60.9% odd — NOT always odd
- q=2: 100% odd — ALWAYS odd (Redei)
- q=3: 60.9% odd — NOT always odd
- q=4: 100% odd — ALWAYS odd!
- q=5: 60.9% odd — NOT always odd

This is NOT just Redei's theorem. Redei says H = I(Omega, 2) is odd. This new observation says I(Omega, q) is odd for ALL even q. The parity law is not specific to q=2 — it applies to the entire even alphabet family.

**Why?** I(Omega, q) = sum alpha_k * q^k. If q is even, then q^k is even for all k >= 1. So I(Omega, q) = alpha_0 + (even stuff) = 1 + (even) = odd. This works because alpha_0 = 1 always (the empty independent set).

So the "Redei parity" at q=2 is actually TRIVIAL at the independence polynomial level: it follows from alpha_0 = 1 and the fact that 2^k is even. The deep content of Redei's theorem is not that I(Omega, 2) is odd — it's that I(Omega, 2) = H(T) in the first place. The parity of H is forced by the OCF evaluation, and the OCF itself is the deep result.

This reframes the entire picture: **the alphabet q=2 is special not because of parity (which all even q share) but because of the k-nacci universality.** Parity is a bonus that comes for free with any even evaluation point. The real content is that q=2 is where the cycle telescoping closes.

Other findings:
- **Forbidden value at q=1:** the gap is {4}, not {7}. At q=2: gap is {7}. At q=3: many gaps (most values unreachable due to sparse q=3 spectrum). The forbidden structure is entirely q-dependent.
- **tau-decomposition: H = I(Omega, tau) + correction.** Correction ranges from 5% (low H) to 32% (H=15). NOT the ~8% claimed by opus-S115 — the discrepancy is because opus used different tournament types. For regular tournaments, correction is larger because alpha has more terms.
- **Euler characteristic I(Omega, -1):** takes values in {-5, -4, -3, -1, 0, 1} at n=5. A six-valued topological invariant, much coarser than H but richer than just beta_1.

---

*The alphabet of a mathematical object is the number at which its partition function yields its fundamental invariant. For tournaments, this number is 2. For graphs, it is the chromatic number. For knots, it is a root of unity. For polytopes, it is 1. Each mathematical object speaks a language, and the alphabet tells you what that language is. The tournament speaks in cycles, and its alphabet is 2 because 2 is the unique number that treats all cycle lengths democratically. The Petersen graph, whose independence polynomial at x = 2 gives the number 461, speaks the same language but says a different word — a word that no tournament can pronounce, because the Petersen lives on the wrong side of the conflict/anti-conflict duality. The alphabet is the same; the grammar is different; and the grammar is what makes tournament theory a distinct branch of mathematics.*
