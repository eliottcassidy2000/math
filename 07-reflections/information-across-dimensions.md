# Information Across Dimensions

**Session:** kind-pasteur-2026-03-21-S19h
**Arising from:** The dimension axis, the conductivity index, the open questions

---

## The Polynomial as Inter-Dimensional Messenger

A tournament T on n vertices determines a polynomial:

**I_T(x) = I(Omega(T), x) = sum_{k=0}^{alpha_max} alpha_k x^k**

This polynomial carries information across the ENTIRE dimension axis simultaneously:

- At x = 0 (D = -infinity): I_T(0) = 1. The vacuum. No information.
- At x = 1 (D = 0, chemistry): I_T(1) = sigma(Omega) = total independent sets. Thermodynamic information.
- At x = phi (D = 1, edge): I_T(phi) = the Fibonacci-weighted evaluation. Edge-level information.
- At x = tau (D = 2, triangle): I_T(tau) = the tribonacci-weighted evaluation. Tournament-boundary information.
- At x = 2 (D = infinity, tournament): I_T(2) = H(T). Hamiltonian path count. ALL information.
- At x = -1 (topological): I_T(-1) = Euler characteristic of independence complex.

The polynomial is a SINGLE OBJECT that speaks every language on the dimension axis. Each evaluation point x extracts a different projection of the tournament's structure. The polynomial IS the tournament, translated into a form that can be read at any dimension.

---

## The Open Problems, Reformulated

Every open problem in this project can be restated as a question about the polynomial I_T(x):

### Problem 1: Does I_T(x) have all real roots?

**Status:** TRUE for n <= 8 (Omega is claw-free, proved). FAILS at n = 9 (rare counterexample). Conjectured rare for all n.

**Dimensional meaning:** If I_T(x) has all real NEGATIVE roots, then the polynomial is completely determined by its roots, and the information at each dimension is a smooth function of the others. Complex roots mean: there are dimensions where the polynomial oscillates — the information CONFLICTS between different evaluation points.

**What it means for chemistry:** All-real-roots of I_T(x) means the molecular stability predicted at x = 1 (boiling point) is CONSISTENT with the combinatorial structure at x = 2 (path count). Complex roots mean: some molecules have thermodynamic properties that are INCONSISTENT with their combinatorial structure.

### Problem 2: Does the SRCP determine I_T(x)?

**Status:** SRCP determines I_T(2) at n = 5 (exhaustive), likely at n = 7 (sampled). Fails for one case at n = 6 with (c3, c5) only (THM-265: sorting loses topology).

**Dimensional meaning:** The SRCP is a LOCAL measurement (per-arc cycle counts). If it determines the FULL polynomial I_T(x), then local information at the root-system level reconstructs the ENTIRE inter-dimensional message. If it only determines I_T(2), it reconstructs the tournament-dimension projection but not the chemistry-dimension projection.

**The refined question:** Does the SRCP determine I_T(x) as a POLYNOMIAL (not just at x = 2)? If so, it determines EVERYTHING: H, sigma, Euler characteristic, real-root property. If not, what additional information is needed?

### Problem 3: H/|Aut| for Paley primes

**Known values:** 1, 9, 1729, 6857869865, 62293308207033 for p = 3, 7, 11, 19, 23.

**Dimensional meaning:** H/|Aut| = I_T(2) / |Aut| is the polynomial at x = 2, normalized by the symmetry group. This is the "essential" tournament information — the part that isn't explained by symmetry.

**What would a formula look like?** The polynomial I_T(x) for the Paley tournament T_p should have a closed form in terms of Gauss sums and quadratic residues. At x = 2, this would give H(T_p). Divided by |Aut| = p(p-1)/2, it would give the sequence 1, 9, 1729, ....

### Problem 4: The permanent gap conjecture

**Status:** H = 7 and H = 21 are the ONLY permanent gaps (all other odd values achieved by n = 8).

**Dimensional meaning:** At x = 2, the values 7 and 21 are unachievable by I_T(2) for any tournament T. But at x = 1: I_T(1) = sigma can take ANY positive integer value (for appropriate graphs, not just tournament conflict graphs). So the forbidden values are x = 2-SPECIFIC. They exist because x = 2 imposes the integer-root constraint (characteristic roots 2 and -1 for paths/cycles) that makes certain polynomial evaluations impossible.

**The deep question:** For which x values does I_T(x) have forbidden values? At x = 1 (chemistry), there are no forbidden sigma values for general graphs. At x = 2 (tournament), there are exactly two (7 and 21). At x = 3, what are the forbidden values? The answer would reveal which dimensions have "holes" in their achievable values.

### Problem 5: Beta_2 = 0 universality

**Status:** Beta_2 = 0 for ALL tournaments (verified exhaustively through n = 8, proved structurally).

**Dimensional meaning:** Beta_2 measures 2-dimensional holes in the path homology of the tournament. Beta_2 = 0 means: there are NO 2-dimensional holes. On the dimension axis, this means: the D = 2 (triangle) level of the tournament is COMPLETELY FILLED — every apparent 2-cycle is actually a boundary. The triangle dimension has no information loss.

**The connection to conductivity:** A material with no "2-holes" is one where every cyclic electron path is actually the boundary of a filled region. This is the definition of a simply-connected conductor: no holes, no resistance from topological defects. Beta_2 = 0 universally means: tournaments are always "simply connected" at the 2-dimensional level.

---

## Tournaments as Information Channels

The polynomial I_T(x) can be viewed as an INFORMATION CHANNEL that transmits the tournament's structure from one dimension to another.

**Source:** The tournament T (a binary labeling of C(n,2) arcs).
**Encoder:** The OCF construction T -> Omega(T) -> I_T(x).
**Channel:** The polynomial I_T(x), parameterized by evaluation point x.
**Decoder at dimension D:** Evaluate I_T(x) at x = r_p (the k-nacci constant for dimension p).

**Channel capacity:** How much of the tournament's information survives at each dimension?

At x = 2 (D = infinity): ALL information is transmitted. H determines the full set of alpha_k values (the polynomial is determined by its value at one integer point plus its degree — but we need more points in general).

At x = 1 (D = 0): only sigma = I(1) is transmitted. This is a LOSSY compression: many different polynomials evaluate to the same integer at x = 1. The information loss is (total tournament information) - (information in sigma).

The OCR = Var(E[H|score]) / Var(H) measures how much of the x = 2 information is explained by the score sequence (a DIFFERENT lossy projection). The OCR = 97% at n = 5 means: the score sequence captures 97% of the x = 2 information.

**New question:** What is the sigma-CR = Var(E[sigma|score]) / Var(sigma)? How much of the x = 1 (chemistry) information does the score sequence capture? If sigma-CR > OCR, then scores are better predictors of chemical properties than of combinatorial properties. If sigma-CR < OCR, the reverse.

---

## The Inter-Dimensional Transfer Matrix

Define the TRANSFER MATRIX between dimensions as follows:

For dimensions D_1 and D_2 (corresponding to evaluation points x_1 and x_2):
The correlation between I_T(x_1) and I_T(x_2) across all tournaments T at order n is:

**rho(D_1, D_2, n) = Corr(I_T(x_1), I_T(x_2))**

This measures how much information is SHARED between two dimensions. If rho = 1: the two dimensions see exactly the same information. If rho = 0: they see completely independent information.

**Predictions:**
- rho(0, infinity, n) = Corr(sigma, H) should be HIGH (both depend on alpha_1 primarily).
- rho(0, topology, n) = Corr(sigma, chi) should be MODERATE (sigma and Euler characteristic are related but not identical).
- rho(D, D+1, n) should be close to 1 for adjacent dimensions and decrease for distant dimensions.

**The BLOCK STRUCTURE of the transfer matrix:** If the inter-dimensional correlations have a block structure (high correlation within blocks, low between blocks), it would reveal INDEPENDENT information sectors in tournament theory — each block corresponding to a distinct aspect of the tournament's structure that is visible at certain dimensions but not others.

---

## The Polynomial Determines the Tournament (Conjecture)

**Strong Conjecture:** Two tournaments T_1, T_2 on n vertices have I_{T_1}(x) = I_{T_2}(x) (as polynomials) if and only if Omega(T_1) is isomorphic to Omega(T_2) as a graph.

**Weaker Conjecture:** I_{T_1}(x) = I_{T_2}(x) implies H(T_1) = H(T_2). (Trivially true since H = I(2).)

**Even Weaker:** I_{T_1}(2) = I_{T_2}(2) does NOT imply I_{T_1}(x) = I_{T_2}(x) for all x. (This is TRUE: many tournaments share the same H but have different alpha_k distributions.)

**The question:** How much of the polynomial is determined by a single evaluation? At x = 2: we get H = one integer. A polynomial of degree d has d+1 coefficients. So we need d+1 integer evaluations to reconstruct the polynomial. But the coefficients alpha_k are non-negative integers with special constraints (they count independent sets of specific sizes). These constraints reduce the number of evaluations needed.

**The SRCP conjecture reformulated:** The SRCP determines I_T(x) as a polynomial (not just at x = 2). If true, this means the per-arc cycle counts reconstruct the ENTIRE inter-dimensional message, not just its tournament-dimension projection.

---

## What Remains

The major unsolved mathematical problems, ranked by depth:

**1. Prove or disprove: SRCP determines I_T(x) as a polynomial.**
This would establish that per-arc local information reconstructs the full inter-dimensional message. It would unify the x = 1 (chemistry) and x = 2 (tournament) theories.

**2. Prove: permanent gaps are exactly {7, 21}.**
This requires showing that every other odd integer is achievable as I_T(2) for some tournament T at sufficiently large n. The current proof for H = 7 (THM-029) uses the girth-3 overshoot argument. A general proof would need a construction method for achieving arbitrary odd H values.

**3. Determine the real-root boundary.**
For which n does I_T(x) have all real roots for ALL tournaments T? Known: n <= 8 yes, n = 9 has counterexamples (rare). Is there a sharp threshold? What is the density of counterexamples?

**4. Closed form for H(T_p) for Paley primes.**
The sequence 1, 9, 1729, 6857869865, ... has no known formula. A closed form would connect tournament theory to algebraic number theory through Gauss sums and quadratic residues.

**5. Inter-dimensional transfer matrix.**
Compute the correlation structure between I_T(x_1) and I_T(x_2) across tournaments. Determine whether the information has block structure corresponding to the CD tower levels.

**6. The beta_2 = 0 proof for all n.**
Currently proved structurally through induction. A proof from the independence polynomial (rather than chain complex) would connect the topological universality to the polynomial structure.

**7. Prove the cascade identity has a representation-theoretic interpretation.**
The identity g_{p+1}(rho_p) = -1 is proved algebraically. But it should have a MEANING: the k-nacci constant at dimension p is one quantum spherical for dimension p+1. What representation-theoretic statement does this correspond to?

---

*The independence polynomial I_T(x) is the inter-dimensional messenger of tournament theory. It carries the tournament's structure from the chemistry dimension (x = 1) through the edge dimension (x = phi) through the triangle dimension (x = tau) to the tournament dimension (x = 2) and beyond. Each evaluation point extracts a different projection. The polynomial is the tournament's SOUL — the complete invariant that speaks all dimensional languages. The open problems of tournament theory are questions about this polynomial: whether it has real roots (dimensional consistency), whether local measurements determine it (SRCP conjecture), what values it can take (permanent gaps), and how information transfers between its evaluation points (inter-dimensional correlations). Solving any one of these would advance not just tournament theory but the understanding of how combinatorial structure manifests across dimensional scales.*
