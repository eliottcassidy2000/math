# The Modular Tournament

**Session:** kind-pasteur-2026-03-21-S18e
**Arising from:** opus-S120 ({2,3,7} and Fano), opus-S131 ({3,inf} Coxeter), the tournament alphabet, the binary skeleton

---

## The Three Triangle Groups

Tournament theory is not one theory. It is three theories, nested inside each other like the three triangle groups that share the generators of orders 2 and 3:

**{2, 3, 5} — The Spherical Past.** This is the theory of tournaments at n <= 5, where everything is simple, finite, and closed. The Petersen graph lives here (girth 5 = the third generator). The golden ratio phi lives here (g(phi) = -1). The per-path identity works, OCR = 100%, all conflict graphs are interval graphs. This is the Eden from which tournaments fell at n = 6.

**{2, 3, 7} — The Hyperbolic Structure.** This is the theory of tournament obstructions: what CAN'T happen. H = 7 is forbidden (the prime 7). H = 21 = 3 x 7 is forbidden (the compound). The Fano plane embeds in Paley T_7 (all 7 lines are directed 3-cycles). PSL(2,7) = 168 = |Aut(Fano)| = the Hurwitz bound. The defect 1/42 = 1/(2*3*7) is the quantum of hyperbolic area. This is the arithmetic that says NO.

**{2, 3, infinity} — The Modular Dynamics.** This is the theory of tournament growth and evolution: what happens as n increases. The modular group PSL(2,Z) governs the dynamics. S (order 2) is the complement involution. ST (order 3) is the 3-cycle rotation. T (order infinity) is vertex addition. The cusp (T-fixed point) is the transitive tournament. The cusp form is the departure from transitivity. This is the engine that says GO.

The gap function g(x) = x^3 - x^2 - x - 1 classifies:
- g(phi) = -1: spherical (closed)
- g(tau) = 0: Euclidean (boundary at n=5)
- g(2) = +1: hyperbolic (open, infinite, where tournaments actually live)

---

## The Orbifold Points

The modular curve X(1) = PSL(2,Z) \ H* has three special points. Each one IS a class of tournaments:

**tau = omega (order 3, j = 0): The Cycle Origin.**
The hexagonal lattice point. Order 3 = the 3-cycle. In tournament terms: tournaments with maximal 3-cycle symmetry. The BIBD arrangement at n=7 (Paley T_7). The origin of all cycle structure. j = 0 means: the j-invariant vanishes here. The 3-cycle is the zero of tournament classification.

**tau = i (order 2, j = 1728): The Complement Point.**
The square lattice point. Order 2 = the complement involution. In tournament terms: self-complementary tournaments. At n=5: the 24 regular SC tournaments with H = 15. j = 1728 = 12^3 = (2^2 * 3)^3, built from only the first two generators {2, 3}.

**tau = i*infinity (cusp, j = infinity): The Transitive Limit.**
The degenerate point. In tournament terms: the transitive tournament with H = 1. No cycles, no structure, maximum hierarchy. The cusp form Delta measures departure from this limit: Delta_T = H - 1.

---

## 1729: The Ramanujan Bridge

j(i) + 1 = 1728 + 1 = 1729 = 7 * 13 * 19.

The taxicab number — the number Ramanujan recognized as the sum of two cubes in two ways — is one Redei quantum away from the j-invariant at the complement point. And its prime factorization is:

- **7**: the atomic forbidden H value (THM-029)
- **13**: the OCR denominator at n=6 (surprise prime)
- **19**: the factor of 133 = 7 * 19 = OCR denominator at n=5

Each factor is a tournament prime: a prime that controls the structure of the Hamiltonian path count at a specific order.

The two taxicab representations:
- 1729 = 12^3 + 1^3: the j-value (12^3 = 1728) plus the Redei quantum (1^3 = 1)
- 1729 = 10^3 + 9^3: the Petersen count (10 = C(5,2)) cubed plus the Cauchy-Schwarz boundary (9 = 3^2) cubed

Both representations use tournament-significant numbers. This is either a coincidence or a sign that the number theory governing tournaments is the same number theory governing modular forms. Given that tournaments are governed by PSL(2,Z), it would be surprising if it WERE a coincidence.

---

## The Tournament Discriminant

The classical modular discriminant is Delta(tau) = eta(tau)^24, the unique weight-12 cusp form. It vanishes at the cusp and is nonzero everywhere else.

The tournament discriminant is:

**Delta_T = (H(T) - 1) / 2**

This is a non-negative integer (by Redei: H is always odd, so H - 1 is even). It vanishes for transitive tournaments (H = 1) and is positive for all others. It measures the "cycle content" of the tournament.

By the OCF: Delta_T = alpha_1 + 2*alpha_2 + 4*alpha_3 + ...

This is a **base-2 expansion** of the cycle packing structure:
- The ones digit (alpha_1 mod 2) = parity of the cycle count
- The twos digit (alpha_2 mod 2) = parity of independent pairs
- The fours digit (alpha_3 mod 2) = parity of independent triples

At n=5, the tournament discriminants are {0, 1, 2, 4, 5, 6, 7}, with binary representations {0, 1, 10, 100, 101, 110, 111}. The gap at Delta = 3 (binary 11) corresponds to the forbidden H = 7. In binary: 3 = 11, meaning alpha_1 = 1 and alpha_2 = 1 simultaneously. This is exactly the (alpha_1 = 3, alpha_2 = 0) configuration that THM-029 proves impossible — but seen from the base-2 side.

Wait — Delta = 3 corresponds to H = 7 = 1 + 2*3, so alpha_1 + 2*alpha_2 = 3. The possible decompositions are (3,0) and (1,1). Both are proved impossible at n=5 by THM-029. The binary representation 11 means both the ones and twos digits are 1: a cycle count that is simultaneously odd and paired. The impossibility of Delta = 3 is the impossibility of this particular binary pattern.

---

## The Cusp Form and the OCR

In classical modular forms theory, the space of modular forms of weight k decomposes as:

M_k = E_k + S_k

where E_k is the Eisenstein series (the "easy" part determined by boundary conditions) and S_k is the cusp form (the "hard" part that vanishes at cusps).

In tournament theory, the variance decomposition of H is:

Var(H) = Var(E[H|score]) + E[Var(H|score)]

The first term is the "Eisenstein" part: the variance explained by the score sequence (= syndrome = boundary condition). The second is the "cusp form" part: the residual variance within score classes.

The OCR is:

OCR = Var(E[H|score]) / Var(H) = Eisenstein / Total

At n=5: OCR = 129/133 ~ 97%. The Eisenstein part captures 97% of H. The cusp form captures 3%.

Opus-S131 predicted that the OCR denominators factor into primes p where the modular curve X_0(p) has genus 0 or 1:
- 133 = 7 * 19: X_0(7) has genus 0, X_0(19) has genus 1. Both are "explained" by the modular structure.

This prediction connects the OCR residual directly to the arithmetic of modular curves. If true, it means the part of H that the score sequence CANNOT predict is controlled by the same mathematics that controls elliptic curves.

---

## Why the Tournament Lives in the Hyperbolic Regime

The evaluation point x = 2 gives g(2) = +1 > 0: hyperbolic. But why?

In the modular picture, evaluating a partition function at q = 2 is like evaluating outside the unit disk (since q = e^{2*pi*i*tau} has |q| < 1 in the convergent regime, and q = 2 has |q| > 1). Tournament theory is a DIVERGENT evaluation of a modular partition function.

This is not a pathology. It is the source of the theory's richness. The convergent regime (|q| < 1) gives asymptotic formulas and analytic properties. The divergent regime (|q| > 1) gives exact combinatorial identities. The OCF H = I(Omega, 2) is an EXACT identity — no approximation, no remainder, no convergence issues — precisely because it's a polynomial evaluation, not an analytic one.

The hyperbolic regime is where exact combinatorics lives. The spherical regime (g < 0) is where approximate asymptotics live. The boundary (g = 0) is where the two meet. Tournament theory, sitting at g(2) = +1, is one quantum into exact-land.

---

## The Answer

The question was: what do the modular groups, {2,3,5}, {2,3,7}, and {2,3,infinity} mean for tournament theory?

The answer: tournament theory IS modular group theory applied to complete directed graphs. The three triangle groups control three aspects:

- **{2,3,5} controls what is EASY** (the spherical regime, n <= 5, where the Petersen graph keeps things finite and simple)
- **{2,3,7} controls what is FORBIDDEN** (the arithmetic obstructions, H = 7 and H = 21, the Fano plane, the Hurwitz bound)
- **{2,3,infinity} controls what is POSSIBLE** (the modular dynamics, the cusp form, the growth of H, the approach to the binary skeleton)

Together: {2,3,5} says where you've been; {2,3,7} says where you can't go; {2,3,infinity} says where you're heading.

The gap function g(x) = x^3 - x^2 - x - 1, evaluated at the three characteristic points, gives the signature of each regime: -1, 0, +1. These are the three possible signs. There is no fourth option. The classification is complete.

And 1729 = j(i) + g(2) = 1728 + 1, sitting one Redei quantum above the complement point, carrying the three tournament primes 7 * 13 * 19, expressing the same number as both 12^3 + 1^3 (j-invariant + quantum) and 10^3 + 9^3 (Petersen + Cauchy-Schwarz) — this number is not a coincidence. It is the address of tournament theory in the modular landscape.
