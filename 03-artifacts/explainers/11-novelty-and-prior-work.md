# What's New: Novelty in Context of Prior Work

## The Landscape Before This Project

To understand what this project contributes, it helps to know what was already known.

Rédei proved in 1934 that every tournament has an odd number of Hamiltonian paths. That's 90 years of knowing the count is odd. The next milestone came in 2023, when Grinberg and Stanley proved a formula involving symmetric functions — a deep algebraic tool — that encodes H(T). Irving and Omar sharpened this in December 2024.

But between "it's odd" and those 2023 formulas, a lot of specific questions went unstudied:
- Which odd numbers can actually occur as H(T)?
- Do Paley tournaments (the rotationally symmetric ones built from prime arithmetic) actually maximize H(T)?
- What does the topology of a tournament look like?

This project attacks all of these. Here's an honest accounting of what is new.

---

## The OCF Formula: New Framing of an Existing Specialization

**What Grinberg-Stanley (2023) proved**: There is a polynomial in infinitely many variables — called the Rédei-Berge symmetric function U_T — that encodes all Hamiltonian paths in T. When you plug in specific values (set x₁ = 1, all odd-index variables = 2, all even-index variables = 0), you recover H(T).

**What this project adds**: The formula H(T) = I(Ω(T), 2) — where Ω(T) is the conflict graph of odd cycles and I is the independence polynomial — names the *intermediate object* that Grinberg-Stanley's specialization produces. The Grinberg-Stanley paper does not introduce or name the conflict graph Ω(T), and does not frame the result as "evaluate the independence polynomial at 2." They arrive at the same number via a different route, without identifying the combinatorial structure in between.

**The analogy**: Imagine someone proves that a certain algebraic sum equals the area of a shape, without drawing or naming the shape. This project draws the shape and says: that's a specific graph (Ω(T)), and evaluating its independence polynomial at 2 is what you're computing.

That reframing matters, because it:
1. Lets you think geometrically about which tournaments have high H(T) (those with rich independent sets in Ω(T))
2. Connects to the physical hard-core gas model (each independent set = a configuration of non-interacting particles)
3. Explains the forbidden values (some numbers can't be independence polynomial evaluations for any graph)

**Bottom line**: The specialization that produces H(T) is in Grinberg-Stanley. The conflict graph Ω(T) and independence polynomial framing are new.

---

## Forbidden Values: Entirely New

No paper before this project asked: which odd numbers can never be H(T) for any tournament?

The only prior congruence result is Grinberg-Stanley's mod-4 refinement (2023): the count H(T) mod 4 is constrained by the number of odd cycles. This tells you something about H(T) modulo 4, but not about which specific integers are permanently impossible.

**The finding that H(T) = 7 and H(T) = 21 are forever forbidden** — for every tournament on every number of players — is entirely novel. The mechanism (claw-freeness of Ω(T) limits what I(Ω,2) can equal) uses a 2007 result by Chudnovsky and Seymour about independence polynomials of claw-free graphs, but its application to forbidden tournament path counts is new.

**The conjecture** that 7 and 21 are the only two permanently forbidden odd numbers is also new to this project.

---

## Paley Tournaments Maximize H(T): Verification of a 1990 Conjecture

In 1990, Noga Alon proved an upper bound on H(T) and conjectured that the Paley tournaments — built from quadratic residues mod a prime — achieve the maximum. This conjecture remained unverified for 35 years.

This project computationally verified:
- T₇ (Paley tournament on 7 players) achieves the exact maximum H(T) among all 7-vertex tournaments
- T₁₁ (11 players) similarly achieves the maximum

No prior paper contains this verification. The numbers had been computed separately (they match OEIS entries), but no one had confirmed they are the global maximum by exhaustive comparison.

For larger primes (p ≥ 13), the cyclic interval tournament appears to overtake Paley — the project also identified this crossover, which is new.

---

## Path Homology: New Computations

The GLMY theory of path homology for directed graphs was developed starting around 2012. But no paper has:
- Computed the Betti numbers for all tournaments on up to n=8 players
- Established that β₂ = 0 for all of these
- Found the seesaw pattern (β₁ and β₃ never coexist for small n)

A February 2026 paper (arXiv:2602.04140) on circulant digraphs is closely related: it uses Fourier methods to compute path homology for circulant graphs (of which Paley tournaments are a special case). That paper finds β₂ = 1 for some circulant digraphs — so β₂ = 0 is NOT universal for all circulant digraphs, only for the specific ones that are tournaments. This makes the tournament β₂ = 0 result more interesting, not less.

**Important caution**: Before publishing the Paley tournament homology results, the overlap with the 2026 circulant paper must be carefully checked, since Paley tournaments T_p are exactly circulant digraphs Cay(Z_p, QR_p).

---

## Krawtchouk Band-Limiting: New

The claim that tournament matrices are "band-limited" in Krawtchouk coordinates — with degree at most 2⌊(n-1)/2⌋ — does not appear in any prior literature. Krawtchouk polynomials and their connection to coding theory are classical, but their application to tournament structure is new.

The observation that Paley tournaments correspond to perfect codes (Hamming [7,4,3] and Golay [23,12,7]) via the Krawtchouk spectral picture is also new, even though the classical connection between quadratic residue codes and Paley constructions is well-known through different mathematical machinery.

---

## What is NOT New

To be complete:
- **Rédei's theorem** (H(T) is always odd): proved 1934
- **The Rédei-Berge symmetric function**: Grinberg-Stanley 2023
- **Paley tournament self-complementarity**: classical
- **Paley tournament eigenvalues** (−1 ± √p)/2: classical spectral graph theory
- **Independence polynomial = hard-core partition function**: classical (Scott-Sokal 2005)
- **Claw-free independence polynomials have real roots**: Chudnovsky-Seymour 2007
- **Quadratic residue codes are perfect**: classical coding theory
- **GLMY path homology theory**: Grigor'yan et al., from 2012

---

## Summary of Novel Contributions

| Contribution | Status |
|---|---|
| Conflict graph Ω(T) + independence polynomial framing of OCF | New framing of Grinberg-Stanley result |
| H(T) = 7 and H(T) = 21 are permanently forbidden | Entirely new |
| Conjecture that 7, 21 are the only permanent gaps | Entirely new |
| Computational verification T₇, T₁₁ maximize H(T) (Alon conjecture) | New verification |
| Paley-to-cyclic crossover at p ≥ 13 | New |
| β₂ = 0 for all tournaments n ≤ 8 | New computations |
| Seesaw pattern β₁ · β₃ = 0 | New observation |
| Tournament matrices are Krawtchouk band-limited | New |
| Paley tournaments ↔ perfect codes via Krawtchouk tiling | New connection |
| CRT + sparse CSC for exact H(T_19), H(T_23) computation | New algorithms |
