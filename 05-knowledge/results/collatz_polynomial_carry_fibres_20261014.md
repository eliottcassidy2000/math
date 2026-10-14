# Polynomial carry fibres: the exponential compression hypothesis is false

Date: 2026-10-14. PROVED: an O(L²) bound on every slope-admissible carry fibre, equality of word and carry exponential growth rates, and equality of the corresponding native-mass rates. FINITE-EXACT: projected counts through L=16. OPEN: finer fibre statistics and any sourcewise descent consequence. No external novelty claim.

## 1. The decisive result

The [carry-frontier note](collatz_survivor_carry_frontier_20261014.md) found genuine collisions and proposed that they might lower the exponential growth rate of distinct ternary carries. That hypothesis is now refuted for its stated slope-admissible word language.

Let W_L count positive valuation words of length L satisfying 2^Ai<3ⁱ at every prefix. Let R_L count their distinct terminal carries B·2^(−A) mod 3ᴸ. Put α=log₂3, bL=floor(αL), and

    C_L=(bL−L+1)ceil(L/3).

Then every carry fibre contains at most C_L words. In particular,

    W_L/C_L≤R_L≤W_L.

The earlier composition argument gives the word growth rate

    β=α^α/(α−1)^(α−1)≈2.8395137305.

Since C_L is polynomial,

    lim R_L^(1/L)=lim W_L^(1/L)=β.

Finite collisions and the stable substitution remain real. They cannot reduce the exponential rate of this language. The ternary carry loses at most O(log L) bits of word identity, despite forgetting the exact affine carry.

## 2. Why a carry fibre is polynomially small

For a length-L word, with A₀=0 and Aᵢ the partial valuation sums,

    B=Σ_{i=0}^{L−1}3^(L−1−i)2^Ai.

Slope admissibility gives 2^Ai≤3ⁱ, with equality only at i=0. Therefore

    0<B≤L·3^(L−1).

Fix the total cost A and the terminal carry r modulo 3ᴸ. The exact B must satisfy

    B≡2ᴬr mod 3ᴸ.

In the displayed interval there are at most ceil(L/3) such B values. Also A ranges over only bL−L+1 possibilities.

It remains to show that (A,B) determines the word, rather than merely a larger class. This has an elementary exact decoder. For L≥2,

    B=3^(L−1)+2^a₁ B_tail,

where B_tail is odd. Hence

    a₁=v₂(B−3^(L−1)),
    B_tail=(B−3^(L−1))/2^a₁.

Repeat on the tail. After recovering the first L−1 letters, the total A determines the final letter. Thus (A,B) is injective on all positive valuation words of a fixed length, independently of the slope restriction. Multiplying the two bounds proves C_L.

The proof does not assert injectivity of (A,r). That stronger property is observed through length 16, but is unnecessary and remains unproved here. Parity of B can improve the numerical fibre bound; the simpler bound above already resolves the exponent question.

## 3. Source weights also preserve the exponential rate

Let

    S_L=Σ_w 2^(−A(w))

be the exact relative odd-source mass of the disjoint native cylinders of slope-admissible words. Let E_L(r) be the minimum cost in a carry fibre, and define the mass of its cheapest representatives by

    S*_L=Σ_r 2^(−E_L(r)).

Each fibre has at least its largest word weight and at most C_L times that weight. Therefore

    S*_L≤S_L≤C_L S*_L.

This does not permit replacing the true mass by S*_L in a finite statement: the preceding note's length-six loss of 9/512 remains a concrete counterexample. It does imply equality of their exponential rates.

For completeness, that common rate can be calculated without fitting data. At fixed total cost bL there are binom(bL−1,L−1) positive compositions. At least one cyclic rotation in every orbit is prefix-admissible, so

    S_L≥binom(bL−1,L−1)/(L·2^bL).

For all sufficiently large L, the unconstrained composition terms

    binom(A−1,L−1)2^(−A),  L≤A≤bL,

increase with A, because their adjacent ratio is A/[2(A−L+1)]>1 before A=2L−2, and bL<2L−2. Thus

    S_L≤(bL−L+1)binom(bL−1,L−1)2^(−bL).

Stirling's formula and 2^α=3 give

    lim S_L^(1/L)=lim (S*_L)^(1/L)=β/3≈0.9465045768.

This is a rate for ensembles of native slope-admissible prefixes. It is not a rate for the stopping time of a particular integer and does not classify finite-source-only exceptions to the slope condition.

## 4. Bounded counts: a finite loss, not a fitted exponent

The independent minimum-cost and weighted-polynomial recursions were continued only to length 16. The integer budget bL is retained explicitly, so its one-or-two-step jumps are not smoothed away.

| L | bL | Words W_L | Carries R_L | R_L/W_L, approximately | Largest fibre |
|---:|---:|---:|---:|---:|---:|
| 6 | 9 | 30 | 24 | 0.800000 | 2 |
| 7 | 11 | 85 | 67 | 0.788235 | 2 |
| 8 | 12 | 173 | 137 | 0.791908 | 2 |
| 10 | 15 | 961 | 768 | 0.799168 | 2 |
| 12 | 19 | 8,045 | 6,464 | 0.803480 | 3 |
| 14 | 22 | 51,033 | 41,192 | 0.807164 | 4 |
| 16 | 25 | 312,455 | 253,264 | 0.810562 | 4 |

Every (carry,total cost) pair in this finite universe has multiplicity one. Thus the weighted polynomials have coefficients zero or one through length 16, although a carry can contain several distinct costs. This is a finite observation, not the mechanism used to prove equal exponential rates.

The theorem resolves the asymptotic question independently of these modest counts. In particular, a ratio near 0.8 is not evidence for an exponentially shrinking ratio.

## 5. A sharper surviving question

Normalize the native word weights 2^(−A) within the surviving language, and regard the word and its terminal carry as random variables. The conditional entropy H(word | carry) measures the ambiguity remaining after observing the carry, with the source weights retained.

For any probability distribution supported on this language, not only the native weighting,

    H(word | carry)≤log₂ C_L=O(log L).

Thus no extensive, linear-in-L entropy loss can occur under this quotient. The finite native-weighted values are approximately 0.387406 bits at L=6, 0.502233 at L=12 and 0.500751 at L=16. These numbers are descriptive outputs, not a claimed limiting constant.

The useful remaining question is finer: are the native-weighted conditional entropy and typical fibre size bounded as L grows, or do they grow slowly? Fixed-cost injectivity, if eventually false, would help identify the first additional lifting mechanism. None of these alternatives would by itself imply descent of one positive source.

The source-preserving boundary from the preceding note still applies. The all-ones Mersenne branch survives every finite depth while its positive witnesses move to infinity. A small carry fibre, or even a fully decoded word, does not supply a fixed positive starting integer for an infinite branch.

## 6. Reproduction and supersession

The standard-library checker computes all projected frontiers through length 16, compares the minimum-cost compiler with the full cost-polynomial frontier, and compares word counts with a separate cost-only recursion. It checks exact B witnesses and the inverse decoder on bounded samples, records every fibre-size histogram, and verifies the polynomial comparison for source mass. Approximate entropy is labelled separately from exact counts and fractions.

11,013 checks pass. Normal and optimized outputs are byte-identical.

    python3 04-computation/experiments/collatz_polynomial_carry_fibres_20261014.py > 05-knowledge/results/collatz_polynomial_carry_fibres_20261014.json
    python3 -O 04-computation/experiments/collatz_polynomial_carry_fibres_20261014.py > 05-knowledge/results/collatz_polynomial_carry_fibres_20261014.optimized.json
    cmp 05-knowledge/results/collatz_polynomial_carry_fibres_20261014.json 05-knowledge/results/collatz_polynomial_carry_fibres_20261014.optimized.json

Canonical output-payload SHA-256, excluding its hash field:

    759cd7bd8afee562dbb723e0b73db9721e846dd2d52e20973847353ebec6ebc7

Superseded: only the strict carry-entropy-gap hypothesis in the preceding note. Its finite counts, source inequalities, stable collision identity, minimum-cost recursion, weighted recursion and sourcewise limitations remain intact. This is a resolved research hypothesis, not a retraction of those proved statements.
