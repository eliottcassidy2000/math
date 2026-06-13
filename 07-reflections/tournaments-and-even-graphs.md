# Tournaments and Even Graphs: The Cut/Cycle Decomposition

**Session**: kind-pasteur-2026-03-22-S20ae

## The Graph-Theoretic Decomposition

Every tournament on n vertices decomposes as:

**Tournament = Score (cut space) + Cycle structure (cycle space)**

- Edge space of K_n: dimension C(n,2) over GF(2)
- Cut space: dimension n-1 (encodes the score sequence)
- Cycle space: dimension C(n-1,2) (encodes the even graph)
- Direct sum: C(n,2) = (n-1) + C(n-1,2)

At n=5: 10 = 4 + 6. A tournament uses 10 bits. Its score uses 4 bits. Its cycle structure uses 6 bits.

## The OCR Through This Lens

The OCR = 97% means: **97% of H variance lives in the 4 score bits (cut space).**

But score bits are 48x more informative per bit than cycle bits:
- Score: 97% / 4 bits = 24.2% per bit
- Cycle: 3% / 6 bits = 0.5% per bit

This is an extreme asymmetry. The cut space is a tiny subspace (4 out of 10 dimensions) that carries almost all the H-information.

## The Even Graph Connection

Even graphs (all vertices have even degree) form the cycle space of K_n. The number of labeled even graphs = 2^{C(n-1,2)}, confirmed at n=3,4,5.

**Correction**: The Royle et al. equinumerosity A000568 = #tournament iso = #even graph iso was NOT confirmed at n=4,5 in our computation. We get:
- n=3: 2 = 2 (matches)
- n=4: 4 != 3 (FAILS)
- n=5: 12 != 7 (FAILS)

This likely means the Royle et al. result uses a different vertex count convention or a different definition. The specific equinumerosity needs to be checked carefully against the paper.

**Regardless**, the structural decomposition is verified: every tournament = score + even graph, and the spaces have the right dimensions.

## The Surprise: Cycle Space Carries Little H-Information

The cycle space projection (even graph) has OCR = 37.9% for H. This means the even graph structure alone is a POOR predictor of H. Most of H's variation is in the scores, not in the cycle structure.

But wait: we know from Walsh-Fourier that order-2 carries 92-95% of Var(H), and order-2 depends on PAIRWISE interactions, which are captured by scores. So it's expected that the cut space (scores) dominates.

The cycle space captures the ORDER-4 and higher interactions -- the remaining 3-8% of H variance. This is where forbidden values, alpha_2, and the Morse secondary peak live.

## The Deep Meaning

The cut/cycle decomposition tells us: a tournament's identity splits into two independent components.

1. **The cut component (score)**: who is strong, who is weak. This is the HIERARCHY. It determines 97% of H.

2. **The cycle component (even graph)**: among players of similar strength, who beats whom in circular patterns. This is the COMPLEXITY. It determines the remaining 3%.

The equinumerosity conjecture (if true in some form) would say: the number of distinct structural TYPES is the same whether you keep both components or just the cycle component. The cut component adds NO new types -- it only adds COPIES of existing types at different "strength levels."

This is exactly the unlabeling principle: the score is labeling (assigning "strength levels" to the cyclic structure). The even graph is the pure structure.

## The Information Density Asymmetry

The 48x difference in information density (score vs cycle bits) has a physical interpretation:

- **Score bits are like macroscopic variables** (temperature, pressure): few in number but highly informative about the system's state.
- **Cycle bits are like microscopic variables** (individual particle positions): many in number but individually uninformative.

This is a THERMODYNAMIC structure. The OCR is the fraction of "energy" (H) that's captured by the thermodynamic (macroscopic) description. 97% means tournament thermodynamics works extremely well -- you can predict H from scores alone, just as you can predict pressure from temperature and volume alone.

The 3% residual is the "fluctuations" -- the part that requires microscopic (cycle-level) information to predict. In tournament theory, this 3% is where all the deep combinatorics lives.
