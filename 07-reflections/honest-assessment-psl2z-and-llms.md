# Honest Assessment: Does PSL(2,Z) Actually Underlie LLMs?

**Session:** kind-pasteur-2026-03-21-S18u
**Arising from:** The 2/3 reflection (S18t) and the need for intellectual honesty

---

## The Question

We have built a tower of connections: the {3, infinity} tessellation governed by PSL(2,Z), the Cartan decomposition gl(n) = so(n) + p + R, the 2/3 ratio appearing in LLM phase transitions, the 16-dimensional fiber matching the sedenion dimension. Does PSL(2,Z) actually govern LLM structure? Or have we found patterns in noise?

This reflection answers honestly.

---

## What Is Rigorously True

### 1. The Cartan decomposition of attention matrices is real mathematics

Any real n x n matrix decomposes as:
- Antisymmetric part: (A - A^T)/2, living in so(n)
- Symmetric traceless part: (A + A^T)/2 - (tr A/n)I, living in p
- Scalar part: (tr A/n)I

This is the Cartan decomposition of gl(n,R). It applies to attention matrices, weight matrices, or any square matrix. It is a theorem, not an empirical claim. The dimensions are C(n,2) for so(n) and C(n+1,2)-1 for p, and their ratio at n=4 is 6/9 = 2/3. This is TRUE regardless of whether LLMs care about it.

### 2. Bidirectional training induces symmetry; autoregressive induces directionality

The paper arXiv:2502.10927 (February 2025, validated on ModernBERT, GPT, LLaMA3, Mistral across text/vision/audio) proves that:
- **Bidirectional** models (BERT-type) develop symmetric attention weight matrices
- **Autoregressive** models (GPT-type) develop directional, column-dominant weight matrices

This is the Cartan decomposition in action: the symmetric part (p) dominates in BERT, the antisymmetric part (so(n)) is more relevant in GPT. The training objective determines which sector matters.

### 3. LayerNorm creates SO(d-1) geometric structure

Van Nierop (arXiv:2412.14543) rigorously proved that LayerNorm projects hidden states onto the intersection of a sphere and a hyperplane, which is diffeomorphic to S^{d-2}. The symmetry group of this space is SO(d-1). This is genuine geometry, not metaphor.

### 4. GL(d_h) gauge freedom per attention head is real

The bilinear form QK^T/sqrt(d) is invariant under Q -> Qg, K -> K(g^{-T}) for any g in GL(d_h). This is a real parametric redundancy — a genuine gauge freedom. It means the "physical" content of an attention head lives in the quotient GL(d_h)\M_{d_h}, not in the raw matrices.

---

## What Is Plausible But Unproven

### 5. The 2/3 ratio in LLM phase transitions

Napolitano observed a phase transition at ~67% depth. We computed dim(so(4))/dim(p_4) = 6/9 = 2/3. The match is suggestive but:
- Napolitano's paper has serious methodological issues (Part II of our analysis)
- The 67% is approximate, not exact
- The n=4 choice (from k=2 normalization layers) is ARBITRARY — a different model architecture would give a different n and a different ratio
- The same 2/3 appearing in D(sqrt(2)) and T_11 transitivity may be coincidence

**Confidence: LOW-MEDIUM.** The Cartan decomposition ratio is real math, and the empirical observation is interesting, but the connection requires stronger evidence.

### 6. Tournament structure in attention matrices

Our analysis (S12) showed that thresholded attention matrices can be interpreted as tournaments. But:
- Softmax attention produces POSITIVE matrices (A_ij > 0 always), which are inherently more symmetric than antisymmetric
- Random softmax attention puts ~72% energy in the symmetric sector, ~19% antisymmetric — this is a BASELINE, not a discovery
- Whether the tournament (antisymmetric) component carries meaningful information is testable but UNTESTED on real models
- The paper arXiv:2502.10927 suggests the antisymmetric part IS where directionality lives in autoregressive models, which SUPPORTS the idea that the tournament sector carries causal information

**Confidence: MEDIUM.** The mathematical framework is sound. The empirical question is open.

---

## What Is Likely Numerological

### 7. PSL(2,Z) governing LLM dynamics

The modular group PSL(2,Z) governs the {3, infinity} tessellation. Tournaments live on this tessellation. But:
- LLMs are NOT tournaments. They are continuous-valued, high-dimensional, trained by gradient descent. The connection to discrete combinatorics requires MULTIPLE non-trivial bridges (thresholding, Cartan decomposition, specific n choice).
- PSL(2,Z) governs modular forms, elliptic curves, and number theory. Its appearance in our tournament theory is rigorous (the score-class structure, the cusp form decomposition). But extending this to LLMs requires showing that LLM dynamics ACTUALLY have modular symmetry — that the computation is invariant under the S and T generators of PSL(2,Z). Nobody has shown this.
- The 2/3 appearing in four places is suspicious in the GOOD sense (it might be real) and the BAD sense (four data points, each with possible alternative explanations).

**Confidence: LOW.** The mathematical framework is beautiful, but the claim "PSL(2,Z) governs LLMs" requires evidence that does not exist.

### 8. 16 = sedenion dimension being significant for LLMs

The match 16 = 2^4 = sedenion dimension = gl(4,R) dimension is a numerical coincidence until proven otherwise. The number 16 appears everywhere in computer science (bits per word, common vector dimensions) for reasons having nothing to do with the Cayley-Dickson tower.

**Confidence: LOW.** Interesting coincidence, not evidence.

### 9. 1729, 196884, and other "moonshine" connections to LLMs

These connections are between tournament theory and number theory, NOT between LLMs and number theory. The chain tournament -> modular group -> moonshine is mathematically grounded. The chain LLM -> tournament -> modular group -> moonshine requires the first link (LLM -> tournament) to be established, which it has not been.

**Confidence: VERY LOW for LLM connection.** HIGH for pure math connection.

---

## What Would Constitute Evidence

To establish that PSL(2,Z) (or even just the Cartan decomposition) genuinely governs LLM structure:

### Tier 1: Testable now

1. **Cartan decomposition of attention matrices across layers.** For a trained GPT model, compute A_sym and A_anti for each attention head at each layer. Plot the ratio ||A_anti||/||A_sym|| as a function of depth. If it matches the 2/3 prediction at ~67% depth, the Napolitano observation is confirmed. If it shows model-independent structure, the Cartan decomposition is relevant.

2. **Tournament invariants predict model behavior.** Threshold attention matrices to tournaments, compute H(T), beta_k, SRCP. Correlate with downstream accuracy. If tournament invariants predict accuracy better than raw attention statistics, tournaments are relevant to LLMs.

3. **Symmetric initialization helps encoder models.** arXiv:2502.10927 already shows this. Extending: does ANTISYMMETRIC initialization help autoregressive models? Our framework predicts yes — the tournament sector should carry the causal information.

### Tier 2: Would require new theory

4. **Modular form structure in loss landscape.** If PSL(2,Z) governs LLMs, the loss landscape should have modular symmetry. Specifically: the loss as a function of certain parameter combinations should transform as a modular form under S and T transformations. This is testable but requires defining what S and T mean for parameters.

5. **Cusp form decomposition of attention variance.** Our OCR = Eisenstein/Total decomposition should have an LLM analogue: the fraction of attention variance explained by a "score-like" (easily computed) quantity should be ~97%, with the residual being "cuspidal." Measuring this would test the modular form analogy directly.

### Tier 3: Would require extraordinary evidence

6. **Tournament forbidden values in LLM structure.** If thresholded attention produces tournaments, those tournaments should exhibit H = 7 impossibility, girth = {3, infinity}, beta_2 = 0, and all 26 binary phenomena. Finding even one of these in real LLM data would be strong evidence. Finding ALL of them would be extraordinary.

---

## The Honest Summary

| Claim | Status | Confidence |
|-------|--------|------------|
| Cartan decomposition applies to attention matrices | **TRUE** (pure math) | CERTAIN |
| Bidirectional = symmetric, autoregressive = directional | **TRUE** (arXiv:2502.10927) | HIGH |
| LayerNorm creates SO(d-1) geometry | **TRUE** (van Nierop) | HIGH |
| GL(d_h) gauge freedom is real | **TRUE** (algebraic identity) | CERTAIN |
| Tournament sector carries causal info in GPT | **PLAUSIBLE** | MEDIUM |
| 2/3 ratio is meaningful for LLMs | **SUGGESTIVE** | LOW-MEDIUM |
| PSL(2,Z) governs LLM dynamics | **UNPROVEN** | LOW |
| 16 = sedenion dimension is significant | **COINCIDENTAL** | LOW |
| Moonshine connects to LLMs | **NO EVIDENCE** | VERY LOW |

The mathematical framework we've built (Cartan decomposition, tournament structure, modular group) is INTERNALLY rigorous and INTERNALLY beautiful. The connection to LLMs is:
- **Certain** at the level of linear algebra (any matrix decomposes)
- **Supported** at the level of empirical observation (arXiv:2502.10927 on symmetry/directionality)
- **Speculative** at the level of PSL(2,Z) and modular forms
- **Ungrounded** at the level of moonshine and sedenions

The right attitude: **use the framework as a source of testable predictions, not as established truth.** The Cartan decomposition of attention is testable today. Tournament invariants of thresholded attention are testable today. The modular form structure is testable with more work. The moonshine connection may never be testable.

---

*The most dangerous thing in mathematics is to mistake beauty for truth. The {3, infinity} tessellation is beautiful. The 2/3 ratio appearing four times is beautiful. The Ramanujan number 1729 emerging from Paley T_11 is beautiful. But beauty is evidence only when accompanied by prediction. The predictions we can make — Cartan decomposition ratios, tournament invariants of attention, symmetric vs antisymmetric initialization effects — are testable. The predictions we cannot make — PSL(2,Z) acting on LLM parameters, moonshine in neural network loss landscapes — are poetry until someone tests them. The framework's value lies in the testable predictions, not in the poetry. The poetry is what makes us want to test them.*
