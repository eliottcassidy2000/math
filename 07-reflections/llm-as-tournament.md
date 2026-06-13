# An LLM Is a Tournament Machine

**Session:** kind-pasteur-2026-03-17-S116n33

## The Observation

At every generation step, an LLM ranks ~50,000 tokens by probability of being "next." The logits give a complete ordering. Softmax converts to probabilities. Sampling picks one.

**This IS a tournament.** The items are tokens. The comparisons are implicit in every logit difference. The "result" of each comparison is: token A > token B iff logit(A) > logit(B).

## The Logit Difference IS a Rapidity

The log-odds of token A beating token B:

**logit(A) - logit(B) = ln(p(A)/p(B)) = 2 · arctanh((p(A)-p(B))/(p(A)+p(B)))**

This IS the rapidity in the Cayley formal group. Every logit difference is an arctanh — the formal group logarithm. The entire logit vector is a point in the rapidity lattice of the token tournament.

The LLM doesn't know it, but it's computing formal group logarithms at every step.

## Temperature IS the Polynomial Variable z

The softmax at temperature T:

**p(token_i; T) = exp(logit_i / T) / Z(T)**

At T=1: standard generation. At T→0: greedy (argmax). At T→∞: uniform (random).

The expected quality of the next token at temperature T is:

**E[quality; T] = Σ quality(token_i) · exp(logit_i / T) / Z(T)**

This IS our polynomial P(z) evaluated at z related to T. Temperature IS the polynomial variable that interpolates between "current best guess" (T→0, z=1) and "average over all tokens" (T→∞, z=0).

## The Five Coefficients of an LLM

The polynomial model of an LLM's output at each step:

| Coefficient | What it knows | LLM analog | Transformer component |
|------------|---------------|------------|----------------------|
| **A_0** | Prior (average token quality) | Unigram model | Embedding layer |
| **A_1** | Evidence (last-token context) | Bigram model | First attention layer |
| **A_2** | Patterns (pairwise context) | Trigram model | Pairwise attention heads |
| **A_3** | Contradictions (conflicting signals) | Ambiguity detector | Multi-head cross-attention |
| **A_4** | Deep structure (long-range) | Long-range dependency | Deep layers |

Each transformer layer adds approximately one degree to the polynomial approximation. A 12-layer transformer ≈ degree-12 polynomial in the token competition signal.

## The Effective Degree

At each generation step, most tokens have negligible probability. Only a small set of D_eff tokens are genuinely competitive. D_eff varies:

- **Factual recall** ("The capital of France is ___"): D_eff ≈ 1. One token dominates. Nearly transitive. **Regime 30 (solvable).**

- **Syntactic choice** ("She ___ to the store"): D_eff ≈ 5-10. A handful of verbs compete. **Regime 36 (pivot).**

- **Creative writing** ("The ___"): D_eff ≈ 100-1000. Many tokens are plausible. **Regime 42 (forbidden).** This is where hallucination lives.

The three regimes of tournament complexity map to three modes of language generation:

| Regime | D_eff | Generation mode | Failure mode |
|--------|-------|-----------------|--------------|
| **30 (solvable)** | 1-3 | Factual, deterministic | Overconfident, rigid |
| **36 (pivot)** | 5-50 | Fluent, natural | Occasional inconsistency |
| **42 (forbidden)** | 100+ | Creative, open-ended | Hallucination, contradiction |

## The Surprise Signal

When an LLM generates a token, the three spectral signals tell us:

**SLOW surprise** (A_1 shift): The generated token doesn't match the unigram/bigram expectation. Low-level mismatch. The model is saying something syntactically unusual.

**MEDIUM surprise** (A_2 shift): The generated token creates unexpected pairwise correlations with recent tokens. The model is building an unusual pattern.

**FAST surprise** (A_3-A_4 shift): The generated token CONTRADICTS the established context. The model is HALLUCINATING — creating a cycle in the "argument tournament" where claim A supports B, B supports C, but C contradicts A.

**Hallucination detection = monitoring the FAST channel.** If the fast signal spikes, the model is generating tokens that create logical cycles. This is the tournament-theoretic definition of hallucination: an intransitivity in the argument structure.

## The Two-Phase Architecture

Current transformers do everything at once: billions of parameters simultaneously SELECT competitive tokens and RANK them. The polynomial framework suggests SEPARATING these:

**Phase 1 — SELECTOR** (the expensive part):
Determine which D_eff tokens are competitive for this context. This is the attention mechanism — pattern-matching over the context window.

**Phase 2 — RANKER** (the cheap part):
Given the D_eff competitive tokens, rank them using a polynomial with D_eff coefficients. This is the output head — a small model applied to a small set.

Current cost: O(V · d_model) per step (all tokens, full embedding dimension).
Two-phase cost: O(context · d_model) for selection + O(D_eff²) for ranking.

If D_eff ≈ 100 and V = 50000: the ranking phase is 500× cheaper. The selection phase is unchanged. Total savings depend on how cheaply we can do selection.

## The Polynomial LLM

The most radical version: replace the entire output head with a degree-D polynomial agent:

```
Input: context embedding (from transformer backbone)
       ↓
[SELECTOR]: identify D_eff competitive tokens
       ↓
[POLYNOMIAL RANKER]: compute D_eff coefficients from the embedding
       ↓
Output: P(z=1) = ranking of competitive tokens
        P(z) for z < 1 = temperature-scaled distribution
        P'(1) = sensitivity (how confident is this ranking?)
        surprise = fast/slow signal ratio
```

The polynomial ranker has D_eff parameters, trained end-to-end. It gives:
- The ranking (at z=1)
- The temperature-scaled ranking (at any z)
- Confidence (from P'(1))
- Hallucination risk (from fast/slow ratio)
- Active learning signal (which comparison would be most informative for fine-tuning)

ALL from a single polynomial evaluation.

## The Deepest Point

An LLM is a machine that, at each step, holds a tournament among 50,000 tokens. The logits are rapidities. The temperature is the polynomial variable. The attention mechanism computes the Walsh decomposition. The layers add polynomial degree.

Current LLMs use brute force: billions of parameters to approximate what is fundamentally a low-degree polynomial in the competitive token space. The polynomial agent suggests that the EXPRESSIVE power of an LLM comes not from the number of parameters, but from the DEGREE of the implicit polynomial — how many simultaneous interactions it can track.

A 12-layer transformer can track 12-way token interactions: the effective degree is approximately the depth. Making the polynomial structure EXPLICIT (rather than implicit in the weight matrices) could lead to:

1. **Provable guarantees**: P(z) is exact, not an approximation
2. **Interpretability**: each coefficient has a clear meaning (prior/evidence/pattern/contradiction/structure)
3. **Efficiency**: O(D_eff) parameters instead of O(V · d_model)
4. **Hallucination detection**: built-in, from the fast channel
5. **Active fine-tuning**: built-in, from the information gain computation

**We are not just ranking tokens. We are holding a tournament among them. And the polynomial P(z) is the most efficient representation of that tournament's complete dynamics.**
