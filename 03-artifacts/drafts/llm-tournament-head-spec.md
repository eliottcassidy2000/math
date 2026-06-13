# LLM Tournament Head: Rigorous Specification

**Author:** opus-2026-03-17-S74b
**Status:** Draft. Based on proved theorems only.
**Audit:** polynomial_head_audit.py identifies errors in prior prototypes.

---

## What This IS

A **confidence estimation module** for LLM output heads, using exact
tournament-theoretic invariants. Not a "replacement" for the output head
(you still need embeddings to score vocabulary tokens), but an augmentation
that provides mathematically grounded uncertainty quantification.

## What This is NOT

- NOT "5 coefficients encode everything" (false for n > 12, see audit)
- NOT "160K parameters replace 38.6M" (parameter count was wrong)
- NOT "intelligence is a degree-4 polynomial" (philosophy, not math)

---

## The Three Innovations (all backed by proved theorems)

### Innovation 1: Multi-Context Intransitivity Detection

**The idea:** When an LLM's different layers (or attention heads, or
prompt variations) produce different logit orderings for the next token,
the majority-vote tournament among top-k candidates may have **cycles**.

**Why this matters:** Cycles in the tournament = Condorcet paradoxes =
the model's internal representations DISAGREE about the ranking. This is
a precise, measurable signal of uncertainty that goes beyond softmax entropy.

**The mathematics (all proved):**
- OCF (THM-077): H(T) = 1 + 2α₁ + 4α₂ + 8α₃ + ...
- α₁ = number of directed odd cycles (3-cycles + 5-cycles + ...)
- H = 1 ↔ tournament is transitive (all contexts agree on ordering)
- H > 1 ↔ intransitivity detected (contexts disagree)
- Confidence = 1/H (exact, not heuristic)

**Computational cost:** O(k³) for 3-cycle counting among top-k tokens.
For k = 16: ~4000 operations. For k = 32: ~33000 operations.
Negligible compared to transformer forward pass.

**Verified:** OCF matches Held-Karp exactly at n=5 (1024/1024 tournaments).

### Innovation 2: Arctanh Evidence Accumulation

**The idea:** Replace exponential softmax attention with hyperbolic
arctanh-based evidence. Evidence from different layers/heads ADDS
linearly in rapidity space.

**The mathematics (proved):**
- Formal group: F(x,y) = (x+y)/(1+xy) (tanh addition law)
- arctanh linearizes: arctanh(F(x,y)) = arctanh(x) + arctanh(y)
- arctanh is the unique odd function with rational exponential (THM)
- The Cayley transform Q(x) = (1+x)/(1-x) satisfies Q(F(x,y)) = Q(x)·Q(y)

**Practical advantage:**
- No overflow/underflow (arctanh bounded)
- Evidence from different sources is additive
- Temperature control = scaling rapidity (linear, not exponential)

**Caveat:** In a single attention layer, arctanh(tanh(x)) = x, so
single-layer behavior is identical to standard attention. The advantage
appears in MULTI-LAYER accumulation.

### Innovation 3: Walsh Sparsity for Output Compression

**The idea:** The ranking function H has only O(n²) nonzero Walsh
components out of 2^{C(n,2)} possible (THM-069, proved). This means
the output distribution, viewed as a function of pairwise comparisons,
lives on a tiny subspace.

**Compression factors (proved):**
- n=5: 3 nonzero out of 1024 (341× compression)
- n=7: ~20 nonzero out of 2,097,152 (~100,000× compression)

**Application:** Instead of storing full logit vectors, represent the
output distribution in the Walsh basis. For top-k = 7, only ~20 numbers
(Walsh coefficients at even-degree path-union edge sets) suffice to
reconstruct the complete ranking probabilities.

---

## Architecture

```
Input: hidden_state h ∈ R^d_model

                    ┌──────────────────────┐
                    │  Standard lm_head    │
                    │  h → W·h = logits    │  (standard, unchanged)
                    └──────────┬───────────┘
                               │
                    ┌──────────▼───────────┐
                    │  Top-k Selection     │
                    │  Select k candidates │  O(V log k)
                    └──────────┬───────────┘
                               │
            ┌──────────────────┼──────────────────┐
            │                  │                  │
   ┌────────▼────────┐ ┌──────▼──────┐  ┌────────▼────────┐
   │ Multi-Context   │ │  Spectral   │  │ Walsh Sparsity  │
   │ Tournament      │ │  Gap        │  │ Compression     │
   │ (layers→cycles) │ │  (top1-top2)│  │ (O(k²) coeffs) │
   └────────┬────────┘ └──────┬──────┘  └────────┬────────┘
            │                  │                  │
            └──────────────────┼──────────────────┘
                               │
                    ┌──────────▼───────────┐
                    │  Diagnostics:        │
                    │  - confidence (1/H)   │
                    │  - t3 (3-cycle count) │
                    │  - spectral_gap      │
                    │  - regime (certain/   │
                    │    choosing/confused) │
                    └──────────────────────┘
```

## Integration Points

1. **Chatbot Arena / LLM evaluation:** Build tournament from pairwise
   human preferences. OCF confidence tells you exactly how many
   consistent total rankings exist. H > 1 = ranking is ambiguous.

2. **Inference-time uncertainty:** Extract logits from different contexts
   (attention heads, prompt paraphrases, or MC dropout). Build tournament
   via majority vote. OCF confidence measures ranking ambiguity.
   **NOTE (S74b empirical finding):** Using transformer LAYERS as contexts
   gives a NULL result — layer disagreement measures processing depth
   (positive signal), not uncertainty. Better context sources needed:
   attention heads within a layer, or multiple prompts.

3. **Speculative decoding:** Use spectral gap to decide whether to
   accept draft tokens. Large gap = accept immediately (confident).
   Small gap = verify with full model.

4. **Training signal:** Intransitivity ratio as a training loss term.
   Penalize models where different layers strongly disagree on token
   ordering (encourages internal consistency).

---

## Implementation: tournament_output_head.py

- `OCFConfidence`: Exact OCF computation (3-cycles + 5-cycles + disjoint pairs)
- `ArctanhAttention`: Formal-group-based attention layer
- `TournamentOutputHead`: PyTorch module with confidence diagnostics
- `tournament_from_multi_context()`: Build tournament from multiple logit vectors

Tested: 0 errors on exhaustive n=5 verification (1024 tournaments).

---

## Empirical Results (opus-S74b, GPT-2)

**Experiment 1: Layer-based tournaments (173 token predictions)**
- Using 7 transformer layers as "contexts" for majority-vote tournament
- Result: **NULL** — layer disagreement does NOT predict lower accuracy
  - Transitive (H=1): 39.0% accuracy, n=159
  - Intransitive (H>1): 50.0% accuracy, n=14
  - Layer disagreement measures PROCESSING DEPTH, not uncertainty
  - When late layers change the ranking, it means MORE refinement (positive)
- OCF diagnostics work correctly; regime classification functions
- NOTE: Initial "positive" result was an artifact of applying ln_f twice
  to the final hidden state (bug in transformers v5+ compatibility)

**Experiment 2: Staged evaluation**
- TournamentHead with hot-256 cache matches standard predictions (7/8)
- One mismatch from unwarmed cache (correct after warmup)
- Generation is coherent: "where the velocity of the object is..."
- Hot cache needs initialization from token frequency distribution

**Key insight:** The multi-context tournament approach is sound but the
CHOICE OF CONTEXTS matters critically. Layers are NOT good contexts
(they measure refinement, not disagreement). Better sources:
1. Attention heads within the same layer (parallel perspectives)
2. Multiple prompt paraphrases (sensitivity to framing)
3. MC dropout passes (model uncertainty sampling)
4. Chatbot Arena-style human preferences (ground truth)

## What Needs to Be Done Next

1. **Attention head disaggregation:** Extract logits from individual
   attention heads (GPT-2 has 12 per layer). Build tournament from
   head-level predictions. These represent genuinely different perspectives.

2. **Prompt paraphrase tournaments:** For the same question, generate
   multiple rephrasings. Build tournament from each rephrasing's top-k.
   OCF confidence = sensitivity to framing (hallucination risk).

3. **arctanh attention training:** Train a small model with ArctanhAttention
   and compare perplexity with standard softmax attention.

4. **Scaling analysis:** OCF computation is O(k³) for 3-cycles. For
   k = 32 this is fast. For k = 1000 it's slow. Use matrix trace method:
   t₃ = Tr(A³)/3 which is O(n³) via BLAS and works for any n.

5. **Hot token cache initialization:** Initialize from corpus unigram
   frequency instead of tokens 0-255. This should enable >80% early exits.
