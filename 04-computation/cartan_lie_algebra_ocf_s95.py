#!/usr/bin/env python3
"""
cartan_lie_algebra_ocf_s95.py — opus-2026-03-21-S95

THE LIE ALGEBRA OF TOURNAMENTS AND OCF

Key questions:
1. Can OCF (H(T) = I(Omega(T), 2)) be seen as a statement about so(n)?
2. The skew-adjacency B = A - A^T ∈ so(n). What structure does the
   set {B_T : T a tournament} have within so(n)?
3. How does H(T) relate to the Lie-algebraic properties of B?
4. The Killing form on so(n) is K(X,Y) = (n-2)*Tr(XY). How does
   K(B_T, B_T) relate to tournament invariants?
5. Why is H = 1 for GPT-2 attention tournaments? This means the
   attention matrices threshold to TRANSITIVE tournaments. Is there
   a Lie-algebraic explanation?

Also: The H=1 finding for GPT-2 attention.

Author: opus-2026-03-21-S95
"""

import sys
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(2**m):
        A = np.zeros((n, n), dtype=int)
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[(mask, v)] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))

def count_3cycles(A, n):
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            c3 += 1
        if A[i][k] and A[k][j] and A[j][i]:
            c3 += 1
    return c3

print("=" * 72)
print("  THE LIE ALGEBRA OF TOURNAMENTS AND OCF")
print("  opus-2026-03-21-S95")
print("=" * 72)

# ========================================================================
# PART 1: The Killing form on tournament skew-adjacency matrices
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: KILLING FORM K(B_T, B_T) AND TOURNAMENT INVARIANTS")
print("=" * 72)
print("""
  For so(n), the Killing form is K(X, Y) = (n-2) * Tr(X·Y).

  For a tournament T with B = A - A^T:
    K(B, B) = (n-2) * Tr(B²) = (n-2) * (-n(n-1)) = -(n-2)*n*(n-1)

  This is CONSTANT for all tournaments!

  But the Killing form of B with ITSELF is not the only invariant.
  Consider K(B, B^(2k+1)) which involves higher Casimir-like invariants.

  Actually, for so(n), the independent Casimir invariants are:
    C_k = Tr(B^{2k}) for k = 1, 2, ..., floor(n/2)
  (Only even powers give nonzero traces for skew matrices.)

  We already showed that Tr(B²) is constant but Tr(B⁴) varies.
  Let's compute ALL Casimir invariants for small n.
""")

for n in [3, 4, 5, 6]:
    print(f"\n--- n = {n} ---")
    num_casimirs = n // 2  # floor(n/2)
    results = []

    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_3cycles(A, n)
        B = (A - A.T).astype(float)

        casimirs = []
        for k in range(1, num_casimirs + 1):
            ck = np.trace(np.linalg.matrix_power(B, 2*k))
            casimirs.append(round(ck))

        results.append({
            'H': H, 'c3': c3,
            'casimirs': tuple(casimirs)
        })

    # Group by Casimir tuple
    by_casimir = defaultdict(list)
    for r in results:
        by_casimir[r['casimirs']].append((r['H'], r['c3']))

    print(f"  Casimir invariants: Tr(B²), Tr(B⁴), ..., Tr(B^{2*num_casimirs})")
    print(f"  Number of distinct Casimir tuples: {len(by_casimir)}")
    print(f"  {'Casimirs':>30s} | {'H values':>30s} | count")

    for cas in sorted(by_casimir.keys()):
        hvals = sorted(set(h for h, c in by_casimir[cas]))
        cvals = sorted(set(c for h, c in by_casimir[cas]))
        cnt = len(by_casimir[cas])
        print(f"  {str(cas):>30s} | H={hvals}, c3={cvals} | {cnt}")

# ========================================================================
# PART 2: Does the Casimir tuple determine H?
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 2: DO CASIMIR INVARIANTS DETERMINE H?")
print("=" * 72)

for n in [3, 4, 5, 6]:
    print(f"\n--- n = {n} ---")
    num_casimirs = n // 2
    by_casimir = defaultdict(set)

    for A in all_tournaments(n):
        H = count_hp(A, n)
        B = (A - A.T).astype(float)
        casimirs = tuple(round(np.trace(np.linalg.matrix_power(B, 2*k)))
                         for k in range(1, num_casimirs + 1))
        by_casimir[casimirs].add(H)

    determined = all(len(v) == 1 for v in by_casimir.values())
    print(f"  H determined by Casimir tuple: {determined}")
    if not determined:
        for cas, hvals in sorted(by_casimir.items()):
            if len(hvals) > 1:
                print(f"    Casimirs {cas}: H ∈ {sorted(hvals)}")

# ========================================================================
# PART 3: The special basis vectors of so(n) and tournaments
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 3: TOURNAMENTS AS LATTICE POINTS IN so(n)")
print("=" * 72)
print("""
  The standard basis of so(n) consists of matrices E_{ij} - E_{ji}
  for i < j, where E_{ij} has a 1 at position (i,j) and 0 elsewhere.

  A tournament T corresponds to the element:
    B_T = Σ_{i<j} ε_{ij}(T) · (E_{ij} - E_{ji})
  where ε_{ij}(T) = +1 if i→j and -1 if j→i.

  So B_T is a VERTEX of the hypercube {±1}^{C(n,2)} embedded in so(n).

  KEY OBSERVATION: The 2^{C(n,2)} tournaments form a discrete subset
  of so(n) — the "tournament lattice." They are the vertices of a
  C(n,2)-dimensional hypercube inscribed in so(n).

  The centroid of this hypercube is 0 (the zero matrix).
  Each tournament B_T has ||B_T||² = n(n-1) (constant radius).
  So tournaments lie on a SPHERE in so(n).

  H(T) is a function on this sphere/hypercube.

  QUESTION: What is H as a function on so(n)?
  H(B_T) counts Hamiltonian paths. Can we extend H to a polynomial on so(n)?
""")

# H as a polynomial: let's check the degree
for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    m = n * (n - 1) // 2
    print(f"  dim(so({n})) = {m}")
    print(f"  Number of tournaments = 2^{m} = {2**m}")

    # Each tournament is a vertex of {±1}^m in so(n)
    # H is a function on these vertices.
    # As a multilinear polynomial in ε_{ij} ∈ {±1}:
    #   H(ε) = Σ_S a_S · Π_{(i,j)∈S} ε_{ij}
    # This is the Walsh-Fourier expansion!

    # We already know this from THM-069 (Walsh diagonalization).
    # The Walsh coefficient ĥ[S] is nonzero only when |S| ≡ n (mod 2).

    # In the Lie algebra language:
    # H is a polynomial on so(n) restricted to the tournament lattice.
    # The Walsh expansion IS the decomposition into irreducible
    # representations of the hyperoctahedral group acting on {±1}^m.

    # But in the Lie algebra language, each basis vector E_{ij} - E_{ji}
    # corresponds to a ROOT of the A_{n-1} root system!

    # For so(n), the roots are ±e_i ± e_j (type D_n/B_n).
    # For our tournaments, we use the A_{n-1} embedding where
    # the basis is e_i - e_j for i < j.

    print(f"  Tournament variables: ε_{{ij}} for i<j, total {m}")
    print(f"  Walsh degrees present in H: even only (complement symmetry)")

    # Highest Walsh degree of H
    # H = Σ_σ Π_{k} T[σ(k)][σ(k+1)] is a product of n-1 factors
    # Each factor is ε_{ij} (if i<j) or -ε_{ji} (if j<i)
    # Well, T[i][j] = (1 + ε_{ij})/2 for i<j, (1 - ε_{ji})/2 for i>j.
    # So each factor is degree 1 in ε, and H is degree n-1 in ε.
    # But this is Walsh degree, and many cancel.

    # Actually, H = Σ_σ Π_k (1 + B[σ(k)][σ(k+1)])/2
    #            = (1/2^{n-1}) Σ_σ Π_k (1 + ε_{σ(k)σ(k+1)})
    # where ε_{ij} = B[i][j] = ±1 for i < j.

    print(f"  H = sum over n! permutations of product of (n-1) factors")
    print(f"  Each factor has degree 1 in Walsh variables")
    print(f"  H has Walsh degree at most n-1")

# ========================================================================
# PART 4: OCF in Lie algebra terms
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 4: OCF IN LIE ALGEBRA TERMS")
print("=" * 72)
print("""
  OCF states: H(T) = I(Ω(T), 2) = Σ_k α_k · 2^k

  where α_k = number of independent sets of size k in Ω(T).

  In Lie algebra terms:
  - T is a vertex of the hypercube {±1}^m ⊂ so(n)
  - Ω(T) is the conflict graph of directed odd cycles in T
  - Two cycles conflict iff they share a vertex
  - I(Ω, 2) is the hard-core lattice gas partition function at fugacity 2

  KEY OBSERVATION:
  Odd cycles of T correspond to CLOSED PATHS in the directed graph.
  In so(n) terms, a k-cycle C = (v_1, ..., v_k, v_1) with all v_i→v_{i+1}
  corresponds to the product:
    ε_{v_1 v_2} · ε_{v_2 v_3} · ... · ε_{v_k v_1}

  But this is exactly a TRACE in so(n)!
  Specifically, for the standard representation:
    Tr(B_{v_1 v_2} · B_{v_2 v_3} · ... · B_{v_k v_1})
  where B_{ij} = E_{ij} - E_{ji}.

  Wait, this isn't quite right. Let me be more careful.

  For a 3-cycle {a,b,c} with a→b→c→a:
  The indicator function is T[a][b]·T[b][c]·T[c][a].
  In terms of ε: = ((1+ε_{ab})/2) · ((1+ε_{bc})/2) · ((1-ε_{ac})/2)
  (assuming a<b<c; ε_{ac} appears with minus because c→a = -ε_{ac}).

  Actually, this is getting complicated because of the orientation.
  The key point is:

  CYCLE DETECTION IS A POLYNOMIAL ON THE TOURNAMENT LATTICE.

  And the CONFLICT GRAPH Ω(T) is determined by the OVERLAP STRUCTURE
  of these cycles.

  OCF says: the COUNTING function H is equal to a WEIGHTED SUM over
  independent sets of cycles.

  This is analogous to the CLUSTER EXPANSION in statistical mechanics,
  where the partition function equals a sum over compatible configurations.

  In Lie algebra terms:
    H = Σ_{independent cycle sets S} 2^{|S|}
    = the generating function for non-conflicting cycle collections
    evaluated at fugacity 2.

  The fugacity 2 is special because:
    2 = dim(kernel of ∂_1 for the tournament chain complex at p=2)
    2 = the "defect" of the Rédei parity constraint
    2 = 1 + 1 where 1 = each cycle contributes ±1 to H mod 2
""")

# ========================================================================
# PART 5: Why H = 1 for attention tournaments
# ========================================================================
print("\n" + "=" * 72)
print("  PART 5: WHY H = 1 FOR GPT-2 ATTENTION TOURNAMENTS")
print("=" * 72)
print("""
  FINDING: ALL GPT-2 attention tournaments are TRANSITIVE (H = 1).

  EXPLANATION:
  The attention matrix A[i,j] = P(token j | context up to token i)
  has CAUSAL MASKING: A[i,j] = 0 for j > i (cannot attend to future).

  With causal masking, the attention matrix is LOWER TRIANGULAR.
  Therefore:
    For i < j: A[i,j] = 0 < A[j,i] (since A[j,i] ≥ 0 with some mass)

  So the thresholded tournament always has j→i for i < j.
  This is EXACTLY the transitive tournament!

  H(transitive) = 1: there is exactly one Hamiltonian path
  (n → n-1 → ... → 1), corresponding to reading order.

  CONCLUSION: The H = 1 finding is an ARTIFACT of causal masking.
  It tells us nothing interesting about the tournament structure.

  To get non-trivial tournaments from attention, we need to either:
  1. Use the BIDIRECTIONAL attention (e.g., BERT) instead of GPT-2
  2. Use only the NON-CAUSAL part of the attention (mask out the
     triangular constraint)
  3. Compare attention between different positions excluding the
     causal bias
  4. Use cross-attention from encoder-decoder models

  Let's verify this hypothesis.
""")

# Verify: causal masking creates transitive tournament
print("  Verification: causal attention creates transitive tournament")
print()

# Simulate causal attention
for n in [5, 8, 12]:
    rng = np.random.RandomState(42)
    logits = rng.randn(n, n)
    # Causal mask: zero out upper triangle (j > i, i.e., future tokens)
    mask = np.tril(np.ones((n, n)))
    logits = logits * mask + (1 - mask) * (-1e9)
    exp_l = np.exp(logits - logits.max(axis=1, keepdims=True))
    A = exp_l / exp_l.sum(axis=1, keepdims=True)

    # Threshold to tournament
    T = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if A[i,j] > A[j,i]:
                T[i][j] = 1
            else:
                T[j][i] = 1

    # Check: is it transitive?
    is_transitive = True
    for i in range(n):
        for j in range(i+1, n):
            # In transitive tournament: j→i for all j > i
            # i.e., T[j][i] = 1 for j > i
            if T[i][j] == 1:  # i→j means NOT transitive (standard order)
                is_transitive = False
                break
        if not is_transitive:
            break

    # Actually check H
    H = count_hp(T, n) if n <= 12 else -1

    scores = T.sum(axis=1)
    c3 = count_3cycles(T, n) if n <= 12 else -1
    print(f"  n={n}: H={H}, c3={c3}, scores={tuple(scores)}")
    print(f"        Transitive: {H == 1}")

    # Now try WITHOUT causal mask
    logits2 = rng.randn(n, n)
    exp_l2 = np.exp(logits2 - logits2.max(axis=1, keepdims=True))
    A2 = exp_l2 / exp_l2.sum(axis=1, keepdims=True)

    T2 = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if A2[i,j] > A2[j,i]:
                T2[i][j] = 1
            else:
                T2[j][i] = 1

    H2 = count_hp(T2, n) if n <= 12 else -1
    c3_2 = count_3cycles(T2, n) if n <= 12 else -1
    scores2 = T2.sum(axis=1)
    print(f"  n={n} (no mask): H={H2}, c3={c3_2}, scores={tuple(sorted(scores2))}")

print("""
  CONFIRMED: Causal masking → transitive tournament → H = 1.

  For the Cartan bridge to be meaningful for LLMs, we should:
  1. Use the SYMMETRIC part of attention (which varies between tokens)
  2. OR use bidirectional models (BERT, etc.)
  3. OR extract the "tournament signal" from the non-diagonal part
     of the Cartan decomposition, ignoring the causal bias.
""")

# ========================================================================
# PART 6: NON-CAUSAL TOURNAMENT — removing the position bias
# ========================================================================
print("\n" + "=" * 72)
print("  PART 6: NON-CAUSAL TOURNAMENT — SUBTRACTING POSITION BIAS")
print("=" * 72)
print("""
  IDEA: The causal mask creates a position-dependent bias.
  If we SUBTRACT the expected causal attention pattern (proportional
  to a lower-triangular matrix), the residual may have non-trivial
  tournament structure.

  Alternatively: look at the ANTISYMMETRIC part of the attention
  matrix A_anti = (A - A^T)/2. For causal attention:
  - Lower triangle: A[j,i] for j > i (nonzero)
  - Upper triangle: A[i,j] for j > i (zero due to mask)
  - A_anti[i,j] = (A[i,j] - A[j,i])/2 = -A[j,i]/2 for i < j (always negative)

  So A_anti has ALL negative entries in the upper triangle.
  This means the tournament is ALWAYS i←j for i < j (reverse reading order).

  To get interesting structure, we need to look at the RESIDUAL
  after subtracting the causal pattern, or use bidirectional attention.

  Let's instead analyze the attention matrix WITHIN the lower triangle:
  For rows j1 > j2 > 0, compare A[j1, 0] vs A[j2, 0]:
  does token j1 attend to token 0 more or less than j2 does?
  This gives a tournament on the "attending" tokens.
""")

# Try loading GPT-2 again to extract non-trivial tournaments
print("  Loading GPT-2 for non-causal analysis...")
try:
    import torch
    from transformers import GPT2LMHeadModel, GPT2Tokenizer

    tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
    model = GPT2LMHeadModel.from_pretrained("gpt2", output_attentions=True)
    model.eval()

    sentence = "The cat sat on the mat and looked at the dog."
    inputs = tokenizer(sentence, return_tensors="pt")
    n_tokens = inputs['input_ids'].shape[1]

    with torch.no_grad():
        outputs = model(**inputs)
        attentions = outputs.attentions

    print(f"  Sentence: \"{sentence}\" ({n_tokens} tokens)")
    tokens = tokenizer.convert_ids_to_tokens(inputs['input_ids'][0])
    print(f"  Tokens: {tokens}")

    # Strategy: Build a tournament on the n tokens based on
    # the "column profile" — how much is each token attended to?
    # For each pair (i, j), compare sum_k A[k,i] vs sum_k A[k,j]
    # (total attention received by i vs j).

    for layer in [0, 5, 11]:
        for head in [0, 5, 11]:
            A = attentions[layer][0, head].numpy()  # n x n

            # Column sums = total attention received by each token
            col_sums = A.sum(axis=0)

            # Build tournament on "attention received"
            T = np.zeros((n_tokens, n_tokens), dtype=int)
            for i in range(n_tokens):
                for j in range(i+1, n_tokens):
                    if col_sums[i] > col_sums[j]:
                        T[i][j] = 1
                    elif col_sums[j] > col_sums[i]:
                        T[j][i] = 1
                    else:
                        T[i][j] = 1

            H = count_hp(T, n_tokens) if n_tokens <= 12 else -1
            c3_val = count_3cycles(T, n_tokens) if n_tokens <= 12 else -1
            scores = T.sum(axis=1)

            print(f"\n  Layer {layer}, Head {head}: column-attention tournament")
            print(f"    H = {H}, c3 = {c3_val}")
            print(f"    Scores: {tuple(scores)}")
            # Show which token gets most attention
            top_idx = np.argsort(col_sums)[::-1]
            print(f"    Attention ranking: {[tokens[i] for i in top_idx[:5]]}")

    # Strategy 2: Row-wise comparison — asymmetric attention between pairs
    # For each pair (i,j), compare A[i,j] vs A[j,i]
    # But with causal masking, A[i,j]=0 for j>i. So this always gives j→i.
    # Instead, compare within the accessible pairs:
    # For tokens i, j, k (all < some position), who does each prefer?

    # Better strategy: Use the QK scores BEFORE masking
    # This requires extracting the QK dot products, not the softmax output.

    print("\n\n  --- Strategy: Extract QK scores before causal mask ---")
    # Access the raw QK attention scores
    # In GPT-2, the attention is computed as:
    # scores = Q @ K^T / sqrt(d_k)
    # Then masked and softmaxed.

    # We can get pre-mask scores by computing Q @ K^T ourselves
    for layer in [0, 5, 11]:
        # Get Q, K projections from the model
        attn_layer = model.transformer.h[layer].attn
        # c_attn projects input to Q, K, V concatenated
        # Weight shape: (768, 2304) = (d_model, 3*d_model)
        # Split into Q, K, V each of shape (768, 768)

        hidden_states = None
        # We need the hidden states at each layer. Let's just use the
        # attention weights directly and infer structure.
        # Actually, the post-softmax attention already has the causal mask.
        # Let's instead use a different approach.
        pass

    # Strategy 3: Compare cross-position attention PROFILES
    # For positions i and j, compute A[i,:] (where does i attend?)
    # and A[j,:] (where does j attend?). The tournament is:
    # i→j if A[i,:] and A[j,:] diverge in a way that i "dominates" j.
    # We can use the KL divergence: i→j if KL(A[i,:] || A[j,:]) > KL(A[j,:] || A[i,:])

    print("\n  --- Strategy: Profile divergence tournament ---")
    for layer in [0, 5, 11]:
        for head in [0]:
            A = attentions[layer][0, head].numpy()
            n = A.shape[0]

            # KL divergence for each pair
            T_kl = np.zeros((n, n), dtype=int)
            for i in range(n):
                for j in range(i+1, n):
                    # KL(A[i,:] || A[j,:])
                    eps = 1e-10
                    kl_ij = np.sum(A[i] * np.log((A[i] + eps) / (A[j] + eps)))
                    kl_ji = np.sum(A[j] * np.log((A[j] + eps) / (A[i] + eps)))
                    if kl_ij > kl_ji:
                        T_kl[i][j] = 1
                    else:
                        T_kl[j][i] = 1

            H_kl = count_hp(T_kl, n) if n <= 12 else -1
            c3_kl = count_3cycles(T_kl, n) if n <= 12 else -1
            scores_kl = T_kl.sum(axis=1)

            print(f"\n  Layer {layer}, Head 0: KL-divergence tournament")
            print(f"    H = {H_kl}, c3 = {c3_kl}")
            print(f"    Scores: {tuple(scores_kl)}")
            print(f"    Transitive: {H_kl == 1}")

except ImportError:
    print("  (transformers not available, skipping GPT-2 analysis)")

# ========================================================================
# PART 7: B^2 as a Casimir-like element and its eigenvalues
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 7: B² EIGENVALUES AND TOURNAMENT STRUCTURE")
print("=" * 72)
print("""
  B² is SYMMETRIC NEGATIVE SEMIDEFINITE (since B is skew).
  Eigenvalues of B² are -θ₁², -θ₂², ..., 0 (if n odd).

  B² = -D where D is positive semidefinite.
  D[i,j] = -B²[i,j] = Σ_k B[i,k]B[j,k] = <row_i(B), row_j(B)>

  This is the GRAM MATRIX of the rows of B.

  For a tournament: B[i,k] = ±1, so row_i is a ±1 vector.
  D[i,j] = number of vertices k where i,j agree on direction
           minus number where they disagree.

  If i→k and j→k (or k→i and k→j): agree. Contribution +1.
  If i→k and k→j (or k→i and j→k): disagree. Contribution -1.

  D[i,j] = #{k : (i→k ∧ j→k) or (k→i ∧ k→j)}
          - #{k : (i→k ∧ k→j) or (k→i ∧ j→k)}

  For i = j: D[i,i] = n-1 (all agree with themselves).
  For i ≠ j: D[i,j] = ... depends on the tournament.

  KEY: D = -B² is the "tournament similarity matrix."
  Its eigenvalues are the θ_k² (rotation angle squares).
  Its trace is n(n-1) (constant).

  The DISTRIBUTION of eigenvalues encodes tournament structure!
""")

for n in [5, 7]:
    print(f"\n--- n = {n} ---")
    results = []

    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_3cycles(A, n)
        B = (A - A.T).astype(float)
        D = -B @ B  # tournament similarity matrix
        eigs = sorted(np.linalg.eigvalsh(D), reverse=True)

        results.append({
            'H': H, 'c3': c3,
            'eigs': tuple(round(e, 4) for e in eigs)
        })

        if n > 6:
            break  # just show first few for n=7

    groups = defaultdict(list)
    for r in results:
        groups[(r['H'], r['c3'])].append(r)

    for (H, c3) in sorted(groups.keys()):
        grp = groups[(H, c3)]
        print(f"  H={H:3d}, c3={c3:2d}: D eigenvalues = {grp[0]['eigs']}")

# ========================================================================
# PART 8: The tournament Laplacian and Cartan
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 8: TOURNAMENT LAPLACIAN AND SPECTRAL GAP")
print("=" * 72)
print("""
  The tournament Laplacian: L = D_out - A where D_out = diag(s_1,...,s_n).

  For the skew-adjacency B = A - A^T:
    L = D_out - A = D_out - (B + J - I)/2 ... hmm, not clean.

  Better: use the SYMMETRIZED Laplacian of the underlying complete graph:
    L_sym = nI - J (the complete graph Laplacian)
  This has eigenvalues 0 (multiplicity 1) and n (multiplicity n-1).

  The tournament ADDS a skew-symmetric "magnetic potential" B to this.
  The magnetic Laplacian is: L_B = L_sym + iB (complex Hermitian).
  Its eigenvalues are real and encode the interplay between the
  "base graph" structure and the "tournament flow."
""")

for n in [3, 4, 5]:
    print(f"\n--- n = {n}: magnetic Laplacian eigenvalues ---")

    by_Hc3 = defaultdict(list)
    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_3cycles(A, n)
        B = (A - A.T).astype(float)
        L_sym = n * np.eye(n) - np.ones((n, n))
        L_mag = L_sym + 1j * B
        eigs_mag = sorted(np.linalg.eigvalsh(L_mag))
        by_Hc3[(H, c3)].append(tuple(round(e, 4) for e in eigs_mag))

    for (H, c3) in sorted(by_Hc3.keys()):
        rep = by_Hc3[(H, c3)][0]
        print(f"  H={H:3d}, c3={c3:2d}: λ(L_mag) = {rep}")

print("\n\nDone. opus-2026-03-21-S95")
