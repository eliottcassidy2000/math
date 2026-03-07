# THM-068: Position Character Decomposition

**Type:** Theorem
**Certainty:** 5 -- PROVED (algebraic, all degrees, all odd n)
**Status:** PROVED
**Added by:** opus-2026-03-07-S35
**Tags:** #Walsh-Hadamard #transfer-matrix #position #diagonal #character

---

## Statement

### (i) Vertex decomposition of transfer matrix diagonal

For odd n, each diagonal entry M[v,v] of the transfer matrix, viewed as a function on the tournament hypercube {0,1}^m (m = C(n,2)), satisfies:

**M[v,v]_hat[S] = (-1)^{[v in N(S)]} * H_hat[S] / (n - 2k)**

for each nonzero Walsh monomial S of degree 2k, where:
- H_hat[S] is the Walsh coefficient of H(T) = tr(M)
- N(S) is a "negative set" of k vertices determined by the monomial S
- The formula sum_v M[v,v]_hat[S] = H_hat[S] follows since |N(S)| = k

### (ii) Even-length component rule (PROVED)

**H_hat[S] ≠ 0 only if every path component of the subgraph S has even length.**

The subgraph induced by |S| = 2k edges (with max degree ≤ 2, required to be embeddable in a Hamiltonian path) is a union of vertex-disjoint paths. Components of odd length L give desc d and L-d from the two traversal directions; since L is odd, these have opposite parity, causing exact cancellation.

For even-length L = 2a: both directions give desc with the same parity, so (-1)^desc is constant.

### (iii) Negative set structure

The negative set N(S) has |N(S)| = k for any monomial S of Walsh degree 2k.

**Degree 0 (k=0):** N(S) = {} (empty). All vertices contribute equally: M[v,v]_hat = H_hat/n.

**Degree 2 (k=1):** N(S) = {shared vertex}. For edge pair S = {(a,b),(a,c)}, the shared vertex a is the unique negative vertex. M[v,v]_hat = ±H_hat/(n-2).

**Degree 4 (k=2):** Two types: P4 (path on 5 vertices) and P2+P2 (two disjoint wedges). For P4: N = {odd-position vertices of the path}. For P2+P2: N = {center vertex of each wedge}.

**Degree n-1 (k=(n-1)/2, top degree):** Only type is P_{n-1} (Hamiltonian path). N(S) = {odd-position vertices}. M[v,v]_hat = (-1)^{pos(v,path)} * H_hat[S]. No scaling (1/(n-2k) = 1).

### (iv) Pattern counts

At Walsh degree 2k, there are exactly C(n,k) distinct position character patterns, corresponding to the C(n,k) choices of k-element negative set from n vertices.

At n=5: C(5,0) = 1, C(5,1) = 5, C(5,2) = 10 patterns at degrees 0, 2, 4.

### (v) Energy identity

The Parseval energy at each Walsh degree level satisfies:

**sum_v sum_{deg(S)=2k} M[v,v]_hat[S]^2 = n/(n-2k)^2 * sum_{deg(S)=2k} H_hat[S]^2**

This follows from sum_v (-1)^{2*[v in N(S)]} = n.

### (vi) Connection to THM-053 and THM-071

This theorem bridges:
- **THM-053** (M[v,v] = sum_P (-1)^{pos(v,P)}): the pos function is encoded in the highest Walsh degree (n-1), where the position character equals (-1)^{pos(v,path)} with no scaling.
- **THM-071** (Walsh-Fourier diagonalization): H decomposes into pure Walsh degree components D_{2k}. The PCD refines this by showing M[v,v] decomposes identically, with each D_{2k} split into vertex-dependent signed components.

The scaling 1/(n-2k) increases from 1/n (degree 0) to 1 (top degree), meaning the highest Walsh degree carries the MOST per-vertex information about position parity, while lower degrees carry progressively averaged information.

---

## Complete algebraic proof (all degrees, all odd n)

### Setup

Let S be a Walsh monomial of degree 2k whose subgraph has r vertex-disjoint
even-length path components P_{2a_1}, ..., P_{2a_r} with sum a_i = k.

Define: valid permutation sigma = one whose Hamiltonian path contains all edges of S.

The "master formula" is:

**M[v,v]_hat[S] = (1/2^{n-1}) sum_sigma (-1)^{pos(v,sigma) + desc(sigma,S)}**

### Step 1: Valid permutation structure

Each component P_{2a_i} must appear as a consecutive block in the HP (it's a connected
subpath). A valid permutation is determined by:
- A **macro-permutation** of n-2k items: r blocks (of widths 2a_1+1, ..., 2a_r+1)
  and n-2k-r free vertices (of width 1)
- An **internal orientation** for each block (2 choices per block)

Total valid perms = **2^r * (n-2k)!**

### Step 2: Constant descent sign

For an even-length path P_{2a} on vertices v_0,...,v_{2a}:
- Forward traversal gives desc = d
- Reverse gives desc = 2a - d
- (-1)^d = (-1)^{2a-d} since 2a is even

So both orientations give the same (-1)^{desc_i}. The macro-permutation does not
affect which S-edges are traversed ascending vs descending. Therefore:

**epsilon := (-1)^{desc(sigma,S)} is constant across all valid sigma.**

(The value of epsilon depends on vertex labels but cancels in all ratios.)

### Step 3: Odd-width parity lemma

Every macro-item has **odd width**: blocks have width 2a_i+1 (odd), free vertices
have width 1 (odd). Therefore the start position of the item at macro-position j is:

**start(j) = sum of widths of items at positions 0,...,j-1 ≡ j (mod 2)**

since each summand contributes 1 mod 2.

For vertex v in block B_i at internal offset d (0-indexed from one endpoint):

**pos(v) ≡ macro_pos(B_i) + d (mod 2)**

For a free vertex w at macro-position j: pos(w) ≡ j (mod 2).

### Step 4: Signed position sum

Fix a vertex v in component i at internal offset d. For each macro-permutation,
B_i is at some macro-position j among n-2k positions. There are (n-2k-1)!
macro-perms with B_i at position j. Internal orientations contribute factor 2^r
(and both orientations of B_i give the same parity since the path length is even).

Raw signed sum for v:

**R(v) = epsilon * 2^r * (n-2k-1)! * sum_{j=0}^{n-2k-1} (-1)^{j+d}**
**     = epsilon * 2^r * (n-2k-1)! * (-1)^d * sum_{j=0}^{n-2k-1} (-1)^j**

For odd n: n-2k is odd, so the alternating sum over n-2k terms equals 1.

**R(v) = epsilon * 2^r * (n-2k-1)! * (-1)^d**

For free vertices (d=0): R(w) = epsilon * 2^r * (n-2k-1)!.

### Step 5: PCD formula

The total H_hat * 2^{n-1} = sum_v R(v). Count: k negative vertices (odd d) contribute
(-1) each, (n-k) positive/free vertices contribute (+1) each.

**H_hat * 2^{n-1} = epsilon * 2^r * (n-2k-1)! * (n-2k) = epsilon * 2^r * (n-2k)!**

Ratio: **R(v) / (H_hat * 2^{n-1}) = (-1)^d / (n-2k)**

This gives: **M[v,v]_hat[S] = (-1)^{[v in N(S)]} * H_hat[S] / (n-2k)**

where **N(S) = {vertices at odd internal offsets}**, with |N(S)| = sum a_i = k. **QED.**

### Degree-specific instances

**Degree 0 (k=0):** Trivially M[v,v]_hat = H_hat/n by S_n symmetry.

**Degree 2 (k=1, P_2):** Single block of width 3. Center vertex has offset 1 (negative).
R(center) = -epsilon*2*(n-3)!, R(other) = +epsilon*2*(n-3)!. Ratio = ∓1/(n-2).

**Degree 4 (k=2):**
- P_4: Single block of width 5. N = {offset 1, offset 3} = odd-position vertices.
  Total = 2*(n-4)!. Ratio = ±1/(n-4).
- P_2+P_2: Two blocks of width 3. N = {center of each wedge}.
  Total = 4*(n-4)!. Ratio = ±1/(n-4).

**Degree n-1 (k=(n-1)/2, top):** Single Hamiltonian path block.
N = {odd-position vertices}. n-2k = 1, so ratio = ±1. No averaging.

### Even-length component rule (proof)

If S has an odd-length component P_L (L odd), then its two traversal directions
give desc d and L-d with (-1)^d = -(-1)^{L-d}. Since both orientations appear
equally in the sum, the contributions cancel exactly. Hence H_hat[S] = 0.

## Exhaustive verification (n=3, n=5)

Computed M[v,v] for all 2^m tournaments and all n vertices. Walsh-Hadamard transformed each M[v,v] function. Verified:
1. All odd Walsh degrees vanish (only degrees 0, 2, ..., n-1 are nonzero)
2. At each degree 2k, all coefficients are ±H_hat[S]/(n-2k)
3. The negative set has cardinality k
4. The pattern counts match C(n,k)

## Numerical verification of general proof

Verified the algebraic formulas against permutation enumeration:
- n=5: P2 (deg 2), P4 (deg 4) -- all pass
- n=7: P2 (deg 2), P4 (deg 4), P2+P2 (deg 4), P6 (deg 6) -- all pass
- n=9: P4+P2 (deg 6), P2+P2+P2 (deg 6), P6 (deg 6) -- all pass

---

## Significance

This result reveals that the signed position function pos(v,P) from Rédei's theorem is not merely a counting device — it is the TOP LEVEL of a hierarchical Walsh decomposition of the transfer matrix diagonal. Each level captures progressively finer vertex-edge correlations:

- Level 0: uniform average (E[H])
- Level 1: vertex-edge pair correlations (3-cycle structure)
- Level k: vertex-vs-2k-edge-subgraph correlations

The hierarchy terminates at k = (n-1)/2 where it recovers the full signed position function.

---

## Scripts

- `04-computation/pcd_verification.py` -- Comprehensive verification at n=3,5,7
- `/tmp/pos_walsh_bridge.py` -- Initial discovery of even-degree-only property for M[v,v]
- `/tmp/pos_walsh_formula.py` -- Discovery of degree-2 and degree-4 formulas
- `/tmp/pcd_algebraic.py` -- Algebraic proof of degree-2 case
- `/tmp/pcd_degree4_other_types.py` -- Even-length component rule discovery
- `/tmp/pcd_even_length_proof.py` -- Proof of descent cancellation for odd-length components
