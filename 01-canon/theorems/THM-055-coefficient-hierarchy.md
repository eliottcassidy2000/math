# THM-055: Coefficient Hierarchy Theorem

**Type:** Theorem (algebraic proof + computational verification)
**Certainty:** 4 -- PROVED algebraically for k=0,1; verified computationally for k=2
**Status:** PROVED (k=0,1), VERIFIED (k=2)
**Added by:** opus-2026-03-06-S11b (continued^4)
**Tags:** #transfer-matrix #r-polynomial #moments #hierarchy

---

## Statement

**Theorem.** For a tournament T on n vertices (n odd), let M(r) = sum_{k=0}^{(n-1)/2} c_{2k} r^{2k} be the even r-polynomial of the transfer matrix. Then:

**(a)** tr(c_{n-1-2k}) = sum_P e_{2k}(s_P) over all n! permutations P, where s_i = A[p_i, p_{i+1}] - 1/2.

**(b)** Each e_{2k}(s_P) is a polynomial of degree 2k in f_P = (number of forward arcs in P), with coefficients depending only on n. Explicitly, via Newton's identities and the key simplification that s_i^{2j} = 1/4^j:

At n=7:
- e_0 = 1
- e_2(f) = f^2/2 - 3f + 15/4
- e_4(f) = f^4/24 - f^3/2 + 47f^2/24 - 11f/4 + 15/16
- e_6(f) = f^6/720 - f^5/40 + 49f^4/288 - 13f^3/24 + 287f^2/360 - 13f/30 + 1/64

**(c)** Consequently: tr(c_{n-1-2k}) is a polynomial in (sum_P f^0, sum_P f^1, ..., sum_P f^{2k}).

**(d)** The moment hierarchy:
- sum_P f^j for j=0,1 is universal (depends only on n).
- sum_P f^j for j=2,3 depends only on t_3 (the 3-cycle count).
- sum_P f^j for j=4 depends on MORE than t_3 (verified at n=5,7).

Therefore:
- **k=0:** tr(c_{n-1}) = (n-1)! (universal)
- **k=1:** tr(c_{n-3}) = 2(n-2)!(t_3 - C(n,3)/4) (depends only on t_3)
- **k=2:** tr(c_{n-5}) depends on t_3 AND sum_P f^4 (the 4th moment of forward-arc count)

---

## Proof of (a): the trace-permutation identity

By THM-053 (diagonal signed position formula): M[v,v](r) = sum_P (-1)^{pos(v,P)} prod_{edges} (r + s_{ij}).

Taking trace at odd n: tr(M(r)) = sum_P [sum_v (-1)^{pos(v,P)}] prod(r+s) = sum_P prod(r+s), since the alternating sum equals 1 at odd n.

Expanding the product: prod_{i=0}^{n-2} (r + s_i) = sum_{k=0}^{n-1} e_k(s_P) r^{n-1-k}.

So tr(c_{n-1-2k}) = sum_P e_{2k}(s_P). QED.

---

## Proof of (b): Newton's identity simplification

Since s_i in {-1/2, +1/2}:
- p_{2l} := sum s_i^{2l} = (n-1)/4^l (constant, independent of P)
- p_{2l+1} := sum s_i^{2l+1} = p_1/4^l where p_1 = f - (n-1)/2

By Newton's identities e_k is determined by p_1, ..., p_k. Since all power sums reduce to polynomials in p_1, each e_k is a polynomial in p_1 (= f - (n-1)/2) of degree k. QED.

---

## Proof of (d): moment hierarchy

**j=0,1 universal:** sum_P 1 = n! and sum_P f = C(n,2)*(n-1)! (each arc appears at (n-1)! permutations times C(n,2) possible positions... actually by symmetry of arc selection).

**j=2 depends on t_3:** This is proved in THM-054 / algebraic_proof_cn3.py. The key: sum_P f^2 = sum_P f_(1) + sum_P f_(2), and f_(2) decomposes into consecutive (involves score squares -> t_3) and non-consecutive (universal) parts.

**j=3 depends only on t_3:** Verified computationally at n=5 (all iso classes with same t_3 have same sum_P f^3).

**j=4 goes beyond t_3:** At n=5 with t_3=4 (score seq (1,2,2,2,3)): sum_P f^4 takes values {6084, 6132, 6180} across non-isomorphic tournaments.

---

## The 4th moment decomposition

sum_P f_(4) = sum_P sum_{i<j<k<l} T_i T_j T_k T_l decomposes by position pattern at n=7:

| Pattern | Subsets | #Vertices | Count | Interpretation |
|---------|---------|-----------|-------|----------------|
| 4 consecutive | {0,1,2,3} etc. | 5 | 3 | (n-5)! * sum_{5-subsets} H(T[S]) |
| 3+1 | {0,1,2,4} etc. | 6 | 6 | 6-vertex path-arc correlations |
| 2+2 | {0,1,3,4} etc. | 6 | 3 | Disjoint 2-path pair correlations |
| 2+1+1 | {0,1,3,5} etc. | 7 | 3 | Full tournament invariant |

The "2+1+1" pattern involves all n=7 vertices and captures tournament structure that NO local invariant (cycle counts, degree sequence) can encode.

---

## Cross-scale universality boundary

| n | c_{n-1} | c_{n-3} | c_{n-5} | c_{n-7} |
|---|---------|---------|---------|---------|
| 5 | (n-1)! (univ) | f(t_3) | f(t_3, H) | - |
| 7 | (n-1)! (univ) | f(t_3) | f(t_3, sum_P f^4) | f(t_3, f^4, f^6) |
| 9 | (n-1)! (univ) | f(t_3) | f(t_3, sum_P f^4) | f(t_3, f^4, f^6) |

At n=5: sum_P f_(4) = H (since C(4,4)=1 position subset = all edges = Hamiltonian path count).

---

## Singleton Cancellation Lemma

**Lemma.** For centered edge deviations s_i = A[p_i, p_{i+1}] - 1/2, if a position subset has an isolated position (not adjacent to any other selected position), then sum_P prod_{selected positions} s_i = 0.

**Proof sketch:** The isolated position involves two vertices not constrained by other selections. Summing over all placements of these two vertices gives sum (A[a,b] - 1/2) = 0.

Verified computationally at n=7: (3,1) and (2,1,1) patterns both give exactly 0.

Consequence: At n=7, only 6 of 15 position 4-subsets contribute to tr(c_2).

---

## Exact Formula at n=7 (VERIFIED)

**Theorem.** At n=7:

tr(c_2) = 24*bc - 60*t_3 + 12*t_5 + 231

where bc = sum over 6-vertex subsets S of #{unordered partitions (T,T') of S into two triples that are both directed 3-cycles}.

**Derivation:**
- (4,) contribution: 126 - 36*t_3 + 12*t_5 (uses OCF at n=5 recursively)
- (2,2) contribution: 24*bc - 24*t_3 + 105 (complementary-triple structure)
- All other patterns: 0 (singleton cancellation)
- Total: 24*bc - 60*t_3 + 12*t_5 + 231

**Consequence:** tr(c_0) = H - 6*bc - 3*t_5 + 249/4

**Complete coefficient table at n=7:**

| Coefficient | Formula | Depends on |
|-------------|---------|------------|
| tr(c_6) | 720 | universal |
| tr(c_4) | 240*t_3 - 2100 | t_3 |
| tr(c_2) | 24*bc - 60*t_3 + 12*t_5 + 231 | t_3, t_5, bc |
| tr(c_0) | H - 6*bc - 3*t_5 + 249/4 | H, t_5, bc |

Verified: max error = 0.000000 over 30 random tournaments.
Verified: c_0 + c_2/4 + c_4/16 + c_6/64 = H (exact).

---

## Scripts

- `04-computation/algebraic_proof_cn3.py` (k=1 algebraic proof)
- `04-computation/coefficient_hierarchy_proof.py` (Newton identity + moment hierarchy)
- `04-computation/trc2_exact_formula.py` (complete n=7 formula verification)
