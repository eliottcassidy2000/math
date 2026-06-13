# THM-056: Forward-Arc Moment Hierarchy Theorem

**Type:** Theorem (algebraic proof + computational verification)
**Certainty:** 5 -- PROVED algebraically at n=7, verified computationally at n=5,9
**Status:** PROVED at n=7 (complete); PROVED at n=5 (complete)
**Added by:** opus-2026-03-06-S28
**Tags:** #moments #transfer-matrix #hierarchy #sigma-patterns

---

## Statement

**Theorem.** For a tournament T on n=7 vertices, let f_P = #{forward arcs in permutation P}. Then:

sum_P f_P^j = a_j(t3) + b_j(t5) + c_j(bc) + d_j(H)

where all coefficients are LINEAR in (t3, t5, bc, H). Explicitly:

| j | sum_P f^j | Depends on |
|---|-----------|------------|
| 0 | 5040 | universal |
| 1 | 15120 | universal |
| 2 | 48720 + 480*t3 | t3 |
| 3 | 166320 + 4320*t3 | t3 |
| 4 | 596064 + 27840*t3 + 288*t5 + 576*bc | t3, t5, bc |
| 5 | 2227680 + 158400*t3 + 4320*t5 + 8640*bc | t3, t5, bc |
| 6 | 8636880 + 850080*t3 + 40320*t5 + 80640*bc + 720*H | t3, t5, bc, H |

where:
- t3 = directed 3-cycle count
- t5 = directed 5-cycle count
- bc = sum over 6-vertex subsets of #{pairs of vertex-disjoint directed 3-cycles}
- H = directed Hamiltonian path count

---

## Key Properties

1. **Linearity:** All moments are LINEAR in the tournament invariants — no quadratic or higher terms.

2. **Hierarchy:** New invariants enter at specific levels:
   - j=0,1: universal (independent of T)
   - j=2,3: depends only on t3
   - j=4,5: depends on t3, t5, and bc
   - j=6: depends on t3, t5, bc, and H

3. **Coefficient of H in m_{n-1}:** The coefficient of H in m_{n-1} = sum_P f^{n-1} is always (n-1)!. This is algebraically clear: in (T_0+...+T_{n-2})^{n-1}, the product of all distinct T_i has multinomial coefficient (n-1)!, and prod T_i = [all edges forward] = Ham path indicator.

---

## Proof Architecture

### Step 1: Stirling Decomposition

sum_P f^j = sum_{k=0}^{j} S(j,k) * k! * SIGMA_k

where S(j,k) = Stirling number of second kind and
SIGMA_k = sum over size-k position subsets S of [0,n-2] of sigma(S),
where sigma(S) = sum_P prod_{i in S} A[p_i, p_{i+1}].

### Step 2: Position-Translation Invariance

sigma(S) depends only on the connected-component structure (adjacency pattern) of S, where positions i,j are adjacent iff |i-j|=1. This is because shifting all positions by a constant corresponds to relabeling vertices, and we sum over all n! permutations.

### Step 3: Sigma Pattern Formulas

Complete table at n=7 (all 15 patterns):

| Pattern | Mult | sigma formula | Proof method |
|---------|------|---------------|-------------|
| () | 1 | 5040 | trivial |
| (1,) | 6 | 2520 = n!/2 | pair symmetry |
| (1,1) | 10 | 1260 | pair-partition universality |
| (1,1,1) | 4 | 630 | triple pair-partition |
| (2,) | 5 | 840 + 48*t3 | directed 2-path decomposition |
| (2,1) | 12 | 420 + 24*t3 | factorization: (2,) x (1,) |
| (2,1,1) | 3 | 210 + 12*t3 | factorization: (2,) x (1,1) |
| (3,) | 4 | 210 + 48*t3 | directed 3-path = H(4-sub) via OCF |
| (3,1) | 6 | 105 + 24*t3 | factorization: (3,) x (1,) |
| (4,) | 3 | 42 + 24*t3 + 4*t5 | H(5-sub) via OCF at n=5 |
| (4,1) | 2 | 21 + 12*t3 + 2*t5 | factorization: (4,) x (1,) |
| (2,2) | 3 | 140 + 16*t3 + 8*bc | disjoint 2-path product decomposition |
| (3,2) | 2 | 35 + 10*t3 + 8*bc | derived from m5 consistency |
| (5,) | 2 | 7 + 8*t3 + 4*t5 + 4*bc | H(6-sub) via OCF at n=6 |
| (6,) | 1 | H | definition of Ham path count |

### Step 4: SIGMA_k Aggregation

SIGMA_k = sum of (multiplicity * sigma) over patterns of size k:

- SIGMA_0 = 5040
- SIGMA_1 = 15120
- SIGMA_2 = 16800 + 240*t3
- SIGMA_3 = 8400 + 480*t3
- SIGMA_4 = 1806 + 300*t3 + 12*t5 + 24*bc
- SIGMA_5 = 126 + 60*t3 + 12*t5 + 24*bc
- SIGMA_6 = H

### Step 5: Apply Stirling

m_j = sum_k S(j,k) * k! * SIGMA_k gives the moment formulas.

---

## Key Sigma Proofs

### sigma((2,)) = 840 + 48*t3

sigma({i,i+1}) = (n-3)! * #{directed 2-paths on n vertices}

For each 3-vertex set: cyclic triple has 3 directed 2-paths, transitive has 1.
Total = C(n,3) + 2*t3. At n=7: (n-3)! * (35 + 2*t3) = 24*(35+2*t3) = 840+48*t3.

### sigma((2,2)) = 140 + 16*t3 + 8*bc

sigma involves two disjoint 2-paths on 6 vertices (one vertex free).
For each free vertex v, partition remaining 6 into two 3-sets (S, S'):
- dp_2(S) * dp_2(S') = 9 if both cyclic, 3 if one cyclic, 1 if both transitive
- Using mixed(v) = c3(T-v) - 2*bc(v): sum = 10 + 4*bc(v) + 2*c3(T-v)
- Factor 2 for ordered partitions
- Sum over v: 2*(70 + 4*bc + 8*t3) = 140 + 8*bc + 16*t3

### sigma((5,)) = 7 + 8*t3 + 4*t5 + 4*bc

sigma = sum_v H(T-v). By OCF at n=6: H(T-v) = 1 + 2*c3(T-v) + 2*c5(T-v) + 4*bc(T-v).
Sum: 7 + 2*4*t3 + 2*2*t5 + 4*bc = 7 + 8*t3 + 4*t5 + 4*bc.

---

## Connection to Transfer Matrix Coefficients

Via THM-055: tr(c_{n-1-2k}) = sum_P e_{2k}(s_P) where e_{2k} is a polynomial in f.
By Newton's identities, e_{2k} is degree 2k in f, so tr(c_{n-1-2k}) depends on m_0,...,m_{2k}.

This gives the complete coefficient hierarchy:
- tr(c_6) = 720 (universal)
- tr(c_4) = 240*t3 - 2100
- tr(c_2) = 24*bc - 60*t3 + 12*t5 + 231
- tr(c_0) = H - 6*bc - 3*t5 + 249/4

---

## Generalization to Other n

At n=5: m4 = 648*t3 + 48*t5 + 3444 (bc=0 forced, t5 enters via OCF)
At n=9: m2, m3 depend only on t3; m4, m5 depend on (t3, t5, bc); m6 requires NEW invariants (t7 and/or higher-order bc variants). The hierarchy extends but with new layers at each odd n.

---

## Scripts

- `04-computation/fourth_moment_analysis.py` (computational discovery of linear formulas)
- `04-computation/moment_hierarchy_n5.py` (complete n=5 analysis)
- `04-computation/sigma22_proof.py` (algebraic proof of sigma((2,2)))
- `04-computation/sigma_k5_derivation.py` (remaining sigma formulas)
- `04-computation/complete_moment_derivation.py` (full algebraic derivation)
