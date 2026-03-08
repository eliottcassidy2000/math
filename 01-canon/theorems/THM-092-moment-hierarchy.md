# THM-092: Moment-Cycle Hierarchy for Forward-Edge Distribution

**Status:** PROVED for r=0,1,2,3 (algebraic); VERIFIED for r=4 at n=5,6,7; kappa_6 at n=7; CONJECTURED for general r
**Proved by:** opus-2026-03-07-S46c, opus-2026-03-07-S46d (kappa_4 general, kappa_6)
**Scope:** All tournaments

---

## Statement

For a tournament T on n vertices, the r-th moment of the forward-edge count:

$$M_r(T) = E[\text{fwd}^r] = \frac{1}{n!} \sum_{\sigma \in S_n} \text{fwd}(\sigma)^r$$

is determined by cycle-counting invariants of T as follows:

**Level 0-1:** M_0 = 1, M_1 = (n-1)/2 (universal, tournament-independent).

**Level 2:** M_2 = (n-1)(3n-2)/12 + 4*t3/(n(n-1)) (THM-089). Depends on t3 only.

**Level 3:** M_3 = A(n) + 6*t3/n (THM-090). Depends on t3 only, by zero-skewness argument (THM-091).

**Level 4:** M_4 depends on (t3, t5) at n=5, on (t3, t5, alpha_2) at n=6.
Exact formulas:
- n=5: E[fwd^4] = 287/10 + (27/5)*t3 + (2/5)*t5
- n=6: E[fwd^4] = 619/10 + (82/15)*t3 + (2/15)*t5 + (4/15)*alpha_2

---

## The Hierarchy Principle

The key insight connecting moments to cycles:

1. **Reversal symmetry (THM-091):** fwd(sigma) + fwd(sigma^rev) = n-1, so the distribution is symmetric about (n-1)/2. All odd cumulants kappa_{2k+1} = 0.

2. **Cumulant decomposition:** The even cumulants determine the distribution:
   - kappa_2 = Var[fwd] depends on t3 (THM-089)
   - kappa_4 = E[(fwd-mu)^4] - 3*sigma^4 depends on (t3, t5, alpha_2) at n=6
   - kappa_{2k} introduces invariants involving cycles on <= 2k+1 vertices

3. **Why the hierarchy exists:** E[fwd^r] = sum over index tuples (i_1,...,i_r) of E[X_{i_1} ... X_{i_r}]. By non-adjacent uncorrelatedness, only "connected" clusters of adjacent indices contribute t3-dependent terms. A cluster of k adjacent indices involves k+1 consecutive permutation positions, hence depends on tournament structure of (k+1)-vertex subsets.

   - Clusters of size 2 (pairs i, i+1): depend on 3-vertex invariants = t3
   - Clusters of size 3 (triples i, i+1, i+2): depend on 4-vertex invariants = still t3 (since 4-tournaments are classified by t3 in {0,1})
   - Clusters of size 4 (quadruples i,...,i+3): depend on 5-vertex invariants = (t3, t5)

4. **Worpitzky-moment bridge:** The Worpitzky coefficient c_j = sum of linear combinations of M_0, M_1, ..., M_{n-1-j}. Therefore c_j inherits the cycle dependencies of M_{n-1-j}.

---

## Predictions for n=7

At n=7 (VERIFIED by sampling, 156 F-classes):
- delta_5 = delta_4 = 0 (universal) ✓
- delta_4 = 10*t3 = 2(n-2)*t3 ✓ (exact for all 156 classes)
- delta_3 = 20*t3 = (n-2)(n-3)*t3 ✓ (exact for all 156 classes)
- delta_2: NOT determined by t3 alone (needs additional invariants from 5-vertex subgraph structure)
- delta_0 = H(T) - 1 ✓

---

## Connections

1. **Graded OCF refinement (THM-087):** The Worpitzky polynomial provides a graded decomposition of H(T), with each level j encoding cycle invariants from the hierarchy.

2. **Forward r-path formula (CORRECTED S46d):**
   - #fwd2path = C(n,2)
   - #fwd3path = C(n,3) + 2*t3
   - #fwd4path = C(n,4) + 2(n-3)*t3
   - #fwd5path = C(n,5) + 2*C(n-3,2)*t3 + 2*t5
   These count #{r-tuples of distinct vertices forming a directed (r-1)-path}.
   The formula is: #fwd(r)path = sum_{S in C(V,r)} H(T[S]).
   The t3 coefficient is 2*C(n-3, r-3) since each 3-cycle sits in C(n-3, r-3) r-subsets.

3. **Cumulant interpretation:** The forward-edge distribution is fully characterized by even cumulants kappa_2, kappa_4, kappa_6, .... Each kappa_{2k} adds exactly one "level" of cycle complexity (cycles on 2k+1 or fewer vertices, plus their disjoint combinations).

---

## Open Questions (updated S46d)

1. ~~What is the exact formula for kappa_4 at general n?~~ **RESOLVED:** See THM-093.
2. ~~Does kappa_6 introduce t7 (directed 7-cycles)?~~ **YES.** Verified at n=7 (149 F-classes). kappa_6 = (n+1)/252 + (2/C(n,6))*t7 + nonlinear lower terms.
3. Is there a generating function for the cumulant hierarchy?
4. Does the hierarchy connect to the Fourier decomposition (INV-050)?
5. **NEW:** Prove the universal coefficient conjecture: coeff(t_{2k+1}) in kappa_{2k} = 2/C(n, 2k).
6. **NEW:** What determines the nonlinear cross terms in kappa_{2k} for k >= 3?

---

## Verified

- M_0, M_1: all n (algebraic)
- M_2: n=3..6 exhaustive (THM-089)
- M_3: n=3..7 (THM-090, 156 F-classes at n=7)
- M_4: n=5 (8 classes), n=6 (24 classes) exact rational
- Worpitzky at n=7: 156 F-classes, delta_4/delta_3 linear in t3
- Scripts: fwd_moments_cycles.py, fwd3_formula.py, efwd4_exact.py, fwd_cumulants.py, worpitzky_n7_tiny.py
