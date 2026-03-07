# THM-057: General n Moment Closed Forms

**Type:** Theorem (algebraic derivation + computational verification)
**Certainty:** 5 -- PROVED
**Status:** PROVED for m0-m4 at all odd n
**Added by:** opus-2026-03-06-S28
**Tags:** #moments #general-n #closed-form #sigma-patterns

---

## Statement

**Theorem.** For a tournament T on n vertices, let f_P = #{forward arcs in permutation P}. The moments m_j = sum_{P in S_n} f_P^j satisfy:

### Universal moments (independent of T)

  m0 = n!
  m1 = n!(n-1)/2

### Moments depending only on t3

  m2 = n!(3n^2 - 5n + 4)/12  +  4(n-2)! * t3
  m3 = n!(n-1)(n^2 - n + 2)/8  +  6(n-1)! * t3

### Moments depending on (t3, t5, bc)

  m4 = n!(15n^4 - 30n^3 + 65n^2 - 82n + 48)/240
       + 2(n-2)!(3n^2 - 5n + 4) * t3
       + 48(n-4)! * t5
       + 96(n-4)! * bc

where:
- t3 = #{directed 3-cycles in T}
- t5 = #{directed 5-cycles in T}
- bc = #{pairs of vertex-disjoint directed 3-cycles in T}

---

## Key Structural Properties

1. **The polynomial 3n^2 - 5n + 4** appears in both m2 (as the constant term normalizer) and m4 (as the t3 coefficient normalizer). Specifically:
   - m2_const = n! * (3n^2-5n+4) / 12
   - m4_t3 = 2(n-2)! * (3n^2-5n+4)
   - Cross-relation: m4_t3 = 24 * m2_const / (n(n-1))

2. **bc/t5 ratio in m4:** bc_coeff = 2 * t5_coeff = 96(n-4)!.

3. **Hierarchy:** New invariants enter at even moments:
   - j=0,1: universal
   - j=2,3: t3 only
   - j=4,5: (t3, t5, bc)
   - j=6,...,n-1: progressively more invariants from OCF hierarchy

4. **At n=5:** bc is structurally zero (needs >= 6 vertices). The formula m4 = 3444 + 648*t3 + 48*t5 agrees with prior computation.

---

## Proof Architecture

### Step 1: Stirling Decomposition

  m_j = sum_{k=0}^{j} S(j,k) * k! * SIGMA_k

where SIGMA_k = sum over size-k position subsets of sigma(S).

### Step 2: General Sigma Formula (THM-056)

For position subset S with connected components of sizes (c_1,...,c_m):

  sigma(S) = (n - sum c_i - m)! * sum_{ordered partitions} prod_j H(T[G_j])

### Step 3: Sub-tournament H-counts

| Component size | Vertex group | H formula |
|---------------|-------------|-----------|
| 1 | 2 vertices | H = 1 (universal) |
| 2 | 3 vertices | H = 1 + 2*cyc (depends on t3) |
| 3 | 4 vertices | H = 1 + 2*c3 (depends on t3) |
| 4 | 5 vertices | H = 1 + 2*(c3+c5) (depends on t3, t5) |
| (2,2) | 3+3 vertices | Product involves bc |

### Step 4: Aggregation

Combining pattern counts with sigma formulas and Stirling weights gives the closed forms.

---

## Computational Verification

- **m2, m3:** Verified at n=3,5,7,9,11,13,15 (closed form matches sigma computation exactly)
- **m4:** Verified at n=5,...,21 (closed form matches sigma computation)
- **Brute-force:** 26 random tournaments at n=5,7,9 — all pass with zero error

---

## Scripts

- `04-computation/general_n_closed_forms_v3.py` (derivation and fitting)
- `04-computation/verify_m4_closed_form.py` (comprehensive verification)
- `04-computation/general_n_moments.py` (initial exploration)

---

## Open Directions

- m5: requires sigma((3,2)) formula at general n, which involves cross-correlations of 3-cycles
- Complete hierarchy at general n: when does t7 enter? (Likely at m6, via 8-vertex sub-tournament H-counts)
- Can the polynomial 3n^2-5n+4 be given a combinatorial interpretation?
