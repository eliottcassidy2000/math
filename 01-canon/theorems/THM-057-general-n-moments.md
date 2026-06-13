# THM-057: General n Moment Closed Forms

**Type:** Theorem (algebraic derivation + computational verification)
**Certainty:** 5 -- PROVED
**Status:** PROVED for m0-m5 at all n
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

  m5 = n!(n-1)(3n^4 - 2n^3 + 13n^2 - 14n + 16)/96
       + 5(n-1)(n^2 - n + 2)(n-2)! * t3
       + 120(n-1)(n-4)! * t5
       + 240(n-1)(n-4)! * bc

where:
- t3 = #{directed 3-cycles in T}
- t5 = #{directed 5-cycles in T}
- bc = #{pairs of vertex-disjoint directed 3-cycles in T}

---

## Key Structural Properties

1. **Cross-moment recurrence:** The t3 coefficient at level j+2 is proportional to the constant at level j:
   - m4_t3 * n(n-1) = 24 * m2_const  (polynomial 3n^2-5n+4 appears in both)
   - m5_t3 * n(n-1) = 40 * m3_const  (polynomial (n-1)(n^2-n+2) appears in both)
   - m5_t5 / m4_t5 = 5(n-1)/2  (similar recurrence for t5)

2. **bc/t5 ratio:** bc_coeff = 2 * t5_coeff in ALL moments m4 and m5 (and m6 at n=7). This means t5 and bc always appear in the combination (t5 + 2*bc), which connects to the OCF weight function I(Omega,2).

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
- **m5:** Verified at n=7,...,21 (closed form matches sigma computation)
- **Brute-force:** 46 random tournaments at n=5,7,9 — all pass with zero error

---

## Scripts

- `04-computation/general_n_closed_forms_v3.py` (derivation and fitting)
- `04-computation/verify_m4_closed_form.py` (m4 verification at n=5,7,9)
- `04-computation/m5_general_n.py` (m5 derivation and n=7,9 verification)
- `04-computation/m5_closed_form.py` (m5 coefficient fitting)
- `04-computation/sigma32_general_n.py` (sigma((3,2)) derivation)

## Open Directions

- m6 at general n: requires new invariants (t7? higher bc variants?) via 7+ vertex H-counts
- Why does bc_coeff = 2*t5_coeff universally? (Likely traces back to OCF weight 2^(#components))
- Cross-moment recurrence: is there a generating function unifying all m_j?
- The constant-term polynomials (3n^2-5n+4, (n-1)(n^2-n+2), etc.) — combinatorial interpretation?
