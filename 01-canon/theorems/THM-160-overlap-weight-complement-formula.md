# THM-160: Overlap Weight Complement Formula

**Status:** PROVED (trivial) + VERIFIED (p=7,11)
**Session:** kind-pasteur-2026-03-13-S60

## Statement

For any tournament T on n vertices, let Omega(T) be the conflict graph of directed odd cycles. For any directed odd cycle C on vertex set V, its overlap weight (degree in Omega) is:

    w(C) = N_total - 1 - C_odd(T[comp(V)])

where:
- N_total = sum_{k odd} c_k(T) = total directed odd cycles in T
- C_odd(T[W]) = sum_{k odd, k <= |W|} c_k(T[W]) = total odd cycles in induced subtournament on W
- comp(V) = [n] \ V

## Proof

w(C) counts cycles C' != C with V(C') ∩ V(C) ≠ ∅.

By complementary counting:
  w(C) = (N_total - 1) - |{C' : V(C') ⊆ comp(V)}|
       = N_total - 1 - C_odd(T[comp(V)])    □

## Consequences

### 1. Trivial constancy when 2k > n
When k > n/2, the complement comp(V) has fewer than k vertices and cannot support any k-cycles. But it may still support shorter odd cycles. When |comp(V)| < 3 (i.e., k ≥ n-2), no odd cycles fit at all, so w(C) = N_total - 1 is constant.

### 2. Circulant tournaments: C_odd(comp) is orbit-invariant
For a circulant tournament on Z_p (p prime), C_odd(comp(V+t)) = C_odd(comp(V)) for all t. So w(C) depends only on the Z_p-orbit of V.

### 3. Paley p=11 overlap weight classification
- k=3: C_odd(8-vertex comp) = 172 for ALL 55 orbits → w = 20996 CONSTANT
- k=5: C_odd(6-vertex comp) ∈ {7, 10, 12, 13, 16} → w ∈ {21152,...,21161}, 5 values
- k=7: C_odd(4-vertex comp) ∈ {0, 1, 2} → w ∈ {21166, 21167, 21168}, 3 values
- k=9,11: comp too small → w = 21168 CONSTANT

### 4. Complete constancy at p=7
For both Paley and Interval at p=7, C_odd(comp) is constant for ALL cycle lengths → uniform overlap weight within each length. This implies the conflict graph has a very regular structure.

## Connection to disj(k1,k2)

The formula immediately gives:
    sum_{V active for k} n(V) * C_odd(comp(V))
    = sum_{k'} sum_{V1,V2 disjoint} n_k(V1) * n_{k'}(V2)
    = sum_{k'} disj(k, k')

So the overlap weight decomposes the disjoint pair counts:
    sum_C w(C) = N_total * c_k - c_k - sum_{k'} disj(k, k')_ordered

## Verification

- overlap_weight_formula.out: p=11 Paley, w formula verified for all 5 types
- master_overlap_formula.out: p=7,11 Paley and Interval, all cycle types verified
- cross_length_overlap.out: full decomposition by target cycle length
