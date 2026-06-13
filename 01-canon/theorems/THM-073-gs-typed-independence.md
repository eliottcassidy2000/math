# THM-073: GS Generator Expansion = Typed Independence Polynomial

**Type:** Theorem (consequence of Grinberg-Stanley)
**Certainty:** 5 -- PROVED
**Status:** PROVED
**Added by:** opus-2026-03-07-S37
**Tags:** #grinberg-stanley #typed-independence #symmetric-function #ocf

---

## Statement

For any tournament T on n vertices, the Rédei-Berge symmetric function U_T can be written as:

**U_T = Σ_{S} α_S · (2p₃)^{s₃} · (2p₅)^{s₅} · (2p₇)^{s₇} · ... · p₁^{n - 3s₃ - 5s₅ - 7s₇ - ...}**

where:
- S = (s₃, s₅, s₇, ...) ranges over all tuples of non-negative integers with 3s₃ + 5s₅ + 7s₇ + ... ≤ n
- α_S = number of collections of s₃ vertex-disjoint directed 3-cycles, s₅ vertex-disjoint directed 5-cycles, s₇ vertex-disjoint directed 7-cycles, ..., all pairwise vertex-disjoint, in T
- All α_S are non-negative integers

The principal specialization at m=1 gives:

**ps₁(U_T)(1) = Σ_{S} α_S · 2^{s₃ + s₅ + s₇ + ...} = I(Ω(T), 2) = H(T)**

which is exactly the OCF.

---

## Proof

The fact that U_T is a polynomial in p₁, 2p₃, 2p₅, ... with non-negative integer coefficients is Theorem 1.39 (combined with Lemma 6.5) of Grinberg & Stanley (arXiv:2307.05569).

The identification of the coefficients α_S with typed cycle collection counts follows from the permutation-theoretic definition of U_T: each valid permutation σ ∈ S_V(T, T^op) with all nontrivial cycles being odd directed T-cycles of lengths l₁, l₂, ..., l_k contributes to the monomial (2p_{l₁})(2p_{l₂})...(2p_{l_k}) · p₁^{n-l₁-...-l_k} with sign (-1)^φ = +1 (since each l_i is odd, l_i - 1 is even).

The specialization ps₁(p_k)(1) = 1 for all k means:
ps₁(U_T)(1) = Σ_S α_S · 2^{|S|} · 1^{n-Σl_i s_i} = Σ_S α_S · 2^{|S|} = I(Ω(T), 2).

This equals H(T) by Grinberg-Stanley's Corollary (restated as Corollary 20 in Irving-Omar arXiv:2412.10572).

---

## Verification

| Tournament | α_S coefficients | ps₁(1) | H(T) |
|---|---|---|---|
| Transitive T₅ | α_∅ = 1 | 1 | 1 |
| Cyclic C₅ | α_∅=1, α_{(1,0)}=5, α_{(0,1)}=2 | 1+10+4=15 | 15 |
| Paley T₇ | α_∅=1, α_{(1,0,0)}=14, α_{(2,0,0)}=7, α_{(0,1,0)}=42, α_{(0,0,1)}=24 | 1+28+28+84+48=189 | 189 |

---

## Significance

This theorem reveals that the GS symmetric function U_T **encodes the full typed independence polynomial** of Ω(T). The non-negative integer coefficients in the generator expansion directly count cycle collections by type. This is a more refined invariant than I(Ω(T), x) alone, as it tracks cycle lengths separately.

The specialization p₁→1, 2p_k→x for all odd k≥3 gives:
ps(U_T) = Σ_S α_S · x^{|S|} · 1^{n-...} = I(Ω(T), x)

recovering the full independence polynomial, not just its evaluation at x=2.

---

## Bug note

A bug was found in `04-computation/bags_of_sticks.py` (line 99): the opposite edge set was computed incorrectly as T itself instead of T^op. This caused the script to compute a "T-only" version of the Rédei-Berge sum, which is NOT U_T. The bug was fixed in this session. The correct U_T for tournaments has support only on all-odd partitions.

---

## References

- Grinberg & Stanley, arXiv:2307.05569, Theorem 1.39 + Lemma 6.5
- Irving & Omar, arXiv:2412.10572, Corollary 20

## Scripts

- `04-computation/bags_of_sticks.py` (with bug fix)
- `04-computation/omega_constraint.py`
