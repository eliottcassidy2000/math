# THM-053: Diagonal Signed Position Formula

**Type:** Theorem (PROVED)
**Certainty:** 5 -- PROVED with clean bijective argument
**Status:** PROVED
**Added by:** opus-2026-03-06-S11b (continued³)
**Tags:** #transfer-matrix #diagonal #position #bijection

---

## Statement

**Theorem.** For any tournament T on n vertices, the diagonal entries of the transfer matrix satisfy:

M[a,a] = sum_{P} (-1)^{pos(a,P)}

where the sum is over all Hamiltonian paths P of T, and pos(a,P) is the 0-indexed position of vertex a in path P.

Equivalently: M[a,a] = #(paths with a at even position) - #(paths with a at odd position).

---

## Proof

The IE formula gives:

M[a,a] = sum_{S subset U} (-1)^|S| * E_a(S+{a}) * B_a(R+{a})

where U = V\{a}, R = U\S, E_a(W) = #Ham paths in T[W] ending at a, B_a(W) = #Ham paths in T[W] starting at a.

**Claim:** The map from triples (S, pi_S, pi_R) to Hamiltonian paths P is a bijection, where:
- pi_S is a Hamiltonian path of T[S+{a}] ending at a
- pi_R is a Hamiltonian path of T[R+{a}] starting at a

**Bijection construction:** Given (S, pi_S, pi_R), define P = pi_S followed by pi_R (concatenation, identifying the shared endpoint a). This gives P = (first(pi_S), ..., a, ..., last(pi_R)), a Hamiltonian path of T visiting all n vertices exactly once.

**Well-defined:** P visits S ∪ {a} ∪ R = V, and S ∩ R = empty (since S, R partition U), so all vertices are distinct. The path is valid because pi_S is a valid path in T[S+{a}] and pi_R in T[R+{a}].

**Injective:** Different triples give different paths because P uniquely determines (S, pi_S, pi_R): S = {vertices before a in P}, pi_S = P restricted to positions 0,...,pos(a,P), pi_R = P restricted to positions pos(a,P),...,n-1.

**Surjective:** Every Hamiltonian path P determines a unique triple by the same decomposition.

Under this bijection, |S| = pos(a,P), so:

M[a,a] = sum_{(S, pi_S, pi_R)} (-1)^|S| = sum_P (-1)^{pos(a,P)}.  QED.

---

## Consequences

1. **Trace formula reproof:** sum_a M[a,a] = sum_P sum_a (-1)^{pos(a,P)} = sum_P (sum_{k=0}^{n-1} (-1)^k) = sum_P 1 at odd n (since alternating sum = 1), giving tr(M) = H(T). At even n: sum = 0, giving tr(M) = 0. This reproves THM-027.

2. **VT diagonal:** For vertex-transitive T at odd n: every vertex appears equally often at each position, so M[a,a] = H/n * (alternating sum over position frequencies) = H/n.

3. **Defect vertex:** A vertex a with M[a,a] significantly different from H/n has an unequal distribution across positions — it appears disproportionately at even or odd positions.

---

## Verification

| n | Method | Result |
|---|--------|--------|
| 5 | All 12 iso classes | All match exactly |
| 7 | Regular (H=171, 175, 189) | All match exactly |
| 9 | 5 random regular | All match exactly |

**Scripts:** `04-computation/diagonal_signed_position_theorem.py`
