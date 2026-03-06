# THM-016: Hamiltonian Path Alternating Sum Identity

**Type:** Theorem
**Certainty:** 5 — PROVED
**Status:** PROVED (induction on m)
**Added by:** opus-2026-03-05-S4
**Tags:** #hamiltonian-path #alternating-sum #induction #ocf #claim-b-path

---

## Statement

For any weighted tournament on m vertices W with distinguished vertex v:

```
sum_{S⊆W\{v}} (-1)^|S| H(S) h_start(W\S, v) = (-1)^{m+1} h_end(W, v)
```

where:
- H(S) = total Hamiltonian path weight on S (all orderings; H(empty) = 1)
- h_start(R, v) = Hamiltonian path weight on R starting at v
- h_end(W, v) = Hamiltonian path weight on W ending at v

This holds as a **polynomial identity** in the arc variables (verified with real-valued arcs in [-2,3]).

---

## Proof

**By induction on m = |W|.**

**Base case m=1:** W = {v}. LHS = H(empty) * h_start({v}, v) = 1. RHS = (-1)^2 * 1 = 1.

**Inductive step:** Assume the identity for all sizes < m. Let W' = W\{v}.

**Step 1.** Separate the S = W' boundary term:

  Phi(v) = (-1)^{m-1} H(W') + sum_{S ⊊ W'} (-1)^|S| H(S) h_start(W\S, v)

**Step 2.** For S ⊊ W', the set R = W\S has |R| >= 2. Apply the first-step decomposition:

  h_start(R, v) = sum_{u in R\{v}} T(v,u) h_start(R\{v}, u)

**Step 3.** Exchange summation order:

  Phi(v) = (-1)^{m-1} H(W') + sum_{u in W'} T(v,u) [sum_{S ⊆ W'\{u}} (-1)^|S| H(S) h_start(W'\S, u)]

The inner sum is exactly Phi_{W'}(u) — this identity applied to the (m-1)-vertex set W' with distinguished vertex u.

**Step 4.** By the induction hypothesis:

  sum_{S ⊆ W'\{u}} (-1)^|S| H(S) h_start(W'\S, u) = (-1)^m h_end(W', u)

**Step 5.** Substitute and use T(v,u) = 1 - T(u,v):

  sum_u T(v,u) h_end(W', u) = sum_u h_end(W', u) - sum_u T(u,v) h_end(W', u)
                              = H(W') - h_end(W, v)

The second equality uses:
- sum_u h_end(W', u) = H(W')  (partition by ending vertex)
- sum_u T(u,v) h_end(W', u) = h_end(W, v)  (last-step decomposition of h_end(W,v))

**Step 6.** Combine:

  Phi(v) = (-1)^{m-1} H(W') + (-1)^m [H(W') - h_end(W, v)]
         = [(-1)^{m-1} + (-1)^m] H(W') + (-1)^{m+1} h_end(W, v)
         = (-1)^{m+1} h_end(W, v)

since (-1)^{m-1} + (-1)^m = 0 always. QED.

---

## Verification Record

| m | Trials | Max error | Status |
|---|--------|-----------|--------|
| 1 | 100 | 0 | PASS |
| 2 | 100 | 0 | PASS |
| 3 | 100 | 1.8e-15 | PASS |
| 4 | 100 | 6.2e-15 | PASS |
| 5 | 100 | 2.8e-14 | PASS |
| 6 | 30 | 2.3e-13 | PASS |
| 7 | 30 | 2.4e-12 | PASS |

Each proof step verified independently at every m. Code: `04-computation/q009_claim_b_proof.py`.

---

## Significance

This identity is the KEY LEMMA for proving B(Li,Rj) = B(Lj,Ri) (the even-odd split / signed adjacency identity) by induction. The proved chain:

1. **THM-016** (this theorem) for all m
2. → B(Li, Rj) = B(Lj, Ri) for all n (via induction on |W| = n-2, THM-017)
3. → Even-odd split for all n

**CAVEAT (MISTAKE-008):** Steps 1-3 do NOT prove OCF. The even-odd split is a necessary
condition for OCF, not sufficient. The chain B(Li,Rj)=B(Lj,Ri) → delta_H=delta_I was
previously claimed but is FALSE (see MISTAKE-008, signed-adjacency-identity.md).
OCF additionally requires proving that delta_H = delta_I (unsigned), which is proved
only at n<=8 (THM-015, THM-018).

The inductive proof of B(Li,Rj)=B(Lj,Ri) uses THM-016 for the d-independent part, and the induction hypothesis on sub-tournaments of size |W|-2 for the d-dependent part.
