# THM-030: Transfer Matrix Symmetry (Even r-Powers)

**Status:** PROVED
**Proved by:** kind-pasteur-2026-03-06-S25, building on opus-2026-03-06-S24 (Key Identity discovery + column recurrence)
**Verification:** Symbolic, all (a,b) pairs, m = 2,...,6
**Script:** `04-computation/complete_even_r_proof.py`

---

## Statement

**Theorem.** Let T be a c-tournament on vertex set W with edge weights t_{ij} = r + s_{ij}, where s_{ij} = -s_{ji} and r = c/2. Define the transfer matrix:

M_W[a,b] = sum_{S subset U} (-1)^|S| E_a(S + {a}) B_b(R + {b})

where U = W \ {a,b}, R = U \ S, E_a counts Hamiltonian paths ending at a, and B_b counts Hamiltonian paths starting at b.

Then M_W[a,b] has only even powers of r. Equivalently, M_W[a,b] = M_W[b,a].

---

## Proof

By strong induction on m = |W|, we prove the **Key Identity**:

> odd_r(B_b(W)) = r * col_sum_W(b)

where col_sum_W(b) = sum_{v != b} M_W[v,b].

### Base cases
- m = 1: B_b({b}) = 1, odd(1) = 0 = r * 0.
- m = 2: B_b({b,v}) = r + s_{bv}, odd = r, col_sum = M[v,b] = 1, r * 1 = r.

### Inductive step (m >= 3)

Assume the Key Identity for all vertex sets of size < m.

**Ingredient I: Decomposition of B_b.**
B_b(W) = sum_{v in W'} t(b,v) B_v(W') = r * T(W') + sum_v s_{bv} B_v(W')
where W' = W \ {b}, T(W') = total Hamiltonian path weight.

**Ingredient II: Column recurrence.**
M_W[a,b] = (-1)^{m-2} E_a(W') + sum_{w in W' \ {a}} t(b,w) M_{W'}[a,w]

Summing over a in W':
col_sum_W(b) = (-1)^{m-2} T(W') + r * Sigma_{W'} + sum_w s_{bw} cs_{W'}(w)
where Sigma_{W'} = sum_{a != w in W'} M_{W'}[a,w].

**Ingredient III: Inductive evaluation of Sigma.**
By induction at size m-1: odd(B_w(W')) = r * cs_{W'}(w) for each w.
Summing over w: odd(T(W')) = r * Sigma_{W'}.

Since T(W') has definite r-parity (-1)^{m-2}:
(-1)^{m-2} T(W') + odd(T(W')) = even_r(T(W'))

Therefore: col_sum_W(b) = even_r(T(W')) + sum_w s_{bw} cs_{W'}(w)  ...(*)

**Ingredient IV: Inductive decomposition of B_v.**
By induction, M_{W'} has even r-powers. Therefore:
- cs_{W'}(v) has even r-powers
- B_v(W') = even_r(B_v) + r * cs_{W'}(v)

**Proof of Key Identity:**
odd(B_b(W)) = odd(r * T + sum_v s_{bv} B_v)
= r * even_r(T) + odd(sum_v s_{bv} B_v)

Substituting B_v = even_r(B_v) + r * cs_v:
sum s_{bv} B_v = sum s_{bv} even_r(B_v) + r * sum s_{bv} cs_v

- odd(sum s * even_r(B_v)) = 0  (product of even-r terms)
- odd(r * sum s * cs) = r * sum s * cs  (sum s*cs is already even in r)

Therefore: odd(B_b) = r * even_r(T) + r * sum s_{bv} cs_v = r * col_sum_W(b) by (*).

### From Key Identity to even r-powers

Row recurrence: M_W[a,b] = B_b(W \ {a}) - sum_u t(u,a) M_{W\{a}}[u,b]

By induction, M_{W\{a}} has even r-powers. So:
odd(M_W[a,b]) = odd(B_b(W \ {a})) - r * col_sum_{W\{a}}(b) = 0

by the Key Identity at size m-1. QED.

---

## Dependencies
- Column recurrence for M (algebraic identity, proved in column_recurrence_proof.py)
- T(W) has definite r-parity (-1)^{|W|-1} (proved via path reversal)
- M[a,b](-r) = M[b,a](r) (proved algebraically in kind-pasteur-S23b)

## Consequences
- M_W[a,b] = M_W[b,a] for all c-tournaments (transfer matrix is symmetric)
- Combined with THM-027 (trace formula): tr(M) = H(T) for odd n, 0 for even n
- Off-diagonal sum: sum_{a!=b} M[a,b] = 0 (odd n), 2*H(T) (even n)
