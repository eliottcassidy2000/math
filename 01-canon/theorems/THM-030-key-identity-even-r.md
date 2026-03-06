# THM-030: Key Identity and Even r-Powers of the Transfer Matrix

**Status:** PROVED (inductive proof, computationally verified m=2..6)
**Added:** opus-2026-03-06-S25
**Proved by:** opus-2026-03-06-S25

## Statement

**Theorem (Key Identity).** Let W be a vertex set with |W| = m >= 2 in a c-tournament with edge weights t(i,j) = r + s_{ij}, s_{ij} = -s_{ji}. For any b in W,

    B_b(W) + (-1)^m E_b(W) = 2r * col_sum_W(b)

where:
- B_b(W) = sum of Hamiltonian path weights starting at b through W
- E_b(W) = sum of Hamiltonian path weights ending at b through W
- col_sum_W(b) = sum_{a != b} M_W[a,b], the column sum of the transfer matrix

## Corollaries

**Corollary 1 (Sigma Identity).** Summing the Key Identity over all b in W:

    T(W) * [1 + (-1)^m] = 2r * Sigma(W)

where T(W) = total Hamiltonian path weight, Sigma(W) = sum_{a!=b} M[a,b].
- Even m: T(W) = r * Sigma(W)
- Odd m: Sigma(W) = 0

**Corollary 2 (Even r-powers).** Every entry M[a,b] of the transfer matrix contains only even powers of r.

**Corollary 3 (Symmetry).** M[a,b] = M[b,a] for all a, b.

## Proof

By strong induction on m = |W|.

### Base case (m = 2)

W = {a, b}. Then M[a,b] = 1 (trivial paths), col_sum(b) = 1.
- B_b(W) = t(b,a) = r - s_{ab}
- E_b(W) = t(a,b) = r + s_{ab}
- LHS = (r - s_{ab}) + (r + s_{ab}) = 2r
- RHS = 2r * 1 = 2r. CHECK.

### Inductive step (m >= 3)

Assume the Key Identity holds for all vertex sets of size < m.

**Step 1: Three recurrences** (algebraic identities from definitions, W' = W\{b}):

(R1) B_b(W) = r * T(W') + sum_{w in W'} s_{bw} * B_w(W')
     [decompose by first edge: t(b,w) = r + s_{bw}]

(R2) E_b(W) = r * T(W') - sum_{w in W'} s_{bw} * E_w(W')
     [decompose by last edge: t(w,b) = r + s_{wb} = r - s_{bw}]

(R3) col_sum_W(b) = (-1)^{m-2} * T(W') + r * Sigma(W') + sum_{w in W'} s_{bw} * cs_{W'}(w)
     [from the column recurrence of M, summing over a]

**Step 2: Compute LHS = B_b + (-1)^m E_b.**

Using (R1) and (R2):

    LHS = r*T(W') + sum_w s_{bw} B_w(W') + (-1)^m [r*T(W') - sum_w s_{bw} E_w(W')]
        = r*T(W')[1 + (-1)^m] + sum_w s_{bw} [B_w(W') + (-1)^{m-1} E_w(W')]

Note: (-1)^m * (-1) = (-1)^{m+1} = (-1)^{m-1}.

By the **inductive hypothesis** (Key Identity at size m-1 = |W'|):

    B_w(W') + (-1)^{m-1} E_w(W') = 2r * cs_{W'}(w)

Therefore:

    LHS = r*T(W')[1 + (-1)^m] + 2r * sum_w s_{bw} * cs_{W'}(w)

**Step 3: Compute RHS = 2r * col_sum_W(b).**

Using (R3):

    RHS = 2r*(-1)^{m-2}*T(W') + 2r^2*Sigma(W') + 2r * sum_w s_{bw} * cs_{W'}(w)

**Step 4: Equate LHS = RHS.**

The terms 2r * sum_w s_{bw} * cs_{W'}(w) cancel from both sides. What remains:

    r*T(W')[1 + (-1)^m] = 2r*(-1)^{m-2}*T(W') + 2r^2*Sigma(W')

Since (-1)^{m-2} = (-1)^m:

    r*T(W')[1 + (-1)^m] = 2r*(-1)^m*T(W') + 2r^2*Sigma(W')
    r*T(W')[1 - (-1)^m] = 2r^2*Sigma(W')      ... (*)

**Step 5: Verify (*) using the Sigma identity at size m-1.**

From the Key Identity at size m-1, summing over all vertices of W':

    T(W')[1 + (-1)^{m-1}] = 2r * Sigma(W')

**Case m even:** m-1 is odd, so 1 + (-1)^{m-1} = 0, giving Sigma(W') = 0.
Equation (*): LHS = r*T(W')[1 - 1] = 0, RHS = 0. CHECK.

**Case m odd:** m-1 is even, so 1 + (-1)^{m-1} = 2, giving T(W') = r*Sigma(W').
Equation (*): LHS = 2r*T(W'), RHS = 2r^2*Sigma(W') = 2r*T(W'). CHECK.

QED.

### Proof of Corollary 2 (Even r-powers)

By the path reversal identity: B_b(W; -r) = (-1)^{m-1} E_b(W; r).
Therefore: odd_r(B_b) = [B_b(r) - (-1)^{m-1} E_b(r)] / 2 = [B_b + (-1)^m E_b] / 2 = r * col_sum(b).

The row recurrence gives: M[a,b] = B_b(V\{a}) - sum_v t(v,a) M_{smaller}[v,b].
By induction, the smaller M entries have even r-powers.
So odd_r(M[a,b]) = odd_r(B_b) - r * col_sum(b) = 0.

### Proof of Corollary 3 (Symmetry)

The re-indexing identity (proved separately): M_T[b,a] = (-1)^{m-2} M_{T^op}[a,b].
T^op corresponds to s_{ij} -> -s_{ij} with r unchanged.
Even r-powers means every monomial r^k * (s-product) has k even.
Since total degree is m-2, the s-degree has parity (m-2) mod 2 = parity of m.
So M[a,b]|_{s->-s} = (-1)^{m-2} M[a,b], giving M[b,a] = M[a,b].

## Dependencies

- Transfer matrix definition (definitions.md)
- Path reversal identity: B_b(W; -r) = (-1)^{m-1} E_b(W; r)
- Re-indexing identity: M_T[b,a] = (-1)^{m-2} M_{T^op}[a,b] (THM-027 or standalone)
- Column recurrence of M (algebraic identity from definition)

## Computational verification

File: `04-computation/key_identity_complete_proof.py`
Verified: m = 2, 3, 4, 5, 6 (all vertices, all recurrences, all corollaries)
