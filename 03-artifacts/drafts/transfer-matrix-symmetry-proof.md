# Transfer Matrix Symmetry: Complete Proof

**Author:** opus-2026-03-06-S25
**Status:** COMPLETE PROOF (THM-030)

## Setup

Let $W$ be a set of $m \geq 2$ vertices in a c-tournament, where each edge has weight $t(i,j) = r + s_{ij}$ with $s_{ij} = -s_{ji}$ and $r = c/2$.

**Definitions:**
- $B_b(W)$: total weight of Hamiltonian paths through $W$ starting at $b$
- $E_b(W)$: total weight of Hamiltonian paths through $W$ ending at $b$
- $T(W) = \sum_b B_b(W) = \sum_b E_b(W)$: total Hamiltonian path weight
- $M_W[a,b] = \sum_{S \subseteq W \setminus \{a,b\}} (-1)^{|S|} E_a(S \cup \{a\}) \cdot B_b(R \cup \{b\})$ where $R = (W \setminus \{a,b\}) \setminus S$
- $\text{cs}_W(b) = \sum_{a \neq b} M_W[a,b]$ (column sum)
- $\Sigma(W) = \sum_{a \neq b} M_W[a,b]$ (total off-diagonal sum)

## Main Theorem

**Theorem (Key Identity).** For all $|W| = m \geq 2$ and $b \in W$:
$$B_b(W) + (-1)^m E_b(W) = 2r \cdot \text{cs}_W(b)$$

**Proof.** By strong induction on $m$.

### Base case: $m = 2$

Let $W = \{a, b\}$. Then $M_W[a,b] = E_a(\{a\}) \cdot B_b(\{b\}) = 1$, so $\text{cs}_W(b) = 1$.

$$B_b(W) + E_b(W) = t(b,a) + t(a,b) = (r - s_{ab}) + (r + s_{ab}) = 2r = 2r \cdot 1. \quad \checkmark$$

### Inductive step: $m \geq 3$

Assume the Key Identity holds for all vertex sets of size $< m$. Write $W' = W \setminus \{b\}$.

**Three recurrences** (algebraic identities from definitions):

**(R1)** First-edge decomposition of $B_b$:
$$B_b(W) = \sum_{w \in W'} t(b,w) \cdot B_w(W') = r \cdot T(W') + \sum_{w \in W'} s_{bw} \cdot B_w(W')$$

**(R2)** Last-edge decomposition of $E_b$:
$$E_b(W) = \sum_{w \in W'} E_w(W') \cdot t(w,b) = r \cdot T(W') - \sum_{w \in W'} s_{bw} \cdot E_w(W')$$

(The sign flip uses $t(w,b) = r + s_{wb} = r - s_{bw}$.)

**(R3)** Column sum from the column recurrence of $M$:
$$\text{cs}_W(b) = (-1)^{m-2} T(W') + r \cdot \Sigma(W') + \sum_{w \in W'} s_{bw} \cdot \text{cs}_{W'}(w)$$

(Derived by summing $M_W[a,b] = (-1)^{m-2} E_a(W') + \sum_w t(b,w) M_{W'}[a,w]$ over $a \neq b$.)

**Combining (R1) and (R2):**

$$B_b + (-1)^m E_b = r T(W')[1 + (-1)^m] + \sum_w s_{bw}\bigl[B_w(W') + (-1)^{m-1} E_w(W')\bigr]$$

where we used $(-1)^m \cdot (-1) = (-1)^{m-1}$.

By the **inductive hypothesis** at size $|W'| = m - 1$:

$$B_w(W') + (-1)^{m-1} E_w(W') = 2r \cdot \text{cs}_{W'}(w)$$

Therefore:

$$\text{LHS} = r T(W')[1 + (-1)^m] + 2r \sum_w s_{bw} \cdot \text{cs}_{W'}(w) \quad \cdots (*)$$

**From (R3):**

$$\text{RHS} = 2r \cdot \text{cs}_W(b) = 2r(-1)^{m-2} T(W') + 2r^2 \Sigma(W') + 2r \sum_w s_{bw} \cdot \text{cs}_{W'}(w) \quad \cdots (**)$$

**Setting $(*) = (**)$**, the $s_{bw} \cdot \text{cs}_{W'}(w)$ terms cancel, leaving:

$$r T(W')[1 + (-1)^m] = 2r(-1)^{m-2} T(W') + 2r^2 \Sigma(W')$$

Since $(-1)^{m-2} = (-1)^m$:

$$r T(W')[1 - (-1)^m] = 2r^2 \Sigma(W') \quad \cdots (\dagger)$$

**Closing the induction via the Sigma identity.**

Summing the Key Identity at size $m-1$ over all vertices of $W'$:

$$T(W')[1 + (-1)^{m-1}] = 2r \cdot \Sigma(W')$$

- **$m$ even:** $m - 1$ odd, so $1 + (-1)^{m-1} = 0$, giving $\Sigma(W') = 0$. Then $(\dagger)$: $0 = 0$. $\checkmark$
- **$m$ odd:** $m - 1$ even, so $1 + (-1)^{m-1} = 2$, giving $T(W') = r \Sigma(W')$. Then $(\dagger)$: $2rT(W') = 2r^2\Sigma(W') = 2rT(W')$. $\checkmark$

$\blacksquare$

## Corollaries

**Corollary 1 (Sigma Identity).** Summing the Key Identity over $b \in W$:
$$T(W)[1 + (-1)^m] = 2r \cdot \Sigma(W)$$
- Even $m$: $T(W) = r \cdot \Sigma(W)$
- Odd $m$: $\Sigma(W) = 0$

**Corollary 2 (Even r-powers).** $M_W[a,b]$ contains only even powers of $r$.

*Proof.* Path reversal: $B_b(W; -r) = (-1)^{m-1} E_b(W; r)$, so $\text{odd}_r(B_b) = [B_b + (-1)^m E_b]/2 = r \cdot \text{cs}_W(b)$.

Row recurrence: $M_W[a,b] = B_b(W \setminus \{a\}) - \sum_v t(v,a) M'[v,b]$ where $M'$ acts on size $m-1$.

By induction, $M'$ has even r-powers, so $\text{odd}_r(M[a,b]) = \text{odd}_r(B_b) - r \cdot \text{cs}(b) = 0$. $\blacksquare$

**Corollary 3 (Symmetry).** $M_W[a,b] = M_W[b,a]$.

*Proof.* The re-indexing identity (path reversal + S-complement swap): $M_T[b,a] = (-1)^{m-2} M_{T^{\text{op}}}[a,b]$, where $T^{\text{op}}$ replaces $s_{ij} \to -s_{ij}$.

Even r-powers $\Rightarrow$ each monomial $r^k \prod s_{ij}$ in $M[a,b]$ has $k$ even. Since total degree is $m-2$, the number of $s$-factors has parity $m-2 \pmod{2}$, so $M[a,b]\big|_{s \to -s} = (-1)^{m-2} M[a,b]$.

Therefore $M[b,a] = (-1)^{m-2} \cdot (-1)^{m-2} M[a,b] = M[a,b]$. $\blacksquare$

## Ingredients used

1. **Path reversal identity:** $B_b(W; -r) = (-1)^{m-1} E_b(W; r)$ — proved by substituting $-r$ for $r$, which negates $t(i,j)$ to $-t(j,i)$, reversing all paths.

2. **Re-indexing identity:** $M_T[b,a] = (-1)^{m-2} M_{T^{\text{op}}}[a,b]$ — proved by the $S \leftrightarrow R$ substitution in the definition of $M$.

3. **Column recurrence:** $M_W[a,b] = (-1)^{m-2} E_a(W \setminus \{b\}) + \sum_{w} t(b,w) M_{W \setminus \{b\}}[a,w]$ — derived by separating $b$'s contribution in the inclusion-exclusion.

All three are algebraic identities from the definitions, requiring no assumptions on the tournament.

## Computational verification

Verified symbolically (as polynomial identities in $r$ and the $s_{ij}$) for $m = 2, 3, 4, 5, 6$ with all vertices. See `04-computation/key_identity_complete_proof.py`.
