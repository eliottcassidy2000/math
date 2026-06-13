# THM-049: Transfer Matrix r-Swap Identity

**Status:** PROVED
**Instance:** opus-2026-03-06-S11 (corrected proof)
**Dependencies:** Path reversal, B_v = (-1)^{|V|-1} E_v(-r)

## Statement

For the c-tournament transfer matrix M[a,b] defined by
$$M[a,b](r) = \sum_{S \subset U} (-1)^{|S|} E_a(S \cup \{a\}; r) \cdot B_b((U \setminus S) \cup \{b\}; r)$$
where $U = V \setminus \{a,b\}$, we have:

$$M[b,a](r) = M[a,b](-r)$$

## Proof

**Step 1: Path reversal identity.**
$B_v(V; r) = (-1)^{|V|-1} E_v(V; -r)$

A Hamiltonian path $v_1 \to \cdots \to v_k \to v$ ending at $v$ has weight
$\prod (r + s_{v_i, v_{i+1}})$. The reversed path $v \to v_k \to \cdots \to v_1$
has weight $\prod (r - s_{v_i, v_{i+1}}) = (-1)^k \prod (-r + s_{v_i, v_{i+1}})$
$= (-1)^{|V|-1} \cdot [\text{original weight at } r \to -r]$.

**Step 2: Rewrite M using path reversal.**
Since $|R \cup \{b\}| - 1 = |U| - |S|$:

$$M[a,b](r) = (-1)^{n-2} \sum_S E_a(S \cup \{a\}; r) \cdot E_b((U \setminus S) \cup \{b\}; -r)$$

(The $(-1)^{|S|}$ sign and $(-1)^{|U|-|S|}$ from B→E conversion combine to give $(-1)^{|U|} = (-1)^{n-2}$ with no $|S|$-dependence.)

**Step 3: Compute M[a,b](-r).**
$$M[a,b](-r) = (-1)^{n-2} \sum_S E_a(S \cup \{a\}; -r) \cdot E_b((U \setminus S) \cup \{b\}; r)$$

**Step 4: Compute M[b,a](r) and compare.**
$$M[b,a](r) = (-1)^{n-2} \sum_S E_b(S \cup \{b\}; r) \cdot E_a((U \setminus S) \cup \{a\}; -r)$$

Substituting $T = U \setminus S$:
$$= (-1)^{n-2} \sum_T E_b((U \setminus T) \cup \{b\}; r) \cdot E_a(T \cup \{a\}; -r)$$

This is identical to M[a,b](-r) from Step 3 (rename $S \to T$, commute factors). **QED.**

## Corollary

The following three properties are **equivalent**:
1. $M[a,b] = M[b,a]$ (transfer matrix symmetry)
2. $M[a,b](r) = M[a,b](-r)$ (only even powers of r)
3. Odd-power r-coefficients of $M[a,b]$ vanish

Proving any one implies the other two.

## Note on Sign

An earlier claim (kind-pasteur-S23b) stated $M[b,a] = (-1)^{n-2} M[a,b](-r)$.
This is INCORRECT. The correct identity has **no sign factor**. The error arose from
omitting the $(-1)^{|R|}$ factor in the $B \to E$ conversion (Step 2). The sign factors
$(-1)^{|S|}$ and $(-1)^{|U|-|S|}$ combine to give $(-1)^{|U|}$, which is independent
of $S$ and absorbs into the global $(-1)^{n-2}$ prefactor that cancels between
M[b,a] and M[a,b](-r).

## Verification

| n | M[a,b]=M[b,a] | M[b,a]=M[a,b](-r) | Even powers only |
|---|---|---|---|
| 3 | YES | YES | YES |
| 4 | YES | YES | YES |
| 5 | YES | YES | YES |
| 6 | YES | YES | YES |
| 7 | YES (numeric) | YES (numeric) | YES (numeric) |
| 8 | YES (numeric) | YES (numeric) | YES (numeric) |

## Additional Results

- **Leading coefficient:** $[r^{n-2}]$ of $M[a,b]$ = $(n-2)!$ if $n$ even, $0$ if $n$ odd.
  Formula: $\sum_{k=0}^{m} (-1)^k \binom{m}{k} k! (m-k)! = m! \cdot [1+(-1)^m]/2$.
- **Key algebraic identity:** $T(-r) = -T^T$ where $T$ is the c-tournament weight matrix.
