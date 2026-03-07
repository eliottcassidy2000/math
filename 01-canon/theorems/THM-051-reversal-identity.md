# THM-051: Transfer Matrix Reversal Identity (M_T vs M_{T^op})

**Status:** PROVED (from path reversal + THM-050, verified n=3,4,5)
**Added:** opus-2026-03-06-S26

## Statement

**Theorem.** For any tournament T on n vertices, let T^op denote the opposite tournament (all arcs reversed). Then:

    M_{T^op}[a,a] = (-1)^{n-1} · M_T[a,a]     (diagonal)
    M_{T^op}[a,b] = (-1)^n · M_T[a,b]           (off-diagonal, a ≠ b)

## Proof

**Path reversal:** If P = (v_0, v_1, ..., v_{n-1}) is a Hamiltonian path in T, then P_rev = (v_{n-1}, ..., v_1, v_0) is a Hamiltonian path in T^op. This gives a bijection between Ham(T) and Ham(T^op).

**Diagonal:** Vertex a at position j in P corresponds to position (n-1-j) in P_rev. So P_T(a,j) = P_{T^op}(a, n-1-j), where P(a,j) = #{paths with a at position j}. Then:

    M_{T^op}[a,a] = Σ_j (-1)^j P_{T^op}(a,j) = Σ_j (-1)^j P_T(a,n-1-j)
                   = Σ_k (-1)^{n-1-k} P_T(a,k) = (-1)^{n-1} M_T[a,a]

**Off-diagonal (using THM-050):** The pair {a,b} at consecutive positions {j,j+1} in P corresponds to {a,b} at positions {n-2-j, n-1-j} in P_rev. So N_T(a,b,j) = N_{T^op}(a,b,n-2-j). Then:

    M_{T^op}[a,b] = Σ_j (-1)^j N_{T^op}(a,b,j) = Σ_j (-1)^j N_T(a,b,n-2-j)
                   = Σ_k (-1)^{n-2-k} N_T(a,b,k) = (-1)^n M_T[a,b]

## Corollaries

**At odd n:**
- M_{T^op}[a,a] = M_T[a,a] (diagonal preserved)
- M_{T^op}[a,b] = -M_T[a,b] for a ≠ b (off-diagonal negated)

**At even n:**
- M_{T^op}[a,a] = -M_T[a,a] (diagonal negated)
- M_{T^op}[a,b] = M_T[a,b] for a ≠ b (off-diagonal preserved)

**Corollary (Palindromic ⟹ Scalar at odd n):** If N(a,b,j) = N(a,b,n-2-j) for all pairs {a,b} and all j, then at odd n, M[a,b] = 0 for all a ≠ b.

## Verification

Exhaustively verified at n=3,4,5 for all tournaments.

## Source Files

- `04-computation/palindromic_n_scalar_m.py`
