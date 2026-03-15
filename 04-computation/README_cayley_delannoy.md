# cayley_delannoy.py — Cayley-Delannoy Tournament Library

A clean Python library implementing the theory from "The Cayley Transform, Delannoy Paths, and the Fourier Analysis of Hamiltonian Paths in Random Tournaments."

## Quick Start

```python
from cayley_delannoy import Q, gk, cv2, W, T_matrix, address, golden_shadow, rodrigues

# The Cayley transform
Q(Fraction(1, 3))  # = 2

# Delannoy weight
gk(3, 5)  # = 85

# CV^2 of Hamiltonian path count for n-vertex random tournament
cv2(10)  # = Fraction(3121, 16200) ≈ 0.1927

# W(n): the NUD-weighted permutation count (NEW, not in OEIS)
W(7)  # = 6350

# Symmetric Delannoy matrix (diagonal = OEIS A108666)
T_matrix(3, 5)  # = 255 = 3*gk(3,5) = 5*gk(5,3)

# Cayley address of natural number n
address(7)  # = Fraction(3, 4)

# Golden shadow (quadratic irrational with CF [n-1; n,n,n,...])
golden_shadow(4)  # = 3.236... = phi^2
```

## Core Functions

| Function | Description |
|---|---|
| `Q(x)` | Cayley transform (1+x)/(1-x) |
| `gk(k, m)` | Delannoy weight (exact, Fraction) |
| `cv2(n)` | Exact CV^2 for n-vertex random tournament |
| `W(n)` | NUD-weighted permutation count (bitmask DP) |
| `T_matrix(k, m)` | Symmetric Delannoy matrix k*gk(k,m) |
| `address(n)` | Cayley address (n-1)/(n+1) |
| `golden_shadow(n)` | (n-2+sqrt(n^2+4))/2 |
| `rodrigues(k, m)` | Rodrigues formula for gk |
| `asymptotic_cv2(n, terms)` | Fast CV^2 approximation |

## Key Identities (all verified by self-test)

- **Master GF**: Q(x)^m = 1 + 2*sum gk(k,m)*x^k
- **Duality**: k*gk(k,m) = m*gk(m,k)
- **Functional equation**: gk(2k,m) = sum (-1)^{j+1} gk(j,m)*gk(2k-j,m)
- **Rodrigues**: gk(k,m) = (1/k!) [d^k/du^k u^m(u+1)^{k-1}]_{u=1}
- **CV^2 vs W**: W(n)/n! = 1 + cv2(n)

## Running Tests

```bash
python3 cayley_delannoy.py
```

All 8 tests should pass: Q, gk, duality, parity, Rodrigues, functional equation, CV2 vs W, W values.

## References

- Paper: `03-artifacts/drafts/cayley-delannoy-tournaments-v2.tex`
- OEIS A108666: diagonal of T_matrix
- W(n) sequence: not yet in OEIS (submission prepared)
