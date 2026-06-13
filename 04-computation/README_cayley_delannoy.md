# Cayley-Delannoy Tournament Toolkit

A Python library + CLI tools implementing the theory from "The Cayley Transform, Delannoy Paths, and the Fourier Analysis of Hamiltonian Paths in Random Tournaments."

## Quick Start

```python
from cayley_delannoy import Q, gk, cv2, W, T_matrix, address, golden_shadow, rodrigues

Q(Fraction(1, 3))       # = 2 (the binary alphabet lives at coupling 1/3)
gk(3, 5)                # = 85 (Delannoy weight)
cv2(10)                  # = Fraction(3121, 16200) ~ 0.193
W(7)                     # = 6350 (NUD-weighted permutation count)
T_matrix(3, 5)           # = 255 (symmetric Delannoy matrix, diagonal = OEIS A108666)
address(7)               # = Fraction(3, 4) (Cayley address)
golden_shadow(4)         # = 3.236... = phi^2 (metallic mean - 1)
```

## CLI Tools

```bash
# Is a tournament ranking statistically significant?
python3 tournament_test.py --n 7 --h 189
# => Z = +2.62, SIGNIFICANT

# Rank items from pairwise A/B test data
python3 ab_test_ranker.py --demo
# => Ranked list + significance score

# Just need the quick formula?
# CV = sqrt(2/n), Z = (H - n!/2^{n-1}) / (n!/2^{n-1} * sqrt(2/n))
```

## Core Library Functions

| Function | Description |
|---|---|
| `Q(x)` | Cayley transform (1+x)/(1-x) = exp(2 arctanh x) |
| `gk(k, m)` | Delannoy weight (exact Fraction arithmetic) |
| `cv2(n)` | Exact CV^2 for n-vertex random tournament |
| `W(n)` | NUD-weighted permutation count (bitmask DP, n <= 20) |
| `T_matrix(k, m)` | Symmetric Delannoy matrix k*gk(k,m) |
| `address(n)` | Cayley address x_n = (n-1)/(n+1) |
| `golden_shadow(n)` | f_n = (n-2+sqrt(n^2+4))/2 = metallic mean - 1 |
| `rodrigues(k, m)` | Rodrigues formula for gk |
| `asymptotic_cv2(n, terms)` | Fast CV^2 with k terms |

## Key Mathematical Results

All verified by the self-test suite (10/10 tests):

- **Master GF**: Q(x)^m = 1 + 2 sum gk(k,m) x^k
- **Delannoy identity**: k gk(k,m) = sum j C(k,j) C(m,j) 2^{j-1} = OEIS A108666 diagonal
- **Duality**: k gk(k,m) = m gk(m,k)
- **Parity**: gk(k,-m) = (-1)^k gk(k,m)
- **Functional equation**: even gk determined by odd gk via Q^m Q(-x)^m = 1
- **Rodrigues**: gk(k,m) = (1/k!) [d^k/du^k u^m(u+1)^{k-1}]_{u=1}
- **Wick rotation**: arctanh(i) = i pi/4 (pi emerges at imaginary coupling)
- **CV^2 = 2/n + 0/n^2 - 14/(3n^3)** (1/n^2 cancellation proved)
- **Golden EPs**: discriminant 4x(x^2-11x-1), EP eigenvalues = 1/phi, -phi (THM-224)

## Running Tests

```bash
python3 cayley_delannoy.py   # 10/10 tests should pass
```

## W(n) Sequence (NEW — not in OEIS)

```
W(n) = 1, 2, 8, 32, 158, 928, 6350, 49752, 439670, 4327904,
       46963358, 556953448, 7166360054, 99428495088,
       1479600188798, 23506712352248, 397095175477430,
       7107209383674112, 134345623603516190, 2674426516381764744
```

OEIS submission prepared at `03-artifacts/oeis/W_n_submission.txt`.

## File Index

| File | Description |
|---|---|
| `cayley_delannoy.py` | Core library (10 tests) |
| `tournament_test.py` | CLI ranking significance tool |
| `ab_test_ranker.py` | A/B test ranking from CSV |
| `setup_cayley.py` | pip install setup |
| Paper | `03-artifacts/drafts/cayley-delannoy-tournaments-v2.tex` |
| OEIS | `03-artifacts/oeis/W_n_submission.txt` |
| Cheat sheet | `03-artifacts/substack/cheat-sheet.tex` |
| Hooks A-I | `03-artifacts/substack/hook-{A..I}.tex` |

## The One-Sentence Summary

The Fourier energy of tournament Hamiltonian path counts is governed by ((1+x)/(1-x))^m, whose coefficients are Delannoy lattice path weights, whose exceptional points encode the golden ratio.
