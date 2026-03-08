# Variable: M[a,b] — Transfer Matrix

**Symbol:** `M[a,b]` or `M(T)` or `M(r)`
**Type:** n x n integer matrix (at r=1/2); polynomial matrix in r
**Defined in:** THM-053, THM-054, THM-M

## Definition
M[a,b] = sum over S subset of V\{a,b} of (-1)^|S| * E_a(S) * B_b(V\{a,b}\S)

where:
- E_a(S) = number of Hamiltonian paths in T[S union {a}] ending at a
- B_b(R) = number of Hamiltonian paths in T[{b} union R] starting at b

For the r-parameterized version M(r):
M(r)[a,b] = sum_S (-1)^|S| * E_a(S,r) * B_b(R,r)
where weights use prod(r + s_e).

## Key properties (ALL VERIFIED computationally)

| Property | Status | Verified through | Source |
|----------|--------|-----------------|--------|
| M is symmetric: M[a,b] = M[b,a] | VERIFIED | n=8 (7500+ tests) | THM-M, opus-S35 |
| M(r) symmetric for ALL r | VERIFIED | n=5 | opus-S35c11 |
| tr(M) = H for odd n | VERIFIED | n=7 | THM-053 |
| M[a,a] = sum_P (-1)^{pos(a,P)} | PROVED | all n | THM-053 |
| M(r) = M(-r) at odd n | VERIFIED | n=5 | opus-S35c11 |
| Each M(r)[a,b] is even function of r at odd n | VERIFIED | n=5 | opus-S35c11 |

## What does NOT hold
- M(T) != M(T^op) in general (reversal_proof_attempt.py)
- M(T) != M(T^op)^T in general
- M is NOT a linear combination of I, A, A^T, A^2 etc. (simple_pattern_search.py)
- M does NOT commute with A or B (M_functional_equation.py)
- M(T\e) is NOT symmetric for deletion (dc_symmetry_path.py)
- M(T/e) is NOT symmetric for contraction

## Values at n=3

For C_3 (directed 3-cycle 0->1->2->0):
```
M = [[1, 0, 0],
     [0, 1, 0],
     [0, 0, 1]]
```
H = tr(M) = 3.

## Diagonal structure at n=7
All diagonal entries are ODD (verified by M_n7_structure.py).

## Equations
- tr(M) = H (odd n) — [hamiltonian-paths.md](hamiltonian-paths.md)
- tr(c_{n-3}) = 2*(n-2)!*(t_3 - C(n,3)/4) — THM-054
- Coefficient hierarchy: each level adds OCF invariants — THM-055

## Highest priority
PROVING M is symmetric would prove OCF. This is INV-001.

## Tags
#transfer-matrix #symmetry #core #proof-target #INV-001
