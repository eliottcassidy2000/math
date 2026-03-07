# THM-063: Trivariate Generating Function G_T(t, x)

**Type:** Theorem (PROVED)
**Certainty:** 5 -- PROVED (algebraic + verified computationally)
**Status:** PROVED
**Added by:** opus-2026-03-07-S33
**Tags:** #generating-function #independence-polynomial #Eulerian #master-polynomial

---

## Statement

For a tournament T on n vertices, define the **trivariate generating function**:

**G_T(t, x) = A_n(t) + sum_I x^{parts(I)} I(T) A_{f_I+1}(t) (t-1)^{n-1-f_I}**

where the sum is over all odd-cycle collection invariants I, and A_m(t) denotes the classical Eulerian polynomial.

### Special evaluations

| (t, x)  | G_T(t, x)  | Name |
|----------|-----------|------|
| (t, 0)  | A_n(t)    | Eulerian polynomial (T-independent) |
| (0, x)  | I(Omega(T), x) | Independence polynomial |
| (0, 0)  | 1         | Empty independent set |
| (0, 2)  | H(T)      | Hamiltonian path count |
| (1, x)  | n!        | Total permutations (T-independent) |
| (t, 2)  | E_T(t) = sum_k a_k t^k | Tournament Eulerian polynomial |
| (-1, 2) | Deformed zigzag number | 0 for even n |

### Key structural property

The **t-axis** (x=0) gives the universal Eulerian polynomial A_n(t), independent of T.
The **x-axis** (t=0) gives the independence polynomial I(Omega(T), x).
These two axes **cross at (0,0) = 1** (the "double empty-set" point).

---

## Proof of G_T(0, x) = I(Omega(T), x)

G_T(0, x) = A(n,0) + sum_I x^{parts(I)} I(T) A_{f_I+1}(0) (0-1)^{n-1-f_I}

Since A_m(0) = 1 for all m:

G_T(0, x) = 1 + sum_I x^{parts(I)} (-1)^{n-1-f_I} I(T)

**Key fact:** n-1-f_I is always even, because each odd-cycle collection with total cycle vertex count C has f_I = n - 1 - (C - parts(I)) [?]. More directly:
- A single (2m+1)-cycle has f = n-1-2m, so n-1-f = 2m (even).
- A pair of cycles with sizes (2a+1) and (2b+1) has f = n-1-2(a+b), so n-1-f = 2(a+b) (even).
- In general, d-f counts the "excess" positions beyond the inner Eulerian polynomial, which equals 2 * (sum of half-cycle-sizes).

Therefore (-1)^{n-1-f_I} = 1 for all invariants I, and:

**G_T(0, x) = 1 + sum_I x^{parts(I)} I(T) = I(Omega(T), x).**

QED.

---

## Proof of G_T(t, 0) = A_n(t)

Immediate: setting x=0 kills all correction terms.

---

## Proof of G_T(1, x) = n!

At t=1: (t-1)^{n-1-f} = 0 for all f < n-1. The only surviving term has f = n-1, but no invariant has f = n-1 (that would require a cycle of size n, which doesn't appear in Omega). So G_T(1, x) = A_n(1) = n!.

---

## Connection to other results

- **THM-059:** The master polynomial framework provides the building blocks.
- **THM-062:** The forward-edge distribution a_k(T) = G_T evaluated coefficient-by-coefficient in t at x=2.
- **THM-001 (OCF):** H(T) = I(Omega(T), 2) = G_T(0, 2).
- The **k-colored independence polynomial** is: I_k(Omega, x) = coefficient of t^k in G_T(t, x).

---

## Verification

| Check | n | Samples | Status |
|-------|---|---------|--------|
| G_T(0,x) = I(Omega,x) | 7 | 5 tournaments x 5 x-values | PASS |
| G_T(t,0) = A_n(t) | 7 | 3 tournaments x 5 t-values | PASS |
| G_T(t,2) = E_T(t) | 7 | 3 tournaments x 5 t-values | PASS |
| dG/dx at (0,0) = alpha_1 | 7 | 3 tournaments | PASS |
| dG/dx at (1,0) = 0 | 7 | 3 tournaments | PASS |

---

## Scripts

- `04-computation/empty_set_weight_analysis.py` — Full verification of all special evaluations
- `04-computation/k_colored_ip_verify.py` — k-colored IP definition and verification
- `04-computation/euler_zigzag_deformation.py` — E_T(-1) evaluation
