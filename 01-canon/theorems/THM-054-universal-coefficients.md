# THM-054: Universal Coefficient Formula

**Type:** Theorem (verified computationally, proof sketch available)
**Certainty:** 4 -- VERIFIED computationally n=5,7,9
**Status:** VERIFIED
**Added by:** opus-2026-03-06-S11b (continued³), building on opus-S26
**Tags:** #transfer-matrix #r-polynomial #3-cycles #universal

---

## Statement

**Theorem.** Let M(r) = sum_{j=0}^{floor((n-1)/2)} c_{2j} * r^{2j} be the even r-power decomposition of the transfer matrix. Then:

**(a)** c_{n-1} = (n-1)! * I for ALL tournaments on n vertices.

**(b)** tr(c_{n-3}) = 2 * (n-2)! * (t_3 - C(n,3)/4)

where t_3 = #directed 3-cycles in T and C(n,3)/4 is the expected 3-cycle count in a random tournament.

**(c)** For regular tournaments at odd n: c_{n-3} = (n-1)!/4 * I (scalar, universal).

---

## Proof of (a)

The coefficient of r^{n-1} in M[a,a] comes from the term where ALL edge weights contribute their r component (no s_{ij} contributions). By THM-053:

c_{n-1}[a,a] = sum_{pi in S_n} (-1)^{pos(a,pi)} * [product of r's]

Since each path has n-1 edges, each contributing r, the r^{n-1} coefficient counts: for every ordering of n vertices, (-1)^{pos(a)}. There are (n-1)! orderings for each position of a, so:

c_{n-1}[a,a] = (n-1)! * sum_{k=0}^{n-1} (-1)^k = (n-1)! at odd n, 0 at even n.

For off-diagonal M[a,b] (a != b): total degree is n-2, so c_{n-1} contributes only on the diagonal, giving c_{n-1} = (n-1)! * I.

---

## Verification of (b)

| n | Formula | Tested | Max error |
|---|---------|--------|-----------|
| 5 | tr(c_2) = 12*t_3 - 30 | All 12 iso classes | 0 |
| 7 | tr(c_4) = 240*t_3 - 2100 | 50 iso classes (sampled) | 0 |
| 9 | tr(c_6) = 10080*t_3 - 211680 | 8 random tournaments | < 1e-10 |

General pattern: coefficient of t_3 is 2*(n-2)! (verified n=5,7,9).
Constant is -(n-2)!*C(n,3)/2 (verified n=5,7,9).

---

## Proof of (c)

For regular tournaments: all out-degrees equal (n-1)/2, so t_3 = C(n,3) - n*C((n-1)/2, 2). This is fixed for all regular tournaments at given n. Combined with (b) and the fact that c_{n-3} is scalar for regular (verified computationally), we get c_{n-3}[v,v] = 2*(n-2)!*(t_3(reg) - C(n,3)/4)/n = (n-1)!/4.

---

## Cross-scale universality boundary

The number of coefficients c_{2k} that are universal scalars for regular tournaments:

| n | Universal | Non-universal | Boundary |
|---|-----------|---------------|----------|
| 5 | c_4 | c_0 | c_2 is boundary (function of t_3 only) |
| 7 | c_2, c_4, c_6 | c_0 | c_0 encodes iso class |
| 9 | c_4, c_6, c_8 | c_0, c_2 | c_2 varies within regular |

---

## Dependencies

- THM-030 (transfer matrix symmetry / even r-powers)
- THM-053 (diagonal signed position formula)
- THM-027 (trace formula)

## Scripts

- `04-computation/universal_coefficient_formula.py`
- `04-computation/c0_concentration_theorem.py`
- `04-computation/c0_concentration_n9.py`
