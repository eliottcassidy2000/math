# THM-018: Coefficient Identity alpha_w^H = alpha_w^I

**Type:** Theorem (proved symbolically at n<=8; numerical at n<=9)
**Certainty:** 5 -- PROVED at n=4,...,8 (symbolic); OPEN for general n
**Status:** PROVED at n<=8; OPEN for general n
**Added by:** opus-2026-03-05-S5/S5b
**Tags:** #ocf #arc-reversal #coefficient-matching #open-q-009 #claim-a

---

## Statement

For n = 4, 5, 6, 7, 8 and any tournament T on n vertices with s_w = 0 for all w in W = V\{i,j}:

**alpha_w^H = alpha_w^I**

as a polynomial identity in the arc variables, where:

- alpha_w^H = d(delta_H)/d(s_w)|_{s=0} is the coefficient of s_w in the Hamiltonian path difference delta_H
- alpha_w^I = d(delta_I)/d(s_w)|_{s=0} is the coefficient of s_w in the OCF expression delta_I (using the FULL inductive formula)

Since delta_H and delta_I are both linear in the s-variables (no cross terms s_a*s_b), and the coefficients match for each w, this proves **delta_H = delta_I** and hence **OCF** at each verified n.

**CRITICAL DISCOVERY (S5b):** The full inductive formula for alpha_w^I, which uses
H(comp) factors from inductive OCF, produces the EXACT same polynomial as alpha_w^H.
This means the OCF proof reduces to proving a single polynomial identity that holds
at every n checked (4,...,8 symbolic, 4,...,9 numerical).

---

## Key Formulas

### alpha_w^H (Insertion Decomposition)

For each permutation pi = (u_1, ..., u_{m-1}) of B_w = W\{w}, define:

**S(pi) = [T(w,u_1) + q_{u_1}] + sum_{k=1}^{m-2} [T(u_k,w)*q_{u_{k+1}} + p_{u_k}*T(w,u_{k+1})]/T(u_k,u_{k+1}) + [T(u_{m-1},w) + p_{u_{m-1}}]**

where p_v = T(v,I), q_v = T(J,v), with p_v + q_v = 1 at s=0.

Then: **alpha_w^H = -sum_pi B_wt(pi) * S(pi)**

where B_wt(pi) = T(u_1,u_2)*...*T(u_{m-2},u_{m-1}) is the Hamiltonian path weight of pi in B_w.

### alpha_w^I (Full Inductive Formula — ALL n)

**alpha_w^I = -2*H(B_w) + sum_{chains} 2*(-p or -q)*internal*H(comp)**

Specifically:

**alpha_w^I = -2*H(B_w) + sum over chains (w_1,...,w_k) of odd length k >= 3 with w at endpoint:**
- **If w = w_1:** contributes 2*(-p_{w_k})*T(w_1,w_2)*...*T(w_{k-1},w_k)*H(W\{w_1,...,w_k})
- **If w = w_k:** contributes 2*(-q_{w_1})*T(w_1,w_2)*...*T(w_{k-1},w_k)*H(W\{w_1,...,w_k})

where H(comp) = H(W\{chain vertices}) is the Hamiltonian path weight of the complement,
computed by INDUCTIVE OCF on sub-tournaments.

**Derivation:** From THM-013's general formula delta_I = sum_{k>=1} 2^k*Delta(alpha_k):
- The 3-cycle terms across ALL k levels telescope: sum_k 2^k*(-s_w)*alpha_{k-1}(B_w)
  = -2*s_w * sum_k 2^{k-1}*alpha_{k-1}(B_w) = -2*s_w * I(Omega(B_w),2) = -2*s_w*H(B_w)
- Similarly, L-cycle terms telescope to give 2*bracket*H(comp) for each chain.
- The H(comp) factor absorbs all higher-order cycle interactions via inductive OCF.

### The Defect Identity

The identity alpha_w^H = alpha_w^I is equivalent to:

**sum_pi B_wt(pi) * (S(pi) - 2) = -cycle_derivs(w)**

At n=4 (|B_w|=1): S(pi) = 2 always (trivial). No cycles. Identity is S=2.
At n>=5: S(pi) != 2 in general; the weighted defect gives cycle derivatives.

---

## Key Structural Insight

At s=0: **T(I,v) = T(J,v) = q_v** for all v in W.

This means i and j are interchangeable in the interface arcs at the s=0 specialization. This symmetry is crucial for the cancellations in the proof.

---

## Proof Method

Direct symbolic computation using SymPy. For each n:

1. Express alpha_w^H as a polynomial in the free variables (arcs T(a,b) for a,b in W, and interface arcs p_v)
2. Express alpha_w^I using the same variables
3. Compute expand(alpha_H - alpha_I) and verify it equals 0

### Variable counts by n:
| n | |B_w| | # internal arcs | # w-arcs | # interface | # perms | # terms (alpha_H) |
|---|-------|----------------|----------|-------------|---------|-------------------|
| 4 | 1 | 0 | 1 | 1 | 1 | 0 |
| 5 | 2 | 1 | 2 | 2 | 2 | ~10 |
| 6 | 3 | 3 | 3 | 3 | 6 | ~50 |
| 7 | 4 | 6 | 4 | 4 | 24 | 256 |

All cases verified: diff = 0 exactly (no simplification needed for n <= 7).

---

## Extension to n >= 8 (SOLVED)

The earlier n<=7 formula used a simplified delta_I that was incomplete at n>=8.
The FULL inductive formula (using H(comp) from inductive OCF) works at ALL n tested:

**n=8 PROVED symbolically** (1338 terms, diff = 0 exactly)
**n=9 confirmed numerically** (max error 1.05e-12)

The key insight: by applying inductive OCF to sub-tournaments, all Delta(alpha_k) terms
collapse into a single formula involving H(comp) factors. This eliminates the need to
separately handle VD cycle pairs, triples, etc.

---

## Verification Record

| n | Method | alpha_H terms | alpha_I terms | diff | Result |
|---|--------|-------------|-------------|------|--------|
| 4 | symbolic | 0 | 0 | 0 | PROVED |
| 5 | symbolic | ~10 | ~10 | 0 | PROVED |
| 6 | symbolic | ~50 | ~50 | 0 | PROVED |
| 7 | symbolic | 256 | 256 | 0 | PROVED |
| 8 | symbolic (full) | 1338 | 1338 | 0 | PROVED |
| 4-9 | numerical (s=0) | - | - | < 1e-12 | CONFIRMED |

---

## Code

- `04-computation/q009_symbolic_proof.py` — symbolic verification at n=4,5,6
- `04-computation/q009_symbolic_n7.py` — symbolic verification at n=7
- `04-computation/q009_alpha_identity.py` — insertion decomposition and factor analysis
- `04-computation/q009_insertion_decomp.py` — defect vs cycle derivatives verification
- `04-computation/q009_symbolic_n8.py` — symbolic verification at n=8 (full inductive formula)
- `04-computation/q009_inductive_proof.py` — numerical verification of full formula at n=4,...,9
- `04-computation/q009_coefficient_proof.py` — analytical formulas (earlier version)
