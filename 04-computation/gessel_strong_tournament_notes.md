# Gessel et al. "Counting acyclic and strong digraphs by descents" (2020)

**Source:** Archer-Gessel-Graves-Liang, Discrete Mathematics 343 (2020)
**Reviewed by:** opus-2026-03-07-S45 (background agent)

## Key definitions

- **Descent of a digraph:** edge (s,t) with s > t (higher to lower label)
- **t_n(u):** descent polynomial for strong tournaments on [n]
  - t_1=1, t_2=0, t_3=u+u^2, t_4=u+6u^2+10u^3+6u^4+u^5

## Key results

1. **Recurrence:** t_n(u) = (1+u)^{C(n,2)} - sum_{k=1}^{n-1} C(n,k)_u * (1+u)^{C(n-k,2)} * t_k(u)
   where C(n,k)_u is the q-binomial (Gaussian binomial) with q=u.

2. **GF:** U(x) = 1/(1-T(x)) where T(x) = sum t_n(u) x^n / n!_u (Eulerian GF)
   Meaning: every tournament = sequence of strong components.

3. **Palindromicity:** t_n(u) = u^{C(n,2)} t_n(1/u) — from reversing all edges.

4. **Divisibility:** t_n(u) divisible by (1+u)^{floor(n/2)}.

## Connection to F(T,x)

- **Different statistics:** Gessel's des(T) counts backward edges globally (label-dependent).
  F(T,x) counts forward edges along Hamiltonian paths (path-dependent).
- **NOT directly equivalent**, but both palindromic for related reasons.
- The strong component decomposition U(x)=1/(1-T(x)) could give structural insights
  into H(T) since strong tournaments have richest Ham path structure.

## "Eulerian graphic GF" (new framework)

Combines descent tracking (u) and edge tracking (y) via q=(1+uy)/(1+y).
- Setting u=1: ordinary graphic GF
- Setting y=0: Eulerian GF
- Connects Wright's strong digraph counting to strong tournament counting:
  eta_n(y^{-2},y) = (1+y)^{C(n,2)} t_n(y^{-1})

## Relevance to our project

- The (1+u)^{floor(n/2)} divisibility could relate to OCF's factor-2 structures
- Strong tournament enumeration by descents may connect to F(T,x) via relabeling
- The framework is the natural algebraic home for descent-type statistics on tournaments
