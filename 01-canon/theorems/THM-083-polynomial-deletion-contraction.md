# THM-083: Polynomial Deletion-Contraction for F(T,x)

**Status:** PROVED (algebraic) + VERIFIED n=4,5 exhaustive
**Proved by:** kind-pasteur-2026-03-07-S35
**Scope:** All tournaments (some parts extend to all digraphs)

---

## Statement

Let T be a tournament on n vertices with edge e = (u → v). Define:

- **F_T(x)** = sum over all permutations P of {1,...,n} of x^{asc_T(P)}, where asc_T(P) = #{i : T[P_i, P_{i+1}] = 1} counts steps following the tournament direction.
- **T \ e** = T with edge u→v removed (a non-tournament digraph on n vertices)
- **T / e** = contraction merging u,v into w (convention: w gets IN from u, OUT from v)

Then:

**(A) Polynomial DC identity:**
$$F_T(x) = F_{T \setminus e}(x) + (x-1) \cdot F(T/e, x)$$

**(B) Key identification:** F(T/e, x) = G_{u,v}(x), where
$$G_{u,v}(x) = \sum_{\substack{P \in S_n \\ P_k = u, P_{k+1} = v \text{ for some } k}} x^{\text{asc}_T(P) - 1}$$

**(C) Arc-flip formula:** For T' = T with u→v flipped to v→u:
$$F_T(x) - F_{T'}(x) = (x-1) \cdot D(x)$$
where D(x) = F(T/e, x) - F(T'/e', x) is **anti-palindromic**: D(x) = -x^{n-2} D(1/x).

---

## Proof of (A)

Partition the set of all n! permutations into two classes:

**Class 1: P does not have u immediately before v.** These permutations have the same asc count in T and in T\e (the only difference is at the u→v position, which is not traversed). Contribution: F_{T\e}(x) minus the contribution of Class 2 permutations to F_{T\e}.

Actually, more directly:

For any permutation P:
- If P does NOT have u immediately before v: asc_T(P) = asc_{T\e}(P) (same edges traversed)
- If P HAS u immediately before v at position k: asc_T(P) = asc_{T\e}(P) + 1 (the u→v step contributes 1 to asc_T but 0 to asc_{T\e})

So:
$$F_T(x) = \sum_{P: \text{no } u \text{ before } v} x^{\text{asc}_T(P)} + \sum_{P: u \text{ before } v} x^{\text{asc}_T(P)}$$
$$= \sum_{P: \text{no } u \text{ before } v} x^{\text{asc}_{T\setminus e}(P)} + \sum_{P: u \text{ before } v} x^{\text{asc}_{T\setminus e}(P) + 1}$$
$$= F_{T\setminus e}(x) + \sum_{P: u \text{ before } v} [x^{\text{asc}_{T\setminus e}(P)+1} - x^{\text{asc}_{T\setminus e}(P)}]$$
$$= F_{T\setminus e}(x) + (x-1) \sum_{P: u \text{ before } v} x^{\text{asc}_{T\setminus e}(P)}$$

The remaining sum is $\sum x^{\text{asc}_{T\setminus e}(P)} = \sum x^{\text{asc}_T(P) - 1} = G_{u,v}(x)$.

## Proof of (B)

There is a bijection between permutations P of {1,...,n} with u at position k, v at position k+1, and permutations P' of V(T/e) = {w} ∪ (V\{u,v}) with w at position k:

φ(P) = (P_0, ..., P_{k-1}, w, P_{k+2}, ..., P_{n-1})

Under this bijection, the edge traversals are preserved:
- Step i < k-1 or i > k+1: same vertices, same adjacency
- Step k-1 to k: P_{k-1} → u in T maps to P_{k-1} → w in T/e, both valid iff T[P_{k-1}][u] = T/e[P_{k-1}][w] = 1 ✓
- Step k+1 to k+2: v → P_{k+2} in T maps to w → P_{k+2} in T/e, both valid iff T[v][P_{k+2}] = T/e[w][P_{k+2}] = 1 ✓

The u→v step (always forward, contributing 1) is removed, so:
asc_T(P) - 1 = asc_{T/e}(P')

Therefore G_{u,v}(x) = F(T/e, x). □

## Proof of (C)

Since T\e = T'\e' (same digraph: T without the u-v edge in either direction):

F_T(x) - F_{T'}(x) = [F_{T\e} + (x-1)F(T/e)] - [F_{T'\e'} + (x-1)F(T'/e')]
                    = (x-1) · [F(T/e) - F(T'/e')]

The anti-palindromicity of D = F(T/e) - F(T'/e') follows from the palindrome property of tournaments:

Both F_T and F_{T'} satisfy F(x) = x^{n-1} F(1/x) (palindrome for tournaments, via path reversal).

So f(x) := F_T(x) - F_{T'}(x) also satisfies f(x) = x^{n-1} f(1/x).

From f(x) = (x-1) · D(x):
(x-1)D(x) = x^{n-1}(1/x - 1)D(1/x) = -x^{n-2}(x-1)D(1/x)

Dividing by (x-1):
**D(x) = -x^{n-2} D(1/x)** (anti-palindrome) □

---

## Corollaries

1. **H(T) - H(T') = D_{n-2}** where D_{n-2} = leading coefficient of D(x). Since D is anti-palindromic with deg ≤ n-2, the constant term D_0 = -D_{n-2} = -(H(T)-H(T')).

2. **At x=1:** (1-1)·D(1) = 0, consistent with F_T(1) = F_{T'}(1) = n!.

3. **Sum formula:** F(T/e, x) + F(T'/e', x) is palindromic of degree n-2 (follows from F_T + F_{T'} palindromic).

4. **THM-082 recovery:** At x^{n-1} coefficient, F_T = H(T), and (x-1)·F(T/e) contributes F(T/e)_{n-2} = H(T/e) at degree n-1. So H(T) = H(T\e) + H(T/e).

---

## Verified

- Polynomial DC (A): n=4 exhaustive (384/384), n=5 exhaustive (10240/10240)
- Key identification (B): n=4 (384/384), n=5 (10240/10240)
- Anti-palindrome (C): n=4 (384/384), n=5 (10240/10240)
