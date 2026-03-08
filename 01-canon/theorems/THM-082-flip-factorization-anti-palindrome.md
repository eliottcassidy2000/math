# THM-082: Arc-Flip Factorization and Anti-Palindromicity of D(x)

**Type:** Theorem (proved algebraically)
**Certainty:** 5 -- PROVED
**Status:** PROVED (algebraic proof from definitions; verified exhaustively at n=4,5)
**Added by:** opus-2026-03-07
**Tags:** #f-polynomial #arc-flip #palindrome #reversal-involution

---

## Setup

Let T be a tournament on n vertices with arc u->v. Let T' be the tournament
obtained by flipping this arc (replacing u->v with v->u). Define:

- **F(T, x)** = sum over Hamiltonian paths P of x^{fwd(P)}, where fwd(P) counts
  forward edges along P.

Partition Hamiltonian paths into three types:

- **Type (a):** P contains ...u,v... consecutively. Since u->v is an arc of T,
  this contributes one forward step. Let fwd(P) denote the total forward count in T.

- **Type (b):** P contains ...v,u... consecutively. Since v->u is NOT an arc of T
  (u->v is), this step is backward. Let fwd(P) denote the total forward count in T.

- **Type (c):** P does not contain u,v in either order consecutively. The flip
  does not change the forward count of P.

Define the generating functions (polynomials in x):

- **G_uv(x) = sum_{P of type (a)} x^{fwd(P)-1}** (shifted down by 1, removing the
  u->v forward step)

- **G_vu(x) = sum_{P of type (b)} x^{fwd(P)}**

- **D(x) = G_uv(x) - G_vu(x)**

Both G_uv and G_vu have degree at most n-2.

---

## Theorem 1: Arc-Flip Factorization

**F(T, x) - F(T', x) = (x - 1) * D(x)**

### Proof

For type (c) paths: fwd_T(P) = fwd_{T'}(P), so the contribution to F(T)-F(T') is 0.

For type (a) paths (containing ...u,v... consecutively):
- In T: the step u->v is forward, so contribution is x^{fwd(P)}.
- In T': the arc is now v->u, so the step u->v is backward. The path P still exists
  but now has fwd_{T'}(P) = fwd_T(P) - 1. Contribution to F(T') is x^{fwd(P)-1}.
- Difference: x^{fwd(P)} - x^{fwd(P)-1} = x^{fwd(P)-1}(x - 1).

For type (b) paths (containing ...v,u... consecutively):
- In T: the step v->u is backward (since u->v is the arc), so contribution is x^{fwd(P)}.
- In T': the arc is now v->u, so this step is forward. fwd_{T'}(P) = fwd(P) + 1.
  Contribution to F(T') is x^{fwd(P)+1}.
- Difference: x^{fwd(P)} - x^{fwd(P)+1} = -x^{fwd(P)}(x - 1).

Summing:

    F(T,x) - F(T',x) = (x-1) * [sum_{type (a)} x^{fwd(P)-1} - sum_{type (b)} x^{fwd(P)}]
                      = (x-1) * [G_uv(x) - G_vu(x)]
                      = (x-1) * D(x).                                              QED

---

## Theorem 2: Anti-Palindromicity of D(x)

**D_k = -D_{n-2-k}** for all 0 <= k <= n-2.

Equivalently, D(x) is anti-palindromic: x^{n-2} D(1/x) = -D(x).

### Proof

Consider the **reversal involution** on Hamiltonian paths:

    P = (v_1, v_2, ..., v_n) |-> P^rev = (v_n, v_{n-1}, ..., v_1).

This is an involution on the set of all Hamiltonian paths. Key properties:

1. **Forward count:** fwd(P^rev) = (n-1) - fwd(P), because each edge is either
   forward or backward in P, and it switches in P^rev.

2. **Consecutive pair mapping:** If P contains ...u,v... at position j (meaning
   P[j]=u, P[j+1]=v), then P^rev contains ...v,u... at position (n-2-j).

Property (2) means reversal maps type-(a) paths bijectively to type-(b) paths
and vice versa.

Now consider G_uv[k], which counts type-(a) paths with fwd(P) = k+1
(the shift accounts for removing the forward step u->v).

Under reversal, such a path P maps to P^rev which:
- Is type (b) (contains ...v,u... instead of ...u,v...)
- Has fwd(P^rev) = (n-1) - (k+1) = n-2-k

So P^rev contributes to G_vu[n-2-k]. Since reversal is a bijection between
type-(a) paths with fwd = k+1 and type-(b) paths with fwd = n-2-k, we have:

    **G_uv[k] = G_vu[n-2-k]**    for all 0 <= k <= n-2.

Therefore:

    D[k] = G_uv[k] - G_vu[k]
         = G_vu[n-2-k] - G_vu[k]
         = -(G_vu[k] - G_vu[n-2-k])
         = -(G_uv[n-2-k] - G_vu[n-2-k])
         = -D[n-2-k].                                                              QED

---

## Corollary 1: D(1) = 0

Since D is anti-palindromic of degree n-2, D(1) = sum_k D_k = 0
(pair D_k with D_{n-2-k} = -D_k; they cancel).

---

## Corollary 2: H(T) = H(T') under single arc flip

H(T) - H(T') = F(T,1) - F(T',1) = (1-1) * D(1) = 0.

Therefore H(T) is invariant under any single arc reversal.

**Remark:** This is well known -- it follows from the fact that every tournament
has an odd number of Hamiltonian paths (Redei's theorem), and flipping one arc
preserves the total count. But the polynomial factorization gives a stronger result:
the entire polynomial difference factors through (x-1).

---

## Corollary 3: (x-1)^2 divides F(T,x) - F(T',x)

Since D(1) = 0, we can write D(x) = (x-1) * E(x) for some polynomial E.
Therefore F(T,x) - F(T',x) = (x-1)^2 * E(x).

---

## Verification Record

| n | Arcs tested | G_uv[k]=G_vu[n-2-k] | Anti-palindromic D | Factorization |
|---|-------------|----------------------|--------------------|---------------|
| 4 | 384 (exhaustive) | 384/384 | 384/384 | 384/384 |
| 5 | 10240 (exhaustive) | 10240/10240 | 10240/10240 | 10240/10240 |

---

## Verification Script

`04-computation/flip_reversal_proof.py`
