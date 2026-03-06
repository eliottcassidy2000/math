# Even-Odd Split Lemma (Signed Position Identity)

**Author:** opus-2026-03-05-S4
**Status:** VERIFIED n=4 (exhaustive), n=5 (exhaustive 10240), n=6 (sampled 750)
**Equivalent to:** OCF / Claim A (via arc-flip induction)

---

## Statement

**Lemma (Even-Odd Split).** For any tournament T on V and any i,j in V with M = V\{i,j}:

sum_{S⊆M} (-1)^|S| [E_i(S)·B_j(M\S) - E_j(S)·B_i(M\S)] = 0

where E_v(S) = h_end(S∪{v}, v) and B_v(R) = h_start({v}∪R, v).

**Equivalent restatement (Signed Position Identity):**

sum_{P: i→j in P} (-1)^{pos_P(i)} = sum_{P': j→i in P'} (-1)^{pos_{P'}(j)}

where pos_P(i) is the 0-indexed position of i in path P.

---

## Key Properties

1. **Tournament-specific:** FAILS for general digraphs (282/500 at n=4). The proof MUST use T[a][b]+T[b][a]=1.

2. **Polynomial identity:** Holds for weighted tournaments (T[a][b] ∈ R with T[a][b]+T[b][a]=1). Verified to machine precision at n=4,5,6.

3. **Each side nonzero individually:** sum (-1)^|S| E_i(S)·B_j(R) ≠ 0 in general. Only the DIFFERENCE vanishes.

4. **Factored dependence:** E_i(S) depends only on arcs within S∪{i}. B_j(R) depends only on arcs within {j}∪R. No overlap (arcs between S and R unused by either factor).

---

## Proof for n=3 (|M|=1)

M = {x}. The alternating sum is:
(B_j({x}) - B_i({x})) - (E_i({x}) - E_j({x}))
= (T[j][x] - T[i][x]) - (T[x][i] - T[x][j])
= (T[j][x] + T[x][j]) - (T[i][x] + T[x][i])
= 1 - 1 = 0

Uses tournament property T[a][b]+T[b][a]=1 directly. QED.

## Proof for n=4 (|M|=2)

M = {x,y}. Using variables p=T[x][i], q=T[j][x], r=T[y][i], s=T[j][y], e=T[x][y]
and tournament substitutions, every monomial cancels. Verified algebraically.

---

## Recursion Structure

E_v(S) = sum_{u∈S} T[u][v]·E_u(S\{u}),  E_v(∅) = 1
B_v(R) = sum_{w∈R} T[v][w]·B_w(R\{w}),  B_v(∅) = 1

Define Ψ_a(b, N) = sum_{S⊆N} (-1)^|S| E_a(S)·B_b(N\S).
Then Φ(i,j,M) = Ψ_i(j,M) - Ψ_j(i,M).

**Induction hypothesis:** Φ(a,b,N) = 0 for |N| < m implies Ψ_a(b,N) = Ψ_b(a,N).

Using B-recursion: Φ(i,j,M) = sum_w [T[j][w]·Ψ_w(i,M_w) - T[i][w]·Ψ_w(j,M_w)].
By induction + reindexing, this reduces back to Ψ_i(j,M) - Ψ_j(i,M) = Φ (circular).

**Obstacle:** The recursion unfolds to itself. A non-recursive insight is needed.

---

## Bracket Structure (from T043)

The "double recursion" (expanding both E and B) gives:
Φ = sum_{u∈S, w∈R} B(u,w)·E_u(S\{u})·B_w(R\{w}) + boundary terms

where B(u,w) = T[u][i]·T[j][w] - T[u][j]·T[i][w].

Bracket by s-type (s_v = 1 - T[v][i] - T[j][v]):
| u\w  | M-(s=-1) | M+(s=+1) | Z1(s=0) | Z0(s=0) |
|-------|----------|----------|---------|---------|
| M-    |    1     |    0     |    0    |    1    |
| M+    |    0     |   -1     |    0    |   -1    |
| Z1    |    1     |   -1     |    0    |    0    |
| Z0    |    0     |    0     |    0    |    0    |

**When all s_v = 0:** ALL brackets vanish AND boundary terms vanish. Φ = 0 trivially.
This could be the base case for induction on #{v : s_v ≠ 0}.

---

## Connection to OCF Proof

Proving this lemma for all n would prove OCF (H(T) = I(Ω(T), 2)) and hence Claim A.
The arc-flip induction: start from transitive tournament (H=1=I), flip arcs one at a time,
show delta_H = delta_I at each step via this lemma (THM-015).

---

## Verification Record

| n | Method | Result |
|---|--------|--------|
| 3 | Hand proof | ✓ |
| 4 | Algebraic (every monomial cancels) | ✓ |
| 4 | Exhaustive 384/384 | ✓ |
| 5 | Exhaustive 10240/10240 | ✓ |
| 6 | Sampled 750/750 | ✓ |
| 4-6 | Weighted (continuous) 500+500+100 | ✓ |
| 4-6 | Digraphs (non-tournament) | FAIL (expected) |
| 5-8 | Sampled (opus-S4) | ✓ |

## Significance

1. delta_H = adj'(j,i) - adj(i,j) = sum_S Delta(S, others\S).
2. By the Even-Odd Split, delta_H = 2 * sum_{|S| odd} Delta(S, R).
3. The odd-S terms correspond to odd cycle lengths (|S|+2 is odd when |S| is odd).
4. The cycle formula (A-clique) gives: delta_I = 2 * sum_{|S| odd} [g(S)-l(S)] * H(R).
5. So OCF reduces to: sum_{|S| odd} Delta(S, R) = sum_{|S| odd} [g(S) - l(S)] * H(R).

Source: opus-2026-03-05-S4, opus-2026-03-05-S4b
