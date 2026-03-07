# THM-061: Anti-Evaluation Identity W(-1/2) = (-1)^{n-1} H(T)

**Type:** Theorem (PROVED)
**Certainty:** 5 — PROVED (elementary + verified exhaustive n<=6, sample n=7)
**Status:** PROVED
**Added by:** kind-pasteur-2026-03-07-S26
**Tags:** #W-polynomial #Hamiltonian #evaluation #master-polynomial

---

## Statement

**Theorem.** For any tournament T on n vertices,

W(-1/2) = (-1)^{n-1} · H(T)

where W(r) is the weighted path polynomial and H(T) is the number of Hamiltonian paths.

**Corollary.** F_f(-1/2) = (-1)^f for all f >= 0, where F_f is the f-th master polynomial of THM-059.

---

## Proof (Elementary)

W(r) = sum_{P in S_n} prod_{i=0}^{n-2} (r + s_i), where s_i = A[P(i), P(i+1)] - 1/2 in {+1/2, -1/2}.

**At r = 1/2:** Each factor (1/2 + s_i) equals 1 if A[P(i),P(i+1)]=1 (forward edge) and 0 otherwise. So W(1/2) = #{permutations P where P(i)->P(i+1) for all i} = H(T).

**At r = -1/2:** Each factor (-1/2 + s_i) equals 0 if A[P(i),P(i+1)]=1, and -1 if A[P(i),P(i+1)]=0 (backward edge). So W(-1/2) = (-1)^{n-1} · #{permutations P where P(i+1)->P(i) for all i}.

The condition P(i+1)->P(i) for all i means P(n-1) -> P(n-2) -> ... -> P(0) is a Hamiltonian path. Reversing any such path gives a Hamiltonian path P(0) -> P(1) -> ... -> P(n-1), and vice versa. This is a bijection, so the count equals H(T).

Therefore W(-1/2) = (-1)^{n-1} · H(T). QED.

---

## Proof of Corollary (via Eulerian numbers)

F_f(r) = sum_{k=0}^f A(f+1, k) · (r+1/2)^{f-k} · (r-1/2)^k (THM-059).

At r = -1/2: (r+1/2) = 0, (r-1/2) = -1. Only the k=f term survives:
F_f(-1/2) = A(f+1, f) · 0^0 · (-1)^f = 1 · (-1)^f = (-1)^f.

(Since A(n, n-1) = 1: the unique permutation with n-1 descents is the decreasing permutation.)

---

## Consistency Check

From THM-059: W(r) = F_{n-1}(r) + sum_I 2^{parts(I)} · F_{f(I)}(r) · I(T).

At r = -1/2:
W(-1/2) = (-1)^{n-1} + sum_I 2^{parts(I)} · (-1)^{f(I)} · I(T)

Since f(I) = (n-1) - 2|pi_I|, we have (-1)^{f(I)} = (-1)^{n-1} · (-1)^{-2|pi_I|} = (-1)^{n-1}.

So W(-1/2) = (-1)^{n-1} · [1 + sum_I 2^{parts(I)} · I(T)] = (-1)^{n-1} · I(Omega(T), 2) = (-1)^{n-1} · H(T). ✓

---

## Parity Consequences

- **Odd n:** W has only even powers of r, so W(-r) = W(r). Hence W(-1/2) = W(1/2) = H(T).
- **Even n:** W has only odd powers of r, so W(-r) = -W(r). Hence W(-1/2) = -W(1/2) = -H(T).

Both are consistent with (-1)^{n-1} · H(T).

---

## Verification

| n | Tested | Failures |
|---|--------|----------|
| 3 | 2/2 (all) | 0 |
| 4 | 8/8 (all) | 0 |
| 5 | 64/64 (all) | 0 |
| 6 | 1024/1024 (all) | 0 |
| 7 | 100 (sample) | 0 |

---

## Connection to THM-060 (Bipartite Skeleton)

At odd n, W(-1/2) = W(1/2) = H(T). Since the GS skeleton bipartition is determined by t3 parity (THM-060), the anti-evaluation W(-1/2) sits on both sides symmetrically. The interesting object is W(0), the "tangent evaluation":

At n=5: W(0) = 1 - t3 + 2t5 (from THM-059 tangent number connection).
Result: W(0) = 0 for ALL odd-t3 GS tilings (100% vanishing on one side of bipartition).

This vanishing follows from parity: W(0) mod 2 = (1 + t3) mod 2, so odd t3 forces W(0) even. At n=5, W(0) happens to always be exactly 0 on the odd side (a stronger statement than just parity).

---

## Scripts

- `04-computation/W_anti_eval_proof.py` — verification and proof
- `04-computation/W_flip_analysis.py` — W(r) flip analysis
- `04-computation/W0_vanishing.py` — W(0) vanishing on bipartition
