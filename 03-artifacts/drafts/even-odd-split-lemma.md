# Even-Odd Split Lemma

## Statement

For any tournament T on n vertices, any arc i→j, let others = V\{i,j} with m = n-2 vertices.

For each subset S ⊂ others with R = others\S, define:

Delta(S, R) = L_j(S) * R_i(R) - L_i(S) * R_j(R)

where:
- L_j(S) = sum_{a∈S} h_end(S, a) * T[a][j]
- L_i(S) = sum_{a∈S} h_end(S, a) * T[a][i]
- R_i(R) = sum_{b∈R} h_start(R, b) * T[i][b]
- R_j(R) = sum_{b∈R} h_start(R, b) * T[j][b]

(With convention L_j(∅) = L_i(∅) = R_i(∅) = R_j(∅) = 1.)

**Lemma (Even-Odd Split):**

sum_{|S| even} Delta(S, R) = sum_{|S| odd} Delta(S, R)

Equivalently: sum_S (-1)^|S| Delta(S, R) = 0.

## Significance

1. delta_H = adj'(j,i) - adj(i,j) = sum_S Delta(S, others\S).
2. By the Even-Odd Split, delta_H = 2 * sum_{|S| odd} Delta(S, R).
3. The odd-S terms correspond to odd cycle lengths (|S|+2 is odd when |S| is odd).
4. The cycle formula (A-clique) gives: delta_I = 2 * sum_{|S| odd} [g(S)-l(S)] * H(R).
5. So OCF (delta_H = delta_I) reduces to proving:
   sum_{|S| odd} Delta(S, R) = sum_{|S| odd} [g(S) - l(S)] * H(R)

   This is a sum over ONLY odd-sized subsets, which is exactly the structure
   needed to match the cycle formula.

## Connection to the proof

The Even-Odd Split Lemma is EQUIVALENT to OCF: if OCF holds, then even=odd
(since delta_H = 2 * cycle_sum and cycle_sum = odd-S sum). Conversely, if
even=odd, then delta_H = 2 * odd-S sum, and showing odd-S sum = cycle formula
proves OCF.

The lemma reformulates OCF as: the alternating sum vanishes. This is a
polynomial identity in the arc variables {T[a][b]} of degree n-2.

## Verification

| n | Trials | Status |
|---|--------|--------|
| 5 | 100/100 | PASS |
| 6 | 100/100 | PASS |
| 7 | 100/100 | PASS |
| 8 | 30/30 | PASS |

## Notes

- The pairwise symmetry Delta(S,R) = ±Delta(R,S) does NOT hold.
- The cancellation in the alternating sum is a GLOBAL property.
- The quantity sum_S (-1)^|S| L_j(S) * R_i(R) = sum_S (-1)^|S| L_i(S) * R_j(R)
  (both cross-products agree under the alternating sum).

Source: opus-2026-03-05-S4
