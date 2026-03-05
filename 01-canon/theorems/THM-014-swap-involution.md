# THM-014: Swap Involution for Arc-Flip Delta

**Type:** Lemma (verified n=4,...,7 computationally)
**Certainty:** 4 -- VERIFIED computationally, proof in progress
**Status:** VERIFIED
**Added by:** opus-2026-03-05-S2
**Tags:** #ocf #arc-reversal #involution #open-q-009

---

## Statement

Let T be a tournament, T' obtained by flipping arc i->j to j->i.

Define the **swap map** on Hamiltonian paths: for a T-path
  pi = (v_1,...,v_{a-1}, i, j, v_{a+2},...,v_n)
(using arc i->j at position (a, a+1)), define
  swap(pi) = (v_1,...,v_{a-1}, j, i, v_{a+2},...,v_n)

**Fact 1:** swap(pi) is a valid T'-path iff:
- (a = 1 OR T[v_{a-1}][j] = 1)  [predecessor of i must also beat j]
- (a+1 = n OR T[i][v_{a+2}] = 1)  [i must beat successor of j]

**Fact 2:** The swap is an involution on matched paths:
if swap(pi) is valid, then reverse-swap of swap(pi) recovers pi.

**Consequence:** adj(i,j) - adj'(j,i) = #U_T - #U_{T'}
where U_T = {unmatched T-paths} and U_{T'} = {unmatched T'-paths}.

---

## Blocking structure

A T-path (..., x, i, j, y, ...) is **unmatched** iff:
- **Pred-blocked:** T[x][j] = 0 (j beats x), equivalently s_x = -1
- **Succ-blocked:** T[i][y] = 0 (y beats i), equivalently s_y = -1

A T'-path (..., x, j, i, y, ...) is **unmatched** iff:
- **Pred-blocked:** T[x][i] = 0 (i beats x), equivalently s_x = +1
- **Succ-blocked:** T[j][y] = 0 (y beats j), equivalently s_y = +1

**Key observation:** Unmatched T-paths are blocked only by vertices
with s = -1, and unmatched T'-paths only by vertices with s = +1.

---

## Inclusion-exclusion decomposition

By inclusion-exclusion:
  #U_T = sum_{x: s=-1} nadj_T(x,i,j) + sum_{y: s=-1} nadj_T(i,j,y)
       - sum_{x,y: s=-1, x!=y} nadj_T(x,i,j,y)

  #U_{T'} = sum_{x: s=+1} nadj_{T'}(x,j,i) + sum_{y: s=+1} nadj_{T'}(j,i,y)
           - sum_{x,y: s=+1, x!=y} nadj_{T'}(x,j,i,y)

where nadj_T(a,b,c,...) = #{T-paths with consecutive subsequence (a,b,c,...)}.

All nadj terms can be expressed using sub-tournament Ham path counts:
  nadj_T(x,i,j) = sum_{S subset B_x} h_end(T[S+{x}], x) * h_start(T[{j}+R], j)
                 where R = B_x \ S.

**Proof goal:** Show this equals -2*sum(s_x*H(B_x)) + 2*sum_{L>=5}(DL-CL).

---

## Connection to 3-cycle terms

For x with s_x = -1: nadj_T(x,i,j) + nadj_T(i,j,x) counts paths
containing the triple (x,i,j) or (i,j,x). These paths "use" both the
arc i->j AND an arc involving x and {i,j}.

For x with s_x = +1: nadj_{T'}(x,j,i) + nadj_{T'}(j,i,x) counts paths
of T' containing (x,j,i) or (j,i,x).

The sum over x of s_x-weighted nadj terms encodes 3-cycle information,
since s_x = -1 iff {i,j,x} forms a 3-cycle using arc i->j.

---

## Verification

| n | Trials | matched_T = matched_{T'} | #U_T - #U_{T'} = delta |
|---|--------|--------------------------|----------------------|
| 5 | 20 | all equal | all match |
| 6 | 10 | all equal | all match |
| 7 | 5 | all equal | all match |
