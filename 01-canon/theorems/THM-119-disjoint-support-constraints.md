# THM-119: Disjoint Support Constraint Property at Omega_2

**Status:** PROVED (pure algebraic, no computation needed)
**Author:** kind-pasteur-S45 (2026-03-09)
**Verified:** Exhaustive n=4,5,6; sampled n=7,8,9,10

## Statement

For any tournament T on n vertices, every 2-path (a,b,c) in A_2(T) has AT MOST ONE non-allowed 1-face. Consequently, the constraint matrix defining Omega_2(T) from A_2(T) has full row rank, and:

    dim(Omega_2) = |A_2| - #NA_faces

where #NA_faces = #{ordered pairs (a,c) : c->a in T, and exists b with a->b->c}.

## Proof

A 2-path (a,b,c) means a->b and b->c in the tournament T. Its boundary faces are:

- d_0(a,b,c) = (b,c): In A_1 because b->c. **Always allowed.**
- d_1(a,b,c) = (a,c): In A_1 iff a->c. **Non-allowed iff c->a.**
- d_2(a,b,c) = (a,b): In A_1 because a->b. **Always allowed.**

So each 2-path has exactly 0 or 1 non-allowed faces, depending on whether its "skip edge" (a,c) is forward or backward in the tournament. QED (disjoint support claim)

**Corollary (Full row rank):** The Omega_2 constraint matrix P has one row per non-allowed 1-face (a,c) with c->a. The support of row (a,c) consists of exactly the 2-paths {(a,b,c) : a->b->c} — the paths whose middle face is (a,c). Since each path has at most one non-allowed face, different rows have DISJOINT column supports. Matrices with disjoint non-empty row supports have full row rank (each row has a private nonzero entry that no other row touches). Therefore rank(P) = #rows = #NA_faces.

**Corollary (Omega_2 dimension):** dim(Omega_2) = |A_2| - rank(P) = |A_2| - #NA_faces.

## Significance

1. This property is UNIQUE to level 2 in tournaments. At level p >= 3, a p-path can have up to p-1 non-allowed faces (verified: at level 3, the distribution is 25%/50%/25% for 0/1/2 NA faces).

2. At level 4, constraint rows can overlap (a single path contributes to multiple constraints), leading to rank deficit. This is WHY beta_4 can be nonzero (first at n=8) while beta_2 cannot.

3. This is a NECESSARY condition for beta_2 = 0, but NOT sufficient. Exactness of the chain complex at Omega_2 (ker d_2 = im d_3) requires additional structure.

## Sharpness

The disjoint support property depends on COMPLETENESS of the tournament:
- Removing even one edge from a tournament can create beta_2 > 0 (verified: 13/500 at n=6).
- The property holds because in a tournament, every endpoint pair (a,b) and (b,c) of a 2-path IS an edge, so faces d_0 and d_2 are always allowed.

## Files

- `04-computation/beta2_disjoint_support_proof.py` — verification and analysis
- `04-computation/beta2_completeness_argument.py` — sharpness tests
- `05-knowledge/results/beta2_disjoint_support_proof.out` — computational output
