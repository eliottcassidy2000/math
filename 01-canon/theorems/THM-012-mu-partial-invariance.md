# THM-012: Partial Invariance of mu Under Arc Flips

**Type:** Theorem
**Certainty:** 5 — PROVED (elementary + verified n=5,6)
**Status:** PROVED
**Added by:** opus-2026-03-05-S2
**Tags:** #mu #arc-reversal #invariance #open-q-009

---

## Statement

Let T be a tournament on n vertices, v a vertex, C a directed odd cycle through v, and i->j an arc of T not used in C. Let T' be obtained from T by flipping i->j to j->i.

**If at least one of i, j is in V(C)\{v}**, then mu_T(C) = mu_{T'}(C).

Equivalently: mu(C) is invariant under arc flips where at least one endpoint belongs to C\{v} (provided the arc is not part of the cycle).

---

## Proof

mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2).

The "available" vertex set is A = V(T-v) \ C\{v} = V\V(C).

The flip changes arc i<->j in T-v (since i,j != v). There are three cases:

**Case 1: Both i,j in C\{v}.** Then i,j are both excluded from A. The sub-tournament T[A] = T'[A] (unchanged). So Omega(T-v)|_{avoid C\{v}} = Omega(T'-v)|_{avoid C\{v}}, and mu is unchanged. (Note: the arc i->j is "not used in C" by hypothesis, but both endpoints are in C.)

**Case 2: Exactly one of i,j in C\{v}, say i in C\{v} and j not in V(C).** Then i is excluded from A but j is in A. The flip changes arc i<->j. Since i is NOT in A, the sub-tournament T[A] is unaffected (the flip involves a vertex outside A). So T[A] = T'[A], and mu is unchanged.

**Case 3: i not in V(C) and j in C\{v}.** Same argument as Case 2 with roles swapped.

In all three cases, the sub-tournament on A is unchanged, so mu(C) = mu'(C). QED.

---

## Non-Invariance When Both Endpoints Outside V(C)

If BOTH i and j are in A (i.e., both outside V(C)), the flip directly changes T[A], which can change the odd cycle structure of A, changing mu.

At n=6 with 3-cycles (|A| = 3): flipping an arc among the 3 available vertices can create or destroy a 3-cycle among them, changing mu between 1 and 3. Verified: 2286 changes out of 5706 cases.

---

## Verification Record

| n | One endpoint in C\{v} | Both outside V(C) |
|---|----------------------|-------------------|
| 5 | 0/30720 changes | 0/7680 changes (trivial: all mu=1) |
| 6 | 0/18932 changes | 2286/5706 changes |

---

## Verification Scripts

- `04-computation/mu_arc_invariance.py` — partial invariance under arc flips verification
- `04-computation/arc_reversal_q009.py` — arc-reversal invariance of μ tests

## Significance for OPEN-Q-009

For the arc-reversal invariance approach to proving Claim A, we need:

D(T,v) = D(T',v) where D(T,v) = H(T) - H(T-v) - 2*sum_{C through v} mu(C).

This theorem shows that for cycles C not using arc i->j, the change in mu(C) depends on whether i,j are both outside V(C). The difficult part of the arc-reversal argument is:
1. Tracking cycles that USE the flipped arc (cycle creation/destruction)
2. Tracking mu changes for cycles where both flipped vertices are in the available set

The mu changes in case (2) are exactly determined by whether the flip creates/destroys a 3-cycle among the available vertices. This is a tractable combinatorial problem.
