# Inbox Processing Report

**Generated:** 2026-03-05 13:41:00  
**Files found:** 2  

---

## `FINAL_FINDINGS.md`

*Size: 7,533 bytes | Extension: .md*

**STATUS: NEW CONTRIBUTION 🆕**

**Content preview:**

```
# Final Findings: The Per-Path Identity, Its Scope, and Its Failure

## The Core Algebraic Fact (proven, unconditional)

For any path P' of Ham(T-v), the quantity (inshat(v,P')-1)/2 equals
the number of Type-II positions in P' -- this is a pure algebraic identity
that holds for any binary string.

## What the Per-Path Identity Actually Says

The per-path identity from the paper states:
  (inshat(v,P')-1)/2 = sum_{Type-II at j} mu(v, P'[j], P'[j+1])

Since LHS = #Type-II positions (algebraic), this is really asking:
  #Type-II positions = sum_{Type-II at j} mu(v, P'[j], P'[j+1])

This can only hold if the weighted sum of mu values equals the count,
which requires careful cancellation.

## Why the Identity Holds for n <= 5 (TRIVIALLY)

For a 3-cycle C = (v, a, b), mu(C) = I(Omega(T-v)|avoid{a,b}, 2).
Omega(T-v)|avoid{a,b} contains only odd cycles of T-v on vertices V\{v,a,b}.
That set V\{v,a,b} has n-3 vertices.

  n=3: V\{v,a,b} has 0 vertices. No cycles possible. mu = 1.
  n=4: V\{v,a,b} has 1 vertex.  No cycles possible. mu = 1.
  n=5: V\{v,a,b} has 2 vertices. No 3-cycle (needs 3 vertices). mu = 1.

So for n <= 5: mu(C) = 1 for ALL 3-cycles C through v.
The identity becomes #TypeII = #TypeII. TRIVIALLY TRUE.
The per-path identity for n<=5 has zero mathematical content beyond the algebraic identity.

## Why the Identity Fails for n = 6 (EXACTLY CHARACTERIZED)

At n=6: V\{v,a,b} has 3 vertices, which CAN form a directed 3-cycle.
If they do: mu(v,a,b) = 3. If not: mu(v,a,b) = 1.

THEOREM (proven computationally, clear theoretically):
The per-path identity holds for path P' if and only if:

... (136 more lines)
```

**Suggested integration:**
- Contains theorem/proof material → extract to `01-canon/theorems/`
- Contains error/bug report → add to `01-canon/MISTAKES.md`
- Contains speculative ideas → add to `00-navigation/TANGENTS.md`

### Claude action required:
- Read the full file (it's been archived; path shown below)
- Extract theorems → `01-canon/theorems/`
- Extract open questions → `00-navigation/OPEN-QUESTIONS.md`
- Extract tangents → `00-navigation/TANGENTS.md`
- Extract mistakes → `01-canon/MISTAKES.md`
- Note any discrepancies with existing canon → open court case if needed
- Update `00-navigation/SESSION-LOG.md`

*Archived to: `inbox/processed/2026-03-05/new/FINAL_FINDINGS.md`*

---

## `file.txt`

*Size: 60,929 bytes | Extension: .txt*

**STATUS: NEW CONTRIBUTION 🆕**

**Content preview:**

```
Looking at these documents, I want to identify what's actually proven versus what's conjectured, then push on genuinely new directions.

---

**What's solid so far:**

The arc-reversal formula H(T) - H(T') = adj_T(i→j) - adj_{T'}(j→i) is clean and provable. The invariance D(T,v) = D(T',v) for flips not involving v is the key structural claim. Let me think carefully about whether the bijection argument actually closes.

---

**The cycle bijection: does it actually work?**

Given flip i→j → j→i with i,j ≠ v. Claim: odd cycles through v containing i→j in T biject with odd cycles through v containing j→i in T', preserving V(C).

A cycle C through v containing i→j looks like: v → P₁ → i → j → P₂ → v, where P₁ is a directed path from v to i (not using arc (i,j)), P₂ is a directed path from j to v (not using arc (i,j)).

The proposed bijection: C ↔ C' where C' = v → P₂⁻¹ → i → ... wait. The reversal sends C to: use the same vertex set but pass through j→i instead. We need C' = ... → j → i → ... as a directed cycle through v. 

The natural candidate: C' traverses the same vertices in the order v → P₂' → j → i → P₁' → v, where P₂' is the path from v to j using the same vertices as P₂ (but now going v to j rather than j to v), and P₁' is the path from i to v using vertices of P₁.

**This requires P₂ to be reversible as a directed path.** But reversing a directed path generally doesn't give a valid directed path in the same tournament! The arcs could go the wrong way.

So the bijection C ↔ C' does NOT obviously work. The vertex sets of C and C' would be the same ({v,i,j} ∪ V(P₁) ∪ V(P₂)), but C' is only valid if the reversed paths are directed paths in T'.

**This is a gap.** Let me think about whether it can be patched.

---

**Alternative: count via "cycle indicator" polynomials.**

Rather than bijecting cycles, count the change in 2∑μ directly.

2∑_{C∋v, T} μ_T(C) - 2∑_{C∋v, T'} μ_{T'}(C) = 2∑_{C∋v: i→j∈C} H(T[V\V(C)]) - 2∑_{C'∋v: j→i∈C'} H(T'[V\V(C')])

Since T and T' agree outside arc (i,j): T[V\V(C)] = T'[V\V(C)] for any cycle C (since V(C) contains both i and j). And T'[V\V(C')] = T[V\V(C')] similarly.

So this equals 2∑_{C∋v: i→j∈C} H(T[V\V(C)]) - 2∑_{C'∋v: j→i∈C'} H(T[V\V(C')]).

Now there IS a natural bijection between:
- Odd cycles through v containing i→j in T

... (828 more lines)
```

**Suggested integration:**
- Contains theorem/proof material → extract to `01-canon/theorems/`
- Contains error/bug report → add to `01-canon/MISTAKES.md`
- Contains open questions → add to `00-navigation/OPEN-QUESTIONS.md`

### Claude action required:
- Read the full file (it's been archived; path shown below)
- Extract theorems → `01-canon/theorems/`
- Extract open questions → `00-navigation/OPEN-QUESTIONS.md`
- Extract tangents → `00-navigation/TANGENTS.md`
- Extract mistakes → `01-canon/MISTAKES.md`
- Note any discrepancies with existing canon → open court case if needed
- Update `00-navigation/SESSION-LOG.md`

*Archived to: `inbox/processed/2026-03-05/new/file.txt`*

---

## Summary

Processed 2 file(s). See sections above for required actions.

When done integrating, add an entry to `00-navigation/SESSION-LOG.md`.
