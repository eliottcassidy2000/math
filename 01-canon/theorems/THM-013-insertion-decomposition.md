# THM-012b: Insertion Decomposition Formula

**Type:** Theorem
**Certainty:** 4 — VERIFIED (exhaustive n≤5, formula holds with 0 failures)
**Status:** VERIFIED (proof pending full rigor)
**Contributed by:** opus-2026-03-05-S3
**NOTE:** This was mislabeled as THM-012 in the original header. The filename
THM-013-insertion-decomposition.md is also misleading. The actual THM-012 is
THM-012-mu-partial-invariance.md and the actual THM-013 is THM-013-arc-flip-delta-I.md.
This file is THM-012b (a standalone result, not referenced by other theorems under this ID).
**Tags:** #insertion #orphan-paths #hamiltonian-paths #claim-a-reformulation

---

## Statement

For every tournament T on n vertices and every vertex v:

```
H(T) - H(T-v) = Σ_{P'∈Ham(T-v)} #{TypeII(P')} + #{orphan Ham paths of T w.r.t. v}
```

where:
- **TypeII(P')** = #{j : s[j]=1 and s[j+1]=0 in the signature of P' w.r.t. v}
  (equivalently, by THM-004: #{TypeII(P')} = (inshat(v,P') - 1)/2)
- **Orphan Ham path** of T w.r.t. v: a directed Hamiltonian path P ∈ Ham(T) such that
  removing v from P gives an invalid directed path in T-v
  (i.e., the two neighbors of v in P have a backward arc: if P[k-1]=a and P[k+1]=b,
  then T(a,v)=1, T(v,b)=1, but T(b,a)=1, so the sequence a→b is backward in T-v)

Equivalently, since #{TypeII(P')} = (inshat(v,P')-1)/2:

```
H(T) - H(T-v) = Σ_{P'∈Ham(T-v)} (inshat(v,P') - 1)/2 + #{orphan Ham paths}
```

---

## Derivation

**Step 1 (from insertion framework):** Every Ham path P ∈ Ham(T) arises in exactly one of two ways:
- *Non-orphan*: removing v from P gives a valid P' ∈ Ham(T-v). This means v is inserted between
  P[k-1] and P[k+1] with T(P[k-1], P[k+1]) = 1 (a valid insertion).
- *Orphan*: removing v gives an invalid sequence (T(P[k-1], P[k+1]) = 0).

So: H(T) = Σ_{P'∈Ham(T-v)} insact(v,P') + #{orphans}.

**Step 2:** insact(v,P') = (valid insertions of v into P') = #{valid TypeI interior positions} + #{valid boundary positions}. And insact(v,P') - 1 = #{TypeII(P')} (derived from inshat = insact + TypeII and inshat = 1 + 2*TypeII).

**Step 3:** H(T) - H(T-v) = Σ(insact-1) + #{orphans} = Σ TypeII + #{orphans}. □

---

## Key Corollaries

**Cor 1 (Orphan structure):** Every orphan path P is associated with exactly one directed
3-cycle through v: if v is between a and b in P with T(a,v)=1, T(v,b)=1, T(b,a)=1, then
a→v→b→a is a directed 3-cycle. Hence #{orphans} = Σ_{3-cycles C through v} N_Orphan(C).

**Cor 2 (Reformulation of Claim A):**
Claim A (H(T)-H(T-v) = 2*Σ mu(C)) is equivalent to:
```
Σ_{P'∈Ham(T-v)} #{TypeII(P')} + #{orphan Ham paths} = 2*Σ_{odd C through v} mu(C)
```
Since Σ TypeII = Σ_{3-cycles C} N_TypeII(C) (by THM-005) and #{orphans} = Σ_{3-cycles C} N_Orphan(C):
```
Σ_{3-cycles C through v} [N_TypeII(C) + N_Orphan(C)] = 2*Σ_{ALL odd C through v} mu(C)
```
This shows: the LHS is determined entirely by 3-cycles (via TypeII positions and orphan paths),
while the RHS includes ALL odd cycles. For n≤5 (only 3-cycles), this becomes:
Σ_{3-cycles C} [N_TypeII(C) + N_Orphan(C)] = 2*Σ_{3-cycles} mu(C).

---

## Verification Record

| n | Pairs (T,v) | Formula failures |
|---|-------------|-----------------|
| 3 | 24 | 0 ✓ |
| 4 | 256 | 0 ✓ |
| 5 | 5,120 | 0 ✓ |

Verified by `04-computation/investigate_inshat.py` (opus-2026-03-05-S3).

---

## Relationship to Paper's Remark (LINE 1261 ERROR)

The paper's Remark at line 1261-1264 claims: "Claim A is equivalent to Σ(inshat-1)/2 = Σ mu(C)."

**This is INCORRECT.** See MISTAKE-006. The correct equivalence requires including the orphan term:
Σ(inshat-1)/2 + #{orphans} = 2*Σ mu(C), not Σ(inshat-1)/2 = Σ mu(C).

The paper's remark implicitly assumes H(T) = Σ inshat(v,P'), which fails whenever there are
orphan paths (96/256 pairs at n=4, 3080/5120 at n=5, while Claim A has 0 failures).

---

## Failed Per-Cycle Sub-Identity

*For 3-cycles at n=4:* The per-cycle identity N_TypeII(C) + N_Orphan(C) = 2*mu(C) holds
exactly (0 failures out of 192 triples). However, this FAILS at n=5 (277/693 failures)
even though mu(3-cycle)=1 at n≤5. The total SUM is correct but individual cycles need not
satisfy the per-cycle identity. This dead end is documented in TANGENTS.
