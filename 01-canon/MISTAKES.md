# Mistakes Log

**Purpose:** Every error that has been made and corrected — with enough context that no Claude instance ever repeats it. Read this before doing any computational or proof work.

Format per entry:
- What was assumed / done
- Why it was wrong
- The correct framing
- Impact on existing results
- Source (who found it, when)

---

## MISTAKE-001: μ Computation Bug in Scripts 6-9

**Date discovered:** ~2026-02-26 (pre-Cowork sessions); logged 2026-03-05
**Found by:** Claude instance (Account unknown, pre-Cowork era), reported in MASTER_FINDINGS.md
**Affects:** Scripts 6-9 (sum_mu() function); does NOT affect scripts 1-5

### What was assumed
`ind_poly_at_2_restricted()` correctly computes μ(C) = I(Ω(T−v)|_{avoid C\{v}}, 2) by iterating over cycles using `T[perm[i]][perm[(i+1)%len]]`.

### Why it was wrong
`perm` is a permutation of T−v vertices that does NOT include v. The cycle-checking code uses the **full T matrix** instead of T−v. Specifically, for a cycle of length L:
- `T[perm[i]][perm[(i+1)%len]]` checks arcs in the full tournament T
- But cycles of T−v (cycles not using v) should be checked in T−v
- The wrap-around `perm[(i+1)%len]` may also have indexing issues for general L

### The correct framing
The restriction to T−v must be done BEFORE checking cycles. Any cycle-finding in `ind_poly_at_2_restricted()` must operate on T−v's adjacency matrix exclusively.

### Impact assessment
- **Scripts 1-5 (inshat analysis):** UNAFFECTED (do not use μ computation)
- **Paper's Claim A verification:** STATED UNAFFECTED by MASTER_FINDINGS (verification runs used a separate code path)
- **Paper's Claim B verification:** Needs confirmation — Claim B involves I(Ω,2) directly
- **The verification table in the paper (0 failures for Claim A at n≤6):** Stated to be valid

⚠️ **UNRESOLVED:** It is not fully clear which code path the paper used for Claim A verification at n=6. If sum_mu() was used, the 0-failure count may not be reliable. This is the subject of DISC-001 in 02-court/active/.

### Lesson
When computing μ(C) for any cycle C in T:
1. First restrict to T−v (remove v and all incident arcs from the adjacency matrix)
2. Then find all odd cycles in T−v that are vertex-disjoint from C\{v}
3. Build Ω(T−v)|_{avoid C\{v}} from these cycles
4. Evaluate I(·, 2) on this restricted conflict graph

Never use the full T matrix when computing anything about T−v.

### Resolution
**Independent verification completed** (opus-2026-03-05-S1): tournament_lib.py implements mu(C) correctly per the above steps. Exhaustive verification of Claim A at n<=6 (196,608 pairs, 0 failures) confirms the paper's results are valid. See DISC-001 for the formal resolution.

---

## MISTAKE-004: OCF is Recursive, Not a Closed-Form Over All Odd Cycles

**Date discovered:** During file.txt exploration (pre-2026-03-05)
**Found by:** Claude instance (Account unknown), via constructing a counterexample
**Affects:** Any attempt to compute H(T) as I(Omega_ALL(T), 2) where Omega_ALL = all odd cycles

### What was assumed
H(T) = I(Omega(T), 2) where Omega(T) is the set of ALL directed odd cycles of T, and I is the independence polynomial evaluated at 2.

### Why it was wrong
Counterexample: T on {1,2,3,4,5,6} with 3-cycle (1->2->3->1), 3-cycle (4->5->6->4), and all arcs from {1,2,3} to {4,5,6}. The only odd cycles are C1 and C2 (no cross-group cycles since all arcs go one way). I(Omega_ALL, 2) = 1 + 2*3 + 2*3 + 4*3*3 = 49. But H(T) = 3*3 = 9 (3 orderings of each group, concatenated).

### The correct framing
OCF is a RECURSIVE formula: H(T) = H(T-v) + 2 * sum_{C through v} H(T[V\V(C)]). The "closed form" H(T) = I(Omega(T), 2) holds ONLY when Omega(T) is defined relative to the recursive vertex-removal order, NOT as the set of all odd cycles. The mu weights mu(C) = H(T[V\V(C)]) are themselves computed recursively.

### Impact
Any proof strategy that tries to express H(T) directly as a sum over all odd-cycle collections with simple weights is using the wrong formula. The recursion is essential.

### Lesson
Always use the recursive formulation Claim A: H(T) - H(T-v) = 2 * sum_{C through v} mu(C), where mu(C) = H(T[V\V(C)]). Do not flatten this into a non-recursive independence polynomial over all cycles.

---

## MISTAKE-005: Cycle Bijection Under Arc Reversal Fails

**Date discovered:** During file.txt exploration
**Found by:** Claude instance (Account unknown)
**Affects:** Proof strategy for arc-reversal invariance D(T,v) = D(T',v)

### What was assumed
When T' is obtained from T by flipping arc i->j to j->i (with i,j != v), odd cycles through v containing i->j in T biject with odd cycles through v containing j->i in T', preserving V(C).

### Why it was wrong
A cycle C through v containing i->j can be written as v ->^{P1} i -> j ->^{P2} v. The "conjugate" C' would need j->i in T', requiring the path segments to be REVERSED. But reversing a directed path P2 does not generally give a valid directed path in the same tournament (arcs may point the wrong way).

### The correct framing
There is NO individual bijection C <-> C' preserving vertex sets. The arc-reversal invariance, if true, must hold as a SUM equality: sum_{C: i->j in C} H(T[V\V(C)]) = sum_{C': j->i in C'} H(T'[V\V(C')]). This is a weaker statement that doesn't require a cycle-by-cycle bijection.

### Impact
The arc-reversal invariance D(T,v) = D(T',v) remains the key unproved step for a general proof of Claim A. The cycle bijection approach is a dead end; a sum-level argument is needed.

---

## MISTAKE-002: Exact Path Formula H(T) = B_v + S_v + R_v

**Date discovered:** During script verification at n=4
**Affects:** Any approach trying to decompose H(T) exactly as B_v + S_v + R_v

### What was assumed
H(T) = B_v + S_v + R_v (exact equality)

### Why it was wrong
Verified FALSE: 96 failures out of 256 pairs at n=4.

### The correct framing
Only the PARITY version holds: H(T) ≡ B_v + S_v + R_v (mod 2). The exact identity is wrong.

### Impact
Similarly, S_v + R_v = 2Σμ(C) is FALSE (144/256 failures at n=4). These exact decompositions are not the right approach.

---

## MISTAKE-003: Per-Path Identity as a Path to Claim A for n≥6

**Date discovered:** During n=6 verification
**Affects:** Any proof strategy that relies on the per-path identity to prove Claim A for general n

### What was assumed
The per-path identity (inshat−1)/2 = Σ_{3-cycle embeddings} μ(C) could serve as an intermediate step toward proving Claim A for all n.

### Why it was wrong
The per-path identity fails for 2,758/9,126 ≈ 30% of (T,v,P') triples at n=6. It is provably a 3-cycle-only formula (by THM-005). Longer cycles contribute to Claim A's RHS but are invisible to the per-path identity.

### The correct framing
The per-path identity is valid and useful for n≤5 (where it implies Claim A). For n≥6, a different per-path formula is needed — one that accounts for contributions from all odd cycles, not just 3-cycles. See OPEN-Q-004.

### Impact
Any proof of Claim A for n≥6 must go beyond the per-path identity. The five proof strategies in the paper are the current best alternatives.
