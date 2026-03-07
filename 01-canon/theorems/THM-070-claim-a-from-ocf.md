# THM-070: Claim A from OCF (Clean Proof Chain)

**Type:** Theorem (PROVED)
**Certainty:** 5 -- PROVED (elementary, 4-step proof)
**Status:** PROVED
**Added by:** opus-2026-03-07-S36
**Tags:** #claim-a #ocf #conflict-graph #clique #independence-polynomial

---

## Statement

**Claim A:** For any tournament T and vertex v: H(T) - H(T-v) = 2 * sum_{C through v} mu(C).

This follows from OCF (H(T) = I(Omega(T), 2)) via a clean 4-step chain.

---

## Proof

### Step 0: OCF (Grinberg-Stanley)
H(T) = I(Omega(T), 2) for all tournaments T.

### Step 1: Through-v Clique Property

**Lemma.** All directed odd cycles through v form a clique in Omega(T).

*Proof.* Any two cycles C1, C2 through v satisfy v in V(C1) cap V(C2). Since adjacency in Omega means sharing a vertex, C1 and C2 are adjacent. QED.

### Step 2: Unique Cycle Decomposition

By Step 1, any independent set S in Omega(T) contains **at most one** cycle through v (since two through-v cycles would be adjacent, contradicting independence).

Therefore:
- I(Omega(T), 2) = [sum over S with NO cycle through v] + [sum over S with EXACTLY ONE cycle through v]
- I(Omega(T), 2) = I(Omega(T-v)*, 2) + sum_{C through v} sum_{S': indep in Omega\N[C]} 2^{|S'|+1}

where Omega(T-v)* is the subgraph of Omega(T) on cycles NOT through v.

**Note:** Omega(T-v)* = Omega(T-v) because:
- A cycle not through v exists in both T and T-v
- Two such cycles conflict (share a vertex) iff they do in Omega(T)

So the first sum = I(Omega(T-v), 2) = H(T-v) by OCF applied to T-v.

### Step 3: Graph Equality (THM-069)

For each C through v:
  sum_{S': indep in Omega\N[C]} 2^{|S'|+1} = 2 * I(Omega\N[C], 2) = f(C) = 2*mu(C)

The equality f(C) = 2*mu(C) follows from the graph equality:
  Omega(T)\N[C] = Omega(T-v)|_{avoid C\{v}}

which holds because both sides have the same vertex set {cycles D : V(D) cap V(C) = empty}
(since v in C implies any cycle disjoint from C also avoids v).

### Step 4: Combining

H(T) = I(Omega(T), 2)
     = H(T-v) + sum_{C through v} 2*mu(C)

Rearranging: H(T) - H(T-v) = 2 * sum_{C through v} mu(C).

QED.

---

## Key Insight

The proof requires NO inclusion-exclusion. The through-v clique property (Step 1) ensures that each independent set in Omega(T) containing a through-v cycle contains EXACTLY ONE such cycle. This allows direct summation rather than IE.

This was already known abstractly (INV-038 noted the clique property), but the consequence that IE is unnecessary was not previously observed.

---

## Dependencies

- OCF: H(T) = I(Omega(T), 2) [Grinberg-Stanley, arXiv:2412.10572]
- THM-069: Graph equality Omega(T)\N[C] = Omega(T-v)|_{avoid C\{v}}
- Claim B (THM-003): Used historically but NOT needed in this proof chain

---

## Significance

This gives the cleanest known proof that OCF implies Claim A. The proof is entirely elementary — no generating functions, no Fourier analysis, no algebraic manipulations beyond the definition of independence polynomial. The key observations are:
1. Through-v cycles are pairwise adjacent (trivial)
2. Removing v from a cycle-disjoint set doesn't affect it (trivial)
3. These two trivial observations plus OCF immediately give Claim A.
