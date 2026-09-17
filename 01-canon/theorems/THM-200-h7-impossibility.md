> **CORRECTION 2026-09-17:** The former gap-tower corollary below is refuted.
> `H=63` is realized by THM-344; adjoining a source preserves H and makes
> the tournament non-strong. The repaired all-order exclusions `{7,21}` are
> [THM-1370-h-spectrum-omits-7-21-all-n](THM-1370-h-spectrum-omits-7-21-all-n.md).
> The H=7 statement below is retained; the withdrawn corollary is historical.

# THM-200: H=7 Impossibility Theorem

**Status:** PROVED
**Proved by:** opus-2026-03-14-S71g
**Dependencies:** OCF (THM-002), directed cycle enumeration

## Statement

For any tournament T on n ≥ 1 vertices, H(T) ≠ 7.

## Proof

By the OCF (THM-002), H(T) = I(Ω(T), 2) where Ω(T) is the directed-odd-cycle conflict graph.

**Step 1.** I(G, 2) = 7 if and only if G = K₃.

*Proof:* I(G, 2) = Σ_k α_k · 2^k = 1 + 2α₁ + 4α₂ + ···. For I = 7: 2α₁ + 4α₂ + ··· = 6, so α₁ = 3 and α_k = 0 for k ≥ 2. This means G has 3 vertices, all pairs adjacent. So G = K₃. Conversely, I(K₃, 2) = 1 + 6 = 7.

**Step 2.** Three pairwise-sharing directed 3-cycles span at most 7 vertices.

*Proof:* Each cycle uses 3 vertices. If all share a common vertex, |V(C₁) ∪ V(C₂) ∪ V(C₃)| ≤ 1 + 2·3 = 7. If no common vertex, the maximum is 6 (each pair shares exactly one distinct vertex).

**Step 3.** No tournament on n ≤ 6 vertices has exactly 3 directed odd cycles.

*Proof:* Verified by exhaustive enumeration of all tournaments at n = 3 (8 total), n = 4 (64), n = 5 (1024), n = 6 (32768). Total = 0 in all cases.

Key sub-facts:
- n = 3, 4: t₃ = 3 is impossible (t₃ ∈ {0,1} at n=3; t₃ ∈ {0,1,2} at n=4).
- n = 5: t₃ = 3 always forces t₅ ≥ 1 (Lemma: on 5 vertices, 3 triangles create a 5-cycle).
- n = 6: t₃ = 3 with t₅ = 0 never occurs (exhaustive).

**Step 4.** At n = 7, every tournament with exactly 3 directed odd cycles has Ω ≠ K₃.

*Proof:* Exhaustive enumeration of all 2^21 = 2,097,152 tournaments on 7 vertices. Found 3360 tournaments with exactly 3 directed odd cycles (all 3-cycles, t₅ = t₇ = 0). In ALL 3360 cases, at least one pair of triangles is vertex-disjoint. The conflict graph is always P₃ (path on 3 vertices), giving I(P₃, 2) = 15 ≠ 7.

**Step 5.** For n ≥ 8, H(T) ≠ 7.

*Proof:* By Step 1, H = 7 requires Ω = K₃, meaning exactly 3 directed odd cycles, all pairwise sharing a vertex. These 3 cycles are 3-cycles (since a 5-cycle on 5 vertices forces ≥3 additional 3-cycles, giving total ≥ 4). Let S = V(C₁) ∪ V(C₂) ∪ V(C₃) with |S| ≤ 7 (Step 2).

Since H(T) = 7 requires total directed odd cycles = 3, and all 3 are within S, the subtournament T|_S also has these 3 cycles. Any additional cycle on T|_S would also be a cycle of T, contributing to the total. So T|_S has at most 3 directed odd cycles (exactly 3, since the 3 planted triangles exist on S).

- If |S| ≤ 6: By Step 3, T|_S cannot have exactly 3 directed odd cycles. Contradiction.
- If |S| = 7: By Step 4, T|_S has exactly 3 directed odd cycles but Ω(T|_S) ≠ K₃ (at least one disjoint pair). Since the 3 cycles of T are exactly the 3 cycles of T|_S, the conflict structure is the same. So Ω(T) ≠ K₃.

In both cases, I(Ω(T), 2) ≠ 7, so H(T) ≠ 7.

**Completeness check:** If total directed odd cycles ≠ 3, then |V(Ω)| ≠ 3, so I(Ω, 2) ≠ 7 (since I = 7 requires exactly 3 vertices). This covers all remaining cases.

∎

## Computational verification

- n ≤ 6: H ≠ 7 verified exhaustively (33,864 tournaments).
- n = 7: H ≠ 7 verified exhaustively (2,097,152 tournaments).
- n = 8: H ≠ 7 verified by sampling (500,000 random tournaments, 0 hits).
- All these match the theorem prediction.

## Withdrawn corollary and repaired consequence

**REFUTED historical claim:** all `7*3^k` were asserted absent for non-strong
tournaments. The SCC product formula does not require a factor7: one strong
component can already have H=63 or189. In particular, take an eight-vertex
H=63 witness from
[THM-344-n8-h63-complete-omega](THM-344-n8-h63-complete-omega.md) and adjoin
a source vertex. Every Hamiltonian path starts at that source and then is
a path of the witness, so H remains63 in a non-strong tournament.

**Survivor:** H=7 and H=21 are excluded at every order by
[THM-1370-h-spectrum-omits-7-21-all-n](THM-1370-h-spectrum-omits-7-21-all-n.md).
Its later proof supersedes the old claim that H=21 had only bounded support.
The missing coordinate was the attainable value of an entire SCC, rather
than the prime factorization of its Hamiltonian count. The elementary first
failed implication is `7 divides a product => some factor equals7`.

## Key scripts

- `04-computation/h7_theorem.py` — Main proof verification
- `04-computation/h7_clean_proof.py` — Supporting lemmas
- `04-computation/h7_n7_check.py` — n=7 exhaustive check
- `05-knowledge/results/h7_theorem.out` — Computation output
