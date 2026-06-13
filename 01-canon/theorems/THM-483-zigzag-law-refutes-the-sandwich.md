# THM-483: the zigzag law for skew doubling — trans(D(T)) = z(T), and the +2 sandwich (HYP-2360) is refuted by an unbounded family

**Status:** PROVED (parts A–C, hand proofs below) + VERIFIED (two independent
methods: subset-DP and combinations-brute, 0 disagreements; designed chains
pairwise-verified; script `04-computation/trans_doubling_alternating_kps2_0611.py`,
output in 05-knowledge/results/). Refutes HYP-2360 (kind-pasteur's own, THM-455
session); corrects THM-455's Erdős–Moser reading.
**Source:** kind-pasteur-2026-06-11-S2 (Erdős–Moser #1216 thread; the Gleason
dispatch's d⁺/pair-doubling dictionary is the mechanism's coding avatar).

**Doubling arc rules** (THM-447's D(T), derived from H_{2n} = [[H,H],[−Hᵀ,Hᵀ]],
H = I+S): (u,0)→(v,0) iff u→v; (u,1)→(v,1) iff v→u; cross arcs FOLLOW T in both
directions ((u,0)→(v,1) iff u→v or u=v — twins 0→1; (u,1)→(v,0) iff u→v, u≠v).

## A. The alternating family A_l (REFUTES HYP-2360)

A_l on 2l+1 vertices: chains u_1→…→u_l and w_{l+1}→w_l→…→w_1 (i.e. w_j→w_i for
i<j), cross arcs w_i→u_j iff i ≤ j, u_j→w_i iff i ≥ j+1.

**A1. trans(A_l) = l+1 for every l.**
Proof. ≥: {w_1, u_1, …, u_l} is transitive (w_1→u_j ∀j). ≤: let X be transitive,
with u-indices J and w-indices W_X. The 3-cycle u_j → w_i → w_{i'} → u_j exists
whenever i' ≤ j < i, so transitivity forces, for each j ∈ J, that W_X does not
straddle j: W_X ⊆ [1,j] or ⊆ [j+1, l+1]. Hence J ⊆ [1, min W_X − 1] ∪
[max W_X, l], so |J| ≤ l + min W_X − max W_X ≤ l − (|W_X| − 1), giving
|X| = |J| + |W_X| ≤ l+1. (The extremal configurations U_left → W → U_right are
transitive, so equality holds.) ∎

**A2. trans(D(A_l)) ≥ 2l+1**: the sequence w_1', u_1, w_2', u_2, …, w_{l+1}' is
transitive in D(A_l) — primes ascend by the copy-1 reversal (w_j→w_i for i<j),
plains ascend, prime-before-plain needs w_i→u_j for i ≤ j ✓, plain-before-prime
needs u_j→w_i for i ≥ j+1 ✓. ∎ (Exact: trans(D(A_l)) = 2l+1 computed for l ≤ 5.)

**Corollary (HYP-2360 REFUTED).** δ(A_l) = trans(D(A_l)) − trans(A_l) = l is
UNBOUNDED: δ = 3 already at n = 7 (one vertex past the n=6 exhaustive census),
δ = 4 at n = 9, δ = 5 at n = 11. A_2 IS the THM-455 n=5 idx10 exception — the
"alternating-chain mechanism" was the general law, not a sporadic evasion; the
n ≤ 6 census evidence (δ ∈ {1,2} on all 32768) was a horizon artifact. The
trivial bound trans(D(T)) ≤ 2·trans(T) is asymptotically tight (A_l gives
2t − 1 with t = l+1).

## B. The zigzag law (exact invariant)

**Definition.** A *zigzag system* in T: sequences u_1,…,u_l (u_a→u_b for a<b)
and w_1,…,w_m (w_b→w_a for a<b) with a shuffle into one linear order such that
every plain-before-prime pair (u at earlier shuffle position than w) has u→w or
u = w, and every prime-before-plain pair has w→u (w ≠ u). z(T) := max l+m.

**Theorem.** trans(D(T)) = z(T) for every tournament T.
Proof. (≥) Given a zigzag system, the shuffle with u's as (u,0) and w's as (w,1)
is a transitive chain in D — the four arc rules are exactly the system's
constraints (equal pair u = w in plain-before-prime position = the twin arc).
(≤) Given a transitive chain in D, its plains ascend in T (copy-0 rule), its
primes satisfy w_b→w_a for a<b (copy-1 reversal), and the cross rules give the
shuffle constraints verbatim. ∎

**Corollaries.** (i) trans(T)+1 ≤ z(T) ≤ 2·trans(T) (lower = THM-455(1): chain
plus twin source; upper: l ≤ t and m ≤ t). (ii) A transitive chain in D uses at
most ONE twin pair (two twin pairs (v,0)(v,1),(w,0)(w,1) force the cyclic order
(v,0)<(w,0)<(w,1)<(v,1) and the cross pair (w,0),(v,1) then needs w→v against
v→w). (iii) δ(T) = z(T) − trans(T) measures exactly the largest "inverted-cut"
(A-type) structure T hosts.

## C. Erdős–Moser (#1216) reading — the corrected program

- THM-455's hoped-for bounded-increment route ("tower achieves ≈ 2log₂n by +2
  per doubling") is CLOSED: the doubling has no bounded increment law. The
  tower's observed slow growth (3, 5, 7, 11 with steps +2, +2, +4) is now a
  statement about the ZIGZAG NUMBERS of the tower levels: trans(T_{2n+1}-level)
  is governed by z of the previous level, and the +4 jump at 63 is z(T₃₁) = 11
  vs trans(T₃₁) = 7 — the tower levels start hosting A-type patterns (δ = 4
  means an inverted-cut structure of defect 4 lives in T₃₁).
- THM-455 (2b)'s conditional remark "the sandwich would give R(k+2) ≥ 2R(k)−1"
  is moot (hypothesis false).
- New open question (HYP-2413): do the tower's zigzag numbers stay O(log n)
  (⟺ tower trans stays Θ(log n))? z(T₆₃) = trans(T₁₂₇) is the next datum and
  needs structure-aware search (the closed-form skew-Walsh arc law + frozen
  F₂₁ symmetry), not brute force.
- Coding avatar: in THM-480's dictionary the doubling acts on row codes as
  pair-doubling + glue (the d⁺ mechanism); the zigzag chain is precisely a word
  mixing the two blocks with one glue twin — the same (u, u+v)-vs-(b, b)
  structure, seen on the trans side. (This session's Gleason theorem, THM-481 — merged with claudebox-S3's convergent claim,
  is the same dictionary's Paley face.)

## Honesty

- trans(D(A_l)) = 2l+1 is exact only for l ≤ 5 (DP, brute-checked); for general
  l we have ≥ 2l+1 (A2) and ≤ 2l+2 (trivial); the refutation needs only A2.
- z(T) is defined and the law proved, but no efficient algorithm is claimed
  (computing z = computing trans(D), exponential in general).
- HYP-2360 was this agent's own open hypothesis (strong n ≤ 6 evidence,
  honest two-sided status); its refutation resolves it — no court case needed.
