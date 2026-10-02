---
id: THM-4525
title: "The edge multiset dimension of the 6-cube is 15 (Allikvere's Open Problem 1). No edge-multiset resolving set of Q_6 has size <= 14; the resolving 15-sets form exactly 229 Aut(Q_6)-orbits, all with trivial stabilizer. edim_m(Q_d) grows superpolynomially, at least exp((0.6215 - o(1)) d^(1/3)), and explicit sets give edim_m(Q_7) <= 19, Q_8 <= 26, Q_9 <= 38, Q_10 <= 48, Q_11 <= 65, Q_12 <= 76, so density-1/2 landmark sets are far from optimal (Open Problem 4)"
status: >
  FINITE-EXACT + INDEPENDENTLY AUDITED: edim_m(Q_6) = 15. Three
  independent exhaustive searches with different symmetry reductions agree:
  the lane's methods A and B, and the orchestrator's own third search. The
  229 orbits were found by A and B and reproduced by the orchestrator's
  search.
  PROVED + audited:
  - L1, a resolving set has trivial stabilizer;
  - L2, antipodal reversal;
  - L3, the Walsh alternating sum;
  - L4, the counting bound edim_m(Q_6) >= 7 (the paper had 6);
  - L5, the entropy bound edim_m(Q_d) >= exp((c - o(1)) d^(1/3)) with
    c = (ln 2/sqrt 2)^(2/3) = 0.6215, so the growth is superpolynomial.
  VERIFIED + audited: the explicit resolving sets for Q_7..Q_12, sizes 19,
  26, 38, 48, 65, 76; the paper had 63, 115, 246, 492 for d = 7..10.
  VERIFIED (lane-only, double precision with a safety margin, not
  interval-certified): sparse random union bounds, e.g. Q_11 <= 511,
  Q_16 <= 1056, Q_32 <= 6638, with ln M_d / d^(1/3) in [2.76, 2.80] for
  d = 11..32.
  FINITE-EXACT (lane-only):
  - at k = 14 the minimum defect is 1, attained by exactly one orbit
    (0x000001810690226d; the orchestrator verified its defect), and 16
    orbits have defect 2;
  - no union of tournament isomorphism classes, H-levels or score classes
    of the n = 5 tiling cube resolves.
  OPEN: HYP-9169 (ln edim_m(Q_d) = Theta(d^(1/3)), the precise form of
  Open Problem 2); edim_m(Q_7) in [8, 19]; Open Problem 3.
  CITED: Allikvere, arXiv:2608.09983v1 (definitions; the projection lemma;
  Thm 6, infinite iff 2 <= d <= 5; Table 1; Lemma 11, the forest lemma);
  THM-474 (Q_6 is the n = 5 tiling cube).
  Definitions: for an edge uv and a vertex s, d(uv, s) = min(d(u,s), d(v,s)).
  S subset V(Q_d) is edge-multiset resolving if the histograms
  H_e(r) = #{s in S : d(e,s) = r} are pairwise distinct over all edges e;
  edim_m(Q_d) is the least such |S|.
source: collatz-procgen-20260922 session, edim lane (2026-10-01), attacking the open problems of arXiv:2608.09983 at the owner's request (2026-10-01 prompt; Q_6 is the n = 5 tiling cube); audited and promoted by the session orchestrator 2026-10-01
depends_on:
  - external: J. Allikvere, The edge multiset dimension of hypercubes, arXiv:2608.09983v1 (2026-08-05)
related:
  - 05-knowledge/hypotheses/HYP-9169-edge-multiset-dimension-of-hypercubes-grows-like-exp-cube-root.md
  - 01-canon/theorems/THM-474-tilings-are-switching-classes.md
  - 01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md (the same owner prompt)
note: 05-knowledge/results/procgen_edim_20261001_edge_multiset_dimension.md
data: 05-knowledge/results/procgen_edim_20261001_q6_resolving15_orbits.txt (sha256 51d8cf31dfca64f04854aebc4ccd2ae3d485df8afae9b8da1ccca6abee9aa717)
scripts: 04-computation/experiments/procgen_edim_20261001_{run,lib}.py, procgen_edim_20261001_{q5orbits,searchA,searchB,canon6,sad}.c
script_audit: 04-computation/experiments/procgen_edim_20261001_orchestrator_check.{py,c}
output: 05-knowledge/results/procgen_edim_20261001.out
output_sha256: 92036688f2b9910f991bbd4fab08915cfbb8408a7616c8e9c2684ddcf0246704
output_audit: 05-knowledge/results/procgen_edim_20261001_orchestrator_check.out
output_audit_sha256: 2036e50deb14a9527f7ca49ee4c5922550a6f4ef518e6b8b15e5ef498f87b762
hash_basis: raw bytes
audit: >
  The orchestrator read the proofs of L1 (an automorphism fixing S fixes
  every edge histogram, hence every edge, hence every vertex), L2, L4 and
  L5. For L5: the injectivity of e -> H_e, subadditivity, the mean
  m C(d-1,r)/2^(d-1), the geometric-law entropy maximum, and the
  asymptotic count of about sqrt(2 d ln m) heavy levels all check out.
  It also read the WLOG reductions of methods A and B.
  Independent code (procgen_edim_20261001_orchestrator_check.{py,c}; the
  lane's code was not read, and only its data, the explicit sets and the
  orbit list, was taken):
  - A third exhaustive search with its own reduction: split by one fixed
    coordinate, |A| >= |B|, A the minimum of its Aut(Q_5)-orbit (the
    orchestrator's own orbit representatives, whose counts match its own
    Burnside computation for a <= 15), and every B enumerated, with no
    imbalance normal form. For k <= 14 it enumerates 14,168,149,784 leaves, every count equal to sum_a #reps(a) C(32, k-a), and finds no resolving set. At k = 15 it finds 1678 resolving leaves (26,655,385,802 leaves in all), which reduce to exactly the deposited 229 orbits. The search took 875 s with peak RSS 235 MB.
  - All 229 deposited representatives are resolving, are orbit minima
    under the orchestrator's 46080 automorphisms, are pairwise
    inequivalent and have trivial stabilizer. The paper's set lies in one
    of them.
  - The k = 14 near miss has defect exactly 1. Its one collision is
    {24,26} ~ {37,39}, with palindromic histogram (1,1,5,5,1,1).
  - The explicit Q_7..Q_12 sets resolve all 448..24576 edges.
  - L4 holds exactly for m <= 6, and the L5 table (4, 5, 5, 6, 8, 11, 15,
    28, 81, 275, 1159, 50116 for d = 6..1024) was reproduced.
  - L2 was checked on random sets.
  The lane's runner was re-run in --quick mode (k <= 12; 100 checks, ALL
  CHECKS PASSED, 98 s). Of its 100 [OK] lines, 99 occur verbatim in the
  deposited certification output; the other differs only in the
  quick-mode random-set sizes. The 35-minute certification run was not
  repeated, because the orchestrator's own exhaustive search covers it.
  Scope notes:
  - The sparse union bounds and the tournament-structured negatives were
    not re-derived.
  - Minimality of the defect-1 orbit is lane-only.
  - No literature search was made beyond the paper. As of the paper's v1,
    Open Problem 1 was open.
---

# THM-4525 — the edge multiset dimension of Q_6 is 15

**FINITE-EXACT + INDEPENDENTLY AUDITED.** Full note: [procgen_edim_20261001_edge_multiset_dimension](../../05-knowledge/results/procgen_edim_20261001_edge_multiset_dimension.md).

## 1. The problem

Allikvere (arXiv:2608.09983, Aug 2026) proved that `edim_m(Q_d)` is infinite exactly for `2 <= d <= 5`. The paper's Open Problem 1 asks for the exact value at the first finite dimension, where the known bounds were `6 <= edim_m(Q_6) <= 15`. The paper's 96 annealing restarts at size 14 found nothing, but it notes that settling the value "would require an exhaustive or SAT-based lower-bound computation".

The owner pointed the session at this paper together with the selfie-tournament prompt. By THM-474, `Q_6` is the `n = 5` tiling cube (6 tiles).

## 2. The answer

**`edim_m(Q_6) = 15`.** No set of at most 14 vertices gives all 192 edges distinct distance-multisets.

Multiset representations do not refine as landmarks are added, because two histograms can merge. So no splitting-type pruning is valid, and the proof is a complete enumeration up to symmetry. The space is cut by `Aut(Q_6)` (order 46080) through a layer split `S = A x {0} + B x {1}`, with `A` a canonical `Aut(Q_5)`-orbit representative.

Three searches with different normal forms all find nothing at `k <= 14`:
- maximal imbalance;
- minimal imbalance;
- the orchestrator's plain `|A| >= |B|`.

At `k = 15` they all find the same 229 orbits. Every orbit has trivial stabilizer, as it must (L1), so there are `229 * 46080 = 10,552,320` resolving 15-sets. The paper's set is one of them.

The near miss at `k = 14` fails on a single antipodal pair of edges whose common histogram `(1,1,5,5,1,1)` is a palindrome. This is exactly the self-paired collision that L2 allows.

## 3. Growth, and the other open problems

**Open Problem 2 (growth).** An entropy argument (L5) gives `edim_m(Q_d) >= exp((0.6215 - o(1)) d^(1/3))`. Hence the growth is superpolynomial in `d`, although far below the `2^(d-1)` scale of density-1/2 sets. HYP-9169 conjectures `ln edim_m(Q_d) = Theta(d^(1/3))`.

**Open Problem 4 (density 1/2).** It is far from optimal:
- annealed explicit sets give `Q_7 <= 19` (the paper had 63), `Q_8 <= 26` (115), `Q_9 <= 38` (246) and `Q_10 <= 48` (492);
- sparse random union bounds give, for example, `Q_32 <= 6638`, against about `2^31` at density 1/2.

So `8 <= edim_m(Q_7) <= 19`.

**Tournament structure.** It does not help: resolving sets have trivial stabilizer, so no union of tournament classes in the tiling cube resolves.

**FORMALIZED 2026-10-01 (Lean 4.30 core; package [`04-computation/lean/ProcgenSelfieEdim/`](../../04-computation/lean/ProcgenSelfieEdim/README.md); orchestrator-audited).**
- The paper's 15-set resolves Q_6, so edim_m(Q_6) <= 15: `paperSet_resolving`.
- L1: `trivial_stabilizer`.
- L2: `hist_antipode`.
- L3: `alt_sum`.
- L4: `not_resolving_of_length_le_six`, i.e. edim_m(Q_6) >= 7.
- Explicit sets: edim_m(Q_7) <= 19, Q_8 <= 26, Q_9 <= 38.
- Not formalized: the lower bound edim_m(Q_6) >= 15 (a 1.4e10-leaf search) and L5.

**Later work (2026-10-02; independent audit owed).**
- [`edge_multiset_dimension_q7_20261002.md`](../../05-knowledge/results/edge_multiset_dimension_q7_20261002.md):
  Q_7 has no resolving set of size <= 12 (exhaustive search, FINITE-EXACT), so 13 <= edim_m(Q_7) <= 19.
- [`edge_multiset_dimension_growth_20261002.md`](../../05-knowledge/results/edge_multiset_dimension_growth_20261002.md):
  ln edim_m(Q_d) = Theta(d^(1/3)) (HYP-9169; Open Problem 2), with L5's constant 0.6216 improved to 0.8146.
