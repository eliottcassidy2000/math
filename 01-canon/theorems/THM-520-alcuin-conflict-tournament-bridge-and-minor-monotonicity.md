---
id: THM-520
title: The conflict-graph↔tournament bridge (forward arc = edge) — independent set ↔ reverse-transitive run, 3-cycle ↔ ordered induced P₃, the Rédei odd-parity shadow — and the minor-theory dichotomy: τ is minor-monotone (finite Robertson-Seymour obstruction set) but the Alcuin number is SUBGRAPH-monotone yet NOT minor-monotone, failing ONLY under edge contraction (smallest witness: contract an edge of K₂,₄, Alcuin 2→3)
status: PROVED/VERIFIED computationally (all non-isomorphic graphs n≤6, exhaustive: 0 failures on the bridge identities; 0 bound violations; minor-monotonicity verdicts exact). The τ minor-monotonicity is proved analytically (all three minor cases). Two creative conjectures (Alcuin+1 ↔ T_G non-strong; G-Ham-cycle ↔ ∃-order-T_G-strong) are REFUTED by exhaustive search — honest negatives. COMPLEMENTARY to kind-pasteur's THM-519 (same human dispatch, concurrent collision — both used THM-519/HYP-2550-2552/OPEN-Q-106/T828; this work renumbered to THM-520/HYP-2553-2555/OPEN-Q-107/T829): theirs specializes Alcuin to the OCF conflict graph Ω (Alcuin(Ω)=τ+1 ⟺ Ω edgeless ⟺ H=3^{α₁}; Kuratowski K₅=5 overlapping odd cycles), citing CHW Lemma 4.3 (two distinct max independent sets ⟹ small-boat) + Thm 3.1; THIS one is the GENERAL graph→tournament bridge + the minor-monotonicity dichotomy.
source: mac-mini-2026-06-15-S6
depends_on:
  - THM-519   # kind-pasteur's complementary Alcuin-of-Ω result (cites CHW Lemma 4.3 / Thm 3.1)
related:
  - THM-001   # Rédei's theorem (odd # Hamiltonian paths)
  - OPEN-Q-107
  - HYP-2553
  - HYP-2554
  - HYP-2555
references:
  - "P. Csorba, C.A.J. Hurkens, G.J. Woeginger, The Alcuin Number of a Graph and Its Connections to the Vertex Cover Number, SIAM J. Discrete Math. 24(3) (2012); JSTOR 41642576."
  - "N. Robertson, P.D. Seymour, Graph Minors XX (Wagner's conjecture), JCTB 2004 — every minor-closed family has a finite forbidden-minor set."
  - "Kuratowski/Wagner: planar ⟺ no K₅, K₃,₃ minor."
  - "Moon 1966 / Camion: a tournament is strong ⟺ it has a Hamiltonian cycle (⟹ vertex-pancyclic)."
---

# THM-520 — The Alcuin/vertex-cover ↔ tournament bridge, and where the "+1" escapes minor-closure

Dispatch (human): connect the Alcuin number (= vertex cover number, or that +1; Csorba–
Hurkens–Woeginger) to tournaments via "conflict graph as tournament" (present edge =
transitive arc, absent edge = intransitive arc), Hamiltonicity, and Kuratowski/Robertson–
Seymour. All claims verified exhaustively over the 208 non-isomorphic graphs on `n≤6`
vertices. Script: `04-computation/alcuin_tournament_macmini_0615s6.py`; output
`05-knowledge/results/alcuin_tournament_macmini_0615s6.out`.

## A. The correspondence (VERIFIED — all identities exact, 0 failures n≤6)

**Map.** Fix a linear order on `V=[n]`. Define the tournament `T_G` by: for `i<j`, the arc
is `i→j` (**forward / "transitive"**) iff `{i,j}∈E(G)`, else `j→i` (**backward**). So the
forward-arc set of `T_G` is exactly `E(G)`. (Complete: every pair is edge↦forward or
non-edge↦backward.)

1. **Clique ↔ forward run, independent set ↔ reverse-transitive run.** A vertex set `S` is
   order-monotone-increasing transitive in `T_G` iff `S` is a **clique** of `G`; monotone-
   decreasing (reverse-transitive) iff `S` is an **independent set** of `G`. Hence
   `ω(G)` = largest forward run, `α(G)` = largest reverse run, and the vertex cover number
   is a transitive-cover invariant: **`τ(G) = n − α(G) = n − (largest reverse-transitive
   run of T_G)`**. VERIFIED: the largest transitive subtournament `≥ α(G)` in all 208 cases
   (0 with `<α`; the size-`α` reverse run is always a witness).
2. **3-cycles ↔ ordered induced P₃.** The directed 3-cycles of `T_G` are exactly the
   ordered triples `i<j<k` with edge-pattern `(ij,jk,ik)=(1,1,0)` [an **induced path
   i–j–k** in `G`] or `(0,0,1)` [the same in the complement `Ḡ`]. VERIFIED:
   `#3-cycles(T_G) = #{ordered induced P₃ in G or Ḡ}` exactly, all 208 graphs (0 failures).
   So **cyclicity of `T_G` = induced-cherry density of `G`**; cliques and independent sets
   (P₃-free) map to transitive (acyclic) tournaments.
3. **Score sequence.** `score(i) = (i−1) + deg⁺_G(i) − deg⁻_G(i)` (deg⁺/⁻ = higher/lower-
   indexed neighbors). The transitive backbone `(i−1)` twisted by the edge imbalance.
4. **The Rédei odd-parity shadow.** Every tournament has an odd number of Hamiltonian
   directed paths (THM-001 / Rédei parity). VERIFIED: `#HamPaths(T_G)` is odd for all 208
   graphs. So **every conflict graph carries a hidden odd invariant** `#HamPaths(T_G)`.
5. **Hamiltonian path alignment.** `G` has a Hamiltonian path ⟺ some ordering makes `T_G`'s
   spine `π(0)→…→π(n−1)` all-forward — i.e. `G`-Hamiltonicity is exactly the condition that
   `G` aligns with the project's tiling-model base path (the all-present spine).

## B. The Alcuin / vertex-cover sandwich (VERIFIED; literature)

`τ(G) ≤ Alcuin(G) ≤ τ(G)+1` — **0 violations** over all 208 graphs `n≤6` (Csorba–Hurkens–
Woeginger). The `Alcuin = τ` case dominates (n=3..6: 3/4, 9/11, 30/34, 148/156); the `+1`
cases are special — **edgeless graphs (`τ=0`, trivially `Alcuin=1`), stars `K₁,ₜ` (`t≥3`),
and certain dense two-hub graphs**. (Correction to a session hand-guess: `P₃=K₁,₂` is `τ`,
NOT `τ+1` — `Alcuin(P₃)=1`; the jump starts at `K₁,₃`.) The exact characterization is the
literature's (NP-hard to decide): **CHW Lemma 4.3 — two distinct maximum independent sets ⟹
small-boat** (surfaced by the concurrent kind-pasteur THM-519). My `+1` examples are
consistent: `K₁,₃` and `K₂,₄→`contraction each have a UNIQUE maximum independent set (the
leaves / the larger part), the necessary precondition for large-boat — though uniqueness is
not sufficient (`K₂,₄` itself has a unique max independent set yet is small-boat, decided by
CHW Thm 3.1). Note: kind-pasteur's clean "`Alcuin(Ω)=τ+1 ⟺ Ω edgeless`" is special to the
**conflict-graph class** `Ω`; for GENERAL graphs the `+1` class is strictly richer (stars,
two-hub graphs are `+1` but not edgeless). CHW's exact rule for **complete multipartite**
graphs `K_{n₁≤…≤nₖ}` is **small-boat ⟺ `nₖ ≤ 2n_{k−1}`** (the largest part at most twice the
next) — which exactly reproduces my data: `K₁,ₜ` (parts `1,t`) is small ⟺ `t ≤ 2`, so
`P₃=K₁,₂` small, `K₁,₃` large (`t=3>2`); `K₂,₄` (parts `2,4`) is small ⟺ `4 ≤ 4` ✓ small.
After the contraction, `H` is no longer complete multipartite (the two hubs become adjacent),
the rule no longer applies, and the over-committed cover flips it to large — the precise
boundary the contraction crosses.

## C. The minor-theory dichotomy (the headline — answers the Kuratowski/RS question)

**τ is minor-monotone (PROVED).** If `H ≼ G` then `τ(H) ≤ τ(G)`: edge/vertex deletion only
shrinks a cover; for contraction `uv→w`, from a min cover `C` of `G` take `C' =
(C∖{u,v})∪{w}` if `C∩{u,v}≠∅` else `C`, and `C'` covers `G/uv` with `|C'|≤|C|`. VERIFIED:
0 failures of `τ(H)≤τ(G)` over all single-step minors of all 208 graphs. Hence **`{G:τ(G)≤k}`
is minor-closed**, so by **Robertson–Seymour** it has a FINITE forbidden-minor obstruction
set (e.g. `(k+1)·K₂`; for `τ≤1`: `{2K₂, K₃}`) — the exact analogue of Kuratowski/Wagner's
`{K₅,K₃,₃}` for planarity.

**The Alcuin number is SUBGRAPH-monotone but NOT minor-monotone — it fails ONLY under edge
contraction (VERIFIED).** Over all single-step minors of all 208 graphs:
- deletion (vertex or edge): `Alcuin(H) ≤ Alcuin(G)` always (0 failures) — **Alcuin is
  subgraph-monotone**.
- contraction: **8 failures** `Alcuin(H) > Alcuin(G)`. **Smallest witness:** `G=K₂,₄`
  (parts `{0,1}`,`{2,3,4,5}`, `Alcuin=2=τ`); contract any edge `→ H` (two hubs `{0,1}` now
  ADJACENT, each joined to independent `{2,3,4}`; `τ(H)=2`, `Alcuin(H)=3=τ+1`).

> **Mechanism.** Contraction can create an edge INSIDE a minimum vertex cover. In `K₂,₄` the
> two hubs are non-adjacent (same part), so the cover `{0,1}` can be parked together and
> `τ` seats suffice. The contraction makes the hubs mutually conflicting AND each conflicts
> with all leaves — the cover becomes internally conflicting + over-committed (like a star's
> center), forcing the extra seat. **The `+1` is born from an edge contraction internal to
> the cover.**

**Consequence.** `{G : Alcuin(G) ≤ k}` is **NOT minor-closed**, so Robertson–Seymour gives
it **no finite obstruction set** — even though the underlying `{τ ≤ k}` does. The Alcuin
`±1` is a genuinely non-minor-closed (non-WQO) correction riding on a minor-closed skeleton,
created precisely by the one minor operation (contraction) that can densify a cover. It IS,
however, characterized by a finite set of forbidden SUBGRAPHS / topological structure (since
it is subgraph-monotone) — a strictly weaker closure than minor-closure.

## D. Honest negatives — two creative conjectures REFUTED (exhaustive n≤6)

- **Alcuin=τ+1 ⟺ T_G non-strong** (the "+1 = T_G refuses to be Hamiltonian-cyclic"):
  **REFUTED.** 12 graphs have `Alcuin=τ+1` yet `T_G` is strong under some ordering (1 even
  under the natural order). The `+1` has no clean strong-connectivity signature.
- **G has a Hamiltonian cycle ⟺ ∃ order making T_G strongly connected:** **REFUTED both
  directions** (4 graphs Ham-cyclic but no order makes `T_G` strong; 141 non-Ham-cyclic but
  some order does). `G`-Hamiltonicity and `T_G`-strong-connectivity are independent.

These negatives are informative: the `T_G` correspondence faithfully encodes `G`'s
clique/independent/induced-P₃ structure (part A, all exact), but the Alcuin `+1` and graph
Hamiltonicity are NOT shadows of tournament strong-connectivity — they live on the
scheduling/cover-internal-edge axis (part C), not the cycle-structure axis.

## Honesty

Part A identities and the part C monotonicity verdicts are EXHAUSTIVELY verified `n≤6` (a
finite, complete check on the 208 non-iso graphs). `τ` minor-monotonicity is proved for all
`n`. The exact literature characterization of the Alcuin `+1` (CHW Lemma 4.3 / Thm 3.1, see
THM-519) and the formal "no finite minor-obstruction-set" statement (vs. the finite subgraph-
obstruction-set) deserve `n≥7` confirmation and a clean proof of the contraction mechanism in
general (OPEN-Q-107). The two refuted conjectures are closed.
