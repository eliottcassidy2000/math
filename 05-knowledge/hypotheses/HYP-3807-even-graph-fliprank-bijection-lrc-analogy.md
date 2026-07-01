---
id: HYP-3807
title: THE TOURNAMENT<->EVEN-GRAPH BIJECTION IS A CUBE-ISO BUT NOT A FLIP-RANK-ISO (the extremal structure lives in the QUOTIENT, not the cube) + the LRC relation-lattice analogy. Two novel obligations constructed together via the cycle-space bijection (tile-vector in GF(2)^m, m=C(n-1,2), -> XOR of fundamental cycles of the base path -> even graph). The SAME cube Q_m carries TWO S_n-quotients: tournaments G_n (A000568=2,4,12,56) and even graphs E_n (A002854=2,3,7,16). OBLIGATION 1 (challenge the S72 assumption that flip-rank is cube-intrinsic): VERIFIED (exhaustive n=3..6) the EVEN-GRAPH flip-rank rho_E(n)=1,2,6,9 and rainbow R_E(n)=1,1,2,3, vs the TOURNAMENT rho_G=1,2,4,7 / R_G=1,2,3,5. They DIFFER sharply: even-graph covering-excess rho_E-ceil(log2|E_n|)=0,0,3,5 (vs tournament 0,0,0,1), and R_E < floor(log2|E_n|) at n=6 (3<4) -- so the S72 law "rainbow R(n)=floor" is NOT self-dual; it was SPECIFIC to the tournament quotient (the balanced-cut shape). CONCLUSION: the tournament<->even-graph bijection is a CUBE-isomorphism but NOT a flip-rank/packing-covering isomorphism -- these invariants depend on the QUOTIENT (the S_n-action), not the cube. (Diagnosis: even graphs have tiny orbits -- the empty even graph is a size-1 orbit -- and iso classes cut across the fund-cycle coordinates, so axis-aligned subcubes are inefficient.) OBLIGATION 2 (transport to LRC): the LRC analogue of "cycle space / even graph" is the RELATION LATTICE Lambda={t: sum t_i v_i=0} (THM-515); the # 3-term relations (a+b=c, the leading ADDITIVE ENERGY) is the analogue of the tournament 3-cycle count c3; the LONELY MEASURE is the even-graph theta; the PARITY LEMMA (S55, odd D => #lonely even) is the analogue of REDEI (# Ham paths odd). Verified: AP {1..n-1} has MORE 3-term relations than the construction (15 vs 10 at n=7, 78 vs 66 at n=14) = high-energy-tight vs lower-energy-loose (THM-515). LESSON for the LRC push: transporting results via a bijection requires preserving the ARITHMETIC quotient (Phi6/phase-residue), not just the abstract lattice -- mirroring "flip-rank depends on the quotient"
status: MIXED (exhaustive small-n + established bijection + analogy). VERIFIED (exhaustive n=3..6, fund-cycle cube): rho_E=1,2,6,9, R_E=1,1,2,3, |E_n|=2,3,7,16; even-graph orbit sizes computed (size-1 empty-graph orbit + spread classes). vs tournament rho_G=1,2,4,7/R_G=1,2,3,5 (HYP-3803/3804). R_E<floor at n=6 is exact. LRC: #3-term relations AP vs construction verified (15/10, 78/66). HONEST: the even-graph flip-rank is a clean NEW computation challenging the S72 self-duality assumption; the LRC analogy (relation lattice = even graph, additive energy = c3, parity lemma = Redei) is a structural DICTIONARY (THM-515 + S55 are established), a synthesis/transport frame, NOT a new proof. Exploratory + assumption-challenging.
source: klein-2026-07-01-S74
depends_on:
  - HYP-3803   # tournament flip-rank rho_G (S71); this is its even-graph dual
  - HYP-3804   # packing/covering asymmetry + R_G=floor (S72); this shows R=floor is NOT self-dual
related:
  - THM-515    # lonely measure = theta over the relation lattice Lambda (the LRC even graph); additive energy governs it
  - HYP-2873   # additive energy = spectral 4th moment (the LRC c3 analogue)
  - HYP-3806   # covering excess (this compares it across the two quotients)
external: A002854 (even/Eulerian graphs); the cycle-space (fundamental cycle) bijection; CLAUDE.md even-graph metagraph E_n; Redei's theorem
results:
  - 04-computation/even_graph_fliprank_lrc_analogy_klein.py
  - 05-knowledge/results/even_graph_fliprank_lrc_analogy_klein.out
---

# HYP-3807 — the bijection is a cube-iso, not a flip-rank-iso; and the LRC relation-lattice analogy

## The setup (the natural bijection)
The cycle-space bijection (CLAUDE.md): a tile-vector in `GF(2)^m` (`m = C(n-1,2)`) maps, via XOR of the
**fundamental cycles** of the base path `0-1-...-(n-1)`, to an **even graph** (all degrees even). The SAME
cube `Q_m` carries two `S_n`-quotients: **tournaments** `G_n` (A000568 `= 2,4,12,56`) and **even graphs**
`E_n` (A002854 `= 2,3,7,16`). Two novel obligations, constructed on this shared cube.

## Obligation 1 — the even-graph flip-rank (challenging the S72 self-duality assumption)
I had (implicitly, S72) assumed the packing/covering structure was a property of the cube-with-`S_n`-action.
Computing the EVEN-GRAPH flip-rank `rho_E` and rainbow `R_E` (exhaustive `n=3..6`) refutes that:

| `n` | `m` | `\|E_n\|` | `rho_E` | `ceil(log2 E)` | `R_E` | `floor(log2 E)` | | `\|G_n\|` | `rho_G` | `R_G` |
|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 1 | 2 | 1 | 1 | 1 | 1 | | 2 | 1 | 1 |
| 4 | 3 | 3 | 2 | 2 | 1 | 1 | | 4 | 2 | 2 |
| 5 | 6 | 7 | **6** | 3 | 2 | 2 | | 12 | 4 | 3 |
| 6 | 10 | 16 | **9** | 4 | **3** | 4 | | 56 | 7 | 5 |

- **Even-graph covering-excess** `rho_E - ceil(log2|E_n|) = 0,0,3,5` — vastly larger than the tournament's
  `0,0,0,1`. Even graphs are *much* harder to cover by subcubes.
- **`R_E < floor(log2|E_n|)` at `n=6`** (`3 < 4`) — the S72 law "**rainbow `R(n)` = floor**" is **NOT
  self-dual**; it was specific to the tournament quotient (the balanced-max-cut shape, which has no
  even-graph analogue in these coordinates).

> **The tournament<->even-graph bijection is a CUBE-isomorphism but NOT a flip-rank / packing-covering
> isomorphism.** These extremal invariants live in the **quotient** (the `S_n`-action), not the cube.
Diagnosis (orbit sizes): even graphs have **tiny orbits** — the empty even graph is a size-1 orbit
(`tile-vector 0`) — and their iso classes cut *across* the fundamental-cycle coordinates, so axis-aligned
subcubes are inefficient. The tournament's efficient "balanced cut" is a tournament-specific shape.

## Obligation 2 — transport to LRC: the relation lattice is the LRC's "even graph"
The LRC analogue of the tournament cycle space (even graphs) is the **relation lattice**
`Lambda = { t : sum_i t_i v_i = 0 }` (THM-515). The dictionary:

| tournament (cycle space / even graph) | LRC (relation lattice `Lambda`) |
|---|---|
| # 3-cycles `c3` (leading OCF term) | # 3-term relations `a+b=c` (leading **additive energy**) |
| the even graph itself | the lattice `Lambda` |
| even-graph generating function | the **lonely measure** `L = theta over Lambda` (THM-515) |
| **Redei**: # Hamiltonian paths is odd | **Parity lemma** (S55): odd `D` => # lonely times even |
| even-graph metagraph `E_n` | the moment/singular-series structure of `L` |

VERIFIED: the tight AP `{1..n-1}` has MORE 3-term relations than the construction (`15` vs `10` at `n=7`;
`78` vs `66` at `n=14`) — high additive energy (AP, tight, `L` small) vs lower (construction, loose), exactly
THM-515's "high additive energy <=> low lonely measure."

## The lesson for pushing LRC (challenge assumptions)
Obligation 1's moral transports: **a bijection of the underlying set need not preserve the extremal
structure — the quotient carries the content.** So when pushing LRC via any "dual" or bijection (to even
graphs, to a lattice, to a tournament), one must preserve the **arithmetic quotient** (the modulus `Phi6`,
the phase-residue `p(w)=nw mod Phi6`, S68), not merely the abstract lattice. The covering-min value depends
on the arithmetic of the speeds, just as the flip-rank depends on the `S_n`-action — this is why the naive
lattice/geometry-of-numbers predictor (`lambda_1`) failed (THM-515) and the additive-energy/phase-residue
predictors succeeded: they see the quotient. The relation-lattice/even-graph frame is the right *object*;
the arithmetic quotient is the required *decoration*.

## Net
Two obligations built on the tournament<->even-graph bijection: (1) the even-graph flip-rank `rho_E =
1,2,6,9` and rainbow `R_E = 1,1,2,3` differ sharply from the tournament's, so the bijection is a cube-iso
but not a flip-rank-iso — the S72 "rainbow = floor" law is quotient-specific (assumption challenged); (2)
the LRC relation lattice is the analogue "even graph," with additive energy = `c3`, lonely measure = the
even-graph theta, and the parity lemma = Redei. The transferable lesson: the quotient/arithmetic, not the
underlying set, carries the extremal content — a caution and a compass for LRC transport. Exploratory +
assumption-challenging; not a new proof.
