# Tournaments and graphs with one fewer edge than vertices: the Euler-characteristic trijection, its tiers, and procedural counterexample families

kind-pasteur-2026-10-27-S8. This builds on the owner's trijection (kind-pasteur S3: a tournament with a marked line,
cycle+solo or star is a tournament on one vertex fewer, modulo 1, `Z_{n-1}` or `S_{n-1}`) and on the symmetry calculus
of session80 (composition, cycle indices, explicit collisions). It was written concurrently with session81, which read
the same prompt as an encoding-capacity question and already posed the pair-coordinate question answered here (§6).

- **Code and data:** `04-computation/tournament_carriers/`.
- **Pre-registration:** `PREREGISTRATION.md` (P1–P3), pushed in efdf0b0a5f before the confirmatory runs; the referee
  verified the chronology.
  - P1 PASSED.
  - P2 FAILED and is recorded.
  - P3 PASSED.
- **Audit:** independent referee, `audit_20261027/REPORT.md` (MAJOR 1, MINOR 5, NOTE 8). All repairs are applied
  (MISTAKE-689). The referee also proved full necessity (now THM-4835 part 4) and refuted HYP-9435(3).
- **Canon:** THM-4835, THM-4836, HYP-9435.

## 1. Inheritance, board, niche, wildcard

- **Closest proved mechanism.** S3 Theorem 1 (three marked bijections); LEM-003.
- **Hostile.** Session80's collision: equal counts, different groups. Here the bigon `C_4 + K_1` has the same marked
  count as `Tour(4)/Z_2`, but admits no pair-level trijection.
- **Corrected near miss.** S3's `8 = 4 + 4` (n = 4 only).
- **Least-used sidecar.** The Burnside ring / table of marks.
- **Niche.** Non-functional orientations of disconnected excess -1 graphs.
- **Wildcard.** Maxwell rigidity counts.
- **META-PATTERNS cards used.**
  - "The same representation is not the same carrier": permutation characters, and even F2-affine structure, do not
    decide the signed A-set (audit m1).
  - "Controlled forgetting and unlabeled quotients require a sidecar".

## 2. The answer in one paragraph

The owner's "graphs with one fewer edge than vertices" are the Euler-characteristic-1 graphs: trees, and disconnected
graphs whose extra cycles are paid for by extra components. For each such graph the complement has `C(n-1,2)` pairs,
so a tournament with a marked oriented copy D is, as a set, a tournament on `n-1` vertices.

The identification is equivariant and pair-level (linear), with the reduced vertex set the edge set of D, **exactly
when** every symmetry of D fixes a sub-carrier of Euler characteristic 1. Equivalently, `V ≅ E ⊔ {pt}` as
Aut(D)-sets (THM-4835).
- Sufficiency is the session's proof.
- Full necessity is the referee's (lemmas L1–L4).
- Trees in any orientation, and functional orientations, always qualify, with explicit maps. The S3 trijection and
  the vertex-deleted deck are instances.
- Symmetric bigons, alternating even cycles, dense cores with twin isolates and odd thetas do not.

The tiers behave as follows.
- Involution defects destroy only the pair-level structure; counts and A-sets never see them.
- When A fixes a vertex, destroying the A-set structure requires an odd-order defect.
- An odd-order defect does not always destroy it: `TT_4 + 3K_1` survives (audit M1).

As "each tournament maps to a pair or a triple" of carrier kinds, the statement is exact precisely for `n <= 5`
(THM-4836); the first singleton is `C_3 ⇒ C_3`.

## 3. Results

### 3.1 The criterion (THM-4835 parts 2, 4)

- **Reduction to marks.** View the ordered non-edges as an `(A x Z2)`-set. Burnside's marks reduce linear
  trijectivity to two families of equations over subgroups K:
  - `C(v_K,2) - e_K = C(x_K,2)`;
  - `v_{K0} - v_K = x_{K0} - x_K` for `K0` of index 2 in K.
- **Sufficiency.** `X = V \ {r}`.
- **Necessity (referee).**
  - Twinned components have `χ = 0` (a transposition argument).
  - `δ_K >= 0` (relative homology).
  - `φ_K = x_K - v_K + 1 ∈ {-1, 0, 1}` is constant along prime-index steps, so solvable K have no defect.
  - Fixed sets of elements contain the 2-core, so cyclic suffices.

### 3.2 Explicit trijections (THM-4835 part 3)

- **Trees.** `Tour(V) ≅ Or(G) x Tour(E(G))` via the end-edge map (cf. the line-graph Wiener identity of Buckley).
  For the path this is S3's staircase; for the star it is `T -> T - v`.
- **Functional orientations.** The root-swap map `{v, f(v)} -> {r, v}`. For the directed cycle plus solo it is S3's
  map.
- **Verification.**
  - Author: 600/600 random cases at n = 8, 9.
  - Referee: all trees with `n <= 9`, all functional digraphs with `n <= 6`, and random samples up to n = 10.

### 3.3 The tiers (THM-4835 parts 5–8)

| tier | sees | what destroys it | witnesses |
|---|---|---|---|
| counts `\|{T ⊇ D}/A\|` | odd-order subgroups only | some odd-order defect (if A fixes a vertex; necessary, not sufficient) | `θ_3^k + 2K_1` (proved for all k); `(p, k) = (5, 2)` checked; `K_4(C_3+apex) + 3K_1`; `K_{3,2} + 2K_1` |
| A-set `{T ⊇ D}` (weak) | odd-order subgroups only | same condition | `θ_p^k + (p-1)K_1` (all odd primes p, `k >= 2`, proved); the n = 7 sporadics. Survives an odd defect in `TT_4 + 3K_1` |
| F2-affine | permutation module plus cocycle | strictly weaker than linear | rigid diamond + 2K_1 (n = 6): affine, not linear (audit m1) |
| signed A-set (linear) | all subgroups | any defect | bigons, alternating cycles, dense cores + twin isolates |

**Census** (linear / weak only / neither):

| n | carriers | linear | weak only | neither | source |
|---|---|---|---|---|---|
| 3 | 3 | 3 | 0 | 0 | author and referee |
| 4 | 10 | 10 | 0 | 0 | author and referee |
| 5 | 41 | 39 | 2 | 0 | author and referee |
| 6 | 187 | 167 | 20 | 0 | author and referee |
| 7 | 1009 | 816 | 188 | 5 | author and referee |
| 8 | 5781 | 4155 | 1469 | 157 | referee, out of sample |

- Linear iff criterion iff cyclic criterion in all 7031 carriers.
- Group dependence: `TT_4 + 3K_1` is linear for `Z_3` but not for the full `S_3`. That is exactly why lemma L1
  needs the full group.

### 3.4 Pairs and triples (THM-4836)

- **The three classes.** `Σ(T) ⊆ {L, C, S}`: L by Rédei, C for strong T or strong `(n-1)` plus a source/sink
  (Camion, Moon), S for a source or sink.
- **Class sizes.**

  | `Σ(T)` | types |
  |---|---|
  | `{L,C}` | `s_n` |
  | `{L,C,S}` | `2s_{n-1}` |
  | `{L,S}` | `2t_{n-1} - t_{n-2} - 2s_{n-1}` |
  | `{L}` | `[x^n] S_{>=3}^2 T` |

- **Singleton counts** by n:

  | n | 6 | 7 | 8 | 9 | 10 |
  |---|---|---|---|---|---|
  | `{L}` types | 1 | 3 | 16 | 99 | 915 |

  The values for n = 9, 10 are the referee's.
- **The pair-or-triple law** holds iff `n <= 5`.

### 3.5 The Maxwell ladder (HYP-9435(2); pre-registered P2 FAILED)

- **Setup.** Marking `jn - C(j+1,2)` edges reduces by j vertices ("pairs" at j = 2, "triples" at j = 3).
- **The naive count fails** (j = 2, n = 6: 6/3817; j = 3: 14/2790) at subgroups of order divisible by 3, with an
  alternative fixed count x.
- **j = 1, the trijection, is the rigid rung** (THM-4835 part 4).

### 3.6 A free identity (THM-4835 corollary)

- `Σ_T Σ_r ∏_{v≠r} s_v(T) = F_n 2^{C(n-1,2)}`.
- `F_n = 9, 72, 815, 11700, 201747` for n = 3..7; the referee checked n = 7.

## 4. Procedural generators of counterexamples

Each comes with a proof.

1. **Involution gadgets** (destroy linearity only, when A fixes a vertex and the odd part is clean): symmetric
   bigons, alternating even cycles, twin isolates next to a dense core.
2. **Odd gadgets**: p identical directed u-w paths plus `p-1` isolates. The `Z_p` mark forces `x = p+1 ≢ 0 (mod p)`,
   so there is no A-set model. For p = 3 there is no count model either.
3. **Chirality.** Remove all involutions so that only odd defects remain. This gives non-linear carriers with no
   involution witness: the referee's n = 10 carrier `K_4(C_3 + source) + K_1 + K_2 + vec P_3`, and the session's
   n = 19 chiral theta.

## 5. Unifying proposal (typed)

**PROVED.**
- Marking an oriented excess -1 carrier is the same as passing to its edge set. The price is one vertex, paid
  equivariantly iff the carrier is a point in the Burnside ring (`[V] - [E] = 1`).
- One statement covers:
  - S3's trijection;
  - the deck;
  - Hamiltonian-path staircases;
  - functional digraphs (`∏ s_v`);
  - arborescences.

**PROVED (corrected tier law).**
- Counts and A-sets see only odd-order subgroups (odd-order rule).
- Involution defects are visible only to the signed structure.
- When A fixes a vertex, failure of the A-set tier requires an odd-order defect.

**ANALOGY (no implication).**
- The Maxwell ladder.
- Session80's information hierarchy.
- Session81's "judge a representation against its future operations".

## 6. Concurrent and prior work

**Session81** (`session81_counterexample_families_20261028`, pushed before this session's first commit).
- **Mark law and the open question.** `GRAPH_CARRIERS.md` §5 "General mark law" states two things. Marked classes are
  the A-orbits of the affine completion cube, and each element fixes 0 or `2^{o}` completions. And a correspondence
  with tournaments on `n-1` vertices "requires a compatible pair-coordinate action". THM-4835 answers when that
  action exists. The counting half of the tier law (5(a)) is session81's mark law.
- **Capacity.** It also proves the capacity obstruction for independent pairs/triples of `(n-1)`-edge graph types.
- **`K_4 + 3K_1`.** It finds this graph as the first non-forest odd-cycle conflict graph with `E = V - 1`. That is a
  different role from our non-weak carrier `K_4(C_3+apex) + 3K_1`; the shared graph is a coincidence, noted and not
  identified.

**Royle–Praeger–Glasby–Freedman–Devillers** (arXiv:2204.01947).
- Tournaments and *even graphs* are equinumerous. The orientation independence of the parity is due to von Brömssen,
  and "sign" is session81's term for that parity.
- No implication is claimed.

**Classical inputs.**
- Rédei; Camion; Moon and Harary–Moser (pancyclicity).
- Burnside's marks; tom Dieck's equivariant Euler characteristic; Serre (fixed subtrees).
- Buckley 1981 (the line-graph Wiener identity, cf. the end-edge map).

## 7. Reproduction

| directory | command | covers |
|---|---|---|
| `04-computation/tournament_carriers/` | `python carrier_census.py 7` | census |
| | `python families_and_controls.py` | P1 |
| | `python maxwell_ladder.py` | P2 |
| | `python coherent_carrier_sets.py` | P3 |
| | `python tiers_and_sporadics.py` | tiers and sporadic carriers |
| `audit_20261027/` | see `REPORT.md` §Reproduction | independent referee |
