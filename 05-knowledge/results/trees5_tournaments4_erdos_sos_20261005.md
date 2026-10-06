# The three trees on five vertices and the three merged classes of four-vertex tournaments: one Klein group, two notions of reversal; the Erdős–Sós half/whole alternation at `t = 5` as a mode pair; what reaches Collatz

**Session:** mac-mini (claude) 2026-10-05, worktree `codex/session-trees-tournaments-20261005`.
**Owner's seed (verbatim gist):** the 3 unlabeled trees on 5 vertices correspond to the isomorphism
classes of 4-vertex tournaments once the two mutually converse classes are merged (3 nodes); trees
are "the structures that create a cycle when an edge is added"; this "deeply aligns" with the
fractal ternary trees studied here; the Erdős–Sós bound `(t-2)n/2` alternates whole/half and this
relates to the repo's modes of natural recursion; weave into the Collatz work.

**Status.** PROVED (elementary, author-audited; independent audit OWED): Propositions 1–5.
FINITE-EXACT: the censuses in §1, §3, §6 (script below, 1,054,454 checks). CITED: Erdős–Gallai 1959,
Faudree–Schelp 1975, McLennan 2005, Otter 1948, OEIS A000055 / A000568 / A002785 / A059735.
Typing of the owner's four claims is in §7; nothing here changes the status of any Collatz
statement (Collatz OPEN). Script: `04-computation/experiments/trees5_tournaments4_erdos_sos_20261005.py`;
output: `trees5_tournaments4_erdos_sos_20261005.out` (this directory).

## 0. Inheritance and concept board

- Closest proved mechanism: THM-584 (`THM-584-complement-is-antipodal-map-level-parity-spectrum`),
  Klein-`V_4` base case: the four 4-classes are `(Z_2)^2` under "source destroyed" / "sink destroyed";
  the merged metagraph `G_4/Z_2 = {T, [±], S}` is labelled by `u = #boundary defects ∈ {0,1,2}`.
- The owner's earlier four-vertex claim, already adjudicated: THM-4472
  (`four-vertex-reading-of-3n-plus-minus-1`): the converse pair `(0,2,2,2)/(1,1,1,3)` is forward map
  vs inverse tree (time reversal) on the AM-fair quadruple `{i, 2i-1, 2i, 3i-1}` of THM-4470; converse
  = negation ∘ time reversal; REFUTED as the sheet swap.
- Canonical hostile example: MISTAKE-049/081 (`SC(n) = A000568(n-1)`, a tournament-count identity
  built on two small matches) and MISTAKE-554/559 (small-number coincidences typed STRUCTURAL /
  DICTIONARY). Rule applied: NUMEROLOGY unless a map carries the mechanism, and a structural claim
  must extend beyond the instance that prompted it.
- Least-used sidecar: the repo's "recursion modes" are cyclotomic factors `(x-1)^depth · Φ_d`
  (`07-reflections/the-recursion-modes-are-cyclotomic-factors-depth-times-phi-d.md`; THM-1955/1960):
  EVEN mode `(x-1)^2`, FULL `(x-1)^3`, ODD `(x-1)^3 Φ_2`; and Mode A (`n -> n-1`) vs Mode B
  (`n -> n-2`, "the natural recursion", `the-two-staircases.md:72`).
- Trees inside tournaments already in canon: THM-4533 / HYP-9173 (Rédei graphs; a tree is good iff
  hereditarily symmetric). Erdős–Sós and Sumner are greenfield (`PROBLEM-LEDGER.md:379`; the S10 note
  of 2026-10-05 names Erdős–Sós once, as a heuristic rhyme).
- Board: {end-defect Klein group; converse vs reflection; `K_4` as the hinge clique; handshake parity;
  `Φ_2` and `Φ_4`; the AM-fair quadruple}.

## 1. The census, and where the count coincidence stops (FINITE-EXACT)

| `m` (tournament order) | iso classes A000568 | self-converse A002785 | converse-merged A059735 | trees on `m+1` vertices A000055 |
|---|---|---|---|---|
| 1 | 1 | 1 | 1 | 1 |
| 2 | 1 | 1 | 1 | 1 |
| 3 | 2 | 2 | 2 | 2 |
| 4 | 4 | 2 | **3** | **3** |
| 5 | 12 | 8 | 10 | 6 |
| 6 | 56 | 12 | 34 | 11 |
| 7 | 456 | 88 | 272 | 23 |

Recomputed here for `m ≤ 6` by brute force over all labelled tournaments (canonical forms), `m = 7`
from OEIS. The owner's `3 = 3` holds, and the equality of the two sequences holds exactly for
`m ≤ 4`. The repo already knows the merged count as `V_merged = (A000568 + SC)/2 = A059735`
(THM-283; `sc_vmerged_exact_burnside_s560.out`). As a *count* identity the correspondence is a
small-number coincidence of the MISTAKE-049 kind. §2 says what is structural.

The three trees on 5 vertices: path `P_5` (degrees 2,2,2,1,1; diameter 4; |Aut| = 2), fork
`S(2,1,1)` (3,2,1,1,1; diameter 3; |Aut| = 2), star `K_{1,4}` (4,1,1,1,1; diameter 2; |Aut| = 24).
The four 4-tournaments: `TT_4 (0,1,2,3)`, `H = 1`; `(0,2,2,2)` 3-cycle over a sink, `H = 3`;
`(1,1,1,3)` source over a 3-cycle, `H = 3`; `S_4 (1,1,2,2)` strong, `H = 5` (Rédei: all odd).

## 2. Proposition 1 (PROVED): one Klein group, two notions of reversal

Let `TT_m` be the transitive tournament on `1 < 2 < … < m` (arc `i -> j` for `i < j`) and `P_{m+1}`
the path `1 – 2 – … – m+1`. Define

- on `TT_m`: `g_L` = reverse the arc `(1,3)`; `g_R` = reverse the arc `(m-2, m)`;
- on `P_{m+1}`: `f_L` = reattach the leaf `1` from `2` to `3`; `f_R` = reattach the leaf `m+1` from `m` to `m-1`.

**(a)** Each pair consists of two commuting involutions, so each generates an action of
`V_4 = (Z_2)^2` ("end defects": left end, right end).

**(b)** The four tournaments `TT_m, g_L TT_m, g_R TT_m, g_L g_R TT_m` lie in four distinct isomorphism
classes; `g_L TT_m` and `g_R TT_m` are converse to each other; `TT_m` and `g_L g_R TT_m` are
self-converse. At `m = 4` they are exactly `TT_4, (0,2,2,2), (1,1,1,3), S_4`, i.e. THM-584's Klein
base case with `x = g_L` (source destroyed) and `y = g_R` (sink destroyed).

**(c)** The four trees `P_{m+1}, f_L P_{m+1}, f_R P_{m+1}, f_L f_R P_{m+1}` lie in three isomorphism
classes: `f_L P_{m+1} ≅ f_R P_{m+1}`. At `m = 4` they are exactly `P_5`, fork, fork, `K_{1,4}`.

**(d)** The order reversal `ρ(i) = m+1-i` (resp. `m+2-i`) conjugates `g_L <-> g_R` and `f_L <-> f_R`.
On the path, `ρ` is an **automorphism**; on `TT_m` it is an **anti-automorphism** (an isomorphism onto
the converse). Hence the identification of the two one-defect objects is automatic for trees and
requires the converse quotient for tournaments. The merged label is the same on both sides:

    u = #end defects = (THM-584's boundary-defect label) = m - diameter(tree) ∈ {0, 1, 2}.

For the 5-vertex trees `u = 0, 1, 2` is path, fork, star; `H = 2u + 1 = 1, 3, 5` and
`#leaves = u + 2 = 2, 3, 4`.

*Proof.* (a) The two arcs `(1,3)` and `(m-2,m)` are distinct for `m ≥ 4`, so the reversals commute;
likewise the two reattachments touch disjoint edge sets for `m ≥ 4` (`{1,2},{1,3}` vs
`{m,m+1},{m-1,m+1}`). (b) Scores: `g_L TT_m` is the ordinal sum `C_3 ⊕ TT_{m-3}` (the 3-cycle
`1 -> 2 -> 3 -> 1` dominating `4..m`) and `g_R TT_m = TT_{m-3} ⊕ C_3`; their score sequences are
`(0,…,m-4, m-2,m-2,m-2)` and `(1,1,1, 4,…,m-1)`, distinct from each other and from the transitive
`(0,…,m-1)` and from the two-defect sequence; converse maps score `s` to `m-1-s`, which exchanges the two
one-defect sequences and fixes the other two; the tournaments `TT_m` and `g_L g_R TT_m` are carried onto
their converses by `ρ` (so they are self-converse), and `ρ(g_L TT_m) = converse(g_R TT_m)`. (c) `f_L P`
is the spider with legs `(m-2, 1, 1)` at vertex `3`; `f_R P` is the same spider at vertex `m-1`; `ρ`
maps one to the other. `f_L f_R P` has two vertices of degree 3 (`3` and `m-1`) for `m ≥ 5`,
adjacent at `m = 5` (the H-shaped tree) and at distance `m-4 ≥ 2` for `m ≥ 6`; at `m = 4` both folds
land on the same vertex `3`, giving `K_{1,4}`. Degree sequences separate the path, the fork and the
double fold. Diameters `m, m-1, m-2`. (d) `ρ` reverses the linear order: it maps
edges of `P_{m+1}` to edges (automorphism) and arcs `i -> j` of `TT_m` to arcs `ρ(j) -> ρ(i)`, i.e.
onto the converse. ∎ (FINITE-EXACT confirmation for `m = 4..7`, including Rédei counts
`H = 1, 3, 3, 5` at `m = 4` and `1, 3, 3, 9` at `m = 5, 6, 7`.)

**Scope.** At `m = 4` the Klein orbit exhausts all four classes and all three trees. For `m ≥ 5` it is
the boundary layer only: 10 merged classes vs 6 trees at `m = 5`. So the owner's correspondence is
STRUCTURAL as the statement "the end-defect `V_4` of the linear order has four tournament orbits (one
converse pair) and three tree orbits, matched by `u`", and the equality of the *full* censuses at
`m = 4` is an EXPLAINED COINCIDENCE (everything on four points is boundary).

**Weighted metagraph check.** The arc-reversal metagraph on 4 vertices (rows = class of a
representative, entries = number of its 6 arcs whose reversal lands in the column class) is

```
            TT4   A=(1,1,1,3)   B=(0,2,2,2)   S4
TT4          3        1             1          1
A            3        0             0          3
B            3        0             0          3
S4           1        1             1          3
```

(labelled flow symmetric: `24·1 = 8·3`). Merging `A, B` into `M` gives rows `TT4: (3,2,1)`,
`M: (3,0,3)`, `S4: (1,2,3)`, a triangle with loops at `TT4` and `S4`, none at `M`, eigenvalues
`{-2, 2, 6}` = THM-584's observed `V_+` spectrum at `n = 4` (reproduced). The edge-rotation graph of
the three trees is a *path* `P_5 – fork – star` with loops at all three (weights 8/8, 8/5/1, 12). So
the two "metagraphs" agree on nodes and on the label `u`, not on adjacency: the tournament side has
the `TT_4 – S_4` edge (reverse the source–sink arc; THM-220's near-transitive construction), the tree
side does not. Destroyed by the dictionary: adjacency and loop structure. Sidecar restoring the
tournament side: the second coordinate `w = x - y ∈ {-1, 0, 1}` of THM-584.

## 3. Proposition 2 (PROVED): "add an edge, get a cycle" — the linear-order dictionary

Trees are the maximal acyclic graphs; transitive tournaments are the acyclic tournaments. The one-step
disturbances correspond, indexed by the distance `d` along the order:

- adding the edge `{u, v}` (distance `d = v - u ≥ 2`) to `P_{m+1}` creates **exactly one** cycle, of
  length `d + 1`;
- reversing the arc `(i, j)` (distance `d = j - i`) in `TT_m` creates **exactly `2^{d-1} - 1`** cycles,
  `C(d-1, L-2)` of each length `L = 3, …, d+1`, and none for `d = 1`.

*Proof.* A cycle through the reversed arc `j -> i` continues by an increasing path from `i` to `j`
through a nonempty subset of `{i+1, …, j-1}` (the direct arc `i -> j` no longer exists); the subset of
size `L-2` gives a cycle of length `L`. The path case is the unique cycle of the unicyclic graph. ∎
(FINITE-EXACT: all arcs/non-edges for `m = 3..6`.)

So the owner's phrase is exact on both sides and the dictionary is "one added edge = one reversed
arc = one distance `d`"; the tournament side carries the binomial multiplicity (its cycles are the
subsets of the interval), the tree side the single cycle. THM-220 is the extreme case `d = m-1`
(`2^{m-2}` cycles; `H = 2^{m-2} + 1`).

## 4. Proposition 3 (PROVED): the hinge is `K_4`, and the half is the handshake

For trees on `t` vertices the Erdős–Sós threshold is "more than `(t-2)n/2` edges". The extremal
configuration is the disjoint union of cliques `K_{t-1}` (one vertex too small to host the tree); its
edge density `(t-2)/2 = e(K_{t-1})/(t-1)` is the **mean score of a tournament on `t-1` vertices**,
because tournaments on `t-1` vertices are the orientations of the same clique. At `t = 5` the hinge
is `K_4`: the owner's two objects are the trees that do not fit into `K_4` and the orientations of
`K_4`. Both "halves" are the same parity fact:

- a `(t-2)`-regular graph on `n` vertices exists iff `(t-2)n` is even (handshake), so `(t-2)n/2` is a
  whole number exactly when the degree-extremal graph exists;
- a regular tournament on `m = t-1` vertices exists iff `m` is odd, i.e. iff `(m-1)/2 = (t-2)/2` is whole.

At `m = 4` the mean score `3/2` is a half: there is no regular 4-tournament, the unique near-regular one
is `S_4 = (1,1,2,2)`, the unique acyclic one is `TT_4`, and the converse pair sits between.

**The dictionary swaps whole and half between its two sides.** A reflection-fixed vertex of `TT_m`
(score exactly `(m-1)/2`) exists iff `m` is odd; a central *vertex* of `P_{m+1}` exists iff `m` is even.
So the tournament has a central vertex exactly when the tree has a central edge. (At `m = 4` the
reflection `ρ` of THM-4472, `ρ(x) = s_o + s_e - x`, is the point reflection about the half-integer
`2i - 1/2`, the mean of the AM-fair pair of THM-4470: the quadruple is in the "half" mode.)

## 5. Proposition 4 (PROVED): Erdős–Sós at `t = 5` with exact extremal numbers, and the two mode periods

For `n ≥ 4`:

    ex(n, K_{1,4}) = ⌊3n/2⌋,
    ex(n, fork)    = ex(n, P_5) = 6⌊n/4⌋ + C(n mod 4, 2).

*Proof.* Star: a graph has no `K_{1,4}` iff its maximum degree is `≤ 3`; near-regular graphs of degree
3 exist for all `n ≥ 4`. Path: Erdős–Gallai 1959 (bound) and Faudree–Schelp 1975 (exact value
`(l-2)n/2 - r(l-1-r)/2`, `n ≡ r mod (l-1)`, here `l = 5`). Fork (**structure lemma**): let `G` be
fork-free and `c` a vertex of degree `≥ 3`. If a neighbour `a` of `c` had a neighbour `e ∉ N[c]`, then
`c` with `a` and two further neighbours, plus `e`, is a fork; so every neighbour of `c` has all its
neighbours in `N[c]`, and `N[c]` is the whole component of `c`. If `deg c ≥ 4`, an edge `ab` inside
`N(c)` gives the fork `c; a, d, e; b`, so `N(c)` is independent and the component is a star. If
`deg c = 3` the component has 4 vertices. Hence every component of a fork-free graph is a star, a path,
a cycle, or a graph on `≤ 4` vertices, and the maximum edge count is attained by `⌊n/4⌋` copies of
`K_4` plus a clique on the remainder (a component of size `s ≠ 4` has at most `s` edges, at most `1`
if `s = 2`, so trading a `K_4` for anything else loses). ∎ (FINITE-EXACT: the lemma agrees with direct
subgraph search on all 1252 graphs with `≤ 7` vertices and on 3000 random graphs with 8–12
vertices; the values `7/6/6`, `9/7/7`, `10/9/9` at `n = 5, 6, 7` by exhaustive search over the atlas.)
In particular the Erdős–Sós bound holds for all three trees on 5 vertices (fork: also McLennan 2005,
diameter `≤ 4`), with equality `ex = 3n/2` iff `n` even (star) resp. iff `4 | n` (path, fork).

**The two modes.** Write the deficit `δ_T(n) = 3n/2 - ex(n, T)`:

| `n mod 4` | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| `δ_star` | 0 | 1/2 | 0 | 1/2 |
| `δ_path = δ_fork` | 0 | 3/2 | 2 | 3/2 |

- `δ_star` has period 2: `ex(n, star)` satisfies the recurrence with characteristic polynomial
  `(x-1)^2 Φ_2` — the owner's whole/half alternation is exactly the `Φ_2 = x+1` factor of the
  repo's recursion-mode grading (the ODD mode's parity character); for even `t` the threshold
  `(t-2)n/2 + 1` is in the EVEN mode `(x-1)^2` with no `Φ_2`.
- `δ_path = δ_fork` has period 4 = `t - 1`: `ex(n, P_5)` lives in the cyclotomic factorization of
  `x^4 - 1`, characteristic polynomial `(x-1)^2 Φ_2 Φ_4`.
- In the repo's Mode A / Mode B language: the forcing threshold `E*(n) = ⌊3n/2⌋ + 1` is **Mode-B
  exact** (`E*(n+2) = E*(n) + 3`) and **Mode-A alternating** (`E*(n+1) - E*(n) = 1, 2, 1, 2, …`): the
  "+2 is the natural recursion" of `the-two-staircases.md` is the statement that the half disappears
  under `n -> n+2`.

**Otter's recursion** is the tree-side instance of the same alternation (CITED Otter 1948; VERIFIED
`n ≤ 10`): `t_n = r_n - (1/2) Σ_{i} r_i r_{n-i} + [n even] r_{n/2}/2`, the half-term present exactly for
even `n` (a bicentroid needs two halves of `n/2` vertices). For the three 5-vertex trees the fork is
the one bicentral tree (diameter 3), which is the tree-side shadow of "the converse pair is the
one-defect class".

## 6. What reaches Collatz (PROVED extension of THM-4472; typed)

**Proposition 5 (PROVED, elementary).** On the AM-fair quadruple `V = {i < 2i-1 < 2i < 3i-1}` of
THM-4470/4472 (`i ≥ 2`, sheet `3n+1`), orient every non-Collatz pair smaller -> larger and orient the
two Collatz-graph edges `{i, 2i}` (the `×2` edge) and `{2i-1, 3i-1}` (the `×3+1` edge) independently
"up" (smaller -> larger) or "down". The four orientation regimes are the four classes:

| `×2` edge | `×3+1` edge | tournament | `u` | reading |
|---|---|---|---|---|
| up (`i -> 2i`) | up (`2i-1 -> 3i-1`) | `TT_4` | 0 | both moves expanding |
| down (`2i -> i`) | up | `(0,2,2,2)` | 1 | the forward map `T` (THM-4472 forward) |
| up | down (`3i-1 -> 2i-1`) | `(1,1,1,3)` | 1 | the inverse tree (THM-4472 backward) |
| down | down | `S_4` | 2 | both moves contracting |

The Klein group is `(Z_2)^2 =` {reverse the time of the `×2` move} × {reverse the time of the `×3+1`
move}; THM-4472's converse (= negation ∘ total time reversal) is the diagonal element; the merged label
`u` is the **number of size-contracting moves**, which is the only datum surviving the loss of time's
arrow. *Proof:* the `×2` edge joins positions 1 and 3 of the order, the `×3+1` edge positions 2 and 4,
so the regimes are exactly `TT_4, g_L TT_4, g_R TT_4, g_L g_R TT_4` of Proposition 1. ∎

Consequences and non-consequences.
- THM-4472's "the forward class is a 3-cycle over a sink, the inverse tree a source over a 3-cycle"
  is the `u = 1` layer; `S_4` never appears there because it needs the two moves to run in opposite
  time directions (one forward step `2i -> i` and one backward step `3i-1 -> 2i-1`): it is the
  "descent" orientation, both arcs toward the smaller number.
- The merged node `[±]` is therefore "the dynamics with the arrow of time forgotten"; the tree it
  corresponds to (the fork, `u = 1`, the one bicentral 5-vertex tree) is the picture of one expanding
  and one contracting move — the `log 2` against `log 3` tension of THM-4482 — but this is a
  dictionary of labels, not a transfer of any inequality. Nothing about Collatz orbits follows.
- THM-4472's "every odd multiplier `q ≥ 3` gives the same diamond" persists: the regime table does
  not see `q`.

**The forcing threshold as a Collatz composite (SYNTAX-LEVEL, all `n ≤ 10^6` checked).**
`E*(n) = ⌊3n/2⌋ + 1 = (3n+1)/2` for odd `n` and `= 3(n/2) + 1` for even `n`: the two orders of
composing `x -> 3x+1` and `x -> x/2`. Both sides are "`3/2` and rounding"; no structure crosses
(MISTAKE-230–235: shared syntax is not a bridge). Recorded because the owner asked; typed as an
elementary identity.

**Two period-2 alternations with different mechanisms (DISTINCT MECHANISMS).** The Erdős–Sós half is
`2 ∤ 3n`, i.e. the parity of `n`. The Collatz inverse step `n = (2^k m - 1)/3` is whole iff
`2^k m ≡ 1 (mod 3)`: `k` even for `m ≡ 1`, `k` odd for `m ≡ 2 (mod 3)`, never for `3 | m` — period
`2 = ord_3(2)`, a multiplicative order, not a divisibility. Same period, different prime doing the work.

**The "fractal ternary tree" (NUMEROLOGY, with a decisive obstruction).** The repo's ternary object is
the inverse Collatz tree's three-row automaton (`collatz_mod6_20260917_reverse_tree_pieces.md`,
Theorem 1.1): children of `u` rotate through the rows `1 -> 5 -> 3 -> 1 (mod 6)` with exact period
`3 = ord_9(4)`, every third child a leaf. Its transition structure is a directed 3-cycle (eigenvalues
the cube roots of unity). The merged 4-metagraph is a symmetric triangle with two loops (eigenvalues
`-2, 2, 6`). No graph map carries one onto the other; the shared datum is the numeral 3. The honest
tree-side object with a genuine `3` is Proposition 1's `u ∈ {0,1,2}` (and `H = 1, 3, 5`), which is a
`(Z_2)^2`-orbit count, not a ternary branching.

## 7. Verdicts on the owner's four claims

| claim | verdict | carrier |
|---|---|---|
| 3 trees on 5 vertices ↔ 3 converse-merged 4-classes | STRUCTURAL as the end-defect `V_4` dictionary (Prop. 1, all `m`); EXPLAINED COINCIDENCE as a census equality (exact for `m ≤ 4`, fails `10 ≠ 6` at `m = 5`) | `u = #end defects = m - diameter`; reflection = automorphism (trees) vs anti-automorphism (tournaments) |
| "structures that create a cycle when an edge is added" | EXACT on both sides (Prop. 2) | one added edge / one reversed arc at distance `d`; one cycle vs `2^{d-1}-1` |
| "deeply aligns with fractal ternary trees" | NUMEROLOGY (spectral obstruction, §6) | none |
| `(t-2)n/2` whole/half ↔ recursion modes | DICTIONARY, exact (Prop. 3–4): `Φ_2` factor; Mode-B exact / Mode-A alternating; handshake = regular-tournament parity; the hinge `K_{t-1}` | `(t-2)/2 = (m-1)/2 = e(K_m)/m` |
| weave into Collatz | PROVED extension of THM-4472 (Prop. 5: `u` = number of contracting moves; `S_4` = both-contracting mixed-time regime); everything else SYNTAX-LEVEL or DISTINCT MECHANISMS | the AM-fair quadruple |

## 8. Open items and cheapest next tests

- **Sumner for `t = 5`** (greenfield in the ledger): does every tournament on 8 vertices contain every
  oriented tree on 5 vertices (27 oriented trees, 6880 host classes)? The regular 7-tournaments are
  the extremal hosts for the out-star, by the same handshake mechanism as §4. A nauty-free class
  generator for `n = 8` is the only cost.
- **Deficit periods for larger `t`:** is `δ_T(n)` eventually periodic with period dividing
  `lcm(2, t-1)` for every tree `T` with `t` vertices where Erdős–Sós is known (stars: 2; paths:
  `t-1` by Faudree–Schelp)? The 6 trees on 6 vertices are the first test (`t - 1 = 5` odd, so no
  half ever appears: the alternation is a `t`-odd phenomenon).
- **Boundary layer of `M_m`:** the three merged nodes `{TT_m, [C_3 ⊕ TT_{m-3}], C_3 ⊕ TT_{m-6} ⊕ C_3}`
  of Proposition 1 are 3 of the 10 merged 5-classes and 3 of the 34 merged 6-classes; whether the
  remaining merged classes admit any natural tree label is not asked here (no candidate map).
- Audit owed: Propositions 1–5 (elementary; the finite checks reproduce THM-584's `n = 4` spectrum).

## 9. Reproduction

```
python3 04-computation/experiments/trees5_tournaments4_erdos_sos_20261005.py > 05-knowledge/results/trees5_tournaments4_erdos_sos_20261005.out
```

~40 s; networkx 3.4.2, Python 3.10. Universe: all labelled tournaments on `≤ 6` vertices; all
graphs on `≤ 7` vertices (networkx atlas); `nonisomorphic_trees(n)`, `n ≤ 10`; random `G(n,p)`
controls `n = 8..12`; `n ≤ 10^6` for the threshold identity; odd `m < 10^5` for the inverse-step
parity. Positive controls: OEIS A000055, A000568, A002785, A000081 values; THM-584's spectrum
`{-2,2,6}`; Rédei parity. Hostile controls: the census split at `m = 5`; the tree-rotation graph vs
the merged metagraph; the spectral obstruction. SHA-256 (raw bytes): script
`59c48bac1f6f1b8f11a06224b461e9ece291e4583211c01d741154a462e1a45d`, output
`40c8539e288dae15afe65617768d1376765cd4d8f25404ed373a8fb9c8450fde`.

References (CITED). P. Erdős, T. Gallai, *On maximal paths and circuits of graphs*, Acta Math.
Hungar. 10 (1959). R. J. Faudree, R. H. Schelp, *Path Ramsey numbers in multicolorings*, J. Combin.
Theory B 19 (1975) (exact `ex(n, P_l)`). A. McLennan, *The Erdős–Sós conjecture for trees of
diameter four*, J. Graph Theory 49 (2005) 291–301. R. Otter, *The number of trees*, Ann. Math. 49
(1948). OEIS A000055, A000568, A002785, A059735 (complementary pairs of tournaments).
