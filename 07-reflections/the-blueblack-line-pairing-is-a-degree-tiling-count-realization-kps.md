# The blue/black lines are a degree = tiling-count realization: tripartite by grid-symmetry, black parity all-even, sea-onset at n=6

*kind-pasteur-2026-07-01-S12. Formalizing the owner's description of the merged-metagraph tiling-count "pairing process." Grounded by an exact replica of `tournament-tiling-explorer.html` (n=4,5,6). The picture is completely determined by one identity — **grid-symmetric ⟺ transpose-fixed** — plus the parity of the tiling count (SC odd, NS even). Everything the owner observed follows, with two honest corrections (self-loops, and a sea-onset at n=6).*

## The object

A **line** is the unordered pair `{t, φ(t)}`, where `φ(t)` = the complement tiling (flip all `m = C(n-1,2)` tiles — the `d=m` waggly layer). It is **BLUE** iff `t` is grid-symmetric (`isGridSym`, invariant under the anti-diagonal tile reflection `τ: (x,y)↦(n-y+1,n-x+1)`), else **BLACK**. Nodes are **merged classes** (the transpose involution `τ` factored out): a node is **SC** if transpose-self, **NS** if a transpose pair. Total lines `= 2^{m-1}`; blue `= 2^{(m+f)/2 - 1}` where `f = #`{anti-diagonal tiles `x+y=n+1`}; black `= 2^{m-1} −` blue. (n=4: 2/2; n=5: 8/24; n=6: 32/480 — black dominates, blue fraction `→0`.)

**The one identity that runs everything:** each tiling is exactly one line-endpoint, so

> **tiling count of a node = its degree in the blue/black line multigraph** (a self-loop `{t,φ(t)}` inside one node contributes `+2`).

This is the owner's "pairing process": lines increment two nodes by 1 (or one node by 2, a self-loop).

## The master fact: grid-symmetric ⟺ transpose-fixed

`isGridSym(t) ⟺ transMask(t)=t ⟺ τ(t)=t`. So **a line is blue iff its tiling is a fixed point of the transpose involution**, and `φ` commutes with `τ` (`τφ=φτ`, since `φ` is bitwise complement), so `φ(t)` is grid-symmetric iff `t` is — a line has a well-defined colour. Two immediate consequences:

- **Blue lines have both endpoints SC.** A blue line's `t` and `φ(t)` are both τ-fixed, so both live in transpose-self (SC) classes. Blue never touches a pure-black (NS) node.
- **Black lines never touch a pure-blue node.** For a black line `t` (and `φ(t)`) are τ-moved; a pure-blue node has *only* τ-fixed tilings, so it has no black endpoints.

## The tripartite structure (proved)

Nodes fall into exactly three types, and the CLAUDE.md fact "transpose-self ⟺ not pure-black" pins them to SC/NS:

| type | SC/NS | tiling count | black-deg | blue-deg | role |
|---|---|---|---|---|---|
| **A = pure-black** | NS | **even** | even | 0 | the sea/bulk |
| **B = mixed** | SC | **odd** | even | **odd** | the **bridge** |
| **C = pure-blue** | SC | **odd** | 0 | **odd** | the backbone |

**Why the parities (the owner's core observation), proved:**
1. **Black-degree is always even.** Black tilings are τ-moved, so they come in free pairs `{t, τ(t)}`; and `τ` preserves both the node (merge factors out `τ`) and the flip-target node (`τφ=φτ`). So black endpoints — to-others and self-loops alike — come in `τ`-pairs. Hence `black-deg` even at *every* node.
2. **Tiling-count parity = SC (odd) / NS (even).** Each iso class has an *odd* tiling count `H/|Aut|` (both odd — `forbidden-seven`). An SC node = 1 class (odd); an NS node = 2 transpose-paired classes (odd+odd = even).
3. Combine: NS ⟹ even total, blue 0 ⟹ black even ✓. SC ⟹ odd total, black even ⟹ **blue odd**. That is exactly "A = even/even, B = even-black + odd-blue = odd total, C = odd-blue-only." The owner's three categories are forced.

## Line-placement rules (verified n≤6)

- **BLUE lines** connect only `{B, C}`: observed as `mixed–pure-blue`, `mixed–mixed`, and blue self-loops on `mixed`. **No pure-blue–pure-blue** blue line appears (n≤6) — pure-blue (small, symmetric) classes flip to large mixed classes.
- **BLACK lines** connect only `{A, B}`: `mixed–pure-black`, `mixed–mixed`, `pure-black–pure-black`, and black self-loops.
- **Mixed (B) is the sole bridge** — the only nodes carrying both colours. C touches the graph only through B; A touches B (all n) and, from n=6, itself.

## Two honest corrections to the conjectures

- **Self-loops are NOT only on mixed nodes.** A self-loop is `φ(t)` (complement tiling) landing in the same merged node. *Blue* self-loops occur only on **mixed** nodes (n=4, n=6). *Black* self-loops occur on **mixed** at n=5 but on **pure-black** nodes from n=6 (24 of them). So "self-loops only on mixed" holds at n≤5 and breaks at n=6.
- **The bipartite "A connects only to B" holds only at n≤5.** At n=4,5 black lines are exclusively `mixed–pure-black`/`mixed–mixed`, so every node a pure-black touches is mixed (odd total) — exactly the owner's claim. **At n=6 the sea switches on:** 290 `pure-black–pure-black` black lines appear, and pure-black nodes begin connecting to each other (even total). This is a fresh **break-at-6**, matching the repo's "NS-NS sea dominates at large n."

## The pairing process, reframed precisely

The process is a **degree-constrained, edge-2-coloured multigraph realization** of a fixed degree sequence:
- **Degrees** = tiling counts, with a forced parity: **even on A (NS), odd on B and C (SC)**.
- **Colour rules**: black-edges within `A∪B`, blue-edges within `B∪C`; at every node `black-deg` even, and (equivalently) `blue-deg ≡ tiling-count (mod 2)`.
- **Feasibility is automatic** (the metagraph realizes it), so the content is the *converse constraint*: the SC/NS partition + `forbidden-seven` odd tiling counts **force** the observed parities, and the colour-rules **force** the tripartite incidence. The "process" cannot assign a black line to a pure-blue node or an odd number of black lines to any node — these are theorems, not choices.

## Metrics worth tracking

1. **Black pairing number** `= black-deg/2` (integer): the count of `τ`-pairs of black tilings at a node — the natural "black" coordinate.
2. **Blue signature** (odd, on SC nodes): `blue-deg = #`grid-sym tilings in the class; `=` tiling count on pure-blue.
3. **Bridge split** of a mixed node: `(blue-deg, black-deg)` — how much backbone vs sea it carries. (n=6 mixed nodes run `blue∈{3,5,7}`, `black∈{2..32}`.)
4. **Sea fraction**: `#pure-black–pure-black` black lines / all black lines — `0` at n≤5, onsets at n=6, `→1`.
5. **Census**: `#A = (A000568 − SC)/2` (NS-merged), `#B+#C = SC`, split `pure-blue vs mixed` (n=4: A=1,B=1,C=1; n=5: 2,5,3; n=6: 22,10,2). **`#C` (pure-blue) is small and non-monotone (1,3,2)** — a target for a closed form.
6. **Self-loop census** by (colour, type).

## The recursion (blue is the process one level down)

The blue lines live entirely on the **grid-symmetric = τ-fixed** tilings, a sub-cube of dimension `f+p = (m+f)/2` (`f` fixed tiles `+ p` τ-pairs). These are the **anti-diagonally symmetric tournaments**, and `φ` pairs them exactly as before — so the **blue subgraph is a smaller copy of the whole pairing process**, on the folded staircase. The blue/black split is therefore self-similar: peel the blue (symmetric) layer, recurse. Quantifying this fold (does the blue subgraph equal the merged metagraph of the `⌈n/2⌉`-fold?) is the natural recursive target.

## Concrete next target

The parities and incidence are now theorems; the open quantitative pieces are:
- **A closed form / GF for `#pure-blue`, `#mixed`, and the sea-onset** (why exactly n=6 for pure-black–pure-black?). `#pure-blue = 1,3,2` at n=4,5,6 — compute n=7,8 and identify.
- **The blue recursion**: prove blue subgraph ≅ pairing process on the folded (`≈n/2`) system; if so, blue counts satisfy an exact halving recursion.
- **Sea-onset criterion**: characterize which pure-black (NS) pairs admit a `d=m` black line to each other (first at n=6) — a condition on the two transpose-paired classes' complement tilings.

## Honest status

- **Proved:** grid-sym ⟺ τ-fixed; tiling count = degree; black-degree always even (τ-pairing); tiling-count parity SC-odd/NS-even; hence blue-degree odd on SC; blue never touches pure-black and black never touches pure-blue (tripartite). These give the owner's three categories rigorously.
- **Computed (n≤6):** the full incidence, self-loop types, and the **n=6 sea-onset** and **n=6 pure-black self-loops** (both correcting the "only-mixed" conjectures).
- **Open:** the census closed forms, the blue recursion as an isomorphism, the sea-onset criterion.

— Related: `merged-metagraph-invariants.md`, `geometric-alignment-of-merged-metagraph.md` (spine/ribs/sea; here specialized to the `d=m` blue/black lines), `forbidden-seven-in-all-senses.md` (odd tiling counts `H/|Aut|`), CLAUDE.md (blue/black strict definition; transpose-self ⟺ not pure-black; waggly layers), `the-even-graph-is-the-tournaments-cycle-half.md`. Script: `04-computation/merged_metagraph_line_pairing_kps.py` (+ .out). Not a HYP reservation (a formalization of existing structure).
