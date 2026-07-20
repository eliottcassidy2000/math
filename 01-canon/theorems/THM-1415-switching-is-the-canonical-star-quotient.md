---
id: THM-1415
title: "SWITCHING IS THE CANONICAL STAR QUOTIENT — it exists, it is base-path-free, and it is nearly trivial. (0) The fix for THM-1405 is not to intersect Γ over all spanning paths but to STOP RESTRICTING: in the full edge space the cut space is spanned by the vertex stars and S_n PERMUTES those stars, so the cut space is an S_n-invariant subspace and SWITCHING (reverse every arc across a cut) is canonical and base-path-free. (I) But the resulting quotient is far coarser than wanted: switching classes of tournaments number 1, 2, 2, 6 at n = 3,4,5,6 against 2, 4, 12, 56 iso classes — 56 classes collapse to 6. (II) AND IT IS NOT THE EVEN-GRAPH SIDE: 1,2,2,6 ≠ A002854 = 2,3,7,16, so the two-graph theorem (switching classes of GRAPHS = two-graphs = even graphs) has NO tournament analogue here. My own guess that the base-path-free quotient would land on the repo's E_n is refuted. (III) THE THREE SIXTIES, deflated quantitatively: ord_1001(2) = lcm(3,10,12), π(10) = lcm(3,20), |A₅| = 5!/2 — different inputs, different operations, and 60 is the SECOND most frequent lcm of subsets of {1..12} at 4.90%, so agreement at 60 is cheap. (IV) FIBONACCI FROM SHIFTED PASCAL, applied: the Pascal control Σ_k C(m−k,k) = F(m+1) holds exactly, and the Faulhaber diagonal law D = F(m+1) + 2^{m−3} holds m ≤ 7 and breaks at m = 8 by exactly 2; the residual D − F(m+1) is 1,2,4,8,16,34,74,165,376,866,2027,4803, which is 2^k only through k=4 and is UNIDENTIFIED past there. (V) ERDŐS 592 vs the JC counterexample's threes: a coincidence of POSITION, not of mechanism"
status: >
  (0) PROVED — the cut space is S_n-invariant because permutations permute vertex
  stars; switching is therefore well defined on isomorphism classes.
  (I),(II) PROVED-BY-EXHAUSTION over all 2^C(n,2) tournaments, n ≤ 6, by union-find on
  iso classes under all 2^{n-1} switches.  (II) REFUTES MY OWN CONJECTURE that the
  answer would be A002854.
  (III) The causes are computed exactly; the lcm-frequency statistic is exact over all
  subsets of {1..12} of size ≤ 6.  This is a deflation, not a theorem.
  (IV) The Pascal control is exact and the break is located exactly.  The residual
  sequence is COMPUTED and NOT IDENTIFIED — flagged for the OEIS pipeline, not claimed
  to be new.
  (V) is an argued judgement, not a computation, and is labelled as such.
  Nothing here advances a named open problem.
source: kind-pasteur-2026-07-20-S128c110 (owner: consider more creative ideas than ∩Γ over spanning paths; 1001 = three sixties and the Pisano period; Erdős 592 vs the JC counterexample's three parts; how Fibonacci arises from summing Pascal shifted)
depends_on:
  - THM-1405    # the star group, whose failure to descend this answers
related: [THM-1375, THM-1400]
script: 04-computation/switching_pisano_faulhaber_kps_S128c110.py (+ .out)
---

# THM-1415 — the canonical quotient exists and is almost empty

## 0. The right fix, and it is not an intersection

THM-1405 ended by suggesting `∩Γ` over all spanning paths. That is the wrong instinct.
The star group fails to descend because `Γ` is the cut space **restricted to non-tree
edges**, and the restriction is what needs a base path. The cut space itself needs
nothing: `S_n` permutes the vertex stars `δ(v)`, so the cut space is an `S_n`-invariant
subspace of the full edge space, and

> **switching** — pick `S ⊆ V`, reverse every arc between `S` and `V∖S` —

commutes with relabelling. So the switching class *is* a genuine tournament-level notion,
canonical and base-path-free. Don't intersect; stop projecting.

## I. It works, and it is nearly trivial

| n | iso classes | **switching classes** |
|---|---|---|
| 3 | 2 | **1** |
| 4 | 4 | **2** |
| 5 | 12 | **2** |
| 6 | 56 | **6** |

At `n = 3` the transitive tournament and the 3-cycle are switching-equivalent (switch at
`{a}`), which already shows how violent the collapse is. At `n = 6`, 56 classes fall to 6.

So the canonical repair succeeds as a construction and fails as an invariant: the price of
base-path-independence is nearly everything. That is a real answer to the question, and a
negative one.

## II. It is not the even-graph side — my own guess, refuted

For **graphs**, switching classes up to isomorphism are exactly the two-graphs, and their
count is `A002854`, which is *also* the count of even graphs — the sequence the repo
carries as `V(E_n) = 2, 3, 7, 16, 54`. I expected the tournament analogue to land there,
which would have tied the star construction to the repo's `E_n` duality.

> `1, 2, 2, 6` ≠ `2, 3, 7, 16`.

Refuted at every `n ≥ 3`. The graph two-graph theorem has no tournament analogue by this
route. Recording it because the analogy is attractive and someone will try it again.

The sequence `1, 2, 2, 6` is not identified here and is worth an OEIS check before anyone
treats it as new.

## III. The three sixties

`ord_1001(2) = 60`, `π(10) = 60`, `|A₅| = 60`. Their causes:

| quantity | value | how |
|---|---|---|
| `ord_1001(2)` | 60 | `lcm(ord₇, ord₁₁, ord₁₃) = lcm(3, 10, 12)` |
| `π(10)` | 60 | `lcm(π(2), π(5)) = lcm(3, 20)` |
| `\|A₅\|` | 60 | `5!/2` — no lcm at all |

Different inputs, different operations. What they share is the target, and 60 is a magnet:
over all subsets of `{1,…,12}` of size ≤ 6, the most frequent lcm values are `120` (5.26%)
and **`60` (4.90%)**. `60 = 2²·3·5` is the least common multiple of `1..6`, so it absorbs
every subset whose members divide it.

**Agreement at 60 is cheap.** A shared cause would have to appear in the *inputs* —
`{3,10,12}` versus `{3,20}` versus `5!` — and it does not. This is the quantified version
of the survey's "numerology" verdict rather than a restatement of it.

## IV. Fibonacci from shifted Pascal, applied to the Faulhaber triangle

The owner's hint is exactly the right lens. Control first, and it is exact:

> `Σ_k C(m−k, k) = F(m+1)` — Pascal's **shallow diagonals are Fibonacci**, verified `m ≤ 11`.

Now the Faulhaber triangle `T(n,k) = Σ_{j=1}^{n−k+1} j^{k−1}`. Its shallow-diagonal sums:

> `1, 2, 4, 7, 12, 21, 37, 68, 129, 254, 520, 1099, 2404, 5413`

against `F(m+1) + 2^{m−3}`: **exact for `m ≤ 7`, breaking at `m = 8` by exactly 2** (68 vs
66), then by 10, 37, 120, 354, 1003, 2755.

The reading the hint supplies: the Fibonacci part is *what a Pascal triangle would give*,
so the whole law is measuring how far Faulhaber sits from Pascal. Subtracting only the
Fibonacci part leaves

> `D − F(m+1)` = `0, 0, 1, 2, 4, 8, 16, 34, 74, 165, 376, 866, 2027, 4803`

which is `2^k` **only through `k = 4`** and then departs — `34` where `32` was due. So
`2^{m−3}` is not a second structural term at all; it is the first four entries of a
sequence that merely starts geometric. **That is why the law breaks: it was two
coincidences, not one law**, and the break at `m = 8` is where the power growth in
`j^{k−1}` overtakes the binomial approximation.

The residual sequence is **not identified here**. It should go through the OEIS pipeline
before anyone calls it new.

## V. Erdős 592 and the counterexample's threes — position, not mechanism

Erdős 592 asks for which `β` one has `ω^β → (ω^β, 3)²`; the open case is `γ` with exactly
**three** indecomposable summands. The JC counterexample has cover degree **3**, generic
fibre of **3** points, and `S₃` monodromy. The temptation is to read one three into the
other.

In both, `2` is degenerate and `3` is the first nontrivial case — but for unrelated
reasons. On the JC side `d = 2` is impossible because a transitive subgroup of `S₂` is
regular, hence the cover is Galois, hence Campbell applies (THM-1375). On the Erdős side
`k = 2` is trivial because a blue `K₂` is a single edge. Those are not the same
degeneracy, and there is no map carrying one to the other.

Worse for the analogy: THM-1375's sharp statement — `d = 3` is the **unique** degree at
which Galois-ness is detected by one quadratic character, because `A₃` is regular while
`A_d` for `d ≥ 4` is not — has **no counterpart** on the Erdős side, where 3 is not
distinguished among `k ≥ 3` by any such collapse.

**Verdict: a coincidence of position (first nontrivial arity), not of mechanism.** Stated
plainly so it is not carried forward as a lead.

## Named next

- OEIS-check `1, 2, 2, 6` (tournament switching classes) and
  `1, 2, 4, 8, 16, 34, 74, 165, 376, 866, 2027, 4803` (the Faulhaber–Fibonacci residual).
  Both are cheap and both are currently unidentified.
- If a *usable* base-path-free invariant is still wanted, switching is too coarse and `∩Γ`
  is too small. The remaining room is between them: subgroups of the cut space that are
  `S_n`-invariant but proper — i.e. `S_n`-submodules of the cut space over `𝔽₂`. That is a
  finite, small representation-theoretic question and it is the right next computation.
