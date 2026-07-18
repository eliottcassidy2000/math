---
id: THM-1121
title: The r=6 finite horn collapses to a 35-point universal weighted atlas — total witness weight 505, while every candidate killer has capacity at most 84, so six killers have capacity at most 504 and cannot cover the atlas
status: PROVED for the complete finite branch 92 <= k_i < 333 by a dependency-free exact integer verifier; this removes THM-1102's estimated 3.64e12-sextuple enumeration wall. It does NOT by itself prove the unbounded r=6 branch: the max-T calculation in THM-1102 scanned a width-16 near-bottom window, so a uniform tail theorem is still required before declaring the whole r=6 clustered case closed
source: codex-2026-07-18-S67 (r6 structural-compression subtask)
depends_on:
  - THM-1041   # rational small-modulus witness criterion
  - THM-1102   # proposed finite range and the enumeration wall replaced here
  - THM-1098   # why the fixed atlas must remain height-bounded
related:
  - THM-1111   # MST overlap prune/dedupe negative; its named-next asks for this dual relaxation
related:
  - THM-1111   # pairwise MST prune is strong but nonterminal; this weighted dual closes its gap
script: 04-computation/r6_universal_weighted_atlas_codex_S67.py
output: 05-knowledge/results/r6_universal_weighted_atlas_codex_S67.out
---

# THM-1121 — a universal weighted atlas closes the r=6 finite horn

## The exact statement

Put

\[
 h(q)=\left\lceil\frac q{14}\right\rceil,
 \qquad
 \ell_q(x)=\min(x\bmod q,\;q-(x\bmod q)).
\]

There is an explicit weighted set `A` of 35 rational-time obligations `(q,a)`, with
`26 <= q <= 112`, having the following three properties.

1. Its total integer weight is
   \[
   w(A)=505.
   \]
2. Every obligation is safe for every possible core speed:
   \[
   \ell_q(pa)\ge h(q)
   \quad(1\le p\le12,\ (q,a)\in A).
   \]
3. Every integer candidate killer in the complete finite range has capacity at most 84:
   \[
   w\{(q,a)\in A:\ell_q(ka)<h(q)\}\le84
   \quad(92\le k\le332).
   \]

Consequently, for **any** six integers `k_1,...,k_6` in `[92,332]`, repetition allowed,

\[
 w\left(\bigcup_{i=1}^6
 \{(q,a)\in A:\ell_q(k_i a)<h(q)\}\right)
 \le \sum_{i=1}^6 84=504<505.
\]

Some `(q,a)` in `A` is therefore safe for all six killers and all speeds `1,...,12`.
At the rational time `t=a/q`, every one of those speeds has

\[
 \|vt\|=\frac{\ell_q(va)}q
 \ge\frac{h(q)}q\ge\frac1{14}.
\]

This proves, more strongly than needed:

> If `P` is any subset of `{1,...,12}` and `K` is any multiset of at most six integers
> in `[92,332]`, then `P union K` has a rational lonely time with denominator at most 112.

For THM-1102's seven-speed core, `max(P)>=7`; its killer convention
`k>13 max(P)` therefore gives `k>=92`. Thus this is exactly the entire proposed
`KB=333` finite horn. Distinctness and the divisibility-covering conditions `2,...,14`
are not used.

## The 35 weighted obligations

The common denominator of the discovery LP dual was 84. After clearing it, the exact
integer table is:

```text
(26,2):81, (27,2):7, (28,2):3, (28,13):7,
(40,3):10, (41,3):3, (41,19):5, (42,13):1,
(53,4):16, (55,4):2, (56,13):21,
(68,5):3, (68,21):21, (69,5):1, (69,16):12, (70,27):20,
(79,6):3, (81,25):5, (82,19):6, (83,6):1, (83,32):5,
(84,13):8, (93,7):9, (94,7):13, (97,7):1, (97,30):2,
(98,15):6, (105,8):4, (106,49):46, (107,33):10,
(109,42):49, (110,17):25, (111,8):1, (111,17):64,
(112,43):34.
```

The proof does not trust a floating-point optimizer. The optimizer was only a discovery
device. The frozen verifier uses integer modular arithmetic to check all
`35*12` core inequalities and all `241*35` killer incidences, recomputes the total 505,
and asserts the exact maximum capacity 84. Its source and frozen-output SHA-256 hashes are

```text
2232e1b115a0f5d4c2c8dab23aff2aa89965a358257480fb6263f242f6deff6d
3058a7b59ab0a6da62da27a3090fae3b7581d372a47b9141723872c56cc3bc3f
```

respectively.

## Structural reading: change the vertices

THM-1102 used one universe of safe `(q,a)` pairs for each of 792 seven-speed cores and
then enumerated sextuples of killer masks. The compression comes from challenging both
choices:

- **Vertices are proof obligations, not runners.** Use rational-time obligations `(q,a)`
  as the left vertices and killers as incidence hyperedges.
- **Make the atlas universal.** Requiring a point to be safe for all twelve possible core
  speeds makes it valid for every subcore, removing the 792-way branch completely.
- **Enlarge the chart before enumerating families.** Moving from the old ad hoc window
  `q<=40` to a still-small `q<=112` creates enough independent obligations that a weighted
  dual certificate exists.
- **Weight residue ownership.** Cardinality was the wrong observable: THM-1102's
  unweighted `sum frac >= 1` prune nearly dissolved at six killers. The dual weights make
  every killer have capacity 84 while retaining total demand 505. The one-unit inequality
  `6*84<505` is the exact obstruction.

Incoming THM-1111 independently sharpened the old cardinality prune by subtracting a
maximum-spanning-tree of pair intersections. It found that random candidates almost all
die, but adversarial per-core masks still leave the MST upper bound 36 bits above coverage;
it explicitly named a **fractional-relaxation bound** as the next object. The certificate
here is that object in dual form, after one further change: optimize over a universal
`q<=112` obligation atlas instead of the core-dependent `q<=40` universe. Thus the two
results are complementary, not competing. THM-1111 locates the higher-overlap gap; this
weighted dual separates all six-edge unions without enumerating that gap.

The incoming THM-1111 independently shows why stopping at pairwise overlap is not enough:
its maximum-spanning-tree correction is extremely strong on random sextuples but leaves an
adversarial margin of `+36`.  The atlas is the proposed fractional/second-order move made
exact.  It does not approximate those overlaps; it exhibits a positive dual functional
that separates **every** single-killer incidence vector from sixfold coverability by one
integer unit.

The quotient preserves precisely the predicate needed here — existence of a rational
time safe for the core and all six killers. It destroys the geometry of the continuous
safe intervals, the identity of the seven-speed subcore, and killer-killer intersection
data. None of those are required by this finite certificate.

## Tournament Analysis

On the 241 candidate-killer vertices, take the pairwise observable
`load(k)-load(l)` and break equal-load ties by the integer label. This gives a transitive
tournament: zero directed cycles, 241 singleton SCCs, score histogram `0,...,240`, and one
Hamiltonian path; reversing the tie gauge flips 1,414 edges. This tournament faithfully
records the scalar capacity ordering used in the proof, but it destroys obligation owners
and all higher overlaps. The **weighted bipartite incidence hypergraph**, not the runner
tournament, is therefore the underlying proof object. This is a useful negative result
about the tournament lens: pair orientation adds no certificate beyond the dual load.

## What remains open

This theorem eliminates the stated `3.64*10^12` finite enumeration rather than merely
accelerating it. It does **not** yet justify a global r=6 closure. THM-1102 obtained
`max T=308.4` from five removed killers in a width-16 near-bottom window. The fact that the
window maximizer is interior is strong evidence about that scan, but is not a proof that
arbitrary larger or differently spaced quintuples have `T<333`. The remaining obligation
is now sharply isolated:

> Prove a uniform measure-tail lemma for the five-killer residual, or give a scaling/cluster
> normal form reducing every `k_6>=333` case to a verified bounded chart.

Once that tail statement is rigorous, THM-1121 supplies the missing finite half with a
small exact certificate.
