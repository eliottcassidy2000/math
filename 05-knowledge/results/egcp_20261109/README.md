# Sharp even-graph Chan-Pak (opus-2026-11-09-S89)

**Result.** [THM-5040](../../../01-canon/theorems/THM-5040-sharp-even-graph-chan-pak.md) proves HYP-9575 (B), and in
sharp form.
- For every n, the dominating-set counts of even graphs on at most n vertices are exactly the odd numbers up to
  M(n) = 2^n - n (n odd) or 2^n - n - 1 (n even), except 3, 7, 13 and, for odd n >= 5, X(n) = 2^n - n - 2.
- Evenness costs at most one vertex compared with all graphs.
- The extra vertex is needed exactly on the top bands (M(m - 1), 2^(m - 1)) and at the values X(m - 1).

Scripts and outputs: `04-computation/experiments/egcp_20261109/` ([README](../../../04-computation/experiments/egcp_20261109/README.md)).

## How the proof goes

1. **Complement formula.** d(G) = 2^n - D(comp G), where D counts the vertex sets lying inside an open neighbourhood.
   - D is additive over components, so joins add deficits.
   - In C_4-free graphs, D = 1 + #(non-isolated vertices) + sum over v of (2^deg - 1 - deg).
2. **Exclusions** (all n):
   - an even graph has at most one universal vertex, so d <= 2^n - n;
   - equality forces n odd;
   - d = 2^n - n - 2 forces the complement to be one of P_3, P_4, C_4 or 2P_3 plus a matching, none co-even.
3. **Four constructions:**
   - E_2-multi-joins, with d = 2^n - 2^g - 2t + d(G); the evenness follows from a multi-join sign rule;
   - cocktail-party graphs, which attain M(n);
   - the leaf stable recursion SR, d -> 3*2^(n-1) + d, which preserves "even, has a leaf, no universal vertex";
   - the jump family J(H, z), with d = 2^(m-1) - D(H), which is even iff the stabiliser of z in Aut(H) is co-even.
4. **Top band at every order.** SR carries each j two orders up. Three explicit sparse H handle the one or two new j
   per order: a path plus K_1 with z an end, the path plus K_1 with z the second path vertex, and a spider. The same
   recursion gives 2^(m-1) + {1, 3, 5, 7, 11, 13} from the census at orders 8 and 9.
5. **Induction** on n for all odd x <= 2^n - n - 3, in three cases:
   - x below 2^(n-1): the top band;
   - D = 2^n - x <= 2^(n-2) + 1: an E_2-multi-join, sidestepping 3, 7 and 13 by moving g up by 2, or else a cone over
     an L-graph;
   - otherwise: a cone over a padded graph, an L-graph or a top-band graph.

## Companion results

- **[THM-5042](../../../01-canon/theorems/THM-5042-odd-graph-companion.md)** is a second, independent proof of mac-mini's THM-5056, which was first in canon. Odd graphs on <= n vertices give
  exactly the odd numbers up to 2^n - 1 except 1, 5 and 37, so oddness costs no vertex. The proof:
  - a cone over a cone is odd (two adjacent universal vertices);
  - padded cones over odd graphs are odd;
  - six census seeds (2^(m-1) + {1, 5, 37} at m = 8, 9, each with an odd automorphism fixing a leaf) carried upward by
    SR.
- **[THM-5041](../../../01-canon/theorems/THM-5041-matrix-tiling-pair-dictionary.md).** The owner's matrix picture
  made exact:
  - an n x n matrix is a diagonal plus an ordered pair (U, L) of (n+1)-tilings;
  - tournaments are lines, semicomplete digraphs are covering pairs, undirected graphs are equal pairs;
  - H ignores loops and is tile-wise additive, so a defect tile carries the same total as the arc-flip pair through it;
  - the 7/21 conjecture is restated in these coordinates (refuted by Session234; corrected form DG2).
- **[HYP-9610](../../hypotheses/HYP-9610-even-graph-domination-no-universal-and-dictionary.md).**
  - The no-universal-vertex sets E0(n). The bound is proved; the sets are exact for n <= 9, and every closure miss at
    n = 10, 11 is realised.
  - The binary/ternary dictionary of exceptional sets, typed ANALOGY.

## What this changes

- **HYP-9575 (B):** OPEN becomes PROVED, together with the exact minimal order V(k).
- **THM-4996 Proposition B5.** The unreached values near 2^m - 2^k are explained by two mechanisms: top bands (which
  need graphs that are neither unions nor joins), and cone/join inputs falling on 3, 7 or 13.
- **Prior art.** Completeness over all graphs is S. Alikhani, Graphs Combin. 29 (2013) 1175-1181. THM-4996 (1)
  reported no prior art; only its conciseness part remains new.

## Open

- The exact no-universal-vertex value sets beyond order 9. They are exact for n <= 9, and cover.py's closure misses
  about 2n values per order, which the theorem never needs.
