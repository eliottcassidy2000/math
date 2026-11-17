# Neutral choices of the core: count-preserving arc flips, fibre by fibre, and the fibres that have none

Owner: *"Work to change the core so it gives more neutral choices without creating another free path run."* A **neutral
choice** here is a balanced arc `a → b` of a tournament, `M(a, b) = M(b, a)`, whose reversal preserves the count
(THM-5267 (A)); it is *trivial* when the reversed tournament is isomorphic to the original and *nontrivial* when it moves
to another class of the same count — a change of the core at fixed `h`. A **free path run** is transitive padding
`M_a ⇒ S ⇒ M_b`: it doubles every response and adds only trivial choices (Lemma below). The *fibre* `F_h` is the set of
classes with count `h`; its *neutral graph* joins classes related by a nontrivial balanced flip.

## 1. Padding is transparent (PROVED)

**Lemma.** Let `T = M_a ⇒ S ⇒ M_b` with `S` strong of order ≥ 3. The balanced arcs of `T` are the adjacent-rank arcs of
the paddings (trivial: the flip is again transitive padding) and the balanced arcs of `S` (trivial or nontrivial exactly
as in `S`). *Proof.* A non-adjacent flip inside a padding creates a strong block of count ≥ 3 and multiplies the count.
Flipping `x → y` with `x` a padding vertex and `y ∈ S`: if `x` is not the last vertex of `M_a`, every old path survives
and the path (`M_a` before `x`, a Hamiltonian path of `S` ending at `y`, `x`, the rest of `M_a`, `M_b`) is new; if `x` is
the last, the `s_y(S)` old paths through `x → y` are lost and the pairs (path of `S` ending at `y`, rest of `S`) give
`E_y(S) ≥ h(S) > s_y(S)` new ones (`y` is not a source of the strong `S`); symmetrically for `M_b`. Flipping `x → z` with
`x ∈ M_a`, `z ∈ M_b`: no old path uses it and (`M_a − x`, `S`, `M_b` up to `z`, `x`, the rest of `M_b`) is new. Arcs inside
`S` flip to `M_a ⇒ S′ ⇒ M_b`, isomorphic to `T` iff `S′ ≅ S` (the condensation is canonical). ∎

So the neutral structure of a fibre is that of its strong realisers, and since a count `h < 125` has strong realisers
only at orders ≤ 9 (Busch's floor `f(10) = 125`), the census below is **complete for every `h < 125`**.

## 2. Census ([receipts/neutral.out](receipts/neutral.out) orders ≤ 8 from the pair matrices of `profiles8.npz`; [receipts/neutral9.out](receipts/neutral9.out) order 9, `neutral9.c` + labelg)

| order | classes | balanced arcs | trivial | nontrivial flips | classes with a nontrivial move |
|---|---|---|---|---|---|
| 3 | 2 | 2 | 2 | 0 | 0 |
| 4 | 4 | 6 | 6 | 0 | 0 |
| 5 | 12 | 20 | 16 | 4 | 2 |
| 6 | 56 | 106 | 62 | 44 | 20 |
| 7 | 456 | 652 | 266 | 386 | 204 |
| 8 | 6880 | 6004 | 2060 | 3944 | 2066 |
| 9 | 191536 | 82384 | 54928 | 27456 | — |

**The fibres with no nontrivial neutral move, for every order** (`h < 125`, complete by §1):

    1, 3, 5, 11, 13, 17, 19, 31, 35, 39, 47, 53, 55, 57, 63, 67, 85, 95, 113, 115, 121.

Only 75 among the fibres rigid at orders ≤ 8 gains moves at order 9 (four edges); 85, 95, 113, 115, 121 do not (N-7's
"half of them" refuted). The list contains the near-join exceptions 11, 13, 47, 67 (THM-5263), the +2-window failures 11,
17, 19, 35 of THM-5267 (27 is not rigid: eight edges), the atom 63, the fertile-only 3, 5, 17 (THM-5264), and THM-5260's
tower-blocked 53, 67, 113. Their response sets are small for their size: `|𝔅(h)|` at orders ≤ 8 is 26, 18, 28, 68 for
11, 13, 47, 67 against 47–94 for the neighbouring flexible fibres; 63 has 16.

Selected fibres (classes of orders ≤ 8, strong among them, balanced arcs, nontrivial edges, neutral-graph components,
`|𝔅(h)|`): `9: 30, 2, 79, 10, 10, 35` · `15: 40, 4, 68, 6, 6, 47` · `29: 30, 6, 68, 18, 6, 68` · `37: 36, 7, 148, 36, 6, 76`
· `45: 36, 12, 53, 5, 1, 94` · `59: 28, 12, 58, 12, 6, 107` · `89: 49, 21, 71, 24, 12, 172` · `109: 50, 22, 84, 24, 18, 169`.
The fibre with the most nontrivial edges at orders ≤ 8 is 37 (36 edges); at order 9 it is 547 (58 edges), then 657,
523, 513, 509, 425, 415, 345.

**Neutral-richest cores.** Order 6: the class of count 37 with scores `1, 2, 2, 3, 3, 4` (7 balanced arcs); order 7: counts
51 and 33 (6 each; the Paley tournament has none); order 8: count 657, scores `3^4 4^4`, next to the maximum 661 (7);
order 9: count 495, scores `1, 1, 4, 4, 4, 4, 4, 7, 7` (10 balanced arcs), then 195 and 129 (8). Neutral richness is not
maximal count and not regularity.

**Neutral choices enrich the responses.** Over the 61 fibres with odd `h ≤ 125`, Spearman's correlation between the
number of nontrivial neutral edges and `|𝔅(h)|` is 0.613 (N-6 confirmed): the cores with more count-preserving
re-arrangements offer more responses `B` at the same `h`, which is what the near-join target equation `p = h + 2^k B`
consumes.
