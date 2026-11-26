# Hall distinctness and exact pair-mask synchronization

2026-11-26. Continues [the tight sheet-synchronization test](TIGHT-SHEET-SYNCHRONIZATION.md). These are finite exclusions for two distinct full charts sharing the specified three frozen stars, not a Moore existence result.

## 1. Four spare high incidences are impossible

Let $L$ be the $l$ closed branches of changed degree four, and $Q_a$ the four changed neighbors of each $a\in L$. Outside low branches have only high changed neighbors. Thus the projected sets $Q_a\cap L$ are the rows of a symmetric simple graph on $L$. The deficit $d=\sum_{a\in L}(4-|Q_a\cap L|)$ uses that many high incidences. In the three four-slack rows, $d\le4$ and $p\ge6$, so restrict to the 148 catalogued quartets with at least three antipodal pairs at each target.

If distinct lows $a,b$ choose the same quartet containing any $t\in L$, they each use the entire active set at $t$. At any fixed active sheet there, the replayed oriented trade for that quartet and target sheet is unique. The two old neighbors therefore have another old contact in common, creating a 4-cycle. Thus the low vertices with nonempty internal rows must receive distinct whole quartets.

For $l=4$, $d\le4$ requires a $K_4$. For each of 316,251 label sets, check that every vertex can find its three other labels in an eligible quartet whose fourth label lies outside the set. None succeeds. For $l=5,6,7$, at least $l-4$ vertices have full internal degree, so seed label sets with the union of $l-4$ distinct eligible quartets contained in them. Enumerate all symmetric projection masks for each set, respecting the deficit cap. The replay tests Hall on the lists of compatible whole quartets.

| $l$ | Label sets | Sets with a symmetric graph | Symmetric graphs | Passed Hall |
|--:|--:|--:|--:|--:|
| 5 | 7,362 | 15 | 15 | 0 |
| 6 | 2,169 | 19 | 20 | 0 |
| 7 | 786 | 82 | 91 | 0 |

In increasing order of Hall witness size, **all 126** failures already have a witness of size two: two nonempty rows are forced to the same sole quartet. This eliminates $(h,s)=(38,7),(39,6),(40,5)$.

## 2. One eight-slack row fails exact synchronization

In $(h,s,l)=(40,6,5)$ the cap is $d\le8$ and $p\in\{8,10\}$. The internal graph has at least six edges, hence at least one vertex of degree three or four. Seed from any triple in an eligible quartet, then pad and deduplicate. The replay examines 605,011 label sets and 257,570 symmetric graphs, of which 19,755 on 8,673 label sets pass Hall. No such graph can have an isolated vertex: this would leave a $K_4$, already excluded.

For each surviving graph, choose the actual whole quartets and intersect availability and disjointness masks at every common closed target. Only 48 graphs permit these local masks at either $p$. For every quartet use the rooted catalogue's partial bijections connecting the pair labels at its closed targets. Shrink the masks until each allowed label maps to an allowed label at every other target of that quartet. This necessary fixed-point test leaves 45 graphs (and 45 specific quartet assignments) at each $p$.

The fixed point is not enough. Lift every partial bijection to edges between vertices $(a,x)$, with branch $a$ and allowed antipodal pair label $x$. In any exact choice of $p/2$ pairs per branch, an edge forces its endpoints to be selected together, hence the selected labels are a union of whole connected components. Conversely any such union obeys these specific mask equalities. Use a finite coefficient test: if a component has layer-size vector $w$, multiply by $(1+z^w)$ and truncate each exponent at $p/2$. The simultaneous size is possible exactly if the coefficient at $(p/2,\ldots,p/2)$ is nonzero. **No assignment remains for either $p=8$ or $p=10$.** This excludes the entire $(40,6,5)$ subrow before filling any actual unknown matchings.

The surviving 128-matching rows are $(39,7)$ with $l=6,7$; $(40,6)$ with $l=6$; and $(40,7)$ with $l=6,7$. The lower bound remains 128, and their survival is only a failure of these filters to exclude them.

Replay: `node hall_closed.js ../moore57_quadratic_tensors_20261125/period_fourteen_site.txt --eight`. Omitting `--eight` runs the four-slack certificate only.

## 3. Relation to exact search and integer selection

Paul Stankovski Wagner's [*Less is Moore*, §4.3–4.6](https://arxiv.org/html/2610.04815v1) applies AllDifferent to the rows of a permutation block and uses replayable search frontiers. Here the same matching principle applies a level higher, to *whole* neighbor quartets. Its residual sheet orbit cuts cannot be copied to the marked fixture while holding branch labels fixed: the frozen transport group is $S_{56}$, and the common sheet permutation in the identity gauge would have to centralize it.

The step beyond domain propagation has a direct analogue on the Chan–Pak side: a projected residue or nonempty domain does not select an integer witness. Here every survivor had locally sufficient pair counts, but the exact cardinality vector could not be assembled from the forced transport components. The appropriate inverse object is the component weight vector, not the single-sheet phase.
