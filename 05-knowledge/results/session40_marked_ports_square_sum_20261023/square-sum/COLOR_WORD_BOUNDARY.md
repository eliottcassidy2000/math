# Birth-color counts, ordered words, and an explicit order80 hostile

**FINITE-EXACT** companion to [README.md](README.md). The first packet is unchanged. This note tests a natural extension of the user's finite three-color partition, rather than assigning the user an unstated all-order rule.

## 1. Extend the birth interpretation, not an arbitrary coloring

Color an edge by its birth vertex `b=max(a,b)`:

1. first color if `Q_b` has a Hamiltonian path;
2. second color if `Q_b` has no Hamiltonian path and `b` arrives with two neighbors;
3. third color otherwise.

Through25 this equals the user's specified first and second classes exactly. The current note provides direct Hamiltonian path witnesses for **each** birth order25 through80. Thus all edges whose larger endpoint lies in25..80 have first color, without needing an infinite theorem as a dependency.

## 2. Marginal colors already miss a port ambiguity at26

Orient each unoriented path from its smaller endpoint and retain the actual endpoint pair. At26 the two distinct paths with endpoint pair `{13,18}` displayed in the receipt have the same color histogram `(8,11,6)`.

Hence endpoint markings plus all three color totals do not reconstruct the path. The single actual edge-membership bit `3–6` does separate this pair and the other two endpoint collisions at26, as established by the main packet. Color names and totals are not a substitute for shared incidence.

The counts of endpoint-plus-histogram labels at orders25..32 are

`10,11,25,38,19,20,292,339`,

against actual path counts

`10,12,35,52,19,20,349,392`.

## 3. Ordered color words work in the finite census, but fail at80

Retaining the **ordered** color sequence along the oriented path, instead of just its histogram, separates every Hamiltonian path with its actual endpoint pair at every tested order25..32. This is a finite positive result, not an all-size theorem.

The file [COLOR_WORD_HOSTILE_80.json](COLOR_WORD_HOSTILE_80.json) supplies two different Hamiltonian paths of `Q_80` with the same actual endpoints1 and2. Every edge in both paths has larger endpoint at least25. Therefore both paths have the **identical all-first-color word of length79**. The existence witnesses for every relevant birth order are in [birth_witnesses_25_80.json](birth_witnesses_25_80.json).

Consequently even endpoint marks plus the complete ordered three-color word are not an all-size decoder for this natural color extension. No claim is made that80 is the first counterexample. The construction search was heuristic; the final counterexample is an exact certificate checked by integer square roots, vertex coverage, endpoint equality, and color equality. Its validity does not depend on the search heuristic.

The failure mechanism is concrete: at later sizes the positive-birth color occupies a large subgraph. Two paths can avoid all old exceptional-color edges entirely. Recording the order of a single repeated color cannot recover their different vertex incidences.

## 4. Strongest survivor

A full **square-sum word** does reconstruct a path from its first vertex: if the next edge has sum `s`, then the next vertex is exactly `s-x`. This is the inherited reflection transport, with domain and no-repeated-vertex checks retained. The three birth colors discard these numerical shell labels. The endpoint-plus-cotree incidence decoder in the main packet is another rigorous all-size survivor.

There is therefore a hierarchy, with explicit tested losses:

`color totals < ordered color word < actual edge incidence / exact shell transport`.

This notation describes information retention, not a claim that each level is sufficient. The first fails already at26; the second has the explicit80 hostile; the last two have the stated reconstruction proofs and their validity gates.

## 5. Reproduce

Run `python3 check_color_words.py > color-word-receipt.json` and repeat with `python3 -O`. The receipts are identical. The verifier reads the main census and the two explicit witness files. It verifies56 independent birth-existence witnesses, every finite color-word label through32, the26 histogram collision, and both80 paths. It neither enumerates every path at80 nor makes any all-large-order conclusion.
