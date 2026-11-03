# Independent audit: explicit projective puzzle coordinates

**Disposition: PASS.** This additive audit does not modify earlier Session 143
checkpoint bundles. No repository changes or publication were performed here.

## Independent construction

The thirteen displayed triples are exactly one representative from each nonzero
scalar class of F3^3. Rather than reuse the producer's cross-product/span code,
this audit constructs all lines as kernels of normalized dual linear forms.
It verifies 13 four-point lines, unique lines through point pairs and unique
intersections of distinct lines.

The legal words are replayed by moving labeled tokens on an actual position
array, then inverting the final board to recover the position permutation.
This independently fixes the chronology convention: later moves compose on the
left. Both words return the hole to zero and give exactly

- `g=[0,1,3,2,4,0] = (1 4 2 12 3)(6 11 10 7 8)`;
- `h=[0,1,7,12,9,4,0] = (1 4 12 2)(6 10 11 8)`.

The orders are five and four, and `h g h^-1=g^2` holds on **all thirteen
positions**, not only on the selected orbit. On `V=(3,1,4,2,12)`, g acts by
index translation and h by multiplication by two. Independent breadth-first
generator closure gives exactly twenty permutations. Their restrictions to V
are twenty distinct affine permutations; hence this is a faithful `5:4` action.

## Tournament and scope checks

The relation with steps {1,2} has exactly one direction per pair and outdegree
two. Its row adjacency is `P+P^2`; independent Newton traces give
`x^5-5x^2-5x-2`. The generated subgroup has a four-object orbit on this relation,
a five-element stabilizer of one tournament, and a five-element kernel on the
family. The multiplier changes the orientation and its square takes converse.
There is no false twenty-element single-tournament automorphism claim.

The direct line image `{0,1,2,12} -> {0,3,4,12}` is not a line. Therefore g is
not a projective collineation. The separate `|PGL3(3)|=5616` order obstruction
agrees. Projective coordinates govern the puzzle moves, not every puzzle
permutation's incidence action.

The underlying puzzle rule and its closed-path identification with M12 are
**CITED**, not rederived from the twenty-element subgroup:
[Conway–Elkies–Martin, section 2.2](https://arxiv.org/html/math/0508630v1).
That primary source was independently inspected during the parent audit;
its move is the required double transposition, and composition retains the
hole endpoint. The explicit labels, words and subgroup witnesses here are
independently finite-verified. No claim maps arbitrary M12/M13 elements to
this fixed five-set; no Baby Monster embedding, Hecke map or canonical
correspondence is established.

## Reproduction

`python3 verify_manifest.py` hashes this bundle, replays the independent checker
under ordinary Python and `-O`, and verifies that all three deliberate
corruptions fail in both modes: omitting the second move transposition,
replacing h by its inverse, and substituting reciprocal edges. Every check
uses an explicit exception; no `assert` node is present. The checker requires
only the standard library and runs from an unrelated working directory.

`reviewed_source_pins.json` freezes the producer proof/checker/output/manifest
reviewed here. The producer's own normal/optimized and hostile tests also
passed, but are separate from the independent implementation.
