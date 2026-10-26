# Session 69 independent prime-action / literal-response audit

**Verdict: ACCEPTED in the stated scope.** The field-action classification, ordered multi-rooted literal-response convolution, its all-ordered-cover recursive extension, and the Paley28 endpoint/flux calculations pass proof review and independent exact checks. This is not an unrestricted prime census, a bounded/minimal-state theorem, a polynomial-time algorithm, or a converse-merged metagraph degree formula.

Input research pin: `61d7cca5d8df80169f939268269239c2f67ea724`. The repository was not modified by this auditor. Producer source hashes are frozen in `MANIFEST.json`; later publication changes are outside this mathematical audit.

## 1. Accepted mechanisms and convention checks

### Field action and fibre

The polynomial uniqueness proof is valid: after reducing exponents modulo seven on F8*, the difference of two claimed formulas has degree at most six and seven distinct roots. Hence all seven coefficients vanish. Scaling gives `b -> bu`, `c -> cu^3`; Frobenius gives `(b,c) -> (b^4,c^4)`. For `b != 0`, the ratio `r=c/b^3` is invariant under scaling and is classified up to Frobenius. The two `b=0` types and four `b!=0` types per parity give exactly twelve classes. The stabilizer strata are 2/6/4 classes of orders 21/3/1.

Independent exhaustive checks confirm all 7! root permutations, the 21-element full tournament group, the 42-element cyclic-triple-hypergraph group, no proper nontrivial modules, all fourteen cyclic triples as two disjoint Fano systems, H=189, all 128 coefficient triples, converse, all 896 pointed flips, 48 block-/quotient-vertex-pointed classes, and the 8/40 point-stabilizer split. The two distinct stabilizer-three weight-three representatives are `0001011` and `0001101`; the third weight-three class is `0000111`, with stabilizer one.

The complete cyclic-triple-realization fibre uses cited structural input, not a new elementary reconstruction claim. Independently checked primary statements: Theorem 10 gives common strong modules; Theorem 24 gives only a prime tournament and its converse; Proposition 51 gives the independent node-choice bijection; Corollary 57 counts choices. See [Boussaïri–Chergui–Ille–Zaidi, arXiv:1805.04607](https://arxiv.org/html/1805.04607), checked 2026-10-26. This supports the stated `2^15` labelled fibre relative to the identified tree. The existing twelve-class census is recovery, not a novelty claim.

### P2: ordinary source-presentation incidence

The proof is a genuine bijection: split an actual Hamilton word at every block change; order each child's runs by occurrence; retain the directed quotient walk. Conversely, these data concatenate uniquely. Reversing all source arcs and restoring internal run edges gives the local targets; restoring cross-block path edges gives exactly the hinges.

Multi-rooted ordinary target isomorphisms preserve every ordered endpoint pair, including repeated endpoints for singleton runs, so independently chosen representatives glue to isomorphic targets. Aggregated local response multiplicities count Aut(child)-orbits, not bare target states. The local action is free on the actual ordered spanning covers, because every position of every component is fixed. Consequently raw counts may be divided by the child automorphism order separately in each response category.

The top child-type stabilizer `K` acts freely on the retained labelled quotient walk, since every slot occurs. Equivalently the underlying Hamilton-word action is free. This does **not** license dividing the number of response marginals after the walk is forgotten. Ordinary target categories are invariant, so the aggregate division by `|K|` is valid. Ordinary source presentations, labelled words, target classes, and merged-converse edges remain different objects.

### P3: recursive closure for every ordered-cover response

The same proof works for all ordered d-path covers, not just Hamilton words. Use ordered d-tuples of nonempty quotient walks, assigning each child's run components in concatenated tuple order. Toggle hinges only inside each walk, never across different global paths. Retain the first and last endpoints of each global path, in order. Divide weighted response totals by `|K|`; there is **no `d!` factor**, because global components are ordered.

Splitting and gluing are inverse, preserve the entire output relation and all ordered marks, and commute with the source automorphism action. This recursively computes the full finite response family from child responses and the actual quotient relation/action. Target modules may shatter; no persistence of the source partition is assumed. Linear canonical nodes have the analogous formula with trivial order-preserving quotient stabilizer. A singleton has its unique one-path response. Thus this is exact finite recursive composition given the canonical input tree, not a catalogue or efficient classification of arbitrary prime nodes.

## 2. Independent numerical method

The producer uses bounded-content last-slot dynamic programming. This audit independently computes endpoint walk series by Cramer's rule:

`W_ij(x) = x_i Cofactor_(j,i)(I-A_Q X) / det(I-A_Q X)`.

Determinants and cofactors are expanded directly over permutations; no producer routine is imported. The reciprocal denominator is computed by integer coefficient recurrence through exponent four in each of seven variables. The Paley denominator has 65 nonzero squarefree monomials. Local endpoint tensors are independently generated by recursive ordered partitions of vertex subsets into actual directed paths, not the producer's permutation-and-cut loop.

Every diagonal/cross endpoint channel and all eight raw flux histograms agree exactly. In particular,

`H(Paley7[A_w]) = 5307704039306577787515`.

This calculation uses `5^7` coefficient cells, not a 28-vertex subset DP. Full transfer-to-graph endpoint histograms are also checked against actual small Hamilton words.

Regularity makes each cut balanced: the directed mixed-colour arc counts are `k(7-k)/2`. With the two ordered-pair orbitals and equal diagonal/same-colour response channels, the full histogram depends only on `k(7-k)`. The standard literal-partner identity is independently verified on every small control path. For distinct global endpoints it is `delta c3 = score(first)-score(last)-1`; the isolated one-vertex case is not in the prime-root scope.

Hence the two stabilizer-three weight-three Fano classes have equal H, full automorphism order, and complete raw flux histogram while remaining nonisomorphic and not converses of one another. The endpoint/flux observer is decisively incomplete.

## 3. Hostiles that locate the lost information

All bitcodes use upper-triangle pairs `(0,1),(0,2),...` in lexicographic order, low bit first, with bit 1 meaning `i -> j`. Canonical codes are minima over all ordinary relabellings, never converse merging.

### Scalar profiles do not determine target incidence

`C3[A0,1,1]` and `C3[A1,1,1]` both have H=27, Aut order3, nine ordinary ports, and identical complete child path-cover profiles. Their nine ordinary target classes differ:

- A0: `{8,19,22,83,145,149,154,179,281}`;
- A1: `{4,17,18,73,83,145,146,152,154}`.

### Even target plus two outer endpoints is insufficient under revisits

Let the four-vertex child S have bitcode7 (it is transitive), and take

- `C=((0,1),(3,2))`;
- `D=((0,2,1),(3,))`.

Both are ordered two-path covers. Their local complement targets have the **same ordinary class and the same global-outer-endpoint rooted class**, represented by `((0,1),4)`. But the quotient walk `(0,1,2,0)` in `C3[S,1,1]` produces targets of canonical codes **81 and 83**. Their ordered internal component ports differ; retaining only the first endpoint of the first component and last endpoint of the last component loses exactly the hinge information.

### Response multiplicities cannot be discarded

For the five-vertex source S of bitcode10, Aut(S)=1, and some distinct ordered covers give the same multi-rooted response. Examples of response multiplicity2 are:

- c=1, target code278, marks `((0,1),)`;
- c=2, target code447, marks `((0,1),(2,3))`;
- c=2, target code393, the same ordered marks.

`C3[S,1,1]` has H=93 and 93 ordinary ports. The correct weighted P2 total is93; replacing every positive local response weight by1 gives **88**. The missing five are three occurrences of the multiplicity-two c=1 response and one occurrence of each of the two multiplicity-two c=2 responses. This is an actual undercount, not merely a formal warning.

### Bare response-marginal orbit division is invalid

For `C3[1,1,1]`, there is one equal child-response tuple but `|K|=3`. Dividing that marginal count gives `1/3`, whereas the three quotient walks give one ordinary Hamilton presentation after the lawful division. The walk is part of the free carrier.

## 4. Finite universe and reproducibility

`audit.py` imports no producer code. It checks 31 target-incidence substitutions: every C3 child tuple from `{1,T2,C3,A0,A1}` of total order at most6, five prime-five substitutions with one edge child, all-singleton prime-five and Paley-seven controls, and the seven-vertex multiplicity hostile. It directly verifies 1,003 Hamilton words and their literal flux identities.

The all-d P3 check independently compares complete multi-rooted response histograms on seven sources:

`C3[1,1,1]`, `C3[T2,1,1]`, `C3[C3,1,1]`, `C3[T2,T2,1]`, `C3[T2,T2,T2]`, `C3 -> 1`, and `T2 -> T2`.

This covers **7,523 ordered spanning covers and 3,044 response categories**, with all d through the source order. The all-edge-child C3 case tests nontrivial top stabilizer3 with non-singleton children; its ordinary ordered-cover counts by d are `15,114,372,648,600,240`.

Reproduce from any directory after preserving the packet's internal relative paths:

```bash
python3 -I -B path/to/session69-tournament-audit/audit.py > normal.out
python3 -I -B -O path/to/session69-tournament-audit/audit.py > optimized.out
cmp normal.out optimized.out
```

Both checked transcripts are byte-identical; checks are explicit exceptions and remain active under `-O`. `audit_results.json` records the exact 4,464 checks and complete target-incidence tables. The sibling input JSON is only a frozen producer comparison reference used after the independent endpoint arithmetic; its provenance is recorded in the manifest.

## 5. Remaining limits

No all-order untouchable-prime uniqueness, ABC result, interpretation of the raw-H gaps as group orders, unrestricted prime enumeration, minimal closed compression, or converse-merged/reflection-coloured incidence is inferred. The original 1,2,4 orbit is literal field arithmetic, whereas the function module decomposition is 1+3+3. Publication and repository checks belong to the owning publisher, not this read-only research audit.
