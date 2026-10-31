# Bounded-backedge tournament counts: an exact carrier and a width barrier

**Restricted-family theorem and finite positive witnesses, not a proof of the unrestricted Hamiltonian-path spectrum conjecture.**

A width-k presentation orders vertices0,...,n-1 so i→j whenever j-i>k. On appending a vertex, its cut mask is indexed by the retained frontier from oldest to newest; bit1 means new→old. Counts are directed Hamiltonian vertex words, without automorphism division.

## Results

1. An exact endpoint-pairing path-cover transfer has2,6,25,116,567 states for widths0,1,2,3,4. For fixed width it counts any supplied presentation in O(n) arithmetic operations, with growing integer bit costs. This is an implementation claim, not a claim of priority for bounded-width dynamic programming.
2. **No width-two tournament of any order has exactly47 Hamiltonian paths.**
3. Width three does attain47: the seven-vertex cut word is `(0,0,0,1,3,2,1)`.
4. `width3_witnesses.tsv` contains825 explicit witnesses, one for every positive odd H≤1653 except7 and21. Each has at most nine vertices. This finite prefix does not imply arbitrary-value coverage.

## Exact carrier proof

For each processed prefix, choose a directed path cover, including singleton paths. Keep the last min(n,k) vertices as frontier. Record each component's ordered (start,end) pair, retaining frontier labels and replacing forgotten endpoints by−1. Keep the component pairing, not just counts of starts and ends.

A forgotten start cannot receive an edge from a later vertex because their gap exceeds k. Thus two forgotten starts can never be merged into one final path. Discard such states. Every retained component has either a distinct frontier start or the unique forgotten start, giving at most k+1 components and a finite state space.

When x arrives, choose at most one incoming component and at most one outgoing component. Incoming edges are allowed from any forgotten end or a frontier end pointing to x. Outgoing edges are allowed only to a frontier start beaten by x. Choosing the same component for both creates a cycle and is forbidden. Join the chosen components through x, retain all others, forget the oldest frontier vertex when needed, discard two-forgotten-start states, and sum multiplicities of coincident endpoint records. A one-component state is accepted.

Every Hamiltonian path selects each edge exactly when its later endpoint arrives, so it has exactly one transfer trajectory. Its partial edge sets form path covers and never have two permanently unjoinable starts. Conversely an accepted trajectory produces a connected acyclic spanning directed graph with indegrees/outdegrees at most one, hence a Hamiltonian path. Forgotten vertices have uniform future adjacencies, so forgetting their identities loses no future extension information. Keeping endpoint pairing excludes spurious cycle joins.

## All-order width-two obstruction

Tournament strongly connected components are totally ordered; every Hamiltonian path traverses them in that order. Consequently H is the product of their path counts. For prime47, exactly one nontrivial strong component can contribute47, and all other factors are1.

For a strong tournament S of order n≥3, Busch's Theorem1 implies `H(S)^3 >= 5^(n-1)`. Since `47^3=103823 < 5^8=390625`, a strong component with H=47 has at most eight vertices. This uses only the weaker uniform consequence of the cited piecewise bound. Source: Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, Electronic Journal of Combinatorics13 (2006), N3, [Theorem1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf).

Width is hereditary for induced subgraphs: a compressed index gap greater than k implies the original gap was greater than k. Thus the strong component would have a width-two presentation on at most eight vertices. The exhaustive10,923 physical cut words of orders1–8 have no H=47, as checked independently by adjacency DP and the transfer. This is the required finite contradiction, not an assumption that small-order data extrapolate.

## Reproduction

Python3 standard library only. No external data, network, C++ executable, numerical solver, or repository checkout is needed to replay the computations. The all-order theorem additionally invokes the cited mathematical theorem.

```sh
python /path/to/package/run_checks.py
```

This runs `transfer.py`, `verify.py`, and `independent_check.py`, each normally and with `-O`, checking their output against frozen receipt hashes and checking mode equality. Files resolve relative to their source paths, so execution from another working directory is supported.

The independent checker constructs adjacency directly from physical masks, reconstructs every finite carrier by materializing endpoint components as short directed paths, compares281 complete small edge-subset path-cover distributions, repeats all10,923 width-two counts and825 witness counts, and counts47 by literal permutations. Explicit exceptions, not `assert`, enforce checks.

`AUDIT.md` records the independent scope audit; `PROVENANCE.json` identifies the producer inputs. `receipts/` holds exact outputs and relocated-package evidence. `MANIFEST.json` and `SHA256SUMS` freeze the package; publication transport metadata is outside it.

## Boundaries

No new globally forbidden H-value is claimed. The width-three47 witness demonstrates that this is a carrier obstruction. The global claim that7 and21 are the only missing positive odd values remains open here. Unused discovery censuses, C++ programs, spectral diagnostics, and unverified rank claims are deliberately excluded from this proof packet.
