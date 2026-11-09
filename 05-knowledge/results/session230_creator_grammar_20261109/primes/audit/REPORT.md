# Independent audit: prime kernels and Fibonacci fans

**Verdict: PASS.** No mathematical defect found in the requested claims. The finite census is now certified independently through order **eight**, stronger than simply inheriting the generator's completeness.

## Deliverables and independence

`check.py`, `check.out`, `check_factors.py`, and `check_factors.out` are in this directory. In the published `primes/audit/` layout, run `python3 check.py > check.out` and `python3 check_factors.py > check_factors.out`. In this scratch layout use `python3 check.py ../session230-prime-cores > check.out`. The checker imports no producer code. Its only producer inputs are the candidate codes and declared Hamiltonian counts in the generated `../classes_1.txt` through `../classes_8.txt`.

The main checker uses Boolean adjacency matrices, recursion counting paths by their **last** endpoint, explicit combinations for every proper nontrivial subset, and independently constructed degree-partition canonical keys. This differs from the producer's outgoing-mask extension DP and bitmask module test.

### Independent completeness certificate

For each candidate representative, vertices are grouped by outdegree, groups ordered by degree, and every permutation within those groups is enumerated. The minimum upper-triangle Boolean word is a canonical key. Isomorphisms must preserve outdegree, so this enumerates all possible canonical labelings and cannot miss an isomorphism. The number of permutations giving any one resulting word is exactly the automorphism-group order. The checker establishes distinct keys for all supplied representatives and sums `n!/|Aut(T)|`.

The resulting exact orbit masses are

| n | classes | orbit mass | required labeled tournaments |
|---|---:|---:|---:|
| 1 | 1 | 1 | 1 |
| 2 | 1 | 2 | 2 |
| 3 | 2 | 8 | 8 |
| 4 | 4 | 64 | 64 |
| 5 | 12 | 1024 | 1024 |
| 6 | 56 | 32768 | 32768 |
| 7 | 456 | 2097152 | 2097152 |
| 8 | 6880 | 268435456 | 268435456 |

Pairwise nonisomorphic valid representatives whose orbit sizes sum to `2^(n choose 2)` exhaust all tournaments. Thus **completeness is independently certified**, though production of the candidate representative list remains inherited. No trust in the producer's class count, canonical-labeling implementation, or automorphism fields is needed. All declared H values are independently reproduced.

## Finite conclusions

The prime / noncreator-prime / critical-noncreator-prime counts for n=3,...,8 are respectively:

- n=3: 1 / 0 / 0.
- n=4: 0 / 0 / 0.
- n=5: 3 / 0 / 0.
- n=6: 15 / 0 / 0.
- n=7: 197 / 8 / 1.
- n=8: 4008 / 214 / 0.

Noncreator means its H occurs at a smaller order, checked against the union of all earlier independently verified spectra. The eight order-seven noncreator primes have H multiset `[27,31,31,41,41,41,43,43]`. This proves the first substitution-prime obstruction to the all-first-creators seed collection occurs at order seven.

The six-vertex code 6241 has row masks `[2,12,33,21,7,27]` and H=31. The all-order proof below makes W7 prime and critical and gives H(W7)=31, certifying the claimed smallest obstruction without relying solely on its position in the census.

Exactly two order-eight representatives have H=63. Their complete proper nontrivial module lists are:

- 135528708: `(0,1), (2,4), (5,7)`.
- 36175876: `(1,2), (4,6), (5,7)`.

Their deletion-H multisets are both `[1,11,25,25,33,33,37,37]`. Contracting the three pairs yields a quotient canonically identical to W5. The separate explicit-factor checker verifies every cross-block arc in the producer's claimed W5 order, giving sizes `(1,2,2,2,1)` and `(2,2,2,1,1)` respectively. These are transitive-child inflations, not new prime seeds.

## All-order proof audit

### W-family

The chain-plus-alternating-hub definition is consistent. A module avoiding the hub must be a chain interval, and a nonsingleton such interval contains both parities, hence is distinguished by the hub. A module containing the hub and a chain vertex must contain both endpoints: omission of the first or last endpoint distinguishes the hub from the existing chain member. With both endpoints present, any omitted even or odd chain vertex distinguishes the hub from one endpoint. Thus the module must be the whole tournament.

The deletion argument is correct: deleting an interior chain vertex makes its neighbors a homogeneous pair; deleting the first or last vertex leaves a source or sink; deleting the hub leaves a transitive tournament. For odd order at least five this proves criticality. Order three has no prime deletion under the report's order-at-least-three convention but is appropriately described as the seed exception.

For Hamiltonian paths, both sides of the hub must increase along the chain. Endpoint placements contribute two. An interior hub must have odd predecessor a and even successor b. The a<b case requires b=a+1 and contributes k-1; for a-b=2d+1>0, the intervening 2d vertices choose either side freely, giving `(k-d)4^d`. Hence

`H(W_(2k+1)) = k+1 + sum_(d=0)^(k-1) (k-d)4^d = (4^(k+1)+6k+5)/9`.

Consequently `H(W_(n+2))=4H(W_n)-n`. The standalone checker verifies primality, all deletion cards, and these counts for orders 3,5,...,13: 3,9,31,117,459,1825. The proof, not these tests, supplies the all-order assertion.

### Substitution obstruction

Every nonempty child of a substitution is a module. A prime output of a genuinely nontrivial substitution therefore has singleton children, so its quotient is the same prime. Removing identity wrappers from an expression shows every prime output must already be a seed. The W-family supplies unbounded prime orders. Conversely, recursively contracting a proper module, with singleton and directed-pair bases, proves that all primes suffice. These statements assume the usual nonempty-child substitution convention. They do not assert the failure of the richer attachment grammar.

### Fibonacci fan kernels

I inspected THM-4748 at `01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md`. It defines F_k exactly as claimed and already proves its Mersenne count. The source's construction and count are inherited prior art, not new claims here.

For even k=2r, the hub losses are r,r+1,2r+1. Maximal constant-orientation chain runs have lengths `(r-1,2,r-1,1)` and alternating hub orientation, so their quotient plus hub is W5. For odd k=2r+1, losses are r,r+2,2r+2; the run lengths are `(r-1,1,1,1,r-1,1)`, yielding W7. Each run is a module by the chain-interval property. This proves the claimed isomorphisms for every r>=2:

`F_(2r)=W5[TT_(r-1),TT2,TT_(r-1),1,1]`;

`F_(2r+1)=W7[TT_(r-1),1,1,1,TT_(r-1),1,1]`.

The checker directly constructs both adjacency matrices and checks exact equality after the stated chain-first/hub-last relabeling for k=4,...,14, and independently confirms H=2^k-1 throughout. In particular F5=W7, and F6 is the W5 inflation with three directed pairs. No assertion that W5 and W7 generate all tournaments follows or is made.

## Boundaries

This audit does not independently reprove the cited Schmerl–Trotter classification, H-neutral insertion theorem, or the richer five-operation grammar. Those ancillary claims remain citations, separate from the elementary substitution obstruction and the independently certified finite census. No all-order first-creator status for W_n (n>7) has been established here. The kernel reinterpretation is a consequence of the existing explicit fan construction, not evidence that the construction or Mersenne realization is new.

## Publication hardening

Both checkers use explicit `require` / `ValueError` checks, with no Python `assert` statements. All verification runs under `python -O`. Both have main guards, so importing either performs no census; the factor checker does not import the census checker at all.

`check.py` resolves input from an optional positional directory argument, then `PRIME_AUDIT_SRC`, then defaults to `Path(__file__).resolve().parents[1]`. That default is the correct `primes/` input location when this script is published as `primes/audit/check.py`.

Fresh ordinary and optimized runs byte-match both each other and the original frozen outputs. They are retained as `check.normal.out`, `check.optimized.out`, `check_factors.normal.out`, and `check_factors.optimized.out`. An optimized negative control explicitly called `require(False, ...)` in each imported module and received the expected `ValueError`. AST inspection confirms no `assert` statements remain. `SHA256SUMS` records all finalized checker and output hashes.
