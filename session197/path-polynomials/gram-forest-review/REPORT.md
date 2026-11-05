# Independent review: cyclic Gram kernels and rooted forests

**PASS, no correction.** This review binds the frozen `session197-odd-cycle-path-polynomials/` files:

- CYCLIC_GRAM_FORESTS.md:2fe61a1fcbd4133fbd9bddcdd8a337dfe521d10ca41514f69ecb25f358327c83.
- verify_cyclic_gram_forests.py:f56c9feed4d93b757a35bf7e847fc03a38d623d25214bae59fc9fefa5c34f6f0.
- cyclic_gram_normal.stdout and optimized counterpart:226c04befca6252e61eda03b2534a079c346dd61603ccac47e723e46fd69921c.

## Mathematical audit

A Hamilton word in a cyclic three-block tournament has maximal run colors following the directed cyclic order. With all blocks nonempty, all run counts are positive. Equal counts admit exactly three starting colors; counts differing by one admit exactly one appropriate start. This yields the stated diagonal weights1,3,1 and neighboring weights1,1. Counts of T/U runs differing by more than one are impossible. Since V has at most m runs, the other coordinates stop exactly at L=m+1.

For fixed unordered local covers, every ordering of the distinct fragments into their color occurrences is valid and unique. This gives the i!j!k! factor, and hence J on ordered coordinates but D J D on unordered coordinates. No cyclic rotation is factored out: these are paths, with a distinguished first vertex. No concatenation joins two same-color fragments. The normalization is correct.

The weighted-square expansion is exact, including m1 and the last coordinate. Its first singleton term and positive edge weights force every coordinate to vanish if the quadratic form vanishes. Thus J is genuinely positive definite on all real vectors, not just physical cover vectors. Its positive off-diagonal entries and positive LDL pivots give a positive lower-bidiagonal/diagonal/upper-bidiagonal factorization. These factors are TN and finite Cauchy-Binet proves TN for J and its positive factorial scaling K. This claim is limited to this finite kernel, not arbitrary context-production matrices or Gram arrays in arbitrary input order.

Alternating-sign conjugation turns the off-diagonal entries into the negative weighted path incidence terms. The diagonal residual is exactly lambda1=2f1, lambdar=fr for2<=r<=m, lambdaL=0. For a selected incidence-column set, a squared determinant is nonzero precisely with one coordinate/root column in each tree component; its value is1. The zero last root weight removes those weighted terms rather than adding a fictitious positive root. This proves the forest determinant with no extra orientation or root multiplicity. The full path rooted at1 supplies a positive term. Factorial rescaling gives the stated determinant of K.

Cauchy-Schwarz and its proportional-vector equality condition follow from positive definiteness and nonzero first coordinates. They do not require isomorphism of T and U. The actual inputs Tr1 through TrL have ordered-cover matrix r!S(n,r), lower triangular with diagonal n!. Thus they span the full L-dimensional space and their physical Gram determinant is exactly det K. This is a fixed-third-block rank statement, compatible with unbounded rank as third blocks grow.

Converse preserves all local cover counts. Swapping T,T^op gives an actual converse symmetry only when the remaining block V also has its own reversing map, as the source explicitly requires. Every displayed cyclic graph is strong because the three blocks are nonempty, irrespective of internal strongness.

For V=K1, the unordered quadratic form is exactly2h^2+(h+2u)^2. The nonresidue-2 modulo7 forces both terms' linear factors to be divisible by7 whenever the output is. The only remaining7p possibility would be49, excluded by h>=7 and2h^2>=98. The exceptional prime p7 is handled explicitly. For V=Tr2 the second kernel and its determinants are exact. The mentioned37050461 witness belongs to the separately reviewed THREE_COVER_CYCLIC_GATE.md at its cited frozen hash; it is not re-certified by this determinant audit.

## Independent controls

Two standalone checkers import no source research code.

`check.cpp` counts ordered covers by literal vertex permutations and all legal cuts. It counts Hamilton paths in actual cyclic graphs by inclusion-exclusion over unrestricted length-(n-1) walks, independently of the source subset-path DP. It checks:

- 1331 physical triples formed from every labelled tournament through order3.
- 150 actual entries of the full transitive-input Gram matrices, with graphs up to order11.
- 64 additional physical triples with arbitrary deterministic third blocks of orders4 and5.
- Every physical Cauchy-Schwarz inequality and its exact equality condition in the complete small triple bank.

`check_forests.py` uses exact fraction-free determinants. For24 arbitrary positive coefficient kernels with m1 through6, it explicitly builds the weighted incidence factor, sums4312 squared incidence minors, and separately enumerates individual rooted forests including zero-weight terminal-root cases. It checks all18792 finite TN minors of these matrices, exact rational positive LDL pivots, literal legal color words, both first kernels and all49 norm residue pairs modulo7. Both source and independent determinant calculations use exact integers.

All optimized/debug C++ outputs agree. Normal/Python -O forest outputs agree. The source and its helper were copied to `source/` before normal/-O replay; both outputs reproduce the frozen source receipt. Original files remain unchanged.

## Reproduction and hashes

    g++ -O3 -std=c++17 check.cpp -o check
    ./check
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 check.cpp -o check-debug
    ./check-debug
    python3 check_forests.py
    python3 -O check_forests.py

- check.cpp:0d7f8ce55c3df29d5b1e998139028009354be7f7c1aac87f904f7c759933dda4.
- independent-physical.txt and independent-physical-debug.txt:7454968bc2583609e41efec41fa7130fc7f42fa12695c6b8b9816fbb58e01062.
- check_forests.py:1aadf627a4b742b53ec05b7854acf299129511d017e1acdaf045eab160dbd8b9.
- independent-forests.json and independent-forests-optimized.json:e357a4dad43b1d88d152002fd8ba21897549cc6b88401ffebf35f48235602012.
- source-normal.json and source-optimized.json:226c04befca6252e61eda03b2534a079c346dd61603ccac47e723e46fd69921c.

The all-order theorem is established by the fragment bijection, positive factorization and incidence proof. The finite checks independently verify its physical interpretation and normalization. Exploratory content/gcd claims discussed separately are not part of this frozen review.
