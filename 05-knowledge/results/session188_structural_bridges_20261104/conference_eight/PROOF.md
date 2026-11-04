# The eight-name cubic kernel is exactly the skew-conference class

**FINITE-EXACT STRUCTURAL CHARACTERIZATION.** For a tournament Q on eight fixed names, let h be the complete shared-cubic parity matrix from `../eight_pair_quotient/PROOF.md`. Define its skew Seidel matrix by

    S_ii=0,  S_ij=+1 if i->j,  S_ij=-1 if j->i.

Then

    h(Q)=0  if and only if  S S^t=7I.                    (1)

There are exactly30,720 labelled tournaments satisfying these equivalent conditions. The complete count and equivalence use exhaustive exact finite certificates, with the normalization coverage proof below. No equivalence at general orders is claimed.

The independently audited eight-quotient packet already proves that all these quotients have Hamilton start/end parity XOR equal to one at every name, and hence their pair inflations are separated by a single block-uniform endpoint query. This packet explains the finite residual's matrix structure; it is not a new response collision or an unrestricted reconstruction theorem.

## 1. Unique root normalization

A diagonal matrix D with entries in{+1,-1} switches a Seidel matrix to DSD. It preserves skew symmetry, every off-diagonal sign and the conference identity SS^t=7I. Fix D_00=1 and choose D_jj=S_0j for j>0. This is the unique such switching that makes all entries S_0j equal+1, that is, name0 is a source.

The resulting matrix has block form

    [ 0    1^t ]
    [ -1    C  ],

where C is a7-by7 skew Seidel matrix. Multiplication gives the exact equivalent conditions

    C1=0,             CC^t=7I-J.                        (2)

The first says that every vertex of the seven-name tournament has outdegree three. The second has diagonal six and off-diagonal entries minus one.

Conversely, every such C yields a normalized skew-conference matrix of order eight. Applying any of the2^7 diagonal switchings with D_00=1 gives another labelled skew-conference matrix. These128 outputs are distinct. Furthermore different normalized C cannot have overlapping switching outputs, because normalizing any output at name0 recovers C uniquely.

Thus normalization is a bijection between labelled order-eight skew-conference tournaments and pairs consisting of a normalized C satisfying(2) and one of128 switchings. This is the completeness argument, not a heuristic use of sampled switchings.

## 2. Exhaustive normalized count and full-set comparison

The executable checks all2^21 labelled seven-name tournaments in the same fixed edge-code convention as the eight-quotient packet. Exactly2,640 have outdegree three everywhere. Testing every off-diagonal integer inner product in(2) leaves exactly240 normalized blocks. Their complete ordered code list is retained in `census.json`.

The normalization bijection therefore yields exactly

    240*128=30,720

labelled skew-conference tournaments at order eight. The checker constructs all30,720 codes, confirms uniqueness and directly verifies SS^t=7I on every generated matrix. It then compares the entire ordered integer code list with the complete h-zero bank retained by `eight_pair_quotient`, not merely its size or a digest. The two lists are identical, proving both directions of(1).

The Python driver independently reconstructs every switching as a fixed edge-code XOR mask and again compares all30,720 exact codes with the complete bank. This is separate from the C++ row-based generation. The bank's SHA256 binds the dependency, but set equality is checked using full integer values.

## 3. Relation to the four quotient orbits

The eight-quotient packet already records four relabelling orbits of the shared-cubic kernel, with sizes1920,13440,13440,1920. Consequently these are also exactly the four labelled skew-conference relabelling orbits. The middle two are prime, giving26,880 prime examples. The other two have a universal source or sink. Whole converse pairs the two middle orbits and the two outer orbits.

These orbit descriptions are finite exact conclusions, not claims about skew-conference matrices in other dimensions. In particular the equivalence is false if the order-eight hypothesis is simply removed: the transitive four-name tournament has zero shared-cubic signature, but two distinct Seidel rows have inner product two, so SS^t is not3I. The verifier retains this explicit boundary.

## 4. Reproduction and scope

Run `python verify.py > normal.json` and `python -O verify.py > optimized.json`, then compare. The driver uses Python's standard library and an available C++17 compiler, installs nothing, reconstructs the complete normalized list and switching class, and binds the complete dependency bank. All checks remain active under optimization.

The sibling `eight_pair_quotient/blind8.bin` is required. This packet depends on that packet's independently audited proof that the bank contains every h-zero eight-name quotient. The new normalization argument and complete conference enumeration independently characterize all order-eight skew-conference tournaments. There is no general-order extrapolation and no assertion that zero shared cubics imply a full response collision.
