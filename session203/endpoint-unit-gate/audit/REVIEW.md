# Independent audit of the fixed rooted unit-correction obstruction

6 November 2026. **PASS in the stated scope**, including the separate tail addendum and the parity-compatible positive identities. This is an obstruction to a physical interface, not a missing Hamilton-path count or a spectrum threshold.

## Frozen scope and inputs

The clarified candidate PROOF.md has SHA-256 00acb8e5db1ab86f8547d998f20b936b84b1905180f2be19683551d14160d1bd. Its sole change from the initial e50639bdbbc35201f6c3f08f72febcb7c44e47f3e54e2b7a10e2ce92bae2d5fe candidate says explicitly that contraction occurs in W before adding the apex. Both versions are retained here; no mathematical change was required.

The reviewed tail addendum has hash ffa6f030bbae46fd85a9b5e29cfe15f9f4ef1c81cccad89db58f42be3691aba4. The family-C companion has hash 43dc39855ddadd2971d06caa92c03f6a0e3b380269ccb20f3db616e7c7bc077d. The producer verifier has hash 20bff3e4b4deea8cd7fa58fb0df84d71269db24f72f2046fa425649f87c2d242.

Each context contains exactly one copy of the distinguished input root, identified with one fixed core vertex. All other fixed core vertices treat every nonroot input vertex uniformly. The two cores need not be related, need not be small, and may differ in order, masks and edges. The input's nonroot graph is retained literally. The optional global reversal of copied-root incidences is also covered.

This definition must not be confused with deleting the input root and retaining only U-v. In that zero-retained-root case, all four test inputs become the identical transitive chain, so their output mixed differences vanish. The present two-one-root theorem does not establish the corresponding two-zero-root exclusion. Conversely an ordinary uniform-module construction retaining the entire U is a special case of the one-retained-root class: absorb its retained input root into the fixed core. Merely not distinguishing that retained root does not create an escape.

## Twin curvature and the all-order contradiction

Only Hamilton words with q internally adjacent to both x and y can contain both free bit variables. Their contributions are E_xy(1-b_x)b_y and E_yx(1-b_y)b_x. Thus their mixed difference is -(E_xy+E_yx).

The pair x,y is twin outside its mutual edge. Contract its exposed consecutive block in W, not in the apex-augmented graph. Each of the two possible exposed block orders is bijective with a Hamilton word of the same actual quotient tournament Q. The exposed internal step is exempt from the directed-path condition, so its original orientation does not matter. Hence E_xy=E_yx=H(Q) and the mixed difference is exactly -2H(Q).

There is no realizability relaxation: the common outside incidences define a genuine quotient tournament. Complementing both copied bits preserves the mixed difference. The proof does not require a bound on the quotient order or its path count.

The four five-vertex sources have private chain x1->x2->x3->x4 and root bits 1,a,b,0. Their spanning cycle proves strongness. Deleting the root leaves a transitive tournament, so any root-start Hamilton path must use its unique order; the first bit makes this path valid. Their counts are exactly 9,9,11,9, and all four have start count 1. These are physical source graphs, not four freely chosen scalar profiles.

After deleting the copied root in either output context, the adjacent pair x2,x3 is twin. By classical Rédei Hamilton-path parity, both actual quotient counts are odd. For every even integer k, including zero and negative even integers,

    Delta(H(child)-kH(parent))
      = -2(H(Q_child)-kH(Q_parent)) = 2 modulo 4.

A constant or affine function of a,b has mixed difference zero. Therefore the proposed fixed-context identity fails on at least one of these four inputs. This uses the classical all-order parity theorem, not a parity extrapolation from the finite controls. No output strongness is assumed or required.

## Tail addendum

The longer root-incidence words 1,0^t,a,b,0 preserve both strongness and the unique root-start Hamilton path. The varied vertices remain adjacent twins after deleting the copied root. The common induced prefix 1,0^t is a cyclic triangle with one vertex replaced by a transitive t-chain. Its three one-run possibilities and 2^t-2 two-run assignments give H=2^t+1 for t>=1. Induced insertion monotonicity therefore makes all four full source counts exceed any chosen cutoff.

The same contradiction rules out a fixed even-multiplier constant gate on all sufficiently large strong unit-start sources. The addendum also correctly allows cores and masks chosen separately for each length, provided each chosen pair is fixed across that length's four assignments. No bounded core size is needed. It does not allow choosing a different context after seeing the free bits or target.

## Positive identities and the larger family-C companion

The small edge and triangle contexts give H(parent)=h+L-s and H(child)=h+2L, where L is the actual ordered two-path-cover count with the designated nonempty path ending at the root. The inherited ordered-root-amalgamation formula gives these identities directly from the core profiles (1,0,1,2,1) and (3,1,1,3,1); its singleton boundary is valid. Their difference is 2s-h, not a constant on the strong unit-start class.

The larger family-C tensor was independently re-enumerated by selecting mandatory and optional occupied gaps, rather than multiplying the producer's gap polynomials. Every one of the 7!+8!=45,360 core orders was examined. The only surviving coefficients of child minus twice parent are

    A_(0,0,0) + A_(1,0,1) - A_(1,1,0) - A_(2,1,1).

Here the first root-adjacency flag constrains an incoming run, and the second an outgoing run; when both are present they concern distinct runs. Fixing the core order fixes their ordered positions, so there is no extra bin-permutation factor. For nonempty V=U-v, the root-start, root-end and internal-root path decompositions give h=A_(1,0,1)+A_(1,1,0)+A_(2,1,1) and s=A_(1,0,1). For empty V the surviving A_(0,0,0)=1 equals 2s-h because h=s=1. Thus the larger physical replacement also has correction 2s-h for every rooted input.

The original core vertex 0 remains a universal sink after the stated birth, including toward all input nonroot vertices. Its deletion preserves H and verifies the same-order interpretation. Singleton U gives 25->51; a directed triangle gives 195->389 and correction -1. Neither identity is a unit extractor.

## Independent verification

The independent scripts import no producer code and use explicit exceptions rather than removable assertions.

check.py uses nonmemoized literal directed-path enumeration and permutation/cut counts, not the producer's subset DP. It verifies:

- 549 labelled contexts with core order 1 through 4 and root labelled 0.
- 1,098 curvature identities including both global copied-root polarities.
- 74 additional curvature identities with the twin pair's mutual edge reversed.
- 1,099 labelled rooted inputs through order 5 for both small positive identities; fixing root label 0 covers every rooted isomorphism class by relabeling.
- 8,058 output/quotient graph calculations plus the four physical sources, with 6,208,430 literal path-prefix visits.
- The zero-retained-root scope control: deleting v erases both free bits in each input square.

check_tail_companion.py verifies the prefix formula and all four strong/unit-start sources for t=1,...,7, independently enumerates the family-C tensors, and checks 75 labelled rooted inputs through order 4 in the larger parent, intermediate child and sink-deleted child. These controls visit 4,323,211 path prefixes. Both scripts' normal and optimized receipts are byte-identical.

The observed quotient H values 1,3,9,17,23 are finite controls, not the scope of the theorem. Its arbitrary-order conclusion follows from literal twin contraction and all-order parity.

## Reproduction and limits

    python3 check.py > NORMAL.json
    python3 -O check.py > OPTIMIZED.json
    python3 check_tail_companion.py > TAIL_COMPANION_NORMAL.json
    python3 -O check_tail_companion.py > TAIL_COMPANION_OPTIMIZED.json

The theorem does not exclude target-dependent context selection, multiple sensitive copied vertices, position-dependent incidences, modifications to the input's internal arcs, nonlinear encodings, or smaller input families avoiding the four-input squares. It does not show that a zero-retained-root method succeeds. No such exit is asserted to solve occupancy or to prove an eventual odd-count threshold.
