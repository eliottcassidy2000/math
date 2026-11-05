# Independent review: retained-root sewing

**PASS. No correction.** Bound to RETAINED_ROOT_SEWING.md SHA256 c5a7a18611363985b49d4ab455ad42160150ab012157c7d4c01aa6f838332655; verify_retained_root_sewing.py 62d824eb95905c6c408b38079eb7f56a3c7e8226539f9bc7588f3d1c3855beed; helper verify.py 10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb; and both source receipts 2b1b3ae04c03b8ff8a3985e0a6bb4227acbb1ec0e544083bb60013fa819111eb. Originals were preserved. Source replays used the copy in `source/`.

## Mathematical review

The operation retains the old vertex a. Deleting it removes the only reverse cross edge, so d'=Hd. A Hamilton path ending at a must have b->a as its final edge, yielding e'=Ed even when A is a singleton. The scalar row uses the inherited forbidden-boundary/one-return split correctly.

For R', a is first in its reserved path and absent from the other, so neither bin can use b->a. Both bins have A before B. The only invalid joining boundary is a->b, and it requires the A part of the reserved bin to be exactly the singleton a. The unreserved A part contributes d and the B allocation contributes Q. Thus R'=CR-Qd, without an omitted second reverse-edge case. Converting to x=R-d and y=d-e gives exactly the displayed positive four-state law.

C>=2Q follows from the bin-swap bijection; Q>=H and L>=H follow from cutting at the marked vertex. For a non-singleton B, all |B|+1 cuts of each Hamilton word are distinct two-bin covers, hence C>2H and C-Q-H>0. This proves strict activation, not merely a finite observation. The singleton coefficient is zero. The stated four physical inputs and continuations give determinant16, so the rational linear lower bound has the correct realization scope.

The pruning condition uses h_next>=h+x. Terminal outputs are recorded before their states cease to be expanded. Exact four-state merging is valid for every later bank word, and the physical inequalities bound all expanded coordinates. Consequently the finite closure really decides the full arbitrary-depth bank language. Initial singleton letters have exactly the seed's whole four-state, so removing them preserves every later continuation. The proof does not strip internal padding. After the first non-singleton, x>=1 and doubles at every further step; this gives h_m>=2^(m-1) and the claimed count-equivalent witness depth/order bound. Count1 and all-singleton words are treated separately.

The bank enumerates every labelled orientation through order5 and every root before deduplication by the full six-constant signature. The order-at-most4 subbank is exact because representatives are chosen in increasing order. The bounded coverage statement and H173 bank separation are appropriately finite-target, all-depth conclusions. Neither the proof nor the receipt establishes cofinite coverage, a repeatable digit compiler, or the unrestricted tournament spectrum.

## Independent numerical evidence

`check.cpp` has no author-code dependency. It computes all path counts by reverse subset endpoint recurrences, forms both reserved-bin responses by an explicit subset convolution, and regenerates all1099 labelled graphs and5405 roots. The complete set of65 signatures, including the19 order-at-most4 signatures, matches exactly.

Its search uses the original coordinates (h,e,R,d), with the signed physical row Hh-Se+LR and R'=CR-Qd, rather than the author's positive-coordinate implementation. It independently reproduces the entire reachable state sets, not just their scalar images:

- Order-at-most4, bound2000:20389 states, maximum representative depth10, missing odds7,21,173. Sorted positive-coordinate image hash869d18e3ef524da6e204e2abc10834067d808fe98a58fd61ccf21662359a74e0.
- Order-at-most5, bound5000:152044 states, maximum depth12, missing odds7,21. Image hashac84c8a13f6505432cd89c52cb36718178d38c99f790d322bceddafb8980fd9d.

All841 actual rooted source/gadget pairs through order3 satisfy the full four-state transfer. The sixteen rank-matrix entries are independently counted by literal Hamilton-word permutations, and a separate determinant expansion gives16. All eleven supplied witness graphs have their complete four-states independently counted, and a second walk inclusion-exclusion method verifies their Hamilton counts. The comparison script reconstructs every witness word with the old root fixed at label0 and checks all476 final arcs.

Optimized C++ and an unoptimized build with standard-library assertions and undefined-behavior sanitizer give identical receipts. Python comparison scripts agree in normal and optimized modes. Source normal/-O replays reproduce the frozen source receipt exactly.

## Reproduction and hashes

From this audit directory:

    g++ -std=c++17 -O2 check.cpp -o /tmp/check-retained
    /tmp/check-retained . > independent-normal.json
    python compare.py > comparison-normal.json

For debug controls, compile with `-O0 -D_GLIBCXX_ASSERTIONS -fsanitize=undefined`, rerun, then use `python -O compare.py`. The checker writes its independently regenerated bank and sorted state CSVs; these are outputs, not required inputs. It reads only witnesses.txt. compare.py additionally reads source/retained_sewing_normal.json to compare complete sets and literal recipes. Copy the named source files before replaying their verifier. No external packages are needed.

Exact audit hashes:

- check.cpp:52f34043daf0eac4933bfd8234da400d4e6c6d081d8831aa45e485501539b517.
- compare.py:767d46033e4b8116825a46c7ee7fd0ebaa03440cc0fceb1243bbe39110d7379e.
- witnesses.txt:fb645d107a2a10eef03e4b14d029c92d0878bce1932b35fc217f4f9d3f78b440.
- independent-bank.csv:1726317a5b4ccaf0a1193a712d5881ef9f7085e70c23054d5e3279b3a7ad5260.
- independent-normal.json and independent-debug.json:66b9eeb9137eddb093480b3bbeecaa64be064d4b692ac35cdbee4c10473b0f88.
- comparison-normal.json and comparison-optimized.json:27f7f74c0866aad3d0b92cb5bdd7845116607f412ba83376edd973918dcc0f1f.
- source-normal.json and source-optimized.json:2b1b3ae04c03b8ff8a3985e0a6bb4227acbb1ec0e544083bb60013fa819111eb.

No correction, denied action or unresolved blocker remains.
