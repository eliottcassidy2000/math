# Session40 actual endpoint and shared-incidence certificates

The main result is a sharp marked reconstruction rule: at most six unmarked vertices require at most one shared-selector chirality bit; at most five require none. At order eight the actual two Hamiltonian endpoints meet the bound. Full endpoint enumeration finds18 surviving marked-refinement collision groups, all repaired by the signed rectangle J=±3, retained modulo8 and even under converse. A regular-seven module gives a sharp next-order obstruction. This authenticates target classes, not individual literal partner ports.

- ENDPOINT_CHIRALITY.md: proof, finite census, ownership and converse scopes.
- PATH_TRANSPORT_AND_FOREST_RULE.md: exact source/path transport, all-size linear marking-budget rule, cheaper conditional forest rule, sharp order-nine boundary.
- scan_endpoints.cpp and endpoint_scan.txt: complete order-eight endpoint census using exact strings and subset DP.
- analyze_collisions.py and collision_analysis.json: all18 groups, actual Hamiltonian witnesses and matching signs.
- check.py and receipt.json: independent matrix-level consumer, old literal-switch repairs, four-point controls, exhaustive fixed-path transport through order six, regular-five uniqueness and order-nine hostile.
- check.normal.out and check.optimized.out: frozen consumer replay logs.

## Reproduce

The C++ scan uses the Session38 certified ordinary-class inventory:

    g++ -O3 -std=c++17 scan_endpoints.cpp -o /tmp/scan_endpoints
    /tmp/scan_endpoints PATH_TO_SESSION38_WIRING/classes_input_n8.tsv > /tmp/endpoint_scan.txt
    cmp /tmp/endpoint_scan.txt endpoint_scan.txt

The Python consumer uses the inherited binary-counterfeit receipt:

    python -B check.py --counterfeits PATH_TO_SESSION38_BINARY_COUNTERFEITS/receipt.json
    python -O -B check.py --counterfeits PATH_TO_SESSION38_BINARY_COUNTERFEITS/receipt.json

The default input path is the original workspace sibling session38-binary-counterfeits. The explicit option is portable; no network or nonstandard Python packages are needed. The C++ scan and consumer serve different purposes: the scan enumerates every actual feasible marked pair, while the consumer independently checks all surviving collision witnesses and mathematical transport controls. The normal/optimized check count is71,202. No removable assertions are used.

## Pinned inheritance

Session38 full wiring and binary-counterfeit packets were published at e5d18e55f4c321f7c34fc67924782d04a3b5fed2. Session39 equitable-triple phase decoder was published at a538127eefd3f7575f2b54b629b43624d2e64ab2. The preceding marked-cell route was incoming e41a2c82b34fea4dffc6c1ee5dab07bfb8c56ded. These are inherited inputs, not new claims of priority. The all-size affine holonomy classification itself is inherited; this packet's new finite carrier is the actual endpoint census and bounded residual-cell budget.

## Full literal-wiring replay

The stronger replay builds the marked-signature lookup only from representative adjacency and feasible endpoint pairs. It then independently decodes all2^21 literal fixed-path tournaments and compares every ordinary target class to the inherited mask ledger. Both source and target classes in the regenerated wiring are taken from the new decoder. Exact ordinary port and merged edge counters are compared, including loops and blue/black multiplicities, and every mask's converse reflection is checked.

    python -B replay_all_masks.py --wiring PATH_TO_SESSION38_WIRING
    python -O -B replay_all_masks.py --wiring PATH_TO_SESSION38_WIRING

This requires g++ supporting C++17. All decompressed inputs and the compiled binary stay in a temporary directory. It does not build classIDs into the signatures and does not use the reference mask ledger to infer a missing class. The reference is only a comparison oracle after decoding. The literal complement mask is retained as actual path ownership; this is not reconstruction of literal partner identity from the compact signature alone.
