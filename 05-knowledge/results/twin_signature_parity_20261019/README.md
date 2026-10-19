# Twin signatures, oriented-tree automorphisms, and a minimal parity counterexample

Read [the proof](TWIN_SIGNATURE_PARITY.md) and [independent review](review/INDEPENDENT_REVIEW.md). A nine-vertex tree refutes the converse of hereditary symmetry implying goodness; the sufficient direction survives. Twin signatures give exact subgroup divisibility and rooted orientation-capacity recurrences. The minimality certificate through eight vertices and the five-vertex tournament hostile retain their stated finite scopes.

Run from this directory:

    python3 -B verify_twin_signature.py
    python3 -B -O verify_twin_signature.py
    python3 -B review/independent_tree_capacity.py
    python3 -B -O review/independent_tree_capacity.py
    python3 -B review/independent_embedding_witnesses.py
    python3 -B -O review/independent_embedding_witnesses.py

Match author outputs with twin_signature_receipt*.json, capacity outputs with review/capacity_receipt.json, and embedding outputs with review/embedding_receipt*.json. The witness certificate is beside each checker that consumes it. The independent capacity audit has 5,464 explicit checks through order eight; embedding checks cover all eight minimality witness pairs and all 1,024 labelled five-tournaments. The author's larger orientation census remains separately identified. SOURCE_PROVENANCE.json pins reviewed and published bytes.
