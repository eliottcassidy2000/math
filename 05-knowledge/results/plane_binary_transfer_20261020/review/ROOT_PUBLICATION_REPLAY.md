# Publication replay and robust verification gates

Root reviewed the circle compatibility, stochastic rounding, Boolean degree boundary, finite norm identities, eight-point transport and finite certificate scope. The independently reviewed base results are unchanged. For portable checking under Python optimization, removable assert statements in the publication copy were mechanically replaced by explicit exceptions; original and new script hashes are recorded in ROOT_PUBLICATION_REPLAY.json. The q17 refutation certificate itself is byte-identical to the reviewed version.

Five exact replay scripts were executed both normally and with -O, from this copied packet; each pair produced byte-identical output and passed all retained gates. The output files are included. This is not a Lean build or a validation of the source blueprint's major claims. The proof-of-search generator is preserved as evidence; replaying the refutation verifier does not rerun the search.

The q17 verifier additionally rejects negative/out-of-range and Boolean indices and malformed tree records, with an optional --certificate argument. Three corrupted-index variants were rejected in both normal and optimized modes. Valid-certificate outputs remained unchanged.
