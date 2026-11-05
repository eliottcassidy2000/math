# Replaying the independent graph-census audit

The independent sources and receipts here are unchanged from the audit bound in ../REVIEW.md. This is a different negative-check route from the main proof's all-length port closure. The main proof does not need the strong-order bound or census.

From this directory:

    c++ -std=c++17 -O2 probe_c3_order13.cpp -o census
    ./census records.csv > inventory-replay.jsonl
    cmp inventory-replay.jsonl c3-order13-inventory.jsonl
    c++ -std=c++17 -O2 probe_lowcycle_821.cpp -o filter
    ./filter records.csv > hits-replay.csv
    cmp hits-replay.csv lowcycle-821-order13.csv
    python check.py
    python -O check.py

Compare the Python stdout with controls-normal.json. The large records.csv is regenerated rather than stored; its expected SHA256 is 71aafd0283d5fe4f93406ce65c3f546ded04ed2dacb951da9532e92259579779. Assertion-enabled unoptimized census/filter reference receipts are included. Both C++ sources and the Python checker are standalone and use no external packages.

C3_COMPLETE.md is the unchanged mathematical census-completeness dependency (SHA a919deb88bc78a84b8d0cf40199fc7966051541f1e7027be3666bb9f0713b444). Its audited dependency bundle is also available at the immutable checkpoint:
https://github.com/eliottcassidy2000/math/tree/231bda1d4b2f57506bfcd6a2b2b3decb5f5f089d/session197/lowcycle/feedback533/source/session197-lowcycle

For the main capped-transfer proof, run the parent directory's verify.py instead. It regenerates the unchanged physical table and reproduces the complete all-length closure without the census data.
