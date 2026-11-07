# Q/Z graph commutation and counting limit

Read COMMUTING.md. Q and Z commute up to rooted graph isomorphism, including transient sources. Words collapse to two nonnegative counters, and valid compiler outputs from any fixed finite source bank have counting function O((log N)^2). The rooted directed triangle attains the sharp order Theta((log N)^2), with an explicit formula for every output. This does not resolve whether a full synchronized paired source exists, nor the unrestricted odd Hamiltonian-path spectrum.

Run the self-contained checks:

    python3 verify_commuting.py
    python3 -O verify_commuting.py
    python3 verify_sharpness.py
    python3 -O verify_sharpness.py

Expected receipts are VERIFY.json and VERIFY_OPT.json. The checker includes the singleton for graph isomorphism, but uses order>=2 for the inherited transfer formulas. It performs 4,800 adjacency-isomorphism checks, 3,552 algebra/word checks and 592 literal endpoint-DP checks.

Only the files listed in PUBLICATION_FILES.txt form the publication packet; other local exploration files are not proof dependencies. Independent bridge-audit material is supplied separately for the bridges packet.
