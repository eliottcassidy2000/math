# Fortunate / square-sum pattern audit

Read [THEOREM.md](THEOREM.md) for the result and its scope.

- The prime pattern is the run-length word of F_n=p_(n+1), not F_n=p_n+1 and not raw numerical differences. Its equality indices are known as OEIS A035346.
- On square-sum graphs10..26, ordinary Hamiltonicity misses the first three bits. Componentwise traceability, equivalently pc(G)=c(G), matches all17 bits.
- The match fails immediately at n18↔G27. Both ends of the selected graph run window are truncated.
- A general endpoint-aware path-cover birth law explains the graph transitions without implying any arithmetic bridge.

Reproduce with Python3 standard library only:

    python3 check.py
    python3 -O check.py
    python3 test_rejections.py
    python3 -O test_rejections.py

The artifact retains alternative interpretations, all paths through28, graph births, components, minimum-cover witnesses, independent no-path certificates, and prime certificates. No all-order square-sum existence theorem is needed. `MANIFEST.json` records hashes; `SOURCES.json` pins inherited material and primary sequence sources.
