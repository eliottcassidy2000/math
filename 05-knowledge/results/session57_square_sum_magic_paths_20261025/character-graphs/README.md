# Character graphs versus magic incidence

Main results: neither direction of magic ↔ square-sum Hamiltonian-path existence holds (F59 and F17). In characteristic3 every magic entry set has square-sum graph K9 or empty according to χ(−1), although exactly432 of its9! arrangements are magic. The exact positive encoding is the eight-edge grid hypergraph embedded in the triple-sum hypergraph H_T.

Read THEOREM.md. Reproduce with:

    python analyze.py
    python check.py

The producer needs SymPy; the independent checker uses the standard library only. results.json exhausts all normalized magic grids at odd prime-power orders≤101. witness-729.json and witness-2187.json add explicit extension-field controls. sources/ contains pinned prior arithmetic and the prior census, not copyrighted paper text. SHA256SUMS records local inputs and outputs.
