# Reproduce the independent Seidel-energy review

The exact scope and frozen producer pins are in REVIEW.md and producer_receipt.json. This packet includes the structural all-order theorem, the substitution/profile-dominance extensions, and the explicit odd-carousel family, but no future appendix.

Use standard Python and the existing C++17 compiler; install nothing. Run from this directory:

    python3 audit.py > normal.stdout
    python3 -O audit.py > optimized.stdout
    cmp normal.stdout optimized.stdout
    python3 physical_audit.py > physical.normal.stdout
    python3 -O physical_audit.py > physical.optimized.stdout
    cmp physical.normal.stdout physical.optimized.stdout
    python3 ambient_audit.py > ambient.normal.stdout
    python3 -O ambient_audit.py > ambient.optimized.stdout
    cmp ambient.normal.stdout ambient.optimized.stdout
    python3 carousel_audit.py > carousel.normal.stdout
    python3 -O carousel_audit.py > carousel.optimized.stdout
    cmp carousel.normal.stdout carousel.optimized.stdout
    python3 extensions_audit.py > extensions.normal.stdout
    python3 -O extensions_audit.py > extensions.optimized.stdout
    cmp extensions.normal.stdout extensions.optimized.stdout
    python3 producer_audit.py ../seidel_energy > producer-audit.normal.stdout
    python3 -O producer_audit.py ../seidel_energy > producer-audit.optimized.stdout
    cmp producer-audit.normal.stdout producer-audit.optimized.stdout

Adjust the explicit producer directory in the last two commands to the actual frozen packet. The local default is the sibling session189-seidel-energy. Run physical_audit.py before ambient_audit.py and extensions_audit.py, which reuse its freshly built physical executable. Producer replay occurs in a temporary copy and never changes the producer packet.

PUBLICATION_FILES.json is the exact allowlist. Binaries, caches, input reading copies, and the unlisted empty normal.json are excluded. All complete formal banks, all six complete physical banks, literal masks and query values, normal/optimized receipts, and source pins are retained. The formal child-uniform banks are not called full physical banks for large inflated tournaments.
