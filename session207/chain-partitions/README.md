# All-order Moon-chain spectra and recognition

PROOF.md classifies the Hamiltonian-path count spectrum of Moon interval chains
at every order, all scalar collisions, and the distinct family after source
padding. Its cardinalities are explicit differences of partition numbers.
The same proof gives polynomial-time recognition and the unique chain order.
These are results for a physical subfamily, not the full tournament spectrum.

Dependencies are the known Moon product law and cited Bang–Zsigmondy theorem;
prior repository constructions and the checked primary source are identified
in the proof. REVIEW.md records independent audit.

Replay from this directory (standard Python only):

    python3 verify.py > /tmp/session207-chains-result.json
    cmp RESULT.json /tmp/session207-chains-result.json
    python3 independent_audit.py > /tmp/session207-chains-independent.out
    cmp independent_audit.out /tmp/session207-chains-independent.out
    sha256sum -c SHA256SUMS

The scripts check physical Hamiltonian paths, partition generating functions,
scalar-product injectivity, the exceptional relation 9=3^2, and greedy
recognition against separately generated finite monoids. Infinite conclusions
come from the proof, not finite extrapolation.
