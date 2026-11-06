# Session 206: modularly weighted O/Z parameter capacity

Status: independently audited all-order theorem for fixed physical sources and fixed odd moduli. The count is of tail parameters with multiplicity, not distinct output values.

For a source (h,R,e,d), R>d>0, the even-input and odd-input capacities are respectively

c0 = log(R/(R-d/2))/d,
c1 = log((R-d/2)/(R-d))/d.

For odd p and g=gcd(R,d,p), parity b contributes capacity (g/p)cb to each scalar residue congruent to h-eb modulo g, and zero elsewhere. The proof combines parity-specific residue counting with the primitive binary digit chain and endpoint-cylinder argument.

The resulting weighted finite-bank criterion excludes the physical bank containing the root, JJ=(5,6,2,3), and OZJ=(11,12,4,5). It misses at least X/252-o(X) odd targets at most X in the class 1 modulo 6. This exclusion is not obtained by the earlier unweighted capacity bound or by an omitted modular class: OZJ itself has every parity-compatible modular residue.

This does not exclude every finite physical bank, supply a covering bank, prove an output-density limit, or provide a uniform estimate for growing moduli or banks. No general even-modulus extension is claimed in this checkpoint.

## Contents and reproduction

- WEIGHTED_CAPACITY.md: complete proof, unchanged from the independently reviewed draft.
- verify.py and RESULT.json: 2,688 exact finite homomorphism-fiber checks over 21 physical source states; exploratory parity/residue input counts.
- audit/REVIEW.md and audit/REPLAY.json: frozen independent review and optimized verifier replay.
- SHA256SUMS: frozen payload manifest.

From this directory, using Python 3 and its standard library:

    python verify.py > /tmp/weighted-capacity-replay.json
    cmp RESULT.json /tmp/weighted-capacity-replay.json
    sha256sum -c SHA256SUMS

The finite checks are regression controls, not proofs of asymptotics. Convergence rates are not asserted uniformly in the source or modulus.

## Audited prerequisites

- [Exact total parameter capacity](https://github.com/eliottcassidy2000/math/blob/ea15da8016b2623c3f10365202301f492a7e9a84/session206/parameter-capacity/PARAMETER_CAPACITY.md).
- [Root residue capacity and the primitive digit chain](https://github.com/eliottcassidy2000/math/blob/16e5d8ee8194eeda2ab79ca22dbf6e96eb6e5a26/session206/root-residue-capacity/ROOT_RESIDUE_CAPACITY.md).
- [Root collision and d=1 collapse](https://github.com/eliottcassidy2000/math/blob/e5686802c8998e01db36cf714d98793ef2b13bdd/session206/root-oz-noncofiniteness/ROOT_OBSTRUCTION.md).
- [Exact modular images of physical sources](https://github.com/eliottcassidy2000/math/blob/1635b0bfdb762eeea7f0c93df273861e722d48bc/session206/exact-modular-image/DYADIC_SURJECTIVITY.md).
