# Finite pure-blue atoms and exact rigid-strong classification

Status: PROVED from CITED published inputs, with FINITE-EXACT certificate coverage through31. See [the full proof](FINITE_ATOMS_AND_RIGID_CLOSURE.md) and [independent audit](AUDIT.md).

## Results

- Every nontrivial strong pure-blue tournament has odd order at most131. The orders121,127,129 are separately excluded by exact count bounds.
- The only nontrivial rigid strong tournament is C3, in every order, closing the earlier above62 gap.
- Every skew circulant presentation at odd orders7 through31 has an explicit non-automorphic Hamiltonian-cycle successor and a reflection-failing cut path. Together with the strong-cycle reduction, the strong pure-blue alphabet through31 is exactly {1,C3,C5}.
- If B(z) counts the remaining finite pure-blue strong alphabet, including the singleton, then the ordinary generating function for all pure-blue tournament isomorphism classes, including the empty object, is (1+B(z))/(1-z²-z⁶).
- At every even order2m, the number of pure-blue classes is c(m), with c(0)=c(1)=c(2)=1 and c(m)=c(m-1)+c(m-3). This is an all-order theorem, not a census conjecture.
- For every n>=132, the full pure-blue counts obey p_n=p_(n-2)+p_(n-6). The odd recurrence begins at133. The unknown numerator has possible extra support only at33,35,...,119,123,125,131.
- Every carousel tournament of odd order at least7 has a symbolic non-automorphic-cycle certificate.

Completeness of the pure-blue alphabet {1,C3,C5} beyond31 is not resolved by this packet. The finite upper bound does not claim that the residual47 possible sizes have been checked. Vertex transitivity, regularity and circulancy are necessary but not sufficient for purity; nonrigidity does not imply nonpurity.

## Reproduction

Run from this directory using standard-library Python3:

    python3 independent_cycle_certificates.py cycle_certificates.json
    python3 independent_rigidity_cutoff.py
    python3 probe_cycles.py
    python3 extend_circulant_certificates.py --max-order 25 --seconds 45
    python3 extend_circulant_certificates.py --min-order 27 --max-order 29 --seconds 45
    python3 extend_circulant_certificates.py --min-order 31 --max-order 31 --seconds 45 --node-cap 2000000

The time and node guards yield INCOMPLETE, never a false completed exclusion. Runtime depends on hardware; the supplied final files cover every mask in each claimed order. Certificate verification is independent of generation and much cheaper than searching. The initial guarded31 attempt in extend_n31.out is incomplete; extend_n31_retry.out records the completed retry. Timings in output are descriptive, not deterministic invariants.

## Dependencies and sources

See [DEPENDENCIES.json](DEPENDENCIES.json) for exact source roles. The new rigid-strong classification uses Benhocine–Wojda1983's published tournament bypass theorem, avoiding the exceptions of later general digraph bypass theorems. The order131 bound uses Thomassen1985 plus Dixon1967 and an inverse-pair counting lemma; it is not inferred numerically from Cuckler2007's asymptotic. Publisher abstract access and full-text limits are disclosed in the proof.
