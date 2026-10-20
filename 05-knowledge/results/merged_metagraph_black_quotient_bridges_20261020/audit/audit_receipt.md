# Independent audit: paired-black bridges and blue-fibre hostiles

Status: general proofs audited; independent exact graph reconstruction through order six; direct order-eight witness verification. This audit does not claim an independent second exhaustive census at orders seven or eight.

## General statements

The following arguments pass:

- Reflection-paired black multigraph Q has degree (merged fibre mass minus blue mass)/2. A diagonal entry of the black half-edge matrix contributes one Q loop per four, so loop conventions and all cut-parity formulas are correct.
- Every non-self-converse node has odd Q degree. The charge/cut laws and the cycle-space identity beta(K)=|E(Q)|+beta(Q) hold with loops and isolates included.
- The two explicit leaf families give, for every n>=5,
  bridges(Q_n) >= (a_n-p_n)/2 + [n even] a_((n-4)/2),
  where a_0=1, a_j=0 for j<0, a_n=a_(n-1)+a_(n-3), p_(2r)=a_r, and p_(2r+1)=a_r+a_(r-1).
- Unique ordered SCC decomposition proves distinctness within and between the two families. The triangle-score separation proves that no two constructed leaves share their edge when n>=5. The n=4 exception is essential: the two constructed leaves share one bridge.
- No all-order equality, connectivity, or classification of all Q bridges follows. Those assertions remain separate from the proved constructive bound.

## Independent finite reconstruction

`audit_structural.py`, containing the separate brute-permutation classifier, rebuilds the full fixed-path cube and graph without importing producer code or an atlas. It obtains Q-bridge counts 0,1,1,3 for n=3,4,5,6. All loop divisibilities pass. The n=5 regular pure-blue node has fibre size three, and the independently reconstructed complement incidence has three distinct neighbours. Thus unrestricted pure-blue-implies-pendant prose is false; the proved parity results do not depend on that prose.

The scripts also reproduce the small pure-blue inventories. They do not claim completeness of the strong pure-blue alphabet at arbitrary orders.

## Independent order-eight cap refutation

`independent_n8_witness.py` directly examines every one of the 8! vertex permutations for each of two explicitly decoded tournaments. It enumerates Hamiltonian paths, tests automorphisms, and tests path-reversal anti-automorphisms. No canonicalizer, stored class census, or producer helper is used.

In the declared explorer mask gauge:

- Mask 40077: H=637, Aut=1, fibre=637, blue Hamiltonian paths=39, blue fibre=39.
- Mask 40109: H=621, Aut=3, fibre=207, blue Hamiltonian paths=117, blue fibre=39.

Both have score multiset (3,3,3,3,4,4,4,4) and c3=20. A blue path directly certifies self-converseness. Therefore each is an exact counterexample to the order-eight value 27 predicted by C1 in `THM-643-blue-black-line-parity-structure.md`. This is a refutation of that conjectural cap, not of THM-643's proved parity theorem. The audit establishes blue fibre 39 for these witnesses; it does not independently establish the global order-eight maximum.

The output JSON supplies explicit adjacency lists, so the witness is not dependent on remembering the mask convention.

## Reproduction

Keep the audit Python files in one directory. Python 3 with the standard library is sufficient.

    python3 audit_structural.py
    python3 independent_n8_witness.py

Both scripts are self-contained and run without producer files. Their `.json` and `.out` outputs record the exact checks. The structural script consolidates the audited classifier and graph checks without changing their mathematics.
