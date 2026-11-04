# S190 planar Hamiltonian interfaces

Read **PROOF.md** for statements, proofs, scope, sources, and the connection ledger.

Key additions to this repository:

1. Exact small-interface specialization of standard matching-connectivity theory. Three/four/five ordered planar ports require respectively 6/26/100 live degree-and-pairing coordinates, with a unique complementary exterior state.
2. Every coordinate has a realizable sparse planar unit probe. Full profile equality is equivalent to equality in every admissible exterior; scalar degree counts are insufficient.
3. The first nonpermutation block occurs at six ports: K_2 plus K_3, determinant -2. Full six-port ranks are 365 over Q and 364 modulo two. This is response-operator blindness, not a claim that cycle parity tests existence.
4. Complete face-stellation gives infinitely many triangulations with linearly many pairwise non-touching reducible degree-three vertices but robust Hamiltonian and traceability failure.
5. Exact square-sum deletion sets: Q18–Q21 restore only on deleting18; Q22 on deleting2 or18; Q24 on deleting11,18,22,24. All124 deletion cards were replayed by two independent search methods.

These are independently derived and reproducible elementary results, not asserted literature-priority breakthroughs. The general interface and meander frameworks are established literature, including a newly checked October2026 paper. Current 4CT v3 hypotheses are distinguished from its older abstract formulation. Cited all-order square-sum existence is not independently audited; the nontrivial no-hypo corollary retains the Q2 convention exception.

Reproduce in this directory:

    python3 verify.py > RESULTS.json
    python3 -O verify.py > RESULTS.optimized.json
    cmp RESULTS.json RESULTS.optimized.json
    python3 independent_kernel.py > KERNEL_REPLAY.json
    python3 stellation.py > STELLATION.json

All programs are Python-standard-library only, with explicit validation exceptions. RESULTS contains all189 physical graph cases and all positive square-sum deletion path witnesses. This packet does not mutate Git or install software.
