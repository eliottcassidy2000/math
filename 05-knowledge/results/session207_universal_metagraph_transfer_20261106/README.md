# Universal tribonacci transport of metagraph incidence

Fix any tournament core T and any directed Hamilton path P. Pad T by a transitive block of size a before it and another of size b after it, and complement every arc outside the extended path. The resulting neighbour count satisfies an exact tribonacci recurrence independently in a and b. Nine corner counts, with 0<=a,b<=2, determine the entire infinite array.

The proof works for arbitrary pointed cores. It classifies each odd cycle crossing an attached tail as a tail suffix followed by a simple directed core path. The inherited odd-cycle collection formula then gives the rational generating function. This is not a fitted recurrence.

Further exact consequences:

- For a,b>=3 the target automorphism group is the original neighbour's automorphism group fixing both old path endpoints.
- One-sided padding stabilizes the group fixing just the first endpoint. Thus a fixed source count and ordinary fibre degree coexist with neighbour degrees following an exact tribonacci sequence at every sufficiently large consecutive size.
- Nine counts and one endpoint-stabilizer order for each original Hamilton-path orbit determine the degree-labelled incidence measure along the two-sided padding family. Parallel incidences and loop flags retain their proper multiplicities.
- Each one-variable sequence has minimal rational recurrence order exactly three. Its growth is asymptotic to a positive constant times the tribonacci root to the padding length.

PROOF.md contains full statements, boundaries and dependencies. REVIEW.md records independent mathematical review. This does not classify all cores, prove global spectrum completeness or connectivity, identify targets by count, or assert count-only lumpability.

Reproduce from this directory:

    python3 check.py > replay.json
    cmp replay.json RESULT.json
    python3 -O check.py > replay-opt.json
    cmp replay-opt.json RESULT.json
    python3 independent_audit.py > independent-replay.out
    cmp independent-replay.out independent_audit.out
    sha256sum -c SHA256SUMS

Main controls cover1710 pointed-tail counts,300 two-sided counts,285 one-sided and48 two-sided stabilizer checks, and four nontrivial order-three stabilizer examples. Its explicit checks survive optimization. The independent checker covers203 pointed-tail counts using a separate recurrence implementation; it uses assertions and must run without -O. All scripts require only Python3. Finite checks supplement the all-order proofs.
