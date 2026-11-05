# Independent all-order H67 audit

2026-11-05. **PASS for the primary H67 classification packet.**

Frozen hashes:
- PROOF.md: `ae7fde9e30af17d8f97836db4302a4afa629e35ca49b099fd428438c79d2a163`
- verify.py: `f2d2060bf9e1651ca5f677146d6a020279619a94f0f8696320c80a2b1726a0b0`
- certificate.json: `0301a1ef4ac75bac492553cbef941fcde2a5c49292725c842549fe45408126c8`

This signoff excludes the separate ROOTED.md supplement.

## Published order bound

Busch's official paper was independently opened from
https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf.
Theorem 1 on PDF page 3 gives the three residue-class lower bounds in terms
of beta=5^(1/3). Each implies `H(T)^3 >= 5^(n-1)` for strong order n>=3.
The arithmetic `67^3=300763 <390625=5^8` is exact and implies n<=8.
The theorem is applicable to directed Hamiltonian vertex-word counts, with
the same convention as the packet. The paper also gives the exact strong
order-six minimum as 15; no bound of 17 is used here.

## Independent complete census

`census.cpp` generates the complete tournament universe inductively by
adjoining each possible incidence pattern. It uses a separate exact
canonicalization: vertices are partitioned only by outdegree, and every
permutation within those degree classes is tried. There is no color
refinement and no dependency on the source canonicalizer. Hamiltonian
paths and cycles are counted by predecessor-removal subset programs;
strongness is checked by repeated reachability expansion.

It independently reproduces the inventory
`1,1,2,4,12,56,456,6880` at orders 1 through 8, all C1/C2 class counts,
exactly six strong H67 classes at order 7, and exactly twelve weak H67
classes at order 8. No strong H67 class occurs at any other order <=8.

`check.py` separately calculates automorphism orders and sums `n!/Aut(T)`
over every independent representative. The resulting labelled orbit
masses are exactly

`1,2,8,64,1024,32768,2097152,268435456`,

which equal `2^(n choose 2)` at every order through 8. Together with exact
canonicalization, this is a strong independent check of full coverage and
absence of duplicate isomorphism classes.

For each of the six certified cores, literal permutation enumeration
confirms H=67 and the claimed cycle count. The cycle distribution is
3,4,3,4,4,3 in certificate order. Every core has automorphism order one,
no transitive single-vertex deletion, and at least one transitive
two-vertex deletion. All saved two-apex presentations were reconstructed
using topological source stripping of the base, rather than the source's
degree-sort procedure. Their numbers are 2,2,3,1,1,2. The converse pairing
is exactly 1<->6, 4<->5, with 2 and 3 self-converse.

The twelve independent order-8 targets are exactly the six source paddings
and six sink paddings. The source verifier was also rerun normally and
under `-O`; stdout and generated certificates are byte-identical.

## All-order proof checks

Hamiltonian paths traverse the tournament SCCs in their unique quotient
order, in contiguous blocks, so the path counts multiply. Since 67 is
prime and nontrivial strong components contain a triangle, exactly one
component has count 67 and every other SCC is a singleton. The latter fact
uses an actual injection of the three triangle paths into full Hamilton
paths, not an unsupported small-order extrapolation.

The strong component is therefore one of the six exhaustively classified
cores. Every tournament has unique form `Tr_p -> Q_i -> Tr_q`. For fixed
n>=7, there are n-6 choices of (p,q) and six choices of core. Directed
SCC order is preserved by isomorphism, so different choices do not merge,
even when the core has a converse partner. This proves `6(n-6)` unlabelled
classes and, by rigidity, `6(n-6)n!` labelled tournaments.

Deleting a certified two-vertex core feedback set makes the whole ordered
join transitive. A one-vertex deletion cannot make the induced core
transitive, so the feedback vertex number is exactly two at every order.
The minimum cycle count is three only in the strong H67 subclass; adding
any source or sink destroys strongness and leaves zero Hamiltonian cycles.
The packet keeps this distinction explicit.

No mathematical or implementation issue was found. The source
canonicalizer's stop rule is also sound: inclusion of the old color in
each signature prevents merging, so unchanged class count means unchanged
partition, and its final `nc` labels remain isomorphism invariant.

## Reproduction

In this folder, run `g++ -O3 -std=c++17 census.cpp -o census`, then
`./census > independent_census.json 2> census.log`, then `python check.py`.
Results are in `AUDIT.json`, `independent_census.json`, and `classes*.txt`.
All author source files were left untouched.
