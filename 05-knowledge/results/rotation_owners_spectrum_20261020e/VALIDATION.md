# Validation and scope

The seven executable checks use explicit exceptions, so optimization does not
remove their predicates. Normal and `python -O` runs must both reproduce the
saved stdout exactly after CRLF normalization. `verify.py --replay` enforces
this and authenticates the recorded source and inherited data hashes.

| Check | Explicit universe and independent or hostile controls |
|---|---|
| black_check.py | All2046 marked circulant step sets at odd orders5..23;2034 interval-rotation controls;3303 replays over arbitrary fixed-HP middles of order at most6; threshold orders7..101; Paley11 extraction hostile and order3 boundary |
| owners_check.py | Every2244 primitive record at kernel budget2 and34698 at3, before pair-free/owner filtering;7480 one/two-edge bulk completions, accepted and rejected categories;469 independent large-light substitutions;168 leading kernels and100 exact difference identities |
| owners_independent.py | Independently written adjacency-matrix and direct-permutation implementation of all2244 budget-two kernels; no repository imports; reproduces C2 with positive triangle and hostile transitive-pair controls |
| regular_check.py | All3 regular classes at7 and15 at9; full endpoint matrices by subset DP and literal DFS; independent labelled degree-demand counts versus orbit masses2640 and3230080; pairwise nonisomorphism and full automorphism enumeration; non-Hamiltonian-connected orders3/5 |
| spectrum_check.py | All33867 labelled tournaments at orders1..6, independent literal path counts through5; all14 carousel7 triangle reversals; all576 ordered labelled regular5 pairs for cycle-flip transport; exact cutoff and arithmetic checks |
| arithmetic_check.py | All7920 SCC words and parity-path matchings through18; every preserved triangle weight;64 literal Hamiltonian controls through8; marked-edge273 partitions; operation fibres through200;3600 Euclidean inverses; braid matrix kernel and linking hostile |
| purity_check.py | Every65534 skew presentation at odd orders3..31; explicit nonblue certificates in all65528 negative cases; all Hamiltonian words of the six positive order3/5 presentations; symmetric SCC grammar through32 |

The main owner checker imports the already published exact ordered-template and
local-owner implementation. Its new kernel extraction, completion and symbolic
weight controls are separate; its counts at budgets2/3 agree with the earlier
full capped-class census. The all-order claim rests on the proof in owners.md,
not just these sampled completion controls. The remaining C4,C5,C6 coefficient
banks have not been enumerated. The independent checker reconstructs all2244
budget-two primitive kernels and their owners without any author-code imports;
it separately reproduces `C2(y)=756+792y+270y^2+30y^3`.

The regular checker alone needs NetworkX. Its order-nine input is the
authenticated inherited census bank; an independent labelled count and exact
orbit-mass equality certify completeness of the selected regular classes.
The large-order Hamiltonian-connectivity result uses Thomassen's cited theorem.
The publisher abstract was checked; direct original article retrieval returned
403, so no audit of that paper's original proof is claimed.

Separate-agent mathematical reviews covered both black macros, the proper-
interval lemma, the owner proof's lexicographic and weight invariance, and the
matching bijection/63+210 partition. Parent review checked the kernel bounds,
leading coefficient, witness cutoff and scope. No Lean formalization is claimed.

Reproduce from the repository root:

```text
python 05-knowledge/results/rotation_owners_spectrum_20261020e/verify.py --replay
```

The manifest hashes packet sources, notes and receipts, as well as inherited
Python/data dependencies. `--record` refreshes it intentionally after edits;
ordinary verification never rewrites it. The manifest is an integrity record,
not evidence of mathematical correctness.
