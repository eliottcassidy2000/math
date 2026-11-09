# Session234: Hamiltonian-path counts 7 and 21 beyond tournaments

**Status: REFUTED (the proposed matrix condition); PROVED (complete all-order
semicomplete classification and fixed-count finite grammar); FINITE-EXACT
(semicomplete census through six, unrestricted digraphs through five);
independently audited.** Date: 2026-11-09.

The conjecture is false. A binary adjacency matrix can have exactly seven
or twenty-one directed Hamiltonian paths while its diagonal is entirely
zero and no mirrored off-diagonal pair consists of two zeros.

The distinction is **tournament versus semicomplete digraph**. A tournament
has exactly one arc between each pair of different vertices. A semicomplete
digraph has at least one, so it may contain reciprocal pairs. The classical
tournament omissions 7 and 21 do not survive that extra freedom. Self-loops
are irrelevant to vertex-distinct Hamiltonian paths.

## Smallest counterexamples

Rows of a four-vertex matrix with exactly seven paths:

```text
0111
0011
1001
1000
```

On vertices 0,1,2,3 its paths are exactly
`0123, 1203, 1230, 1302, 2013, 2301, 3012`.

Rows of a five-vertex matrix with exactly twenty-one paths:

```text
01111
00110
00011
10001
11000
```

These orders are minimal even among arbitrary digraphs: at most `n!`
directed Hamiltonian paths exist; three vertices allow at most six, and a
four-vertex graph missing any one directed nonloop arc has at most
`24-6=18`, whereas the complete bidirected four-vertex graph has 24.

## Complete finite catalogue

| order | labelled semicomplete matrices | H=7 | H=21 | H=7 classes | H=21 classes |
|---:|---:|---:|---:|---:|---:|
|1|1|0|0|0|0|
|2|3|0|0|0|0|
|3|27|0|0|0|0|
|4|729|72|0|3|0|
|5|59,049|720|1,760|6|16|
|6|14,348,907|6,480|32,640|9|48|

Exactly three strong seven-path cores occur, all of order four. There are
sixteen strong twenty-one-path cores of order five and sixteen of order
six. The latter include two classes with only **one** reciprocal pair.
Matrices, paths, reciprocal pairs, tournament resolutions and orbit sizes
are saved, not merely aggregate counts.

- [Finite census report and reproduction commands](census/REPORT.md)
- [Order 1–5 catalogue](census/results.json)
- [Order-six catalogue](census/classes6.json)
- [The two strong three-path cores](census/strong_h3.json)
- [Complete proofs and counting formulas](PROOFS.md)
- [General fixed-count grammar theorem](GENERAL_GRAMMAR.md)
- [Independent adversarial audit](audit/AUDIT.md)
- [Standalone matrix generator](census/generator.py)
- [Validation and full replay instructions](VALIDATION.md)

**Important hostile control:** a directed triangle with one reverse arc
added still has three Hamiltonian paths. The added arc is unused. Thus
strong connectivity does not make arc addition strictly count-increasing.
Both three-path cores must be retained in a twenty-one-path product.

## Exact all-order generation

For seven paths, choose one of the three strong four-vertex cores and put
any number of singleton vertices before and after it in a total component
order. For twenty-one paths, either choose one of the 32 strong cores and
pad it in the same way, or choose one seven-path core and either of the two
three-path cores, put them in either order, and place singleton chains
before, between and after them. Every cross-component arc points forward.

**These rules generate every loopless semicomplete example, at every
order, exactly once up to isomorphism.** Relabel for labelled matrices.
Choose diagonal bits freely for loop-allowed versions. The classification
is complete because strong semicomplete graphs contain spanning strong
tournaments, whose proved lower bounds restrict the core orders.

For order `n>=4`, seven has exactly `3(n-3)` loopless isomorphism classes
and `3(n-3)n!` labelled matrices. The twenty-one class counts are 16 at
order five, 48 at six, 92 at seven and 148 at eight, with an explicit
quadratic formula thereafter. [All formulas, including loops](PROOFS.md#5-exact-counts-and-loop-symmetry).

Run, from `census/`:

```bash
python3 generator.py --target 21 --n 8 --limit 3
python3 generator.py --target 7 --n 5 --with-loops --limit 3
python3 generator.py --validate
```

## What survives of the conjecture

The complete classification gives these sharp loop-independent necessary
conditions, writing `b` for the number of reciprocal pairs:

- Seven paths force an absent pair **or exactly two reciprocal pairs**.
- Twenty-one paths force an absent pair **or `1<=b<=3`**.

Exactly one reciprocal pair suffices for 21 at order six. Its two
tournament resolutions have 17 and 19 paths; their union has 21.
Tournament resolution preserves a *weighted* path count, not the original
count, explaining why the tournament omissions do not transfer.

## General theorem and unrestricted matrices

For **every fixed positive path total**, semicomplete graphs have a finite
strong-core grammar. A strong core has at least as many Hamiltonian paths
as vertices, so a core contributing to total `h` has at most `h` vertices.
Singleton gaps then give rational generating functions: loopless
isomorphism counts are eventually polynomial in the order; labelled counts
are `n!` times a polynomial. Allowing loops gives `2^n` times the respective
quantities, except that isomorphism coefficients must first use the actual
loop-decoration orbits. Exact degree is the maximum realizable number of
nontrivial components, at most the number of prime factors of `h` counted
with multiplicity. This is not a universal-spectrum claim.

For **arbitrary digraphs with missing pairs**, the complete scan through
order five gives 120 seven-path graphs at order four, and 49,440
seven-path / 7,760 twenty-one-path graphs at order five. An exact
path-support/hypergraph construction is given in [DG8](PROOFS.md#7-arbitrary-matrices-and-a-complete-if-expensive-generator).
The smaller strong-core catalogue is not claimed for this unrestricted
sector.

## Verification and boundaries

The semicomplete scan is independently checked by a different algorithm
that scatters each permutation into every compatible matrix. Every
histogram entry through order six agrees; all 2,097,152 order-seven
tournaments independently verify the required floor. Unrestricted target
strata agree through order five. Loop decorations and every generated
loopless representative through order eight are verified separately.

No completeness theorem for the remaining tournament spectrum is claimed:
the Chan–Pak conjecture remains **OPEN**. Classical floor, parity and
Hamiltonicity mechanisms are credited, not claimed as new.

Inherited reference:
[THM-1370 / h-spectrum-omits-7-21-all-n](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
The similarly numbered elliptic-Jacobian file is unrelated.
