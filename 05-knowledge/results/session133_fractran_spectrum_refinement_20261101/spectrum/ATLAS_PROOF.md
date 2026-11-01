# Exact low Hamiltonian-count strong atlas

**Status: FINITE-EXACT fresh independent reproduction of inherited results;
CITED all-order cutoff. No novelty claim for the exclusions 7 and 21.**

**Correction:** checkpoint 5811bb4 incorrectly used H15 as a packing-polynomial
collision by comparing trailing zeros. The H9 control below repairs this;
[CORRECTION.md](CORRECTION.md) preserves the exact failure mechanism.

`H(T)` counts ordered directed Hamiltonian vertex words, not path orbits.
`build_atlas.py` produces every isomorphism class through order six, starting
with the empty tournament and adding every orientation of one new vertex.
Every tournament is reached because deleting its last vertex gives a previous
class; extending an isomorphism of that deletion transports its incident arcs.
Each generated adjacency is canonicalized by checking every vertex permutation.
No sample or unverified catalogue is used.

For each of the 76 resulting classes, the program compares subset-DP path
counts with direct enumeration of all Hamiltonian vertex permutations, and
with the literal odd-cycle collection formula. Directed cycles are quotiented
by cyclic rotation only. Distinct cycles on the same support stay distinct.
It independently counts automorphisms and verifies the orbit masses
`1,2,8,64,1024,32768` at orders one through six. Converse is tested by reversing
all arcs and canonicalizing again.

The complete strong inventory with H at most 21 is:

| Order | H | Ordinary strong classes | Self-converse strong classes |
|---:|---:|---:|---:|
| 3 | 3 | 1 | 1 |
| 4 | 5 | 1 | 1 |
| 5 | 9 | 2 | 2 |
| 5 | 11 | 1 | 1 |
| 5 | 13 | 1 | 1 |
| 5 | 15 | 2 | 2 |
| 6 | 15 | 2 | 0 |
| 6 | 17 | 2 | 2 |
| 6 | 19 | 2 | 0 |

Arthur H. Busch's published bound gives `H(C)^3 >= 5^(|C|-1)` for every
nontrivial strong tournament. Since `5^6 > 21^3`, no strong tournament of
order at least seven has H at most 21. Thus the finite inventory above is
complete **at all orders in this bounded count range**, not by extrapolation.
[Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*,
EJC 13 (2006), N3, Theorem 1 and Corollary 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf).

Strong components have a unique dominance order and every Hamiltonian path
concatenates paths in that order, so H is their product. The table misses both
7 and 21; a product equal to 7 or 21 needs a factor 7 or 21. This replays the
inherited all-order exclusions. It does **not** show every other odd integer
is attainable.

## Inheritance and boundaries

* `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md` already excludes
  both values at all orders and records the 2016 source/priority correction.
* `05-knowledge/results/rotation_owners_spectrum_20261020e/spectrum.md`
  already gives this small strong spectrum, the general singleton-SCC kernel
  bound, and the open completeness boundary.
* `05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md`
  already derives fixed-H ordinary, self-converse and converse-merged rational
  series and their eventual polynomial/parity-quasipolynomial behavior.
* Session 124's cycle-packing packet and Sessions 130/131's reflection/carry
  packets retain the tournament-realizability obligation. Congruence classes
  do not settle exact-value exclusions.
* OCF is attributed to Grinberg--Stanley, and appears as Corollary 20 in
  [Irving--Omar v3](https://arxiv.org/html/2412.10572v3), versus Corollary 19
  in v1. This packet uses literal directed cycles, not supports.

A hostile checked in the atlas is that H=9 does not determine the packing
polynomial: a strong order-five witness has `1+4x`, while the ordered sum of
two directed triangles has `1+2x+x^2`. Trailing zero coefficients are removed
before comparison. Self-converse flags also cannot be inferred from H alone: the two
strong order-five H15 classes are self-converse; the two strong order-six
H15 classes form a converse pair.

## Reproduction

```
python3 build_atlas.py > atlas.new.json
python3 -O build_atlas.py > atlas.new.optimized.json
cmp atlas.new.json atlas.new.optimized.json
cmp atlas.json atlas.new.json
```

All validity checks are explicit exceptions; no correctness test disappears
under Python optimization. No nonstandard dependencies or downloaded data.
