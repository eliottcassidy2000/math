# Thirty Fano systems under the Paley affine group

**Status: PROVED IN NOTES after independent audit; FINITE-EXACT replay.** The proposed orbit profile `1+1+7+21` is **REFUTED**.
The exact profile is

```
30 = 1 + 1 + 7 + 7 + 7 + 7.
```

This is a new comparison of the already recovered objects, not a claim that
the Fano plane, its affine group, or Burnside's lemma is new. No tournament
Hamiltonian count or Hadwiger–Nelson bound is inferred.

## Objects, action, and inheritance

Let V=F7, and let F be the set of all labelled Fano line systems on V: each
has seven three-element lines and exactly one line through each pair. The
group is `G={x↦ax+b : a∈{1,2,4}, b∈F7}`, of order21, acting by point
permutation. This is the same labelled G-action as the recovered
[Paley-blocker augmentation bridge](../../session54_reframed_numerology_20261025/sixtythree-bridge/README.md).
The source here has30 line systems, not63 blockers or35 individual triples.

The [HN palette/deck packet](../../session61_hadwiger_six_seven_20261025/palette-deck/RESULT.md)
supplies the precise matching-to-Fano completion map. Fixing a point0, its
three incident lines partition the other six points into pairs. There are15
pairings; each has exactly two Fano completions, by the parity of transversal
endpoint choices. Therefore |F|=30. These are also the30 labelled binary
simplex codes recovered there. The marking, not the number30 alone, identifies
the objects.

In the field frame `F8=F2[alpha]/(alpha³+alpha+1)`, the exponent triple
`D={1,2,4}` has `alpha+alpha²+alpha⁴=0`. Define

```
F_+ = {b+D : b∈F7},
F_- = {b−D : b∈F7},       −D={3,5,6}.
```

These are two distinct Fano systems, fixed individually by G. They are the
two Fourier-axis geometries in [RESULT.md, B3](RESULT.md#b3-the-f8f4-comparison-exists-only-after-restricting-symmetry).
Passing from an exponent triple to the corresponding additive F8 line keeps
the selected alpha frame; it does not make all Fano systems simultaneous
Euclidean unit-distance realizations.

## Structural orbit calculation

**Order-seven elements fix exactly two systems.** A nontrivial translation
generates all translations. No three-element subset is translation-fixed,
so an invariant Fano system must consist of the seven translates of one
triple B. This is a Fano system iff its six ordered nonzero differences are
exactly1,…,6, each once. There are five translation orbits of triples:

| Base triple | Cyclic gaps | Difference-set condition |
|---|---|---|
|012|1,1,5|fails|
|013|1,2,4|holds: F_+|
|023|2,1,4|holds: F_-|
|014|1,3,3|fails|
|024|2,2,3|fails|

Completeness can be seen from the15 positive ordered compositions of7 into
three parts, divided into five cyclic classes. Each valid difference system
is fixed by doubling, so these are also exactly the two G-fixed systems.

**Order-three elements fix exactly six systems.** All14 nontranslation
elements have order three and are conjugate to a power of `sigma:x↦2x`.
It is enough to count fixed systems of sigma, whose point orbits are
`{0}`, `{1,2,4}`, `{3,5,6}`. The pencil through0 of an invariant Fano system
must be an invariant perfect matching on the other six points. Every pair
joins the two length-three cycles: an orbit of pairs lying within one cycle
would give degree two, not a matching. An invariant cross-cycle matching is
an equivariant bijection between the cycles, giving exactly three choices.
Each matching has two Fano completions. An order-three action on that two-set
must be trivial, so both completions are fixed. The count is therefore6.

Burnside now gives

```
#orbits = (30 + 6·2 + 14·6)/21 = 6.
```

An orbit size divides21. A size-three orbit would have a stabilizer containing
the unique C7 subgroup, contradicting the fact that both C7-fixed systems
are already globally fixed. After removing the two fixed systems, four
orbits remain, of size at least seven, on28 objects. Hence each has size
seven; every nonfixed stabilizer is a C3 complement.

**Failure anatomy of the proposed `1+1+7+21`.** The scalar30 admits that sum,
but it would give four orbits, whereas the actual order-three fixed-point
count forces six. Replacing a count of whole Fano systems by a count of their
individual triples is the first invalid transfer. On the35 triples the same
group really does have orbit sizes `7+7+21`; it is a different G-set.

## A lossless intrinsic address for the remaining28

For a nonfixed system L, its C3 stabilizer fixes a unique point c(L) of F7.
Put

```
type(L) = (|L∩F_+|, |L∩F_-|).
```

The four possible types are `(0,1),(1,0),(1,3),(3,1)`. The following four
systems have centre0 and those respective types. A string123 denotes the
unordered triple {1,2,3}; the records retain all seven lines.

| Type | Lines of the centre-zero system |
|---|---|
|(0,1)|123,034,025,145,016,246,356|
|(1,0)|124,034,025,135,016,236,456|
|(1,3)|023,124,015,345,136,046,256|
|(3,1)|013,234,125,045,026,146,356|

Each row is directly a Fano system invariant under doubling, with the stated
intersections. These four and F_+,F_- exhaust the six systems fixed by the
centre-zero C3. Each displayed system has a seven-element orbit; different
types cannot share an orbit. The structural count above therefore makes the
map

```
L ↦ (type(L), c(L))
```

a bijection from the remaining28 to four copies of F7. It is G-equivariant:
the type is fixed because F_+ and F_- are fixed, and stabilizers conjugate,
so `c(gL)=g c(L)`. Explicit inversion translates the appropriate centre-zero
row by c. Dropping the type makes the centre map exactly four-to-one.

This supplies the required source, target, map and preserved predicate:
labelled Fano systems map to two fixed symbols plus four labelled affine
seven-sets; the full incidence and G-action are recovered from the retained
type/centre address. It does not identify their native operations with arc
XOR, field addition, or a tournament path-complement operation.

## The larger-group and incoming-object boundaries

The full automorphism group of F_+ has168 elements. The point permutation
`(2 4)(5 6)` fixes F_+ but sends F_- into the remaining28. Thus removing the
two G-fixed systems does **not** produce a subset invariant under that larger
168-element action. No identification with the recovered Coxeter28 or a
PSL(2,7)-equivariant conic construction follows merely from the number28.
This is a literal failed invariance test, not a claim that no differently
marked or differently defined map could ever exist.

The incoming
[Session321 FANO_METAGRAPH.md](../../session321_selective_returns_and_fano_wiring_20261113/FANO_METAGRAPH.md)
uses **128 subsets of the seven quotient points** as independent internal
triangle-flip masks of a specific21-vertex tournament. Its12 ordinary and6
merged classes are a different universe and quotient. That packet's three
types of size-three masks are the `7+7+21` triple orbits; they are not the
orbits of the30 complete Fano line systems here. The common exponent labels
and G21 action give a legitimate comparison while retaining both source
objects. The six orbit count here is not a second derivation of its six
merged tournament classes or their Hamiltonian spectra.

## Exact replay

[check_fano_orbits.py](check_fano_orbits.py) uses two independent generators:
all relabellings of the seven nonzero vectors of F2³, and an exact-cover search
partitioning the21 unordered pairs into three-pair triples. Both give the
same30 systems. It checks every affine group composition and every action
on all30 systems, stabilizers, fixed-point counts, all35 triples, invariant
pencils, the full four-copy address and the larger-group hostile.

```
python3 -B 05-knowledge/results/integer32_twoedge_20261113/bridges/check_fano_orbits.py
python3 -B -O 05-knowledge/results/integer32_twoedge_20261113/bridges/check_fano_orbits.py
```

Both saved receipts, [normal](fano.normal.json) and
[optimized](fano.optimized.json), agree at **1,798 always-active gates**.
They list representatives, all stabilizers, all28 addresses, and all21
fixed-point counts. Checker SHA256:
`8effd5070d1b6d79718a205d969fa22b3ae54ea408d8d1c3d8d5606da6a504cd`;
each receipt SHA256:
`2d0525d93012115e5ce251e76d9e73f8702e29b219f9943eeab47f7ecd4518e7`.
The finite object is completely enumerated; no all-order or external
geometric realization claim is attached to that completeness.
