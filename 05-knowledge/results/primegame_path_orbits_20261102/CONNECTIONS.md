# The numerical prompts, with their actual maps

**CITED inheritance / VERIFIED small calculations / PROVED elementary
incidence consequences.** No identification is inferred from equal counts.

## 1. From the n5 metagraph to480 and288

The owner's-gauge table in
`groups_of_five_primegame_merged_metagraph_20261031.md`, section5, has
12 ordinary classes, eight self-converse classes and ten converse-merged
nodes. Of the eight self-converse classes:

| Kind | H values |
|---|---|
| Three pure-blue classes | 1,3,15 |
| Five mixed classes | 9,9,11,13,15 |

Two additional nodes each represent a non-self-converse pair. Thus3+5=8
does not count all merged nodes, and the three pure-blue classes are not
three converse-pair mergers.

The two H15 cores give the actual bridge to the trap computation:

| Five-vertex core | Aut | Marked-path orbits H/Aut | Six-vertex source/sink traps |
|---|---:|---:|---:|
| Mixed, scores1,2,2,2,3 | 3 | 5 | 2*6!/3=480 |
| Pure-blue regular, scores2,2,2,2,2 | 5 | 3 | 2*6!/5=288 |

**Map:** adjoin a universal source or a universal sink. H and the core's
automorphism group are preserved; deleting that unique added extremal
vertex recovers the core. The two versions are distinct labelled orbits.
An automorphism fixing a Hamilton word fixes every vertex, so path action is
free and H/Aut counts its orbits. Hence the same ratio480:288=5:3 is also
the ratio of the two core path-orbit counts. This is a structural explanation,
not an unexplained reuse of the numbers3 and5 from the class census.

The inherited `prime_inverse_ports_20261102/MOVES.md` and
`two_edge_transport_20261102/ESCAPES_AND_CARRIERS.md` establish respectively
two-edge versus three-edge minimum repair distance to H17, and the sharp
common three-step energy barrier. We do not claim to rediscover those results.
Their equivariant-selector obstruction remains:3 and5 do not divide17.
Our new cores A,B instead have trivial groups: their hidden bits show that
attachment memory can survive even after all automorphism ambiguity is gone.

## 2. A literal19=12+7 inside PRIMEGAME

For the fourteen fractions in [Conway1987, Theorem1](https://www.cs.unc.edu/~stotts/COMP210-s23/madMath/Conway87.pdf),
let N be the set of distinct numerators and D the distinct denominators.
This session's literal-entry count gives

```text
N = {1,11,13,15,17,19,23,29,55,77,78,95},                    |N|=12
D\N = {2,7,33,38,51,85,91},                                |D\N|=7
N∩D = {1,11,13,17,19,23,29},                               |N∩D|=7
N\D = {15,55,77,78,95},                                    |N\D|=5.
```

Thus |N∪D|=19=12+7=5+7+7, while |D|=14. There are28 occurrences of entries,
19 distinct entry values,14 instructions, and10 supported prime coordinates.
These are different objects. The older note already lists the nineteen
entries; the numerator/denominator partition is the present explicit reading.
The bounded recovery also found the earlier equality19=7+12 in THM-2422 /
`operation-fibres-summand-closure-and-twin-center-ancestry`, equation22,
as a distinct-summand ancestry witness. No operation map identifies that
summand rule with tournament H:12 is even and7 is forbidden as a full H.

Draw one directed entry edge d→n per instruction n/d. The multigraph has
19 vertices,14 edges and six components of sizes2,2,2,3,5,5. Its undirected
cycle rank is14-19+6=1, supplied by the reciprocal pair11/13,13/11; deleting
one of those parallel underlying edges leaves a forest. This graph records
fraction entries, **not** actual execution on positive integers.

A precise comparison is the integer matrix diagram

```text
Z^14 --B--> Z^19 --V--> Z^10,
B(e_instruction)=e_numerator-e_denominator,
V(e_a)=(nu_p(a))_p,
Delta=VB.
```

B has rank13. V has rank10: the eight prime entries provide eight basis
directions, while55 and15 supply5 and3 after subtraction. The inherited
Session129 kernel calculation proves rank(Delta)=10, kernel rank4. So the
valuation map creates three additional formal circulation directions beyond
the single entry-graph cycle. This explains exactly how factorization can
create arithmetic relations invisible in the ordinary graph. It does not
make all formal circulations executable.

The decisive hostile is the **inherited, fully guarded** cycle

```text
4301→17765→16302→9867→8671→23023→4301,
instruction word GBDEFA.
```

Each step is the first enabled fraction, not a jump between bare entry
vertices. The standard-library verifier replays all six decisions. The
complete periodic-family classification belongs to
`session129_primegame_tournament_partition_20261031/partition/README.md`,
not this packet. It retains simultaneous tokens, nonnegative exponents and
priority. The entry graph lacks all three.

## 3. What the pair plus three quartets really preserves

Recover the user's exact partition, with letters in Conway's original order:

```text
X0=JK={11/13,13/11}
X1=AIMN={17/91,1/17,1/7,55/1}
X2=CEFH={19/51,29/33,77/29,77/19}
X3=BDGL={78/85,23/38,95/23,15/2}.
```

Session129 already proves the unique primitive additive character under the
stated block conditions. In prime order(2,3,5,7,11,13,17,19,23,29) it is

```text
Phi=(-4,0,-4,2,2,2,2,3,-1,3) dot valuations.
Delta Phi is0,-2,+1,0 on X0,X1,X2,X3 respectively.
```

Every executed return word therefore obeys #X2=2#X1. Conversely that scalar
law alone does not imply execution. The source, residual data counters and
first-enabled guard must accompany the height. THM-4854 /
`primegame-four-blocks-register-decomposition-bipartite-flags-triangular-semiprimes-metagraph-quotient-of-complement-graph`
also records why the three quartets' shared-flag triangle is not the
order-four metagraph, which is a star with a loop.

**Transferable mechanism:** both PRIMEGAME and actual tournament re-marking
have short command alphabets with conjunctive source guards. A coarse height
or current guard signature may preserve one observation while losing future
tests. In our tournament family, actual reversible words implement all
Boolean conjunctions and invert the full response law. No simulation of
PRIMEGAME by these tournaments, or full-spectrum theorem, follows.

## 4. The proposed radicals require a different operator

The first proposed multiset
0, +/-i sqrt(2), (-1+/-sqrt(3))/2 has characteristic polynomial

```text
t^5+t^4+(3/2)t^3+2t^2-t.
```

It has trace-1 and a noninteger coefficient. More strongly, the last two
roots satisfy the irreducible polynomial2t^2+2t-1, and are not algebraic
integers. They cannot be eigenvalues of **any integer matrix**, even after
adding other eigenvalues. The other multiset
2,+/-i,(-1+/-sqrt(5))/2 has polynomial

```text
t^5-t^4-2t^3+t^2-3t+2
```

and trace1. Neither is a loopless tournament adjacency spectrum, whose
matrix has zero diagonal and integer entries. An independently defined
rational operator might have the first spectrum, but none is supplied here.

The corrected positive cyclic-operator bridge is inherited from
`session143_magic_square_cyclic_bridges_20261103/tournaments/RESEARCH.md`:
the regular T5 adjacency is P+P^2 for a five-cycle permutation matrix P.
Its characteristic polynomial is t^5-5t^2-5t-2. The representation-theoretic
bridge there retains this actual operator and its multiplicities. We do not
replace it with the supplied radical lists.

## 5. What this contributes toward the major target

The spectrum statement is [Chan–Pak, Conjecture6.2](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf),
credited there to bof and Gordon Royle. Our response inverse is a constructive
theorem about a particular operation language, not a solution of that
existence conjecture. It identifies a concrete missing coordinate: the
ordered joint source word, observable after actual complement and relocation.

The next useful positive target is to combine a source-entry construction
for the prime and7p lanes with a declared operation family whose preserved
joint law is proved sufficient. The existing four-coordinate arithmetic
compiler remains useful inside its language. The new obstruction forbids
extending it merely by identifying sources with equal shallow guards or
equal total H. Integer invertibility and realization are separate tests.
