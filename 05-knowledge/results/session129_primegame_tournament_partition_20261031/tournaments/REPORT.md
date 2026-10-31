# Session129 — a genuine three-coordinate model of the ten merged nodes

**Status:** independently reconstructed FINITE-EXACT census through order 6; inherited incidence mechanism from **THM-796 §6, `01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md`**. The arithmetic encoding below is a chosen, proved-injective dictionary, not a proof that PRIMEGAME generates the metagraph. No priority claim for the incidence result.

## Results first

1. The actual merged-node census is **2, 3, 10, 34** at orders 3, 4, 5, 6. Order 3 has **one complement line**, not one node. This repairs the user's first proposed identification without discarding the reciprocal pair.
2. The three order-4 nodes have Hamiltonian-path counts **1, 3, 5**. Endpoint deletion from order 5 gives a three-coordinate integer vector. Its primitive form distinguishes **all ten** order-5 nodes.
3. The three pure-blue order-5 nodes have precisely the three coordinate-axis primitive vectors. Thus there is an intrinsic correspondence between the three order-4 nodes and the three order-5 symmetry points.
4. Choosing multiplicatively independent labels **7, 17, 55** turns these vectors into ten distinct integers. The three symmetry points then encode as exactly **7, 17, 55**. This is a faithful arithmetic encoding, **not a tiling-mask interpretation** and not an F=U+S classification.

## 1. Exact carriers and the order-3 correction

Fix the path `n -> n-1 -> ... -> 1`. Its free nonpath arcs form the cube

`X_n = {0,1}^{binom(n-1,2)}`.

Bits are ordered by tiles `(x,y)` with `y=1,...,n-2`, then `x=n,...,y+2`. Bit zero directs `x -> y`; bit one reverses it. Define:

- `a(t)=t+all_ones`: complement **only the free arcs**;
- `r(t)(x,y)=t(n-y+1,n-x+1)`: grid reflection, realizing full tournament converse after reversing vertex names;
- `M_n`: tournament isomorphism classes modulo full converse;
- `L_n=X_n/<a>`: literal complement lines, including loops and parallel lines after projecting endpoints to `M_n`.

These are different quotients. In particular, full converse does not identify a transitive triangle with a directed 3-cycle. Both are self-converse and their score multisets differ. At order 3, the two masks `0,1` produce these two nodes; complement pairs them into a single blue line.

The fraction pair `{11/13,13/11}` is likewise a two-element set with a fixed-point-free reciprocal involution. Choosing either endpoint assignment gives an involution-equivariant bijection with `X_3`, and identifies its reciprocal pair with the single element of `L_3`. It does **not** identify it with a single merged node. The endpoint assignment itself is a choice, not forced by factorization.

## 2. The complete order-4 incidence structure

Use the canonical coordinate order determined by distinct Hamiltonian-path counts:

| Coordinate | Merged node | H | Tiling fibre | Deletion profile over `(TT3,C3)` |
|---|---|---:|---|---|
| A | TT4 | 1 | `{0}` | `(1,0)` |
| B | `1 -> C3` merged with `C3 -> 1` | 3 | `{2,4}` | `(1,1)` |
| C | unique strong order-4 tournament | 5 | `{1,3,5,6,7}` | `(2,3)` |

There are four literal complement lines: one blue `A--C`, one blue loop at C, and two black `B--C` lines. In particular, the three node fibres have sizes **1,2,5**, not 4,4,4. Thus the three four-instruction blocks can be assigned as *names of nodes*, but cannot be bijected instruction-by-instruction to those node fibres. That proposed stronger map fails by cardinality before any arithmetic is used.

## 3. Three coordinates really do distinguish the ten nodes

For a merged node `u` at order n, let `F_n(u)` be its full fixed-path tiling fibre. Delete the bottom path endpoint 1 and relabel, obtaining `d:X_n -> X_(n-1)`. In the chosen bit order, `d(t)=t >> (n-2)`. Define

`R_n(u,v) = #{t in F_n(u): d(t) in F_(n-1)(v)}`.

Deleting the other endpoint gives exactly the same count vector: reflection bijects the two counts and is invisible after converse merging. Divide the vector by the gcd of its entries to obtain `Rbar_n(u)`.

The following table is independently computed, with owner display-node labels and owner ordinary-class labels (lexicographic adjacency canonization). Except for the last row, raw and primitive vectors agree.

| Owner node | Ordinary classes merged | H | Type | Raw `(A,B,C)` | Primitive `(a,b,c)` | Chosen code |
|---:|---|---:|---|---|---|---|
| 0 | #0 | 1 | pure-blue | (1,0,0) | (1,0,0) | 7 |
| 1 | #1,#7 | 3 | pure-black | (1,1,0) | (1,1,0) | `7*17` |
| 2 | #2,#5 | 5 | pure-black | (2,3,5) | (2,3,5) | `7^2*17^3*55^5` |
| 3 | #3 | 3 | pure-blue | (0,1,0) | (0,1,0) | 17 |
| 4 | #4 | 9 | mixed | (2,0,7) | (2,0,7) | `7^2*55^7` |
| 5 | #6 | 9 | mixed | (1,4,4) | (1,4,4) | `7*17^4*55^4` |
| 6 | #8 | 11 | mixed | (1,0,10) | (1,0,10) | `7*55^10` |
| 7 | #9 | 15 | mixed | (0,2,3) | (0,2,3) | `17^2*55^3` |
| 8 | #10 | 13 | mixed | (0,5,8) | (0,5,8) | `17^5*55^8` |
| 9 | #11 | 15 | pure-blue | (0,0,3) | (0,0,1) | 55 |

There are 12 ordinary classes. Eight are self-converse; the other four form exactly the two displayed pairs. Burnside's count is `(12+8)/2=10`.

**Finite completeness proof.** Exhaustively enumerate every fixed-path mask, canonize its tournament over all vertex permutations, and merge the canonical forms of the tournament and its converse. Every tournament has a Hamiltonian path; independently, the checker also enumerates all labelled tournaments through order 5 and confirms the same class counts. The table contains all resulting merged nodes, and its ten primitive vectors are pairwise different. Hence equality of primitive vectors is equivalent to equality of merged nodes **at order 5**.

The same test gives 34 distinct primitive vectors for 34 order-6 nodes. THM-796 also reports 272 distinct vectors at order 7; that latter run is inherited here, not rerun. An all-order injectivity theorem is not asserted.

**Sharp order-6 boundary for the symmetry interpretation.** The axis criterion is not an all-order pure-blue characterization. In the checker's order-6 node order there are five primitive axis profiles: nodes 0 and 17 are pure-blue, but node 5 (H=5, fibre size 5, profile `e_2`) and nodes 32 and 33 (H=45, fibre size 15, profiles `e_6,e_8`) are mixed. Indices on these ten-dimensional basis vectors start at zero, matching the order-5 node table. Thus the useful order-5 correspondence between coordinate primes and symmetry points fails already at the next order, even though completeness of the whole primitive profile survives.

The smallest counterexample has a general proof. If S is self-converse, put `T=1 -> S -> 1`. Every Hamiltonian path of T begins at its unique source and ends at its unique sink. Deleting either endpoint yields `1 -> S` or `S -> 1`, which have the same merged class because S is self-converse. Hence T has an axis deletion profile. Positional reversal of a path is an anti-automorphism of T exactly when reversal of its middle path is an anti-automorphism of S; therefore T is pure-blue **iff** S is pure-blue. Taking the unique strong S4 gives the mixed order-6 axis node 5 above, with five tilings of which three are blue. No smaller counterexample exists, by the exhaustive order-3 through order-5 census.

## 4. Faithful arithmetic encoding, including the composite label 55

Define the selected code

`E(u)=7^a * 17^b * 55^c`, where `(a,b,c)=Rbar_5(u)`.

Its prime factorization is `7^a * 17^b * 5^c * 11^c`. Therefore

`a=nu_7(E), b=nu_17(E), c=nu_5(E)=nu_11(E)`.

Unique prime factorization proves injectivity of this encoding on **all** nonnegative triples, hence on the ten-node set. The coordinate labels need not all be prime: multiplicative independence is sufficient, and pairwise coprimeness of 7,17,55 establishes it. The pure-blue nodes TT5, `1 -> C3 -> 1`, and T5 have primitive axis profiles, so they encode as 7,17,55 respectively.

This recovers the requested triple in a legitimate arithmetic role while respecting earlier obstructions:

- It is **not** the claim that the integer bit masks 7,17,55 represent the pure-blue nodes. They do not, and the previous Hamming-distance obstruction remains valid.
- It is **not** the claim that 55 is F=U+S-balanced. It is not.
- The assignment `A->7, B->17, C->55` is selected. The tournament structure determines the three axes and their canonical H-order, not these numerical labels. Many other multiplicatively independent triples work equally well.
- The prime factorizations decode a ten-node invariant **after the incidence table has been constructed**. They do not independently derive the number ten from Conway's fourteen fractions.
- The three scalar labels 7,17,55 all occur in the user's *same* reset block. Accordingly this encoding does not, by itself, justify identifying the user's three four-instruction blocks with A,B,C.

### Using the actual three four-instruction blocks instead

The companion [block-product codec audit](../bridge_audit/README.md) proves that the three literal block products are multiplicatively independent and admit an integral valuation decoder. Therefore the same ten vectors may be encoded by their products `q1^a q2^b q3^c`, using the user's actual factorization rather than the selected labels 7,17,55. This is a second faithful dictionary; its assignment of the three blocks to the canonically H-ordered coordinates is chosen. The count ten still comes from the independently constructed incidence profiles, not from unique factorization alone. See that audit for the full codec proof, exact formulas and limits.

## 5. The exact information-loss boundary

A complete static fingerprint need not be a sufficient recursive state. The two order-4 masks **2 and 4** lie in the same merged node B. The pure-blue order-5 sandwich node has sole mask **16**. Since `16 >> 3 = 2`, lower masks 2 and 4 have respectively **one and zero** lifts into that node. Thus growth is not strongly lumpable over merged nodes. The missing coordinate is the marked Hamiltonian-path presentation (already the converse sheet matters in this smallest witness).

Consequently no generic Markov recurrence on bare node labels or their prime encodings follows. Keep the incidence span `M_n <- X_n -> M_(n-1)`, not an invented parent map. This is the inherited obstruction in THM-796 §§6,8, now supplied with a literal minimal-size witness.

## 6. Triangular and Fibonacci contacts, typed precisely

**Universal triangular structure.** Total arc count is `binom(n,2)`; fixing the path removes `n-1` forced arcs, leaving `binom(n-1,2)` free bits. At orders 3,4,5 these dimensions are **1,3,6**, not the merged-node counts **2,3,10**. At order 5 the total arc count ten happens also to equal the merged-node count; these are different carriers. The next order already gives 15 arcs versus 34 nodes.

**Exact labelled arithmetic.** `55=5*11=T_10=F_10`, `91=7*13=T_13`, `78=2*3*13=T_12`, and `15=3*5=T_5`. The earlier file `05-knowledge/results/triangular_fibonacci_heron_boat_kps.out` explicitly recorded the first identity and the Fibonacci ranks of the ten primes. Its bounded computed output is not used here as a proof that a finite list of Fibonacci-triangular intersections is globally complete.

**Exact incidence contacts.** The strong order-4 node has deletion counts `(2,3)`. Order-5 rows include `(2,3,5)`, `(0,2,3)`, and `(0,5,8)`, so consecutive Fibonacci windows really occur in the table. These are finite exact identities, **not a proved Fibonacci recurrence** for all nodes or all orders. The canonical transitive node's complement partner has an inherited **tribonacci**, not Fibonacci, Hamiltonian-path law (THM-4731 §6).

This gives two honest directions for exploration: Fibonacci arithmetic in Conway's chosen labels, and Fibonacci-sized finite deletion fibres in the metagraph. A map preserving a stated operation between them is still missing.

## 7. Reproduction and audit scope

Run from any directory, Python standard library only:

```sh
python3 check_tournaments.py > census.json
python3 -O check_tournaments.py > census_optimized.json
cmp census.json census_optimized.json
```

The checker imports no inherited producer or data file. It enumerates fixed-path cubes through order 6; canonizes over all vertex permutations; checks converse Burnside counts, Hamiltonian-fibre counts, left/right deletion equality, black-line parity, the injectivity of primitive signatures through order 6, literal non-lumpability, and prime-valuation decoding. An independent full-labelled enumeration through order 5 checks that fixed-path canonization missed no class. Failure gates are explicit exceptions, active under `-O`. To avoid a large trace, the order-6 output retains every node's full primitive signature and invariants but omits individual masks and literal lines; they are recomputed internally.
