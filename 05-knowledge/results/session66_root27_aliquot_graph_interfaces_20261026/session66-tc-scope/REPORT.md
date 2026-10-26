# Session 66 appendix — ordinary Hamiltonian counts versus presentation-orbit counts

**Status:** elementary proof restatements and exact codec consequences; independent finite checks only through order five. No novelty claim. The incoming order-ten output is **REPORTED / NOT INDEPENDENTLY REPLAYED OR CERTIFIED HERE**. The all-order presentation-count exclusion remains **OPEN** in the pinned source.

## 1. The scope correction

The sealed `session66-graph-pairs/REPORT.md` addresses the classical gaps
`H(T) not in {7,21}`, where `H` counts directed Hamiltonian paths as vertex
words. It does **not** settle the broader prior question about

```
tc(T) = H(T)/|Aut(T)|,
```

the number of automorphism orbits of Hamiltonian paths, also called
presentation counts in the repository. Read its sentence “{7,21} are forbidden
Hamiltonian-path counts” as describing **that theorem**, not denying the
existence of the separate orbit-count question.

The exact lineage is:

- **LEM-003**, `01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md`:
  the elementary free-action lemma below is inherited, not newly discovered.
- **THM-1370**, `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`:
  the all-order ordinary-`H` exclusion, recovered in the sealed graph report.
- **THM-4741**,
  `01-canon/theorems/THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md`:
  separately reports finite `tc` exclusions through order nine.
- **HYP-9352(2)**,
  `05-knowledge/hypotheses/HYP-9352-presentation-count-laws-strong-spectrum-forbidden-fibre-sizes-and-the-three-move-law-at-n-11.md`:
  the all-order `tc not in {7,21}` claim is an open forbidden-fibre question,
  not a consequence of THM-1370 for non-rigid tournaments.
- Incoming **T-P1**,
  `05-knowledge/results/s57_untouchables_root27_20261026/PREREGISTRATION.md`:
  predicts no **strong order-ten tournament with a nontrivial automorphism**
  has `tc` equal to seven or twenty-one.

LEM-003, THM-4741, HYP-9352 and T-P1, together with the incoming
script/output, are frozen at commit
`67e6ca529b9ca594a7ddc042abe9800a6fdb0276` in `SOURCES.json`.
The sealed graph report and its older input-pin manifest are copied without
modification; their inherited source pin is
`63f52936f8dd8592d4f2fee869748baad2e1997e`.
This bounded clarification neither modifies that report nor promotes any
hypothesis or incoming computation to canon.

## 2. Why the quotient is an integer (PROVED, inherited elementary lemma)

Let `P(T)` be the set of vertex words `(v1,...,vn)` that use every vertex once
and satisfy `vi -> v(i+1)`. An automorphism `g` acts coordinatewise. It
preserves membership in `P(T)` because it preserves arcs. If it fixes one
word, then `g(vi)=vi` for every position; because the word contains all
vertices, `g` is the identity. Thus every stabilizer is trivial and every
orbit has exactly `|Aut(T)|` members. Partitioning `P(T)` into its orbits gives

```
H(T) = |Aut(T)| * |P(T)/Aut(T)|,
tc(T) = |P(T)/Aut(T)| in Z_{≥0}.
```

For nonempty tournaments a Hamiltonian path exists, so the quotient is
positive. With the empty-word convention the empty tournament has
`H=|Aut|=tc=1` as well. The free-action proof actually works for any finite
digraph. It does **not** work for directed Hamiltonian cycles modulo rotation:
rotations can stabilize a cycle.

For a rigid tournament `|Aut(T)|=1`, so `tc=H` and the ordinary-`H` theorem
settles both exclusions. For a non-rigid tournament, `tc=7` would instead
mean `H=7|Aut(T)|`, and similarly for twenty-one. Excluding the two values
`H=7,21` does not exclude these larger multiples. This is a logical scope
boundary, not a claimed counterexample to HYP-9352(2).

## 3. What the faithful graph codecs preserve (PROVED)

### Joint seven-vertex-gadget pair

For each tournament vertex use the tree with a degree-three center and arms
of lengths `1,2,3`. Its tail port is the length-one endpoint and its head port
the length-two endpoint. The first graph is the disjoint union of these
asymmetric trees; the second joins tail `u` to head `v` exactly when `u -> v`.

A **joint** pair automorphism must preserve the first-layer components and
each intrinsically distinct position within them. It therefore induces a
permutation of the tournament vertices with no residual gadget choices. The
second layer says precisely that this permutation preserves every arc.
Conversely every tournament automorphism has a unique such lift. Hence

```
Aut(G1(T),G2(T)) ≅ Aut(T),
tc(T) = H(decoded(G1,G2))/|Aut(G1,G2)|.
```

The reconstruction also bijects Hamiltonian words and intertwines the two
group actions. This transports the **orbit-count invariant itself**, not just
integrality. It is about the common-bijection joint pair; forgetting the
cross-layer matching does not retain this theorem.

### Independent graph-class pair and the leaf kernel

In the graph `Γ(T)`, each original vertex has center, tail, and head cores;
the center is joined to its two ports. Attach respectively `1,2,3` leaves to
these cores, and encode each arc as a tail-to-head edge. Leaves have degree
one and cores degree at least three. The number of leaf neighbors intrinsically
recognizes the three core types; a center recognizes its own two ports.
Thus the graph reconstructs the tournament and induces a surjection

```
Aut(Γ(T)) -> Aut(T)
```

via its action on recovered centers. In the kernel every center and hence
every port is fixed. Its only freedom is to permute, independently, the two
tail leaves and three head leaves in each gadget. All such permutations
really are automorphisms, so the kernel is exactly `(S2 x S3)^n`, of order
`12^n`. Every tournament automorphism lifts, by matching local leaf slots;
therefore the sequence splits (after choosing those slots), and

```
Aut(Γ(T)) ≅ (S2 x S3)^n ⋊ Aut(T),
|Aut(Γ(T))| = 12^n |Aut(T)|.
```

Let `Ω(T)` be the directed-odd-cycle conflict graph from the sealed report,
using literal directed cycles modulo rotation, not only their supports. The
inherited odd-cycle formula says `H(T)=I(Ω(T),2)`. Consequently the faithful
independent pair `([Γ(T)],[Ω(T)])` determines

```
tc(T) = I(Ω(T),2) / |image(Aut(Γ(T)) acting on centers)|
      = 12^n I(Ω(T),2) / |Aut(Γ(T))|.
```

The naive quotient `I(Ω(T),2)/|Aut(Γ(T))|` is smaller by `12^n` and is
wrong for every nonempty tournament. For `K1` it is `1/12` instead of one;
for the directed triangle it is `3/5184=1/1728` instead of one. The graph
kernel acts trivially on decoded path words, so the **full graph action**
on those words is not free. Passing to its induced center action is the
necessary repair. None of these formulas concerns the unrelated number of
Hamiltonian paths in the larger gadget graph itself. Also, `H=I(Ω,2)` is
used only as a cardinality identity here: no equivariant bijection between
weighted independent sets and Hamiltonian paths is asserted, and no
division by `|Aut(Ω)|` is licensed.

This is why fidelity on isomorphism classes alone is insufficient to divide
by a graph's full automorphism order. A faithful codec makes the orbit-count
question expressible in graph language; it does not prove its exclusions.
Nor is the informative pair claimed polynomial-size: the cycle-conflict
coordinate can be large, as already scoped in the graph report.

## 4. Incoming finite claim and exact verification boundary

The pinned `04-computation/untouchables_s57/tc_fibres_n10.json` reports empty
hit lists in six symmetry families:

| permutation type | invariant assignments | strong assignments | candidate checks | reported hits |
|---|---:|---:|---:|---:|
| `(3)(3)+4` | 524288 | 391104 | 101568 | 0 |
| `(3)(3)(3)+1` | 32768 | 23520 | 8784 | 0 |
| `(5)+5` | 131072 | 89280 | 28800 | 0 |
| `(5)(5)` | 512 | 480 | 40 | 0 |
| `(7)+3` | 512 | 192 | 240 | 0 |
| `(9)+1` | 32 | 0 | 0 | 0 |

These are the stored output's counters, **not** counts of disjoint
isomorphism classes. The source script can perform separate candidate
checks for both target values, explaining why a candidate counter can exceed
a strong-assignment counter. Its reduction and coverage argument are only
read and pinned here, not independently certified. Neither its order-ten
execution nor the inherited through-nine census was replayed in this appendix.
No all-order conclusion follows from either finite bound.

Run the independent, standard-library appendix controls with:

```
python3 check.py > normal.json
python3 -O check.py > optimized.json
cmp normal.json optimized.json
```

Both outputs are identical and all gates are explicit exceptions, not
optimization-removable assertions. The controls cover:

- all **1,100 labelled tournaments on orders 0 through 5**: enumerate all
  Hamiltonian words and automorphisms, explicitly partition into **6,924 path
  orbits**, check free orbit sizes and the integer quotient;
- both codec decodings on every such tournament; independently evaluate the
  cycle-conflict independence polynomial and compare with the path count;
- all **21 isomorphism classes** on those orders: count both graph
  automorphism groups with a generic edge-colored graph backtracker after
  false-twin compression (no decoder or tournament automorphism routine is
  used by that graph counter); check exact group-order formulas and repeat
  decoding after graph-label reversal;
- the `C3` wrong-denominator hostile, empty boundary, small `tc` spectra, and
  SHA-256/Git-blob integrity of the copied sources. Source integrity is not a
  certificate of the incoming computation's mathematical correctness.

The all-order statements in Sections 2–3 rest on the supplied elementary
proofs, not extrapolation from the finite tests. The tests are bounded
independent implementation checks, not a new large census, global forbidden-
fibre theorem, or classification of graph-code images.
