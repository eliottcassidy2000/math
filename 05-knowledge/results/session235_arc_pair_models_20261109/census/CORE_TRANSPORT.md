# All-order transport of the semicomplete H7/H21 grammar

**Status: PROVED from Session234's complete core classification, with
FINITE-EXACT independent verification of every core.** All graphs in counting
formulas are loopless. No novelty claim is made for SCC decomposition or
ordinal-sum multiplication.

Let D be semicomplete and let E(D) delete **both** arcs of every reciprocal
pair, leaving each single arc unchanged. This is precisely decoding the same
raw `10/01/11` pair word with `11=empty`. It is a bijection from semicomplete
graphs to oriented graphs, with inverse replacing every missing pair by both
arcs. Both maps commute with every vertex permutation. Consequently they
preserve isomorphism distinctions and the entire automorphism group, even
though they need not preserve SCCs or Hamiltonian paths.

## 1. Why the finite core table proves an all-order statement

Use [Session234 DG5 and classification](../../session234_digraphs_20261109/PROOFS.md)
and its [complete generator](../../session234_digraphs_20261109/census/generator.py).
Every semicomplete H7 graph is an ordinal sum of one of three strong order-four
H7 cores and singleton gaps. Every semicomplete H21 graph is either one of 32
strong H21 cores (16 order five and 16 order six), padded with singleton gaps,
or one H7 core and one of two H3 cores in either order, with three singleton gaps.

All cross-SCC pairs in this source decomposition are single forward arcs, so
E leaves them unchanged. Thus

`E(C1 ⊕ ... ⊕ Ck) = E(C1) ⊕ ... ⊕ E(Ck)`.

Hamiltonian-path counts multiply under ordinal sum even if the individual
summands are not strongly connected: a path cannot return to an earlier
summand, so it must traverse each summand consecutively, and arbitrary
summand paths concatenate uniquely. This proves

`H(E(D)) = product_i H(E(Ci))`.

**Do not assert that decoded summands remain SCCs.** The source SCC signature
still classifies the decoded object because E has an isomorphism-equivariant
inverse. This is the sidecar needed to retain uniqueness and automorphism
weights.

## 2. Independently verified complete core table

| Source H | Order | Decoded H | Automorphism order | Classes | Sum 1/Aut |
|---:|---:|---:|---:|---:|---:|
|3|3|1|1|1|1|
|3|3|3|3|1|1/3|
|7|4|1|1|3|3|
|21|5|1|1|8|8|
|21|5|3|1|4|4|
|21|5|3|3|2|2/3|
|21|5|5|1|1|1|
|21|5|7|1|1|1|
|21|6|1|1|8|8|
|21|6|3|1|6|6|
|21|6|15|1|2|2|

The machine-readable [core_transport.json](core_transport.json) includes all
37 original matrices, decoded matrices, full path lists in both models, and
complete automorphism lists. [core_transport.py](core_transport.py) independently
enumerates every permutation for path counts and automorphisms, cross-checks
subset DP, and verifies strong connectivity of each source core. Source core
completeness is inherited from Session234, not inferred from this short script.

## 3. Exact all-order conclusions

- If H(D)=7, then H(E(D))=1, at every order.
- Among semicomplete graphs with H(D)=21, the complete set of decoded path
  counts over all orders is **{1,3,5,7,15}**. The first four occur at order five;
  15 first occurs at order six. Padding by singletons preserves them.
- Therefore the same raw array can never have H=7 under both semantics, nor
  H=21 under both semantics, at any order.
- Raw arrays with reciprocal-semantics H21 and empty-semantics H7 exist at
  every order n>=5, in exactly n-4 isomorphism classes and (n-4)n! labelled
  forms. This is a complete joint fibre, not a count of all oriented H7 graphs.

These conclusions concern fibres over the semicomplete targets 7 and 21.
They do not restrict the entire oriented spectrum. In particular the main
Session235 theory establishes its separate unrestricted all-order result.

## 4. Joint generating functions

Let U_(h,r)(x) count isomorphism classes of source semicomplete graphs D with
H(D)=h and H(E(D))=r by order. It equally counts their distinct decoded graphs.
Let L_(h,r)(x)=sum_n labelled_count(n)x^n/n! be the corresponding labelled EGF.
Put `B_k=x^k/(1-x)^2` and `C_7=x^7/(1-x)^3`.

The nonzero H7 fibre is

`U_(7,1)=L_(7,1)=3B_4`.

The H21 fibres are exactly:

| r | U_(21,r) | L_(21,r) |
|---:|---|---|
|1|8B_5+8B_6+6C_7|8B_5+8B_6+6C_7|
|3|6B_5+6B_6+6C_7|(14/3)B_5+6B_6+2C_7|
|5|B_5|B_5|
|7|B_5|B_5|
|15|2B_6|2B_6|

All other fibres vanish. Each single core has two independent singleton gaps.
For a mixed core word, there are three gaps and two possible core orders.
The three H7 cores all have decoded count one and trivial automorphism group;
pairing them with the H3 core of decoded count one gives coefficient six in
both OGFs/EGFs. Pairing with the directed triangle gives coefficient six in
the OGF and two in the EGF because its automorphism group has order three.

Equivalently, coefficients of B_k are max(n-k+1,0), and coefficients of C_7
are binomial(n-5,2) when n>=7, zero otherwise. Multiply EGF coefficients by n!
for labelled counts. All formulas are exact for every n, not merely eventual.
Summing rows recovers Session234's complete H21 class and labelled formulas.

## 5. Reproduction and controls

From the repository root:

```sh
python3 05-knowledge/results/session235_arc_pair_models_20261109/census/core_transport.py
```

The script writes core_transport.json. Besides full permutation checks on all
37 cores, it verifies the generating-function coefficients against every
matrix produced by the inherited grammar through order ten. For orders one
through six, its joint labelled counts also match the entirely independent
complete raw-array census in jointN.csv. The algebraic argument above, not
agreement through ten, supplies the all-order extension.
