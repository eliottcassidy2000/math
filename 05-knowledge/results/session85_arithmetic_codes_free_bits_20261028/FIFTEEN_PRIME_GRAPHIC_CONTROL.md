# Two different fifteens: a graphic supersingular bank and the89 Fano obstruction

Date:2026-10-28. **PROVED finite square-class representation; FINITE-EXACT checks.** The supersingular-prime list is CITED classical input. No claim identifies Monster-group structure or modular-curve geometry with the graph constructed below.

## 1. Use the actual old fifteen as a control

The classical supersingular/Ogg prime set is

\[
P=\{2,3,5,7,11,13,17,19,23,29,31,41,47,59,71\}.
\]

This list is stated in Gregory Sankaran's *A supersingular coincidence*, Ramanujan Journal59(2022),609–613; see the [author's institutional publication record](https://researchportal.bath.ac.uk/en/publications/a-supersingular-coincidence/) and [primary preprint](https://arxiv.org/abs/2009.11379). Only the list and its classical name are imported here.

For a selected subset A of these distinct primes, square sigma of the squarefree product is equivalent to zero sum of the parity vectors of p+1. In the valuation basis(2,3,5,7), the masks are

|p|2|3|5|7|11|13|17|19|23|29|31|41|47|59|71|
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|mask|2|0|3|1|2|9|1|4|3|7|1|11|2|6|1|

The matrix has rank4, so the square-sigma subset code has dimension11 and2048 words. These are the same length and dimension as Hamming[15,11,3], **but it is not that code**: selecting the prime3 alone gives a weight-one word, and equal columns give ten weight-two words. Equal dimensions do not identify dependency structures.

## 2. A literal graph-cycle representation on five vertices

Change the square-class basis to

\[
(3,2,15,14).
\]

These four classes are independent: their vectors in the original basis are(2,1,6,9). Give graph vertices0,1,2,3,4 the potential classes1,3,2,15,14. The following labelled multigraph has one edge per prime:

|Endpoints|Prime edge labels|Square class of each p+1|
|---|---|---|
|0–1|2,11,47|3|
|0–2|7,17,31,71|2|
|0–3|59|15|
|1–2|5,23|6|
|1–3|19|5|
|2–3|29|30|
|0–4|13|14|
|1–4|41|42|
|loop at0|3|1|

Each edge's square-class column is exactly the sum of the two endpoint potential vectors. Therefore its binary square-class matrix is row-equivalent to the graph's reduced incidence matrix.

**Theorem.** A subset of the fifteen primes has square sigma iff its labelled graph edges form an even subgraph. The arithmetic binary matroid is graphic, with a loop and parallel edges retained.

Indeed independence of the four nonroot potentials makes zero total column equivalent to even degree at every nonroot vertex. The remaining degree is even by the handshaking identity. The loop contributes zero parity, just as sigma(3)=4 contributes no square-class obstruction. The graph is connected with15 edges and5 vertices, so its cycle-space dimension is15-5+1=11.

The loop explains the weight-one codeword. Parallel bundles of sizes4,3,2 supply respectively6,3,1 weight-two circuits, explaining all ten. Discarding these labels or simplifying the multigraph would change the selection problem.

The five vertices here are a concrete incidence realization, not the five-dimensional affine-check space in the separate sixteen-prime extended-Hamming construction.

## 3. Why89 is the decisive extra prime in the Fano bank

Six primes of the user's Fano bank already lie in P:

\[
\{2,5,7,19,29,59\}.
\]

In the above multigraph, choosing just those six labelled edges gives a simple K4 on vertices0,1,2,3. Their square classes in the original three-coordinate basis(2,3,5) are all the nonzero vectors **except the class10**.

The missing class is supplied by89, because89+1=90 has squarefree part10. Moreover89 is the least prime of that class: the candidates are10k²-1, and k=1,2 give9,39 while k=3 gives89.

Adding89 completes the seven nonzero vectors of F2³, hence the Fano matroid. A simple graphic matroid of rank3 has at most six elements: a simple graph of total graphic rank3 has at most the six edges of K4 (splitting into components cannot increase that maximum). The seven-column Fano matroid is therefore not graphic.

This is an exact boundary rather than a general rejection of graph analogies:

- the old fifteen-prime square-class model **is graphic** on five vertices;
- six of the Fano primes form its K4 restriction;
- the additional prime89 fills the missing projective point and destroys graphicness.

In particular the sixteen-prime bank P union{89} has a Fano restriction and is nongraphic, although its square-class rank remains4 and its code dimension becomes12. This is **not** the separately constructed extended Hamming[16,11,4] bank, whose rank is5. The same number16 labels different objects with different constraints.

### The missing generator is a four-vertex hyperedge

In the five-vertex representation,89's new column is the incidence vector of the **four-vertex hyperedge{0,1,2,3}**, not an ordinary edge. This gives an explicit successful higher-order model, rather than merely a negative graphicness test. Its boundary is also the XOR of either of the following perfect matchings of the K4:

- edges01 and23, labelled2 and29;
- edges02 and13, labelled7 and19;
- edges03 and12, labelled59 and5.

Consequently three of the Fano circuits are exactly{2,29,89},{7,19,89},{5,59,89}. The other four are the ordinary triangles of K4. This describes all seven three-element Fano circuits without losing the four-point circuits of the full matroid.

The total signature of the original fifteen columns is precisely the missing class10. Thus their full product has nonsquare sigma and has no unused prime left for an addition-only repair within that bank. Adding89 makes the product square-sigma. The finite-bank bound from the complete Hamming bank cannot be exported to an arbitrary bank without its hypotheses.

The same three perfect matchings are the indices in the earlier four-charge pairing identity. For any four vectors q_i in a real inner-product space, let E_children=sum||q_i||² and E_parent=||sum q_i||². If E_pi sums the two squared pair sums in matching pi, then

\[
\frac13\sum_{\pi}E_\pi=(E_{\rm parent}+2E_{\rm children})/3.
\]

Each unordered pair occurs exactly once across the three matchings, which proves the identity by expansion. Their binary boundaries agree, but their energies need not: charges(1,1,-1,-1) give energies8,0,0. Thus this is an exact shared combinatorial indexing object, not a claim that parity retains the energy or that one pairing is universally optimal.

## 4. Different free-bit behavior at the same15/11 dimensions

The old fifteen-prime code has a dual parity check of weight2:

\[
x_{13}+x_{41}=0,
\]

because only13+1 and41+1 contain the prime7 to odd exponent. Thus either selection bit is recoverable from the other, assuming square sigma. There is no weight-one dual check, so individual selection bits remain unbiased under the uniform code law; the first nontrivial correlations occur at order2.

The explicit Hamming15 bank instead has dual minimum weight8 and uniform coset marginals through order7. The old and new fifteen-prime banks therefore have different local recovery and error-detection behavior despite equal length, rank and kernel size.

The complete dual weight enumerator of the old bank is

\[
1+z^2+z^3+z^5+3z^7+3z^8+3z^9+3z^{10}.
\]

Its primal code has minimum distance1, not3. This is a concrete hostile to inferring code structure from dimension alone.

The result provides a real connection to the previous supersingular-prime work through a specified observable, namely the square classes of p+1. It supplies no transport of supersingularity under arbitrary code symmetries and no claim about a canonical Monster action on these graph edges.

## Reproduction

Run `python3 check_fifteen_control.py`, also under `-O`. It checks the prime list, basis change, all32768 squarefree selections by exact sigma multiplication and integer square roots, the explicit even-subgraph equivalence, exact primal/dual enumerators, the13/41 recovery check, the K4-to-Fano hyperedge extension and the three-pairing identity. The graphic representation is proved by the displayed incidence identity, not inferred only from the census.
