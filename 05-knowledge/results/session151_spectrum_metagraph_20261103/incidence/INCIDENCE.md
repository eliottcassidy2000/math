# Colored metagraph incidence by path orbits and double cosets

## Scope and notation

Fix n≥3. Take one labelled representative T_C for each tournament isomorphism class C, with automorphism group G_C and order a_C. For each directed Hamiltonian path p=(v_1,...,v_n) in T_C, let A_p(T_C) reverse every edge except its consecutive path edges. This is the chord complement of THM-643, THM-4729, and the accompanying path-overlap bridge. Its marked path stays fixed and applying it twice returns T_C. Define the reflection r_p(v_i)=v_{n+1-i}.

Write ε_B(C,p)=1 if r_p is an anti-automorphism of T_C, and 0 otherwise, and ε_K=1−ε_B. This agrees with the usual blue condition S(t)=t in p-coordinates. Path reversal here is an anti-automorphism, not an automorphism.

## Path formula

For colors x=B,K define the directed incidence number

N^x_CD = (1/a_C) Σ_{p∈HP(T_C)} ε_x(C,p) · 1[A_p(T_C) ≅ T_D].

Then N^x_CD=N^x_DC. For C≠D, this is the number of individual color-x lines between C and D. For C=D, N^x_CC is twice the number of color-x loops, so the loop count is N^x_CC/2. The row sum is the color-x degree, with loops counted twice. In particular N^x_CC is even.

Proof. Automorphisms act freely on directed Hamiltonian paths, because preserving the sequence fixes every vertex; the orbits correspond bijectively to fixed-path tilings in fiber C. Chord complement bijects tilings in C ending at D with tilings in D ending at C, preserving color. A nonloop contributes one tiling to each endpoint fiber, while a loop has two distinct tilings in its fiber. The two are distinct since n≥3 gives at least one chord. Hence the formulas.

Equivalently, for blue the sum may be partitioned by anti-involutions α of T_C, summing paths with α(v_i)=v_{n+1-i}; each path has a unique such α. Since a tournament anti-automorphism fixes at most one vertex, only involutions with zero fixed vertices (even n) or exactly one (odd n) can occur. This potentially prunes the sum but is not a general complexity bound.

## Individual adjacency via agreement-path double cosets

A bijection φ:V(T_C)→V(T_D) is an agreement-path witness if the pairs on which T_C and the pullback φ* T_D agree are exactly the n−1 edges of one directed Hamiltonian path p in T_C. This path is unique: the agreement edges and their orientation determine its sequence. Equivalently φ is an isomorphism A_p(T_C)→T_D. Classify the witness blue when r_p is an anti-automorphism of T_C and black otherwise. Let W^x_CD be the set of witnesses of color x. Then

N^x_CD = |W^x_CD| / (a_C a_D).

Indeed G_D × G_C acts on W via (b,a)·φ=b φ a^(-1). The agreement path of the transformed witness is a(p). If a witness is fixed, a preserves the oriented spanning path and so is the identity; then b is also the identity. Hence this action is free, and its orbits are the marked paths in the first formula. In particular, after identifying vertex sets with {1,...,n}, N^x_CD equals the number of color-x double cosets G_D φ G_C contained in W^x_CD; each such double coset has exactly a_C a_D elements. This is a concrete filtered double-coset count, not a claim that all double cosets qualify. The inverse witness belongs to W^x_DC, proving symmetry also directly.

This counts individual line multiplicities, not merely whether adjacency exists. The formula requires knowing representative classes and automorphism groups. One can check at most n! bijections per ordered class pair, rather than listing 2^((n−1)(n−2)/2) marked tilings globally, but it is still factorial, has a potentially huge class input, and is not a polynomial-time algorithm or uniformly faster than global tiling enumeration. It is not a closed form for all entries or an enumeration-free classification of tournament classes.

## Converse-merged graph

Let E_U be the one or two original classes in merged node U. For every U,V define the folded directed incidence matrix

M^x_UV = Σ_{C∈E_U} Σ_{D∈E_V} N^x_CD.

It is symmetric. For U≠V it is the number of individual color-x edges joining U and V; for U=V, the number of loops is M^x_UU/2. The diagonal summation intentionally includes both orientations of a line between two distinct converse classes as well as twice each original loop. Thus no original cross-line is lost or counted twice as a merged loop. This fold and the previously established free reflection pairing immediately retain the even-black-bundle theorem; they do not claim a new theorem for it.

## Check and provenance

`python check.py` compares the agreement-path formula, for every pair of original classes and both colors, to direct fixed-path tiling enumeration for n=3,4,5. It checks symmetry, divisibility by the two automorphism orders, and the diagonal factor. Output: 2,4,12 classes and 1,4,32 lines respectively, all comparisons pass. The check is finite and the proof applies to all n≥3.

The original S,F definitions, fiber interpretation, color parity, and black reflection-pair theorem are recorded in THM-643 and THM-4731 (also THM-796); THM-4729 establishes the line convention. The Hamiltonian-path orbit weighting is already present in the repository's metagraph spectrum work. The path and witness formulas are deductions from those structures; no literature novelty is asserted. See [`STRUCTURE.md`](../metagraph/STRUCTURE.md) and [`OVERLAP_AND_PATH_WALK.md`](../bridges/OVERLAP_AND_PATH_WALK.md) for definitions, cautions and citations.
