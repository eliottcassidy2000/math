# The full order-eight electrical network and its dyadic filtration

October 23, 2026. **PROVED:** the general symmetric integral filtration
below. **FINITE-EXACT:** its application to the complete order-eight
tournament line metagraph. **OPEN:** a formula in the tournament order for
these groups, and an all-order characterization using small incidence data.

## What the graph and its weights mean

Use the converse-merged line multigraph, with a vertex for each converse
orbit of tournament classes. A blue line is fixed by the staircase
reflection. The remaining black lines come in reflected pairs with the
same merged endpoints. Write the literal graph as

    G = L + 2 B,

where L is the blue multigraph and B contains one edge per paired black
line. The input `black` column counts literal lines, hence must be even;
it is divided by two only to construct B. Parallel lines are conductances.
Self-loops remain in the export but contribute zero to the Laplacian.

The closest proved mechanism is the half-energy form in the incoming
[Session36 electrical packet](../session36_information_curl_contexts_20261023/metagraph/DYADIC-ELECTRICAL-METAGRAPH.md).
That packet computed the full graph through order seven and the blue graph
at order eight. This packet computes the **full** order-eight graph from the
independently reconstructed [tiling/port export](../metagraph_order8_20261023/README.md).
The earlier [dyadic counterfeit](../metagraph_bucket_obstruction_20261023/README.md)
is the hostile example: matching binary Laplacians and diagonal energies
need not preserve the bilinear half-energy. The least-used sidecar is a
literal integer lift of each surviving harmonic mode.

The invariant transferred here is the 2-primary part of the critical
group, equivalently the even Smith factors of a root-reduced Laplacian A.
It discards the odd factors and does not reconstruct the graph. Thus a
different group disproves equality of networks, but an equal group is not
an isomorphism or tournament-realizability certificate.

## An all-size lifting rule

Let A be any symmetric nonsingular integer r-by-r matrix. For j >= 1 put

    L_j = {x in Z^r : A x = 0 modulo 2^j},
    V_j = L_j modulo 2, as a subspace of F_2^r.

For u,v in V_j, choose any integer lifts x,y in L_j and define

    q_j(u,v) = (x^T A y)/2^j modulo 2.

**Filtration theorem.** This is a well-defined symmetric bilinear form and

    radical(q_j) = V_(j+1).

If s_1,...,s_r are the Smith factors of A, then

    dim V_j = #{i : 2^j divides s_i},
    rank q_j = #{i : v_2(s_i) = j}.

Consequently the entire 2-primary group can be recovered by successively
retaining the radical, lifting it one binary digit, and evaluating the next
divided energy. The process terminates for each nonsingular finite A.

**Proof.** Replacing x by x+2z changes the quotient by
`2 z^T(Ay/2^j)`, an even integer; symmetry gives the same result for y.
The form is therefore independent of lifts. If Ax=2^j b, the class of x
admits a lift to L_(j+1) exactly when there is z with

    A z = -2^(j-1) b modulo 2^j.

Over Z/(2^j), the image of a symmetric matrix is the annihilator of its
kernel under the perfect dot-product pairing. This follows directly by
Smith reduction, or by containment and equality of their finite orders.
The displayed congruence is soluble precisely when `w^T b` is even for
every w in L_j. That is exactly `q_j(w mod 2, x mod 2)=0`. It proves the
radical identity. Finally Smith reduction is invertible modulo two: a
coordinate of its transformed x can have arbitrary parity precisely when
the corresponding Smith factor is divisible by 2^j. This proves the
dimension and rank formulas. Nonsingularity bounds every valuation. QED.

For a connected graph, any reduced Laplacian is positive definite and
nonsingular, so the theorem applies without additional hypotheses. It is
an all-size electrical rule, not an all-order formula for tournament
metagraphs. For G=L+2B the initial kernel is that of L modulo two; the first
divided form is where B first becomes visible. Subsequent layers retain
the information that reducing modulo two had hidden.

## Order eight: the persistent mode dies at level eight

There are 3,528 merged vertices. Blue support uses 176 vertices and is
connected on that support. Thus the full blue graph has 3,353 components;
its binary bicycle space has dimension 11. The full graph is connected.
Its reduced binary kernel therefore has dimension

    (3353 - 1) + 11 = 3363.

The exact divided-energy computation is:

| j | dim V_j | rank q_j | dim V_(j+1) |
|---|---:|---:|---:|
| 1 | 3363 | 3358 | 5 |
| 2 | 5 | 4 | 1 |
| 3 | 1 | 0 | 1 |
| 4 | 1 | 0 | 1 |
| 5 | 1 | 0 | 1 |
| 6 | 1 | 0 | 1 |
| 7 | 1 | 0 | 1 |
| 8 | 1 | 1 | 0 |

Hence the 2-primary critical group is exactly

    (Z/2)^3358 direct-sum (Z/4)^4 direct-sum Z/256.

The spanning-tree count has 2-adic valuation `3358 + 2*4 + 8 = 3374`.
The actual tree count and odd primary factors are not computed here.
Through order seven the inherited full groups only had factors of order
two and four. Thus order eight is a concrete obstruction to extrapolating
that small-order exponent bound. The remaining one-dimensional mode is
not a free choice: its authenticated lifts survive five further zero
forms before a nonzero form finally removes it.

## Odd conductance scaling: a small-order pattern breaks

Keep the exact same support and blue conductances, and replace G by
`G_b=L+2bB` for positive odd b. This is a weighted-network probe, not a
claim that the rescaled graph is a different native tournament metagraph.
The first form is unchanged: the blue half-energy is fixed, and the black
contribution is multiplied by b, which is one modulo two. Thus V_1 and
q_1, in particular the number of factors divisible by four, are identical
for every odd b. Higher integer lifts need not be identical.

The exact `odd_scaling.py` sidecar finds:

| b | Common factors | Remaining factor | Tree-count valuation |
|---|---|---|---:|
| 1 | (Z/2)^3358 + (Z/4)^4 | Z/256 | 3374 |
| 3 | (Z/2)^3358 + (Z/4)^4 | Z/8 | 3369 |
| 5 | (Z/2)^3358 + (Z/4)^4 | Z/16 | 3370 |
| 7 | (Z/2)^3358 + (Z/4)^4 | Z/8 | 3369 |

All four filtrations terminate. They agree even at the second rank in
this sample, but that extra agreement is not inferred from the first-form
argument. This provides a finite counterexample to extrapolating the
through-order-seven invariance of the entire 2-primary group under odd
color rescaling. It retains the valid all-order first-form invariance.
The higher exponents measure an arithmetic interaction between the two
colored networks that their parity reduction cannot recover.

## Reproduction and independent checks

```text
python -B 05-knowledge/results/metagraph_order8_electrical_20261023/check.py
python -B -O 05-knowledge/results/metagraph_order8_electrical_20261023/check.py
python -B 05-knowledge/results/metagraph_order8_electrical_20261023/odd_scaling.py
```

The consumer reads the full merged node/edge exports and checks their
normalized hashes. It reconstructs the blue binary kernel, forms the full
bilinear half-energy, and stores literal mod-4, mod-8, ..., mod-256 lifts.
Every lift is checked by applying the actual sparse integer Laplacian;
every divided energy is computed by an integer dot product. The small
order-three through order-seven networks provide inherited controls.
No diagonal-only test or floating-point eigensolver is used.

The saved receipt has 439,932 generating checks; ordinary/optimized replay
has 439,933 including equality to the receipt. Its LF SHA256 is
`417e4cc3af2688df39bbb71701789a6967677ae82ea8d0f0badb8ff50d739b27`.
The [separate audit](../metagraph_allsizes_20261023/DYADIC-AUDIT.md) supplies
independent small symmetric-matrix controls and checks the higher-lift
certificate. In particular scaled hyperbolic blocks have zero diagonal
but nonzero bilinear energy, so they reject a diagonal-only implementation.

The discovery program has a level-13 safety cap. It reports whether the
filtration is complete and never treats a surviving radical at the cap as
zero. The present graph terminates at level eight. A general-purpose use
may increase the cap, or use a determinant bound; the theorem itself has
no fixed-depth restriction.

## Next precise boundary

Port signatures and electrical forms remove different ambiguities. The
former retain local directed incidence and resolve the finite port
transport problem. The latter retain global cycle/cut interactions and
can reject rewiring that local margins permit. Electrical forms provide
an independent rejection test; they do not reconstruct the ports. Once
full port wiring is known, these forms are determined by it. The useful
combined test retains local port constraints and bilinear energies when
proposing smaller realizability certificates.
Proving that some bounded signature determines every order remains open.
