# Changing the operation changes the exact quotient: all regular order-seven triangle reversals

Status: FINITE-EXACT exhaustive finite classification and transition certificate; PROVED consequences for its finite Markov chain. This is not a new classification of regular order-seven tournaments. The inherited three classes and Hamiltonian counts are rechecked, then used to expose the boundary of the fixed-Fano operation model.

## Why this is the decisive next hostile

`FANO_CONTINUATION_SPECTRUM.md` uses seven specified Fano-line reversals. Every tournament in that chart has14 directed triangles, so restricting to7 lines omits half its regularity-preserving elementary reversals. An apparent universal eight-state rule must be tested against that missing operation before it is transferred to tournaments generally.

Reversing a directed triangle changes one outgoing and one incoming arc at each of its vertices, preserving every score. For regular order7, every vertex has score3, and the number of directed triangles is C(7,3)−7C(3,2)=14. Thus uniformly choosing a directed triangle is a well-defined fourteen-neighbor random walk on the full labelled regular-tournament class.

## Exact three-class transition

In the order (Paley P, Fano-defect B, rotational R), Hamiltonian counts are (189,171,175), labelled orbit sizes are (240,1680,720), and the number of elementary reversals leading to each class is

Q = [[0,14,0], [2,9,3], [0,7,7]].

These are neighbor counts, not probabilities; divide by14 for the transition matrix. Diagonal entries count different labelled tournaments isomorphic to the starting tournament, not identity moves.

The three representatives are encoded by masks1956555,1956430,2079175 on lexicographically ordered pairs of {0,...,6}, where bit1 directs the lower label toward the upper. P is the circulant with steps{1,2,4}; R has steps{1,2,3}; B is one directed-triangle reversal of P.

### Certificate of completeness, rather than a guessed three-state model

The script constructs each full permutation orbit, of the sizes above, and checks they are disjoint. Independently a score-pruned recursive generator constructs **every** labelled tournament with all scores3, without assuming the orbit classification; its set equals their union of2640 masks. Every one of the2640 masks has all14 reversals inspected and the same class-neighbor counts as its representative. Thus no unrepresented regular class is silently assumed absent.

Hamiltonian counts of the three representatives are computed independently by subset dynamic programming and by checking all5040 vertex words. Isomorphism transports these counts throughout their orbits.

The fixed-Fano chart contains only P and B. For instance the certificate gives a chart state with mask848974 and a directed-triangle reversal to849223 of typeR, outside the chart. All three such R-neighbors of any B state are omitted from its fixed-Fano moves. The syndrome chart is therefore not closed under all triangle reversals.

## Exact stochastic consequences and what they do not imply

The quotient is reversible with stationary distribution (1,7,3)/11, agreeing with labelled orbit sizes. Detailed balance includes240·14=1680·2 and1680·3=720·7.

Its characteristic polynomial is (x−14)(x²−2x−14), so its nonconstant transition eigenvalues are(1±sqrt15)/14. These are quotient eigenvalues; this is not the full2640-state spectrum or a claim about its spectral gap.

Mean first hitting times of the Paley **class** are10 fromB and12 fromR. Indeed

14 t_B=14+9t_B+3t_R,
14 t_R=14+7t_B+7t_R.

The mean first return to the Paley class starting there is11, not the8 of the restricted Fano walk. The stationary mean Hamiltonian count is1911/11, not693/4. The selected operation and its stationary measure have changed; there is no contradiction.

At the isomorphism-class level H itself determines the three states, so H is sufficient for this randomized walk. It still does not supply a complete labelled controlled-operation state. Nor is triangle reversal the fixed-path nonpath-complement operation used in the merged metagraph.

## Transferable conclusion

The same object can support several exact, incompatible-looking reductions because they preserve different continuation languages. Before using a symmetry quotient for enumeration, dynamics or wiring, check **closure under the actual allowed moves** and **the measure used for averaging**. A chart can be exact internally and still omit a new isomorphism type accessible by one lawful ambient operation.

Reproduction: `python3 check_regular7_triangle_walk.py > regular7_triangle_walk.json` (also `python3 -O`). 5294 explicit checks; every directed triangle at every labelled regular order7 tournament is inspected. The finite census is small and structurally targeted, not used as evidence for an all-order classification.
