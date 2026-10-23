# Independent proof review of assembly and transport lifts

23 October 2026. Reviewer: native electrical-code lane.

## Binary full-record assembly

**PASS.** Reviewed the full proof in
[`assembly/BINARY_ASSEMBLY.md`](../assembly/BINARY_ASSEMBLY.md), with the
following potentially delicate implications checked independently.

- Literal owner equality and distinct finite resource inventories really
  reduce to unary and binary forbidden state pairs when the supplied
  records are complete. The additional higher-arity guard exclusion is
  essential and explicitly retained.
- A complementary implication-component pair contains each represented
  packet with a fixed sign; different pairs partition the indices.
  Mutual implication makes simultaneous change of all those indices
  necessary. Reachability to the complementary component correctly
  characterizes a forced truth value.
- In the constructive path proof, a current-true/target-false component
  with no predecessor among such components has no retained true
  predecessor. Adding its complement cannot force a still-false distinct
  component, by contraposition. The potential implication from its
  complement back to itself is separately excluded by target feasibility.
  Thus every step is accepted, flips one component pair, and never revisits
  a changed packet. This proves both the exact weighted threshold and the
  weighted-Hamming total cost, without assuming free components are
  independent choices.
- For the native family, the endpoint-fixed complement isomorphism
  survives transitive prefix/suffix padding. The ordered strong-component
  profile distinguishes the chosen sources and their converses under the
  printed bound `n>=2k+7`. The complement partners remain strong while
  every padded source is decomposable, so source and inverse endpoint
  inventories are disjoint. Complementation of normalized words and the
  trivial source automorphism group supply the claimed `8k` reservoir.
  Shifted actual clocks and the separate cost of identified shared
  variables remain explicit.

Normal and optimized receipt replays both pass **180,007** checks on
4,784 declared Boolean controls and 36 native padded hosts. This review
audits the all-size proof; the finite host range is not its replacement.

## Two-sided dyadic transport

**PASS.** Reviewed
[`dyadic/TRANSPORT_LIFTS.md`](../dyadic/TRANSPORT_LIFTS.md).
The zero-margin factorization `M=R C Q^t` is unique and retains the whole
integer matrix. Both source and target equations are necessary. For the
binary correction map `D -> (Du,D^t v)`, the common-corner condition and
the two zero-vector guards are necessary and sufficient, with the stated
counts. The reciprocity identity between divided residuals holds over the
integers before reduction, so it removes the shared-corner obstruction
when the residuals arise from one matrix. The even-vector hostile and the
bounded-box hostile correctly retain the other obstructions.

The scope is unrestricted coefficient lifting at levels `h>=1`, with
literal capacities intersected afterward. Extra binary coordinates are
unchanged by `2^h D`; multiple independent cyclic coordinates would need
a stacked correction problem. No general port allocation or bounded
integer correction follows from the dyadic solution alone.
