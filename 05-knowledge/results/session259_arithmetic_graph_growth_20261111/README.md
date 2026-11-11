# Session259: growing the two arithmetic graphs

**2026-11-11. Status: COMPLETE.**

Start with the [complete synthesis](SYNTHESIS.md) and [one-node birth atlas](BIRTH_ATLAS.csv). This session follows the owner's
summand/multiplicand question, prioritizing new constructions and connections
rather than repeating audits of inherited work. New statements carry proofs;
finite experiments remain separate. No external novelty claim is made.

## Objects and conventions

The ground nodes are positive integers, inserted in increasing numerical order.
An operation event remembers **both parents**, not merely either incoming arc.
For a target n, write

- A_n = {{a,b}: 1 <= a < b, a+b=n};
- M_n = {{a,b}: 2 <= a < b, ab=n}.

The weak variants permit a=b. Restoring that diagonal adds one event precisely
at a double n=2a or a square n=a². Multiplication excludes the unit parent1:
otherwise {1,n} is a same-time event and its unpaired shadow includes loops.
Node1 may be retained as an isolated ground node in the unit-deleted version.

These definitions inherit [THM-2422](../../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md).
That theorem already proves the basic fibre counts, the special {2,3} closure,
and the swap-fixed explanation of parity versus squares. This session does
not claim those as new.

Three distinct clocks must remain visible:

1. **Parent arrival:** {a,b} becomes available at max(a,b).
2. **Target arrival:** its sum event enters at a+b, its proper product event at ab.
3. **Synchronous generation:** start from chosen seeds and apply all currently
   legal operations once per generation. This is not the numerical-node clock.

A fourth object, the **existential shadow**, forgets the pair relation. With
equal parents allowed its additive version is ordinary order, whereas its
multiplicative version is divisibility. [THM-362](../../../01-canon/theorems/THM-362-natural-operation-graph-shadows.md)
is inherited for exactly that weak/simple convention. The strict shadow deletes
a->2a or a->a² respectively; transitive closure can repair some but not all such
arcs. Never silently identify these objects.

## Research board

| Lane | Exact object | Main question / retained coordinate |
|---|---|---|
| Summand | paired-parent birth fibres and two-seed closure | numerical birth vs synchronous frontier, holes and diagonal restoration |
| Multiplicand | proper factor fibres and divisibility complexes | prime exponents, square vs squarefree births, actual persistence maps |
| Dichotomies | swap action and joint (sum,product) packet | fixed points, higher arity, Frobenius, lossless transformations |
| Transport | edge a+b -> ab, 2<=a<b | birth delay, affine doubling vs square boundary, directed vs undirected geometry |

**Inheritance:** closest mechanisms THM-2422, THM-361 and THM-362; canonical
hostiles4 and9 (square atoms disappear if diagonals are deleted); corrected near
miss the mistaken identification of integer addition with addition on log(N);
least-used sidecar the actual parent pair and its two different target clocks.
The prior bridge atlases HYP-3003 and HYP-3004 are inspiration, not universal
transfer theorems. Historical SC-spine enumeration claims are not imported.

The anchor is arithmetic graph growth. The niche is topology of divisor births;
the wildcard is causal sum-product transport and its affine dynamics. Any
Chan–Pak discussion remains a comparison of operation carriers, not a new
completeness claim.

## Files

- [Transport](transport/RESEARCH.md): a causal joint graph and growth laws.
- [Summand](summand/RESEARCH.md): exact seed closure, holes and delayed parity frontier.
- [Multiplicand](multiplicand/RESEARCH.md): full weak and strict reachability barcodes.
- [Dichotomies](dichotomies/RESEARCH.md): joint charts, higher arity and insertion laws.
- [Focused review](review/REVIEW.md), [final replay](FINAL_REPLAY.json), and SHA256SUMS.
- [Synthesis](SYNTHESIS.md): integrated conclusions, boundaries and next questions.
