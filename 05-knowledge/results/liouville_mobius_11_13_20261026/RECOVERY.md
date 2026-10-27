# Recovering the different roles of eleven and thirteen

The earlier [random recovery](../tournament_composition_classification_20261026/RANDOM_CONNECTIONS.md)
already froze a 24-file random sample. This continuation targets occurrences
of 11/13 and Liouville/Möbius, then follows their proof sources. Numerical
equality alone is not used as a map between subjects. Source hashes and
the reviewed base revision are in `source-pins.json`.

## The same two five vertex tournament classes

[THM-466 / weighted odd-cycle digits](../../../01-canon/theorems/THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts.md)
contains two tournaments with score sequence `(1,2,2,2,3)`, four cyclic
triples, and H=11/13. In row-bit convention (bit j in row i means i beats
j), the historical representatives are

```text
H11: (24,1,3,6,14)      c3=4, c4=4, c5=1
H13: (18,8,3,5,14)     c3=4, c4=3, c5=2.
```

The earlier [metagraph result M4–M5](../tournament_composition_classification_20261026/metagraph/RESULT.md)
obtained H=11/13 by adjoining a singleton to the two half-edge presentations
of the strong four-tournament and taking path complement. In its free-bit
order `(0,2),(0,3),(1,3)`, these are masks2 and5. Exhaustive relabelling
now proves that the resulting classes are exactly the historical pair.
This is a recovery of the same classes, not a new pair with matching counts.

At five vertices no two directed odd cycles can be vertex-disjoint, so
the OCF formula reduces to `H=1+2(c3+c5)`. The missing five-cycle changes
H by2 while preserving the score sequence and c3 count. This does **not**
say that the pair has the same labelled cyclic triples: differing H already
rules that out by the inherited C3 path-count invariance. Counts of triples
and the labelled triple relation are different observers.

The mod4 character distinguishes these values (`chi4(11)=-1`,
`chi4(13)=1`), whereas both Liouville and Möbius give -1 to each. The mod4
sign here reads the parity of the total number of directed odd cycles.
It is not Liouville's prime-multiplicity parity. [A4](ARITHMETIC.md#a4-what-multiplication-contexts-recover)
supplies the actual join contexts that distinguish the two Möbius responses.

## Other typed occurrences

| Earlier source | Meaning of 11 or 13 | Preserved lesson and scope |
|---|---|---|
| [THM-519 / Alcuin conflict graph](../../../01-canon/theorems/THM-519-alcuin-number-of-the-conflict-graph.md), section B | H=11 first realizes K5 as an odd-cycle conflict graph | At five vertices all directed odd cycles pairwise intersect; the pair above gives K5 and K6. Both are nonplanar. No general Alcuin conjecture is needed. |
| [THM-383 / half-turn boundary compactification](../../../01-canon/theorems/THM-383-lrc-half-turn-boundary-compactification-audit.md) | H=11 appears in the order5 open-cell clock menu; H=13 does not | This is the inherited finite menu with moving speeds<=8. It is not a global exclusion of H13 or a proof of LRC. |
| [Boundary support metagraph](../boundary_support_20261020_metagraph.md), section4 | Class label11 is transitive; a switched partner has H=13 instead of9 | Here label11 is an index, not H=11. The degree/fibre-size data survive a blue-edge swap while actual wiring fails. Retain the presentation and its complement. |
| [Carousel classification](../cyclictriple_composition_20261026/README.md), Theorem C; [THM-4804 / module shattering and Frobenius](../../../01-canon/theorems/THM-4804-complement-lines-shatter-every-module-and-the-carousel-fibre-as-frobenius-orbits-of-a-normal-basis.md), normal-basis part | 11 and13 are decoration periods, giving186 and630 primitive binary necklaces | The divisor sum is Möbius inversion. The associated four-block tournaments have44 and52 vertices, respectively; period is not tournament size. The new rigid-block sharpness controls instead have55 and65 vertices. |
| [Prime-power bias and the Collatz collector](../prime_power_bias_collatz_collector_20261014.md), sections1–2 | 11 and13 occur as prime moduli in a finite collector test | The relevant recovered arithmetic is the exact square lift of lambda from mu. The collector is an ensemble statement; no single-orbit or global cancellation conclusion transfers. |
| [Fresh re-marking quotient and divisor arrangements](../fresh_remarking_quotient_and_divisor_arrangements_20261026.md), sections4–5 | 11/13 is one of several context-separation witnesses; prime coordinates index a divisor cube | A Boolean cube remembers which primes occur only if its prime labels are retained. The signed disjointness matrix in A4 gives a precise multiplication-response map with an integral inverse. |
| [Euler-product ghosts](../../../07-reflections/euler-product-ghosts-and-zeta-side-cancellation-gates.md) | Coefficient indices11 and13 in the new replay | Exponents, ghosts and ordinary coefficients are different carriers. A3 proves all-order ghost nonvanishing and the product square lift, extending the earlier finite pilot. |

The [CORE-PAPERS](../../reference/CORE-PAPERS.md) route was read before
literature claims. Classical arithmetic is cited to DLMF in A1–A2.
The finite-field correspondence is inherited from the named normal-basis
source; the new computation enumerates binary words and does not claim
to have constructed finite fields or independently reaudited the earlier
million-presentation complement-line census.

## A repaired inversion claim in the recovered theorem

THM-466's congruence is correct, but its former phrase that H's digits
"read off" the collection counts overstates it. Two literal six-vertex
controls are

```text
rows (58,60,57,16,32,8):  alpha=(1,2,1), H=9
rows (48,33,35,39,46,0):  alpha=(1,4,0), H=9.
```

The first is `C3 -> C3`; the second is an H9 five-tournament joined to a
singleton. Directed cycles are enumerated with distinct vertices and a
fixed minimum starting vertex; alpha2 counts disjoint pairs. Independent
subset path counts give9 for both. Size6 is minimal because at size<=5
the only nonconstant coefficient is alpha1. The repaired theorem states
only the forward implication; the full independence polynomial is the
missing coordinate. [MISTAKES.md](../../../01-canon/MISTAKES.md) preserves
the error mechanism and correction lineage. The theorem's H11/13 example
and forward OCF identities remain intact.

## Connection contracts and the concept board

1. **Prime multiplicity to operation responses.** Source: integer prime
   exponents. Target: the observer's response under multiplication. Map:
   `m -> (f(mk))_k`. Lambda preserves parity but loses prime support;
   mu recovers squarefree support and loses every squareful input. Sidecar:
   the actual allowed contexts. Test:121/143/169. For tournaments use only
   attained H factors and lawful order-joins.
2. **Divisor inversion to primitive periods.** Source: words grouped by
   least period. Target: total word counts. Map: divisor summation.
   Preserved: the full sequence when mu inversion is retained. Lost by
   lambda substitution: separation of square-divisor layers. Sidecar:
   period labels, and a chosen normal basis for the field interpretation.
   Test:prime-square periods121/169; false numerator remainder2.
3. **Composition to symmetry.** Source: child isomorphism colours.
   Target: cyclic triples. Map: forget allowed orientations and orders.
   Preserved: path-cover counts; lost: some quotient stabilizers. Sidecar:
   oriented child colours and their internal groups. Test:81/27 at12,
   then3/1 at15. Squarefree sensitivity changes4p to5p.
4. **Cycles to H and ports.** Source: the directed odd-cycle structure
   and actual Hamiltonian presentations. Target: H or coarser score data.
   Evaluation at2 preserves H but loses cycle coefficients; dropping the
   rooted presentation loses composition response. Tests:the recovered
   11/13 pair and the minimal equal-H unequal-profile pair at6.
5. **Euler exponents to ghosts.** Source: arithmetic exponent sequences.
   Target: logarithmic derivative coefficients. Map: the weighted divisor
   transform. Preserved: the full sequence with its weighted inverse;
   not preserved: sparsity merely inferred from another carrier. Sidecar:
   which transform and which coefficient index. Test:two independent
   expansions through degree80, then local-factor proofs for all n.

After the cycle recovery, the board gained a concrete operation map rather
than another occurrence of the same primes. After the ghost calculation,
square dilation became a common exact transform across three different
carriers. After the rigid-child construction, repeated factors became a
size obstruction with a sharp boundary. These are the three productive
cross-connections; no numerical coincidence between H, period and modulus
is promoted to a theorem. The next open problem is to classify observer
responses under the actual path-complement and rooted-gluing operations,
whose allowed contexts are richer than the multiplicative monoid of H.
