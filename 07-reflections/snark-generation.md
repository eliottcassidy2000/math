> **HISTORICAL SPECULATION / SUPERSEDED (2026-08-25).**  No unique
> factorization, Euclid operation, four-level chain, or shared `{2,3,7}`
> obstruction was proved.  A repeated numeral is not a transfer map.
> MISTAKE-507 withdraws those claims and this file records the repaired
> research direction.

# From “Snarks as Primes” to Boundary Obstruction Compilers

## Why the arithmetic analogy stopped

Graph dot products, 2-sums, and other gluing operations depend on marked
interfaces and need not give unique factors.  “Prime snark” also varies with
the chosen reduction operation.  Consequently, Petersen, the Blanusa
snarks, and flower snarks cannot be assigned arithmetic prime/composite
values without first defining and proving a factorization theory.  Local
modifications such as Loupekine constructions do not play the role of
`product+1`, and no four-level termination law is known.

## Exact replacement: an interface relation

Cut a graph along an ordered edge boundary.  For each piece `X`, retain

```text
f_X(sigma) = number of internal proper colorings extending boundary word sigma.
```

THM-4116 proves the exact gluing law

```text
#Col_3(G) = sum_sigma f_X(sigma) f_(V-X)(sigma).
```

The Petersen hostile has two nonempty boundary supports with zero overlap;
the analogous `K_4` control has dot product six.  Thus noncolorability is not
the absence of local states but an incompatibility between two extension
relations.  This is the reusable obstruction compiler that the prime analogy
was trying to locate.

## Past snark data that remains useful

- THM-261 identifies Petersen exactly as the disjoint-support graph on
  positive `A_4` roots.  It preserves support orthogonality, but loses arc
  orientation, signs, and the distinction between the `so(5)` and `sl(5)`
  vectors sharing the same pair index.
- Full cycle-length profiles show that high girth is not cycle poverty:
  Petersen has `(c_5,c_6,c_8,c_9)=(12,10,15,20)`.  Witness cycles and their
  cut incidences are needed in addition to the counts.
- Resistance, oddness, criticality, and boundary-extension rank measure
  different repair costs.  A single “is a snark” bit destroys all of them.

## Exact adjacent-edge atlas and next experiment

THM-4116 now builds the adjacent-edge atlas for the colorable `K_4` control,
Petersen, `J_5`, and both Sage-convention Blanusa snarks.  Parity gives four
boundary-word orbits; for every uncolorable finite simple cubic graph, a
Kempe involution compresses them to `(m_e,m_e,0,0)` with even `m_e`.  The
exact histograms are
`{2:15}`, `{4:5,6:10,10:10,12:5}`, `{4:23,6:4}`, and `{4:25,6:2}`.

The next experiment should compare matched six-pole cuts and record:

1. support size and extension multiplicities beyond the forced four-pole law;
2. compatibility rank after quotienting the global `S_3` color gauge;
3. the smallest vertex or edge deletion making the pairing nonzero; and
4. cycle witnesses meeting the cut.

The experiment can reveal reusable obstruction types without asserting an
arithmetic factorization that the graph operations do not preserve.
