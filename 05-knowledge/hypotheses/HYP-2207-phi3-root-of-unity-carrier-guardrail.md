# HYP-2207: Phi3 Is a Root-of-Unity Carrier Guardrail, Not a Raw Equality

**Status:** OPEN, supported by S630 synthesis and earlier HYP-2181/HYP-2184/
HYP-2204 packets.

## Claim

The shared values

```text
Phi3(2) = 2^2 + 2 + 1 = 7
3*Phi3(2) = 21
```

are real root-of-unity signals, but the signal is a carrier guardrail rather
than a raw equality between problems.

The primitive cube-root angle gives a common vocabulary:

- tournament `H` forbids `7` and `21` as Hamiltonian-path evaluations;
- triangular/Eisenstein unit-distance carriers can show `7` and `21` as
  decomposed edge-count scalars;
- self-complementary tournaments realize the conjugation/complement side by
  anti-cosets `Anti(T)=Iso(T,T^op)`;
- LRC rows see the same prime-3/cube-root burden through `C=2n-1` strata.

So `Phi3` should be used as a test for whether a quotient has retained the
right side channels. If it maps a decomposed edge carrier directly to a
forbidden tournament `H` value, it has collapsed the wrong object.

## Evidence

S630 adds `04-computation/sc_complement_perspective_s630.py` and stores
`05-knowledge/results/sc_complement_perspective_s630.out`.

The script prints the explicit guardrail:

```text
Phi3(2)=7
3*Phi3(2)=21
centered_hex(3)=37
```

The same run verifies the self-complementary side through `n=7`: `H=7` and
`H=21` are absent from exact SC classes through that range, while the
complementing operation is organized by anti-cosets rather than by a single
canonical node-swap map.

On the unit-distance side, the relevant `n=21` row is

```text
57 = 20 spine + 37 bulk,
```

where `37` is the centered-hex shell `3*3*(3+1)+1`. This makes the larger
`n=21` signal a spine/bulk carrier statement, not a forbidden `H=21`
realization.

## Relation To The Claude Snippet

The snippet correctly points at the primitive cube root of unity, but S630
separates three layers that should not be identified:

1. `Phi3(2)` as an arithmetic/cyclotomic scalar;
2. tournament `H` as an OCF/Hamiltonian-path evaluation;
3. unit-distance `u(n)` as an additive geometric edge carrier with spine,
   shell, embedding, and exact-vs-lattice side channels.

The root-of-unity bridge is meaningful precisely because it survives as a
guardrail across those layers. It is not evidence that the layers are
equinumerous.

## Tournament Analysis

S630's proof-lens tournament ranks the anti-coset torsor first, then vertex
perspective orbits, the OCF H-gap guardrail, and the Eisenstein `Phi3` norm.
Raw scalar `7/21` is last. The majority tournament has one directed
`3`-cycle and `3` Hamiltonian paths, showing that SC-duality, unit-distance
carriers, and H-gap guardrails are compatible but not reducible to one total
order.

## Assumption Challenge

Alternate vertices considered: roots of unity, tournament vertices, SC
classes, anti-maps, unit-distance points, unit edges, shell events, OCF cycle
packets, LRC residues, and proof obligations.

Chosen vertices: proof/carrier lenses. They preserve the predicate
"root-of-unity scalar appears only with its carrier side channels retained."
They destroy exact geometric embeddings, complete higher-`n` SC spectra, and
the raw equality relation among visible integers.

The challenged assumption is that the same integer appearing in two problems
means the same object has appeared. The carrier decides the object.

## Next Tests

1. Measure whether long anti-cycles `(6)` and `(6,1)` in SC classes correlate
   with `F(T, omega)` or other cube-root evaluations.
2. Add centered-hex bulk statistics to the Moser beam ledger and compare them
   with Eisenstein lattice shell counts.
3. Look for a genuine carrier-entropy quantity that can compare the
   unit-distance `n^1.014` side with OCF/proof-obligation growth without raw
   exponent numerology.
