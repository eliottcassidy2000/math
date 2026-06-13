# HYP-2211: Rationality Is a Two-Shadow Carrier Obstruction, Parallel to Parity and Add/Multiply Quotients

**Status:** OPEN synthesis, with S635 finite-field shadow audit and a standard
algebraic lemma.

## Correction

The working fact is not:

```text
exactly one of pi + e and pi*e is rational.
```

The known fact is:

```text
pi + e and pi*e cannot both be algebraic.
```

Therefore at least one of them is transcendental and hence irrational.  It is
not known which one, and it is not known whether both are transcendental.

This correction matters because the useful structure is not an exclusive
choice.  It is a side-channel obstruction: two tame shadows would reconstruct
an object known not to be tame.

## Claim

Rationality/irrationality should be treated as a field-of-definition shadow in
the same family as:

- odd/even parity as descent to `F_2`;
- sum/product as complementary elementary symmetric coordinates;
- LRC reset periods plus relation lattices;
- unit-distance edge count plus spine/bulk carrier;
- tournament fixed/merged/nonfixed layers under complement/converse.

For a field `K`, if

```text
S = x + y in K,
P = xy in K,
```

then `x` and `y` are roots of

```text
T^2 - S*T + P in K[T].
```

Thus the unordered pair `{x,y}` descends to a quadratic carrier over `K`.
Contrapositively, if the hidden pair contains a coordinate known
transcendental over `K`, then the two shadows `S` and `P` cannot both be
`K`-algebraic.

With `K = Qbar` and `{x,y} = {e, pi}`, this proves `e+pi` and `e*pi` cannot
both be algebraic.  The scalar uncertainty is not noise; it is exactly the
signal that the missing field-of-definition side channel has to live in at
least one shadow.

## S635 Evidence

S635 adds `04-computation/rational_shadow_carrier_s635.py` and stores
`05-knowledge/results/rational_shadow_carrier_s635.out`.

The finite-field microscope checks unordered pairs over `F_p`.  Sum alone and
product alone have nontrivial fibers, but the joint `(sum, product)` carrier
has max fiber `1` for every checked prime:

```text
p=5:  sum max fiber 3, product max fiber 5, joint max fiber 1
p=7:  sum max fiber 4, product max fiber 7, joint max fiber 1
p=11: sum max fiber 6, product max fiber 11, joint max fiber 1
p=17: sum max fiber 9, product max fiber 17, joint max fiber 1
```

This is the finite version of Vieta's reconstruction: the pair is not seen by
either one-shadow, but it is seen by the two-shadow carrier.

The S635 proof-lens Tournament Analysis chooses carrier lenses as vertices,
not constants, runners, or tournament vertices.  Its majority tournament is
transitive with one Hamiltonian path and ranks:

```text
joint_sum_product_carrier
field_descent_obstruction
lrc_relation_lattice
anti_coset_transporter
sequence_shadow_recursion
unit_spine_bulk_split
parity_redei_shadow
product_shadow_only
sum_shadow_only
```

This ranking is not a theorem about all problems.  It records the immediate
algorithmic lesson: one-shadows are cheap but dangerous; the joint carrier is
what prevents false scalar collapse.

## LRC Reading

For LRC speed sets, the analogous first split is:

```text
reset period / commensurability  +  relation lattice / short circuits.
```

If speeds are rationally commensurable, the orbit has a finite reset period
and "hit the safe box once" is a finite periodic question after quotienting.
If the speed ratios have irrational rank, the orbit is dense in a torus
modulo its relation lattice; then the relevant quotient is the relation module,
not a literal all-runner reset clock.

This suggests a safer preprocessing rule:

1. factor out common scale;
2. compute rational rank / reset period when it exists;
3. compute short additive circuits in the relation lattice;
4. only then choose the time-search model.

The pi/e obstruction is the warning sign: if two very compressed shadows would
force the hidden state into a tame carrier, but the hidden state is known not
tame, at least one shadow must be carrying the hard part.

## Unit-Distance And Tournament Reading

For unit distance, `57` at `n=21` is not a bare scalar.  HYP-2210 records the
carrier split:

```text
57 = 20 spine edges + 37 bulk edges = 20 + C_hex(3).
```

For self-converse tournaments, the raw complement merge is similarly too
coarse.  The retained data are rooted perspectives, fixed layers, merged
layers, nonfixed pairs, and anti-coset transporters.

These are the same pattern as `(e+pi, e*pi)`: one scalar may be unknown or
misleading, while a two-shadow carrier explains what information is forced to
survive.

## Assumption Challenge

Alternate vertices considered: constants, runners, speed ratios, relation
lattice vectors, reset clocks, proof obligations, field descents, parity
quotients, unit spines, bulk shells, rooted perspectives, anti-cosets, and
sequence-shadow methods.

S635 chooses proof obligations and carrier choices as vertices.

Preserved predicate: whether a quotient can still certify that a target
obstruction is impossible or must survive in at least one shadow.

Destroyed data: representatives, embeddings, exact times, individual labels,
and which of the two shadows contains the hard coordinate.

The challenged assumption is that rational/irrational is a terminal
classification.  Here it is a diagnostic that a quotient has reached a
field-of-definition boundary.

## Next Tests

1. Add a reusable `two_shadow_obstruction` helper: if two shadows together
   reconstruct a forbidden/tame hidden carrier, then they cannot both be tame.
2. Apply the helper to LRC preprocessing: `(reset period, relation lattice,
   circuit support)` before dense time checking.
3. Extend HYP-2209 sequence-shadow work from class counts to strong-`H`
   spectra: fixed/merged/nonfixed layers around the strong-component norm.
4. Add `(spine, bulk, transporter)` fields to unit-distance beams so raw edge
   counts are never the only state key.

## See Also

`04-computation/rational_shadow_carrier_s635.py`;
`05-knowledge/results/rational_shadow_carrier_s635.out`;
`07-reflections/rational-shadow-carriers-s635.md`;
HYP-2210; HYP-2209; HYP-2208; HYP-2207; HYP-2154; HYP-2155;
`07-reflections/natural-operation-digraphs-and-product-sum-s365.md`;
`07-reflections/the-summand-graph-node-is-the-lrc-pinch-denominator-addition-and-multiplication-shadows-s559o.md`.
