# Full complement-word costs reduce the surviving clock ceiling to 16,704

**Status: PROVED analytic transfer and FINITE-EXACT complete certificate;
INDEPENDENTLY AUDITED.**
Restoring the full seven-state complement words removes
99 more scales from the audited finite grid baseline. The refined necessary
set has **8,202 scales**, with maximum **16,704**. Every scale is retained
explicitly in the output. LRC(14) remains open; membership in this finite
set does not assert a realizable or unsafe physical row.

The previous [audited bootstrap](third_20260906_grid_bootstrap.md) remains a
separate frozen baseline. This note consumes the independently computed
[full-word maximum-cost theorem](third_20260906_grid_full_words.md); it
does not replace its optimization proof with selected maximizing examples.

## 1. Inheritance and the new coordinate

The parent [translated-grid theorem](third_20260906_grid.md) bounds the
total seven-tail ceiling excess and recovers actual overlap credit. Its
lower-dimensional lonely-runner phase input remains **CITED** as documented
there. The [small sheet-gcd edge theorem](third_20260906_decoder_profile_graph.md)
supplies some actual edge with `e<=6`, while the actual primitive ratio of
that edge remains in the entire5,855-pair **THM-3818** atlas,
[scaled inert cube-class support-two packet](../../01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md).

The audited bootstrap preserves each ratio's complete clipped component
lengths, every divisor `e|t`, and the exact forced margins

```text
d_p=e*gcd(t/e,p),       d_q=e*gcd(t/e,q).
```

Its marginal cost upper bound still maximized five remaining values using
separate multiplicity capacities. The full complement words are the
coordinate restored here. The new source-to-target map is a complete
finite profile optimization to an upper bound on actual bad-incidence
excess, followed by the already audited physical phase-lifting consumer.

The board is: finite clock divisors; whole complement words; maximum total
ceiling excess; a selected actual pair; separate component rounding; and
compatibility of the maximizing word with that pair. The first probe was
the former largest clock23760: its relaxed excess252 falls to a full-word
maximum168. This positive signal led to the complete8,301-clock pass.

## 2. Exactly which full-profile optimization is consumed

For each baseline clock t, let `F_t` contain all seven-state multisets
`d=(d_1,...,d_7)` satisfying:

* Every d_i is a divisor of t in the inherited42-value seven-body set.
* `gcd(d)=1`.
* For every nonempty subset I of k positions, `1<=k<=6`, put
  `c=gcd(d_i:i in I)`. The pair

  ```text
  (c, sort(gcd(c,d_j):j not in I))
  ```

  belongs to inherited profile level `7-k`.

These are the complete projected words from the six physical V labels
together with selected U labels. They do not assert every other physical
profile, the actual ratios among all seven U labels, or an actual decoder
realization. Repeated states are allowed even when physical speeds are
distinct.

The supplier computes the exact finite maximum

```text
M(t)=max_(d in F_t) [sum_i d_i*ceil(t/(7d_i))-t].     (1)
```

Every actual hypothetical failure has its seven sheet states in F_t,
so its true ceiling excess E satisfies `E<=M(t)`. The exact upper bounds
come from the supplier's complete prefix-profile optimization; every
retained maximizing owner proves attainment in F_t. The consumer directly
rechecks all126 nonempty proper-subset words of every owner, all seven
clock divisibilities, and the literal ceiling cost. Such attainment checks
alone are not substituted for the separate optimization proof.

The input is exactly the audited baseline's8,301 clocks, whose compact
JSON array has SHA-256

```text
a25d83f0eeb630bb82e84cdfac4e3cf7312f892f6c426d6affd5239a064e4b58.
```

The compact array of all ordered `[t,M(t)]` rows is pinned by

```text
ca6b6f562db1fc3632f8b7570b89a16020a981ae8aa130be200dc1bdcb4264ca.
```

Owner order and search-node counts are diagnostic data, not part of the
maximum-cost semantics hash. The consumer requires an exact clock match
and a successfully completed supplier replay.

## 3. The complete refined consumer

For each t, every divisor `e|t` with `e<=6`, and every actual coprime
atlas pair `(p,q)`, retain the baseline's full component credit C and
forced-pair cost bound `E_pair`. Reject pairs with invalid forced margins
or exhausted multiplicity capacity. The new necessary relaxed test is

```text
C(t,e,p,q) <= min(M(t), E_pair(t,e,p,q)).             (2)
```

The two upper bounds in (2) need not have the same maximizing word. Their
minimum is still a valid upper bound on the actual excess, so this is a
sound relaxation. It does not solve the joint full-profile optimization
conditioned on every selected pair.

The exact two-profile aggregate envelope from the baseline is used only
as a cheap exclusion filter. All required ratios are retained for the
individual-component computation. The three complete output sets are:

| Consumer on the baseline8,301 clocks | Count | Maximum |
|---|---:|---:|
| Full-word M(t), aggregate component ceiling | 8,288 | 21,600 |
| Full-word M(t), individual component ceilings | **8,202** | **16,704** |
| Individual ceilings against `min(M(t),E_pair)` | **8,202** | **16,704** |

In this finite universe the last additional minimum removes no further
clock. The output records the full arrays and all99 exclusions from the
baseline; no sampled-clock argument is used. The only survivor exceeding
14,904 is **16,704**.

For physical transfer, suppose a balanced actual decoder row fails weak
safety. The inherited theorem puts t in the baseline set, and the
small-edge theorem supplies some actual edge with e<=6. Its sheet states
have excess at most both bounds in (2). If t is absent from the refined
set, every candidate for that actual edge has overlap credit strictly
exceeding the total ceiling excess. The inherited translated-grid count
then produces a lift of a six-body safe phase that is weak-safe for all
seven remaining labels, a contradiction. Thus every hypothetical failure
in the stated actual decoder domain satisfies

```text
t in S_refined,    |S_refined|=8202,    max S_refined=16704,
1<=g<=90.                                           (3)
```

This is an exact finite necessary scale set, not a finite enumeration of
all physical component shapes. Strict physical clearance is not inferred
from the weak grid endpoints.

**General connected-complement corollary.** Let n be any primitive row of
thirteen distinct positive integer speeds. Choose any six-label subset A,
and let B be its seven-label complement. Put `t=gcd(A)`, `g=gcd(B)`,
`V=A/t`, and `U=B/g`. Suppose the graph on U is connected when an edge
means that its coprime positive ratio `(p,q)` belongs to the strict
5,855-pair atlas. Then the row is weakly safe whenever

```text
t not in S_refined; in particular, whenever gcd(A)>16704.   (4)
```

Indeed primitiveness gives `gcd(V)=gcd(U)=gcd(t,g)=1`. Under hypothetical
failure, the inherited joint-shadow profiles apply globally to the
thirteen-speed row, so their projected seven-state words, `g<=90`, and
the initial grid ceiling all hold. Connectedness of the stated U graph
is precisely the graph hypothesis used by the small-edge certificate;
it supplies some edge with e<=6. The rest of the cost and lift argument
is unchanged. This corollary requires neither decoder-span equality, a
physical height box, nor connectedness of the chosen six-label set.
It claims existence of a weak-safe phase, not safety on every grid lift.

## 4. A failed final exclusion supplies an honest stopping witness

The first retained pair found at t16704 had forced margins `(6,12)`.
Conditioning the full-profile optimizer on those values lowers its
maximum to133, so that particular candidate disappears. It is invalid
to infer that every actual pair disappears.

The alternative pair

```text
t=16704, e=4, (p,q)=(3,308),       p+q=311,
(d_p,d_q)=(12,16)
```

is in the actual inert atlas and is compatible with the full valid word

```text
d=(12,16,72,58,64,9,9).
```

All seven values divide t; every inherited projected complement word
passes. Its literal excess is188, exactly the unconditional maximum
`M(16704)`. Since it contains the two forced margins, the pair-conditioned
full maximum is also188. The exact component credit is only172. Thus

```text
C=172 < 188=M(16704)
```

even after this pair compatibility is imposed. The supplier and consumer
both retain the explicit word, forced margins, atlas membership, and
component calculation. This refutes the proposed implication that
conditioning the first candidate suffices to lower the ceiling to14,904.
The strongest surviving result is (3).

This witness is a compatible word and selected pair in the declared
finite quotient. It is not a complete connected seven-component
realization, a strict failure, or a claim that the actual scale bound
cannot improve. The unresolved information includes the coherent actual
ratios among the other labels and simultaneous attainability of the
separate arc counts.

## 5. Reproduction and verification

The [consumer source](../../04-computation/third_20260906_grid_refined.py)
and [retained output](third_20260906_grid_refined.out) reuse the pinned,
independently audited baseline pair compiler and consume the separately
audited maximum table. The baseline source bytes are pinned; the maximum
table is pinned by its full mathematical `[t,M(t)]` sequence. Independent
verification uses a separate literal-interval and bounded-knapsack engine,
not these producer imports.

Each `SCALE_SET` line is a JSON object with `name`, `count`, `maximum`,
the full sorted `survivors` array, and its compact-JSON SHA-256. Three
complete refined sets are recorded. The output also gives all99 removed
clocks and the positive boundary above.

```text
python3 -B 04-computation/third_20260906_grid_full_words.py
python3 -B 04-computation/third_20260906_grid_refined.py
python3 -B -O 04-computation/third_20260906_grid_refined.py
```

Normal and optimized consumer outputs agree byte for byte after
**1,467,330 active requirements**. Raw LF-byte SHA-256 values are:

```text
source 10af6948ed59efe5e7cc06953410d73839dc3bdfe280519208a9e0099df912db
output 4113f3180f27702e23d1dbe3529ac1fb47f8a6b4795f7003034048d589c03447
```

The final `full_words_coupled` array has compact-JSON SHA-256
`4f29481d984ead40d0144556ce1c45dce210e30b964bb65835a7904ca6353e59`.
The [independent consumer audit](third_20260906_grid_refined_audit.md)
reconstructs all three arrays by literal interval intersections and a
separate bounded-knapsack engine, accepting the transfer relative to the
pinned maximum table. The
[separate optimizer audit](third_20260906_grid_full_words_audit.md)
independently verifies all8,301 maxima, the global residue envelope and
both conditioned boundary controls. The optimizer's status is not inferred
from the downstream count.
