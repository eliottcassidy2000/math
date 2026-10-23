# Residual counterfeits and a converse repair at order eight

23 October 2026. **FINITE-EXACT:** explicit alternative transport tables and
literal endpoint-port matchings for the complete order-eight universe.
**PROVED:** the singleton-or-converse sufficient criterion below at every
order. No alternative tournament realization or universal finite signature
classification is claimed.

## 1. The question is projected wiring, not class separation

A row retains source ordinary class, closed-path endpoint scores, color,
and the computed target signature S. A column retains target ordinary
class, open-path endpoint scores and color. The actual complement operation
provides a feasible table. A residual alternating cycle changes a relaxed
table while preserving its row and column margins; it proves a wiring
ambiguity only if its merged colored-edge projection is nonzero.

The packet `PROJECTED_TRANSPORT_UNIQUENESS.md` gives the exact residual-cycle
criterion. This consumer constructs positive-capacity cycles, authenticates
every changed cell, and tests their projection rather than interpreting
stalled capacity peeling as a counterexample.

There are two distinct models:

1. **Ordinary transport:** the stated source-signature buckets and target
   endpoint capacities are preserved.
2. **Converse-authenticated transport:** additionally preserve the actual
   class-level reflection equation
   x(i,j,a,b,color)=x(i*,j*,7-b,7-a,color). Literal black multiplicities must
   remain even, including reflection-fixed cells. In our negative witnesses
   an even stronger, literal equivariant matching of actual ports is given.

For this packet's signatures, let c_v be the number of cyclic triples
through vertex v and e_vw the number through the unordered edge {v,w}.
Write

    E_(p,q)(T) = multiset_v multiset_(v→w) (c_w mod p, e_vw mod q).

An exact coordinate instead of a modulus is indicated by “exact.” The
nested row boundaries and tuple association are retained. The inherited
J4 uses vertex score, cyclic-triple incidence and strong-four incidence,
then their outgoing-neighbor records.

## 2. Four actual residual failures, three surviving converse

| Target menu | Ordinary transport | Converse-authenticated transport |
|---|---|---|
| J4 | Nonunique | Nonunique |
| E_(exact,2) | Nonunique | Nonunique |
| E_(2,exact) | Nonunique | Nonunique |
| E_(3,exact) | Nonunique | Unique ordinary target-class table |
| E_(3,3) | Nonunique | Unique ordinary target-class table |

Every negative entry is certified by a literal swap, not only a collision
between class signatures. Each ordinary witness is a four-cell unit cycle.
Each converse-symmetric witness adds the reflected swap, producing eight
changed ordinary cells and a nonzero change of ±2 on four merged black
edge counts. The symmetric J4 and E_(2,exact) witnesses can be the same
swap. Their directed merged changes are

    (2345,1005): -2; (2345,1006): +2;
    (2348,1005): +2; (2348,1006): -2.

The E_(exact,2) symmetric witness changes

    (1058,1751): -2; (1058,2752): +2;
    (1721,1751): +2; (1721,2752): -2.

The undirected changes are nonzero as well; the receipt lists them
separately. These are the sorted-canonical node IDs of the full n8 wiring
packet, not historical graph IDs.

The exact cells, their before/delta/after counts, class canonical codes,
residual cycles and literal masks are in `residual_counterfeits.json`.
Every witness is checked to preserve:

- every source endpoint/color/target-signature bucket;
- every target endpoint/color capacity;
- all blue wiring and total line mass;
- the pointwise triangle transport equation
  c3(target)-c3(source)=start_score-end_score-1;
- the complete merged c3 electrical-current vector.

The symmetric witnesses additionally preserve full converse reflection,
every literal merged black parity, and evenness of reflection-fixed black
cells. Their literal versions take two real closed masks m1,m2 and exchange
the open partners m1^F,m2^F, then do the reflected exchange. The checker
verifies a bijection of real open ports, endpoint/color compatibility,
equality of the actual transported signatures, and commutation with
reflection. Every other literal port keeps its native partner. The resulting table
delta is exactly the recorded one.

These alternate matchings deliberately do **not** send each closed mask m
to its actual complement m^F. Thus they are rigorous counterfeits of a
lossy port relaxation, not a second native tournament metagraph.

## 3. An all-size sufficient rule using an unused orientation

Let C be any class set with converse involution j→j*. Let source categories
r carry an involution r→r*, and suppose the unknown nonnegative table
satisfies x(r,j)=x(r*,j*). Here r omits the target signature: it retains
source class, endpoint type and color; the signature labels a separate
bucket. The supplied source buckets retain

    R(r,s) = sum_(S(j)=s) x(r,j).

**Singleton-or-converse criterion.** If every j in C has either a globally
singleton S-signature or a converse j* with a globally singleton
S-signature, then all entries x(r,j) are uniquely determined by these
bucket margins and reflection.

Proof: in the first case x(r,j)=R(r,S(j)). In the second,
x(r,j)=x(r*,j*)=R(r*,S(j*)). These equations determine every entry.
Native feasibility supplies consistency. Endpoint/color restrictions and
target capacities may be retained but are not needed for this uniqueness
argument. This proves uniqueness of aggregated target-class counts, not
individual pairings of indistinguishable ports.

A weaker claim would be false: injectivity of the pair
(S(j),S(j*)) alone does not necessarily recover joint counts from separate
margins. A four-class rectangle can have distinct signature pairs and an
alternating circulation preserving both marginals. The criterion above
avoids that lost-association inference.

## 4. Why E_(3,exact) is repaired by converse

There are exactly eight classes whose E_(3,exact) signature is not
singleton, forming four pairs. Their sorted canonical IDs are

    804,814; 3022,3031; 3140,3143; 5092,5101.

Their converse classes are respectively

    2910,3095; 1307,1283; 584,605; 953,907.

Every one of these eight converse signatures is globally singleton.
Consequently 6,872 classes are directly forced and the remaining eight
are forced by their reflected singleton buckets. This certifies the
entire converse-symmetric ordinary target-class table, and hence its
colored merged projection, without an exhaustive integer-program search.

The ordinary four-cell witness still exists. Simply appending its reflected
copy fails: it changes source-signature buckets whose target signature
identifies a unique converse class. The receipt explicitly records these
violated reflected bucket margins. This is a genuine symmetry repair,
not an unsuccessful switch search misreported as uniqueness.

Independent recomputation on all 6,880 classes shows that E_(3,3) and
E_(3,exact) induce precisely the same partition. The checker compares
class sets in every block, not just the two signature counts. Therefore
the ordinary counterfeit and the converse uniqueness proof both carry
over to E_(3,3).

## 5. The exact nine-symbol minimum, with its scope

Within this rectangular residue family and with full converse constraints,
the symmetric E_(2,exact) counterfeit excludes every p=2 menu, regardless
of q, and also p=1 by coarsening. The symmetric E_(exact,2) counterfeit
excludes every q=2 menu, regardless of p, and also q=1. Therefore p≥3 and
q≥3 are necessary. The proved E_(3,3) success attains the product pq=9.

Thus **nine symbols attain the minimum alphabet product in this specified
family for converse-authenticated order-eight wiring**. This differs from
the unreflected ordinary transport problem. It does not assert that every
larger pair of moduli succeeds: residue moduli are not monotone refinements.
It does not claim a universal smallest invariant or any all-order success
of this nine-symbol menu.

## 6. Reproduction and dependencies

From this directory in the session workspace:

    python -B check_residual_counterfeits.py
    python -O -B check_residual_counterfeits.py

For relocated published packets, supply the sibling packet directories:

    python -B check_residual_counterfeits.py --incidence PATH --wiring PATH

The incidence directory must contain `check.py` and
`residual_J4.json`, `residual_E_edge2.json`, `residual_E_neighbor2.json`,
`residual_E_neighbor3.json`. The wiring directory must contain
`classes.tsv.gz`, `ports.tsv.gz`, `endpoint_capacities.tsv.gz`, and
`class_of_mask.u16le.gz`. No uncompressed binary or external solver is
required. All input hashes use stable logical filenames independent of
where these directories are installed.

Default runs reconstruct the witnesses and compare the frozen receipt.
Use `--write` only to regenerate it intentionally. The literal mask decoder,
canonical numbering, full class-coverage proof and compressed-ledger
certificate are documented in the full wiring packet. Independent review
checks the negative witnesses and the singleton-converse repair separately.
