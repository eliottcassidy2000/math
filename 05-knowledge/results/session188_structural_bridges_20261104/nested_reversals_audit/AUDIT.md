# Independent audit of nested and exposed-laminar noncancellation

**Verdict: PASS.** I independently derived the weighted averaging identity
and odd/even anchor battery before reading `simultaneous_response/ANCHORED.md`,
then checked the written proof, its operation semantics, and its checker.
The result proves nested-pair and exposed-laminar separation, not unrestricted
laminar reconstruction or arbitrary simultaneous reversal separation.

## 1. Exact weighted averaging

The partially tested Hamilton-word involution is valid: reverse every maximal
contiguous run in the actual module A, leaving all exterior names and both
ports in place. Every tested old-to-module boundary remains valid because
the old exterior incidences and the second port's incidences are uniform.
The x-A contacts are the only deliberately untested edges.

For one-hot x attachments, the number of allowed choices of the marked
physical name is m, 1, m-1, or 1 according as x has no A neighbor, only an
incoming A neighbor, only an outgoing A neighbor, or both. In the last case
the names are distinct even if the A runs are singletons: a Hamilton word
cannot revisit a name. Thus reversal preserves both partial validity and
weight. This proves the integer equality of sums, with no truncation of
cubic or quartic response terms.

The result does not say that the individual nonuniform terms are invariant.
The common outside mask across the entire one-hot block is essential.
The module condition ensures that this mask is indeed common in the anchor
battery because a fixed exterior anchor sees all of A in one direction.

## 2. Battery algebra at every core order

With f(a)=delta_B(a,e_v), substitution of the incoming finite differences
into E1 leaves:

    sum_(u->v) f(e_u) + sum_(v->u) f(e_u+e_v)
       + (d mod 2) f(e_v)
       + ((n-1+d) mod 2) f(0) = 1+b_v.

The b contribution is `sum_(u!=v)b_u=b_v`, since total starting and ending
path counts coincide. All masks in the displayed battery are distinct.
For even n exactly one final query occurs; for odd n either both or neither
occurs. Thus the written n, n-1, and n+1 battery sizes are correct.

The singleton query `(e_v,0)` handles b_v=1. The n=2 boundary is not an
application of E1: it is handled solely by this endpoint query, and both
endpoints have b_v=1. Nested nontrivial A properly inside B automatically
has |B|>=3, so the main nested proof does not use the exceptional order.

## 3. Nested pairs and larger laminar families

For A properly inside B, choose v outside A. The A-indexed part of the
battery is exactly an averaging block. Every other battery term is uniform
on A and hence individually invariant under its internal converse. This
continues to hold after B is conversed: the chosen query masks stay fixed;
their common outside values need not match the new root incidences for the
averaging lemma to apply. Only uniformity of the other port on A is needed.

Lifting the battery to the ambient tournament uses the actual unchanged
exterior star of B. The old B plus the ports is a module on both sides.
The incoming modulo-four substitution law multiplies the core contrast by
the odd Hamilton count of the quotient. Its parity cannot erase the core
certificate. No exact integer cancellation of the quotient is assumed.

For an exposed laminar family, each proper selected module avoids v and
remains a module after every earlier laminar reversal. A containing reversal
reverses its uniform exterior directions; a contained or disjoint reversal
does not destroy uniformity. The anchor's incidence to each proper selected
module stays uniform, so its indexed one-hot block continues to use a common
outside mask. Other maximal selected modules are disjoint from B and are
uniformly seen by both lifted ports; all their reversals are individually
invisible. The successive exact averaging argument is therefore valid.

The exposed-anchor assumption is not cosmetic. A maximal selected module
covered by proper selected modules has no eligible anchor, and the proof
does not silently select one. The six-name whole-plus-three-pairs example
is outside the theorem's coverage despite its separately known separator.

## 4. Independent finite controls

`verify_independent.py` imports no producer module. It verifies the actual
run-reversal bijection term by term in 6,602 configurations: all tournaments
through order three, the explicitly declared ten-code slice at order four,
every nonempty module including singleton controls, and all allowed exterior
and uniform-second-port masks. It checks 172,733 valid partially tested
Hamilton words with their multiplicities.

It independently builds 24 literal nested-pair examples of core orders
three through eight, arbitrary inner tournaments and quotients, and optional
one-vertex ambient quotients. Actual Hamilton-path counts of the complete
two-port attachments verify the exact aggregate invariance, the odd/even
residue, the endpoint case, and an actual separating query. Physical rows,
names encoded by their fixed positions, masks, and counts are retained.
No response coefficient of any degree is omitted in these literal counts.

Normal and optimized executions match byte-for-byte. The receipt pins the
audited `ANCHORED.md`. The finite controls corroborate the universal proof;
they are not a new unrestricted reconstruction census.
